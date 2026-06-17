#!/usr/bin/env python3
"""Plot PETSc KSP convergence histories from this example's log.txt.

The script is intentionally dependency-free: it writes a standalone SVG rather
than requiring matplotlib in the active Python environment.
"""

from __future__ import annotations

import html
import math
import re
import sys
from pathlib import Path, PurePath


LOG_PATH = Path(__file__).with_name("log.txt")
OUTPUT_PATH = Path(__file__).with_name("convergence_rates.svg")


def parse_log(path: Path) -> list[dict[str, object]]:
    """Extract completed mesh solves and their residual histories."""
    cases: list[dict[str, object]] = []
    current: dict[str, object] | None = None

    for line in path.read_text(encoding="utf-8").splitlines():
        mesh_match = re.match(r"Mesh:\s*(\S+)", line)
        if mesh_match:
            current = {
                "mesh": mesh_match.group(1),
                "hmax": None,
                "iterations": [],
                "residuals": [],
                "converged_iterations": None,
            }
            cases.append(current)
            continue

        if current is None:
            continue

        hmax_match = re.match(r"hmax:\s*([0-9.eE+-]+)", line)
        if hmax_match:
            current["hmax"] = float(hmax_match.group(1))
            continue

        residual_match = re.match(
            r"\s*(\d+)\s+KSP Residual norm\s+([0-9.eE+-]+)", line
        )
        if residual_match:
            iterations = current["iterations"]
            residuals = current["residuals"]
            assert isinstance(iterations, list)
            assert isinstance(residuals, list)
            iterations.append(int(residual_match.group(1)))
            residuals.append(float(residual_match.group(2)))
            continue

        converged_match = re.match(
            r"Linear solve converged due to \S+ iterations (\d+)", line
        )
        if converged_match:
            current["converged_iterations"] = int(converged_match.group(1))

    completed_cases = []
    for case in cases:
        iterations = case["iterations"]
        residuals = case["residuals"]
        if (
            isinstance(iterations, list)
            and isinstance(residuals, list)
            and iterations
            and residuals
            and case["converged_iterations"] is not None
        ):
            completed_cases.append(case)
    return completed_cases


def decade_label(exponent: int) -> str:
    if exponent == 0:
        return "1"
    return f"1e{exponent}"


def polyline_path(points: list[tuple[float, float]]) -> str:
    return " ".join(
        f"{'M' if index == 0 else 'L'} {x:.2f} {y:.2f}"
        for index, (x, y) in enumerate(points)
    )


def add_text(
    parts: list[str],
    x: float,
    y: float,
    text: str,
    *,
    size: int = 13,
    fill: str = "#172033",
    anchor: str = "start",
    weight: str = "400",
    rotate: float | None = None,
) -> None:
    transform = f' transform="rotate({rotate:.1f} {x:.1f} {y:.1f})"' if rotate else ""
    parts.append(
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-family="Arial, Helvetica, sans-serif" fill="{fill}" '
        f'text-anchor="{anchor}" font-weight="{weight}"{transform}>'
        f"{html.escape(text)}</text>"
    )


def build_svg(cases: list[dict[str, object]]) -> str:
    if not cases:
        raise ValueError(f"No completed KSP residual histories found in {LOG_PATH}")

    width = 1120
    height = 700
    left = 92
    top = 82
    plot_width = 710
    plot_height = 500
    bottom = top + plot_height
    right = left + plot_width
    legend_x = 840
    legend_y = 126

    max_iteration = max(
        int(case["iterations"][-1])  # type: ignore[index]
        for case in cases
    )
    x_step = 50
    x_max = int(math.ceil(max_iteration / x_step) * x_step)

    relative_histories = []
    for case in cases:
        residuals = case["residuals"]
        assert isinstance(residuals, list)
        initial_residual = float(residuals[0])
        relative_histories.append([float(value) / initial_residual for value in residuals])

    min_relative = min(min(history) for history in relative_histories)
    y_min_exp = min(-1, math.floor(math.log10(min_relative)))
    y_max_exp = 0
    if y_min_exp > -8:
        y_min_exp = -8

    def x_coord(iteration: float) -> float:
        return left + (iteration / x_max) * plot_width

    def y_coord(relative_residual: float) -> float:
        log_value = math.log10(max(relative_residual, 1.0e-300))
        fraction = (log_value - y_min_exp) / (y_max_exp - y_min_exp)
        return bottom - fraction * plot_height

    colors = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#9467bd",
        "#ff7f0e",
        "#17becf",
        "#8c564b",
    ]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
    ]

    add_text(
        parts,
        left,
        38,
        "MINRES convergence from log.txt",
        size=24,
        weight="700",
    )
    add_text(
        parts,
        left,
        62,
        "Relative preconditioned residual norm; rho = (final / initial)^(1 / N)",
        size=13,
        fill="#5b6472",
    )

    parts.append(
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" '
        'fill="#fbfcfe" stroke="#c8ced8" stroke-width="1"/>'
    )

    for tick in range(0, x_max + 1, x_step):
        x = x_coord(tick)
        parts.append(
            f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{bottom}" '
            'stroke="#e6e9ef" stroke-width="1"/>'
        )
        parts.append(
            f'<line x1="{x:.1f}" y1="{bottom}" x2="{x:.1f}" y2="{bottom + 6}" '
            'stroke="#6b7280" stroke-width="1"/>'
        )
        add_text(parts, x, bottom + 24, str(tick), size=12, fill="#4b5563", anchor="middle")

    for exponent in range(y_max_exp, y_min_exp - 1, -1):
        y = y_coord(10.0**exponent)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" '
            'stroke="#e6e9ef" stroke-width="1"/>'
        )
        parts.append(
            f'<line x1="{left - 6}" y1="{y:.1f}" x2="{left}" y2="{y:.1f}" '
            'stroke="#6b7280" stroke-width="1"/>'
        )
        add_text(
            parts,
            left - 12,
            y + 4,
            decade_label(exponent),
            size=12,
            fill="#4b5563",
            anchor="end",
        )

    parts.append(
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" '
        'stroke="#111827" stroke-width="1.2"/>'
    )
    parts.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" '
        'stroke="#111827" stroke-width="1.2"/>'
    )
    add_text(parts, (left + right) / 2, height - 48, "MINRES iteration", size=14, anchor="middle")
    add_text(
        parts,
        24,
        (top + bottom) / 2,
        "relative residual norm",
        size=14,
        anchor="middle",
        rotate=-90,
    )

    legend_lines = []
    for index, case in enumerate(cases):
        iterations = case["iterations"]
        residuals = case["residuals"]
        assert isinstance(iterations, list)
        assert isinstance(residuals, list)
        initial_residual = float(residuals[0])
        relative = [float(value) / initial_residual for value in residuals]
        last_iteration = int(iterations[-1])
        convergence_factor = relative[-1] ** (1.0 / last_iteration)
        orders_per_iteration = -math.log10(relative[-1]) / last_iteration
        color = colors[index % len(colors)]
        points = [
            (x_coord(float(iteration)), y_coord(value))
            for iteration, value in zip(iterations, relative)
        ]
        parts.append(
            f'<path d="{polyline_path(points)}" fill="none" stroke="{color}" '
            'stroke-width="2.2" stroke-linejoin="round" stroke-linecap="round"/>'
        )
        parts.append(
            f'<circle cx="{points[-1][0]:.1f}" cy="{points[-1][1]:.1f}" r="3.8" '
            f'fill="{color}" stroke="#ffffff" stroke-width="1.2"/>'
        )

        mesh_name = PurePath(str(case["mesh"])).name
        hmax = case["hmax"]
        hmax_text = f"{float(hmax):.4g}" if hmax is not None else "n/a"
        legend_lines.append(
            (
                color,
                mesh_name,
                f"h={hmax_text}, N={last_iteration}, rho={convergence_factor:.5f}",
                f"final rel={relative[-1]:.2e}, dec/it={orders_per_iteration:.5f}",
            )
        )

    parts.append(
        f'<rect x="{legend_x - 18}" y="{legend_y - 42}" width="252" height="{82 + 66 * len(legend_lines)}" '
        'rx="6" fill="#ffffff" stroke="#d5dae3" stroke-width="1"/>'
    )
    add_text(parts, legend_x, legend_y - 15, "Completed mesh solves", size=14, weight="700")
    for row, (color, mesh_name, line_one, line_two) in enumerate(legend_lines):
        y = legend_y + 20 + row * 66
        parts.append(
            f'<line x1="{legend_x}" y1="{y:.1f}" x2="{legend_x + 30}" y2="{y:.1f}" '
            f'stroke="{color}" stroke-width="3" stroke-linecap="round"/>'
        )
        add_text(parts, legend_x + 40, y + 4, mesh_name, size=13, weight="700")
        add_text(parts, legend_x + 40, y + 23, line_one, size=12, fill="#4b5563")
        add_text(parts, legend_x + 40, y + 41, line_two, size=12, fill="#4b5563")

    parts.append("</svg>")
    return "\n".join(parts)


def main() -> int:
    log_path = Path(sys.argv[1]) if len(sys.argv) > 1 else LOG_PATH
    output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else OUTPUT_PATH
    cases = parse_log(log_path)
    svg = build_svg(cases)
    output_path.write_text(svg, encoding="utf-8")

    print(f"Wrote {output_path}")
    for case in cases:
        iterations = case["iterations"]
        residuals = case["residuals"]
        assert isinstance(iterations, list)
        assert isinstance(residuals, list)
        relative_final = float(residuals[-1]) / float(residuals[0])
        rho = relative_final ** (1.0 / int(iterations[-1]))
        mesh_name = PurePath(str(case["mesh"])).name
        print(
            f"{mesh_name}: h={float(case['hmax']):.6g}, "
            f"N={int(iterations[-1])}, final_rel={relative_final:.3e}, rho={rho:.5f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
