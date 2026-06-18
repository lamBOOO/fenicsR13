"""
Tests for bubble-enriched finite element construction.
"""

import subprocess

import pytest
import yaml

df = pytest.importorskip("dolfin")

from fenicsR13.input import Input


def test_input_accepts_bubble_enriched_sigma():
    """The paper-style R13 input should pass schema validation."""
    Input(
        "tests/2d_r13/inputs/"
        "r13_1_coeffs_nosources_norot_inflow_p1p2p1p2p2b_nostab.yml"
    )


def test_sigma_bubble_enrichment_builds_larger_space():
    """The sigma element should be enriched by cell bubbles."""
    from fenicsR13.solver import Solver  # pylint: disable=C0415

    solver = Solver.__new__(Solver)
    solver.cell = df.triangle
    solver.nsd = 2
    solver.var_ranks = {
        "theta": 0,
        "s": 1,
        "p": 0,
        "u": 1,
        "sigma": 2,
    }
    solver.params = {
        "elements": {
            "sigma": {
                "shape": "Lagrange",
                "degree": 2,
                "bubble_enriched": True,
            }
        }
    }

    symmetry = solver._Solver__sigma_symmetry()
    enriched = solver._Solver__create_element("sigma", symmetry=symmetry)
    base = df.TensorElement("Lagrange", df.triangle, 2, symmetry=symmetry)

    mesh = df.UnitSquareMesh(1, 1)
    assert solver._Solver__bubble_degree("sigma") == 4
    assert solver._Solver__element_output_degree("sigma") == 4
    assert df.FunctionSpace(mesh, enriched).dim() > df.FunctionSpace(
        mesh, base
    ).dim()


def test_paper_bubble_r13_smoke_solve(tmp_path):
    """The unstabilized paper-style element should solve on a coarse mesh."""
    input_path = (
        "tests/2d_r13/inputs/"
        "r13_1_coeffs_nosources_norot_inflow_p1p2p1p2p2b_nostab.yml"
    )
    with open(input_path, "r") as stream:
        params = yaml.safe_load(stream)

    output_folder = tmp_path / "paper_bubble_r13_smoke"
    params["meshes"] = ["../2d_mesh/ring0.h5"]
    params["output_folder"] = str(output_folder)

    smoke_input = tmp_path / "input.yml"
    with open(smoke_input, "w") as stream:
        yaml.safe_dump(params, stream)

    subprocess.check_call(
        ["fenicsR13", str(smoke_input)],
        cwd="tests/2d_r13"
    )
    assert (output_folder / "errors.csv").exists()
