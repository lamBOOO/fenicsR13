"""
Check the exact solutions themselves, not only the solver.

The convergence tests measure the distance between the discrete solution and
an exact solution. They cannot tell whether that exact solution is the right
one: a wrong reference would still produce clean convergence, just towards the
wrong limit. These tests insert the exported exact solution into the balance
equations and require the residual to follow the working precision, which an
only approximately satisfied equation cannot do.

The tests need ``mpmath`` but not dolfin, so they also run outside the
container.
"""

import math
import os
import re
from fractions import Fraction

import pytest

from esols.verify_exact_solution import (
    EQUATIONS,
    R14AnnulusSolution,
    SAMPLE_POINTS,
    check_precision_scaling,
    read_constants,
    read_layer_rates,
    residuals,
    stress_balance_residual,
)

HERE = os.path.dirname(os.path.abspath(__file__))
R14_ESOL = os.path.join(
    HERE, "esols", "1_coeffs_nosources_norot_inflow_positive.cpp"
)

# Residuals come from numerical differentiation, so they sit a few digits above
# the working precision. The margin is what separates "exact up to rounding"
# from "not a solution".
PRECISION_MARGIN = 1e5

# One interior point, used by the checks that do not need the whole set.
SAMPLE_POINT = (0.7, 0.3)


@pytest.fixture(name="solution", scope="module")
def fixture_solution():
    """Exact R14 solution built from the constants in the exported C++ file."""
    return R14AnnulusSolution(read_constants(R14_ESOL))


@pytest.fixture(name="precision_table", scope="module")
def fixture_precision_table(solution):
    """
    Residuals of every equation at two working precisions.

    Computed once for the whole module. The five cheap equations are evaluated
    on the full set of sample points. The stress balance needs second
    derivatives of the stress, which costs about fifty times as much, so it is
    evaluated on one interior point; the precision scaling is what the test
    looks at, and one point already shows it (4.7e-17 at 15 digits against
    1.7e-31 at 30). The full sweep over all points and three precisions is a
    manual run, see the module docstring of ``verify_exact_solution``.
    """
    cheap = [name for name in EQUATIONS if name != "stress"]
    table = check_precision_scaling(
        solution, precisions=(15, 30), equations=cheap
    )
    expensive = check_precision_scaling(
        solution,
        precisions=(15, 30),
        points=SAMPLE_POINTS[:1],
        equations=("stress",),
    )
    for precision, values in expensive.items():
        table[precision].update(values)
    return table


def test_constants_are_complete(solution):
    """The exported file must define all 28 integration constants."""
    assert len(solution.c) == 28


def test_knudsen_layer_rate_matches_export():
    """
    The decay rate of the new moment follows from its own balance equation.

    The R equation reduces to a modified Helmholtz equation with rate
    ``sqrt(2/7)``. The rate is read from the exported file rather than repeated
    here, so that the comparison still means something if the export changes.
    """
    rates = read_layer_rates(R14_ESOL)
    derived = math.sqrt(2 / 7)
    closest = min(float(value) for value in rates.values())
    assert closest == pytest.approx(derived, rel=1e-15)


def test_all_layer_rates_are_square_roots_of_simple_fractions():
    """
    Every layer rate in the export is the square root of a simple fraction.

    The four rates are ``sqrt(2/7)``, ``sqrt(5/6)``, ``sqrt(3/2)`` and
    ``sqrt(5/9)``; the first is the one R14 adds, the other three already
    appear in the R13 solutions. A rate that is not of this form would mean the
    export carries a rounded or wrong constant.
    """
    rates = read_layer_rates(R14_ESOL)
    assert len(rates) == 4
    for name, value in rates.items():
        fraction = Fraction(float(value) ** 2).limit_denominator(60)
        assert fraction.denominator <= 9, (
            "%s = %s is not a square root of a simple fraction" % (name, value)
        )
        assert float(value) == pytest.approx(
            math.sqrt(fraction.numerator / fraction.denominator), rel=1e-15
        )


def test_rank_three_term_is_needed(solution):
    """
    Third negative control, on the term that is hardest to get right.

    The stress balance closes with the rank-three moment
    ``m = -2 Kn (grad sigma)_stf``. Dropping it has to break the balance;
    otherwise the stress residual would be passing for the wrong reason, and
    an error in the rank-three algebra would go unnoticed.
    """
    radius, angle = SAMPLE_POINT
    with_term, scale = stress_balance_residual(solution, radius, angle)
    without_term, _ = stress_balance_residual(
        solution, radius, angle, with_rank_three=False
    )

    assert float(with_term / scale) < 1e-12
    assert float(without_term / scale) > 1e-3
    assert without_term > with_term * 1e6


@pytest.mark.parametrize(
    "equation",
    ["mass", "energy", "scalar_moment", "momentum", "heat_flux", "stress"],
)
def test_residual_follows_working_precision(precision_table, equation):
    """
    The residual must shrink when the working precision is raised.

    A solution that satisfies the equation exactly leaves only rounding error,
    which tracks the precision. A solution that satisfies it approximately
    leaves a residual that stops shrinking.
    """
    coarse = float(precision_table[15][equation])
    fine = float(precision_table[30][equation])

    if coarse == 0.0:
        # Exactly representable cases (mass, energy here) stay at zero.
        assert fine == 0.0
        return

    assert fine < coarse / 1e10, (
        "residual of %s does not follow the working precision: "
        "%.1e at 15 digits, %.1e at 30 digits" % (equation, coarse, fine)
    )


def test_check_detects_a_wrong_stress(solution):
    """
    Guard against a check that cannot fail.

    A deliberately perturbed stress field must break the momentum balance. If
    this test ever passes silently, the residual computation is not looking at
    the fields it claims to look at.
    """

    class PerturbedSolution(R14AnnulusSolution):
        """Same solution with a one percent error in the radial stress."""

        def stress(self, radius, angle):
            s_rr, s_rp, s_pp = super().stress(radius, angle)
            return 1.01 * s_rr, s_rp, s_pp

    perturbed = PerturbedSolution(solution.constants)
    only = ("momentum",)
    clean = float(residuals(solution, equations=only)["momentum"])
    broken = float(residuals(perturbed, equations=only)["momentum"])

    assert broken > 1e-6
    assert broken > clean * 1e10


def test_check_detects_a_wrong_scalar_moment(solution):
    """
    Second negative control, on the field that is new in R14.

    The scalar moment enters the heat flux balance through ``1/6 grad(R)``.
    Perturbing it must therefore break that balance, which also proves the
    heat flux residual really uses the R14 form and not the R13 one.
    """

    class PerturbedSolution(R14AnnulusSolution):
        """Same solution with a one percent error in the scalar moment."""

        def scalar_moment(self, radius, angle):
            return 1.01 * super().scalar_moment(radius, angle)

    perturbed = PerturbedSolution(solution.constants)
    only = ("heat_flux",)
    clean = float(residuals(solution, equations=only)["heat_flux"])
    broken = float(residuals(perturbed, equations=only)["heat_flux"])

    assert broken > 1e-6
    assert broken > clean * 1e10


# Of the four exported exact solutions in ``esols``, only one is an R14
# solution. The other three are byte-identical copies of the R13 files: they
# carry three Knudsen layers instead of four, no ``sqrt(2/7)`` layer and no
# ``ScriptR`` class. The case below points at one of those copies, which is
# why it aborts during error computation.
INCONSISTENT_CASE = "r14_1_coeffs_sources_rot_noinflow_p1p1p1p1p1_gls"

R14_INPUTS = [
    "r14_1_coeffs_nosources_norot_inflow_p1p1p1p1p1_stab",
    "r14_1_coeffs_nosources_norot_inflow_p2p2p2p2p2_stab",
    pytest.param(
        INCONSISTENT_CASE,
        marks=pytest.mark.xfail(
            reason="exact solution is the R13 one; the R14 solution for this "
                   "parameter set is not in the repository",
            strict=True,
        ),
    ),
]


@pytest.mark.parametrize("case", R14_INPUTS)
def test_r14_input_points_at_an_r14_exact_solution(case):
    """
    An ``r14`` case must reference an exact solution that carries R14.

    The solver only notices the mismatch deep inside the error computation,
    where it fails on a missing ``ScriptR`` class. Checking the pairing here
    names the real cause instead: a reference solution for the other theory.
    """
    with open(os.path.join(HERE, "inputs", case + ".yml")) as handle:
        text = handle.read()
    assert "mode: r14" in text

    # Anchored to the start of the line so that the commented-out
    # documentation line "# - exact_solution: Path to ..." is not picked up.
    match = re.search(
        r"^[ \t]*exact_solution:[ \t]*(\S+)", text, re.MULTILINE
    )
    assert match, "no exact solution configured in %s" % case
    esol = os.path.join(HERE, match.group(1))

    rates = read_layer_rates(esol)
    with open(esol) as handle:
        source = handle.read()

    assert "class ScriptR" in source, (
        "%s has no ScriptR class, so it cannot carry the R14 moment" % esol
    )
    assert any(
        float(rate) == pytest.approx(math.sqrt(2 / 7), rel=1e-15)
        for rate in rates.values()
    ), "%s has no sqrt(2/7) Knudsen layer, so it is not an R14 solution" % esol
