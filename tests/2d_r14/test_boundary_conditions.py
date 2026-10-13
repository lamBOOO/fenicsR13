"""
Check the exact solutions against the wall boundary conditions.

``test_exact_solution.py`` inserts the exported solutions into the balance
equations. That is necessary but not sufficient: the balance equations hold
for any choice of the integration constants, because every homogeneous
solution satisfies them by construction. The constants are fixed by the
boundary conditions alone, and the conditions were not checked anywhere.

These tests close that gap. They insert the same solutions into the wall
conditions of ``R14analytical/R14BoundaryConditions.wls`` (seven equations)
and ``R13BoundaryConditions.wls`` (the same six without the scalar moment
terms), on both walls and at four angles. The negative control shows what the
balance equations cannot show: perturbing one constant by a relative 1e-5
moves these residuals by many orders of magnitude.

The tests need ``mpmath`` but not dolfin, so they also run outside the
container.
"""

import os

import pytest

from esols.verify_boundary_conditions import (
    CONDITIONS,
    R13_CONDITIONS,
    SAMPLE_ANGLES,
    boundary_residuals,
)
from esols.verify_exact_solution import (
    R13AnnulusSolution,
    R14AnnulusSolution,
    read_constants,
)

HERE = os.path.dirname(os.path.abspath(__file__))
R14_ESOL = os.path.join(
    HERE, "esols", "1_coeffs_nosources_norot_inflow_positive.cpp"
)
R13_ESOL = os.path.join(
    HERE, os.pardir, "2d_r13", "esols",
    "1_coeffs_nosources_norot_inflow_positive.cpp",
)

# The exported constants are decimal numbers of double precision, so the
# conditions cannot cancel below about 1e-16 no matter how many digits the
# evaluation uses. Cancellation between terms of different size lifts that
# floor by a few orders; 1e-10 separates "satisfied" from "not satisfied"
# while leaving room for it.
SATISFIED = 1e-10

# A constant that is off by a relative 1e-5 must show up far above the floor.
DETECTION_FACTOR = 1e4


@pytest.fixture(scope="module")
def r14_solution():
    """The exported R14 solution of the inflow case."""
    return R14AnnulusSolution(read_constants(R14_ESOL))


@pytest.fixture(scope="module")
def r13_solution():
    """The exported R13 solution of the same geometry."""
    return R13AnnulusSolution(read_constants(R13_ESOL))


@pytest.fixture(scope="module")
def r14_residuals(r14_solution):
    """Worst relative residual of each R14 condition, over both walls."""
    return boundary_residuals(r14_solution)


@pytest.fixture(scope="module")
def r13_residuals(r13_solution):
    """Worst relative residual of each R13 condition, over both walls."""
    return boundary_residuals(r13_solution, conditions=R13_CONDITIONS)


@pytest.mark.parametrize("condition", CONDITIONS)
def test_r14_solution_satisfies_the_wall_conditions(r14_residuals, condition):
    """Every one of the seven R14 conditions holds on both walls."""
    assert float(r14_residuals[condition]) < SATISFIED, (
        "R14 condition %s: relative residual %.3e"
        % (condition, float(r14_residuals[condition]))
    )


@pytest.mark.parametrize("condition", R13_CONDITIONS)
def test_r13_solution_satisfies_the_wall_conditions(r13_residuals, condition):
    """
    The same for R13, where the scalar moment is closed instead of solved.

    ``Delta = -12 Kn div(s)`` vanishes identically for this case, which is why
    the R14 expressions reduce to the R13 ones when the scalar moment is
    carried along.
    """
    assert float(r13_residuals[condition]) < SATISFIED, (
        "R13 condition %s: relative residual %.3e"
        % (condition, float(r13_residuals[condition]))
    )


def test_check_detects_a_wrong_integration_constant(r14_solution):
    """
    A constant that is off must break the conditions.

    This is the test the balance equations cannot provide: there, the same
    perturbation changes the residual by a factor below two, because a
    different constant is still a solution of the homogeneous system. Here it
    has to be visible, otherwise the conditions would not be what fixes the
    constants.

    One angle is enough to see it and keeps the test short; the perturbation
    is not angle specific.
    """
    constants = r14_solution.constants

    def carries_a_layer(key):
        """``CK`` constants multiply the Knudsen layers, ``C`` the bulk."""
        return key.startswith("CK") and abs(float(constants[key])) > 1e-9

    name = sorted(key for key in constants if carries_a_layer(key))[0]
    disturbed = dict(r14_solution.constants)
    disturbed[name] = float(disturbed[name]) * (1 + 1e-5)

    angles = SAMPLE_ANGLES[:1]
    worst_clean = max(
        float(value) for value in
        boundary_residuals(r14_solution, angles=angles).values()
    )
    worst_disturbed = max(
        float(value) for value in
        boundary_residuals(R14AnnulusSolution(disturbed),
                           angles=angles).values()
    )
    assert worst_disturbed > worst_clean * DETECTION_FACTOR, (
        "perturbing %s moved the worst residual only from %.3e to %.3e"
        % (name, worst_clean, worst_disturbed)
    )
