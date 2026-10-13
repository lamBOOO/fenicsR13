"""
Check the Knudsen layer rates of the exported exact solutions.

The rates sit in the exported files as plain constants, and nothing checks
where they come from. They are the non-zero roots of the characteristic
polynomial of the homogeneous moment system, so they can be derived from the
equations and compared against the export. That turns four hard-coded numbers
into a checkable consequence of the theory.

These tests need ``sympy`` but not dolfin, so they also run outside the
container.
"""

import os

import pytest

sympy = pytest.importorskip("sympy")

from esols.layer_rates import (  # noqa: E402
    characteristic_polynomial,
    layer_rates,
    read_layer_rates,
    read_weak_form_coefficients,
    scalar_moment_coefficient_from_weak_form,
)

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.path.join(HERE, os.pardir, os.pardir, "fenicsR13", "solver.py")
R14_ESOL = os.path.join(
    HERE, "esols", "1_coeffs_nosources_norot_inflow_positive.cpp"
)
R13_ESOL = os.path.join(
    HERE, os.pardir, "2d_r13", "esols",
    "1_coeffs_nosources_norot_inflow_positive.cpp",
)


def test_r13_rates_follow_from_the_equations():
    """The three R13 rates are roots of the R13 system, not free constants."""
    derived, _ = layer_rates("r13")
    exported = set(read_layer_rates(R13_ESOL).values())

    assert set(derived) == exported, (
        "derived %s, exported %s"
        % (sorted(derived, key=float), sorted(exported, key=float))
    )


def test_r14_rates_follow_from_the_equations():
    """The R14 export adds one rate, and the system produces exactly it."""
    derived, _ = layer_rates("r14")
    exported = set(read_layer_rates(R14_ESOL).values())

    assert set(derived) == exported, (
        "derived %s, exported %s"
        % (sorted(derived, key=float), sorted(exported, key=float))
    )
    assert sympy.Rational(2, 7) in derived, (
        "the rate that is new in R14 is missing"
    )


def test_r14_adds_exactly_one_layer():
    """
    R14 leaves the R13 layers alone.

    The R14 characteristic polynomial is the R13 one times a single new
    factor. That is a stronger statement than the rates agreeing one by one:
    it says the new moment does not shift the existing layers.
    """
    r13 = characteristic_polynomial("r13")
    r14 = characteristic_polynomial("r14")

    quotient, remainder = sympy.div(r14.as_expr(), r13.as_expr(),
                                    sympy.symbols("y"))
    assert sympy.simplify(remainder) == 0, (
        "the R14 polynomial is not divisible by the R13 one"
    )
    quotient = sympy.factor(sympy.simplify(quotient))
    roots = sympy.roots(sympy.Poly(quotient, sympy.symbols("y")))
    non_zero = [root for root in roots if root != 0]
    assert non_zero == [sympy.Rational(2, 7)], (
        "the extra factor should contribute only 2/7, got %s" % non_zero
    )


def test_the_truncated_omega_closure_keeps_the_layers():
    r"""
    The shortened closure for Omega costs nothing in the layer rates.

    The R14 derivation closes the fifth-order flux as
    :math:`\Omega = -\frac{7}{3} Kn \nabla \mathcal{R}`, equation (3.22). The
    linear R26 system of Rana, Gupta, Sprittles and Torrilhon, JFM 924 (2021)
    A16, equation (21), carries one term more,
    :math:`-\frac{7}{3} Kn \frac{12}{7} \partial_j R_{ij}`. In the R14 system
    :math:`R_{ij}` is closed through the heat flux, so that term is a higher
    gradient of ``s`` and could in principle move the layers.

    It does not: both closures give the same four rates. The truncation is
    therefore invisible to the Knudsen layer structure, which is worth
    recording because the two systems differ elsewhere.
    """
    thesis, _ = layer_rates("r14", omega_closure="thesis")
    published, _ = layer_rates("r14", omega_closure="lr26")

    assert set(thesis) == set(published), (
        "thesis closure %s, R26 closure %s"
        % (sorted(thesis, key=float), sorted(published, key=float))
    )


def test_the_bulk_root_is_not_counted_as_a_layer():
    """
    The root at zero belongs to the bulk solution.

    It must be there, and it must not be reported as a layer rate: a solution
    that decays nowhere is the polynomial part of the solution, not a Knudsen
    layer.
    """
    for mode in ("r13", "r14"):
        rates, bulk = layer_rates(mode)
        assert bulk > 0, "no bulk root in %s" % mode
        assert all(rate != 0 for rate in rates)


@pytest.mark.parametrize("knudsen", [1, 3])
def test_rates_do_not_depend_on_the_knudsen_number(knudsen):
    """
    A Knudsen layer has to thin out as the Knudsen number falls.

    Its thickness is ``Kn / lambda``, so the rate itself must be independent
    of ``Kn``. This is true for all three R13 rates and, with the
    dimensionally consistent reading of the scalar moment balance, for the
    fourth one as well.
    """
    reference, _ = layer_rates("r14", knudsen=sympy.Integer(1))
    scaled, _ = layer_rates("r14", knudsen=sympy.Integer(knudsen))
    assert set(reference) == set(scaled)


def test_literal_source_reading_breaks_the_knudsen_scaling():
    """
    Negative control, and a question about the Mathematica source.

    ``R14Equations.wls`` writes the last term of the scalar moment balance as
    ``(2/3*Kn)*scriptR``, which reads as ``2/3 Kn``; every other equation uses
    ``2/(3 Kn)``. At ``Kn = 1`` the two agree, and every test case uses
    ``Kn = 1``, so nothing in the suite can tell them apart.

    They are not the same, though: with the literal reading the new rate comes
    out as ``2 Kn^2 / 7`` instead of ``2/7``, so the layer thickness
    ``Kn / lambda`` becomes ``sqrt(7/2)`` regardless of ``Kn`` and the layer
    stops thinning. This test pins that difference down, so the choice is
    recorded rather than implicit.
    """
    consistent_one, _ = layer_rates("r14", knudsen=sympy.Integer(1))
    consistent_three, _ = layer_rates("r14", knudsen=sympy.Integer(3))
    literal_one, _ = layer_rates(
        "r14", knudsen=sympy.Integer(1), literal_source=True
    )
    literal_three, _ = layer_rates(
        "r14", knudsen=sympy.Integer(3), literal_source=True
    )

    # At Kn = 1 the two readings cannot be distinguished.
    assert set(consistent_one) == set(literal_one)

    # Away from Kn = 1 they differ, and only one of them keeps the rate fixed.
    assert set(consistent_three) == set(consistent_one)
    assert set(literal_three) != set(literal_one)
    assert sympy.Rational(2, 7) * 9 in literal_three, (
        "with the literal reading the new rate should scale with Kn^2"
    )


def test_the_solver_row_is_a_multiple_of_the_balance():
    """
    The three numbers in row A[5] are not independent.

    Row ``A[5] = j(s, Q) + i(R, Q)`` is the scalar moment balance tested with
    ``Q``, scaled by some factor. The coupling term and the gradient term each
    determine that factor on their own, so they have to agree. If they do not,
    the row is not a scaled balance and reading a coefficient out of it would
    be meaningless.
    """
    coefficients = read_weak_form_coefficients(SOLVER)

    alpha_from_coupling = sympy.Rational(
        1, coefficients["coupling_denominator"]
    ) / 8
    alpha_from_gradient = sympy.Rational(
        coefficients["gradient_numerator"],
        coefficients["gradient_denominator"],
    ) / sympy.Rational(7, 3)

    assert alpha_from_coupling == alpha_from_gradient, (
        "coupling gives %s, gradient gives %s"
        % (alpha_from_coupling, alpha_from_gradient)
    )


def test_the_solver_implements_the_dimensionally_consistent_reading():
    """
    Which of the two readings is in the code, measured instead of argued.

    With the scale factor fixed by the other two terms, the mass term leaves
    no freedom: it gives the coefficient of the scalar moment outright. The
    answer is ``2 / (3 Kn)``, the spelling every other equation uses, and not
    the literal ``2/3 Kn`` of ``R14Equations.wls``. The two differ by
    ``Kn^2``, so keeping ``Kn`` symbolic is what makes the statement possible.
    """
    knudsen = sympy.symbols("Kn", positive=True)
    recovered = scalar_moment_coefficient_from_weak_form(SOLVER, knudsen)

    consistent = sympy.Rational(2, 3) / knudsen
    literal = sympy.Rational(2, 3) * knudsen

    assert sympy.simplify(recovered - consistent) == 0, (
        "recovered %s, expected %s" % (recovered, consistent)
    )
    assert sympy.simplify(recovered - literal) != 0, (
        "the two readings must be distinguishable with Kn symbolic"
    )


def test_the_implemented_coefficient_reproduces_the_exported_layer():
    """
    Close the loop from the implemented weak form to the exported constant.

    The coefficient read out of the solver is fed back into the homogeneous
    system, and the rate it produces is compared against the number in the
    exported exact solution. This ties three independent places together: the
    weak form in the solver, the layer rate in the cpp export, and the balance
    equations the derivation starts from.
    """
    knudsen = sympy.Integer(3)
    recovered = scalar_moment_coefficient_from_weak_form(SOLVER, knudsen)

    rates, _ = layer_rates(
        "r14", knudsen=knudsen, scalar_moment_coefficient=recovered
    )
    exported = set(read_layer_rates(R14_ESOL).values())

    assert set(rates) == exported, (
        "derived %s, exported %s"
        % (sorted(rates, key=float), sorted(exported, key=float))
    )
