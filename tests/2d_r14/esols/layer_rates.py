r"""
Derive the Knudsen layer rates from the moment equations themselves.

The exported exact solutions carry their layer rates as constants: the R13
files write ``sqrt(5.0/9.0)``, ``sqrt(5.0/6.0)``, ``sqrt(3.0/2.0)``, and the
R14 file adds a fourth one, ``sqrt(2/7)``. Nothing checks where those numbers
come from, so a wrong rate in an export would simply become the reference.

They are not free parameters. They are the eigenvalues of the homogeneous
system. With the plane wave ansatz

.. math::
    \\text{every field} \\sim v \\, e^{\\lambda x / Kn}

the homogeneous equations become :math:`A(\\lambda) v = 0`, which has a
non-trivial solution only where :math:`\\det A(\\lambda) = 0`. The non-zero
roots are the layer rates; the root at zero carries the bulk solution, which
is polynomial in space rather than exponential.

Running this for R13 and for R14 reproduces exactly the constants in the two
exported files, and shows that the R14 characteristic polynomial is the R13
one multiplied by a single new factor: R14 adds one Knudsen layer and leaves
the other three untouched.

The same machinery applies to any linear moment system whose balance
equations and closures are known, which is what makes it useful beyond
checking the present export.

Equations and closures are taken from the Mathematica sources that ship with
the R14 exact solutions, ``R14analytical/R13Equations.wls`` and
``R14Equations.wls``:

.. math::
    \\nabla \\cdot u &= 0 \\\\
    \\nabla p + \\nabla \\cdot \\sigma &= 0 \\\\
    \\nabla \\cdot s + \\nabla \\cdot u &= 0 \\\\
    2 (\\nabla u)_{stf} + \\frac{4}{5} (\\nabla s)_{stf}
      + \\frac{1}{Kn} \\sigma + \\nabla \\cdot m &= 0 \\\\
    \\frac{5}{2} \\nabla \\theta + \\nabla \\cdot \\sigma
      + \\frac{1}{2} \\nabla \\cdot R + \\frac{2}{3 Kn} s
      \\;[+ \\tfrac{1}{6} \\nabla \\mathcal{R}] &= 0 \\\\
    8 \\nabla \\cdot s + \\nabla \\cdot \\Omega
      + \\frac{2}{3 Kn} \\mathcal{R} &= 0 \\quad (\\text{r14 only})

with the closures :math:`m = -2 Kn (\\nabla \\sigma)_{stf}`,
:math:`R = -\\frac{24}{5} Kn (\\nabla s)_{stf}` and
:math:`\\Omega = -\\frac{7}{3} Kn \\nabla \\mathcal{R}`.

Requires ``sympy``; it does not import dolfin.
"""

import re
from fractions import Fraction

import sympy as sp

# All six permutations of three indices, for the rank-three symmetrisation.
_PERMUTATIONS_OF_THREE = (
    (0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0),
)

# Largest denominator accepted when recognising a rate as a simple fraction.
MAX_DENOMINATOR = 60


def _delta(i, j):
    return sp.Integer(1) if i == j else sp.Integer(0)


def symmetric_trace_free_rank_two(tensor):
    """Return the symmetric trace-free part, with the trace over three axes."""
    symmetric = [
        [sp.expand((tensor[i][j] + tensor[j][i]) / 2) for j in range(3)]
        for i in range(3)
    ]
    trace = symmetric[0][0] + symmetric[1][1] + symmetric[2][2]
    return [
        [sp.expand(symmetric[i][j] - trace * _delta(i, j) / 3)
         for j in range(3)]
        for i in range(3)
    ]


def symmetric_trace_free_rank_three(tensor):
    """Return the symmetric trace-free part of a rank-three tensor."""
    symmetric = [[[sp.expand(sum(
        tensor[(i, j, k)[p[0]]][(i, j, k)[p[1]]][(i, j, k)[p[2]]]
        for p in _PERMUTATIONS_OF_THREE) / 6)
        for k in range(3)] for j in range(3)] for i in range(3)]
    trace = [sp.expand(sum(symmetric[i][j][j] for j in range(3)))
             for i in range(3)]
    return [[[sp.expand(
        symmetric[i][j][k]
        - (trace[i] * _delta(j, k) + trace[j] * _delta(i, k)
           + trace[k] * _delta(i, j)) / 5)
        for k in range(3)] for j in range(3)] for i in range(3)]


def read_layer_rates(cpp_file):
    """
    Return the squared layer rates stored in an exported exact solution.

    The R13 files spell them as ``sqrt(a/b)`` and the R14 files as decimals,
    so both are accepted. The result maps the constant name to the exact
    fraction :math:`\\lambda^2`.
    """
    with open(cpp_file) as handle:
        text = handle.read()

    rates = {}
    for name, numerator, denominator in re.findall(
        r"double\s+(lambda_\d+)\s*=\s*sqrt\(\s*([0-9.]+)\s*/\s*([0-9.]+)\s*\)",
        text,
    ):
        fraction = Fraction(
            float(numerator) / float(denominator)
        ).limit_denominator(MAX_DENOMINATOR)
        rates[name] = sp.Rational(fraction.numerator, fraction.denominator)
    for name, value in re.findall(
        r"double\s+(lambda_\d+)\s*=\s*([-0-9.eE+]+)\s*;", text
    ):
        if name in rates:
            continue
        fraction = Fraction(
            float(value) ** 2
        ).limit_denominator(MAX_DENOMINATOR)
        rates[name] = sp.Rational(fraction.numerator, fraction.denominator)
    return rates


def _scalar_moment_coefficient(knudsen, literal_source, override=None):
    """
    Return the coefficient of the scalar moment in its own balance.

    ``override`` takes precedence over both spellings and is what allows the
    coefficient recovered from the solver to be fed back in, so the loop from
    the implemented weak form to the layer rate can be closed.

    ``R14Equations.wls`` writes the term as ``(2/3*Kn)*scriptR``, which reads
    as :math:`\\frac{2}{3} Kn`. Every other equation uses
    :math:`\\frac{2}{3 Kn}`, and only that form keeps the layer rate
    independent of the Knudsen number, as the other three rates are. Both
    spellings agree at ``Kn = 1``, which is the value every test case uses, so
    the difference is invisible there. ``literal_source=True`` selects the
    spelling as written, which is what makes the two distinguishable.
    """
    if override is not None:
        return override
    if literal_source:
        return sp.Rational(2, 3) * knudsen
    return sp.Rational(2, 3) / knudsen


def _functional_body(text, name, signature):
    """Return the source of one weak-form functional of the solver."""
    start = text.index("def %s(%s):" % (name, signature))
    rest = text[start:]
    end = rest.index("\n        def ", 1)
    return rest[:end]


def read_weak_form_coefficients(solver_file):
    r"""
    Return the coefficients the solver uses in the scalar moment row.

    Row ``A[5]`` of the system is ``j(s, Q) + i(R, Q)``, so its volume terms
    carry three numbers: the coupling to the heat flux in ``j``, and the mass
    and gradient terms of :math:`\\mathcal{R}` in ``i``. They are read from the
    source rather than copied here, so a change in the solver shows up as a
    failure instead of being silently outvoted by a constant in a test.
    """
    with open(solver_file) as handle:
        text = handle.read()

    coupling_body = _functional_body(text, "j", "s, Q")
    mass_body = _functional_body(text, "i", "R, Q")

    coupling = re.search(
        r"-\s*\(1/(\d+)\)\s*\*\s*df\.inner\(s,\s*df\.grad\(Q\)\)",
        coupling_body,
    )
    mass = re.search(
        r"\(1\s*/\s*\((\d+)\s*\*\s*regs\[reg\]\[\"kn\"\]\)\)\s*\*\s*R\s*\*\s*Q",
        mass_body,
    )
    gradient = re.search(
        r"\(\((\d+)\s*\*\s*regs\[reg\]\[\"kn\"\]\)\s*/\s*(\d+)\)",
        mass_body,
    )
    if not (coupling and mass and gradient):
        raise ValueError(
            "could not read the scalar moment row from %s" % solver_file
        )
    return {
        "coupling_denominator": int(coupling.group(1)),
        "mass_denominator": int(mass.group(1)),
        "gradient_numerator": int(gradient.group(1)),
        "gradient_denominator": int(gradient.group(2)),
    }


def scalar_moment_coefficient_from_weak_form(solver_file, knudsen):
    r"""
    Recover the coefficient of :math:`\\mathcal{R}` that the solver implements.

    The strong balance is

    .. math::
        8 \\nabla \\cdot s + \\nabla \\cdot \\Omega + C \\mathcal{R} = 0,
        \\qquad \\Omega = -\\tfrac{7}{3} Kn \\nabla \\mathcal{R}

    Testing it with :math:`Q`, integrating the divergence of :math:`\\Omega` by
    parts and scaling the whole row by :math:`\\alpha` leaves three volume
    terms:

    .. math::
        -8 \\alpha \\int s \\cdot \\nabla Q
        + \\tfrac{7}{3} \\alpha Kn \\int \\nabla \\mathcal{R} \\cdot \\nabla Q
        + \\alpha C \\int \\mathcal{R} Q

    The solver supplies all three numbers, so the system is overdetermined:
    the coupling term and the gradient term each fix :math:`\\alpha` on their
    own, and only then does the mass term give :math:`C`. Their agreement is
    the evidence; it is checked here rather than assumed, and the recovered
    :math:`C` is returned with ``Kn`` kept symbolic so the two readings of the
    Mathematica source can be told apart.
    """
    coefficients = read_weak_form_coefficients(solver_file)

    # The strong form carries 8 * div(s), so the coupling term is -8 * alpha.
    alpha_from_coupling = sp.Rational(
        1, coefficients["coupling_denominator"]
    ) / 8
    # The gradient term is (7/3) * alpha * Kn; the Kn in the source cancels.
    alpha_from_gradient = sp.Rational(
        coefficients["gradient_numerator"],
        coefficients["gradient_denominator"],
    ) / sp.Rational(7, 3)

    if sp.simplify(alpha_from_coupling - alpha_from_gradient) != 0:
        raise ValueError(
            "the solver row is not a scalar multiple of the balance: the "
            "coupling term gives alpha = %s, the gradient term %s"
            % (alpha_from_coupling, alpha_from_gradient)
        )

    mass_coefficient = sp.Rational(
        1, coefficients["mass_denominator"]
    ) / knudsen
    return sp.simplify(mass_coefficient / alpha_from_coupling)


def build_homogeneous_system(mode, knudsen, rate, literal_source=False,
                             scalar_moment_coefficient=None,
                             omega_closure="thesis"):
    r"""
    Return the matrix of the homogeneous system under the plane wave ansatz.

    One spatial direction is enough: derivatives along the other two vanish,
    and the layer rates do not depend on the direction chosen.

    ``omega_closure`` selects how the fifth-order flux is closed, which only
    matters in r14 mode:

    ``"thesis"``
        :math:`\Omega = -\tfrac{7}{3} Kn \nabla \mathcal{R}`, equation (3.22)
        of the R14 derivation and what the solver implements.
    ``"lr26"``
        the closure of the linear R26 system, which carries one term more,
        :math:`\Omega_i = -\tfrac{7}{3} Kn (\partial_i \Delta
        + \tfrac{12}{7} \partial_j R_{ij})`.

    The second form is here to measure what the truncation costs, not because
    it belongs to the R14 system.
    """
    theta, pressure = sp.symbols("theta p")
    velocity = list(sp.symbols("u1 u2 u3"))
    heat_flux = list(sp.symbols("s1 s2 s3"))
    s11, s12, s13, s22, s23 = sp.symbols("S11 S12 S13 S22 S23")
    # The stress is symmetric and trace free, so five components are free.
    stress = [[s11, s12, s13],
              [s12, s22, s23],
              [s13, s23, -(s11 + s22)]]

    unknowns = ([theta, pressure] + velocity + heat_flux
                + [s11, s12, s13, s22, s23])
    scalar_moment = None
    if mode == "r14":
        scalar_moment = sp.symbols("scriptR")
        unknowns = unknowns + [scalar_moment]

    wave_number = rate / knudsen

    def derivative(index, expression):
        return wave_number * expression if index == 0 else sp.Integer(0)

    equations = []

    # Mass.
    equations.append(sp.expand(
        sum(derivative(i, velocity[i]) for i in range(3))
    ))

    # Momentum, three components.
    for i in range(3):
        equations.append(sp.expand(
            derivative(i, pressure)
            + sum(derivative(j, stress[i][j]) for j in range(3))
        ))

    # Energy.
    equations.append(sp.expand(
        sum(derivative(i, heat_flux[i]) for i in range(3))
        + sum(derivative(i, velocity[i]) for i in range(3))
    ))

    # Stress balance, five independent components.
    gradient_u = [[derivative(j, velocity[i]) for j in range(3)]
                  for i in range(3)]
    gradient_s = [[derivative(j, heat_flux[i]) for j in range(3)]
                  for i in range(3)]
    stf_u = symmetric_trace_free_rank_two(gradient_u)
    stf_s = symmetric_trace_free_rank_two(gradient_s)

    gradient_stress = [[[derivative(k, stress[i][j]) for k in range(3)]
                        for j in range(3)] for i in range(3)]
    moment = symmetric_trace_free_rank_three(gradient_stress)
    moment = [[[sp.expand(-2 * knudsen * moment[i][j][k]) for k in range(3)]
               for j in range(3)] for i in range(3)]
    divergence_moment = [
        [sp.expand(sum(derivative(k, moment[i][j][k]) for k in range(3)))
         for j in range(3)] for i in range(3)
    ]
    for (i, j) in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2)]:
        equations.append(sp.expand(
            2 * stf_u[i][j] + sp.Rational(4, 5) * stf_s[i][j]
            + stress[i][j] / knudsen + divergence_moment[i][j]
        ))

    # Heat flux balance, three components.
    fourth_moment = symmetric_trace_free_rank_two(gradient_s)
    fourth_moment = [
        [sp.expand(-sp.Rational(24, 5) * knudsen * fourth_moment[i][j])
         for j in range(3)] for i in range(3)
    ]
    divergence_fourth = [
        sp.expand(sum(derivative(j, fourth_moment[i][j]) for j in range(3)))
        for i in range(3)
    ]
    divergence_stress = [
        sp.expand(sum(derivative(j, stress[i][j]) for j in range(3)))
        for i in range(3)
    ]
    for i in range(3):
        expression = (
            sp.Rational(5, 2) * derivative(i, theta)
            + divergence_stress[i]
            + divergence_fourth[i] / 2
            + sp.Rational(2, 3) / knudsen * heat_flux[i]
        )
        if mode == "r14":
            expression += derivative(i, scalar_moment) / 6
        equations.append(sp.expand(expression))

    # Scalar moment balance, r14 only.
    if mode == "r14":
        omega_gradient = [derivative(i, scalar_moment) for i in range(3)]
        if omega_closure == "lr26":
            omega_gradient = [
                omega_gradient[i]
                + sp.Rational(12, 7) * divergence_fourth[i]
                for i in range(3)
            ]
        omega = [
            sp.expand(-sp.Rational(7, 3) * knudsen * omega_gradient[i])
            for i in range(3)
        ]
        equations.append(sp.expand(
            8 * sum(derivative(i, heat_flux[i]) for i in range(3))
            + sum(derivative(i, omega[i]) for i in range(3))
            + _scalar_moment_coefficient(
                knudsen, literal_source,
                override=scalar_moment_coefficient,
            ) * scalar_moment
        ))

    if len(equations) != len(unknowns):
        raise ValueError(
            "the system is not square: %d equations, %d unknowns"
            % (len(equations), len(unknowns))
        )

    matrix = sp.Matrix([
        [sp.expand(sp.diff(equation, unknown)) for unknown in unknowns]
        for equation in equations
    ])
    return matrix, unknowns


def characteristic_polynomial(mode, knudsen=sp.Integer(1),
                              literal_source=False,
                              scalar_moment_coefficient=None,
                              omega_closure="thesis"):
    r"""
    Return the characteristic polynomial in :math:`y = \\lambda^2`.

    The determinant contains only even powers of the rate, so rewriting it in
    :math:`y` halves the degree and keeps the arithmetic exact.
    """
    rate = sp.symbols("lambda_rate")
    matrix, _ = build_homogeneous_system(
        mode, knudsen, rate, literal_source=literal_source,
        scalar_moment_coefficient=scalar_moment_coefficient,
        omega_closure=omega_closure,
    )
    determinant = sp.expand(matrix.det(method="berkowitz"))

    polynomial = sp.Poly(determinant, rate)
    squared = sp.symbols("y")
    terms = {}
    for monomial, coefficient in zip(polynomial.monoms(),
                                     polynomial.coeffs()):
        power = monomial[0]
        if power % 2:
            raise ValueError(
                "the determinant has an odd power of the rate: %d" % power
            )
        terms[power // 2] = terms.get(power // 2, 0) + coefficient
    return sp.Poly(
        sum(coefficient * squared ** power
            for power, coefficient in terms.items()),
        squared,
    )


def layer_rates(mode, knudsen=sp.Integer(1), literal_source=False,
                scalar_moment_coefficient=None, omega_closure="thesis"):
    r"""
    Return the squared layer rates and the multiplicity of the bulk root.

    The returned mapping holds :math:`\\lambda^2` for every non-zero root; the
    root at zero is reported separately because it is the bulk solution, not
    a layer.
    """
    polynomial = characteristic_polynomial(
        mode, knudsen=knudsen, literal_source=literal_source,
        scalar_moment_coefficient=scalar_moment_coefficient,
        omega_closure=omega_closure,
    )
    roots = sp.roots(polynomial)
    rates = {}
    bulk_multiplicity = 0
    for root, multiplicity in roots.items():
        if root == 0:
            bulk_multiplicity = multiplicity
            continue
        rates[sp.simplify(root)] = multiplicity
    return rates, bulk_multiplicity
