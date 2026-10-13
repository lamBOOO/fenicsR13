# pylint: disable=invalid-name

"""
Verify that a shipped exact solution actually solves the governing equations.

The convergence tests compare the finite element solution against an exact
solution that was derived with Mathematica and exported to C++. Nothing in the
test suite checks the exported solution itself, so a mistake in the derivation
or in the export would show up as a converging but wrong reference. This module
closes that gap: the fields of an exact solution are evaluated directly and
inserted into the balance equations, and the remaining residual is reported.

The residuals are expected to be at the level of the working precision. Because
the derivatives are evaluated numerically, an exactly satisfied equation does
not give a zero residual but one that shrinks when the precision is increased.
Running the check at several precisions therefore separates "satisfied exactly"
from "satisfied approximately": only the former scales with the precision.

Equations (steady state, no sources, dimensionless), see [PEL2026]_
(3.14)-(3.22):

.. math::
    \\nabla \\cdot u &= 0 \\\\
    \\nabla p + \\nabla \\cdot \\sigma &= 0 \\\\
    \\nabla \\cdot u + \\nabla \\cdot s &= 0 \\\\
    8 \\nabla \\cdot s + \\nabla \\cdot \\Omega &= - \\frac{2}{3 Kn} R

with the closure :math:`\\Omega = -\\frac{7}{3} Kn \\nabla R`, which turns the
last equation into

.. math::
    \\Delta R = \\frac{2}{7 Kn^2} R ,

a modified Helmholtz equation. Its decay length is the Knudsen layer scale
:math:`Kn / \\lambda` with :math:`\\lambda = \\sqrt{2/7}`, which is the value
hard-coded as ``lambda_1`` in the R14 exact solutions.

The heat flux balance is checked as well; it needs the fourth-order tensor
closure :math:`R_{ij} = -\\frac{24}{5} Kn (\\nabla s)_{stf}` and is the
equation in which R14 differs structurally from R13, through the extra term
:math:`\\frac{1}{6} \\nabla R`.

The stress balance is checked too,

.. math::
    2 (\\nabla u)_{stf} + \\frac{4}{5} (\\nabla s)_{stf}
    + \\frac{1}{Kn} \\sigma + \\nabla \\cdot m = 0 ,

with the rank-three closure :math:`m = -2 Kn (\\nabla \\sigma)_{stf}`. The
symmetric trace-free part of a rank-three tensor is awkward in polar
components, so this one equation is evaluated in Cartesian components, where
the symmetrisation and the divergence are plain index operations. The gas is
three dimensional while the geometry is planar, so derivatives along z vanish
and :math:`\\sigma_{zz} = -(\\sigma_{xx} + \\sigma_{yy})`, which is the same
convention the exported solution uses. Because the divergence of ``m`` needs
second derivatives of the stress, those are taken directly as multivariate
derivatives rather than by nesting one numerical derivative inside another,
which would cost half of the available digits.

Requires ``mpmath``; it does not import dolfin and therefore runs outside the
FEniCS container as well.

.. [PEL2026] Pelzer, J. (2026). Derivation and Implementation of the R14
    equations. Bachelor thesis, RWTH Aachen University.
"""


from mpmath import (
    mp, mpf, atan2, besseli, besselk, cos, sin, log, sqrt, diff
)

# Radii of the annulus used by the convergence test cases, see the mesh
# generator in ``tests/2d_mesh/ring.geo``.
R_INNER = mpf("0.5")
R_OUTER = mpf("2")

# Sample points (radius, angle) inside the annulus. They avoid the boundaries
# so that the numerical derivatives stay inside the domain.
SAMPLE_POINTS = [
    (mpf("0.7"), mpf("0.3")),
    (mpf("1.0"), mpf("0.0")),
    (mpf("1.3"), mpf("1.1")),
    (mpf("1.8"), mpf("2.4")),
]


def _I(order, x):
    return besseli(order, x)


def _K(order, x):
    return besselk(order, x)


def divergence(vector_field, radius, angle):
    r"""
    Return the divergence of a vector field in polar coordinates.

    The field is passed as a callable returning the physical components
    :math:`(v_r, v_\varphi)`.

    .. math::
        \nabla \cdot v = \frac{1}{r} \partial_r (r v_r)
                       + \frac{1}{r} \partial_\varphi v_\varphi
    """
    radial = diff(lambda r: r * vector_field(r, angle)[0], radius)
    angular = diff(lambda phi: vector_field(radius, phi)[1], angle)
    return radial / radius + angular / radius


def laplacian(scalar_field, radius, angle):
    r"""
    Return the Laplacian of a scalar field in polar coordinates.

    .. math::
        \Delta f = \frac{1}{r} \partial_r (r \partial_r f)
                 + \frac{1}{r^2} \partial^2_\varphi f
    """
    radial = diff(
        lambda r: r * diff(lambda rr: scalar_field(rr, angle), r), radius
    )
    angular = diff(lambda phi: scalar_field(radius, phi), angle, 2)
    return radial / radius + angular / radius ** 2


def divergence_tensor(tensor_field, radius, angle):
    r"""
    Return the divergence of a symmetric 2-tensor in polar coordinates.

    The field is passed as a callable returning
    :math:`(\sigma_{rr}, \sigma_{r\varphi}, \sigma_{\varphi\varphi})`. The
    component :math:`\sigma_{zz}` does not enter the in-plane divergence.

    .. math::
        (\nabla \cdot \sigma)_r &= \partial_r \sigma_{rr}
            + \frac{1}{r} \partial_\varphi \sigma_{r\varphi}
            + \frac{\sigma_{rr} - \sigma_{\varphi\varphi}}{r} \\
        (\nabla \cdot \sigma)_\varphi &= \partial_r \sigma_{r\varphi}
            + \frac{1}{r} \partial_\varphi \sigma_{\varphi\varphi}
            + \frac{2 \sigma_{r\varphi}}{r}
    """
    s_rr, s_rp, s_pp = tensor_field(radius, angle)
    d_rr = diff(lambda r: tensor_field(r, angle)[0], radius)
    d_rp_phi = diff(lambda phi: tensor_field(radius, phi)[1], angle)
    d_rp_r = diff(lambda r: tensor_field(r, angle)[1], radius)
    d_pp_phi = diff(lambda phi: tensor_field(radius, phi)[2], angle)
    return (
        d_rr + d_rp_phi / radius + (s_rr - s_pp) / radius,
        d_rp_r + d_pp_phi / radius + 2 * s_rp / radius,
    )


class R14AnnulusSolution:
    """
    Exact solution of the linear R14 system on an annulus, in polar form.

    The expressions are the ones exported to
    ``esols/1_coeffs_nosources_norot_inflow_positive.cpp``; the integration
    constants are read from that file so that this class and the C++ reference
    cannot drift apart. Every field is written as the axisymmetric mode plus
    the first angular mode, which is the ansatz of the derivation.
    """

    # Squares of the Knudsen layer decay rates. They are stored as exact
    # rationals and the square roots are taken at the working precision in use,
    # because a rate rounded at import time would put a floor under every
    # residual computed later at a higher precision.
    LAMBDA_SQUARED = {1: (2, 7), 2: (5, 6), 3: (3, 2), 4: (5, 9)}

    def __init__(self, constants, knudsen=1):
        self.constants = dict(constants)
        self._knudsen = knudsen

    @property
    def c(self):
        """Integration constants at the current working precision."""
        return {key: mpf(value) for key, value in self.constants.items()}

    @property
    def kn(self):
        """Knudsen number at the current working precision."""
        return mpf(self._knudsen)

    def lam(self, index):
        """Return a Knudsen layer decay rate at the current precision."""
        numerator, denominator = self.LAMBDA_SQUARED[index]
        return sqrt(mpf(numerator) / denominator)

    @property
    def LAMBDA_1(self):
        """Decay rate of the layer carried by the fourth-order moment."""
        return self.lam(1)

    @property
    def LAMBDA_2(self):
        return self.lam(2)

    @property
    def LAMBDA_3(self):
        return self.lam(3)

    @property
    def LAMBDA_4(self):
        return self.lam(4)

    def temperature(self, radius, angle):
        """Return the temperature deviation ``theta``."""
        c, kn = self.c, self.kn
        mode_0 = (
            c["C_4"]
            - c["CK_1"] * _I(0, self.LAMBDA_1 * radius / kn) / 15
            - 8 * c["CK_5"] * _I(0, self.LAMBDA_2 * radius / kn) / 15
            - c["CK_2"] * _K(0, self.LAMBDA_1 * radius / kn) / 15
            - 8 * c["CK_6"] * _K(0, self.LAMBDA_2 * radius / kn) / 15
            - 4 * c["C_2"] * log(radius / kn) / 15
        )
        mode_1 = (
            c["C_9"] * kn / radius
            + c["C_10"] * radius / kn
            - c["CK_3"] * _I(1, self.LAMBDA_1 * radius / kn) / 15
            - 16 * c["CK_7"] * _I(1, self.LAMBDA_2 * radius / kn) / 15
            - c["CK_4"] * _K(1, self.LAMBDA_1 * radius / kn) / 15
            - 16 * c["CK_8"] * _K(1, self.LAMBDA_2 * radius / kn) / 15
        )
        return mode_0 + cos(angle) * mode_1

    def scalar_moment(self, radius, angle):
        """Return the fourth-order scalar moment ``R`` (``Delta`` in R13)."""
        c, kn = self.c, self.kn
        mode_0 = (
            c["CK_1"] * _I(0, self.LAMBDA_1 * radius / kn)
            + c["CK_2"] * _K(0, self.LAMBDA_1 * radius / kn)
        )
        mode_1 = (
            c["CK_3"] * _I(1, self.LAMBDA_1 * radius / kn)
            + c["CK_4"] * _K(1, self.LAMBDA_1 * radius / kn)
        )
        return mode_0 + cos(angle) * mode_1

    def heat_flux(self, radius, angle):
        """Return the heat flux components ``(s_r, s_phi)``."""
        c, kn = self.c, self.kn
        lam = self.LAMBDA_4
        radial_0 = c["C_2"] * kn / radius
        radial_1 = (
            -15 * c["C_10"] / 4
            - 3 * c["C_8"]
            - mpf("1.5") * c["C_7"] * kn ** 2 / radius ** 2
            + mpf("3.75") * c["C_9"] * kn ** 2 / radius ** 2
            + c["CK_15"] * (_I(0, lam * radius / kn) - _I(2, lam * radius / kn))
            + c["CK_16"] * (_K(0, lam * radius / kn) - _K(2, lam * radius / kn))
        )
        angular_0 = (
            c["CK_13"] * _I(1, lam * radius / kn)
            + c["CK_14"] * _K(1, lam * radius / kn)
        )
        angular_1 = (
            -15 * c["C_10"] / 4
            - 3 * c["C_8"]
            + mpf("1.5") * c["C_7"] * kn ** 2 / radius ** 2
            - mpf("3.75") * c["C_9"] * kn ** 2 / radius ** 2
            + c["CK_15"] * (_I(0, lam * radius / kn) + _I(2, lam * radius / kn))
            + c["CK_16"] * (_K(0, lam * radius / kn) + _K(2, lam * radius / kn))
        )
        return (
            radial_0 + cos(angle) * radial_1,
            angular_0 - sin(angle) * angular_1,
        )

    def velocity(self, radius, angle):
        """Return the velocity components ``(u_r, u_phi)``."""
        c, kn = self.c, self.kn
        lam = self.LAMBDA_4
        radial_0 = (
            mpf("0.5") * c["C_1"] * kn / radius
            - mpf("0.4") * c["C_2"] * kn / radius
        )
        radial_1 = (
            c["C_11"]
            + mpf("0.25") * c["C_6"] * kn ** 2 / radius ** 2
            - mpf("1.5") * c["C_9"] * kn ** 2 / radius ** 2
            - mpf("0.25") * c["C_8"] * radius ** 2 / kn ** 2
            + c["CK_15"] * (
                -2 * _I(0, lam * radius / kn) / 5
                + 2 * _I(2, lam * radius / kn) / 5
            )
            + c["CK_16"] * (
                -2 * _K(0, lam * radius / kn) / 5
                + 2 * _K(2, lam * radius / kn) / 5
            )
            + c["C_7"] * (
                mpf("0.5")
                + mpf(5) / 3 * kn ** 2 / radius ** 2
                - log(radius / kn) / 2
            )
        )
        angular_0 = (
            mpf("0.5") * c["C_0"] * kn / radius
            + c["C_5"] * radius / kn
            - 2 * c["CK_13"] * _I(1, lam * radius / kn) / 5
            - 2 * c["CK_14"] * _K(1, lam * radius / kn) / 5
        )
        angular_1 = (
            c["C_11"]
            - mpf("0.25") * c["C_6"] * kn ** 2 / radius ** 2
            + mpf("1.5") * c["C_9"] * kn ** 2 / radius ** 2
            - mpf("0.75") * c["C_8"] * radius ** 2 / kn ** 2
            + c["CK_15"] * (
                -2 * _I(0, lam * radius / kn) / 5
                - 2 * _I(2, lam * radius / kn) / 5
            )
            + c["CK_16"] * (
                -2 * _K(0, lam * radius / kn) / 5
                - 2 * _K(2, lam * radius / kn) / 5
            )
            + c["C_7"] * (
                -mpf(5) / 3 * kn ** 2 / radius ** 2 - log(radius / kn) / 2
            )
        )
        return (
            radial_0 + cos(angle) * radial_1,
            angular_0 - sin(angle) * angular_1,
        )

    def pressure(self, radius, angle):
        """Return the pressure deviation ``p``."""
        c, kn = self.c, self.kn
        lam = self.LAMBDA_2
        mode_0 = (
            c["C_3"]
            - 4 * c["CK_5"] * _I(0, lam * radius / kn) / 3
            - 4 * c["CK_6"] * _K(0, lam * radius / kn) / 3
        )
        mode_1 = (
            c["C_7"] * kn / radius
            - 2 * c["C_8"] * radius / kn
            - 8 * c["CK_7"] * _I(1, lam * radius / kn) / 3
            - 8 * c["CK_8"] * _K(1, lam * radius / kn) / 3
        )
        return mode_0 + cos(angle) * mode_1

    def stress(self, radius, angle):
        """Return ``(sigma_rr, sigma_rphi, sigma_phiphi)``."""
        c, kn = self.c, self.kn
        lam_2, lam_3 = self.LAMBDA_2, self.LAMBDA_3
        rr_0 = (
            c["C_1"] * kn ** 2 / radius ** 2
            + c["CK_5"] * (
                _I(0, lam_2 * radius / kn) / 3 + _I(2, lam_2 * radius / kn)
            )
            + c["CK_9"] * (
                -_I(0, lam_3 * radius / kn) + _I(2, lam_3 * radius / kn)
            )
            + c["CK_6"] * (
                _K(0, lam_2 * radius / kn) / 3 + _K(2, lam_2 * radius / kn)
            )
            + c["CK_10"] * (
                -_K(0, lam_3 * radius / kn) + _K(2, lam_3 * radius / kn)
            )
        )
        rr_1 = (
            c["C_6"] * kn ** 3 / radius ** 3
            + c["C_7"] * kn / radius
            + c["C_8"] * radius / kn
            + c["CK_7"] * (
                5 * _I(1, lam_2 * radius / kn) / 3 + _I(3, lam_2 * radius / kn)
            )
            + c["CK_11"] * (
                -_I(1, lam_3 * radius / kn) + _I(3, lam_3 * radius / kn)
            )
            + c["CK_8"] * (
                5 * _K(1, lam_2 * radius / kn) / 3 + _K(3, lam_2 * radius / kn)
            )
            + c["CK_12"] * (
                -_K(1, lam_3 * radius / kn) + _K(3, lam_3 * radius / kn)
            )
        )
        rp_0 = c["C_0"] * kn ** 2 / radius ** 2
        rp_1 = (
            c["C_6"] * kn ** 3 / radius ** 3
            - c["C_8"] * radius / kn
            + c["CK_7"] * (
                -_I(1, lam_2 * radius / kn) + _I(3, lam_2 * radius / kn)
            )
            + c["CK_11"] * (
                -_I(1, lam_3 * radius / kn) + _I(3, lam_3 * radius / kn)
            )
            + c["CK_8"] * (
                -_K(1, lam_2 * radius / kn) + _K(3, lam_2 * radius / kn)
            )
            + c["CK_12"] * (
                -_K(1, lam_3 * radius / kn) + _K(3, lam_3 * radius / kn)
            )
        )
        pp_0 = (
            c["C_1"] * kn ** 2 / radius ** 2
            + c["CK_5"] * (
                -_I(0, lam_2 * radius / kn) / 3 + _I(2, lam_2 * radius / kn)
            )
            + c["CK_9"] * (
                _I(0, lam_3 * radius / kn) + _I(2, lam_3 * radius / kn)
            )
            + c["CK_6"] * (
                -_K(0, lam_2 * radius / kn) / 3 + _K(2, lam_2 * radius / kn)
            )
            + c["CK_10"] * (
                _K(0, lam_3 * radius / kn) + _K(2, lam_3 * radius / kn)
            )
        )
        pp_1 = (
            c["C_6"] * kn ** 3 / radius ** 3
            + c["C_7"] * kn / radius
            + c["C_8"] * radius / kn
            + c["CK_7"] * (
                _I(1, lam_2 * radius / kn) / 3 + _I(3, lam_2 * radius / kn)
            )
            + c["CK_11"] * (
                3 * _I(1, lam_3 * radius / kn) + _I(3, lam_3 * radius / kn)
            )
            + c["CK_8"] * (
                _K(1, lam_2 * radius / kn) / 3 + _K(3, lam_2 * radius / kn)
            )
            + c["CK_12"] * (
                3 * _K(1, lam_3 * radius / kn) + _K(3, lam_3 * radius / kn)
            )
        )
        return (
            rr_0 + cos(angle) * rr_1,
            rp_0 + sin(angle) * rp_1,
            -(pp_0 + cos(angle) * pp_1),
        )


class R13AnnulusSolution:
    """
    Exact solution of the linear R13 system on the same annulus.

    Taken from ``tests/2d_r13/esols/
    1_coeffs_nosources_norot_inflow_positive.cpp``.
    The R13 system has no balance equation for the fourth-order moment; there
    it is closed by :math:`\\Delta = -12 Kn \\nabla \\cdot s`, which vanishes
    for this case because the heat flux is divergence free. Checking this
    solution with the same residuals as the R14 one makes the comparison
    between the two theories rest on two independently verified references.

    The exported file stores one set of constants per inflow weight; the set is
    selected by commenting blocks in and out. ``read_constants`` keeps the
    first uncommented set, which is the one the test cases use.
    """

    LAMBDA_SQUARED = {1: (5, 9), 2: (5, 6), 3: (3, 2)}

    def __init__(self, constants, knudsen=1):
        self.constants = dict(constants)
        self._knudsen = knudsen

    @property
    def c(self):
        return {key: mpf(value) for key, value in self.constants.items()}

    @property
    def kn(self):
        return mpf(self._knudsen)

    def lam(self, index):
        numerator, denominator = self.LAMBDA_SQUARED[index]
        return sqrt(mpf(numerator) / denominator)

    def temperature(self, radius, angle):
        """Return the temperature deviation ``theta``."""
        c, kn, lam_2 = self.c, self.kn, self.lam(2)
        mode_0 = (
            c["C_8"]
            + 2 * c["CK_4"] * _I(0, radius * lam_2 / kn) / 5
            + 2 * c["CK_3"] * _K(0, radius * lam_2 / kn) / 5
            - 4 * c["C_6"] * log(radius / kn) / 15
        )
        mode_1 = (
            -4 * c["C_2"] * radius / (15 * kn)
            + 8 * c["C_4"] * radius / (5 * kn)
            - 2 * c["C_1"] * kn / (15 * radius)
            - 4 * c["C_5"] * kn / (5 * radius)
            + 2 * c["CK_10"] * _I(1, radius * lam_2 / kn) / 5
            + 2 * c["CK_9"] * _K(1, radius * lam_2 / kn) / 5
        )
        return mode_0 + cos(angle) * mode_1

    def heat_flux(self, radius, angle):
        """Return the heat flux components ``(s_r, s_phi)``."""
        c, kn, lam_1 = self.c, self.kn, self.lam(1)
        # Argument of the Bessel functions, and the factor the first-order
        # ones are divided by; both appear in every term below.
        argument = radius * lam_1 / kn
        scale = radius * lam_1
        radial_0 = c["C_6"] * kn / radius
        radial_1 = (
            c["C_2"]
            - c["C_1"] * kn ** 2 / (2 * radius ** 2)
            - 5 * (
                2 * c["CK_1"] * kn * _I(1, argument) / scale
                + 2 * c["CK_2"] * kn * _K(1, argument) / scale
            ) / 2
        )
        angular_1 = (
            c["C_2"]
            + c["C_1"] * kn ** 2 / (2 * radius ** 2)
            - 5 * (
                c["CK_1"] * (_I(0, argument) + _I(2, argument))
                - c["CK_2"] * (_K(0, argument) + _K(2, argument))
            ) / 2
        )
        return (radial_0 + cos(angle) * radial_1, -sin(angle) * angular_1)

    def pressure(self, radius, angle):
        """Return the pressure deviation ``p``."""
        c, kn, lam_2 = self.c, self.kn, self.lam(2)
        mode_0 = (
            c["C_9"]
            + c["CK_4"] * _I(0, radius * lam_2 / kn)
            + c["CK_3"] * _K(0, radius * lam_2 / kn)
        )
        mode_1 = (
            4 * c["C_4"] * radius / kn
            - 2 * c["C_5"] * kn / radius
            + c["CK_10"] * _I(1, radius * lam_2 / kn)
            + c["CK_9"] * _K(1, radius * lam_2 / kn)
        )
        return mode_0 + cos(angle) * mode_1

    def velocity(self, radius, angle):
        """Return the velocity components ``(u_r, u_phi)``."""
        c, kn, lam_1 = self.c, self.kn, self.lam(1)
        radial_0 = c["C_7"] * kn / radius
        radial_1 = (
            c["C_0"]
            + c["C_4"] * radius ** 2 / (2 * kn ** 2)
            - c["C_3"] * kn ** 2 / (2 * radius ** 2)
            + 2 * c["CK_1"] * kn * _I(1, radius * lam_1 / kn) / (radius * lam_1)
            + 2 * c["CK_2"] * kn * _K(1, radius * lam_1 / kn) / (radius * lam_1)
            + c["C_5"] * (-mpf("0.5") + log(radius / kn))
        )
        angular_1 = (
            c["C_0"]
            + 3 * c["C_4"] * radius ** 2 / (2 * kn ** 2)
            + c["C_3"] * kn ** 2 / (2 * radius ** 2)
            + c["CK_1"] * (
                _I(0, radius * lam_1 / kn) + _I(2, radius * lam_1 / kn)
            )
            - c["CK_2"] * (
                _K(0, radius * lam_1 / kn) + _K(2, radius * lam_1 / kn)
            )
            + c["C_5"] * (mpf("0.5") + log(radius / kn))
        )
        return (radial_0 + cos(angle) * radial_1, -sin(angle) * angular_1)

    def stress(self, radius, angle):
        """Return ``(sigma_rr, sigma_rphi, sigma_phiphi)``."""
        c, kn = self.c, self.kn
        lam_2, lam_3 = self.lam(2), self.lam(3)
        # Bessel arguments and the factors the terms are divided by.
        # They repeat in every line below, which is what pushes the
        # expressions past the line limit.
        arg_2, arg_3 = radius * lam_2 / kn, radius * lam_3 / kn
        scale_2, scale_3 = radius * lam_2, radius * lam_3
        rr_0 = (
            4 * c["C_6"] * kn ** 2 / (5 * radius ** 2)
            + 2 * c["C_7"] * kn ** 2 / radius ** 2
            + c["CK_6"] * kn * _I(1, arg_3) / scale_3
            + c["CK_4"] * (
                -kn * _I(1, arg_2) / (2 * scale_2)
                - _I(2, arg_2)
            )
            + c["CK_5"] * kn * _K(1, arg_3) / scale_3
            + c["CK_3"] * (
                kn * _K(1, arg_2) / (2 * scale_2)
                - _K(2, arg_2)
            )
        )
        rr_1 = (
            -2 * c["C_4"] * radius / kn
            - 2 * c["C_5"] * kn / radius
            - 4 * c["C_1"] * kn ** 3 / (5 * radius ** 3)
            - 2 * (c["C_3"] - 64 * c["C_5"] / 15) * kn ** 3 / radius ** 3
            + c["CK_10"] * (
                -_I(1, arg_2)
                + 3 * kn * _I(2, arg_2) / (2 * scale_2)
            )
            + c["CK_7"] * kn * _I(2, arg_3) / scale_3
            + c["CK_9"] * (
                -_K(1, arg_2)
                - 3 * kn * _K(2, arg_2) / (2 * scale_2)
            )
            + c["CK_8"] * kn * _K(2, arg_3) / scale_3
        )
        rp_1 = (
            2 * c["C_4"] * radius / kn
            - 4 * c["C_1"] * kn ** 3 / (5 * radius ** 3)
            - 2 * (c["C_3"] - 64 * c["C_5"] / 15) * kn ** 3 / radius ** 3
            + 3 * c["CK_10"] * kn * _I(2, arg_2) / (2 * scale_2)
            + c["CK_7"] * kn * _I(2, arg_3) / scale_3
            - 3 * c["CK_9"] * kn * _K(2, arg_2) / (2 * scale_2)
            + c["CK_8"] * kn * _K(2, arg_3) / scale_3
        )
        pp_0 = (
            4 * c["C_6"] * kn ** 2 / (5 * radius ** 2)
            + 2 * c["C_7"] * kn ** 2 / radius ** 2
            + c["CK_4"] * (
                -_I(0, arg_2) / 2
                + 3 * kn * _I(1, arg_2) / (2 * scale_2)
            )
            + c["CK_6"] * (
                -_I(0, arg_3)
                + kn * _I(1, arg_3) / scale_3
            )
            + c["CK_3"] * (
                -_K(0, arg_2) / 2
                - 3 * kn * _K(1, arg_2) / (2 * scale_2)
            )
            + c["CK_5"] * (
                _K(0, arg_3)
                + kn * _K(1, arg_3) / scale_3
            )
        )
        pp_1 = (
            -2 * c["C_4"] * radius / kn
            - 2 * c["C_5"] * kn / radius
            - 4 * c["C_1"] * kn ** 3 / (5 * radius ** 3)
            - 2 * (c["C_3"] - 64 * c["C_5"] / 15) * kn ** 3 / radius ** 3
            + c["CK_7"] * (
                -_I(1, arg_3)
                + kn * _I(2, arg_3) / scale_3
            )
            + c["CK_10"] * (
                -kn * _I(2, arg_2) / (2 * scale_2)
                - _I(3, arg_2) / 2
            )
            + c["CK_8"] * (
                _K(1, arg_3)
                + kn * _K(2, arg_3) / scale_3
            )
            + c["CK_9"] * (
                kn * _K(2, arg_2) / (2 * scale_2)
                - _K(3, arg_2) / 2
            )
        )
        return (
            rr_0 + cos(angle) * rr_1,
            sin(angle) * rp_1,
            -(pp_0 + cos(angle) * pp_1),
        )

    def scalar_moment(self, radius, angle):
        r"""
        Return the fourth-order moment, which R13 closes instead of solving.

        :math:`\Delta = -12 Kn \nabla \cdot s`, and the heat flux of this case
        is divergence free, so the moment vanishes identically.
        """
        return -12 * self.kn * divergence(self.heat_flux, radius, angle)


def stf_gradient_of_vector(vector_field, radius, angle):
    r"""
    Return the symmetric trace-free gradient of a vector, in polar components.

    Needed for the closure of the fourth-order tensor
    :math:`R_{ij} = -\frac{24}{5} Kn (\nabla s)_{stf}`, see [PEL2026]_ (3.20).
    The trace is removed over three components because the gas is treated as
    three dimensional while the geometry is planar; the out-of-plane gradient
    vanishes, so :math:`(\nabla s)_{zz} = 0` before the trace is subtracted.

    Returns ``(rr, r_phi, phi_phi, zz)``.
    """
    s_r, s_phi = vector_field(radius, angle)
    d_sr_dr = diff(lambda r: vector_field(r, angle)[0], radius)
    d_sr_dphi = diff(lambda phi: vector_field(radius, phi)[0], angle)
    d_sphi_dr = diff(lambda r: vector_field(r, angle)[1], radius)
    d_sphi_dphi = diff(lambda phi: vector_field(radius, phi)[1], angle)

    # Physical components of the gradient in polar coordinates.
    g_rr = d_sr_dr
    g_rp = d_sr_dphi / radius - s_phi / radius
    g_pr = d_sphi_dr
    g_pp = d_sphi_dphi / radius + s_r / radius

    sym_rr = g_rr
    sym_rp = (g_rp + g_pr) / 2
    sym_pp = g_pp
    trace = sym_rr + sym_pp                     # the zz component is zero here
    return (
        sym_rr - trace / 3,
        sym_rp,
        sym_pp - trace / 3,
        -trace / 3,
    )


def divergence_of_R_tensor(solution, radius, angle):
    r"""
    Return the divergence of the closed fourth-order tensor :math:`R_{ij}`.

    With :math:`R_{ij} = -\frac{24}{5} Kn (\nabla s)_{stf}` the divergence is
    taken the same way as for the stress tensor. The out-of-plane component
    does not contribute to the in-plane divergence but does enter through the
    trace, which is why it is carried along in
    :func:`stf_gradient_of_vector`.
    """
    factor = -mpf(24) / 5 * solution.kn

    def tensor(r, phi):
        rr, rp, pp, _zz = stf_gradient_of_vector(solution.heat_flux, r, phi)
        return (factor * rr, factor * rp, factor * pp)

    return divergence_tensor(tensor, radius, angle)


def _delta(i, j):
    """Kronecker delta over the three spatial directions."""
    return mpf(1) if i == j else mpf(0)


# All six permutations of three indices, used by the rank-three symmetrisation.
_PERMUTATIONS_OF_THREE = (
    (0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0),
)


def symmetric_part_of_rank_three(tensor):
    """Return the fully symmetric part, the mean over all six permutations."""
    result = [[[mpf(0)] * 3 for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                index = (i, j, k)
                total = mpf(0)
                for permutation in _PERMUTATIONS_OF_THREE:
                    total += tensor[index[permutation[0]]][
                        index[permutation[1]]][index[permutation[2]]]
                result[i][j][k] = total / 6
    return result


def stf_rank_three(tensor):
    r"""
    Return the symmetric trace-free part of a rank-three tensor.

    .. math::
        A_{ijk}^{stf} = A_{(ijk)} - \\frac{1}{5}
        \\left( t_i \\delta_{jk} + t_j \\delta_{ik} + t_k \\delta_{ij} \\right),
        \\qquad t_i = A_{(ijj)}

    This is the convention of ``TraceFreeSym3`` in the Mathematica sources
    under ``R14analytical``, over three directions.
    """
    symmetric = symmetric_part_of_rank_three(tensor)
    trace = [sum(symmetric[i][j][j] for j in range(3)) for i in range(3)]
    result = [[[mpf(0)] * 3 for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                result[i][j][k] = symmetric[i][j][k] - (
                    trace[i] * _delta(j, k)
                    + trace[j] * _delta(i, k)
                    + trace[k] * _delta(i, j)
                ) / 5
    return result


def stf_rank_two(tensor):
    """
    Return the symmetric trace-free part of a rank-two tensor.

    The trace runs over three directions, matching ``TraceFreeSym2`` in the
    Mathematica sources.
    """
    symmetric = [
        [(tensor[i][j] + tensor[j][i]) / 2 for j in range(3)] for i in range(3)
    ]
    trace = symmetric[0][0] + symmetric[1][1] + symmetric[2][2]
    return [
        [symmetric[i][j] - trace * _delta(i, j) / 3 for j in range(3)]
        for i in range(3)
    ]


class CartesianView:
    """
    The fields of a solution in Cartesian components.

    The solutions are written in polar components, which is convenient for
    every equation except the stress balance: there the rank-three moment
    would drag in Christoffel symbols. In Cartesian components the same
    operations are plain index arithmetic. Out-of-plane derivatives vanish and
    the out-of-plane stress follows from the trace, exactly as in the exported
    solution.
    """

    def __init__(self, solution):
        self.solution = solution

    @staticmethod
    def polar(x, y):
        """Return the polar coordinates of a Cartesian point."""
        return sqrt(x * x + y * y), atan2(y, x)

    def _rotate_vector(self, field, x, y):
        radius, angle = self.polar(x, y)
        radial, tangential = field(radius, angle)
        c, s = cos(angle), sin(angle)
        return [radial * c - tangential * s, radial * s + tangential * c,
                mpf(0)]

    def velocity(self, x, y):
        """Velocity in Cartesian components."""
        return self._rotate_vector(self.solution.velocity, x, y)

    def heat_flux(self, x, y):
        """Heat flux in Cartesian components."""
        return self._rotate_vector(self.solution.heat_flux, x, y)

    def stress(self, x, y):
        """Stress in Cartesian components, including the out-of-plane entry."""
        radius, angle = self.polar(x, y)
        rr, rp, pp = self.solution.stress(radius, angle)
        c, s = cos(angle), sin(angle)
        xx = c * c * rr - 2 * c * s * rp + s * s * pp
        xy = c * s * (rr - pp) + (c * c - s * s) * rp
        yy = s * s * rr + 2 * c * s * rp + c * c * pp
        zz = -(rr + pp)
        return [[xx, xy, mpf(0)], [xy, yy, mpf(0)], [mpf(0), mpf(0), zz]]

    def gradient_of_vector(self, field, x, y):
        """Return ``(grad V)_ik = d_k V_i``; the z derivative vanishes."""
        result = [[mpf(0)] * 3 for _ in range(3)]
        for i in range(3):
            result[i][0] = diff(lambda value: field(value, y)[i], x)
            result[i][1] = diff(lambda value: field(x, value)[i], y)
        return result

    def second_gradient_of_stress(self, x, y):
        """
        Return ``G[i][j][k][l] = d_l d_k sigma_ij``.

        The second derivatives are taken as multivariate derivatives of one
        scalar component, not as a derivative of a derivative.
        """
        orders = {(0, 0): (2, 0), (0, 1): (1, 1), (1, 0): (1, 1),
                  (1, 1): (0, 2)}
        result = [[[[mpf(0)] * 3 for _ in range(3)] for _ in range(3)]
                  for _ in range(3)]
        for i in range(3):
            for j in range(3):
                def component(value_x, value_y, i=i, j=j):
                    return self.stress(value_x, value_y)[i][j]
                for (k, l), order in orders.items():
                    result[i][j][k][l] = diff(component, (x, y), order)
        return result


def divergence_of_rank_three_moment(view, kn, x, y):
    r"""
    Return :math:`(\\nabla \\cdot m)_{ij} = \\partial_k m_{ijk}`.

    With :math:`m = -2 Kn (\\nabla \\sigma)_{stf}` and a symmetrisation that
    has constant coefficients, the derivative commutes with the symmetrisation:

    .. math::
        \\partial_l m_{ijk}
        = -2 Kn \\left( \\partial_l \\partial_k \\sigma \\right)^{stf}_{ijk}

    so only second derivatives of the stress are needed.
    """
    second = view.second_gradient_of_stress(x, y)
    result = [[mpf(0)] * 3 for _ in range(3)]
    for axis in range(3):
        sliced = [[[second[i][j][k][axis] for k in range(3)]
                   for j in range(3)] for i in range(3)]
        moment_derivative = stf_rank_three(sliced)
        for i in range(3):
            for j in range(3):
                result[i][j] += -2 * kn * moment_derivative[i][j][axis]
    return result


def stress_balance_residual(solution, radius, angle, with_rank_three=True):
    r"""
    Return the residual and the term scale of the stress balance.

    .. math::
        2 (\\nabla u)_{stf} + \\frac{4}{5} (\\nabla s)_{stf}
        + \\frac{1}{Kn} \\sigma + \\nabla \\cdot m

    ``with_rank_three=False`` drops the rank-three term. It is used by a test
    that checks the term matters: if the balance held without it, the check
    would not be testing the closure at all.
    """
    view = CartesianView(solution)
    x, y = radius * cos(angle), radius * sin(angle)

    velocity_term = stf_rank_two(
        view.gradient_of_vector(view.velocity, x, y)
    )
    heat_flux_term = stf_rank_two(
        view.gradient_of_vector(view.heat_flux, x, y)
    )
    stress = view.stress(x, y)
    if with_rank_three:
        moment_term = divergence_of_rank_three_moment(
            view, solution.kn, x, y
        )
    else:
        moment_term = [[mpf(0)] * 3 for _ in range(3)]

    worst = mpf(0)
    scale = mpf(0)
    for i in range(3):
        for j in range(3):
            terms = (
                2 * velocity_term[i][j],
                mpf(4) / 5 * heat_flux_term[i][j],
                stress[i][j] / solution.kn,
                moment_term[i][j],
            )
            worst = max(worst, abs(sum(terms)))
            scale = max(scale, sum(abs(term) for term in terms))
    return worst, scale


EQUATIONS = (
    "mass", "energy", "scalar_moment", "momentum", "heat_flux", "stress",
)


def residuals(solution, points=None, equations=None):
    """
    Return the relative residual of each checked balance equation.

    The residuals are normalised with the size of the terms that make up the
    equation, so that they can be compared across equations and precisions.

    ``equations`` restricts the work to the named equations. The stress
    balance needs second derivatives of the stress and costs more than the
    other five together, so a caller that only looks at one equation should
    say so.
    """
    points = points or SAMPLE_POINTS
    wanted = set(equations or EQUATIONS)
    unknown = wanted - set(EQUATIONS)
    if unknown:
        raise ValueError("unknown equation(s): %s" % ", ".join(sorted(unknown)))
    tiny = mpf(10) ** (-mp.dps - 5)
    worst = {name: mpf(0) for name in EQUATIONS if name in wanted}

    for radius, angle in points:
        velocity = solution.velocity
        heat_flux = solution.heat_flux

        if "mass" in wanted or "energy" in wanted:
            div_u = divergence(velocity, radius, angle)
        if "mass" in wanted:
            scale_u = sum(abs(comp) for comp in velocity(radius, angle)) + tiny
            worst["mass"] = max(worst["mass"], abs(div_u) / scale_u)

        if "energy" in wanted:
            div_s = divergence(heat_flux, radius, angle)
            scale_s = sum(abs(comp) for comp in heat_flux(radius, angle)) + tiny
            worst["energy"] = max(worst["energy"], abs(div_u + div_s) / scale_s)

        # 8 div(s) + div(Omega) = -(2/(3 Kn)) R with Omega = -(7/3) Kn grad(R).
        # The divergence of the heat flux vanishes for this case, which reduces
        # the balance to a modified Helmholtz equation for R.
        if "scalar_moment" in wanted:
            laplace_r = laplacian(solution.scalar_moment, radius, angle)
            expected = (
                2 / (7 * solution.kn ** 2)
                * solution.scalar_moment(radius, angle)
            )
            worst["scalar_moment"] = max(
                worst["scalar_moment"],
                abs(laplace_r - expected) / (abs(expected) + tiny),
            )

        if not ("momentum" in wanted or "heat_flux" in wanted):
            if "stress" in wanted:
                stress_res, stress_scale = stress_balance_residual(
                    solution, radius, angle
                )
                worst["stress"] = max(
                    worst["stress"], stress_res / (stress_scale + tiny)
                )
            continue

        grad_p_r = diff(lambda r: solution.pressure(r, angle), radius)
        grad_p_phi = (
            diff(lambda phi: solution.pressure(radius, phi), angle) / radius
        )
        div_sigma_r, div_sigma_phi = divergence_tensor(
            solution.stress, radius, angle
        )
        scale_m = (
            abs(grad_p_r) + abs(div_sigma_r)
            + abs(grad_p_phi) + abs(div_sigma_phi) + tiny
        )
        if "momentum" in wanted:
            worst["momentum"] = max(
                worst["momentum"],
                (abs(grad_p_r + div_sigma_r) + abs(grad_p_phi + div_sigma_phi))
                / scale_m,
            )

        if "heat_flux" not in wanted:
            if "stress" in wanted:
                stress_res, stress_scale = stress_balance_residual(
                    solution, radius, angle
                )
                worst["stress"] = max(
                    worst["stress"], stress_res / (stress_scale + tiny)
                )
            continue

        # Heat flux balance, [PEL2026]_ (3.18) with the closure (3.20):
        #   5/2 grad(theta) + 1/6 grad(R) + div(sigma) + 1/2 div(R_ij)
        #     = -(2/(3 Kn)) s
        # Only solutions that carry a scalar moment R can satisfy this form;
        # for R13 the term 1/6 grad(R) is absent because R is a closure there.
        grad_th_r = diff(lambda r: solution.temperature(r, angle), radius)
        grad_th_phi = (
            diff(lambda phi: solution.temperature(radius, phi), angle) / radius
        )
        grad_R_r = diff(lambda r: solution.scalar_moment(r, angle), radius)
        grad_R_phi = (
            diff(lambda phi: solution.scalar_moment(radius, phi), angle)
            / radius
        )
        div_R_r, div_R_phi = divergence_of_R_tensor(solution, radius, angle)
        s_r, s_phi = heat_flux(radius, angle)

        res_r = (
            mpf(5) / 2 * grad_th_r
            + grad_R_r / 6
            + div_sigma_r
            + div_R_r / 2
            + 2 / (3 * solution.kn) * s_r
        )
        res_phi = (
            mpf(5) / 2 * grad_th_phi
            + grad_R_phi / 6
            + div_sigma_phi
            + div_R_phi / 2
            + 2 / (3 * solution.kn) * s_phi
        )
        scale_h = (
            abs(mpf(5) / 2 * grad_th_r) + abs(div_sigma_r) + abs(div_R_r / 2)
            + abs(2 / (3 * solution.kn) * s_r)
            + abs(mpf(5) / 2 * grad_th_phi) + abs(div_sigma_phi)
            + abs(div_R_phi / 2) + abs(2 / (3 * solution.kn) * s_phi)
            + tiny
        )
        worst["heat_flux"] = max(
            worst["heat_flux"], (abs(res_r) + abs(res_phi)) / scale_h
        )

        # Stress balance, [PEL2026]_ (3.17) with the closure (3.19). Evaluated
        # in Cartesian components, see the module docstring.
        if "stress" in wanted:
            stress_res, stress_scale = stress_balance_residual(
                solution, radius, angle
            )
            worst["stress"] = max(
                worst["stress"], stress_res / (stress_scale + tiny)
            )

    return worst


def read_constants(cpp_file):
    """
    Read the integration constants from an exported C++ exact solution.

    Only the first definition of each constant is used, so commented-out
    blocks for other parameter sets are ignored.
    """
    import re

    pattern = re.compile(
        r"^double\s+(C_\d+|CK_\d+)\s*=\s*([-0-9.eE+]+)\s*;", re.MULTILINE
    )
    constants = {}
    with open(cpp_file) as handle:
        for match in pattern.finditer(handle.read()):
            constants.setdefault(match.group(1), float(match.group(2)))
    return constants


def read_layer_rates(cpp_file):
    """
    Read the Knudsen layer rates from an exported exact solution.

    The R13 files write them symbolically as ``sqrt(5.0/9.0)``, the R14 files
    as decimals, so both spellings are accepted. Reading them instead of
    repeating them in the test is what makes the comparison meaningful: a
    hard-coded copy would keep matching after the export changed.
    """
    import re

    rates = {}
    with open(cpp_file) as handle:
        text = handle.read()
    for name, a, b in re.findall(
        r"double\s+(lambda_\d+)\s*=\s*sqrt\(\s*([0-9.]+)\s*/\s*([0-9.]+)\s*\)",
        text,
    ):
        rates[name] = sqrt(mpf(a) / mpf(b))
    for name, value in re.findall(
        r"double\s+(lambda_\d+)\s*=\s*([-0-9.eE+]+)\s*;", text
    ):
        rates.setdefault(name, mpf(value))
    return rates


def check_precision_scaling(solution, precisions=(15, 25, 35), points=None,
                            equations=None):
    """
    Return the residuals at several working precisions.

    An equation that holds exactly produces residuals that follow the working
    precision. One that holds only approximately produces residuals that stop
    decreasing.

    ``points`` and ``equations`` are passed through, which lets a caller run
    the expensive stress balance on fewer points than the cheap equations.
    """
    table = {}
    original = mp.dps
    try:
        for precision in precisions:
            mp.dps = precision
            table[precision] = residuals(
                solution, points=points, equations=equations
            )
    finally:
        mp.dps = original
    return table


def main():
    """Report the residuals of the shipped R14 exact solution."""
    import os

    here = os.path.dirname(os.path.abspath(__file__))
    cpp_file = os.path.join(
        here, "1_coeffs_nosources_norot_inflow_positive.cpp"
    )
    solution = R14AnnulusSolution(read_constants(cpp_file))

    table = check_precision_scaling(solution)
    precisions = sorted(table)
    equations = sorted(table[precisions[0]])

    header = "{:<16}".format("equation")
    header += "".join("{:>14}".format("dps=%d" % p) for p in precisions)
    print(header)
    for equation in equations:
        line = "{:<16}".format(equation)
        line += "".join(
            "{:>14.1e}".format(float(table[p][equation])) for p in precisions
        )
        print(line)

    derived = sqrt(mpf(2) / 7)
    rates = read_layer_rates(cpp_file)
    print()
    print("Knudsen layer rate from the R equation: %.15f" % float(derived))
    print("layer rates read from the exported solution:")
    for name in sorted(rates):
        print("  %-9s = %.15f  (difference to sqrt(2/7): %.1e)"
              % (name, float(rates[name]), float(abs(rates[name] - derived))))


if __name__ == "__main__":
    main()
