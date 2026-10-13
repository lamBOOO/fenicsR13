r"""
Insert an exact solution into the wall boundary conditions.

Why this exists next to :mod:`verify_exact_solution`: the balance equations
cannot see the integration constants. Every homogeneous solution satisfies
them by construction, so scaling one of the ``CK`` constants produces another
solution of the same system and leaves the residuals of the balance equations
unchanged. Measured on the exported R14 solution, a relative change of 1e-5 in
a single constant moves the balance residuals by a factor of 1.8, which is
nothing, while it moves the residuals below by up to ten orders of magnitude.

The constants are fixed by the boundary conditions alone, so only the
boundary conditions can check them. Together the two files answer the full
question: the balance equations say that the ansatz and the closures are
right, and this file says that the solution is the one belonging to *this*
boundary value problem.

Provenance, so that nothing here is invented:

* the seven R14 conditions are ``R14analytical/R14BoundaryConditions.wls``;
  the six R13 conditions in ``R13BoundaryConditions.wls`` are the same
  equations without the ``RR`` terms and without the ``Omega`` equation.
  For the case considered here the R13 closure
  :math:`\Delta = -12 Kn \nabla \cdot s` vanishes identically, so dropping
  those terms is exactly what setting ``RR = 0`` does.
* the wall data and the sign convention are ``R14analytical/ExampleMain.wls``:
  ``R0 = 0.5``, ``R1 = 2``, ``chi2 = 1``, ``facBC = 1``, ``eps = 1e-3``,
  wall temperatures 1 and 2, ``vnW = cos(phi)`` and ``vtW = -sin(phi)`` on the
  outer wall, ``pW = -0.27 cos(phi)`` there, and ``nn[R0] = -1``,
  ``nn[R1] = +1`` with the same signs for the tangential direction. The same
  numbers appear in the input files of the convergence cases.

The linear system uses :math:`p = \rho + \theta`, so the first condition is
evaluated with ``p`` in place of ``rho + temperature``; a wrong reading here
would show up immediately as a residual of order one instead of 1e-13.

The module needs ``mpmath`` but not dolfin, so it runs outside the container.
"""

from mpmath import cos, diff, mp, mpf, sin

# The tests import this as ``esols.verify_boundary_conditions``; the module is
# also meant to be runnable on its own, and then the relative form has no
# package to resolve against.
try:
    from .verify_exact_solution import (
        CartesianView,
        R13AnnulusSolution,
        R14AnnulusSolution,
        read_constants,
        stf_gradient_of_vector,
        stf_rank_three,
    )
except ImportError:  # pragma: no cover - direct execution
    from verify_exact_solution import (
        CartesianView,
        R13AnnulusSolution,
        R14AnnulusSolution,
        read_constants,
        stf_gradient_of_vector,
        stf_rank_three,
    )

#: Names of the conditions, in the order of ``BCeqn0``. The seventh exists
#: only in R14; it is the balance for the fifth-order flux ``Omega``.
CONDITIONS = (
    "normal", "slip", "temperature_jump", "mnnn", "mntt", "Rnt", "Omegan",
)

R13_CONDITIONS = tuple(name for name in CONDITIONS if name != "Omegan")

INNER_RADIUS = mpf("0.5")
OUTER_RADIUS = mpf(2)

#: Angles at which the conditions are checked. Four values, none of them a
#: multiple of pi/2, so that no angular mode drops out by accident.
SAMPLE_ANGLES = (mpf(0), mpf("1.1"), mpf("2.7"), mpf("4.4"))

_CHI2 = mpf(1)
_FAC_BC = mpf(1)
_EPS = mpf(10) ** -3

#: Digits used while evaluating the conditions. Well above the precision of
#: the exported constants, so that the residual measures them and not the
#: arithmetic.
WORKING_PRECISION = 30


def wall_data(radius):
    """Return the wall values and the sign convention at ``radius``."""
    if radius == INNER_RADIUS:
        return {
            "normal_sign": mpf(-1),
            "tangential_sign": mpf(-1),
            "temperature": mpf(1),
            "eps_pressure": _EPS,
            "eps_velocity": mpf(1),
            "normal_velocity": lambda angle: mpf(0),
            "tangential_velocity": lambda angle: mpf(0),
            "pressure": lambda angle: mpf(0),
            "name": "inner",
        }
    return {
        "normal_sign": mpf(1),
        "tangential_sign": mpf(1),
        "temperature": mpf(2),
        "eps_pressure": mpf(1),
        "eps_velocity": _EPS,
        "normal_velocity": cos,
        "tangential_velocity": lambda angle: -sin(angle),
        "pressure": lambda angle: mpf("-0.27") * cos(angle),
        "name": "outer",
    }


def rank_three_moment(view, knudsen, x, y):
    r"""
    Return :math:`m_{ijk} = -2 Kn (\nabla \sigma)^{stf}_{ijk}` in Cartesian
    components.

    Cartesian, because the polar form of a rank-three gradient carries
    Christoffel terms. The divergence of the same moment in
    :mod:`verify_exact_solution` makes that choice too, and that route was
    checked against the polar one.
    """
    gradient = [[[mpf(0)] * 3 for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            def component(value_x, value_y, i=i, j=j):
                return view.stress(value_x, value_y)[i][j]
            gradient[i][j][0] = diff(lambda v: component(v, y), x)
            gradient[i][j][1] = diff(lambda v: component(x, v), y)
    traceless = stf_rank_three(gradient)
    return [[[-2 * knudsen * traceless[i][j][k] for k in range(3)]
             for j in range(3)] for i in range(3)]


def project_rank_three(moment, first, second, third):
    """Contract a rank-three tensor with three direction vectors."""
    total = mpf(0)
    for i in range(3):
        if first[i] == 0:
            continue
        for j in range(3):
            if second[j] == 0:
                continue
            for k in range(3):
                if third[k] == 0:
                    continue
                total += moment[i][j][k] * first[i] * second[j] * third[k]
    return total


def _wall_quantities(solution, view, radius, angle):
    """Return every moment the conditions need, in wall components."""
    data = wall_data(radius)
    normal_sign = data["normal_sign"]
    tangential_sign = data["tangential_sign"]
    mixed_sign = normal_sign * tangential_sign
    knudsen = solution.kn

    heat_r, heat_phi = solution.heat_flux(radius, angle)
    velocity_r, velocity_phi = solution.velocity(radius, angle)
    stress_rr, stress_rp, stress_pp = solution.stress(radius, angle)

    gradient_rr, gradient_rp, _pp, _zz = stf_gradient_of_vector(
        solution.heat_flux, radius, angle
    )
    fourth_factor = -mpf(24) / 5 * knudsen

    x, y = radius * cos(angle), radius * sin(angle)
    moment = rank_three_moment(view, knudsen, x, y)
    radial = [cos(angle), sin(angle), mpf(0)]
    azimuthal = [-sin(angle), cos(angle), mpf(0)]

    scalar_moment = solution.scalar_moment(radius, angle)
    omega_factor = -mpf(7) / 3 * knudsen

    return {
        "temperature": solution.temperature(radius, angle),
        "pressure": solution.pressure(radius, angle),
        "vn": normal_sign * velocity_r,
        "vt": tangential_sign * velocity_phi,
        "qn": normal_sign * heat_r,
        "qt": tangential_sign * heat_phi,
        "Sigmann": stress_rr,
        "Sigmant": mixed_sign * stress_rp,
        "Sigmatt": stress_pp,
        "Rnn": fourth_factor * gradient_rr,
        "Rnt": mixed_sign * fourth_factor * gradient_rp,
        "RR": scalar_moment,
        "mnnn": normal_sign * project_rank_three(
            moment, radial, radial, radial),
        "mnnt": tangential_sign * project_rank_three(
            moment, radial, radial, azimuthal),
        "mntt": normal_sign * project_rank_three(
            moment, radial, azimuthal, azimuthal),
        "Omegan": normal_sign * omega_factor * diff(
            lambda r: solution.scalar_moment(r, angle), radius),
    }


def condition_residuals(solution, radius, angle, view=None):
    """
    Return ``(residual, scale)`` for every condition at one wall point.

    The scale is the sum of the magnitudes of the terms that make up the
    condition, so that the relative residual says how much the terms fail to
    cancel rather than how large they happen to be.
    """
    view = view or CartesianView(solution)
    data = wall_data(radius)
    q = _wall_quantities(solution, view, radius, angle)

    delta_temperature = q["temperature"] - data["temperature"]
    delta_velocity_t = q["vt"] - data["tangential_velocity"](angle)
    # The right-hand side of the mnnn condition reappears in the mntt one.
    shared = ((-(mpf(2) / 5)) * delta_temperature + 7 * q["Sigmann"] / 5
              - 2 * q["Rnn"] / 25 - q["RR"] / 75)

    residual = {
        "normal": (
            -data["eps_velocity"] * (q["vn"] - data["normal_velocity"](angle))
            + data["eps_pressure"] * (q["pressure"] + _FAC_BC * q["Sigmann"]
                                      - data["pressure"](angle))
        ),
        "slip": -q["Sigmant"] + _CHI2 * (
            delta_velocity_t + q["qt"] / 5 + q["mnnt"]),
        "temperature_jump": -q["qn"] + _CHI2 * (
            2 * delta_temperature + q["Sigmann"] / 2
            + 2 * q["Rnn"] / 5 + q["RR"] / 15),
        "mnnn": -q["mnnn"] + _CHI2 * shared,
        "mntt": (-q["mntt"] - _CHI2 * shared / 2
                 + _CHI2 * (q["Sigmann"] / 2 + q["Sigmatt"])),
        "Rnt": -q["Rnt"] + _CHI2 * (
            -delta_velocity_t + 11 * q["qt"] / 5 - q["mnnt"]),
        "Omegan": -q["Omegan"] + _CHI2 * (
            -8 * delta_temperature - 2 * q["Sigmann"]
            - 8 * q["Rnn"] / 5 + 4 * q["RR"] / 3),
    }
    scale = {
        "normal": (abs(data["eps_velocity"] * q["vn"])
                   + abs(data["eps_velocity"] * data["normal_velocity"](angle))
                   + abs(data["eps_pressure"] * q["pressure"])
                   + abs(data["eps_pressure"] * q["Sigmann"])
                   + abs(data["eps_pressure"] * data["pressure"](angle))),
        "slip": (abs(q["Sigmant"]) + abs(q["vt"])
                 + abs(data["tangential_velocity"](angle))
                 + abs(q["qt"] / 5) + abs(q["mnnt"])),
        "temperature_jump": (abs(q["qn"]) + abs(2 * delta_temperature)
                             + abs(q["Sigmann"] / 2) + abs(2 * q["Rnn"] / 5)
                             + abs(q["RR"] / 15)),
        "mnnn": abs(q["mnnn"]) + abs(shared),
        "mntt": (abs(q["mntt"]) + abs(shared / 2) + abs(q["Sigmann"] / 2)
                 + abs(q["Sigmatt"])),
        "Rnt": (abs(q["Rnt"]) + abs(q["vt"])
                + abs(data["tangential_velocity"](angle))
                + abs(11 * q["qt"] / 5) + abs(q["mnnt"])),
        "Omegan": (abs(q["Omegan"]) + abs(8 * delta_temperature)
                   + abs(2 * q["Sigmann"]) + abs(8 * q["Rnn"] / 5)
                   + abs(4 * q["RR"] / 3)),
    }
    return residual, scale


def boundary_residuals(solution, conditions=None, angles=None):
    """
    Return the worst relative residual of every condition over both walls.

    ``conditions`` restricts the work; pass :data:`R13_CONDITIONS` for a
    solution of the R13 system, where the ``Omega`` balance does not exist.

    The working precision is fixed here so that what the residual reports is
    the quality of the exported constants and not the arithmetic: those
    constants are decimal numbers of double precision, so nothing can cancel
    below roughly 1e-16, and the evaluation must stay well away from that
    floor to make the statement meaningful.
    """
    wanted = tuple(conditions or CONDITIONS)
    worst = {name: mpf(0) for name in wanted}
    with mp.workdps(WORKING_PRECISION):
        view = CartesianView(solution)
        for radius in (INNER_RADIUS, OUTER_RADIUS):
            for angle in (angles or SAMPLE_ANGLES):
                residual, scale = condition_residuals(solution, radius, angle,
                                                      view=view)
                for name in wanted:
                    if scale[name] > 0:
                        worst[name] = max(worst[name],
                                          abs(residual[name]) / scale[name])
    return worst


def main():
    """Report the boundary residuals of both shipped exact solutions."""
    import os

    here = os.path.dirname(os.path.abspath(__file__))
    cases = (
        ("R14", R14AnnulusSolution,
         os.path.join(here, "1_coeffs_nosources_norot_inflow_positive.cpp"),
         CONDITIONS),
        ("R13", R13AnnulusSolution,
         os.path.join(here, os.pardir, os.pardir, "2d_r13", "esols",
                      "1_coeffs_nosources_norot_inflow_positive.cpp"),
         R13_CONDITIONS),
    )

    for label, factory, path, conditions in cases:
        solution = factory(read_constants(path))
        worst = boundary_residuals(solution, conditions=conditions)
        print("%s, %d constants" % (label, len(solution.constants)))
        for name in conditions:
            print("  %-18s %.3e" % (name, float(worst[name])))
        print()


if __name__ == "__main__":
    main()
