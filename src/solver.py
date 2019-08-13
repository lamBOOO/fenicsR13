# pylint: disable=invalid-name

"solver module"

import copy
import dolfin as df
import ufl
import numpy as np
import tensoroperations as to

class Solver:
    r"""
    Solver class

    ::

        # Example usage:
        params = Input("input.yml").dict
        msh = meshes.H5Mesh("mesh.h5")
        solver = Solver(params, msh, "0") # "0" means time=0

    Assembles and solves the linear system.

    .. math::
        \mathbf{A} \mathbf{x} = \mathbf{b}

    The system results from the two dimensional, linearized R13 equations
    [TOR2003]_.

    .. [TOR2003] H Struchtrup, M Torrilhon (2003). Regularization of Grad’s 13
       moment equations: derivation and linear analysis.

    :ivar params: parameter dict
    :ivar mesh: Dolfin mesh
    :ivar cell: ``ufl_cell()`` for internal usage

    """
    def __init__(self, params, mesh, time):
        """
        Initializes solver
        """

        self.params = params #: Doctest
        self.mesh = mesh.mesh
        self.boundaries = mesh.boundaries
        self.cell = self.mesh.ufl_cell()
        self.time = time
        self.mode = params["mode"]
        self.use_coeffs = params["use_coeffs"]
        self.tau = params["tau"]
        self.xi_tilde = params["xi_tilde"]
        self.use_cip = self.params["stabilization"]["cip"]["enable"]
        self.delta_1 = self.params["stabilization"]["cip"]["delta_1"]
        self.delta_2 = self.params["stabilization"]["cip"]["delta_2"]
        self.delta_3 = self.params["stabilization"]["cip"]["delta_3"]

        # Create boundary field and sources expressions
        self.bcs = copy.deepcopy(self.params["bcs"])
        for edge_id in self.bcs:
            for field in self.bcs[edge_id].keys():
                self.bcs[edge_id][field] = self.createMacroExpr(
                    self.bcs[edge_id][field]
                )
        self.heat_source = self.createMacroExpr(self.params["heat_source"])
        self.mass_source = self.createMacroExpr(self.params["mass_source"])

        self.exact_solution = self.params["convergence_study"]["exact_solution"]
        self.rescale_p = self.params["convergence_study"]["rescale_pressure"]
        self.relative_error = self.params["convergence_study"][
            "relative_error"
        ]
        self.output_folder = self.params["case_name"] + "/"
        self.var_ranks = {
            "theta": 0,
            "s": 1,
            "p": 0,
            "u": 1,
            "sigma": 2,
        }
        self.elems = {
            "theta": None,
            "s": None,
            "p": None,
            "u": None,
            "sigma": None,
        }
        self.fspaces = {
            "theta": None,
            "s": None,
            "p": None,
            "u": None,
            "sigma": None,
        }
        self.mxd_elems = {
            "heat": None,
            "stress": None,
            "r13": None,
        }
        self.mxd_fspaces = {
            "heat": None,
            "stress": None,
            "r13": None,
        }
        self.form_a = None
        self.form_b = None
        self.sol = {
            "theta": None,
            "s": None,
            "p": None,
            "u": None,
            "sigma": None,
        }
        self.esol = {
            "theta": None,
            "s": None,
            "p": None,
            "u": None,
            "sigma": None,
        }
        self.errors = {}

    def createMacroExpr(self, cpp_string):
        """
        Return a DOLFIN expression with predefined macros.
        These macros include:

        ============================ ======= =================================
        Name                         Macro   CPP Replacement
        ============================ ======= =================================
        Radius wrt. to :math:`(0,0)` ``R``   ``sqrt(pow(x[0],2)+pow(x[1],2))``
        Angle wrt. :math:`(0,0)`     ``phi`` ``atan2(x[1],x[0])``
        Knudsen number               ``tau`` ``self.tau``
        ============================ ======= =================================

        The following expressions are therefore equal:

        .. code-block:: python

            # expr1 is equal to expr2
            expr1 = self.createMacroExpr("R*cos(phi)")
            expr2 = dolfin.Expression(
                "R*cos(phi)",
                degree=2,
                R=dolfin.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2),
                phi=dolfin.Expression("atan2(x[1],x[0])", degree=2),
            )
        """
        R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        phi = df.Expression("atan2(x[1],x[0])", degree=2)
        tau = self.tau
        return df.Expression(
            str(cpp_string),
            degree=2,
            tau=tau,
            phi=phi,
            R=R
        )

    def setup_function_spaces(self):
        "Setup function spaces"
        cell = self.cell
        msh = self.mesh
        for var in self.elems:
            e = self.params["elements"][var]["shape"]
            deg = self.params["elements"][var]["degree"]
            if self.var_ranks[var] == 0:
                self.elems[var] = df.FiniteElement(e, cell, deg)
            elif self.var_ranks[var] == 1:
                self.elems[var] = df.VectorElement(e, cell, deg)
            elif self.var_ranks[var] == 2:
                self.elems[var] = df.TensorElement(e, cell, deg, symmetry=True)
            self.fspaces[var] = df.FunctionSpace(msh, self.elems[var])

        heat_elems = [self.elems["theta"], self.elems["s"]]
        self.mxd_elems["heat"] = df.MixedElement(heat_elems)
        self.mxd_fspaces["heat"] = df.FunctionSpace(
            msh, self.mxd_elems["heat"]
        )

        stress_elems = [self.elems["p"], self.elems["u"], self.elems["sigma"]]
        self.mxd_elems["stress"] = df.MixedElement(stress_elems)
        self.mxd_fspaces["stress"] = df.FunctionSpace(
            msh, self.mxd_elems["stress"]
        )

        r13_elems = heat_elems + stress_elems
        self.mxd_elems["r13"] = df.MixedElement(r13_elems)
        self.mxd_fspaces["r13"] = df.FunctionSpace(
            msh, self.mxd_elems["r13"]
        )

    def check_bcs(self):
        """
        Checks if all mesh boundaries have bcs prescribed.
        Raises an exception if something is wrong.
        """
        boundary_ids = self.boundaries.array()
        bcs_specified = list(self.bcs.keys())

        for edge_id in boundary_ids:
            if not edge_id in [0] + bcs_specified: # inner zero allowed
                raise Exception("Mesh edge id {} has no bcs!".format(edge_id))

    def assemble(self):
        r"""
        Assembles the weak forms of the either the decoupled heat system,
        the decoupled stress system or the whole coupled system.

        .. |Rt| mathmacro:: \underline{\underline{R}}
        .. |st| mathmacro:: \underline{s}

        **Heat**:

        .. math::
            -\frac{24}{5} \tau (\nabla \st)_{\mathrm{dev}} - \Rt &= 0 \\
            \frac{1}{2} \nabla \cdot \Rt + \frac{2}{3\tau} \st + \frac{5}{2}
            \nabla \theta &= 0 \\
            \nabla \cdot \st &= f \\

        **Stress**

        Includes i.a. the term

        .. math::

            (\nabla \underline{\underline{\sigma}})_{\mathrm{dev}} :
            \nabla \underline{\underline{\psi}}

        for :math:`\theta` and :math:`\st` with a given heat source :math:`f`
        and the Knudsen number :math:`\tau`.

        """

        # Check if all mesh boundaries have bcs presibed frm input
        self.check_bcs()

        # Get local variables
        mesh = self.mesh
        boundaries = self.boundaries
        bcs = self.bcs
        tau = self.tau
        xi_tilde = self.xi_tilde
        delta_1 = self.delta_1
        delta_2 = self.delta_2
        delta_3 = self.delta_3

        # Normal and tangential components
        # - tangential (tx,ty) = (-ny,nx) = perp(n) only for 2D
        n = df.FacetNormal(mesh)
        t = ufl.perp(n)

        # Define custom measeasure for boundaries
        df.ds = df.Measure("ds", domain=mesh, subdomain_data=boundaries)
        df.dS = df.Measure("dS", domain=mesh, subdomain_data=boundaries)

        h = df.CellDiameter(mesh)
        h_avg = (h("+") + h("-"))/2.0 # pylint: disable=not-callable

        # Setup function spaces
        w_heat = self.mxd_fspaces["heat"]
        w_stress = self.mxd_fspaces["stress"]
        w_r13 = self.mxd_fspaces["r13"]
        if self.mode == "r13":
            (theta, s, p, u, sigma) = df.TrialFunctions(w_r13)
            (kappa, r, q, v, psi) = df.TestFunctions(w_r13)
        else:
            # Pure heat or pure stress: setup all functions..
            (theta, s) = df.TrialFunctions(w_heat)
            (kappa, r) = df.TestFunctions(w_heat)
            (p, u, sigma) = df.TrialFunctions(w_stress)
            (q, v, psi) = df.TestFunctions(w_stress)

        # Setup projections
        s_n = df.dot(s, n)
        r_n = df.dot(r, n)
        s_t = df.dot(s, t)
        r_t = df.dot(r, t)
        sigma_nn = df.dot(sigma*n, n)
        psi_nn = df.dot(psi*n, n)
        sigma_tt = df.dot(sigma*t, t)
        psi_tt = df.dot(psi*t, t)
        sigma_nt = df.dot(sigma*n, t)
        psi_nt = df.dot(psi*n, t)

        # Setup source functions
        f_heat = self.heat_source
        f_mass = self.mass_source

        if self.mode == "r13":
            cpl = 1.0
        else:
            cpl = 0.0

        # Setup both weak forms
        if self.use_coeffs:
            a1 = (
                + 12/5 * tau * df.inner(to.dev3d2(df.grad(s)), df.grad(r))
                + 2/3 * (1/tau) * df.inner(s, r)
                - (5/2) * theta * df.div(r)
                + cpl * df.dot(df.div(sigma), r)
            ) * df.dx + (
                + (
                    + 5/(4*xi_tilde) * s_n
                    - cpl * 5/8 * sigma_nn
                ) * r_n
                + (
                    + 11/10 * xi_tilde * s_t
                    + cpl * 1/10 * xi_tilde * s_t
                    - cpl * 1/2 * sigma_nt
                ) * r_t
            ) * df.ds
            l1 = sum([
                - 5.0/2.0 * r_n * bcs[bc]["theta_w"] * df.ds(bc)
                for bc in bcs.keys()
            ])

            a2 = - (df.div(s) * kappa) * df.dx
            l2 = - (f_heat * kappa) * df.dx

            a3 = (
                + 2 * tau * df.inner(
                    to.dev3dOfSym3(to.grad3dOf2(to.gen3dTF2(sigma))),
                    to.grad3dOf2(to.gen3dTF2(psi))
                )
                + (1/tau) * df.inner(to.gen3dTF2(sigma), to.gen3dTF2(psi))
                - 2 * df.dot(u, df.div(df.sym(psi)))
                + cpl * 4/5 * df.inner(to.dev3d2(df.grad(s)), psi)
            ) * df.dx + (
                + (
                    + 21/10 * xi_tilde * sigma_nn
                    + cpl * 3/20 * xi_tilde * sigma_nn
                    - cpl * 3/10  * s_n
                ) * psi_nn
                + 2 * xi_tilde * (
                    (sigma_tt + (1/2)*sigma_nn)*(psi_tt + (1/2)*psi_nn)
                )
                + (
                    + (2/xi_tilde) * sigma_nt
                    - cpl * 2/5 * s_t
                ) * psi_nt
            ) * df.ds + 2 * sum([
                bcs[bc]["gamma_w"] * (p + sigma_nn) * psi_nn * df.ds(bc)
                for bc in bcs.keys()
            ])
            l3 = sum([
                - 2.0 * psi_nt * bcs[bc]["u_t_w"] * df.ds(bc)
                - 2.0 * (
                    - bcs[bc]["gamma_w"] * bcs[bc]["p_w"]
                    + bcs[bc]["u_n_w"]
                ) * psi_nn * df.ds(bc)
                for bc in bcs.keys()
            ])

            a4 = (
                + df.dot(df.div(sigma), v)
                + df.dot(df.grad(p), v)
            ) * df.dx
            l4 = + df.Constant(0) * df.div(v) * df.dx

            a5 = + (
                df.dot(u, df.grad(q))
            ) * df.dx - sum([
                bcs[bc]["gamma_w"] * (p + sigma_nn) * q * df.ds(bc)
                for bc in bcs.keys()
            ])
            l5 = - (f_mass * q) * df.dx + sum([
                (
                    - bcs[bc]["gamma_w"] * bcs[bc]["p_w"]
                    + bcs[bc]["u_n_w"]
                ) * q * df.ds(bc)
                for bc in bcs.keys()
            ])
        else:
            a1 = (
                tau * df.inner(to.dev3d2(df.grad(s)), df.grad(r))
                + (1/tau) * df.inner(s, r)
                - theta * df.div(r)
            ) * df.dx + (
                + 1/(xi_tilde) * s_n * r_n
                + xi_tilde * s_t * r_t
            ) * df.ds
            a2 = - (df.div(s) * kappa) * df.dx
            l1 = sum([
                - 1 * r_n * bcs[bc]["theta_w"] * df.ds(bc)
                for bc in bcs.keys()
            ])
            l2 = - (f_heat * kappa) * df.dx

        # stabilization
        if self.use_cip:
            stab_heat = - (
                delta_1 * h_avg**3 *
                df.jump(df.grad(theta), n) * df.jump(df.grad(kappa), n)
            ) * df.dS

            stab_stress = (
                + delta_2 * h_avg**3 *
                df.dot(df.jump(df.grad(u), n), df.jump(df.grad(v), n))
                - delta_3 * h_avg *
                df.jump(df.grad(p), n) * df.jump(df.grad(q), n)
            ) * df.dS
        else:
            stab_heat = 0
            stab_stress = 0

        # Combine all equations
        if self.mode == "heat":
            self.form_a = a1 + a2 + stab_heat
            self.form_b = l1 + l2
        elif self.mode == "stress":
            self.form_a = a3 + a4 + a5 + stab_stress
            self.form_b = l3 + l4 + l5
        elif self.mode == "r13":
            self.form_a = a1 + a2 + stab_heat + a3 + a4 + a5 + stab_stress
            self.form_b = l1 + l2 + l3 + l4 + l5

    def solve(self):
        """
        Solves the system.
        Some available solver options:

        .. code-block:: python

            # Some solver params
            solver_parameters={
                'linear_solver': 'gmres', 'preconditioner': 'ilu' # or
                'linear_solver': 'petsc', 'preconditioner': 'ilu' # or
                'linear_solver': 'direct' # or
                'linear_solver': 'mumps'
            }

        """
        if self.mode == "heat":
            w = self.mxd_fspaces["heat"]
        elif self.mode == "stress":
            w = self.mxd_fspaces["stress"]
        elif self.mode == "r13":
            w = self.mxd_fspaces["r13"]

        sol = df.Function(w)
        df.solve(
            self.form_a == self.form_b, sol, [],
            solver_parameters={"linear_solver": "mumps"}
        )

        if self.mode == "heat":
            (self.sol["theta"], self.sol["s"]) = sol.split()
        elif self.mode == "stress":
            (self.sol["p"], self.sol["u"], self.sol["sigma"]) = sol.split()
        elif self.mode == "r13":
            (
                self.sol["theta"], self.sol["s"],
                self.sol["p"], self.sol["u"], self.sol["sigma"]
            ) = sol.split()

        if self.mode == "stress" or self.mode == "r13":
            if self.rescale_p:
                # Scale pressure to have zero mean
                p_i = df.interpolate(self.sol["p"], self.fspaces["p"])
                mean_p_value = self.calc_sf_mean(p_i)
                mean_p_fct = df.Function(self.fspaces["p"])
                mean_p_fct.assign(df.Constant(mean_p_value))
                p_i.assign(p_i - mean_p_fct)
                self.sol["p"] = p_i

    def load_exact_solution(self):
        "Writes exact solution"
        if self.mode == "heat" or self.mode == "r13":

            with open(self.exact_solution, "r") as file:
                exact_solution_cpp_code = file.read()

            esol = df.compile_cpp_code(exact_solution_cpp_code)

            self.esol["theta"] = df.CompiledExpression(
                esol.Temperature(), degree=2
            )

            self.esol["s"] = df.CompiledExpression(
                esol.Heatflux(), degree=2
            )
        if self.mode == "stress" or self.mode == "r13":

            with open(self.exact_solution, "r") as file:
                exact_solution_cpp_code = file.read()

            esol = df.compile_cpp_code(exact_solution_cpp_code)

            self.esol["p"] = df.CompiledExpression(
                esol.Pressure(), degree=2
            )

            self.esol["u"] = df.CompiledExpression(
                esol.Velocity(), degree=2
            )

            self.esol["sigma"] = df.CompiledExpression(
                esol.Stress(), degree=2
            )

    def calc_sf_mean(self, scalar_function):
        """
        Calculates the mean of a scalar function.
        Relative errors is based per field component and the maximum value of
        the analytical solution. If the analytical solution is uniformly zero,
        then the absolute erorrs is used.
        (equivalent to setting the maximum to 1)

        .. code-block:: python

            np.set_printoptions(precision=16)
            print(mean) # precision is not soo nice, only 9 digits
            print(self.calc_sf_mean(self.sol["p"])) # in solve() has m. prec. hm

        .. note::

            The following does not work in parallel because the mean is
            then only local:

            .. code-block:: python

                mean = np.mean(scalar_function.vector()[:], dtype=np.float64)
        """
        #
        v = scalar_function.vector()
        mean = v.sum()/v.size()
        return mean


    def calc_errors(self):
        "Calculate errors"

        def calc_scalarfield_errors(sol_, sol_e_, v_sol, name_):
            "TODO"

            field_e_i = df.interpolate(sol_e_, v_sol)
            field_i = df.interpolate(sol_, v_sol)

            difference = df.project(sol_e_ - sol_, v_sol)
            self.write_xdmf("difference_{}".format(name_), difference)

            err_f_L2 = df.errornorm(field_e_i, field_i, "L2")
            err_v_linf = df.norm(field_e_i.vector()-field_i.vector(), "linf")

            if self.relative_error:
                max_esol = df.norm(field_e_i.vector(), "linf") or 1 # avoid 0
                err_f_L2 /= max_esol
                err_v_linf /= max_esol

            print(str(name_) + " L_2 error:", err_f_L2)
            print(str(name_) + " l_inf error:", err_v_linf)

            self.write_xdmf(name_ + "_e", field_e_i)

            return {
                "L_2": err_f_L2,
                "l_inf": err_v_linf,
            }

        def calc_vectorfield_errors(sol_, sol_e_, v_sol, name_):
            "TODO"

            field_e_i = df.interpolate(sol_e_, v_sol)
            field_i = df.interpolate(sol_, v_sol)

            difference = df.project(sol_e_ - sol_, v_sol)
            self.write_xdmf("difference_{}".format(name_), difference)

            dofs = len(field_e_i.split())
            errs_f_L2 = [df.errornorm(
                field_e_i.split()[i], field_i.split()[i], "L2"
            ) for i in range(dofs)] # ignore warning
            errs_v_linf = [
                np.max(
                    np.abs(
                        field_e_i.split()[i].compute_vertex_values()
                        - field_i.split()[i].compute_vertex_values()
                    )
                )
                for i in range(dofs)
            ]

            if self.relative_error:
                max_esols = [
                    np.max(
                        np.abs(field_e_i.split()[i].compute_vertex_values())
                    )
                    for i in range(dofs)
                ]
                errs_f_L2 = [x/y for x, y in zip(errs_f_L2, max_esols)]
                errs_v_linf = [x/y for x, y in zip(errs_v_linf, max_esols)]

            print(str(name_) + " L_2 error:", errs_f_L2)
            print(str(name_) + " l_inf error:", errs_v_linf)

            self.write_xdmf(name_ + "_e", field_e_i)

            return [{
                "L_2": errs_f_L2[i],
                "l_inf": errs_v_linf[i],
            } for i in range(dofs)]

        def calc_tensorfield_errors(sol_, sol_e_, v_sol, name_):
            "TODO"

            field_e_i = df.interpolate(sol_e_, v_sol)
            field_i = df.interpolate(sol_, v_sol)

            # difference = df.project(sol_e_ - sol_, v_sol) # different outpu
            difference = df.project(field_e_i - field_i, v_sol)
            self.write_xdmf("difference_{}".format(name_), difference)

            dofs = len(field_e_i.split())
            errs_f_L2 = [df.errornorm(
                field_e_i.split()[i], field_i.split()[i], "L2"
            ) for i in range(dofs)] # ignore warning
            errs_v_linf = [
                np.max(
                    np.abs(
                        field_e_i.split()[i].compute_vertex_values()
                        - field_i.split()[i].compute_vertex_values()
                    )
                )
                for i in range(dofs)
            ]

            if self.relative_error:
                max_esols = [
                    np.max(
                        np.abs(field_e_i.split()[i].compute_vertex_values())
                    )
                    for i in range(dofs)
                ]
                errs_f_L2 = [x/y for x, y in zip(errs_f_L2, max_esols)]
                errs_v_linf = [x/y for x, y in zip(errs_v_linf, max_esols)]

            print(str(name_) + " L_2 error:", errs_f_L2)
            print(str(name_) + " l_inf error:", errs_v_linf)

            self.write_xdmf(name_ + "_e", field_e_i)

            return [{
                "L_2": errs_f_L2[i],
                "l_inf": errs_v_linf[i],
            } for i in range(dofs)]

        if self.mode == "heat" or self.mode == "r13":
            se = calc_scalarfield_errors(
                self.sol["theta"], self.esol["theta"],
                self.fspaces["theta"], "theta"
            )
            ve = calc_vectorfield_errors(
                self.sol["s"], self.esol["s"],
                self.fspaces["s"], "s"
            )
            ers = self.errors
            ers["theta"] = se
            ers["sx"] = ve[0]
            ers["sy"] = ve[1]
        if self.mode == "stress" or self.mode == "r13":
            se = calc_scalarfield_errors(
                self.sol["p"], self.esol["p"],
                self.fspaces["p"], "p"
            )
            ve = calc_vectorfield_errors(
                self.sol["u"], self.esol["u"],
                self.fspaces["u"], "u"
            )
            te = calc_tensorfield_errors(
                self.sol["sigma"], self.esol["sigma"],
                self.fspaces["sigma"], "sigma"
            )
            ers = self.errors
            ers["p"] = se
            ers["ux"] = ve[0]
            ers["uy"] = ve[1]
            ers["sigmaxx"] = te[0]
            ers["sigmaxy"] = te[1]
            ers["sigmayy"] = te[2]

    def write_solutions(self):
        "Write Solutions"
        sols = self.sol
        for field in sols:
            if sols[field] is not None:
                self.write_xdmf(field, sols[field])

    def write_parameters(self):
        "Write Parameters: Heat source or Mass Source"

        el = "Lagrange"
        deg = 1

        # Heat source
        f_heat = df.interpolate(
            self.heat_source,
            df.FunctionSpace(
                self.mesh,
                df.FiniteElement(el, degree=deg, cell=self.cell)
            )
        )
        self.write_xdmf("f_heat", f_heat)

        # Mass source
        f_mass = df.interpolate(
            self.mass_source,
            df.FunctionSpace(
                self.mesh,
                df.FiniteElement(el, degree=deg, cell=self.cell)
            )
        )
        self.write_xdmf("f_mass", f_mass)


    def write_xdmf(self, name, field):
        "Writes a renamed field to XDMF format"
        filename = self.output_folder + name + "_" + str(self.time) + ".xdmf"
        with df.XDMFFile(self.mesh.mpi_comm(), filename) as file:

            for degree in range(5): # test until degree five
                # Writing symmetric tensors crashes.
                # Therefore project symmetric tensor in nonsymmetric space
                # This is only a temporary fix, see:
                # https://fenicsproject.discourse.group/t/...
                # ...writing-symmetric-tensor-function-fails/1136
                el_symm = df.TensorElement(
                    df.FiniteElement(
                        "Lagrange", df.triangle, degree+1
                    ), symmetry=True
                ) # symmetric tensor element
                el_sol = field.ufl_function_space().ufl_element()
                if el_sol == el_symm:
                    # Remove symmetry with projection
                    field = df.project(
                        field, df.TensorFunctionSpace(
                            self.mesh, "Lagrange", degree+1
                        )
                    )
                    break

            field.rename(name, name)
            file.write(field, self.time)

    def write_systemmatrix(self):
        """
        Writes system matrix. Can be used to analyze e.g. condition number with
        decreasing mesh sizes.

        Import to Matlab with:

        .. code-block:: matlab

            % First adapt to filename...
            T = readtable("a.txt");
            M = table2array(T);

        """

        np.savetxt(
            self.output_folder + "A_{}.txt".format(self.time),
            df.assemble(self.form_a).array()
        )
