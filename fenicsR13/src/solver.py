# pylint: disable=invalid-name

"solver module"

import dolfin as df
import ufl

class Solver:
    "Solver class"
    def __init__(self, params, mesh, time):
        "Initializes solver"
        self.params = params
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
        self.theta_w_inner = self.params["theta_w_inner"]
        self.theta_w_outer = self.params["theta_w_outer"]
        self.v_t_inner = self.params["v_t_inner"]
        self.v_t_outer = self.params["v_t_outer"]
        self.heat_source = df.Expression(self.params["heat_source"], degree=2)
        self.mass_source = df.Expression(self.params["mass_source"], degree=2)
        self.exact_solution = self.params["convergence_study"]["exact_solution"]
        self.output_folder = self.params["output_folder"]
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
            "coupled": None,
        }
        self.mxd_fspaces = {
            "heat": None,
            "stress": None,
            "coupled": None,
        }
        self.form_a = None
        self.form_b = None
        self.sol = {
            "theta": None,
            "s": None,
        }
        self.esol = {
            "theta": None,
            "s": None,
            "p": None,
            "u": None,
            "sigma": None,
        }
        self.errors = {
            "f": {
                "l2": {
                    "theta": None,
                    "s": None,
                }
            },
            "v": {
                "linf": {
                    "theta": None,
                    "s": None,
                }
            }
        }

    def setup_function_spaces(self):
        "Setup function spaces"
        cell = self.cell
        e = "Lagrange"
        msh = self.mesh
        for var in self.elems:
            if self.var_ranks[var] == 0:
                self.elems[var] = df.FiniteElement(e, cell, 1)
            elif self.var_ranks[var] == 1:
                self.elems[var] = df.VectorElement(e, cell, 1)
            elif self.var_ranks[var] == 2:
                self.elems[var] = df.TensorElement(e, cell, 1, symmetry=True)
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

        self.mxd_elems["coupled"] = df.MixedElement(heat_elems, stress_elems)
        self.mxd_fspaces["coupled"] = df.FunctionSpace(
            msh, self.mxd_elems["coupled"]
        )

    def assemble(self):
        "Assemble system"

        # Special tensor functions for 3D problems on 2D domains
        def dev3d(mat):
            "2d deviatoric part of actually 3d matrix"
            return (
                0.5 * (mat + ufl.transpose(mat))
                - (1/3) * ufl.tr(mat) * ufl.Identity(2)
            )
        def devOfGrad2(rank2):
            "From Henning's book p232"
            i, j, k, r = ufl.indices(4)
            entry_ijk = (
                (1/3) * (
                    rank2[i, j].dx(k) + rank2[i, k].dx(j) + rank2[j, k].dx(i)
                )
                # ufl.sym(ufl.grad(rank2))
                - (1/15) * (
                    + (2 * rank2[i, r].dx(r) + rank2[r, r].dx(i)) * delta[j, k]
                    + (2 * rank2[j, r].dx(r) + rank2[r, r].dx(j)) * delta[i, k]
                    + (2 * rank2[k, r].dx(r) + rank2[r, r].dx(k)) * delta[i, j]
                )
            )
            tensor = ufl.as_tensor(entry_ijk, (i, j, k))
            return tensor

        # Local variables
        mesh = self.mesh
        boundaries = self.boundaries
        tau = self.tau
        xi_tilde = self.xi_tilde
        theta_w_inner = self.theta_w_inner
        theta_w_outer = self.theta_w_outer
        v_t_inner = self.v_t_inner
        v_t_outer = self.v_t_outer
        delta_1 = self.delta_1
        delta_2 = self.delta_2
        delta_3 = self.delta_3

        # Normal and tangential components
        # => tangential (tx,ty) = (-ny,nx) = perp(n) only for 2D
        n = df.FacetNormal(mesh)
        t = ufl.perp(n)

        # d x d identiy matrix to use for Kronecker delta
        delta = df.Identity(2)

        # Define custom measeasure for boundaries
        df.ds = df.Measure("ds", domain=mesh, subdomain_data=boundaries)
        df.dS = df.Measure("dS", domain=mesh, subdomain_data=boundaries)

        h = df.CellDiameter(mesh)
        h_avg = (h("+") + h("-"))/2.0  # pylint: disable=not-callable

        if self.mode == "heat":

            w = self.mxd_fspaces["heat"]

            # Define trial and testfunction
            (theta, s) = df.TrialFunctions(w)
            (kappa, r) = df.TestFunctions(w)

            s_n = df.dot(s, n)
            r_n = df.dot(r, n)
            s_t = df.dot(s, t)
            r_t = df.dot(r, t)

            # Define heat source function
            f = self.heat_source

            if self.use_coeffs:
                a1 = (
                    + 12/5 * tau * df.inner(dev3d(df.grad(s)), df.grad(r))
                    + 2/3 * (1/tau) * df.inner(s, r)
                    - (5/2) * theta * df.div(r)
                ) * df.dx + (
                    + 5/(4*xi_tilde) * s_n * r_n
                    + 11/10 * xi_tilde * s_t * r_t
                ) * df.ds
                a2 = - (df.div(s) * kappa) * df.dx
                l1 = (- 5.0/2.0 * r_n * theta_w_outer * df.ds(3100)
                      - 5.0/2.0 * r_n * theta_w_inner * df.ds(3000))
                l2 = - (f * kappa) * df.dx
            else:
                a1 = (
                    tau * df.inner(dev3d(df.grad(s)), df.grad(r))
                    + (1/tau) * df.inner(s, r)
                    - theta * df.div(r)
                ) * df.dx + (
                    + 1/(xi_tilde) * s_n * r_n
                    + xi_tilde * s_t * r_t
                ) * df.ds
                a2 = - (df.div(s) * kappa) * df.dx
                l1 = (-1 * r_n * theta_w_inner * df.ds(3000)
                      -1 * r_n * theta_w_outer * df.ds(3100))
                l2 = - (f * kappa) * df.dx

            # stabilization
            if self.use_cip:
                stab = - (delta_1 * h_avg**3 * df.jump(df.grad(theta), n)
                          * df.jump(df.grad(kappa), n)) * df.dS
            else:
                stab = 0

            self.form_a = a1 + a2 + stab
            self.form_b = l1 + l2

        elif self.mode == "stress":

            w = self.mxd_fspaces["stress"]

            # Define trial and testfunction
            (p, u, sigma) = df.TrialFunctions(w)
            (q, v, psi) = df.TestFunctions(w)

            sigma_nn = df.dot(sigma*n, n)
            psi_nn = df.dot(psi*n, n)
            sigma_tt = df.dot(sigma*t, t)
            psi_tt = df.dot(psi*t, t)
            sigma_nt = df.dot(sigma*n, t)
            psi_nt = df.dot(psi*n, t)

            # Define mass source function
            f = self.mass_source

            if self.use_coeffs:
                a1 = (
                    + 2 * tau * df.inner(devOfGrad2(sigma), df.grad(psi))
                    # + 2 * tau * d.inner(ufl.tr(d.grad(sigma_)), d.grad(psi_))
                    + (1/tau) * df.inner(sigma, psi)
                    - 2 * df.dot(u, df.div(df.sym(psi)))
                ) * df.dx + (
                    + 21/(10*xi_tilde) * sigma_nn * psi_nn
                    + 2 * xi_tilde * (
                        (sigma_tt + (1/2)*sigma_nn)*(psi_tt + (1/2)*psi_nn)
                    )
                    + (2/xi_tilde) * sigma_nt * psi_nt
                ) * df.ds
                # l1 = (- 2 * v_t_inner * df.sin(phi) * psi_nt * df.ds(3000)
                #       - 2 * v_t_outer * df.sin(phi) * psi_nt * df.ds(3100))
                l1 = (- 2 * v_t_inner * psi_nt * df.ds(3000)
                      - 2 * v_t_outer * psi_nt * df.ds(3100))
                a2 = (df.dot(df.div(sigma), v) + df.dot(df.grad(p), v)) * df.dx
                l2 = df.Constant(0) * df.div(v) * df.dx # dummy
                a3 = df.dot(u, df.grad(q)) * df.dx
                l3 = - (f * q) * df.dx

            if self.use_cip:
                stab = (
                    + delta_2 * h_avg**3
                    * df.dot(df.jump(df.grad(u), n), df.jump(df.grad(v), n))
                    - delta_3 * h_avg *
                    df.jump(df.grad(p), n) * df.jump(df.grad(q), n)
                    ) * df.dS
            else:
                stab = 0

            self.form_a = a1 + a2 + a3 + stab
            self.form_b = l1 + l2 + l3

    def solve(self):
        "Solves the system"
        if self.mode == "heat":

            w = self.mxd_fspaces["heat"]
            sol = df.Function(w)
            df.solve(
                self.form_a == self.form_b, sol, [],
                solver_parameters={"linear_solver": "direct"}
            )

            (self.sol["theta"], self.sol["s"]) = sol.split()
        elif self.mode == "stress":

            w = self.mxd_fspaces["stress"]
            sol = df.Function(w)
            df.solve(
                self.form_a == self.form_b, sol, [],
                solver_parameters={"linear_solver": "direct"}
            )

            (self.sol["p"], self.sol["u"], self.sol["sigma"]) = sol.split()

    def load_exact_solution(self):
        "Writes exact solution"
        if self.mode == "heat":

            with open(self.exact_solution, "r") as file:
                exact_solution_cpp_code = file.read()

            esol = df.compile_cpp_code(exact_solution_cpp_code)

            self.esol["theta"] = df.CompiledExpression(
                esol.Temperature(), degree=2
            )

            self.esol["s"] = df.CompiledExpression(
                esol.Heatflux(), degree=2
            )

    def calc_errors(self):
        "Calculate errors"

        def calc_scalarfield_errors(sol_, sol_e_, v_theta_, name_):
            "TODO"

            field_e_i = df.interpolate(sol_e_, v_theta_)
            field_i = df.interpolate(sol_, v_theta_)

            difference = df.project(sol_e_ - sol_, v_theta_)
            self.write_xdmf("difference_{}".format(name_), difference)

            err_f_L2 = df.errornorm(sol_e_, sol_, "L2")
            err_v_linf = df.norm(field_e_i.vector()-field_i.vector(), "linf")
            print("L_2 error:", err_f_L2)
            print("l_inf error:", err_v_linf)

            self.write_xdmf(name_ + "_e", field_e_i)

            return (err_f_L2, err_v_linf)

        def calc_vectorfield_errors(sol_, sol_e_, v_sol, name_):
            "TODO"

            # Vector values functions interpolated
            field_e_i = df.interpolate(sol_e_, v_sol)
            field_i = df.interpolate(sol_, v_sol)

            difference = df.project(sol_e_ - sol_, v_sol)
            self.write_xdmf("difference_{}".format(name_), difference)

            dim = field_i.geometric_dimension()
            errs_f_L2 = [df.errornorm(
                field_e_i.split()[i], field_i.split()[i], "L2"
            ) for i in range(dim)] # ignore warning
            errs_v_linf = [df.norm(
                field_e_i.split()[i].vector()
                -field_i.split()[i].vector(), "linf"
            ) for i in range(dim)]
            print("L_2 error:", errs_f_L2)
            print("l_inf error:", errs_v_linf)

            self.write_xdmf(name_ + "_e", field_e_i)

            return (errs_f_L2, errs_v_linf)

        if self.mode == "heat":
            ve = calc_vectorfield_errors(
                self.sol["s"], self.esol["s"],
                self.fspaces["s"], "s"
            )
            se = calc_scalarfield_errors(
                self.sol["theta"], self.esol["theta"],
                self.fspaces["theta"], "theta"
            )
            f_l2 = self.errors["f"]["l2"]
            v_linf = self.errors["v"]["linf"]
            (f_l2["s"], v_linf["s"]) = ve
            (f_l2["theta"], v_linf["theta"]) = se

    def write_solutions(self):
        "Write Solutions"
        sols = self.sol
        for field in sols:
            if sols[field] is not None:

                # Writing symmetric tensors crashes.
                # Therefore project symmetric tensor in nonsymmetric space
                # This is only a temporary fix, see:
                # https://fenicsproject.discourse.group/t/...
                # ...writing-symmetric-tensor-function-fails/1136
                el_symm = df.TensorElement(
                    df.FiniteElement('Lagrange', df.triangle, 1), symmetry=True
                ) # check for symmetric tensors
                el_sol = sols[field].ufl_function_space().ufl_element()
                if el_sol == el_symm:
                    # Remove symmetry with projection
                    sols[field] = df.project(
                        sols[field], df.TensorFunctionSpace(
                            self.mesh, "Lagrange", 1
                        )
                    )

                self.write_xdmf(field, sols[field])

    def write_parameters(self):
        "Write Parameters: Heat source or Mass Source"

        # Heat source
        f_heat = df.interpolate(
            self.heat_source,
            df.FunctionSpace(
                self.mesh,
                df.FiniteElement("Lagrange", degree=1, cell=self.cell)
            )
        )
        self.write_xdmf("f_heat", f_heat)

        # Mass source
        f_mass = df.interpolate(
            self.mass_source,
            df.FunctionSpace(
                self.mesh,
                df.FiniteElement("Lagrange", degree=1, cell=self.cell)
            )
        )
        self.write_xdmf("f_mass", f_mass)


    def write_xdmf(self, name, field):
        "Writes a renamed field to XDMF format"
        filename = self.output_folder + name + "_" + str(self.time) + ".xdmf"
        field.rename(name, name)
        with df.XDMFFile(self.mesh.mpi_comm(), filename) as file:
            file.write(field, self.time)

    def extract_matrix(self):
        """
        Extract system matrix after assemble.
        Import to Matlab with:
        $ T = readtable("a.txt");
        $ M=table2array(T);
        """
        raise NotImplementedError()
