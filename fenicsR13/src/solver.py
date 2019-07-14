# pylint: disable=invalid-name

"solver module"

import dolfin as df
import ufl

class Solver:
    "Solver class"
    def __init__(self, p, mesh):
        "Initializes solver"
        self.params = p
        self.mesh = mesh.mesh
        self.boundaries = mesh.boundaries
        self.cell = self.mesh.ufl_cell()
        self.mode = p["mode"]
        self.use_coeffs = p["use_coeffs"]
        self.tau = p["tau"]
        self.xi_tilde = p["xi_tilde"]
        self.use_cip = self.params["stabilization"]["cip"]["enable"]
        self.delta_1 = self.params["stabilization"]["cip"]["delta_1"]
        self.theta_w_inner = self.params["theta_w_inner"]
        self.theta_w_outer = self.params["theta_w_outer"]
        self.heat_source = df.Expression(self.params["heat_source"], degree=2)
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
        e_type = "Lagrange"
        msh = self.mesh
        for var in self.elems:
            if self.var_ranks[var] == 0:
                self.elems[var] = df.FiniteElement(e_type, cell, degree=1)
            elif self.var_ranks[var] == 1:
                self.elems[var] = df.VectorElement(e_type, cell, degree=1)
            elif self.var_ranks[var] == 2:
                self.elems[var] = df.TensorElement(e_type, cell, degree=1, symmetry=True)
            self.fspaces[var] = df.FunctionSpace(msh, self.elems[var])

        heat_elems = [self.elems["theta"], self.elems["s"]]
        self.mxd_elems["heat"] = df.MixedElement([self.elems["theta"], self.elems["s"]])
        self.mxd_fspaces["heat"] = df.FunctionSpace(msh, self.mxd_elems["heat"])

        stress_elems = [self.elems["p"], self.elems["u"], self.elems["sigma"]]
        self.mxd_elems["stress"] = df.MixedElement(stress_elems)
        self.mxd_fspaces["stress"] = df.FunctionSpace(msh, self.mxd_elems["stress"])

        self.mxd_elems["coupled"] = df.MixedElement(heat_elems, stress_elems)
        self.mxd_fspaces["coupled"] = df.FunctionSpace(msh, self.mxd_elems["coupled"])

    def assemble(self):
        "Assemble system"

        if self.mode == "heat":

            w = self.mxd_fspaces["heat"]
            mesh = self.mesh
            boundaries = self.boundaries

            tau = self.tau
            xi_tilde = self.xi_tilde
            theta_w_inner = self.theta_w_inner
            theta_w_outer = self.theta_w_outer

            # Define trial and testfunction
            (theta, s) = df.TrialFunctions(w)
            (kappa, r) = df.TestFunctions(w)

            # Define custom measeasure for boundaries
            df.ds = df.Measure('ds', domain=mesh, subdomain_data=boundaries)
            df.dS = df.Measure('dS', domain=mesh, subdomain_data=boundaries)

            # Normal and tangential components
            # => tangential (tx,ty) = (-ny,nx) = perp(n) only for 2D
            n = df.FacetNormal(mesh)
            t = ufl.perp(n)
            s_n = df.dot(s, n)
            r_n = df.dot(r, n)
            s_t = df.dot(s, t)
            r_t = df.dot(r, t)

            # Define source function
            f = self.heat_source

            def dev3d(mat):
                "2d deviatoric part of actually 3d matrix"
                return (
                    0.5 * (mat + ufl.transpose(mat))
                    - (1/3) * ufl.tr(mat) * ufl.Identity(2)
                )

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
                delta_1 = self.delta_1
                h = df.CellDiameter(mesh)
                h_avg = (h('+') + h('-'))/2.0  # pylint: disable=not-callable
                stab = - (delta_1 * h_avg**3 * df.jump(df.grad(theta), n)
                          * df.jump(df.grad(kappa), n)) * df.dS
            else:
                stab = 0

            self.form_a = a1 + a2 + stab
            self.form_b = l1 + l2

    def solve(self):
        "Solves the system"
        if self.mode == "heat":

            w = self.mxd_fspaces["heat"]
            sol = df.Function(w)
            df.solve(self.form_a == self.form_b, sol, [], solver_parameters={'linear_solver': 'direct'})

            (self.sol["theta"], self.sol["s"]) = sol.split()

    def load_exact_solution(self):
        "Writes exact solution"
        if self.mode == "heat":

            with open(self.exact_solution, 'r') as file:
                exact_solution_cpp_code = file.read()

            exact_solution = df.compile_cpp_code(exact_solution_cpp_code)

            self.esol["theta"] = df.CompiledExpression(exact_solution.Temperature(), degree=2)

            self.esol["s"] = df.CompiledExpression(exact_solution.Heatflux(), degree=2)

    def calc_errors(self):
        "Calculate errors"

        def calc_scalarfield_errors(sol_, sol_e_, v_theta_, name_, p_):
            "TODO"

            of = self.output_folder

            field_e_i = df.interpolate(sol_e_, v_theta_)
            field_i = df.interpolate(sol_, v_theta_)

            difference = df.project(sol_e_ - sol_, v_theta_)
            difference.rename("difference_{}".format(name_),
                              "difference_{}".format(name_))
            file_difference = df.File(of + "difference_{}_{}.pvd".format(name_, p_))
            file_difference.write(difference)

            err_f_L2 = df.errornorm(sol_e_, sol_, 'L2')
            err_v_linf = df.norm(field_e_i.vector()-field_i.vector(), 'linf')
            print("L_2 error:", err_f_L2)
            print("l_inf error:", err_v_linf)

            field_e_i.rename("{}_e_i".format(name_), "{}_e_i".format(name_))
            file_field_e = df.File(of + "{}_e.pvd".format(name_))
            file_field_e.write(field_e_i)

            field_i.rename("{}_i".format(name_), "{}_i".format(name_))
            file_field = df.File(of + "{}_i.pvd".format(name_))
            file_field.write(field_i)

            return (err_f_L2, err_v_linf)

        def calc_vectorfield_errors(sol_, sol_e_, v_sol, mesh_, name_, p_):
            "TODO"

            of = self.output_folder

            # Vector values functions interpolated
            field_e_i = df.interpolate(sol_e_, v_sol)
            field_i = df.interpolate(sol_, v_sol)

            difference = df.project(sol_e_ - sol_, v_sol)
            difference.rename("difference_{}".format(name_),
                              "difference_{}".format(name_))
            file_difference = df.File(of + "difference_{}_{}.pvd".format(name_, p_))
            file_difference.write(difference)


            dim = field_i.geometric_dimension()
            errs_f_L2 = [df.errornorm(
                field_e_i.split()[i], field_i.split()[i], 'L2'
            ) for i in range(dim)] # ignore warning
            errs_v_linf = [df.norm(
                field_e_i.split()[i].vector()-field_i.split()[i].vector(), 'linf'
            ) for i in range(dim)]
            print("L_2 error:", errs_f_L2)
            print("l_inf error:", errs_v_linf)

            field_e_i.rename("{}_e_i".format(name_), "{}_e_i".format(name_))
            file_field_e = df.File(of + "{}_e.pvd".format(name_))
            file_field_e.write(field_e_i)

            field_i.rename("{}_i".format(name_), "{}_i".format(name_))
            file_field = df.File(of + "{}_i.pvd".format(name_))
            file_field.write(field_i)

            return (errs_f_L2, errs_v_linf)

        if self.mode == "heat":
            (self.errors["f"]["l2"]["s"], self.errors["v"]["linf"]["s"]) = calc_vectorfield_errors(
                self.sol["s"], self.esol["s"], self.fspaces["s"], self.mesh,
                "s", 1
            )
            (self.errors["f"]["l2"]["theta"], self.errors["v"]["linf"]["theta"]) = calc_scalarfield_errors(
                self.sol["theta"], self.esol["theta"], self.fspaces["theta"],
                "theta", 1
            )

    def write_solutions(self):
        "Write Solutions"
        sols = self.sol
        for field in sols:
            sols[field].rename(field, field)
            file = df.File(self.output_folder + field + ".pvd")
            file.write(sols[field])

    def write_parameters(self):
        "Write Parameters: Heat source"
        f_heat = df.interpolate(
            self.heat_source,
            df.FunctionSpace(
                self.mesh,
                df.FiniteElement("Lagrange", degree=1, cell=self.cell)
            )
        )
        f_heat.rename("f_heat", "f_heat")
        with df.XDMFFile(self.output_folder + 'f_heat.xdmf') as file:
            file.write(f_heat)

    def write(self):
        "Writes to Paraview format"







    def todo(self):
        "todos"
        # if save_matrix:
        #     np.savetxt("A.mat", df.assemble(a_).array())
        #     # Use in matrix with:
        #     # >> T = readtable("a.txt");
        #     # >> M=table2array(T);
        #     # >> spy(M);
        #     # >> cond(M);
        #     # >> det(M);
        #     # >> svd(M);
