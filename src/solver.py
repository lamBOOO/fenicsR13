# pylint: disable=invalid-name

"""
Solver module, contains the Solver class.

For usage examples, see the :class:`solver.Solver` description.
"""

import os
import copy
import time as time_module
import dolfin as df
import ufl
import numpy as np
import tensoroperations as to


class Solver:
    r"""
    Class to store the actual solver.

    Possible order of methods in context of convergence study
    (see main program):

    .. digraph:: foo

        "START, mesh=meshes[i=0]" ->
        "__init__" ->
        "assemble()" ->
        "solve()" ->
        "write()" ->
        "..." ->
        "mesh=meshes[i+1]" ->
        "__init__";
        "..." -> "END";

    Parameters
    ----------
    params : input.Input
        Input parameters
    mesh : meshes.H5Mesh
        Mesh
    time : string
        Identifier for e.g. convergence study naming

    Returns
    -------
    Solver
        Solver object

    Example
    -------
    >>> # Example usage:
    >>> from input import Input
    >>> from meshes import H5Mesh
    >>> params = Input(
    ...     "tests/heat/inputs/heat_01_coeffs_p1p1_stab.yml"
    ... ) # doctest: +ELLIPSIS
    Input:...
    >>> msh = H5Mesh("tests/mesh/ring0.h5")
    >>> solver = Solver(params.dict, msh, "0") # "0" means time=0

    """

    def __init__(self, params, mesh, time):
        """Initialize solver and setup variables from input parameters."""
        self.params = params #: Doctest
        self.mesh = mesh.mesh
        self.boundaries = mesh.boundaries
        self.cell = self.mesh.ufl_cell()
        self.time = time
        self.mode = params["mode"]
        self.use_coeffs = params["use_coeffs"]
        self.kn = params["kn"]
        self.xi_tilde = params["xi_tilde"]
        self.use_cip = self.params["stabilization"]["cip"]["enable"]
        self.delta_1 = self.params["stabilization"]["cip"]["delta_1"]
        self.delta_2 = self.params["stabilization"]["cip"]["delta_2"]
        self.delta_3 = self.params["stabilization"]["cip"]["delta_3"]

        self.write_pdfs = self.params["postprocessing"]["write_pdfs"]
        self.massflow = self.params["postprocessing"]["massflow"]

        # Create boundary field and sources expressions
        self.bcs = copy.deepcopy(self.params["bcs"])
        for edge_id in self.bcs:
            for field in self.bcs[edge_id].keys():
                self.bcs[edge_id][field] = self.__createMacroExpr(
                    self.bcs[edge_id][field]
                )
        self.heat_source = self.__createMacroExpr(self.params["heat_source"])
        self.mass_source = self.__createMacroExpr(self.params["mass_source"])

        self.exact_solution = self.params["convergence_study"]["exact_solution"]
        self.write_systemmatrix = self.params["convergence_study"][
            "write_systemmatrix"
        ]
        self.rescale_p = self.params["convergence_study"]["rescale_pressure"]
        self.relative_error = self.params["convergence_study"][
            "relative_error"
        ]
        self.output_folder = self.params["output_folder"] + "/"
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

    def __createMacroExpr(self, cpp_string):
        """
        Return a DOLFIN expression with predefined macros.

        These macros include:

        ============================ ======= =================================
        Name                         Macro   CPP Replacement
        ============================ ======= =================================
        Radius wrt. to :math:`(0,0)` ``R``   ``sqrt(pow(x[0],2)+pow(x[1],2))``
        Angle wrt. :math:`(0,0)`     ``phi`` ``atan2(x[1],x[0])``
        Knudsen number               ``kn``  ``self.kn``
        ============================ ======= =================================

        The following expressions are therefore equal:

        .. code-block:: python

            # expr1 is equal to expr2
            expr1 = self.__createMacroExpr("R*cos(phi)")
            expr2 = dolfin.Expression(
                "R*cos(phi)",
                degree=2,
                R=dolfin.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2),
                phi=dolfin.Expression("atan2(x[1],x[0])", degree=2),
            )
        """
        R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        phi = df.Expression("atan2(x[1],x[0])", degree=2)
        kn = self.kn
        return df.Expression(
            str(cpp_string),
            degree=2,
            kn=kn,
            phi=phi,
            R=R
        )

    def __setup_function_spaces(self):
        """
        Set up function spaces for trial and test functions for assembling.

        Depends on the ``mode``.
        Function spaces depend on the choice of the element and its degree
        (see the input file :class:`input.Input`).

        The following DOLFIN functions are used:

        ========= =================
        field     DOLFIN Function
        ========= =================
        ``theta`` ``FiniteElement``
        ``s``     ``VectorElement``
        ``p``     ``FiniteElement``
        ``u``     ``VectorElement``
        ``sigma`` ``TensorElement``
        ========= =================
        """
        # Setup elements for all fields
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

        # Bundle elements per mode into `mxd_elems` dict
        # 1) heat
        heat_elems = [self.elems["theta"], self.elems["s"]]
        self.mxd_elems["heat"] = df.MixedElement(heat_elems)
        self.mxd_fspaces["heat"] = df.FunctionSpace(
            msh, self.mxd_elems["heat"]
        )
        # 2) stress
        stress_elems = [self.elems["p"], self.elems["u"], self.elems["sigma"]]
        self.mxd_elems["stress"] = df.MixedElement(stress_elems)
        self.mxd_fspaces["stress"] = df.FunctionSpace(
            msh, self.mxd_elems["stress"]
        )
        # 3) r13
        r13_elems = heat_elems + stress_elems
        self.mxd_elems["r13"] = df.MixedElement(r13_elems)
        self.mxd_fspaces["r13"] = df.FunctionSpace(
            msh, self.mxd_elems["r13"]
        )

    def __check_bcs(self):
        """
        Check if all boundaries from the input mesh have BCs prescribed.

        Raises an exception if one BC is missing.
        """
        boundary_ids = self.boundaries.array()
        bcs_specified = list(self.bcs.keys())

        for edge_id in boundary_ids:
            if not edge_id in [0] + bcs_specified: # inner zero allowed
                raise Exception("Mesh edge id {} has no bcs!".format(edge_id))

    def assemble(self):
        r"""
        Assemble the weak form of the system, depending on the mode.

        The system results from the two dimensional, linearized R13 equations
        [TOR2003]_.

        .. [TOR2003] H Struchtrup, M Torrilhon (2003). Regularization of
            Grad's 13 moment equations: derivation and linear analysis.

        .. |Rt| mathmacro:: \underline{\underline{R}}
        .. |st| mathmacro:: \underline{s}

        **Heat**:

        .. math::
            -\frac{24}{5} \mathrm{Kn} (\nabla \st)_{\mathrm{STF}} - \Rt &= 0 \\
            \frac{1}{2} \nabla \cdot \Rt + \frac{2}{3\mathrm{Kn}} \st
            + \frac{5}{2}
            \nabla \theta &= 0 \\
            \nabla \cdot \st &= f \\

        **Stress**

        Includes i.a. the term

        .. math::

            (\nabla \underline{\underline{\sigma}})_{\mathrm{STF}} :
            \nabla \underline{\underline{\psi}}

        for :math:`\theta` and :math:`\st` with a given heat source :math:`f`
        and the Knudsen number :math:`\mathrm{Kn}`.

        """
        # Check if all mesh boundaries have bcs presibed frm input
        self.__check_bcs()

        # Setup required function spaces
        self.__setup_function_spaces()

        # Get local variables
        mesh = self.mesh
        boundaries = self.boundaries
        bcs = self.bcs
        kn = df.Constant(self.kn)
        xi_tilde = df.Constant(self.xi_tilde)
        delta_1 = df.Constant(self.delta_1)
        delta_2 = df.Constant(self.delta_2)
        delta_3 = df.Constant(self.delta_3)

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
                + 12/5 * kn * df.inner(to.stf3d2(df.grad(s)), df.grad(r))
                + 2/3 * (1/kn) * df.inner(s, r)
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
                + 2 * kn * df.inner(
                    to.stf3d3(to.grad3dOf2(to.gen3dTF2(sigma))),
                    to.grad3dOf2(to.gen3dTF2(psi))
                )
                + (1/kn) * df.inner(to.gen3dTF2(sigma), to.gen3dTF2(psi))
                - 2 * df.dot(u, df.div(psi))
                + cpl * 4/5 * df.inner(to.stf3d2(df.grad(s)), psi)
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
                bcs[bc]["epsilon_w"] * (p + sigma_nn) * psi_nn * df.ds(bc)
                for bc in bcs.keys()
            ])
            l3 = sum([
                - 2.0 * psi_nt * bcs[bc]["u_t_w"] * df.ds(bc)
                - 2.0 * (
                    - bcs[bc]["epsilon_w"] * bcs[bc]["p_w"]
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
                bcs[bc]["epsilon_w"] * (p + sigma_nn) * q * df.ds(bc)
                for bc in bcs.keys()
            ])
            l5 = - (f_mass * q) * df.dx + sum([
                (
                    - bcs[bc]["epsilon_w"] * bcs[bc]["p_w"]
                    + bcs[bc]["u_n_w"]
                ) * q * df.ds(bc)
                for bc in bcs.keys()
            ])
        else:
            a1 = (
                kn * df.inner(to.stf3d2(df.grad(s)), df.grad(r))
                + (1/kn) * df.inner(s, r)
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
        Solve the previously assembled system.

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

        print("Start solving system..")
        start_t = time_module.time()
        sol = df.Function(w)
        df.solve(
            self.form_a == self.form_b, sol, [],
            solver_parameters={"linear_solver": "mumps"}
        )
        end_t = time_module.time()
        print("Finished solving system in: {}".format(str(end_t - start_t)))

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
                mean_p_value = self.__calc_sf_mean(p_i)
                mean_p_fct = df.Function(self.fspaces["p"])
                mean_p_fct.assign(df.Constant(mean_p_value))
                p_i.assign(p_i - mean_p_fct)
                self.sol["p"] = p_i

        # Calculate mass flows
        for bc_id in self.massflow:
            if bc_id not in self.boundaries.array():
                raise Exception("Massflow: {} is no boundary.".format(bc_id))
            n = df.FacetNormal(self.mesh)
            mass_flow_rate = df.assemble(
                df.inner(self.sol["u"], n)*df.ds(bc_id)
            )
            print("mass flow rate of BC", bc_id, ":", mass_flow_rate)
            self.write_content_to_file("massflow_"+str(bc_id), mass_flow_rate)

    def __load_exact_solution(self):
        """
        Load exact solution from the location given in ``input.yml``.

        The exact solution must be C++ format with a specific syntax.
        The ``esol.cpp`` must contain the classes:

        =============== =====================
        Class           mode
        =============== =====================
        ``Temperature`` ``heat`` or ``r13``
        ``Heatflux``    ``heat`` or ``r13``
        ``Pressure``    ``stress`` or ``r13``
        ``Velocity``    ``stress`` or ``r13``
        ``Stress``      ``stress`` or ``r13``
        =============== =====================

        The file has to follow a specific syntax for DOLFIN.
        An example file could look like:

        .. code-block:: c++

            #include <pybind11/pybind11.h>
            #include <pybind11/eigen.h>
            #include <cmath>
            // additional includes
            #include <boost/math/special_functions/bessel.hpp>
            using namespace std;
            namespace py = pybind11;
            #include <dolfin/function/Expression.h>
            double lambda_3 = sqrt(3.0/2.0); // some own constants
            class Temperature : public dolfin::Expression {
                public:
                Temperature() : dolfin::Expression() {}
                void eval(Eigen::Ref<Eigen::VectorXd> values,
                        Eigen::Ref<const Eigen::VectorXd> x) const override {
                    values[0] = 1; // value
                }
            };
            class Heatflux : public dolfin::Expression {
                public:
                Heatflux() : dolfin::Expression(2) {} // note components=2!
                void eval(Eigen::Ref<Eigen::VectorXd> values,
                        Eigen::Ref<const Eigen::VectorXd> x) const override {
                    values[0] = 42;
                    values[1] = 3.14;
                }
            };
            class Pressure : public dolfin::Expression {
                public:
                Pressure() : dolfin::Expression() {}
                void eval(Eigen::Ref<Eigen::VectorXd> values,
                        Eigen::Ref<const Eigen::VectorXd> x) const override {
                    values[0] = boost::math::cyl_bessel_i(1,2.71); // external
                }
            };
            class Velocity : public dolfin::Expression {
                public:
                Velocity() : dolfin::Expression(2) {}
                void eval(Eigen::Ref<Eigen::VectorXd> values,
                        Eigen::Ref<const Eigen::VectorXd> x) const override {
                    values[0] = lambda_3;
                    values[1] = 2;
                }
            };
            class Stress : public dolfin::Expression {
                public:
                Stress() : dolfin::Expression(2,2) {} // note dim=2, shape=(2,2)
                void eval(Eigen::Ref<Eigen::VectorXd> values,
                        Eigen::Ref<const Eigen::VectorXd> x) const override {
                    double xx_val = 1.23;
                    double xy_val = 1.23;
                    double yy_val = 1.23;
                    values[0] = xx_val;
                    values[1] = xy_val;
                    values[2] = yy_val;
                    // values[3] = xy_val // not used due to symmetry, skip
                }
            };
            PYBIND11_MODULE(SIGNATURE, m) { // needed for DOLFIN
                py::class_<Temperature, std::shared_ptr<Temperature>,
                           dolfin::Expression>
                    (m, "Temperature")
                .def(py::init<>());
                py::class_<Heatflux, std::shared_ptr<Heatflux>,
                           dolfin::Expression>
                    (m, "Heatflux")
                .def(py::init<>());
                py::class_<Pressure, std::shared_ptr<Pressure>,
                           dolfin::Expression>
                    (m, "Pressure")
                .def(py::init<>());
                py::class_<Velocity, std::shared_ptr<Velocity>,
                           dolfin::Expression>
                    (m, "Velocity")
                .def(py::init<>());
                py::class_<Stress, std::shared_ptr<Stress>, dolfin::Expression>
                    (m, "Stress")
                .def(py::init<>());
            }
        """
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

    def __calc_sf_mean(self, scalar_function):
        """
        Calculate the mean of a scalar function.

        .. code-block:: python

            np.set_printoptions(precision=16)
            # Precision is not soo nice, only 9 digits:
            print(mean)
            # In solve() has m. prec. hmmm:
            print(self.__calc_sf_mean(self.sol["p"]))

        .. note::

            The following does not work in parallel because the mean is
            then only local. So convergence studies have to be performed in
            serial:

            .. code-block:: python

                mean = np.mean(scalar_function.compute_vertex_values())
        """
        #
        v = scalar_function.compute_vertex_values()
        mean = np.mean(v)
        return mean

    def __calc_field_errors(self, field_, field_e_, v_field, name_):
        r"""
        Calculate both :math:`L_2` and :math:`l_\infty` errors.

        Works for scalars, vectors and tensors.
        The difference is written to a file.
        The exact solution is written to a file.
        Relative errors are based per field component and the maximum value of
        the analytical solution. If the analytical solution is uniformly zero,
        then the absolute erorrs is used.
        (equivalent to setting the maximum to 1)

        Parameters
        ----------
        field_ : DOLFIN function
            calculated field
        field_e_ : DOLFIN function
            exact solution of field
        v_field : DOLFIN fspace
            function space for error calculation
        name_ : string
            name of the field, used to write difference

        Returns
        -------
        dict
            Dict with an error list for "L_2" and a list for "l_inf"

        Raises
        ------
        Nothing

        See Also
        --------
        calculate_errors: Function to return all field errors

        Notes
        -----
        For other norm types, see DOLFIN documentation [1]_ and search for
        norms.

        References
        ----------
        .. [1] `DOLFIN documentation <https://fenicsproject.org/docs/dolfin/>`_

        Examples
        --------
        Here should be some doctest examples.

        >>> a=1
        >>> b=2
        >>> print(a+b)
        3

        """
        field_e_i = df.interpolate(field_e_, v_field)
        field_i = df.interpolate(field_, v_field)

        difference = df.project(field_e_i - field_i, v_field)
        self.__write_xdmf("difference_{}".format(name_), difference, False)

        dofs = len(field_e_i.split()) or 1

        if dofs == 1:
            # scalar
            errs_f_L2 = [df.errornorm(field_e_i, field_i, "L2")]
            errs_v_linf = [
                np.max(
                    np.abs(
                        field_e_i.compute_vertex_values()
                        - field_i.compute_vertex_values()
                    )
                )
            ]
        else:
            # vector or tensor
            errs_f_L2 = [df.errornorm(
                field_e_i.split()[i], field_i.split()[i], "L2"
            ) for i in range(dofs)]
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
            if dofs == 1:
                # scalar
                max_esols = [
                    np.max(np.abs(field_e_i.compute_vertex_values())) or 1
                ]
            else:
                # vector or tensor
                max_esols = [
                    np.max(
                        np.abs(field_e_i.split()[i].compute_vertex_values())
                    )
                    for i in range(dofs)
                ]
            errs_f_L2 = [x/y for x, y in zip(errs_f_L2, max_esols)]
            errs_v_linf = [x/y for x, y in zip(errs_v_linf, max_esols)]

        print("Error " + str(name_) + " L_2:", errs_f_L2)
        print("Error " + str(name_) + " l_inf:", errs_v_linf)

        self.__write_xdmf(name_ + "_e", field_e_i, False)

        return [{
            "L_2": errs_f_L2[i],
            "l_inf": errs_v_linf[i],
        } for i in range(dofs)]

    def calculate_errors(self):
        """
        Calculate and return the errors of numerical to exact solution.

        This includes all calculated fields.

        .. note::

            Usage of `np.max()` does not work in parallel.
            So convergence studies have to be performed in serial for now.
            Final fields should be right, so MPI can be used for production
            simulations.

        Returns:
            dict -- Errors

        """
        print("Calculate errors..")

        self.__load_exact_solution()

        if self.mode == "heat" or self.mode == "r13":
            se = self.__calc_field_errors(
                self.sol["theta"], self.esol["theta"],
                self.fspaces["theta"], "theta"
            )
            ve = self.__calc_field_errors(
                self.sol["s"], self.esol["s"],
                self.fspaces["s"], "s"
            )
            ers = self.errors
            ers["theta"] = se[0]
            ers["sx"] = ve[0]
            ers["sy"] = ve[1]
        if self.mode == "stress" or self.mode == "r13":
            se = self.__calc_field_errors(
                self.sol["p"], self.esol["p"],
                self.fspaces["p"], "p"
            )
            ve = self.__calc_field_errors(
                self.sol["u"], self.esol["u"],
                self.fspaces["u"], "u"
            )
            te = self.__calc_field_errors(
                self.sol["sigma"], self.esol["sigma"],
                self.fspaces["sigma"], "sigma"
            )
            ers = self.errors
            ers["p"] = se[0]
            ers["ux"] = ve[0]
            ers["uy"] = ve[1]
            ers["sigmaxx"] = te[0]
            ers["sigmaxy"] = te[1]
            ers["sigmayy"] = te[2]

        return self.errors

    def write_content_to_file(self, filename, content):
        """Write content to a file in the output folder."""
        path = self.output_folder + "/" + filename
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, mode='w') as file:
            print("Write: {}".format(path))
            file.write(str(content))

    def write(self):
        """
        Write all solver data to separate folder.

        This includes the writing of:

        (#) Solutions
        (#) Parameter functions
        (#) System matrices if set in input file
        """
        print("Write fields..")

        self.__write_solutions()
        self.__write_parameters()
        if self.write_systemmatrix:
            self.__write_discrete_system()

    def __write_solutions(self):
        """Write all solutions fields."""
        sols = self.sol
        for field in sols:
            if sols[field] is not None:
                self.__write_xdmf(field, sols[field], self.write_pdfs)

    def __write_parameters(self):
        """
        Write parameter functions for debug reasons.

        This includes:

        (#) Heat source as `f_mass`
        (#) Mass Source as `f_heat`

        The parameter functions are internally interpolated into a :math:`P_1`
        space before writing.
        """
        # Interpolation setup
        el_str = "Lagrange"
        deg = 1
        el = df.FiniteElement(el_str, degree=deg, cell=self.cell)
        V = df.FunctionSpace(self.mesh, el)

        # Heat source
        f_heat = df.interpolate(self.heat_source, V)
        self.__write_xdmf("f_heat", f_heat, False)

        # Mass source
        f_mass = df.interpolate(self.mass_source, V)
        self.__write_xdmf("f_mass", f_mass, False)

    def __write_discrete_system(self):
        r"""
        Write the discrete system matrix and the RHS vector.

        Can be used to analyze e.g. condition number.
        Include writing of :math:`\mathbf{A}` and :math:`\mathbf{b}`
        of :math:`\mathbf{A} \mathbf{x} = \mathbf{b}`.

        File-ending is `.mat`.

        Import the matrices/vectors e.g. into Matlab with:

        .. code-block:: matlab

            % Input into MATLAB
            At = readtable("A.mat");
            bt = readtable("b.mat");
            A = table2array(At);
            b = table2array(bt);

        Example
        -------
        >>> # Construct LHS
        >>> from dolfin import *
        >>> mesh = IntervalMesh(2 ,0, 1)
        >>> V = FunctionSpace(mesh, "Lagrange", 1)
        >>> u = TrialFunction(V)
        >>> v = TestFunction(V)
        >>> a = inner(grad(u),grad(v))*dx
        >>> L = df.Constant(1)*v*dx
        >>> lhs = assemble(a)
        >>> print(lhs.array())
        [[ 2. -2.  0.]
         [-2.  4. -2.]
         [ 0. -2.  2.]]
        >>> rhs = assemble(L)
        >>> print(rhs.get_local())
        [ 0.25  0.5   0.25]
        >>> # Assign LHS to solver
        >>> from input import Input
        >>> from meshes import H5Mesh
        >>> params = Input(
        ...     "tests/heat/inputs/heat_01_coeffs_p1p1_stab.yml"
        ... ) # doctest: +ELLIPSIS
        Input:...
        >>> msh = H5Mesh("tests/mesh/ring0.h5")
        >>> solver = Solver(params.dict, msh, "0") # "0" means time=0
        >>> solver.form_a = a
        >>> solver.form_b = L
        >>> solver.output_folder = "./"
        >>> solver._Solver__write_discrete_system()
        >>> print(open("A_0.mat","r").read())
        2.0000...00000e+00 -2.0000...00000e+00 0.000...000000e+00
        -2.0000...00000e+00 4.0000...00000e+00 -2.0000...00000e+00
        0.0000...00000e+00 -2.0000...00000e+00 2.000...000000e+00
        <BLANKLINE>
        >>> print(open("b_0.mat","r").read())
        2.500000000000000000e-01
        5.000000000000000000e-01
        2.500000000000000000e-01
        <BLANKLINE>
        """
        file_ending = ".mat"
        np.savetxt(
            self.output_folder + "A_{}".format(self.time) + file_ending,
            df.assemble(self.form_a).array()
        )
        np.savetxt(
            self.output_folder + "b_{}".format(self.time) + file_ending,
            df.assemble(self.form_b)
        )

    def __write_xdmf(self, name, field, write_pdf):
        """
        Write a given field to a XDMF file in the output folder.

        *Arguments*
            name
                The name to be given to the field. Will be used as filename
                and is the name of the field in e.g. Paraview.
            field
                The field to write.
            write_pdf
                If true, write a simple PDF plot for all solution fields
        """
        fname_xdmf = (
            self.output_folder + name + "_" + str(self.time) + ".xdmf"
        )
        with df.XDMFFile(self.mesh.mpi_comm(), fname_xdmf) as file:
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

            print("Write {}".format(fname_xdmf))
            file.write(field, self.time)

        if write_pdf:
            import matplotlib.pyplot as plt
            plt.switch_backend('agg')
            dimension = len(field.value_shape())

            if dimension < 2:
                # skip tensors
                fname_pdf = (
                    self.output_folder + name + "_" + str(self.time) + ".pdf"
                )
                plot = df.plot(field)
                plt.colorbar(plot)
                plt.xlabel("x")
                plt.ylabel("y")
                plt.title(field)
                print("Write {}".format(fname_pdf))
                plt.savefig(fname_pdf, dpi=150)
                plt.close()

            if dimension > 0:
                # skip scalars
                components = len(field.split())
                indexMap = {
                    1: {
                        1: "x",
                        2: "y"
                    },
                    2: {
                        1: "xx",
                        2: "xy",
                        3: "yx",
                        4: "yy",
                    }
                }
                for i in range(components):
                    fieldname = name + "_" + str(indexMap[dimension][i+1])
                    fname_pdf = (
                        self.output_folder + fieldname
                        + "_" + str(self.time) + ".pdf"
                    )
                    plot = df.plot(field.split()[i])
                    plt.colorbar(plot)
                    plt.xlabel("x")
                    plt.ylabel("y")
                    plt.title(fieldname)
                    print("Write {}".format(fname_pdf))
                    plt.savefig(fname_pdf, dpi=150)
                    plt.close()
