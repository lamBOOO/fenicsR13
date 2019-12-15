# pylint: disable=invalid-name,too-many-lines

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
import fenicsR13.tensoroperations as to


class Solver:
    r"""
    Class to store the actual solver.

    Possible order of methods in context of convergence study
    (see main program):
    "mesh=meshes[i=0]",
    "__init__", "assemble()",
    "solve()",
    "write()",
    "...",
    "mesh=meshes[i+1]",
    "__init__",
    "..."

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
    >>> from fenicsR13.input import Input
    >>> from fenicsR13.meshes import H5Mesh
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
        self.form_lhs = None
        self.form_rhs = None
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
        Assemble the weak form of the system, depending on the mode. This
        routine contains the main implementation of a stable variational
        formulation of the R13 system which allows equal order ansatz
        spaces for all variables and correctly incorporates the R13
        boundary conditions.

        The system of partial differential equations results from the two
        dimensional, linearized R13 equations of [1]_, that have been
        derived from Boltzmann's equation as a continuum model for
        rarefied gas flows. The boundary conditions have been presented
        in [2]_, while the variational formulation and its stabilization was
        introduced in [3]_. See also [4]_ for a general review and [5]_ for
        the modified stable boundary conditions.

        .. [1] H. Struchtrup, M. Torrilhon (2003). Regularization of
            Grad's 13 Moment Equations: Derivation and Linear Analysis.
        .. [2] M. Torrilhon, H. Struchtrup (2007). Boundary Conditions for
            Regularized 13-Moment-Equations for Micro-Channel-Flows
        .. [3] A. Westerkamp, M. Torrilhon (2019). Finite Element Methods
            for the Linear Regularized 13-Moment Equations Describing Slow
            Rarefied Gas Flows
        .. [4] M. Torrilhon (2016). Modeling Nonequilibrium Gas Flow Based
           on Moment Equations
        .. [5] M. Torrilhon, N. Sarna (2017). Hierarchical Boltzmann
           Simulations and Model Error Estimation

        Identities:

        - :math:`\boldsymbol{x}_1 \in \mathbb{R}^3`,
          :math:`\boldsymbol{x}_2 \in \mathbb{R}^{3 \times 3}`,
          ...
        - Test function :math:`\boldsymbol{\psi}` is assumed to be
          symmetric and trace-less:

        .. math::
            \boldsymbol{\psi} : \boldsymbol{I} = 0 \\
            (\boldsymbol{\psi})_{\text{sym}} = \boldsymbol{\psi}

        - Trace of vector gradient is divergence of vector:

        .. math::
            \langle
            \boldsymbol{I},\boldsymbol\nabla \boldsymbol{x}_1
            \rangle
            =
            \boldsymbol{I} : \boldsymbol\nabla \boldsymbol{x}_1
            =
            \text{div}(\boldsymbol{x}_1)

        - Inner product has orthogonality property with respect to the
          additive symmetric/skewsymmetric tensor decomposition:

        .. math::
            \langle
            (\boldsymbol{x}_2)_{\text{sym}}
            ,
            \boldsymbol{y}_2
            \rangle
            &=
            \langle
            (\boldsymbol{x}_2)_{\text{sym}}
            ,
            (\boldsymbol{y}_2)_{\text{sym}}
            +
            (\boldsymbol{y}_2)_{\text{skew}}
            \rangle
            \\
            &=
            \langle
            (\boldsymbol{x}_2)_{\text{sym}}
            ,
            (\boldsymbol{y}_2)_{\text{sym}}
            \rangle
            +
            \langle
            (\boldsymbol{x}_2)_{\text{sym}}
            ,
            (\boldsymbol{y}_2)_{\text{skew}}
            \rangle
            \\
            &=
            \langle
            (\boldsymbol{x}_2)_{\text{sym}}
            ,
            (\boldsymbol{y}_2)_{\text{sym}}
            \rangle

        - Inner product of STF tensor with arbitrary tensor:

        .. math::

            T_{i_{1} i_{2} \ldots i_{n} k_{1} \cdots k_{r}}
            A_{i_{1} i_{2} \ldots i_{n} j_{1} \cdots j_{n}}
            =
            T_{i_{1} i_{2} \ldots i_{n} k_{1} \cdots k_{r}}
            A_{\langle i_{1} i_{2} \ldots i_{n}\rangle j_{1} \cdots j_{n}}

        Tricks of the trade:

        .. math::
            \langle
            {(\boldsymbol\nabla \boldsymbol{s})}_{\text{STF}}
            ,
            \boldsymbol\nabla \boldsymbol{r}
            \rangle
            &=
            \langle
            {(\boldsymbol\nabla \boldsymbol{s})}_{\text{sym}}
            -
            \frac13\text{tr}(\boldsymbol\nabla \boldsymbol{s})\boldsymbol{I}
            ,
            \boldsymbol\nabla \boldsymbol{r}
            \rangle
            \\
            &=
            \langle
            {(\boldsymbol\nabla \boldsymbol{s})}_{\text{sym}}
            ,
            \boldsymbol\nabla \boldsymbol{r}
            \rangle
            -\frac13\text{tr}(\boldsymbol\nabla \boldsymbol{s})
            \langle
            \boldsymbol{I},(\boldsymbol\nabla \boldsymbol{r})
            \rangle
            \\
            &=
            \langle
            {(\boldsymbol\nabla \boldsymbol{s})}_{\text{sym}}
            ,
            {(\boldsymbol\nabla \boldsymbol{r})}_{\text{sym}}
            \rangle
            -\frac13\text{tr}(\boldsymbol\nabla \boldsymbol{s})
            \text{tr}(\boldsymbol\nabla \boldsymbol{r})
            \\
            &=
            \langle
            {(\boldsymbol\nabla \boldsymbol{s})}_{\text{sym}}
            ,
            {(\boldsymbol\nabla \boldsymbol{r})}_{\text{sym}}
            \rangle
            -\frac13\text{div}(\boldsymbol{s})
            \text{div}(\boldsymbol{r})

        .. math::
            \langle
            {(\boldsymbol\nabla \boldsymbol{s})}_{\text{STF}}
            ,
            \boldsymbol{\psi}
            \rangle
            &=
            \langle
            {(\boldsymbol\nabla \boldsymbol{s})}_{\text{sym}}
            -
            \frac13\text{tr}(\boldsymbol\nabla \boldsymbol{s})\boldsymbol{I}
            ,
            \boldsymbol{\psi}
            \rangle
            \\
            &=
            \langle
            {(\boldsymbol\nabla \boldsymbol{s})}_{\text{sym}}
            ,
            \boldsymbol{\psi}
            \rangle
            -\frac13\text{tr}(\boldsymbol\nabla \boldsymbol{s})
            \langle
            \boldsymbol{I},\boldsymbol{\psi}
            \rangle
            \\
            &=
            \langle
            {\boldsymbol\nabla \boldsymbol{s}}
            ,
            {\boldsymbol{\psi}}
            \rangle
            -\frac13\text{div}(\boldsymbol{s})
            \text{tr}(\boldsymbol{\psi})
            \\
            &=
            \langle
            {\boldsymbol\nabla \boldsymbol{s}}
            ,
            {\boldsymbol{\psi}}
            \rangle

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

        # Define custom measeasure for boundaries
        df.ds = df.Measure("ds", domain=mesh, subdomain_data=boundaries)
        df.dS = df.Measure("dS", domain=mesh, subdomain_data=boundaries)

        # Define mesh measuers
        h_msh = df.CellDiameter(mesh)
        h_avg = (h_msh("+") + h_msh("-"))/2.0 # pylint: disable=not-callable
        # TODO: Study this, is it more precise?
        # fa = df.FacetArea(mesh)
        # h_avg_new = (fa("+") + fa("-"))/2.0 # pylint: disable=not-callable

        # Setup trial and test functions
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

        # Setup source functions
        f_heat = self.heat_source
        f_mass = self.mass_source

        # Decouple heat/stress switch
        if self.mode == "r13":
            cpl = 1
        else:
            cpl = 0

        # Stabilization switch
        if self.use_cip:
            cip = 1
        else:
            cip = 0

        # Setup normal/tangential projections
        # => tangential (tx,ty) = (-ny,nx) = perp(n) only for 2D
        n_vec = df.FacetNormal(mesh)
        t_vec = ufl.perp(n_vec)
        def n(rank1):
            return df.dot(rank1, n_vec)
        def t(rank1):
            return df.dot(rank1, t_vec)
        def nn(rank2):
            return df.dot(rank2 * n_vec, n_vec)
        def tt(rank2):
            return df.dot(rank2 * t_vec, t_vec)
        def nt(rank2):
            return df.dot(rank2 * n_vec, t_vec)

        # Sub functionals:
        # 1) Diagonals:
        def a(s_, r_):
            return (
                # => 24/25*stf(grad)*grad
                + 24/25 * kn * df.inner(
                    df.sym(df.grad(s_)), df.sym(df.grad(r_))
                )
                - 24/75 * kn * df.div(s_) * df.div(r_)
                + 4/15 * (1/kn) * df.inner(s_, r_)
            ) * df.dx + (
                + 1/(2*xi_tilde) * n(s_) * n(r_)
                + 12/25 * xi_tilde * t(s_) * t(r_)
                - (1-cpl) * 1/25 * xi_tilde * t(s_) * t(r_)
            ) * df.ds
        def d(sigma_, psi_):
            return (
                kn * df.inner(
                    to.stf3d3(to.grad3dOf2(to.gen3dTF2(sigma_))),
                    to.stf3d3(to.grad3dOf2(to.gen3dTF2(psi_)))
                )
                + (1/(2*kn)) * df.inner(
                    to.gen3dTF2(sigma_), to.gen3dTF2(psi_)
                )
            ) * df.dx + (
                + xi_tilde * 9/8 * nn(sigma_) * nn(psi_)
                - xi_tilde * (1-cpl) * 3/40 * nn(sigma_) * nn(psi_)
                + xi_tilde * (
                    (tt(sigma_) + (1/2) * nn(sigma_)) *
                    (tt(psi_) + (1/2) * nn(psi_))
                )
                + (1/xi_tilde) * nt(sigma_) * nt(psi_)
            ) * df.ds + sum([
                bcs[bc]["epsilon_w"] * nn(sigma_) * nn(psi_) * df.ds(bc)
                for bc in bcs.keys()
            ])
        def h(p, q):
            return sum([
                bcs[bc]["epsilon_w"] * p * q * df.ds(bc)
                for bc in bcs.keys()
            ])
        # 2) Offdiagonals:
        def b(scalar, vector):
            return 1 * scalar * df.div(vector) * df.dx
        def c(vector, tensor):
            return cpl * ((
                2/5 * df.inner(tensor, df.grad(vector))
            ) * df.dx + (
                - 3/20 * nn(tensor) * n(vector)
                - 1/5 * nt(tensor) * t(vector)
            ) * df.ds)
        def e(vector, tensor):
            return 1 * df.dot(df.div(tensor), vector) * df.dx
        def f(scalar, tensor):
            return sum([
                bcs[bc]["epsilon_w"] * scalar * nn(tensor) * df.ds(bc)
                for bc in bcs.keys()
            ])
        def g(scalar, vector):
            return 1 * df.inner(vector, df.grad(scalar)) * df.dx
        # 3) CIP Stabilization:
        def j_theta():
            return (
                + delta_1 * h_avg**3 *
                df.jump(df.grad(theta), n_vec) * df.jump(df.grad(kappa), n_vec)
            ) * df.dS
        def j_u():
            return (
                + delta_2 * h_avg**3 *
                df.dot(df.jump(df.grad(u), n_vec), df.jump(df.grad(v), n_vec))
            ) * df.dS
        def j_p():
            return (
                + delta_3 * h_avg *
                df.jump(df.grad(p), n_vec) * df.jump(df.grad(q), n_vec)
            ) * df.dS

        # Setup all equations
        lhs = [None] * 5
        rhs = [None] * 5
        # 1) Left-hand sides
        lhs[0] = +1*a(s, r)    -b(theta, r)-c(r, sigma)  +0        +0
        lhs[1] = +1*b(kappa, s)+0          +0            +0        +0
        lhs[2] = +1*c(s, psi)  +0          +d(sigma, psi)-e(u, psi)+f(p, psi)
        lhs[3] = +1*0          +0          +e(v, sigma)  +0        +g(p, v)
        lhs[4] = +1*0          +0          +f(q, sigma)  -g(q, u)  +h(p, q)
        # 2) Right-hand sides:
        rhs[0] = sum([
            - 1 * n(r) * bcs[bc]["theta_w"] * df.ds(bc)
            for bc in bcs.keys()
        ])
        rhs[1] = f_heat * kappa * df.dx
        rhs[2] = sum([
            - 1 * nt(psi) * bcs[bc]["u_t_w"] * df.ds(bc)
            - 1 * (
                - bcs[bc]["epsilon_w"] * bcs[bc]["p_w"]
                + bcs[bc]["u_n_w"]
            ) * nn(psi) * df.ds(bc)
            for bc in bcs.keys()
        ])
        rhs[3] = + df.Constant(0) * df.div(v) * df.dx
        rhs[4] = + (f_mass * q) * df.dx - sum([
            (
                - bcs[bc]["epsilon_w"] * bcs[bc]["p_w"]
                + bcs[bc]["u_n_w"]
            ) * q * df.ds(bc)
            for bc in bcs.keys()
        ])

        # Combine all equations to compound weak form and add CIP
        if self.mode == "heat":
            self.form_lhs = sum(lhs[0:2]) + cip * j_theta()
            self.form_rhs = sum(rhs[0:2])
        elif self.mode == "stress":
            self.form_lhs = sum(lhs[2:5]) + cip * (j_u() + j_p())
            self.form_rhs = sum(rhs[2:5])
        elif self.mode == "r13":
            self.form_lhs = sum(lhs) + cip * (j_theta() + j_u() + j_p())
            self.form_rhs = sum(rhs)

    def solve(self):
        """
        Solve the previously assembled system.

        Some available solver options:

        .. code-block:: python

            # Some solver params
            solver_parameters={
                "linear_solver": "gmres", "preconditioner": "ilu" # or
                "linear_solver": "petsc", "preconditioner": "ilu" # or
                "linear_solver": "direct" # or
                "linear_solver": "mumps" # or
                "linear_solver": "mumps"
            }
            # List all available solvers:
            list_linear_solver_methods()
            list_krylov_solver_preconditioners()
            # "direct" means "default" means "lu" of default backend
            print(parameters["linear_algebra_backend"]) # usually PETSc

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
            self.form_lhs == self.form_rhs, sol, [],
            solver_parameters={"linear_solver": "umfpack"}
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
        For other norm types, see DOLFIN documentation [6]_ and search for
        norms.

        References
        ----------
        .. [6] `DOLFIN documentation <https://fenicsproject.org/docs/dolfin/>`_

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
            A = table2array(readtable("A_0.mat","FileType","text"));
            b = table2array(readtable("b_0.mat","FileType","text"));

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
        >>> from fenicsR13.input import Input
        >>> from fenicsR13.meshes import H5Mesh
        >>> params = Input(
        ...     "tests/heat/inputs/heat_01_coeffs_p1p1_stab.yml"
        ... ) # doctest: +ELLIPSIS
        Input:...
        >>> msh = H5Mesh("tests/mesh/ring0.h5")
        >>> solver = Solver(params.dict, msh, "0") # "0" means time=0
        >>> solver.form_lhs = a
        >>> solver.form_rhs = L
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
            df.assemble(self.form_lhs).array()
        )
        np.savetxt(
            self.output_folder + "b_{}".format(self.time) + file_ending,
            df.assemble(self.form_rhs)
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
            import matplotlib.pyplot as plt # pylint: disable=C0415
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