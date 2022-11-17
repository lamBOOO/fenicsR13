# pylint: disable=invalid-name,too-many-lines
# pylint: disable=not-callable
# noqa: E226

"""
Solver module, contains the Solver class.

For usage examples, see the :class:`solver.Solver` description.
"""

import sys
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
        self.params = params  #: Doctest
        self.mesh = mesh.mesh
        self.regions = mesh.subdomains
        self.boundaries = mesh.boundaries
        self.cell = self.mesh.ufl_cell()
        self.time = time

        self.comm = df.MPI.comm_world
        self.rank = df.MPI.rank(self.comm)

        # CIP
        self.delta_theta = self.params["stabilization"]["cip"]["delta_theta"]
        self.delta_u = self.params["stabilization"]["cip"]["delta_u"]
        self.delta_p = self.params["stabilization"]["cip"]["delta_p"]

        # Create region field expressions
        self.regs = copy.deepcopy(self.params["regs"])
        for reg_id in self.regs:
            for field in self.regs[reg_id].keys():
                self.regs[reg_id][field] = df.Expression(str(self.regs[reg_id][field]),degree=2)

        # Create boundary field expressions
        self.bcs = copy.deepcopy(self.params["bcs"])
        for edge_id in self.bcs:
            for field in self.bcs[edge_id].keys():
                self.bcs[edge_id][field] = df.Expression(str(self.bcs[edge_id][field]),degree=2)

        self.heat_source = df.Expression(str(self.params["heat_source"]),degree=2)
        self.mass_source = df.Expression(str(self.params["mass_source"]),degree=2)
        self.body_force = df.Expression([str(i) for i in self.params["body_force"]],degree=2)

        self.output_folder = self.params["output_folder"] + "/"
        self.elems = {
            "theta": None,
            "s": None,
            "p": None,
            "u": None,
            "sigma": None,
        }
        ## fspaces[var] seems only used when projecting pressure to mean-zero. 
        self.fspaces = {
            "theta": None,
            "s": None,
            "p": None,
            "u": None,
            "sigma": None,
        }
        self.r13_elems = None
        self.r13_fspaces = None
        self.form_lhs = None
        self.form_rhs = None
        self.sol = {
            "theta": None,
            "s": None,
            "p": None,
            "u": None,
            "sigma": None,
        }

    def __setup_function_spaces(self):
        """
        Set up function spaces for trial and test functions for assembling.

        """
        # Setup elements for all fields
        cell = self.cell
        msh = self.mesh

        e = "Lagrange"
        deg = 1
        # scalar variables
        self.elems["theta"] = df.FiniteElement(e, cell, deg)
        self.elems["p"] = df.FiniteElement(e, cell, deg)
        # vector variables
        self.elems["s"] = df.VectorElement(e, cell, deg)
        self.elems["u"] = df.VectorElement(e, cell, deg)
        # 2-tensor variables
        self.elems["sigma"] = df.TensorElement(e, cell, deg, symmetry={(0, 1): (1, 0)})

        ## fspaces[var] seems only used when projecting pressure to mean-zero. 
        for var in self.elems:
            self.fspaces[var] = df.FunctionSpace(msh, self.elems[var])

        # Bundle elements per mode into `mxd_elems` dict
        r13_elems = [self.elems["theta"], self.elems["s"], self.elems["p"], self.elems["u"], self.elems["sigma"]]
        self.r13_elems = df.MixedElement(r13_elems)
        self.r13_fspaces = df.FunctionSpace(msh, self.r13_elems)

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
        # Setup required function spaces
        self.__setup_function_spaces()

        # Get local variables
        mesh = self.mesh
        regions = self.regions
        regs = self.regs
        boundaries = self.boundaries
        bcs = self.bcs
        delta_theta = df.Constant(self.delta_theta)
        delta_u = df.Constant(self.delta_u)
        delta_p = df.Constant(self.delta_p)

        # Define custom measeasures for boundary edges and inner edges
        df.dx = df.Measure("dx", domain=mesh, subdomain_data=regions)
        df.ds = df.Measure("ds", domain=mesh, subdomain_data=boundaries)
        df.dS = df.Measure("dS", domain=mesh, subdomain_data=boundaries)

        # Define mesh measuers
        h_msh = df.CellDiameter(mesh)
        h_avg = (h_msh("+") + h_msh("-")) / 2.0
    
        # TODO: Study this, is it more precise?
        # fa = df.FacetArea(mesh)
        # h_avg_new = (fa("+") + fa("-"))/2.0

        # Setup trial and test functions
        w_r13 = self.r13_fspaces
        (theta, s, p, u, sigma) = df.TrialFunctions(w_r13)
        (kappa, r, q, v, psi) = df.TestFunctions(w_r13)

        # Setup source functions
        f_heat = self.heat_source
        f_mass = self.mass_source
        f_body = self.body_force

        # Stabilization
        cip = 1

        # Setup normal/tangential projections
        # => tangential (tx,ty) = (-ny,nx) = perp(n) only for 2D
        n_vec = df.FacetNormal(mesh)
        t_vec = ufl.perp(n_vec)

        def n(rank1): return df.dot(rank1, n_vec)

        def t(rank1): return df.dot(rank1, t_vec)

        def nn(rank2): return df.dot(rank2 * n_vec, n_vec)

        def tt(rank2): return df.dot(rank2 * t_vec, t_vec)

        def nt(rank2): return df.dot(rank2 * n_vec, t_vec)

        # Sub functionals:
        # 1) Diagonals:
        def a(s, r):
            return sum([(
                + 24 / 25 * regs[reg]["kn"] * df.inner(df.sym(df.grad(s)), df.sym(df.grad(r)))
                - 24 / 75 * regs[reg]["kn"] * df.div(s) * df.div(r)
                + 4 / 5 * regs[reg]["kn"] * df.div(s) * df.div(r)
                + 4 / 15 * (1 / regs[reg]["kn"]) * df.inner(s, r)
            ) * df.dx(reg) for reg in regs.keys()]) + sum([(
                + 1 / (2 * bcs[bc]["chi_tilde"]) * n(s) * n(r)
                + 11 / 25 * bcs[bc]["chi_tilde"] * t(s) * t(r)
                + 1 / 25 * bcs[bc]["chi_tilde"] * t(s) * t(r)
            ) * df.ds(bc) for bc in bcs.keys()])

        def d(si, ps):
           return sum([(
                + regs[reg]["kn"] * df.inner(to.stf3d3(to.grad3dOf2(to.gen3dTF2(si))),to.stf3d3(to.grad3dOf2(to.gen3dTF2(ps))))
                + (1 / (2 * regs[reg]["kn"])) * df.inner(to.gen3dTF2(si), to.gen3dTF2(ps))
            ) * df.dx(reg) for reg in regs.keys()]) + sum([(
                + bcs[bc]["chi_tilde"] * 21 / 20 * nn(si) * nn(ps)
                + bcs[bc]["chi_tilde"] * 3 / 40 * nn(si) * nn(ps)
                + bcs[bc]["chi_tilde"] * ( (tt(si) + (1 / 2) * nn(si)) * (tt(ps) + (1 / 2) * nn(ps)) )
                + (1 / bcs[bc]["chi_tilde"]) * nt(si) * nt(ps)
                + bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * nn(si) * nn(ps)
            ) * df.ds(bc) for bc in bcs.keys()])

        def h(p, q):
            return sum([(
                bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * p * q
            ) * df.ds(bc) for bc in bcs.keys()])

        # 2) Offdiagonals:
        def b(th, r):
            return sum([(
                th * df.div(r)
            ) * df.dx(reg) for reg in regs.keys()])

        def c(r, si):
            return sum([(
                2 / 5 * df.inner(si, df.grad(r))
            ) * df.dx(reg) for reg in regs.keys()]) - sum([(
                3 / 20 * nn(si) * n(r) + 1 / 5 * nt(si) * t(r)
            ) * df.ds(bc) for bc in bcs.keys()])

        def e(u, ps):
            return sum([(
                df.dot(df.div(ps), u)
            ) * df.dx(reg) for reg in regs.keys()])

        def f(p, ps):
            return sum([(
                bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * p * nn(ps)
            ) * df.ds(bc) for bc in bcs.keys()])

        def g(p, v):
            return sum([(
                df.inner(v, df.grad(p))
            ) * df.dx(reg) for reg in regs.keys()])

        # 3.1) CIP Stabilization:
        def j_theta(theta, kappa):
            return ( delta_theta * h_avg**3 * df.jump(df.grad(theta), n_vec) * df.jump(df.grad(kappa), n_vec) ) * df.dS

        def j_u(u, v):
            return ( delta_u * h_avg**3 * df.dot(df.jump(df.grad(u), n_vec), df.jump(df.grad(v), n_vec)) ) * df.dS

        def j_p(p, q):
            return ( delta_p * h_avg * df.jump(df.grad(p), n_vec) * df.jump(df.grad(q), n_vec) ) * df.dS

        # Setup all equations
        A = [None] * 5
        L = [None] * 5
        # 1) Left-hand sides, bilinear form A[..]:
        # Changed inflow condition => minus before f(q, sigma)
        A[0] = a(s, r)     - b(theta, r) - c(r, sigma)   + 0         + 0
        A[1] = b(kappa, s) + 0           + 0             + 0         + 0
        A[2] = c(s, psi)   + 0           + d(sigma, psi) - e(u, psi) + f(p, psi)
        A[3] = 0           + 0           + e(v, sigma)   + 0         + g(p, v)
        A[4] = 0           + 0           + f(q, sigma)   - g(q, u)   + h(p, q)

        # 2) Right-hand sides, linear functional L[..]:
        L[0] = - sum([( bcs[bc]["theta_w"] * n(r) ) * df.ds(bc) for bc in bcs.keys()])
        # Use div(u)=f_mass to remain sym. (density-form doesnt need this):
        L[1] = (f_heat - f_mass) * kappa * df.dx
        L[2] = - sum([(
            + bcs[bc]["u_t_w"] * nt(psi)
            + ( bcs[bc]["u_n_w"] - bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * bcs[bc]["p_w"] ) * nn(psi)
        ) * df.ds(bc) for bc in bcs.keys()])
        L[3] = df.dot(f_body, v) * df.dx
        L[4] = (f_mass * q) * df.dx - sum([(
            ( bcs[bc]["u_n_w"] - bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * bcs[bc]["p_w"] ) * q
        ) * df.ds(bc) for bc in bcs.keys()])

        # Combine all equations to compound weak form and add stabilization
        self.form_lhs = sum(A) + cip * (j_theta(theta, kappa) + j_u(u, v) + j_p(p, q))
        self.form_rhs = sum(L)

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

            # List all available solvers and preconditioners:

            list_linear_solver_methods()
            bicgstab  | Biconjugate gradient stabilized method
            cg        | Conjugate gradient method
            default   | default linear solver
            gmres     | Generalized minimal residual method
            minres    | Minimal residual method
            mumps     | MUMPS (MUltifrontal Massively Parallel Sparse Direct)
            petsc     | PETSc built in LU solver
            richardson| Richardson method
            superlu   | SuperLU
            tfqmr     | Transpose-free quasi-minimal residual method
            umfpack   | UMFPACK (Unsymmetric MultiFrontal sparse LU factoriz.)
            # "direct" means "default" means "lu" of default backend
            print(parameters["linear_algebra_backend"]) # usually PETSc

            list_krylov_solver_preconditioners()
            amg              |  Algebraic multigrid
            default          |  default preconditioner
            hypre_amg        |  Hypre algebraic multigrid (BoomerAMG)
            hypre_euclid     |  Hypre parallel incomplete LU factorization
            hypre_parasails  |  Hypre parallel sparse approximate inverse
            icc              |  Incomplete Cholesky factorization
            ilu              |  Incomplete LU factorization
            jacobi           |  Jacobi iteration
            none             |  No preconditioner
            petsc_amg        |  PETSc algebraic multigrid
            sor              |  Successive over-relaxation

        """

        w = self.r13_fspaces

        print("Start assemble")
        sys.stdout.flush()
        start_t = time_module.time()
        AA = df.assemble(self.form_lhs)
        LL = df.assemble(self.form_rhs)
        end_t = time_module.time()
        secs = end_t - start_t
        print("Finish assemble: {}".format(str(secs)))
        sys.stdout.flush()

        print("Start solve")
        sys.stdout.flush()
        start_t = time_module.time()
        sol = df.Function(w)
        df.solve( AA, sol.vector(), LL, "mumps", "none" )
        end_t = time_module.time()
        secs = end_t - start_t
        print("Finished solve: {}".format(str(secs)))
        sys.stdout.flush()

        ( self.sol["theta"],self.sol["s"],self.sol["p"],self.sol["u"],self.sol["sigma"] ) = sol.split()

        # Scale pressure to have zero mean
        p_i = df.interpolate(self.sol["p"], self.fspaces["p"])
        v = p_i.compute_vertex_values()
        mean_p_value = np.mean(v)
        mean_p_fct = df.Function(self.fspaces["p"])
        mean_p_fct.assign(df.Constant(mean_p_value))
        p_i.assign(p_i - mean_p_fct)
        self.sol["p"] = p_i

        # compute average velocity
        vol = df.assemble(df.Constant(1) * df.dx)
        avgvel = df.assemble( abs(df.inner(self.sol["u"], self.sol["u"])) * df.dx ) / vol
        print("avg vel:", avgvel)

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

    def __write_solutions(self):
        """Write all solutions fields."""
        sols = self.sol
        for field in sols:
            if sols[field] is not None:
                self.__write_xdmf(field, sols[field] )

    def __write_parameters(self):
        """
        Write parameter functions for debug reasons.

        This includes:

        (#) Heat source as `f_mass`
        (#) Mass Source as `f_heat`
        (#) Body force as `f_body`

        The parameter functions are internally interpolated into a :math:`P_1`
        space before writing.
        """
        # Interpolation setup
        el_str = "Lagrange"
        deg = 1
        el_s = df.FiniteElement(el_str, degree=deg, cell=self.cell)
        el_v = df.VectorElement(el_str, degree=deg, cell=self.cell)
        V_s = df.FunctionSpace(self.mesh, el_s)
        V_v = df.FunctionSpace(self.mesh, el_v)

        # Heat source
        f_heat = df.interpolate(self.heat_source, V_s)
        self.__write_xdmf("f_heat", f_heat )

        # Mass source
        f_mass = df.interpolate(self.mass_source, V_s)
        self.__write_xdmf("f_mass", f_mass )

        # Body force
        f_body = df.interpolate(self.body_force, V_v)
        self.__write_xdmf("f_body", f_body )

    def __write_xdmf(self, name, field ):
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
            for degree in range(5):  # test until degree five
                # Writing symmetric tensors crashes.
                # Therefore project symmetric tensor in nonsymmetric space
                # This is only a temporary fix, see:
                # https://fenicsproject.discourse.group/t/...
                # ...writing-symmetric-tensor-function-fails/1136
                el_symm = df.TensorElement(
                    df.FiniteElement(
                        "Lagrange", df.triangle, degree + 1
                    ), symmetry={(0, 1): (1, 0)}
                )  # symmetric tensor element
                el_sol = field.ufl_function_space().ufl_element()
                if el_sol == el_symm:
                    # Remove symmetry with projection
                    field = df.project(
                        field, df.TensorFunctionSpace(
                            self.mesh, "Lagrange", degree + 1
                        )
                    )
                    break

            field.rename(name, name)

            print("Write {}".format(fname_xdmf))
            file.write(field, self.time)
