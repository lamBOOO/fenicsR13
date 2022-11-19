# pylint: disable=invalid-name,too-many-lines
# pylint: disable=not-callable
# noqa: E226

"""
Solver module, contains the Solver class.

For usage examples, see the :class:`solver.Solver` description.
"""

import sys
import copy
import time as time_module
import dolfin as df


class Solver:
    r"""
    Class to store the actual solver.


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

    """

    def __init__(self, params, mesh, time):
        """Initialize solver and setup variables from input parameters."""

        self.params = params  
        
        self.mesh = mesh.mesh
        self.regions = mesh.subdomains
        self.boundaries = mesh.boundaries
        self.cell = self.mesh.ufl_cell()
        self.time = time

        self.comm = df.MPI.comm_world
        self.rank = df.MPI.rank(self.comm)

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

        # read source terms
        self.heat_source = df.Expression(str(self.params["heat_source"]),degree=2)
        self.mass_source = df.Expression(str(self.params["mass_source"]),degree=2)

        # initialization
        self.output_folder = self.params["output_folder"] + "/"
        self.elems = {
            "theta": None,
            "p": None
        }
        self.poisson_elems = None
        self.poisson_fspaces = None
        self.form_lhs = None
        self.form_rhs = None
        self.sol = {
            "theta": None,
            "p": None
        }

    def assemble(self):
        r"""
        Assemble the weak form of the system.

        """

        # Setup elements for all fields
        cell = self.cell
        mesh = self.mesh

        e = "Lagrange"
        deg = 1
        # scalar variables
        self.elems["theta"] = df.FiniteElement(e, cell, deg)
        self.elems["p"] = df.FiniteElement(e, cell, deg)

        # Bundle elements per mode into `mxd_elems` dict
        poisson_elems = [ self.elems["theta"], self.elems["p"] ]
        self.poisson_elems = df.MixedElement(poisson_elems)
        self.poisson_fspaces = df.FunctionSpace(mesh, self.poisson_elems)

        # Get local variables
        regions = self.regions
        regs = self.regs
        boundaries = self.boundaries
        bcs = self.bcs

        # Define custom measures for boundary edges and inner edges
        df.dx = df.Measure("dx", domain=mesh, subdomain_data=regions)
        df.ds = df.Measure("ds", domain=mesh, subdomain_data=boundaries)
    
        # Setup trial and test functions
        w_poisson = self.poisson_fspaces
        (theta, p ) = df.TrialFunctions(w_poisson)
        (kappa, q ) = df.TestFunctions(w_poisson)

        # Setup source functions
        f_heat = self.heat_source
        f_mass = self.mass_source

        # Sub functionals:
        def a1(theta, kappa):
            return sum([(
                + regs[reg]["kn"] * df.inner(df.grad(theta), df.grad(kappa))
                + (1 / regs[reg]["kn"]) * df.inner(theta, kappa)
            ) * df.dx(reg) for reg in regs.keys()]) + sum([(
                + bcs[bc]["chi_tilde"] * kappa * theta
            ) * df.ds(bc) for bc in bcs.keys()])

        def a2(p, q):
            return sum([(
                + regs[reg]["kn"] * df.inner(df.grad(p), df.grad(q))
                + (1 / regs[reg]["kn"]) * df.inner(p, q)
            ) * df.dx(reg) for reg in regs.keys()]) + sum([(
                + bcs[bc]["chi_tilde"] * q * p
            ) * df.ds(bc) for bc in bcs.keys()])

        # Setup all equations
        A = [None] * 2
        L = [None] * 2
        # 1) Left-hand sides, bilinear form A[..]:
        A[0] = a1(kappa, theta) + 0
        A[1] = 0                + a2(q, p) 

        # 2) Right-hand sides, linear functional L[..]:
        L[0] = f_heat * kappa * df.dx + (
               + sum([( bcs[bc]["chi_tilde"]*bcs[bc]["theta_w"] * kappa ) * df.ds(bc) for bc in bcs.keys()]) )
        L[1] = f_mass * q * df.dx + (
               + sum([( bcs[bc]["chi_tilde"]*bcs[bc]["p_w"] * q ) * df.ds(bc) for bc in bcs.keys()]) )

        # Combine all equations to compound weak form and add stabilization
        self.form_lhs = sum(A)
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

        w = self.poisson_fspaces

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

        ( self.sol["theta"], self.sol["p"] ) = sol.split()


        # compute average theta
        vol = df.assemble(df.Constant(1) * df.dx)
        avgvel = df.assemble( abs(df.inner(self.sol["theta"], self.sol["theta"])) * df.dx ) / vol
        print("avg theta:", avgvel)

    def write(self):
        """
        Write all solver data to separate folder.

        This includes the writing of:

        (#) Solutions
        (#) Parameter functions
        (#) System matrices if set in input file
        """
        print("Write fields..")

        sols = self.sol
        for field in sols:
            if sols[field] is not None:
                self.__write_xdmf(field, sols[field] )


    def __write_xdmf(self, name, field ):
        """
        Write a given field to a XDMF file in the output folder.

        *Arguments*
            name
                The name to be given to the field. Will be used as filename
                and is the name of the field in e.g. Paraview.
            field
                The field to write.
        """
        fname_xdmf = (
            self.output_folder + name + "_" + str(self.time) + ".xdmf"
        )
        with df.XDMFFile(self.mesh.mpi_comm(), fname_xdmf) as file:

            field.rename(name, name)

            print("Write {}".format(fname_xdmf))
            file.write(field, self.time)
