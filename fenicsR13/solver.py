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
from petsc4py import PETSc
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
    ...     "tests/2d_heat/inputs/heat_01_coeffs_p1p1_stab.yml"
    ... ) # doctest: +ELLIPSIS
    Input:...
    >>> msh = H5Mesh("tests/2d_mesh/ring0.h5")
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
        self.mode = params["mode"]
        self.nsd = params["nsd"]
        self.comm = df.MPI.comm_world
        self.rank = df.MPI.rank(self.comm)

        # CIP
        self.use_cip = self.params["stabilization"]["cip"]["enable"]
        self.delta_theta = self.params["stabilization"]["cip"]["delta_theta"]
        self.delta_u = self.params["stabilization"]["cip"]["delta_u"]
        self.delta_p = self.params["stabilization"]["cip"]["delta_p"]

        # GLS
        self.use_gls = self.params["stabilization"]["gls"]["enable"]
        self.tau_energy = self.params["stabilization"]["gls"]["tau_energy"]
        self.tau_heatflux = self.params["stabilization"]["gls"]["tau_heatflux"]
        self.tau_mass = self.params["stabilization"]["gls"]["tau_mass"]
        self.tau_momentum = self.params["stabilization"]["gls"]["tau_momentum"]
        self.tau_stress = self.params["stabilization"]["gls"]["tau_stress"]

        self.petsc_options = self.params["petsc_options"]

        self.write_pdfs = self.params["postprocessing"]["write_pdfs"]
        self.write_vecs = self.params["postprocessing"]["write_vecs"]
        self.flows = self.params["postprocessing"]["flows"]
        self.line_integrals = self.params["postprocessing"]["line_integrals"]

        # Create region field expressions
        self.regs = copy.deepcopy(self.params["regs"])
        for reg_id in self.regs:
            for field in self.regs[reg_id].keys():
                self.regs[reg_id][field] = self.__createMacroScaExpr(
                    self.regs[reg_id][field]
                )

        # Create boundary field expressions
        self.bcs = copy.deepcopy(self.params["bcs"])
        for edge_id in self.bcs:
            for field in self.bcs[edge_id].keys():
                self.bcs[edge_id][field] = self.__createMacroScaExpr(
                    self.bcs[edge_id][field]
                )

        self.polar_system = self.params["polar_coord_syst"]

        self.heat_source = self.__createMacroScaExpr(self.params["heat_source"])
        self.mass_source = self.__createMacroScaExpr(self.params["mass_source"])
        self.body_force = self.__createMacroVecExpr(self.params["body_force"])
        self.f_s = self.__createMacroVecExpr(self.params["f_s"])
        self.f_sigma = self.__createMacroTenExpr(self.params["f_sigma"])

        self.exact_solution = self.params["convergence_study"]["exact_solution"]
        self.write_systemmatrix = self.params["convergence_study"][
            "write_systemmatrix"
        ]
        self.write_mpi_information = self.params["convergence_study"][
            "write_mpi_information"
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
        self.solver = None
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

    def __createSolMacroScaExpr(self, cpp_string):
        """
        Return a DOLFIN scalar expression with predefined macros after solve.

        These macros include:

        ============================ =========== ===============================
        Name                         Macro   CPP Replacement
        ============================ =========== ===============================
        ``theta``                    ``theta``   ``sol["theta"]``
        ``sx``                       ``sx``      ``sol["s"].split()[0]``
        ``sy``                       ``sy``      ``sol["s"].split()[1]``
        ``p``                        ``p``       ``sol["p"]``
        ``ux``                       ``ux``      ``sol["u"].split()[0]``
        ``uy``                       ``uy``      ``sol["u"].split()[1]``
        ``sigmaxx``                  ``sigmaxx`` ``sol["sigma"].split()[0]``
        ``sigmaxy``                  ``sigmaxy`` ``sol["sigma"].split()[1]``
        ``sigmayx``                  ``sigmayx`` ``sol["sigma"].split()[1]``
        ``sigmayy``                  ``sigmayy`` ``sol["sigma"].split()[2]``
        ============================ =========== ===============================
        """
        if self.nsd == 3:  # TODO: Implement this with proper sigma mapping
            raise Exception("__createSolMacroScaExpr not supported in 3D.")

        return df.Expression(
            str(cpp_string),
            degree=2,
            theta=self.sol["theta"],
            sx=self.sol["s"].split()[0],
            sy=self.sol["s"].split()[1],
            p=self.sol["p"],
            ux=self.sol["u"].split()[0],
            uy=self.sol["u"].split()[1],
            sigmaxx=self.sol["sigma"].split()[0],
            sigmaxy=self.sol["sigma"].split()[1],
            sigmayx=self.sol["sigma"].split()[1],
            sigmayy=self.sol["sigma"].split()[2]
        )

    def __createMacroScaExpr(self, cpp_string):
        """
        Return a DOLFIN scalar expression with predefined macros.

        These macros include:

        ============================ ======= =================================
        Name                         Macro   CPP Replacement
        ============================ ======= =================================
        Radius wrt. to :math:`(0,0)` ``R``   ``sqrt(pow(x[0],2)+pow(x[1],2))``
        Angle wrt. :math:`(0,0)`     ``phi`` ``atan2(x[1],x[0])``
        ============================ ======= =================================

        The following expressions are therefore equal:

        .. code-block:: python

            # expr1 is equal to expr2
            expr1 = self.__createMacroScaExpr("R*cos(phi)")
            expr2 = dolfin.Expression(
                "R*cos(phi)",
                degree=2,
                R=dolfin.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2),
                phi=dolfin.Expression("atan2(x[1],x[0])", degree=2),
            )
        """
        # TODO: Check for constants and then use df.Constant
        R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        phi = df.Expression("atan2(x[1],x[0])", degree=2)
        return df.Expression(
            str(cpp_string),
            degree=2,
            phi=phi,
            R=R
        )

    def __createMacroVecExpr(self, cpp_strings):
        """
        Return a DOLFIN vector expression with predefined macros.

        These macros include:

        ============================ ======= =================================
        Name                         Macro   CPP Replacement
        ============================ ======= =================================
        Radius wrt. to :math:`(0,0)` ``R``   ``sqrt(pow(x[0],2)+pow(x[1],2))``
        Angle wrt. :math:`(0,0)`     ``phi`` ``atan2(x[1],x[0])``
        ============================ ======= =================================

        See Also
        --------
        _Solver__createMacroScaExpr
        """
        R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        phi = df.Expression("atan2(x[1],x[0])", degree=2)
        # TODO Handle R, phi for 3D case => raise error?
        cpp_strings = [str(i) for i in cpp_strings]
        if len(cpp_strings) == 2:
            return df.Expression(
                cpp_strings,  # strange that no cast to list is needed
                degree=2,
                phi=phi,
                R=R
            )
        else:
            return df.Expression(cpp_strings, degree=2)

    def __createMacroTenExpr(self, cpp_strings):
        """
        Return a DOLFIN tensor expression with predefined macros.

        These macros include:

        ============================ ======= =================================
        Name                         Macro   CPP Replacement
        ============================ ======= =================================
        Radius wrt. to :math:`(0,0)` ``R``   ``sqrt(pow(x[0],2)+pow(x[1],2))``
        Angle wrt. :math:`(0,0)`     ``phi`` ``atan2(x[1],x[0])``
        ============================ ======= =================================

        See Also
        --------
        _Solver__createMacroScaExpr
        """
        R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        phi = df.Expression("atan2(x[1],x[0])", degree=2)
        # TODO Handle R, phi for 3D case => raise error?
        cpp_strings = [[str(el) for el in firstdim] for firstdim in cpp_strings]
        if len(cpp_strings) == 2:
            return df.Expression(
                cpp_strings,  # strange that no cast to list is needed
                degree=2,
                phi=phi,
                R=R
            )
        else:
            return df.Expression(cpp_strings, degree=2)

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
                if self.nsd == 2:
                    self.elems[var] = df.TensorElement(
                        e, cell, deg, symmetry={(1, 0): (0, 1)}
                    )
                else:  # nsd==3 in this case
                    self.elems[var] = df.TensorElement(
                        e, cell, deg,
                        symmetry={  # the ordering influences sigma.split()
                            (1, 0): (0, 1),
                            (2, 0): (0, 2),
                            (2, 1): (1, 2),
                            (2, 2): (0, 0)
                            # (needed for sigma treatment below)
                            # sigma.split()[i] map:
                            # [
                            #     [0, 1, 2]
                            #     [1, 3, 4]
                            #     [2, 4, ?not_important?]
                            # ]
                        }
                    )
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

    def __check_regions(self):
        """
        Check if all regions from the input mesh have params prescribed.

        Raises an exception if one region is missing.
        """
        # TODO: Implement this
        region_ids = self.regions.array()
        regs_specified = list(self.regs.keys())

        for reg_id in region_ids:
            if reg_id not in [0] + regs_specified:  # inner zero allowed
                raise Exception("Mesh region {} has no params!".format(reg_id))

    def __check_bcs(self):
        """
        Check if all boundaries from the input mesh have BCs prescribed.

        Raises an exception if one BC is missing.
        """
        boundary_ids = self.boundaries.array()
        bcs_specified = list(self.bcs.keys())

        for edge_id in boundary_ids:
            if edge_id not in [0] + bcs_specified:  # inner zero allowed
                raise Exception("Mesh edge {} has no bcs!".format(edge_id))

    def normal(self):
        mesh1 = self.mesh
        v1 = df.as_vector([1.0, 0.0, 0.0])
        v2 = df.as_vector([0.0, 1.0, 0.0])
        n_vec = df.FacetNormal(mesh1)
        t_vec1 = df.conditional(
            df.gt(abs(n_vec[0]), np.finfo(float).eps),
            df.cross(n_vec, v2 / df.sqrt(n_vec[0] ** 2 + n_vec[2] ** 2)),
            df.cross(n_vec, v1 / df.sqrt(n_vec[1] ** 2 + n_vec[2] ** 2))
        )
        t_vec2 = df.cross(n_vec, t_vec1)
        return n_vec, t_vec1, t_vec2

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
        # Check if all mesh regions have params prescribed
        self.__check_regions()

        # Check if all mesh boundaries have bcs prescribed
        self.__check_bcs()

        # Setup required function spaces
        self.__setup_function_spaces()

        # Get local variables
        mesh = self.mesh
        regions = self.regions
        nsd = self.nsd
        regs = self.regs
        boundaries = self.boundaries
        bcs = self.bcs
        delta_theta = df.Constant(self.delta_theta)
        delta_u = df.Constant(self.delta_u)
        delta_p = df.Constant(self.delta_p)
        tau_energy = df.Constant(self.tau_energy)
        tau_heatflux = df.Constant(self.tau_heatflux)
        tau_mass = df.Constant(self.tau_mass)
        tau_momentum = df.Constant(self.tau_momentum)
        tau_stress = df.Constant(self.tau_stress)

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

        sigma = to.gen3DTFdim3(sigma)
        psi = to.gen3DTFdim3(psi)

        # Setup source functions
        f_heat = self.heat_source
        f_mass = self.mass_source
        f_body = self.body_force
        f_sigma = self.f_sigma
        f_s = self.f_s

        # Decouple heat/stress switch
        if self.mode == "r13":
            cpl = 1
        else:
            cpl = 0

        # Stabilization
        if self.use_cip:
            cip = 1
        else:
            cip = 0
        if self.use_gls:
            gls = 1
        else:
            gls = 0

        # Setup normal/tangential projections
        # => tangential (tx,ty) = (-ny,nx) = perp(n) only for 2D

        if nsd == 2:
            n_vec = df.FacetNormal(mesh)
            t_vec1 = ufl.perp(n_vec)
            t_vec2 = df.as_vector([0, 0])  # the second tangent is zeroed.
        else:  # defnining normals and tangents for 3D case
            n_vec, t_vec1, t_vec2 = self.normal()

        def n(rank1):
            return df.dot(rank1, n_vec)

        def t1(rank1):
            return df.dot(rank1, t_vec1)

        def t2(rank1):
            return df.dot(rank1, t_vec2)

        def nn(rank2):
            return df.dot(rank2 * n_vec, n_vec)

        def t1t1(rank2):
            return df.dot(rank2 * t_vec1, t_vec1)

        def nt1(rank2):
            return df.dot(rank2 * n_vec, t_vec1)

        def nt2(rank2):
            return df.dot(rank2 * n_vec, t_vec2)

        def t1t2(rank2):
            return df.dot(rank2 * t_vec1, t_vec2)

        if self.mode == "heat":
            # Note that the heat system does not include the Δ variable from [1]
            # Reason for this was, that I didn't have an exact solution for that
            # The heat system historically was the first one implemented and is
            # documented in [2].
            # Keep in mind: If f_heat = 0 => delta = 0 anyways...
            # [1]: fenicsR13 paper: https://dx.doi.org/10.1145/3442378
            # [2]: M.Sc thesis: https://doi.org/10.18154/RWTH-2023-12262
            incl_delta = 0  # TODO: Include (adapt Mathemathica scripts & esol)
        else:
            incl_delta = 1

        # Sub functionals:
        # 1) Diagonals:

        def a(s, r):
            # Notes:
            # - We expanded the stf(grad(s))-term, therefore the constants read:
            #   4/5-24/75 = (60-24)/75 = 36/75 = 12/25
            return sum([(
                # => 24/25*stf(grad)*grad
                + 24 / 25 * regs[reg]["kn"] * df.inner(
                    df.sym(df.grad(s)), df.sym(df.grad(r))
                )
                - 24 / 75 * regs[reg]["kn"] * df.div(s) * df.div(r)
                + 4 / 5 * incl_delta * regs[reg]["kn"] * df.div(s) * df.div(r)
                + 4 / 15 * (1 / regs[reg]["kn"]) * df.inner(s, r)
            ) * df.dx(reg) for reg in regs.keys()]) + sum([(
                + 1 / (2 * bcs[bc]["chi_tilde"]) * n(s) * n(r)
                + 11 / 25 * bcs[bc]["chi_tilde"] * t1(s) * t1(r)
                + cpl * 1 / 25 * bcs[bc]["chi_tilde"] * t1(s) * t1(r)
                + 11 / 25 * bcs[bc]["chi_tilde"] * t2(s) * t2(r)
                + cpl * 1 / 25 * bcs[bc]["chi_tilde"] * t2(s) * t2(r)
                # + df.inner(
                #     ufl.transpose(df.as_tensor(
                #         [n_vec, t_vec1]
                #     )) *
                #     df.as_tensor(np.diag(
                #         [
                #             1 / (2 * bcs[bc]["chi_tilde"]),
                #             (11+cpl)/25 * bcs[bc]["chi_tilde"]
                #         ]
                #     )) * df.as_tensor(
                #         [n_vec, t_vec1]
                #     ) * s, r
                # )
            ) * df.ds(bc) for bc in bcs.keys()])

        def d(si, ps):
            # Notes:
            # 21/20+3/40=45/40=9/8
            return sum([(
                + regs[reg]["kn"] * df.inner(
                    to.stf3d3(to.grad3dOf2(to.gen3DTFdim2(si), nsd)),
                    to.stf3d3(to.grad3dOf2(to.gen3DTFdim2(ps), nsd))
                )
                # TODO: this one is equivalent to the above, document it!
                # + regs[reg]["kn"] * df.inner(
                #     to.sym3d3(to.grad3dOf2(to.gen3DTFdim2(si), nsd)),
                #     to.sym3d3(to.grad3dOf2(to.gen3DTFdim2(ps), nsd))
                # )
                # - regs[reg]["kn"] * 4 / 15 * df.inner(
                #     df.div(si),
                #     df.div(ps)
                # )
                + (1 / (2 * regs[reg]["kn"])) * df.inner(
                    to.gen3DTFdim2(si), to.gen3DTFdim2(ps)
                )
            ) * df.dx(reg) for reg in regs.keys()]) + sum([(
                + bcs[bc]["chi_tilde"] * 21 / 20 * nn(si) * nn(ps)
                + bcs[bc]["chi_tilde"] * cpl * 3 / 40 * nn(si) * nn(ps)
                + bcs[bc]["chi_tilde"] * (
                    (t1t1(si) + (1 / 2) * nn(si)) *
                    (t1t1(ps) + (1 / 2) * nn(ps))
                )
                + (1 / bcs[bc]["chi_tilde"]) * nt1(si) * nt1(ps)
                + (1 / bcs[bc]["chi_tilde"]) * nt2(si) * nt2(ps)
                + bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * nn(si) * nn(ps)
                + bcs[bc]["chi_tilde"] * t1t2(si) * t1t2(ps)
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
        #     return sum([(  # form 2
        #         - df.inner(r, df.grad(th))
        #     ) * df.dx(reg) for reg in regs.keys()]) + sum([(
        #         th * n(r)
        #     ) * df.ds(bc) for bc in bcs.keys()])

        def c(r, si):
            return cpl * (sum([(
                2 / 5 * df.inner(si, df.grad(r))
            ) * df.dx(reg) for reg in regs.keys()]) - sum([(
                3 / 20 * nn(si) * n(r)
                + 1 / 5 * nt1(si) * t1(r)
                + 1 / 5 * nt2(si) * t2(r)
            ) * df.ds(bc) for bc in bcs.keys()]))

        def e(u, ps):
            return sum([(
                df.dot(df.div(ps), u)
            ) * df.dx(reg) for reg in regs.keys()])
            # return sum([(  # form 2
            #     - df.inner(ps, df.grad(u))
            # ) * df.dx(reg) for reg in regs.keys()]) + sum([(
            #     df.dot(u, n(ps))
            # ) * df.ds(bc) for bc in bcs.keys()])

        def f(p, ps):
            return sum([(
                bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * p * nn(ps)
            ) * df.ds(bc) for bc in bcs.keys()])

        def g(p, v):
            return sum([(
                df.inner(v, df.grad(p))
            ) * df.dx(reg) for reg in regs.keys()])
        #     return sum([(  # form 2
        #         - df.div(v) * p
        #     ) * df.dx(reg) for reg in regs.keys()]) + sum([(
        #         p * n(v)
        #     ) * df.ds(bc) for bc in bcs.keys()])

        # 3.1) CIP Stabilization:
        def j_theta(theta, kappa):
            return (
                + delta_theta * h_avg**3 *
                df.jump(df.grad(theta), n_vec) * df.jump(df.grad(kappa), n_vec)
            ) * df.dS

        def j_u(u, v):
            return (
                + delta_u * h_avg**3 *
                df.dot(df.jump(df.grad(u), n_vec), df.jump(df.grad(v), n_vec))
            ) * df.dS

        def j_p(p, q):
            return (
                + delta_p * h_avg *
                df.jump(df.grad(p), n_vec) * df.jump(df.grad(q), n_vec)
            ) * df.dS

        # 3.2) GLS Stabilization
        def gls_heat(theta, kappa, s, r):
            return sum([(
                tau_energy * h_msh**1 * (
                    df.inner(
                        df.div(s) + cpl * df.div(u) - f_heat,
                        df.div(r) + cpl * df.div(v)
                    )
                )  # energy
                + tau_heatflux * h_msh**1 *
                df.inner(
                    (5 / 2) * df.grad(theta)
                    - (12 / 5) * regs[reg]["kn"] * df.div(to.stf3d2(df.grad(s)))
                    - (1 / 6) * regs[reg]["kn"] * 12 * df.grad(df.div(s))
                    + 1 / regs[reg]["kn"] * 2 / 3 * s,
                    (5 / 2) * df.grad(kappa)
                    - (12 / 5) * regs[reg]["kn"] * df.div(to.stf3d2(df.grad(r)))
                    - (1 / 6) * regs[reg]["kn"] * 12 * df.grad(df.div(r))
                    + 1 / regs[reg]["kn"] * 2 / 3 * r
                )  # heatflux
            ) * df.dx(reg) for reg in regs.keys()])

        def gls_stress(p, q, u, v, sigma, psi):
            return sum([(
                tau_mass * h_msh**1.5 *
                df.inner(
                    df.div(v), df.div(u) - f_mass
                )  # mass
                + tau_momentum * h_msh**1.5 *
                df.inner(
                    df.grad(q) + df.div(psi),
                    df.grad(p) + df.div(sigma) - f_body
                )  # momentum
                + tau_stress * h_msh**1.5 *
                df.inner(
                    cpl * (4 / 5) * to.gen3DTFdim2(df.grad(r))
                    + 2 * to.stf3d2(to.gen3d2(df.grad(v)))
                    - 2 * regs[reg]["kn"] * to.div3d3(
                        to.stf3d3(to.grad3dOf2(to.gen3DTFdim2(psi), nsd))
                    )
                    + (1 / regs[reg]["kn"]) * to.gen3DTFdim2(psi),
                    cpl * (4 / 5) * to.gen3DTFdim2(df.grad(s))
                    + 2 * to.stf3d2(to.gen3d2(df.grad(u)))
                    - 2 * regs[reg]["kn"] * to.div3d3(
                        to.stf3d3(to.grad3dOf2(to.gen3DTFdim2(sigma), nsd))
                    )
                    + (1 / regs[reg]["kn"]) * to.gen3DTFdim2(sigma)
                )  # stress
            ) * df.dx(reg) for reg in regs.keys()])

        v1 = {}  # boundary velocity vector

        if nsd == 2:
            if self.polar_system is True:  # Polar implementation
                for bc in bcs.keys():
                    v1.update({
                        bc: bcs[bc]["u_n_w"] * n_vec + bcs[bc]["u_t_w"] * t_vec1
                    })
            else:  # Cartesian implementation
                for bc in bcs.keys():
                    v1.update({
                        bc: df.as_vector([bcs[bc]["u_x_w"], bcs[bc]["u_y_w"]])
                    })
        else:
            if self.polar_system is True:
                raise Exception("Polar coordinates are NOT supported in 3D")
            for bc in bcs.keys():
                v1.update({
                    bc: df.as_vector([
                        bcs[bc]["u_x_w"], bcs[bc]["u_y_w"], bcs[bc]["u_z_w"]
                    ])
                })
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
        L[0] = + df.dot(f_s, r) * df.dx - sum([(
            bcs[bc]["theta_w"] * n(r)
        ) * df.ds(bc) for bc in bcs.keys()])
        # Use div(u)=f_mass to remain sym. (density-form doesnt need this):
        L[1] = (f_heat - cpl * f_mass) * kappa * df.dx
        L[2] = + df.inner(
            to.gen3DTFdim2(f_sigma), to.gen3DTFdim2(psi)
        ) * df.dx - sum([(
            + t1(v1[bc]) * nt1(psi)
            + t2(v1[bc]) * nt2(psi)
            + (
                + n(v1[bc])
                - bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * bcs[bc]["p_w"]
            ) * nn(psi)
        ) * df.ds(bc) for bc in bcs.keys()])
        L[3] = + df.dot(f_body, v) * df.dx
        L[4] = + (f_mass * q) * df.dx - sum([(
            (
                + n(v1[bc])
                - bcs[bc]["epsilon_w"] * bcs[bc]["chi_tilde"] * bcs[bc]["p_w"]
            ) * q
        ) * df.ds(bc) for bc in bcs.keys()])

        # Combine all equations to compound weak form and add stabilization
        if self.mode == "heat":
            self.form_lhs = sum(A[0:2]) + (
                cip * (j_theta(theta, kappa))
                + gls * df.lhs(gls_heat(theta, kappa, s, r))
            )
            self.form_rhs = sum(L[0:2]) + (
                gls * df.rhs(gls_heat(theta, kappa, s, r))
            )
        elif self.mode == "stress":
            self.form_lhs = sum(A[2:5]) + (
                cip * (j_u(u, v) + j_p(p, q))
                + gls * df.lhs(gls_stress(p, q, u, v, sigma, psi))
            )
            self.form_rhs = sum(L[2:5]) + (
                gls * df.rhs(gls_stress(p, q, u, v, sigma, psi))
            )
        elif self.mode == "r13":
            self.form_lhs = sum(A) + (
                cip * (j_theta(theta, kappa) + j_u(u, v) + j_p(p, q))
                + gls * df.lhs(
                    gls_heat(theta, kappa, s, r)
                    + gls_stress(p, q, u, v, sigma, psi)
                )
            )
            self.form_rhs = sum(L) + (
                gls * df.rhs(
                    gls_heat(theta, kappa, s, r)
                    + gls_stress(p, q, u, v, sigma, psi)
                )
            )

    def solve(self):
        """
        Solve the previously assembled system.

        Note that setting the mg levels for `gamg` doesn't work (see here_).
        Also see some note here2_.

        _here: https://lists.mcs.anl.gov/mailman/htdig/petsc-users/2019-April/038223.html  # noqa
        _here2: https://www.firedrakeproject.org/demos/saddle_point_systems.py.html  # noqa

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
        if self.mode == "heat":
            w = self.mxd_fspaces["heat"]
        elif self.mode == "stress":
            w = self.mxd_fspaces["stress"]
        elif self.mode == "r13":
            w = self.mxd_fspaces["r13"]

        print("Start assemble")
        sys.stdout.flush()
        start_t = time_module.time()
        AA = df.assemble(self.form_lhs)
        LL = df.assemble(self.form_rhs)
        end_t = time_module.time()
        secs = end_t - start_t
        self.write_content_to_file("assemble", secs)
        print("Finish assemble: {}".format(str(secs)))
        sys.stdout.flush()

        print("Start solve")
        sys.stdout.flush()
        start_t = time_module.time()
        sol = df.Function(w)

        # TODO: Implement nonzero init guess
        # df.parameters['krylov_solver']['nonzero_initial_guess'] = True
        # w.dofmap().set(sol.vector(), 1.0)
        # w.dofmap().set(sol.split()[2].vector(), 0.0)
        # df.solve(
        #     AA, sol.vector(), LL, solver_name, preconditioner
        # )
        # TODO: Test this
        # # Create Krylov solver
        # solver = df.PETScKrylovSolver("bicgstab", "amg")
        # solver.set_operator(AA)
        # # Create vector that spans the null space and normalize
        # null_vec = df.Vector(sol.vector())
        # w.dofmap().set(null_vec, 1.0)
        # null_vec *= 1.0/null_vec.norm("l2")
        # # Create null space basis object and attach to PETSc matrix
        # null_space = df.VectorSpaceBasis([null_vec])
        # df.as_backend_type(AA).set_nullspace(null_space)
        # null_space.orthogonalize(LL);
        # TODO: Test this also (problem infered preconditioner)
        # See: https://fenicsproject.org/olddocs/dolfin/1.6.0/python/demo/...
        # documented/stokes-iterative/python/documentation.html
        # self.form_p = None  # add above...
        # self.form_p = d(sigma, psi) + df.dot(u, v) * df.dx + (
        #     df.inner(df.grad(p), df.grad(q)) * df.dx
        # )  # add above...
        # PP = df.assemble(self.form_p)
        # solver.parameters['relative_tolerance'] = 1E-2
        # solver.solve(sol.vector(), LL)
        # solver = df.PETScKrylovSolver(solver_name, "jacobi")
        # solver.set_operators(AA, PP)
        # solver.parameters['monitor_convergence'] = True
        # solver.parameters['relative_tolerance'] = 1E-2

        # Use PETSc solver with provided arguments
        self.solver = df.PETScKrylovSolver()
        opts = PETSc.Options()
        if "log_view" in self.petsc_options.keys():
            PETSc.Log().begin()
        for key in self.petsc_options:
            opts[key] = self.petsc_options[key]
        print(opts.view())
        self.solver.set_from_options()
        self.solver.set_operator(AA)
        self.solver.solve(sol.vector(), LL)

        # Print matrix (serial only)
        # A_mat = self.solver.ksp().getOperators()[0]
        # print(
        #     A_mat.getValues(
        #         range(A_mat.getSizes()[0][0]), range(A_mat.getSizes()[1][0])
        #     )
        # )

        end_t = time_module.time()
        secs = end_t - start_t
        self.write_content_to_file("solve", secs)
        print("Finished solve: {}".format(str(secs)))
        sys.stdout.flush()

        if self.mode == "heat":
            (self.sol["theta"], self.sol["s"]) = sol.split()
        elif self.mode == "stress":
            (self.sol["p"], self.sol["u"], self.sol["sigma"]) = sol.split()
        elif self.mode == "r13":
            (
                self.sol["theta"], self.sol["s"],
                self.sol["p"], self.sol["u"], self.sol["sigma"]
            ) = sol.split()
        # Special treatment for sigma
        if self.mode == "stress" or self.mode == "r13":
            print("Start sigma projection")
            sys.stdout.flush()
            start_t = time_module.time()

            # This can sometimes fail in parallel & is very expensive
            # self.sol["sigma"] = df.project(
            #     to.gen3DTFdim3(self.sol["sigma"]), df.TensorFunctionSpace(
            #         self.mesh, "Lagrange",
            #         self.params["elements"]["sigma"]["degree"]
            #     ), solver_type=solver_name, preconditioner_type=preconditioner
            # )
            # -> See workaroud below...

            # This is a workaround: Make sigma a proper STF tensor
            sp = df.TensorFunctionSpace(
                self.mesh, self.params["elements"]["sigma"]["shape"],
                self.params["elements"]["sigma"]["degree"]
            )
            if self.nsd == 2:
                sigmaxx, sigmaxy, sigmayy = self.sol["sigma"].split()
                # strangely, the following is needed:
                sigmaxx.set_allow_extrapolation(True)
                sigmaxy.set_allow_extrapolation(True)
                sigmayy.set_allow_extrapolation(True)
                self.sol["sigma"] = df.interpolate(
                    df.Expression(
                        [
                            ["sigmaxx", "sigmaxy"],
                            ["sigmaxy", "sigmayy"],
                        ],
                        degree=2,
                        sigmaxx=sigmaxx,
                        sigmaxy=sigmaxy,
                        sigmayy=sigmayy
                    ), sp
                )
            elif self.nsd == 3:
                sigmaxx, sigmaxy, sigmaxz, sigmayy, sigmayz = (
                    self.sol["sigma"].split()
                )
                sigmaxx.set_allow_extrapolation(True)
                sigmaxy.set_allow_extrapolation(True)
                sigmaxz.set_allow_extrapolation(True)
                sigmayy.set_allow_extrapolation(True)
                sigmayz.set_allow_extrapolation(True)
                self.sol["sigma"] = df.interpolate(
                    df.Expression(
                        [
                            ["sigmaxx", "sigmaxy", "sigmaxz"],
                            ["sigmaxy", "sigmayy", "sigmayz"],
                            ["sigmaxz", "sigmayz", "-sigmaxx-sigmayy"]
                        ],
                        degree=2,
                        sigmaxx=sigmaxx,
                        sigmaxy=sigmaxy,
                        sigmaxz=sigmaxz,
                        sigmayy=sigmayy,
                        sigmayz=sigmayz
                    ), sp
                )
            end_t = time_module.time()
            secs = end_t - start_t
            print("Finished sigma projection: {}".format(str(secs)))
            sys.stdout.flush()

        if self.mode == "stress" or self.mode == "r13":
            if self.rescale_p is False:
                print("No pressure rescaling applied.")
            else:
                p_i = df.interpolate(self.sol["p"], self.fspaces["p"])
                shift = 0.0
                if self.rescale_p == "zeromean":
                    shift = self.__calc_sf_mean(p_i)
                elif self.rescale_p == "zerominimum":
                    shift = self.__calc_sf_min(p_i)
                else:
                    raise RuntimeError("Wrong rescale_pressure option.")
                print("Pressure rescaling with shift=", shift)
                shift_function = df.Function(self.fspaces["p"])
                shift_function.assign(df.Constant(shift))
                p_i.assign(p_i - shift_function)
                self.sol["p"] = p_i

        # Calculate mass flows
        for bc_id in self.flows:
            if bc_id not in self.boundaries.array():
                raise Exception(
                    "Flow calculation: {} is no boundary.".format(bc_id)
                )
            n = df.FacetNormal(self.mesh)
            mass_flow_rate = df.assemble(
                df.inner(self.sol["u"], n) * df.ds(bc_id)
            )
            print("mass flow rate of BC", bc_id, ":", mass_flow_rate)
            self.write_content_to_file("massflow_" + str(bc_id), mass_flow_rate)
            # TODO: Implement generically
            heat_flow_rate = df.assemble(
                df.inner(self.sol["s"], n) * df.ds(bc_id)
            )
            print("heat flow rate of BC", bc_id, ":", heat_flow_rate)
            self.write_content_to_file("heatflow_" + str(bc_id), heat_flow_rate)

        if self.mode == "stress" or self.mode == "r13":
            vol = df.assemble(df.Constant(1) * df.dx)
            avgvel = df.assemble(
                abs(df.inner(self.sol["u"], self.sol["u"])) * df.dx
            ) / vol
            print("avgvel:", avgvel)
            self.write_content_to_file("avgvel", avgvel)

            kinetic_energy_velocity = df.assemble(
                0.5 * df.inner(self.sol["u"], self.sol["u"]) * df.dx
            )
            print("kinetic_energy_velocity:", kinetic_energy_velocity)
            self.write_content_to_file(
                "kinetic_energy_velocity", kinetic_energy_velocity
            )

            u_vals = self.sol["u"].compute_vertex_values()
            max_norm_u = np.max(
                np.linalg.norm(u_vals.reshape(-1, self.nsd), axis=1)
            )
            print("max_norm_u:", max_norm_u)
            self.write_content_to_file("max_norm_u", max_norm_u)

        def __line_integral_average(u, A, B, n):
            """
            Integrate u over segment [A, B] partitioned into n elements.

            See: https://fenicsproject.org/qa/13863/integrate-along-axes-lines/
            """
            assert u.value_rank() == 0
            assert len(A) == len(B) > 1 and np.linalg.norm(A - B) > 0
            assert n > 0

            # Mesh line for integration
            mesh_points = [A + t * (B - A) for t in np.linspace(0, 1, n + 1)]
            tdim, gdim = 1, len(A)

            # Create line mesh
            mesh = df.Mesh(df.MPI.comm_world)
            if self.rank == 0:
                editor = df.MeshEditor()
                editor.open(mesh, "interval", tdim, gdim)
                editor.init_vertices(n + 1)
                editor.init_cells(n)
                for vi, v in enumerate(mesh_points):
                    editor.add_vertex(vi, v)
                for ci in range(n):
                    editor.add_cell(ci, np.array([ci, ci + 1], dtype='uintp'))
                editor.close()
            df.MeshPartitioning.build_distributed_mesh(mesh)

            # Setup function space
            elm = u.function_space().ufl_element()
            family = elm.family()
            degree = elm.degree()
            V = df.FunctionSpace(mesh, family, degree)
            v = df.Function(V)

            # Interpolate to line mesh
            # df.interpolate does not work in parallel from different meshes
            # -> https://fenicsproject.org/docs/dolfin/1.6.0/python/
            #    programmers-reference/cpp/function/LagrangeInterpolator.html
            # Use LagrangeInterpolator, which works in parallel
            # -> https://fenicsproject.org/qa/7758/
            #    function-interpolation-in-parallel/?show=12144#c12144
            df.LagrangeInterpolator.interpolate(v, u)

            custom_dx = df.Measure("dx", domain=mesh)

            return df.assemble(v * custom_dx)

        print("Start line_integrals")
        start_t = time_module.time()
        sys.stdout.flush()
        for li in self.line_integrals:
            name = li["name"]

            # Interpolate in "best" Lagrange space
            maxdeg = max([
                self.params["elements"]["theta"]["degree"],
                self.params["elements"]["s"]["degree"],
                self.params["elements"]["p"]["degree"],
                self.params["elements"]["u"]["degree"],
                self.params["elements"]["sigma"]["degree"]
            ])
            expr = df.interpolate(
                self.__createSolMacroScaExpr(li["expr"]),
                df.FunctionSpace(self.mesh, "Lagrange", maxdeg)
            )
            start = np.array(li["start"])
            end = np.array(li["end"])
            res = li["res"]
            result = __line_integral_average(
                expr, start, end, res
            )
            self.write_content_to_file(
                "line_integral_{}".format(name), result
            )
        end_t = time_module.time()
        secs = end_t - start_t
        print("Finish line_integrals: {}".format(str(secs)))
        sys.stdout.flush()

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
        """
        mean_val = (
            df.assemble(scalar_function * df.dx) / (
                df.assemble(df.Constant(1.0) * df.dx)
            )
        )
        print("Calculated mean value:", mean_val)
        return mean_val

    def __calc_sf_min(self, scalar_function):
        """
        Calculate the minimum of a scalar function (works in parallel).
        """
        v = scalar_function.compute_vertex_values()
        min_val = np.min(v)
        min_val = df.MPI.min(self.comm, min_val)
        print("Calculated minimum value:", min_val)
        return min_val

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

        The following is very slow, I don't know why:
        .. code-block:: python
            err_f_L2_full = df.errornorm(
                field_e_i, field_i, "L2", degree_rise=1
            )

        References
        ----------
        .. [6] `DOLFIN documentation <https://fenicsproject.org/docs/dolfin/>`_

        """
        field_e_i = df.interpolate(field_e_, v_field)
        field_i = df.interpolate(field_, v_field)
        error = df.Function(v_field)
        error.assign(field_e_i - field_i)

        err_f_L2_full = df.norm(error, "L2")
        err_v_linf_full = df.MPI.max(self.comm, np.max(np.abs(
            error.compute_vertex_values()
        )))
        err_f_H1_full = df.norm(error, "H1")

        dofs = len(field_e_i.split()) or 1

        if dofs == 1:
            # scalar
            errs_f_L2 = [df.norm(error, "L2")]
            errs_v_linf = [df.MPI.max(self.comm, np.max(np.abs(
                error.compute_vertex_values()
            )))]
            errs_f_H1 = [df.norm(error, "H1")]
        else:
            # vector or tensor
            errs_f_L2 = [df.norm(error.split()[i], "L2") for i in range(dofs)]
            errs_v_linf = [
                df.MPI.max(self.comm, np.max(np.abs(
                    error.split()[i].compute_vertex_values()
                )))
                for i in range(dofs)
            ]
            errs_f_H1 = [df.norm(error.split()[i], "H1") for i in range(dofs)]

        if self.relative_error:
            def wmsg(n):
                return f"WARN: {name_} has zero {n}-norm, abs. error used"
            norm_f_L2_full = (
                df.norm(field_e_i, "L2") or any([1, print(wmsg("L2"))])
            )
            norm_v_linf_full = (
                df.MPI.max(self.comm, np.max(
                    np.abs(field_e_i.compute_vertex_values())
                )) or any([1, print(wmsg("linf"))])
            )
            norm_f_H1_full = (
                df.norm(field_e_i, "H1") or any([1, print(wmsg("L2"))])
            )
            if dofs == 1:
                # scalar
                norms_f_L2 = [
                    df.norm(field_e_i, "L2") or any([1, print(wmsg("L2"))])
                ]
                norms_v_linf = [
                    df.MPI.max(self.comm, np.max(
                        np.abs(field_e_i.compute_vertex_values())
                    ))
                    or any([1, print(wmsg("linf"))])
                ]
                norms_f_H1 = [
                    df.norm(field_e_i, "H1") or any([1, print(wmsg("H1"))])
                ]
            else:
                # vector or tensor
                norms_f_L2 = [df.norm(
                    field_e_i.split()[i], "L2"
                ) or any([1, print(wmsg("L2"))]) for i in range(dofs)]
                norms_v_linf = [
                    df.MPI.max(self.comm, np.max(
                        np.abs(field_e_i.split()[i].compute_vertex_values())
                    ) or any([1, print(wmsg("linf"))]))
                    for i in range(dofs)
                ]
                norms_f_H1 = [df.norm(
                    field_e_i.split()[i], "H1"
                ) or any([1, print(wmsg("H1"))]) for i in range(dofs)]
            err_f_L2_full = err_f_L2_full / norm_f_L2_full
            err_v_linf_full = err_v_linf_full / norm_v_linf_full
            err_f_H1_full = err_f_H1_full / norm_f_H1_full
            errs_f_L2 = [x / y for x, y in zip(errs_f_L2, norms_f_L2)]
            errs_v_linf = [x / y for x, y in zip(errs_v_linf, norms_v_linf)]
            errs_f_H1 = [x / y for x, y in zip(errs_f_H1, norms_f_H1)]

        print("Calculate errors..")
        print("Full error " + str(name_) + " L_2:", err_f_L2_full)
        print("Full error " + str(name_) + " l_inf:", err_v_linf_full)
        print("Full error " + str(name_) + " H_1:", err_f_H1_full)
        print("Component Error " + str(name_) + " L_2:", errs_f_L2)
        print("Component Error " + str(name_) + " l_inf:", errs_v_linf)
        print("Component Error " + str(name_) + " H_1:", errs_f_H1)

        self.__write_xdmf(name_ + "_e", field_e_i, False)
        self.__write_xdmf(name_ + "_error", error, False)

        return (
            [{
                "L_2": err_f_L2_full,
                "l_inf": err_v_linf_full,
                "H_1": err_f_H1_full,
            }]
            +
            [{
                "L_2": errs_f_L2[i],
                "l_inf": errs_v_linf[i],
                "H_1": errs_f_H1[i],
            } for i in range(dofs)]
        )

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
        print("Start error calculation")
        sys.stdout.flush()
        start_t = time_module.time()

        cell = self.cell
        msh = self.mesh
        e = self.params["elements"]["sigma"]["shape"]
        deg = self.params["elements"]["sigma"]["degree"]
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
            ers["s"] = ve[0]
            ers["sx"] = ve[1]
            ers["sy"] = ve[2]
            if self.nsd > 2:
                ers["sz"] = ve[3]

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
                df.FunctionSpace(msh, df.TensorElement(
                    e, cell, deg)), "sigma"
            )
            ers = self.errors
            if self.nsd == 2:
                ers["p"] = se[0]
                ers["u"] = ve[0]
                ers["ux"] = ve[1]
                ers["uy"] = ve[2]
                ers["sigma"] = te[0]
                ers["sigmaxx"] = te[1]
                ers["sigmaxy"] = te[2]
                ers["sigmayy"] = te[4]
            else:
                ers["p"] = se[0]
                ers["u"] = ve[0]
                ers["ux"] = ve[1]
                ers["uy"] = ve[2]
                ers["uz"] = ve[3]
                ers["sigma"] = te[0]
                ers["sigmaxx"] = te[1]
                ers["sigmaxy"] = te[2]
                ers["sigmaxz"] = te[3]
                ers["sigmayy"] = te[5]
                ers["sigmayz"] = te[6]

        end_t = time_module.time()
        secs = end_t - start_t
        print("Finished error calculation: {}".format(str(secs)))
        sys.stdout.flush()

        return self.errors

    def write_content_to_file(self, filename, content):
        """Write content to a file in the output folder."""
        path = self.output_folder + filename + "_" + str(self.time)
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
        if self.comm.size == 1:
            print("Write fields..")
            sys.stdout.flush()
        else:
            print("Write fields.. (MPI bug: can fail too small meshes)")
            sys.stdout.flush()

        self.__write_solutions()
        if self.write_vecs:
            self.__write_vecs()
        if self.nsd == 2:
            self.__write_parameters()
        if self.write_systemmatrix:
            self.__write_discrete_system()
        if self.write_mpi_information:
            self.__write_mpi_information()

    def __write_solutions(self):
        """Write all solutions fields."""
        sols = self.sol
        for field in sols:
            if sols[field] is not None:
                self.__write_xdmf(field, sols[field], self.write_pdfs)

    def __write_vecs(self):
        """Write all solutions to vectors."""
        sols = self.sol
        for field in sols:
            if sols[field] is not None:
                self.__write_vec(field, sols[field])

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
        el_t = df.TensorElement(el_str, degree=deg, cell=self.cell)
        V_s = df.FunctionSpace(self.mesh, el_s)
        V_v = df.FunctionSpace(self.mesh, el_v)
        V_t = df.FunctionSpace(self.mesh, el_t)

        # Heat source
        f_heat = df.interpolate(self.heat_source, V_s)
        self.__write_xdmf("f_heat", f_heat, False)

        # Mass source
        f_mass = df.interpolate(self.mass_source, V_s)
        self.__write_xdmf("f_mass", f_mass, False)

        # Body force
        f_body = df.interpolate(self.body_force, V_v)
        self.__write_xdmf("f_body", f_body, False)

        # Heatflux force
        f_s = df.interpolate(self.f_s, V_v)
        self.__write_xdmf("f_s", f_s, False)

        # Stress force
        f_sigma = df.interpolate(self.f_sigma, V_t)
        self.__write_xdmf("f_sigma", f_sigma, False)

    def __write_discrete_system(self):
        r"""
        Write the discrete system matrix and the RHS vector.

        Also writes the ownership range of the matrix and the dofmap of the
        function space to the output folder.

        Can be used to analyze e.g. condition number.
        Include writing of :math:`\mathbf{A}` and :math:`\mathbf{b}`
        of :math:`\mathbf{A} \mathbf{x} = \mathbf{b}`.

        File-ending is `.mat`.

        Import the matrices/vectors e.g. into Matlab with:

        .. code-block:: matlab

            % Input into MATLAB
            A = table2array(readtable("A_0.mat","FileType","text"));
            b = table2array(readtable("b_0.mat","FileType","text"));

        Julia:

        .. code-block:: julia

            using DelimitedFiles
            A = readdlm("A_0.mat", ' ', Float64, '\n')

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
        ...     "tests/2d_heat/inputs/heat_01_coeffs_p1p1_stab.yml"
        ... ) # doctest: +ELLIPSIS
        Input:...
        >>> msh = H5Mesh("tests/2d_mesh/ring0.h5")
        >>> solver = Solver(params.dict, msh, "0") # "0" means time=0
        >>> solver.form_lhs = a
        >>> solver.form_rhs = L
        >>> solver.output_folder = "./"
        >>> solver._Solver__write_discrete_system()
        Write ./A_0_0.mat
        Write ./b_0_0.mat
        >>> print(open("A_0_0.mat","r").read())
        2.0000...00000e+00 -2.0000...00000e+00 0.000...000000e+00
        -2.0000...00000e+00 4.0000...00000e+00 -2.0000...00000e+00
        0.0000...00000e+00 -2.0000...00000e+00 2.000...000000e+00
        <BLANKLINE>
        >>> print(open("b_0_0.mat","r").read())
        2.500000000000000000e-01
        5.000000000000000000e-01
        2.500000000000000000e-01
        <BLANKLINE>
        """
        file_ending = ".mat"

        # Write left-hand side matrix
        A_name = (
            self.output_folder + "A_{}_{}".format(self.time, self.rank)
            + file_ending
        )
        print("Write {}".format(A_name))
        np.savetxt(A_name, df.assemble(self.form_lhs).array())

        # Write right-hand side vector
        b_name = (
            self.output_folder + "b_{}_{}".format(self.time, self.rank)
            + file_ending
        )
        print("Write {}".format(b_name))
        np.savetxt(b_name, df.assemble(self.form_rhs))

    def __write_mpi_information(self):
        r"""
        Writes the ownership range of the matrix and the dofmap of the
        function space to the output folder.
        """
        file_ending = ".mat"

        # Write ownership
        ownership_name = (
            self.output_folder + "ownership_{}_{}".format(self.time, self.rank)
            + file_ending
        )
        print("Write {}".format(ownership_name))
        np.savetxt(
            ownership_name,
            self.solver.ksp().getOperators()[0].getOwnershipRange(),
            fmt='%-i'
        )

        # Write dofmap
        if self.mode == "heat":
            w = self.mxd_fspaces["heat"]
            num_fields = 2
        elif self.mode == "stress":
            w = self.mxd_fspaces["stress"]
            num_fields = 3
        elif self.mode == "r13":
            w = self.mxd_fspaces["r13"]
            num_fields = 5
        dofmap_func = df.Function(w)
        for i in range(num_fields):
            dofmap_name = (
                self.output_folder + "dofmap_{}_{}".format(self.time, self.rank)
                + file_ending
            )

            # TODO: Remove that in the future
            # dofmap_func.split()[i].vector()[:] += self.rank
            # print(
            #     self.rank, " ", dofmap_func.split()[i].vector()[:],
            #     len(dofmap_func.split()[i].vector()[:])
            # )

            w.sub(i).dofmap().set(dofmap_func.split()[i].vector(), 1.0 * i)
        print("Write {}".format(dofmap_name))
        np.savetxt(dofmap_name, dofmap_func.vector(), fmt='%-i')

        # TODO: Remove that in the future
        # names = [
        #     "theta", "s", "p", "u", "sigma"
        # ]
        # tmp = dofmap_func.split()
        # print(tmp)
        # for i, field in enumerate(tmp):
        #     if field is not None:
        #         print(field.vector())
        #         print("DEBUG", names[i])
        #         self.__write_xdmf(names[i] + "_dist", field, False)

        for var in self.elems:
            print(
                "[", self.rank, "] ", str(var), " DOFs: ",
                len(self.fspaces[var].dofmap().dofs()), " (blocksize: ",
                self.fspaces[var].dofmap().block_size(), ")", sep=""
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

        *Problems*
            For MPi execution, XDMFFile.write function fails if the mesh
            is too small. I wasn't able to fix this in 2 days.
            Possible reasons:
            - Too much overlap in Schwarz-like decomposition
        """
        fname_xdmf = (
            self.output_folder + name + "_" + str(self.time) + ".xdmf"
        )
        with df.XDMFFile(self.mesh.mpi_comm(), fname_xdmf) as file:
            field.rename(name, name)
            print("Write {}".format(fname_xdmf))
            file.write(field, self.time)

        if write_pdf:
            import matplotlib.pyplot as plt  # pylint: disable=C0415
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
                    fieldname = name + "_" + str(indexMap[dimension][i + 1])
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

    def __write_vec(self, name, field):
        """
        Write a given field to a vector.
        """

        dimension = len(field.value_shape())
        if dimension == 0:
            # scalar
            fname_mat = (
                self.output_folder + name + "_" + str(self.time)
                + ".mat"
            )
            np.savetxt(
                self.output_folder + name + "_" + str(self.time) + ".mat",
                field.compute_vertex_values()
            )
            print("Write {}".format(fname_mat))
        else:
            # vector or tensor
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
                fieldname = name + "_" + str(indexMap[dimension][i + 1])
                fname_mat = (
                    self.output_folder + fieldname + "_" + str(self.time)
                    + ".mat"
                )
                np.savetxt(fname_mat, field.split()[i].compute_vertex_values())
                print("Write {}".format(fname_mat))
