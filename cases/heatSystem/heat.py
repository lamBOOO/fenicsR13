#!/usr/bin/env python3
# pylint: disable=invalid-name


# ------------------------------------------------------------------------------
# PYLINT SETTINGS
# ------------------------------------------------------------------------------
# pylint: disable=unsubscriptable-object
# pylint: disable=unused-import
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# TODOs
# ------------------------------------------------------------------------------
# - Export series of meshes
# - Fix the value for xi_tilde
# - Add edge scaling to CIP term
# ------------------------------------------------------------------------------


"""
Program to solve the decoupled (removed coupling term) heat system of the linearized R13 equations
"""


# ------------------------------------------------------------------------------
# IMPORTS
# ------------------------------------------------------------------------------
import os
import warnings
import matplotlib.pyplot as plt
import dolfin as d
import ufl as u
import mshr as m
import numpy as np
d.set_log_level(1000)  # 1: all logs
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------

# Problem parameters
tau = d.Constant(0.1)
A0 = d.Constant(2)
A1 = d.Constant(0)
A2 = d.Constant(-1)
xi_tilde = d.Constant(1.0)
theta_w_inner = d.Constant(1.0)
theta_w_outer = d.Constant(0.5)

# Define system:
# - 1: Westerkamp2019
# - 2: Coefficientless
# - 3: Westerkamp2014: Mixed poisson
system = 1

# FEM parameters
deg_s = 2
deg_theta = 1
el_s = "Lagrange"
el_theta = "Lagrange"

# Convergence Study Parameters
max_exponent = 5

# Continous Interior Penalty (CIP) Stabilization with parameter delta_1:
stab_cip = False
delta_1 = d.Constant(1)

# Meshing parameters
use_gmsh = True

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# SETUP COMPUTATIONAL DOMAIN BY GENERATING MESH
# ------------------------------------------------------------------------------
def create_mesh(p, plot_mesh_=False, overwrite_=False):
    """
    3000 = inner circle, 3100 = outer circle
    """

    if use_gmsh:
        gmsh_path = "/Applications/gmsh/Gmsh.app/Contents/MacOS/gmsh"
        geo_name = "ring"
        mesh_name = "{}{}".format(geo_name, p)

        mesh_avail = os.path.isfile(mesh_name + ".xml") and os.path.isfile(
            mesh_name + "_facet_region.xml") and os.path.isfile(mesh_name + "_physical_region.xml")

        if not mesh_avail or overwrite_:
            os.system(
                "{} -setnumber p {} -2 -o {}.msh {}.geo".format(gmsh_path, p, mesh_name, geo_name))
            os.system("dolfin-convert {0}.msh {0}.xml".format(mesh_name))

        mesh_ = d.Mesh("{}.xml".format(mesh_name))
        domain_markers_ = d.MeshFunction(
            'size_t', mesh_, "{}_physical_region.xml".format(mesh_name))
        boundary_markers_ = d.MeshFunction(
            'size_t', mesh_, "{}_facet_region.xml".format(mesh_name))

    else:
        raise Exception("FIXME: BCs are not working automatically here")

    if plot_mesh_:
        plt.figure()
        d.plot(mesh_, title="Mesh")
        plt.draw()

    print("Max edge length:", mesh_.hmax())

    return (mesh_, domain_markers_, boundary_markers_)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Setup function spaces and shape functions
# ------------------------------------------------------------------------------
def setup_function_spaces_heat(mesh_, deg_s_, deg_theta_):
    "TODO"
    cell = mesh_.ufl_cell()
    el_s_ = d.VectorElement(el_s, cell, degree=deg_s_)
    el_theta_ = d.FiniteElement(el_theta, cell, degree=deg_theta_)
    el_mxd_ = d.MixedElement([el_s_, el_theta_])
    v_s_ = d.FunctionSpace(mesh_, el_s_)
    v_theta_ = d.FunctionSpace(mesh_, el_theta_)
    w_ = d.FunctionSpace(mesh_, el_mxd_)
    return (w_, v_s_, v_theta_)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Setup problem
# ------------------------------------------------------------------------------
def setup_variational_formulation(w_, v_scalar_, mesh_, mesh_bounds_):
    """
    xi_tilde normally d.sqrt(2/d.pi), but we use 1 that looks right
    Note: Sign of a2 and l2 are correlated to sign of cip stabilization!!
    """

    # Define trial and testfunction
    (s_, theta_) = d.TrialFunctions(w_)
    (r_, kappa_) = d.TestFunctions(w_)

    # Define custom measeasure for boundaries
    d.ds = d.Measure('ds', domain=mesh_, subdomain_data=mesh_bounds_)
    d.dS = d.Measure('dS', domain=mesh_, subdomain_data=mesh_bounds_)

    # Normal and tangential components
    # => tangential (tx,ty) = (-ny,nx) only for 2D
    n = d.FacetNormal(mesh_)
    s_n = s_[0] * n[0] + s_[1] * n[1]
    r_n = r_[0] * n[0] + r_[1] * n[1]
    s_t = - s_[0] * n[1] + s_[1] * n[0]
    r_t = - r_[0] * n[1] + r_[1] * n[0]

    # Define source function
    R = d.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
    phi = d.Expression("atan2(x[1],x[0])", degree=2)
    f = d.Expression("A0 + A2 * pow(R,2) + A1 * cos(phi)", degree=2, R=R, phi=phi, A0=A0, A1=A1, A2=A2)
    f_i = d.interpolate(f, v_scalar_)
    file_f = d.File("f.pvd")
    file_f.write(f_i)

    # # No interpolation of f: No improvement
    # x = d.SpatialCoordinate(mesh_)
    # f = A0 + A2 * (x[0]**2 + x[1]**2)

    if system == 1:
        a1 = (
            # Attention: We actually solve a 3D problem!
            + 12/5 * tau * d.inner(0.5*(d.grad(s_)+u.transpose(d.grad(s_)))-(1/3)*u.tr(d.grad(s_))*u.Identity(2), d.grad(r_))
            # + 12/5 * tau * d.inner(d.dev(d.grad(s_)), d.grad(r_))
            + 2/3 * (1/tau) * d.inner(s_, r_)
            - (5/2) * theta_ * d.div(r_)
        ) * d.dx + (
            + 5/(4*xi_tilde) * s_n * r_n
            + 11/10 * xi_tilde * s_t * r_t
        ) * d.ds
        a2 = - (d.div(s_) * kappa_) * d.dx
        l1 = - 5.0/2.0 * r_n * theta_w_outer * d.ds(3100) - 5.0/2.0 * r_n * theta_w_inner * d.ds(3000)
        l2 = - (f * kappa_) * d.dx
    elif system == 2:
        a1 = (
            tau * d.inner(0.5*(d.grad(s_)+u.transpose(d.grad(s_)))-(1/3)*u.tr(d.grad(s_))*u.Identity(2), d.grad(r_))
            # tau * d.inner(d.dev(d.grad(s_)), d.grad(r_))
            + (1/tau) * d.inner(s_, r_)
            - theta_ * d.div(r_)
        ) * d.dx + (
            + 1/(xi_tilde) * s_n * r_n
            + xi_tilde * s_t * r_t
        ) * d.ds
        a2 = - (d.div(s_) * kappa_) * d.dx
        l1 = - r_n * theta_w_outer * d.ds(3100) - r_n * theta_w_inner * d.ds(3000)
        l2 = - (f * kappa_) * d.dx
    elif system == 3:
        a1 = (
            (1/tau) * d.inner(s_, r_)
            - theta_ * d.div(r_)
        ) * d.dx + (
            + 1/(xi_tilde) * s_n * r_n
        ) * d.ds
        a2 = - (d.div(s_) * kappa_) * d.dx
        l1 = - r_n * theta_w_outer * d.ds(3100) -  r_n * theta_w_inner * d.ds(3000)
        l2 = - (f * kappa_) * d.dx
    else:
        raise Exception('system={} is undefined'.format(system))

    if stab_cip:

        # 1)
        h_avg = mesh_.hmax()

        # 2)
        # h = d.CellDiameter(mesh_)
        # h_avg = (h('+') + h('-'))/2.0

        stab = - (delta_1 * h_avg**3 * d.jump(d.grad(theta_), n) * d.jump(d.grad(kappa_), n)) * d.dS
    else:
        stab = 0

    a_ = a1 + a2 + stab
    l_ = l1 + l2

    # Save matrix
    # np.savetxt("a.txt", d.assemble(a_).array())
    ### >> T = readtable("a.txt");
    ### >> M=table2array(T);
    ### >> spy(M);
    ### >> cond(M);
    ### >> det(M);
    ### >> svd(M);

    return (a_, l_)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# SOLVE THE SYSTEM
# - d.solve(a == l, sol, bcs): PETSc is default but RAM limited in conda
# ------------------------------------------------------------------------------
def solve_variational_formulation(a_, l_, w, bcs_, plot_=False):
    """
    Available solvers:
    solver_parameters={'linear_solver': 'gmres', 'preconditioner': 'ilu'}
    solver_parameters={'linear_solver': 'petsc', 'preconditioner': 'ilu'}
    solver_parameters={'linear_solver': 'direct'}
    solver_parameters={'linear_solver': 'mumps'}
    """

    sol_ = d.Function(w)
    d.solve(a_ == l_, sol_, bcs_, solver_parameters={'linear_solver': 'mumps'})

    (s_, theta_) = sol_.split()

    # Write files
    s_.rename('s', 's')
    file_s = d.File("s.pvd")
    file_s.write(s_)
    theta_.rename('theta', 'theta')
    file_theta = d.File("theta.pvd")
    file_theta.write(theta_)

    if plot_:
        plt.figure()
        d.plot(s_, title="s")
        plt.figure()
        d.plot(theta_, title="theta")
        plt.show()

    return (s_, theta_)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# CREATE EXACT SOLUTION
# => Exact solution and L_2/L_inf errors, high degree for good quadr.
# ------------------------------------------------------------------------------
def get_exact_solution(tau_):
    """
    s_e = (s_R, s_phi)
    TODO: Print parameters
    """

    if system == 1:

        R = d.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        phi = d.Expression("atan2(x[1],x[0])", degree=2)

        if tau_.values() == d.Constant(0.1).values():
            C_1 = d.Expression("-0.40855716127979214", degree=2)
            C_2 = d.Expression("2.4471587630476663", degree=2)
        elif tau_.values() == d.Constant(10.0).values():
            C_1 = d.Expression("-0.23487596350870713", degree=2)
            C_2 = d.Expression("13.827308558560057", degree=2)
        else:
            # raise Exception("No exact solution avail for given tau")
            warnings.warn("No exact solution avail for given tau")
            zero_dummy = d.Expression("0", degree=1)
            return ((zero_dummy, zero_dummy), zero_dummy)

        theta_e = d.Expression("C_2 + (- (20.0*C_1*std::log(R)) + ((5.0/4.0)*std::pow(R, 4)) - (2.0*std::pow(R, 2)*(24.0*std::pow(tau, 2) + 5.0)))/(tau*75.0)", degree=2, R=R, tau=tau, C_1=C_1, C_2=C_2)
        s_R = d.Expression(""" C_1/R + ( pow(R,2) - (pow(R,4)/4)) /R""",
                           degree=2, R=R, C_1=C_1)
        s_e = d.Expression((" s_R * cos(phi) ", " s_R * sin(phi)"), degree=2, phi=phi, s_R=s_R)

    elif system == 2:
        R = d.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        phi = d.Expression("atan2(x[1],x[0])", degree=2)
        if tau_.values() == d.Constant(100).values():
            if A0.values() == d.Constant(2).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(-1).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(0.5).values():
                print("22")
                theta_e = d.Expression("0.000625*(223878.98943374405 - 106674.66666666667*std::pow(R,2) + std::pow(R,4) + 3.972412068014572*std::log(R))", degree=2, R=R)
                s_R = d.Expression("-0.24827575425091075/R + R - std::pow(R,3)/4.", degree=2, R=R)
            else:
                warnings.warn("No exact solution avail for given tau=100")
                zero_dummy = d.Expression("0", degree=1)
                return ((zero_dummy, zero_dummy), zero_dummy)
        elif tau_.values() == d.Constant(0.01).values():
            if A0.values() == d.Constant(1).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(0).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(1.0).values():
                print("2")
                theta_e = d.Expression("77.26748929184684 - 0.25*std::pow(R,2) + 0.49360872123410277*std::log(R)", degree=2, R=R)
                s_R = d.Expression("0. - 0.49360872123410277/R + R/2.", degree=2, R=R)
            else:
                warnings.warn("No exact solution avail for given tau=0.01")
                zero_dummy = d.Expression("0", degree=1)
                return ((zero_dummy, zero_dummy), zero_dummy)
        elif tau_.values() == d.Constant(0.1).values():
            if A0.values() == d.Constant(0).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(0).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(0.5).values():
                theta_e = d.Expression("(-67 + 200*std::log(5) - 400*std::log(20) + 200*std::log(10*R))/ (5.*(-23 + 80*std::log(5) - 80*std::log(20)))", degree=2, R=R)
                s_R = d.Expression(""" -4/(R*(-23 + 80*std::log(5) - 80*std::log(20)))""", degree=2, R=R)
            elif A0.values() == d.Constant(2).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(-1).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(0.5).values():
                theta_e = d.Expression("(-76*std::pow(R,2))/15. + (5*std::pow(R,4))/8. + (11*(-36817 + 148480*std::log(5) - 24880*std::log(20)))/(1920.*(-23 + 80*std::log(5) - 80*std::log(20))) - (5665*std::log(10*R))/(8.*(-23 + 80*std::log(5) - 80*std::log(20)))", degree=2, R=R)
                s_R = d.Expression("""R - std::pow(R,3)/4. + 1133/(16.*R*(-23 + 80*std::log(5) - 80*std::log(20)))""", degree=2, R=R)
            elif A0.values() == d.Constant(2).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(0).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(0.5).values():
                theta_e = d.Expression("-(5273 - 21632*std::log(5) + 60*std::pow(R,2)*(-23 + 80*std::log(5) - 80*std::log(20)) + 1712*std::log(20) + 19920*std::log(10*R))/(12.*(-23 + 80*std::log(5) - 80*std::log(20)))", degree=2, R=R)
                s_R = d.Expression("""R + 166/(R*(-23 + 80*std::log(5) - 80*std::log(20)))""", degree=2, R=R)
            elif A0.values() == d.Constant(1).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(0).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(0.5).values():
                theta_e = d.Expression("7.036395447356292 - 2.5*std::pow(R,2) + 6.049130188983095*std::log(R)", degree=2, R=R)
                s_R = d.Expression("0. - 0.6049130188983095/R + R/2.", degree=2, R=R)
            elif A0.values() == d.Constant(1).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(0).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(1.0).values():
                print("1p")
                theta_e = d.Expression("9.73961676589525 - 0.25*std::pow(R,2) + 0.44589091940037223*std::log(R)", degree=2, R=R)
                s_R = d.Expression("0. - 0.44589091940037223/R + R/2.", degree=2, R=R)
            elif A0.values() == d.Constant(0).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(0).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(1.0).values():
                print("2")
                theta_e = d.Expression("1. + 3.230519938209405e-18*std::log(R)", degree=2, R=R)
                s_R = d.Expression("0. - 3.230519938209405e-18/R", degree=2, R=R)
            elif A0.values() == d.Constant(2).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(-1).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(1.0).values():
                print("2")
                theta_e = d.Expression("6.995305054386943 - 5.066666666666666*std::pow(R,2) + 0.625*std::pow(R,4) + 5.587043855102441*std::log(R)", degree=2, R=R)
                s_R = d.Expression("0. - 0.5587043855102442/R + R - std::pow(R,3)/4.", degree=2, R=R)
            else:
                warnings.warn("No exact solution avail for given tau1")
                zero_dummy = d.Expression("0", degree=1)
                return ((zero_dummy, zero_dummy), zero_dummy)
        elif tau_.values() == d.Constant(1).values():
            if A0.values() == d.Constant(2).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(-1).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(1.0).values():
                theta_e = d.Expression("(-896*std::pow(R,2) + 48*std::pow(R,4) + (5*(12449 + 2692*std::log(2)))/(25 + std::log(256)) + (7020*std::log(R))/(25 + std::log(256)))/768.", degree=2, R=R, A0=A0, A2=A2)
                s_R = d.Expression("R - std::pow(R,3)/4. - 585/(64.*R*(25 + std::log(256)))", degree=2, R=R)
            elif A0.values() == d.Constant(2).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(-0.5).values() and theta_w_inner.values() == d.Constant(1.0).values() and theta_w_outer.values() == d.Constant(1.0).values():
                theta_e = d.Expression(" + 3.9220404147903793 - 0.8333333333333333*std::pow(R,2) + 0.03125*std::pow(R,4) + 0.4360856152894074*std::log(R)", degree=2, R=R, A0=A0, A2=A2)
                s_R = d.Expression("-0.4360856152894074/R + R - 0.125*std::pow(R,3)", degree=2, R=R)
            else:
                warnings.warn("No exact solution avail for given tau=1")
                zero_dummy = d.Expression("0", degree=1)
                return ((zero_dummy, zero_dummy), zero_dummy)
        elif tau_.values() == d.Constant(10.0).values():
            if A0.values() == d.Constant(1).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(-1).values():
                theta_e = d.Expression("14.386999290419796 - 6.716666666666666*std::pow(R,2) + 0.00625*std::pow(R,4) + 0.023497664861756654*std::log(R)", degree=2, R=R)
                s_R = d.Expression("""- 0.23497664861756654/R + R - std::pow(R,3)/4.""", degree=2, R=R)
            elif A0.values() == d.Constant(2).values() and A1.values() == d.Constant(0).values() and A2.values() == d.Constant(-1).values():
                theta_e = d.Expression("14.386999290419796 - 6.716666666666666*std::pow(R,2) + 0.00625*std::pow(R,4) + 0.023497664861756654*std::log(R)", degree=2, R=R)
                s_R = d.Expression("""- 0.23497664861756654/R + R - std::pow(R,3)/4.""", degree=2, R=R)
            else:
                warnings.warn("No exact solution avail for given tau1")
                zero_dummy = d.Expression("0", degree=1)
                return ((zero_dummy, zero_dummy), zero_dummy)
        else:
            warnings.warn("No exact solution avail for given tau2")
            zero_dummy = d.Expression("0", degree=1)
            return ((zero_dummy, zero_dummy), zero_dummy)
        s_e = d.Expression(("s_R * cos(phi)", "s_R * sin(phi)"), degree=2, phi=phi, s_R=s_R)
    elif system == 3:
        R = d.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        phi = d.Expression("atan2(x[1],x[0])", degree=2)
        A0_expr = d.Expression(str(A0.values()[0]), degree=2)
        A1_expr = d.Expression(str(A1.values()[0]), degree=2)
        A2_expr = d.Expression(str(A2.values()[0]), degree=2)
        print(tau_.values())
        print(d.Constant(0.1).values())
        if tau_.values() == d.Constant(0.1).values():
            theta_e = d.Expression("(80*A0*(71 + 364*std::log(2)) + 5*A2*(1229 + 6148*std::log(2)) + 16*(24 + 240*std::log(2) + 3*A1*std::cos(1)*(1 + 20*std::log(2)) + 352*A1*std::cos(0.5)*(1 + std::log(32))))/(640.*(1 + std::log(256))) - (5*(4*A0*std::pow(R,2) + A2*std::pow(R,4) + 4*A1*std::pow(R,2)*std::cos(phi) - ((-64 + 1360*A0 + 1535*A2 + 1408*A1*std::cos(0.5) - 48*A1*std::cos(1))*std::log(R))/(20.*(1 + std::log(256)))))/8.", degree=2, R=R, phi=phi, A0=A0_expr, A1=A1_expr, A2=A2_expr)
            s_R = d.Expression("(8*A0*R + 4*A2*std::pow(R,3) + 8*A1*R*std::cos(phi) - (1360*A0 + 1535*A2 - 16*(4 - 88*A1*std::cos(0.5) + 3*A1*std::cos(1)))/(20.*R*(1 + std::log(256))))/16.", degree=2, R=R, phi=phi, A0=A0_expr, A1=A1_expr, A2=A2_expr)
            s_e = d.Expression(("s_R * cos(phi)", "s_R * sin(phi)"), degree=2, phi=phi, s_R=s_R)
        else:
            warnings.warn("No exact solution avail for given tau3")
            zero_dummy = d.Expression("0", degree=1)
            return ((zero_dummy, zero_dummy), zero_dummy)
    else:
        warnings.warn("No exact solution avail for system={}".format(system))

    return (s_e, theta_e)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# CALCULATE VARIOUS ERRORS BETWEEN NUMERICAL AND EXACT SOLUTION
# ------------------------------------------------------------------------------
def calc_field_errors(theta_, theta_e_, v_theta_, name_, p_):
    "TODO"

    # Theta
    field_e_i = d.interpolate(theta_e_, v_theta_)
    field_i = d.interpolate(theta_, v_theta_)

    difference = d.project(theta_e_ - theta_, v_theta_)
    difference.rename("difference_{}".format(name_),
                      "difference_{}".format(name_))
    file_difference = d.File("difference_{}_{}.pvd".format(name_, p_))
    file_difference.write(difference)

    err_f_L2 = d.errornorm(theta_e_, theta_, 'L2')
    err_v_linf = d.norm(field_e_i.vector()-field_i.vector(), 'linf')
    print("L_2 error:", err_f_L2)
    print("l_inf error:", err_v_linf)

    field_e_i.rename("{}_e_i".format(name_), "{}_e_i".format(name_))
    file_field_e = d.File("{}_e.pvd".format(name_))
    file_field_e.write(field_e_i)

    field_i.rename("{}_i".format(name_), "{}_i".format(name_))
    file_field = d.File("{}_i.pvd".format(name_))
    file_field.write(field_i)

    return (err_f_L2, err_v_linf)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# CALCULATE VARIOUS ERRORS BETWEEN NUMERICAL AND EXACT SOLUTION
# ------------------------------------------------------------------------------
def calc_vectorfield_errors(sol_, sol_e_, v_sol, mesh_, name_, p_):
    "TODO"

    # Vector values functions interpolated
    field_e_i = d.interpolate(sol_e_, v_sol)
    field_i = d.interpolate(sol_, v_sol)

    difference = d.project(sol_e_ - sol_, v_sol)
    difference.rename("difference_{}".format(name_),
                      "difference_{}".format(name_))
    file_difference = d.File("difference_{}_{}.pvd".format(name_, p_))
    file_difference.write(difference)

    dim = field_i.geometric_dimension()
    errs_f_L2 = [d.errornorm(field_e_i.split()[i], field_i.split()[i], 'L2', mesh=mesh_) for i in range(dim)]
    errs_v_linf = [d.norm(field_e_i.split()[i].vector()-field_i.split()[i].vector(), 'linf') for i in range(dim)]
    print("L_2 error:", errs_f_L2)
    print("l_inf error:", errs_v_linf)

    field_e_i.rename("{}_e_i".format(name_), "{}_e_i".format(name_))
    file_field_e = d.File("{}_e.pvd".format(name_))
    file_field_e.write(field_e_i)

    field_i.rename("{}_i".format(name_), "{}_i".format(name_))
    file_field = d.File("{}_i.pvd".format(name_))
    file_field.write(field_i)

    return (errs_f_L2, errs_v_linf)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# CREATE ERROR/CONVERGENCE PLOT
# ------------------------------------------------------------------------------
def plot_errors(data_, title_):
    "TODO"
    plt.loglog(data_["h"], data_["L_2"], "-o", label="L_2")
    plt.loglog(data_["h"], data_["l_inf"], "-o", label="l_inf")
    plt.loglog(data_["h"], np.array(
        2*np.power(data_["h"], 1)), "--", label="1st order")
    plt.loglog(data_["h"], np.array(
        0.02*np.power(data_["h"], 2)), "--", label="2nd order")
    plt.loglog(data_["h"], np.array(
        0.02*np.power(data_["h"], 3)), "--", label="3nd order")
    plt.xlabel("h_max")
    plt.ylabel("Error")
    plt.title(title_)
    plt.legend()
    plt.show()
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# PLOT SINGLE FIELD, I.E. CHANGE OF A FIELD
# ------------------------------------------------------------------------------
def plot_single(data_x_, data_y_, title_, legend_):
    "TODO"
    plt.loglog(data_x_, data_y_, "-o", label=legend_)
    plt.loglog(data_x_, np.array(
        2*np.power(data_x_, 1)), "--", label="1st order")
    plt.loglog(data_x_, np.array(
        0.02*np.power(data_x_, 2)), "--", label="2nd order")
    plt.xlabel("h_max")
    plt.ylabel("norm(theta_i-theta_{i-1})_L2")
    plt.title(title_)
    plt.legend()
    plt.show()
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# SOLVE DECOUPLED HEAT SYSTEM
# ------------------------------------------------------------------------------
def solve_heat_system():
    "TODO"

    max_exponent_ = max_exponent

    data_sx, data_sy, data_theta = ({
        "p": [], "h": [], "L_2": [], "l_2": [], "l_inf": []} for _ in range(3))
    theta_array, theta_fspaces, theta_l2_change = ([] for _ in range(3))

    for p in range(max_exponent_):
        (mesh, _, mesh_bounds) = create_mesh(p)
        (w, v_s, v_theta) = setup_function_spaces_heat(mesh, deg_s, deg_theta)
        (a, l) = setup_variational_formulation(w, v_theta, mesh, mesh_bounds)
        (s, theta) = solve_variational_formulation(a, l, w, [])
        (s_e, theta_e) = get_exact_solution(tau)

        theta_array.append(theta)
        theta_fspaces.append(v_theta)
        if not p == 0:
            d.Function.set_allow_extrapolation(theta_array[p], True)
            theta_l2_change.append(d.errornorm(d.project(theta_array[p], theta_fspaces[p-1]), theta_array[p-1]))

        data_sx["p"].append(p)
        data_sx["h"].append(mesh.hmax())
        data_sy["p"].append(p)
        data_sy["h"].append(mesh.hmax())
        data_theta["p"].append(p)
        data_theta["h"].append(mesh.hmax())

        (err_f_l2_s, err_v_linf_s) = calc_vectorfield_errors(
            s, s_e, v_s, mesh, "s", p)
        data_sx["L_2"].append(err_f_l2_s[0])
        data_sx["l_inf"].append(err_v_linf_s[0])
        data_sy["L_2"].append(err_f_l2_s[1])
        data_sy["l_inf"].append(err_v_linf_s[1])

        (err_f_l2_theta, err_v_linf_theta) = calc_field_errors(
            theta, theta_e, v_theta, "theta", p)
        data_theta["L_2"].append(err_f_l2_theta)
        data_theta["l_inf"].append(err_v_linf_theta)

    plot_errors(data_sx, "sx")
    plot_errors(data_sy, "sy")
    plot_errors(data_theta, "Theta")
    plot_single(data_sx["h"][:-1], theta_l2_change, "norm(theta_i-theta_{i-1})_L2", "theta change")
# ------------------------------------------------------------------------------


if __name__ == '__main__':
    solve_heat_system()
