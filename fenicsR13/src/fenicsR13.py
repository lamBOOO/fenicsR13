#!/usr/bin/env python3


# **************************************************************************** #
# DOCUMENTATION
# **************************************************************************** #
"""
Program to solve the decoupled (removed coupling term) heat system of the
linearized R13 equations
"""
# **************************************************************************** #


# **************************************************************************** #
# SETTINGS
# **************************************************************************** #

# ------------------------------------------------------------------------ #
# PYLINT SETTINGS
# ------------------------------------------------------------------------ #
# pylint: disable=unsubscriptable-object
# pylint: disable=unused-import
# pylint: disable=invalid-name
# ------------------------------------------------------------------------ #

# ------------------------------------------------------------------------ #
# TODOs
# ------------------------------------------------------------------------ #
# - Export series of meshes
# ------------------------------------------------------------------------ #

# **************************************************************************** #


# **************************************************************************** #
# IMPORTS
# **************************************************************************** #
import os
import warnings
import csv
import matplotlib.pyplot as plt
import dolfin as df
import ufl
import mshr as m
import numpy as np

import meshes
from input import Input
# **************************************************************************** #


# **************************************************************************** #
# SETTINGS
# **************************************************************************** #

# Dolfin settings
df.set_log_level(1000)  # 1: all logs

# Problem parameters
system_ = "1"  # 1=westerkamp2019, 2=coefficientless
tau_ = "0.1"  # float
A0_ = "2"  # int
A1_ = "0"  # int
A2_ = "-1"  # int
xi_tilde_ = "1.0"  # float
theta_w_inner_ = "1.0"  # float
theta_w_outer_ = "0.5"  # float
v_t_inner_ = "10.0"  # float
v_t_outer_ = "0.0"  # float

# Model definitions
system = int(system_)
# solve_stress = False
# solve_heat = True
solve_heat = False
solve_stress = True

# UFL vars
tau = df.Constant(float(tau_))
A0 = df.Constant(float(A0_))
A1 = df.Constant(float(A1_))
A2 = df.Constant(float(A2_))
xi_tilde = df.Constant(float(xi_tilde_))
theta_w_inner = df.Constant(float(theta_w_inner_))
theta_w_outer = df.Constant(float(theta_w_outer_))
v_t_inner = df.Constant(float(v_t_inner_))
v_t_outer = df.Constant(float(v_t_outer_))

# FEM parameters
deg_s = 1
deg_theta = 1
deg_sigma = 1
deg_u = 1
deg_p = 1

# Convergence Study Parameters
max_exponent = 2

# Continous Interior Penalty (CIP) Stabilization with parameter delta_1:
stab_cip = True
delta_1 = df.Constant(1)
delta_2 = df.Constant(1)
delta_3 = df.Constant(0.01)

# Meshing parameters
use_gmsh = True

# Etc
save_matrix = False
plot_conv_rates = True
output_folder = "results/"

# Parallel MPI Settings
# -> parameters["ghost_mode"] = "shared_vertex"
# -> parameters["ghost_mode"] = "shared_facet"
df.parameters["ghost_mode"] = "shared_vertex"

# **************************************************************************** #


# **************************************************************************** #
# Setup function spaces and shape functions
# **************************************************************************** #
def setup_function_spaces_heat(mesh_, deg_theta_, deg_s_):
    "TODO"
    c = mesh_.ufl_cell()

    el_theta_ = df.FiniteElement("Lagrange", c, degree=deg_theta_)
    el_s_ = df.VectorElement("Lagrange", c, degree=deg_s_)
    el_mxd_ = df.MixedElement([el_theta_, el_s_])

    v_theta_ = df.FunctionSpace(mesh_, el_theta_)
    v_s_ = df.FunctionSpace(mesh_, el_s_)

    w_ = df.FunctionSpace(mesh_, el_mxd_)
    return (w_, v_theta_, v_s_)

def setup_function_spaces_stress(mesh_, deg_p_, deg_u_, deg_sigma_):
    "TODO"
    c = mesh_.ufl_cell()

    el_theta_ = df.FiniteElement("Lagrange", c, degree=deg_p_)
    el_s_ = df.VectorElement("Lagrange", c, degree=deg_u_)
    el_sigma_ = df.TensorElement("Lagrange", c, degree=deg_sigma_, symmetry=True)
    el_mxd_ = df.MixedElement([el_theta_, el_s_, el_sigma_])

    v_theta_ = df.FunctionSpace(mesh_, el_theta_)
    v_s_ = df.FunctionSpace(mesh_, el_s_)
    v_sigma_ = df.FunctionSpace(mesh_, el_sigma_)

    w_ = df.FunctionSpace(mesh_, el_mxd_)
    return (w_, v_theta_, v_s_, v_sigma_)
# **************************************************************************** #


# **************************************************************************** #
# Setup problem
# **************************************************************************** #
def setup_variational_form_heat(w_, v_scalar_, mesh_, mesh_bounds_):
    """
    xi_tilde normally d.sqrt(2/d.pi), but we use 1 that looks right
    Note: Sign of a2 and l2 are correlated to sign of cip stabilization!!
    Attention: We actually solve a 3D problem! Therefore adjust dev part
    """

    # Define trial and testfunction
    (theta_, s_) = df.TrialFunctions(w_)
    (kappa_, r_) = df.TestFunctions(w_)

    # Define custom measeasure for boundaries
    df.ds = df.Measure('ds', domain=mesh_, subdomain_data=mesh_bounds_)
    df.dS = df.Measure('dS', domain=mesh_, subdomain_data=mesh_bounds_)

    # Normal and tangential components
    # => tangential (tx,ty) = (-ny,nx) = perp(n) only for 2D
    n = df.FacetNormal(mesh_)
    t = ufl.perp(n)
    s_n = df.dot(s_, n)
    r_n = df.dot(r_, n)
    s_t = df.dot(s_, t)
    r_t = df.dot(r_, t)

    # Define source function
    R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
    phi = df.Expression("atan2(x[1],x[0])", degree=2)
    f_str = "A0 + A2 * pow(R,2) + A1 * cos(phi)"
    f = df.Expression(f_str, degree=2, R=R, phi=phi, A0=A0, A1=A1, A2=A2)
    f_i = df.interpolate(f, v_scalar_)
    f_i.rename('f_i', 'f_i')
    file_f = df.File(output_folder + "f.pvd")
    file_f.write(f_i)

    def dev3d(mat):
        "2d deviatoric part of actually 3d matrix"
        return (
            0.5 * (mat + ufl.transpose(mat))
            - (1/3) * ufl.tr(mat) * ufl.Identity(2)
        )

    if system == 1:
        a1 = (
            + 12/5 * tau * df.inner(dev3d(df.grad(s_)), df.grad(r_))
            + 2/3 * (1/tau) * df.inner(s_, r_)
            - (5/2) * theta_ * df.div(r_)
        ) * df.dx + (
            + 5/(4*xi_tilde) * s_n * r_n
            + 11/10 * xi_tilde * s_t * r_t
        ) * df.ds
        a2 = - (df.div(s_) * kappa_) * df.dx
        l1 = (- 5.0/2.0 * r_n * theta_w_outer * df.ds(3100)
              - 5.0/2.0 * r_n * theta_w_inner * df.ds(3000))
        l2 = - (f * kappa_) * df.dx
    elif system == 2:
        a1 = (
            tau * df.inner(dev3d(df.grad(s_)), df.grad(r_))
            + (1/tau) * df.inner(s_, r_)
            - theta_ * df.div(r_)
        ) * df.dx + (
            + 1/(xi_tilde) * s_n * r_n
            + xi_tilde * s_t * r_t
        ) * df.ds
        a2 = - (df.div(s_) * kappa_) * df.dx
        l1 = (-1 * r_n * theta_w_inner * df.ds(3000)
              -1 * r_n * theta_w_outer * df.ds(3100))
        l2 = - (f * kappa_) * df.dx
    else:
        raise Exception('system={} is undefined'.format(system))

    if stab_cip:
        # 1)
        # h_avg = mesh_.hmax()
        # 2)
        h = df.CellDiameter(mesh_)
        h_avg = (h('+') + h('-'))/2.0  # pylint: disable=not-callable
        stab = - (delta_1 * h_avg**3 * df.jump(df.grad(theta_), n)
                  * df.jump(df.grad(kappa_), n)) * df.dS
    else:
        stab = 0

    a_ = a1 + a2 + stab
    l_ = l1 + l2

    if save_matrix:
        np.savetxt("A.mat", df.assemble(a_).array())
        # Use in matrix with:
        # >> T = readtable("a.txt");
        # >> M=table2array(T);
        # >> spy(M);
        # >> cond(M);
        # >> det(M);
        # >> svd(M);

    return (a_, l_)

def setup_variational_form_stress(w_, v_scalar_, mesh_, mesh_bounds_):
    """
    xi_tilde normally d.sqrt(2/d.pi), but we use 1 that looks right
    Note: Sign of a2 and l2 are correlated to sign of cip stabilization!!
    Attention: We actually solve a 3D problem! Therefore adjust dev part
    """

    # Define trial and testfunction
    (p_, u_, sigma_) = df.TrialFunctions(w_)
    (q_, v_, psi_) = df.TestFunctions(w_)

    # Define custom measeasure for boundaries
    df.ds = df.Measure('ds', domain=mesh_, subdomain_data=mesh_bounds_)
    df.dS = df.Measure('dS', domain=mesh_, subdomain_data=mesh_bounds_)

    # Normal and tangential components
    # => tangential (tx,ty) = (-ny,nx) = perp(n) only for 2D
    n = df.FacetNormal(mesh_)
    t = ufl.perp(n)
    sigma_nn = df.dot(sigma_*n, n)
    psi_nn = df.dot(psi_*n, n)
    sigma_tt = df.dot(sigma_*t, t)
    psi_tt = df.dot(psi_*t, t)
    sigma_nt = df.dot(sigma_*n, t)
    psi_nt = df.dot(psi_*n, t)

    # Define source function
    R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
    phi = df.Expression("atan2(x[1],x[0])", degree=2)
    f_str = "2.0/5.0 * (1.0 - (5.0*std::pow(R,2))/(18.0*tau)) * std::cos(phi)"
    # f_str = "0"
    f = df.Expression(f_str, degree=2, R=R, phi=phi, A0=A0, A1=A1, A2=A2,
                      tau=tau)
    f_i = df.project(f, v_scalar_)
    f_i.rename('f', 'f')
    file_f = df.File(output_folder + "f.pvd")
    file_f.write(f_i)

    # d x d identiy matrix to use for Kronecker delta
    delta = df.Identity(2)

    def devOfGrad2(rank2):
        "From Henning's book p232"
        i, j, k, r = ufl.indices(4)
        entry_ijk = (
            (1/3) * (rank2[i, j].dx(k) + rank2[i, k].dx(j) + rank2[j, k].dx(i))
            # ufl.sym(ufl.grad(rank2))
            - (1/15) * (
                + (2 * rank2[i, r].dx(r) + rank2[r, r].dx(i)) * delta[j, k]
                + (2 * rank2[j, r].dx(r) + rank2[r, r].dx(j)) * delta[i, k]
                + (2 * rank2[k, r].dx(r) + rank2[r, r].dx(k)) * delta[i, j]
            )
        )
        tensor = ufl.as_tensor(entry_ijk, (i, j, k))
        return tensor

    if system == 1:
        a1 = (
            + 2 * tau * df.inner(devOfGrad2(sigma_), df.grad(psi_))
            # + 2 * tau * d.inner(ufl.tr(d.grad(sigma_)), d.grad(psi_))
            + (1/tau) * df.inner(sigma_, psi_)
            - 2 * df.dot(u_, df.div(df.sym(psi_)))
        ) * df.dx + (
            + 21/(10*xi_tilde) * sigma_nn * psi_nn
            + 2 * xi_tilde * (
                (sigma_tt + (1/2) * sigma_nn) * (psi_tt + (1/2) * psi_nn)
            )
            + (2/xi_tilde) * sigma_nt * psi_nt
        ) * df.ds
        # l1 = (- 2 * v_t_inner * df.sin(phi) * psi_nt * df.ds(3000)
        #       - 2 * v_t_outer * df.sin(phi) * psi_nt * df.ds(3100))
        l1 = (- 2 * v_t_inner * psi_nt * df.ds(3000)
              - 2 * v_t_outer * psi_nt * df.ds(3100))
        a2 = (df.dot(df.div(sigma_), v_) + df.dot(df.grad(p_), v_)) * df.dx
        l2 = df.Constant(0) * df.div(v_) * df.dx # dummy
        a3 = df.dot(u_, df.grad(q_)) * df.dx
        l3 = - (f * q_) * df.dx
    elif system == 2:
        raise Exception("not avail")
    else:
        raise Exception('system={} is undefined'.format(system))

    if stab_cip:
        # 1)
        # h_avg = mesh_.hmax()
        # 2)
        h = df.CellDiameter(mesh_)
        h_avg = (h('+') + h('-'))/2.0  # pylint: disable=not-callable
        stab = (
            + delta_2 * h_avg**3
            * df.dot(df.jump(df.grad(u_), n), df.jump(df.grad(v_), n))
            - delta_3 * h_avg *
            df.jump(df.grad(p_), n) * df.jump(df.grad(q_), n)
            ) * df.dS
    else:
        stab = 0

    a_ = a1 + a2 + a3 + stab
    l_ = l1 + l2 + l3

    if save_matrix:
        np.savetxt("A.mat", df.assemble(a_).array())
        # Use in matrix with:
        # >> T = readtable("a.txt");
        # >> M=table2array(T);
        # >> spy(M);
        # >> cond(M);
        # >> det(M);
        # >> svd(M);

    return (a_, l_)
# **************************************************************************** #


# **************************************************************************** #
# SOLVE THE SYSTEM
# **************************************************************************** #
# - d.solve(a == l, sol, bcs): PETSc is default but RAM limited in conda
# **************************************************************************** #
def solve_variational_formulation_heat(a_, l_, w, bcs_, plot_=False):
    """
    Available solvers:
    solver_parameters={'linear_solver': 'gmres', 'preconditioner': 'ilu'}
    solver_parameters={'linear_solver': 'petsc', 'preconditioner': 'ilu'}
    solver_parameters={'linear_solver': 'direct'}
    solver_parameters={'linear_solver': 'mumps'}
    """

    sol_ = df.Function(w)
    df.solve(a_ == l_, sol_, bcs_, solver_parameters={'linear_solver': 'mumps'})

    (theta_, s_) = sol_.split()

    # Write files
    s_.rename('s', 's')
    file_s = df.File(output_folder + "s.pvd")
    file_s.write(s_)
    theta_.rename('theta', 'theta')
    file_theta = df.File(output_folder + "theta.pvd")
    file_theta.write(theta_)

    if plot_:
        plt.figure()
        df.plot(s_, title="s")
        plt.figure()
        df.plot(theta_, title="theta")
        plt.show()

    return (s_, theta_)

def solve_variational_formulation_stress(a_, l_, w, bcs_, plot_=False):
    """
    Available solvers:
    solver_parameters={'linear_solver': 'gmres', 'preconditioner': 'ilu'}
    solver_parameters={'linear_solver': 'petsc', 'preconditioner': 'ilu'}
    solver_parameters={'linear_solver': 'direct'}
    solver_parameters={'linear_solver': 'mumps'}
    """

    sol_ = df.Function(w)
    df.solve(a_ == l_, sol_, bcs_, solver_parameters={'linear_solver': 'mumps'})

    (p_, u_, sigma_) = sol_.split()

    # Write files

    # symmetric tensor has three comps
    sigma_.rename('sigma_xx', 'sigma_xx')
    file_sigma_xx = df.File(output_folder + "sigma_xx.pvd")
    file_sigma_xx.write(sigma_.split()[0])
    sigma_.rename('sigma_xy', 'sigma_xy')
    file_sigma_xy = df.File(output_folder + "sigma_xy.pvd")
    file_sigma_xy.write(sigma_.split()[1])
    sigma_.rename('sigma_yy', 'sigma_yy')
    file_sigma_yy = df.File(output_folder + "sigma_yy.pvd")
    file_sigma_yy.write(sigma_.split()[2])

    u_.rename('u', 'u')
    file_u = df.File(output_folder + "u.pvd")
    file_u.write(u_)

    p_.rename('p', 'p')
    file_p = df.File(output_folder + "p.pvd")
    file_p.write(p_)

    if plot_:

        plt.figure()
        df.plot(p_, title="p")

        plt.figure()
        df.plot(u_, title="u")


        plt.figure()
        df.plot(sigma_, title="sigma")

        plt.show()

    return (p_, u_, sigma_)
# **************************************************************************** #


# **************************************************************************** #
# CREATE EXACT SOLUTION
# => Exact solution and L_2/L_inf errors, high degree for good quadr.
# **************************************************************************** #
def get_exact_solution_heat():
    """
    s_e = (s_R, s_phi)
    """

    data = open("exact_solutions_heat.csv")
    csv_dict = csv.DictReader(data, delimiter=",", quotechar="\"")
    for item in csv_dict:
        entry_available = (
            item.get("system") == system_ and
            item.get("tau") == tau_ and
            item.get("A0") == A0_ and
            item.get("A1") == A1_ and
            item.get("A2") == A2_ and
            item.get("theta_w_inner") == theta_w_inner_ and
            item.get("theta_w_outer") == theta_w_outer_
        )
        if entry_available:
            R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
            phi = df.Expression("atan2(x[1],x[0])", degree=2)
            theta_e = df.Expression(item.get("theta_e"), degree=2, R=R, tau=tau)
            s_R = df.Expression(item.get("s_R"), degree=2, R=R, tau=tau)
            s_e = df.Expression(("s_R * cos(phi)", "s_R * sin(phi)"),
                                degree=2, phi=phi, s_R=s_R)
            return (s_e, theta_e)

    # no exact solution is csv file
    warnings.warn("No exact solution avail for given tau3")
    scalar_dummy = df.Expression("0", degree=1)
    vector_dummy = df.Expression(("0", "0"), degree=1)
    return (vector_dummy, scalar_dummy)

def get_exact_solution_stress(space0, space1, space2, mesh_):
    ".."

    if system_ == "1":

        with open('exact_solution.cpp', 'r') as file:
            exact_solution_cpp_code = file.read()

        exact_solution = df.compile_cpp_code(exact_solution_cpp_code)

        p_e = df.CompiledExpression(exact_solution.Pressure(), degree=1)
        p_e_i = df.interpolate(p_e, space0)
        p_e_i.rename("p_e", "p_e")
        file_p = df.File(output_folder + "p_e.pvd")
        file_p.write(p_e_i)

        u_e = df.CompiledExpression(exact_solution.Velocity(), degree=1)
        u_e_i = df.interpolate(u_e, space1)
        u_e_i.rename("u_e", "u_e")
        file_u = df.File(output_folder + "u_e.pvd")
        file_u.write(u_e_i)

        sigma_e = df.CompiledExpression(exact_solution.Stress(), degree=1)
        sigma_e_i = df.interpolate(sigma_e, space2)
        sigma_e_i.rename("sigma_e", "sigma_e")
        # file_sigma = df.File(output_folder + "sigma_e.pvd")
        # file_sigma.write(sigma_e_i)

        return (p_e, u_e, sigma_e)




        # ##### OLD, NOT WORKING WELL
        # R = df.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
        # phi = df.Expression("atan2(x[1],x[0])", degree=5)

        # x = df.Expression("x[0]", degree=1)
        # y = df.Expression("x[1]", degree=1)

        # C_0 = df.Constant(-50.80230139855979)
        # C_1 = df.Constant(0.6015037593984962)
        # C_2 = df.Constant(0)
        # C_3 = df.Constant(-444.7738727200452)
        # C_4 = df.Constant(-0.12443443849461801)
        # C_5 = df.Constant(39.38867688999618)
        # C_6 = df.Constant(-0.6917293233082705)
        # C_7 = df.Constant(0)
        # C_8 = df.Constant(0)
        # C_9 = df.Constant(0)
        # C_10 = df.Constant(0)
        # C_11 = df.Constant(0)
        # C_12 = df.Constant(2.255312046238658E-11)
        # C_13 = df.Constant(407.2248457002586)
        # C_14 = df.Constant(-104.89346597195336)
        # C_15 = df.Constant(4.870715709115059E-7)

        # lambda_1 = df.Constant(df.sqrt(5/9))
        # lambda_2 = df.Constant(df.sqrt(5/6))
        # lambda_3 = df.Constant(df.sqrt(3/2))

        # gamma_0 = 1
        # gamma = 1

        # A_1 = df.Constant(0.4)

        # # # x = df.SpatialCoordinate(mesh_)
        # # class PressureExact(df.UserExpression):
        # #     def eval(self, values, x):
        # #         C_0 = df.Constant(-50.80230139855979)
        # #         values[0] = ufl.atan_2(x[1], x[0])*C_0
        # # p_e_i = df.interpolate(PressureExact(), space)
        # # # p_e_i = df.interpolate(df.Expression("2*bessel", degree=2, bessel=PressureExact()), space)
        # # p_e_i.rename("p_e", "p_e")
        # # file_p = df.File(output_folder + "p_e.pvd")
        # # file_p.write(p_e_i)

        # d_0 = C_9 + C_2*ufl.bessel_K(0, R*lambda_2/tau) + C_8*ufl.bessel_I(0, R*lambda_2/tau)
        # d = - (10*A_1 * R**2)/(27*tau) + (4*C_4*R)/(tau) - (2*C_5*tau)/R + C_14*ufl.bessel_K(1, R*lambda_2/tau) + C_15* ufl.bessel_I(1, R*lambda_2/tau)
        # p = d_0 + ufl.cos(phi) * d
        # # test = ufl.bessel_K(0,1)
        # # p = df.Expression("d_0 + cos(phi) * d", degree=2, R=R, phi=phi, d_0=d_0, d=d)

        # # p_e_i = df.project(p, space)

        # # p_e_i = df.interpolate(phi, space)
        # # p_e_i = df.project(test, space)
        # # p_e_i = df.project(phi, space)

        # x = df.SpatialCoordinate(mesh_)
        # # p_e_i = df.project(ufl.atan_2(x[1], x[0]), space)
    else:
        warnings.warn("No exact solution avail")
        scalar_dummy = df.Expression("0", degree=1)
        vector_dummy = df.Expression(("0", "0"), degree=1)
        return (vector_dummy, scalar_dummy)
# **************************************************************************** #


# **************************************************************************** #
# CALCULATE VARIOUS ERRORS BETWEEN NUMERICAL AND EXACT SOLUTION
# **************************************************************************** #
def calc_scalarfield_errors(theta_, theta_e_, v_theta_, name_, p_):
    "TODO"

    of = output_folder

    # Theta
    field_e_i = df.interpolate(theta_e_, v_theta_)
    field_i = df.interpolate(theta_, v_theta_)

    difference = df.project(theta_e_ - theta_, v_theta_)
    difference.rename("difference_{}".format(name_),
                      "difference_{}".format(name_))
    file_difference = df.File(of + "difference_{}_{}.pvd".format(name_, p_))
    file_difference.write(difference)

    err_f_L2 = df.errornorm(theta_e_, theta_, 'L2')
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
# **************************************************************************** #


# **************************************************************************** #
# CALCULATE VARIOUS ERRORS BETWEEN NUMERICAL AND EXACT SOLUTION
# **************************************************************************** #
def calc_vectorfield_errors(sol_, sol_e_, v_sol, mesh_, name_, p_):
    "TODO"

    of = output_folder

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
        field_e_i.split()[i], field_i.split()[i], 'L2', mesh=mesh_
    ) for i in range(dim)]
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
# **************************************************************************** #


# **************************************************************************** #
# CREATE ERROR/CONVERGENCE PLOT
# **************************************************************************** #
def plot_errrorsNew(data_):
    "TODO"

    h = [df["h"] for df in data_]
    for key in data_[0]:
        if key != "h":

            # LaTeX text fonts:
            # Use with raw strings: r"$\mathcal{O}(h^1)$"
            # plt.rc('text', usetex=True)
            # plt.rc('font', family='serif')

            field = [df[key] for df in data_]
            for etype in field[0]:
                values = [df[etype] for df in field]
                plt.loglog(h, values, "-o", label=etype)
            plt.loglog(h, np.array(2*np.power(h, 1)), "--", label="O(h^1)")
            plt.loglog(h, np.array(0.02*np.power(h, 2)), "--", label="O(h^2)")
            plt.loglog(h, np.array(0.02*np.power(h, 3)), "--", label="O(h^3)")
            plt.xlabel("max(h)")
            plt.ylabel("Error")
            plt.title(key)
            plt.legend()
            plt.show()
# **************************************************************************** #


# **************************************************************************** #
# PLOT SINGLE FIELD, I.E. CHANGE OF A FIELD
# **************************************************************************** #
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
# **************************************************************************** #





# **************************************************************************** #
# SOLVE DECOUPLED HEAT SYSTEM
# **************************************************************************** #
def solve_system_heat():
    "TODO"

    settings = Input("input.yml").dict
    mesh_names = settings["meshes"]

    data = []

    for p, mesh_name in enumerate(mesh_names):

        mesh_name = mesh_names[p]

        # setup and solve problem
        current_mesh = meshes.H5Mesh(mesh_name)
        (w, v_theta, v_s) = setup_function_spaces_heat(current_mesh.mesh, deg_theta, deg_s)
        (a, l) = setup_variational_form_heat(w, v_theta, current_mesh.mesh, current_mesh.boundaries)
        (s, theta) = solve_variational_formulation_heat(a, l, w, [])
        (s_e, theta_e) = get_exact_solution_heat()

        # calc errors
        (f_l2_s, v_linf_s) = calc_vectorfield_errors(
            s, s_e, v_s, current_mesh.mesh, "s", p)
        (f_l2_theta, v_linf_theta) = calc_scalarfield_errors(
            theta, theta_e, v_theta, "theta", p)

        # store errors
        data.append({
            "h": current_mesh.mesh.hmax(),
            "theta": {"L_2": f_l2_theta, "l_inf": v_linf_theta},
            "sx": {"L_2": f_l2_s[0], "l_inf": v_linf_s[0]},
            "sy": {"L_2": f_l2_s[1], "l_inf": v_linf_s[1]}
        })

    # plot errors
    if plot_conv_rates:
        plot_errrorsNew(data)
# ------------------------------------------------------------------------------
def solve_system_stress():
    "TODO"

    settings = Input("input.yml").dict
    mesh_names = settings["meshes"]

    data = []

    for p, mesh_name in enumerate(mesh_names):

        mesh_name = mesh_names[p]

        # setup and solve problem
        current_mesh = meshes.H5Mesh(mesh_name)
        (w, v_p, v_u, v_sigma) = setup_function_spaces_stress(
            current_mesh.mesh, deg_p, deg_u, deg_sigma
        )
        (a, l) = setup_variational_form_stress(w, v_p, current_mesh.mesh, current_mesh.boundaries)
        (p, u, sigma) = solve_variational_formulation_stress(a, l, w, [])
        (p_e, u_e, sigma_e) = get_exact_solution_stress(v_p, v_u, v_sigma, current_mesh.mesh)

    #     # calc errors
        (f_l2_u, v_linf_u) = calc_vectorfield_errors(
            u, u_e, v_u, current_mesh.mesh, "u", p)
        (f_l2_p, v_linf_p) = calc_scalarfield_errors(
            p, p_e, v_p, "p", p)

    #     # store errors
        data.append({
            "h": current_mesh.mesh.hmax(),
            "p": {"L_2": f_l2_p, "l_inf": v_linf_p},
            "ux": {"L_2": f_l2_u[0], "l_inf": v_linf_u[0]},
            "uy": {"L_2": f_l2_u[1], "l_inf": v_linf_u[1]}
        })

    # plot errors
    if plot_conv_rates:
        plot_errrorsNew(data)
# **************************************************************************** #



# **************************************************************************** #
# MAIN
# **************************************************************************** #
if __name__ == '__main__':
    if solve_heat:
        solve_system_heat()
    if solve_stress:
        solve_system_stress()
# **************************************************************************** #
