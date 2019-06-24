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
import dolfin as d
import ufl as u
import mshr as m
import numpy as np
d.set_log_level(1000)  # 1: all logs
# **************************************************************************** #




# **************************************************************************** #
# SETTINGS
# **************************************************************************** #

# Problem parameters
system_ = "1" # 1=westerkamp2019, 2=coefficientless
tau_ = "0.1" # float
A0_ = "2" # int
A1_ = "0" # int
A2_ = "-1" # int
xi_tilde_ = "1.0" # float
theta_w_inner_ = "1.0" # float
theta_w_outer_ = "0.5" # float

# UFL vars
system = int(system_)
tau = d.Constant(float(tau_))
A0 = d.Constant(float(A0_))
A1 = d.Constant(float(A1_))
A2 = d.Constant(float(A2_))
xi_tilde = d.Constant(float(xi_tilde_))
theta_w_inner = d.Constant(float(theta_w_inner_))
theta_w_outer = d.Constant(float(theta_w_outer_))

# FEM parameters
deg_s = 1
deg_theta = 1
el_s = "Lagrange"
el_theta = "Lagrange"

# Convergence Study Parameters
max_exponent = 5

# Continous Interior Penalty (CIP) Stabilization with parameter delta_1:
stab_cip = True
delta_1 = d.Constant(1)

# Meshing parameters
use_gmsh = True

# Etc
save_matrix = False

# **************************************************************************** #



# **************************************************************************** #
# SETUP COMPUTATIONAL DOMAIN BY GENERATING MESH
# **************************************************************************** #
def create_mesh(p, plot_mesh_=False, overwrite_=False):
    """
    3000 = inner circle, 3100 = outer circle
    """

    if use_gmsh:
        gmsh_path = "/Applications/gmsh/Gmsh.app/Contents/MacOS/gmsh"
        geo_name = "ring"
        mesh_name = "{}{}".format(geo_name, p)

        mesh_avail = (
            os.path.isfile(mesh_name + ".xml") and
            os.path.isfile(mesh_name + "_facet_region.xml") and
            os.path.isfile(mesh_name + "_physical_region.xml")
        )

        if not mesh_avail or overwrite_:
            os.system(
                "{} -setnumber p {} -2 -o {}.msh {}.geo".format(
                    gmsh_path, p, mesh_name, geo_name))
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
# **************************************************************************** #



# **************************************************************************** #
# Setup function spaces and shape functions
# **************************************************************************** #
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
    f_str = "A0 + A2 * pow(R,2) + A1 * cos(phi)"
    f = d.Expression(f_str, degree=2, R=R, phi=phi, A0=A0, A1=A1, A2=A2)
    f_i = d.interpolate(f, v_scalar_)
    file_f = d.File("f.pvd")
    file_f.write(f_i)

    dev_grad_s = (0.5*(d.grad(s_)+u.transpose(d.grad(s_)))
                  - (1/3)*u.tr(d.grad(s_))*u.Identity(2))

    if system == 1:
        a1 = (
            + 12/5 * tau * d.inner(dev_grad_s, d.grad(r_))
            + 2/3 * (1/tau) * d.inner(s_, r_)
            - (5/2) * theta_ * d.div(r_)
        ) * d.dx + (
            + 5/(4*xi_tilde) * s_n * r_n
            + 11/10 * xi_tilde * s_t * r_t
        ) * d.ds
        a2 = - (d.div(s_) * kappa_) * d.dx
        l1 = (- 5.0/2.0 * r_n * theta_w_outer * d.ds(3100)
              - 5.0/2.0 * r_n * theta_w_inner * d.ds(3000))
        l2 = - (f * kappa_) * d.dx
    elif system == 2:
        a1 = (
            tau * d.inner(dev_grad_s, d.grad(r_))
            + (1/tau) * d.inner(s_, r_)
            - theta_ * d.div(r_)
        ) * d.dx + (
            + 1/(xi_tilde) * s_n * r_n
            + xi_tilde * s_t * r_t
        ) * d.ds
        a2 = - (d.div(s_) * kappa_) * d.dx
        l1 = (- r_n * theta_w_outer * d.ds(3100)
              - r_n * theta_w_inner * d.ds(3000))
        l2 = - (f * kappa_) * d.dx
    else:
        raise Exception('system={} is undefined'.format(system))

    if stab_cip:

        # 1)
        h_avg = mesh_.hmax()

        # 2)
        # TODO: Test this
        # h = d.CellDiameter(mesh_)
        # h_avg = (h('+') + h('-'))/2.0

        stab = - (delta_1 * h_avg**3 * d.jump(d.grad(theta_), n)
                  * d.jump(d.grad(kappa_), n)) * d.dS
    else:
        stab = 0

    a_ = a1 + a2 + stab
    l_ = l1 + l2

    if save_matrix:
        np.savetxt("A.mat", d.assemble(a_).array())
        ### Use in matrix with:
        ### >> T = readtable("a.txt");
        ### >> M=table2array(T);
        ### >> spy(M);
        ### >> cond(M);
        ### >> det(M);
        ### >> svd(M);

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
# **************************************************************************** #



# **************************************************************************** #
# CREATE EXACT SOLUTION
# => Exact solution and L_2/L_inf errors, high degree for good quadr.
# **************************************************************************** #
def get_exact_solution_heat():
    """
    s_e = (s_R, s_phi)
    TODO: Print parameters
    """

    data = open("exact_solutions_heat.csv")
    csv_dict = csv.DictReader(data, delimiter=",", quotechar="\"")
    for item in csv_dict:
        if (item.get("system") == system_ and
                item.get("tau") == tau_ and
                item.get("A0") == A0_ and
                item.get("A1") == A1_ and
                item.get("A2") == A2_ and
                item.get("theta_w_inner") == theta_w_inner_ and
                item.get("theta_w_outer") == theta_w_outer_):
            R = d.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2)
            phi = d.Expression("atan2(x[1],x[0])", degree=2)
            theta_e = d.Expression(item.get("theta_e"), degree=2, R=R, tau=tau)
            s_R = d.Expression(item.get("s_R"), degree=2, R=R, tau=tau)
            s_e = d.Expression(("s_R * cos(phi)", "s_R * sin(phi)"),
                               degree=2, phi=phi, s_R=s_R)
            return (s_e, theta_e)

    # no exact solution is csv file
    warnings.warn("No exact solution avail for given tau3")
    zero_dummy = d.Expression("0", degree=1)
    return ((zero_dummy, zero_dummy), zero_dummy)
# **************************************************************************** #


# **************************************************************************** #
# CALCULATE VARIOUS ERRORS BETWEEN NUMERICAL AND EXACT SOLUTION
# **************************************************************************** #
def calc_scalarfield_errors(theta_, theta_e_, v_theta_, name_, p_):
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
# **************************************************************************** #


# **************************************************************************** #
# CALCULATE VARIOUS ERRORS BETWEEN NUMERICAL AND EXACT SOLUTION
# **************************************************************************** #
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
    errs_f_L2 = [d.errornorm(
        field_e_i.split()[i], field_i.split()[i], 'L2', mesh=mesh_
    ) for i in range(dim)]
    errs_v_linf = [d.norm(
        field_e_i.split()[i].vector()-field_i.split()[i].vector(), 'linf'
    ) for i in range(dim)]
    print("L_2 error:", errs_f_L2)
    print("l_inf error:", errs_v_linf)

    field_e_i.rename("{}_e_i".format(name_), "{}_e_i".format(name_))
    file_field_e = d.File("{}_e.pvd".format(name_))
    file_field_e.write(field_e_i)

    field_i.rename("{}_i".format(name_), "{}_i".format(name_))
    file_field = d.File("{}_i.pvd".format(name_))
    file_field.write(field_i)

    return (errs_f_L2, errs_v_linf)
# **************************************************************************** #


# **************************************************************************** #
# CREATE ERROR/CONVERGENCE PLOT
# **************************************************************************** #
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

    data_sx, data_sy, data_theta = ({
        "p": [], "h": [], "L_2": [], "l_2": [], "l_inf": []} for _ in range(3))
    theta_array, theta_fspaces, theta_l2_change = ([] for _ in range(3))

    for p in range(max_exponent):
        (mesh, _, mesh_bounds) = create_mesh(p)
        (w, v_s, v_theta) = setup_function_spaces_heat(mesh, deg_s, deg_theta)
        (a, l) = setup_variational_form_heat(w, v_theta, mesh, mesh_bounds)
        (s, theta) = solve_variational_formulation_heat(a, l, w, [])
        (s_e, theta_e) = get_exact_solution_heat()

        theta_array.append(theta)
        theta_fspaces.append(v_theta)
        if not p == 0:
            d.Function.set_allow_extrapolation(theta_array[p], True)
            theta_l2_change.append(d.errornorm(
                d.project(theta_array[p], theta_fspaces[p-1]), theta_array[p-1]
            ))

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

        (err_f_l2_theta, err_v_linf_theta) = calc_scalarfield_errors(
            theta, theta_e, v_theta, "theta", p)
        data_theta["L_2"].append(err_f_l2_theta)
        data_theta["l_inf"].append(err_v_linf_theta)

    plot_errors(data_sx, "sx")
    plot_errors(data_sy, "sy")
    plot_errors(data_theta, "Theta")
    plot_single(data_sx["h"][:-1], theta_l2_change,
                "norm(theta_i-theta_{i-1})_L2", "theta change")
# **************************************************************************** #



# **************************************************************************** #
# MAIN
# **************************************************************************** #
if __name__ == '__main__':
    solve_system_heat()
# **************************************************************************** #
