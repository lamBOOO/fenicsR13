#!/usr/bin/env python3
# pylint: disable=invalid-name


# ------------------------------------------------------------------------------
# PYLINT SETTINGS
# ------------------------------------------------------------------------------
# none
# ------------------------------------------------------------------------------


"""
Program to solve the decoupled (removed coupling term) heat system of the linearized R13 equations
"""


# ------------------------------------------------------------------------------
# IMPORTS
# ------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import dolfin as d
d.set_log_level(1)  # all logs
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------

# Continous Interior Penalty (CIP) Stabilization with parameter delta_1:
stab_cip = True
delta_1 = d.Constant(1.0)

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# SETUP COMPUTATIONAL DOMAIN
# ------------------------------------------------------------------------------
mesh = d.Mesh("ring.xml") # TODO: Include internal mesher
mesh_regions = d.MeshFunction('size_t', mesh, "ring_physical_region.xml")
mesh_bounds = d.MeshFunction('size_t', mesh, "ring_facet_region.xml")

plt.figure()
d.plot(mesh, title="Mesh")
plt.draw()
print("Max edge length:", mesh.hmax())
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Setup function spaces and shape functions
# ------------------------------------------------------------------------------
el_p2 = d.VectorElement("Lagrange", mesh.ufl_cell(), 2)
el_p1 = d.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
el_mxd = d.MixedElement([el_p2, el_p1])
v_p2 = d.FunctionSpace(mesh, el_p2)
v_p1 = d.FunctionSpace(mesh, el_p1)
w = d.FunctionSpace(mesh, el_mxd)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Boundary conditions
# ------------------------------------------------------------------------------
# FIXME: Insert the right BC
theta_inner = d.Constant(1.0)
theta_outer = d.Constant(0.5)

# TODO: Check removing s BCs
s_inner = d.Constant((0.0, 0.0))
s_outer = d.Constant((0.0, 0.0))

# Inner=3000
bc_theta_inner = d.DirichletBC(w.sub(1), theta_inner, mesh_bounds, 3000)
bc_s_inner = d.DirichletBC(w.sub(0), s_inner, mesh_bounds, 3000)

# Outer=3100
bc_theta_outer = d.DirichletBC(w.sub(1), theta_outer, mesh_bounds, 3100)
bc_s_outer = d.DirichletBC(w.sub(0), s_outer, mesh_bounds, 3100)

bcs_theta = [bc_theta_inner, bc_theta_outer]
bcs_s = [bc_s_inner, bc_s_outer]
# bcs_full = bcs_theta + bcs_s
# bcs = bcs_theta
bcs = []
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Setup problem definition
# ------------------------------------------------------------------------------
tau = d.Constant(0.1)
xi_tilde = d.Constant(1.0) # FIXME: Lookup right value
theta_w_inner = d.Constant(1.0)
theta_w_outer = d.Constant(0.5)

# Define trial and testfunction
(s, theta) = d.TrialFunctions(w)
(r, kappa) = d.TestFunctions(w)

# Define source function
f = d.Expression("2 - pow(x[0],2) - pow(x[1],2)", degree=2)

print("Setup variational formulation")

# Normal and tangential components
# tangential (tx,ty) = (-ny,nx) only for 2D

# Define custom measeasure for boundaries
ds = d.Measure('ds', domain=mesh, subdomain_data=mesh_bounds)

n = d.FacetNormal(mesh)
s_n = s[0] * n[0] + s[1] * n[1]
r_n = r[0] * n[0] + r[1] * n[1]
s_t = - s[0] * n[1] + s[1] * n[0]
r_t = - r[0] * n[1] + r[1] * n[0]

a1 = (
    + 12/5 * tau * d.inner(d.dev(d.grad(s)), d.grad(r))
    + 2/3 * 1/tau * d.dot(s, r)
    - 5/2 * theta * d.div(r)
    ) * d.dx + (
        + 5/(4*xi_tilde) * s_n * r_n
        + 11/10 * xi_tilde * s_t * r_t
    ) * d.ds
a2 = - (d.div(s) * kappa) * d.dx
l1 = - 5/2 * r_n * theta_w_inner * ds(3000) - 5/2 * r_n * theta_w_outer * ds(3100)
l2 = - (f * kappa) * d.dx

if stab_cip:
    stab = - (delta_1 * d.jump(d.grad(theta), n)
              * d.jump(d.grad(kappa), n)) * d.dS
else:
    stab = 0

a = a1 + a2 + stab
l = l1 + l2

# FIXME: Add edge scaling to CIP term
# FIXME: Think about how to treat boudary integrals:
#   div(S)*KAPPA) * dx + (5/(4*XI_TILDE) * S_N * R_N
#   + (11*XI_TILDE)/10 * S_T * R_T)*ds
# FIXME: Include jump term, see MS thesis from sweden: "jump"
# FIXME: L += (5/2 * THETA_W) * R_N * ds
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# SOLVE THE SYSTEM
# ------------------------------------------------------------------------------
print("Solving system...")
sol = d.Function(w)
d.solve(a == l, sol, bcs)
(s, theta) = sol.split()

print("Writing output files...")

s.rename('s', 's')
file_s = d.File("s.pvd")
file_s.write(s)

theta.rename('theta', 'theta')
file_theta = d.File("theta.pvd")
file_theta.write(theta)

print("Plotting solution...")
plt.figure()
d.plot(s, title="s")
plt.figure()
d.plot(theta, title="theta")
plt.show()
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# COMPARE WITH EXACT SOLUTION
# ------------------------------------------------------------------------------
# Exact solution and L_2/L_inf errors, high degree for good quadr.
# X0 = d.FunctionSpace(mesh, "Lagrange", 1)
R = d.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=50)
C_1 = d.Expression("-0.40855716127979214", degree=5)
C_2 = d.Expression("2.4471587630476663", degree=5)
theta_e = d.Expression(""" C_2 + (1.0/75.0)*(-20*C_1*std::log(R) + (5.0/4.0)*std::pow(R, 4) - 2*std::pow(R, 2)*(24*std::pow(tau, 2) + 5))/tau """, degree=5, R=R, tau=tau, C_1=C_1, C_2=C_2)
theta_e = d.interpolate(theta_e, v_p1)
err_l2 = d.errornorm(theta_e, theta, 'L2')
print("L2 error:", err_l2)
plt.figure()
d.plot(theta_e, title="theta_e")
plt.show()

theta_e.rename('theta_e', 'theta_e')
file_theta_e = d.File("theta_e.pvd")
file_theta_e.write(theta_e)
# ------------------------------------------------------------------------------
