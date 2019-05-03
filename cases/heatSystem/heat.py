#!/usr/bin/env python3
# pylint: disable=invalid-name

"""
Program to solve the heat system
"""
__author__ = "L. Theisen"
__copyright__ = "2019 %s" % __author__

import matplotlib.pyplot as plt
import dolfin as d
d.set_log_level(1)  # all logs

# Load mesh and setup spaces
mesh = d.Mesh("ring.xml") # TODO: Include internal mesher
mesh_regions = d.MeshFunction('size_t', mesh, "ring_physical_region.xml")
mesh_bounds = d.MeshFunction('size_t', mesh, "ring_facet_region.xml")

plt.figure()
d.plot(mesh, title="Mesh")
plt.draw()
print("Max edge length:", mesh.hmax())

# Setup function spaces and shape functions
el_p2 = d.VectorElement("Lagrange", mesh.ufl_cell(), 2)
el_p1 = d.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
el_mxd = d.MixedElement([el_p2, el_p1])
v_p2 = d.FunctionSpace(mesh, el_p2)
v_p1 = d.FunctionSpace(mesh, el_p1)
w = d.FunctionSpace(mesh, el_mxd)

# Boundary conditions
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
bcs = bcs_theta

# Setup problem definition
tau = d.Constant(0.1)
xi_tilde = d.Constant(0.1) # FIXME: Lookup right values
# THETA_W = d.Constant(0.5) # FIXME: Implement BC
delta_1 = d.Constant(1.0)

# Define trial and testfunction
(s, theta) = d.TrialFunctions(w)
(r, kappa) = d.TestFunctions(w)

# Define source function
f = d.Expression("2 - pow(x[0],2) - pow(x[1],2)", degree=2)
# F = d.Expression("0", degree=0)

print("Setup variational formulation")

n = d.FacetNormal(mesh)
s_n = d.dot(s, n)
r_n = d.dot(r, n)
s_t = d.sqrt(d.dot(s, s) - s_n)
r_t = d.sqrt(d.dot(r, r) - r_n)

a1 = (
    + 12/5 * tau * d.inner(d.dev(d.grad(s)), d.grad(r))
    + 2/3 * 1/tau * d.dot(s, r)
    - 5/2 * theta * d.div(r)
    ) * d.dx
a2 = - (d.div(s) * kappa) * d.dx
stab_cip = - (delta_1 * d.jump(d.grad(theta), n)
              * d.jump(d.grad(kappa), n)) * d.dS
l1 = 0
l2 = - (f * kappa) * d.dx

a = a1 + a2 + stab_cip
l = l1 + l2

# FIXME: Add edge scaling to CIP term
# FIXME: Think about how to treat boudary integrals:
#   div(S)*KAPPA) * dx + (5/(4*XI_TILDE) * S_N * R_N
#   + (11*XI_TILDE)/10 * S_T * R_T)*ds
# FIXME: Include jump term, see MS thesis from sweden: "jump"
# FIXME: L += (5/2 * THETA_W) * R_N * ds

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

# Exact solution and L_2/L_inf errors, high degree for good quadr.
# X0 = d.FunctionSpace(mesh, "Lagrange", 1)
r = d.Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=5)
c1 = d.Expression("-0.40855716127979214", degree=5)
c2 = d.Expression("2.4471587630476663", degree=5)
theta_e = d.Expression("""
    (-20*c1*log2(r) + (5*pow(r,4))/4 - 2*pow(r,2)*(24*pow(tau,2)+5))/(75*tau)
    + c2 """, degree=5, r=r, tau=tau, c1=c1, c2=c2)
theta_e = d.interpolate(theta_e, v_p1)
err_l2 = d.errornorm(theta_e, theta, 'L2')
print("L2 error:", err_l2)
plt.figure()
d.plot(theta_e, title="theta_e")
plt.show()

theta_e.rename('theta_e', 'theta_e')
file_theta_e = d.File("theta_e.pvd")
file_theta_e.write(theta_e)
