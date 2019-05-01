"""
Program to solve the heat system
Analytical solution for Paraview approx:
(-20*(-0.40855)*ln(sqrt(coordsX^2+coordsY^2))+(5*(coordsX^2+coordsY^2)^2)/4
-2*(coordsX^2+coordsY^2)*(24*(1.0)^2+5))/(75*1.0) + 2.44715
"""
__author__ = "L. Theisen"
__copyright__ = "2019 %s" % __author__

# %%
import matplotlib.pyplot as plt
import dolfin as d
d.set_log_level(1)  # all logs

# %% Load mesh and setup spaces
MESH = d.Mesh("ring.xml")
SUBDOMAINS = d.MeshFunction('size_t', MESH, "ring_physical_region.xml")
BOUNDARIES = d.MeshFunction('size_t', MESH, "ring_facet_region.xml")

plt.figure()
d.plot(MESH, title="GMSH Mesh")
plt.draw()

# Setup function spaces and shape functions
P2 = d.VectorElement("Lagrange", MESH.ufl_cell(), 2)
P1 = d.FiniteElement("Lagrange", MESH.ufl_cell(), 1)
ME = P2 * P1  # mixed element space
W = d.FunctionSpace(MESH, ME)

# Boundary conditions
# FIXME: Insert the right BC
INNER_THETA = d.Constant(1.0)
OUTER_THETA = d.Constant(0.5)
INNER_S = d.Constant((0.0, 0.0))
OUTER_S = d.Constant((0.0, 0.0))

# Inner=3000
BC_INNER_THETA = d.DirichletBC(W.sub(1), INNER_THETA, BOUNDARIES, 3000)
BC_INNER_S = d.DirichletBC(W.sub(0), OUTER_S, BOUNDARIES, 3000)

# Outer=3100
BC_OUTER_THETA = d.DirichletBC(W.sub(1), OUTER_THETA, BOUNDARIES, 3100)
BC_OUTER_S = d.DirichletBC(W.sub(0), OUTER_S, BOUNDARIES, 3100)

BCS = [BC_INNER_THETA, BC_INNER_S, BC_OUTER_THETA, BC_OUTER_S]

# %% Setup problem definition
TAU = d.Constant(0.1)
XI_TILDE = d.Constant(0.1) # FIXME: Lookup right values
# THETA_W = d.Constant(0.5) # FIXME: Implement BC
DELTA_1 = d.Constant(0.1)

# Define trial and testfunction
(S, THETA) = d.TrialFunctions(W)
(R, KAPPA) = d.TestFunctions(W)

# Define source function
F = d.Expression("2 - pow(x[0],2) - pow(x[1],2)", degree=2)
# F = d.Expression("0", degree=0)

print("Setup variational formulation")

N = d.FacetNormal(MESH)
S_N = d.dot(S, N)
R_N = d.dot(R, N)
S_T = d.sqrt(d.dot(S, S) - S_N)
R_T = d.sqrt(d.dot(R, R) - R_N)

A1 = (
    + 12/5 * TAU * d.inner(d.dev(d.grad(S)), d.grad(R))
    + 2/3 * 1/TAU * d.dot(S, R)
    - 5/2 * THETA * d.div(R)
    ) * d.dx
A2 = - (d.div(S) * KAPPA) * d.dx
STAB = - (DELTA_1 * d.jump(d.grad(THETA), N) * d.jump(d.grad(KAPPA), N)) * d.dS
L1 = 0
L2 = - (F * KAPPA) * d.dx

A = A1 + A2 + STAB
L = L1 + L2

# FIXME: Think about how to treat boudary integrals:
#   div(S)*KAPPA) * dx + (5/(4*XI_TILDE) * S_N * R_N
#   + (11*XI_TILDE)/10 * S_T * R_T)*ds
# FIXME: Include jump term, see MS thesis from sweden: "jump"
# FIXME: L += (5/2 * THETA_W) * R_N * ds

# %%
print("Solving system...")
SOL = d.Function(W)
d.solve(A == L, SOL, BCS)
(S, THETA) = SOL.split()

print("Writing output files...")
S.rename('s', 's')
S_FILE_PVD = d.File("s.pvd")
S_FILE_PVD << S
THETA.rename('theta', 'theta')
THETA_FILE_PVD = d.File("theta.pvd")
THETA_FILE_PVD << THETA

print("Plotting solution...")
plt.figure()
d.plot(S, title="s")
plt.figure()
d.plot(THETA, title="Theta")
plt.show()

# TODO: Plot and implement analytical solution