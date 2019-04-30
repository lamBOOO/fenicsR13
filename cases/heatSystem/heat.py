"""
Program to solve the heat system
"""
__author__ = "L. Theisen"
__copyright__ = "2019 %s" % __author__

# %%
# from mshr import Circle, generate_mesh
from dolfin import plot, set_log_level, FiniteElement, VectorElement, \
    FunctionSpace, Mesh, MeshFunction, TrialFunctions, TestFunctions, Expression, \
    dx, dot, grad, dev, div, ds, FacetNormal, inner, DirichletBC, Constant, sqrt, solve, Function
import dolfin
import matplotlib.pyplot as plt

# Settings
set_log_level(1)  # all logs

# %%
# Load gmsh mesh
MESH = Mesh("ring.xml")
SUBDOMAINS = MeshFunction('size_t', MESH, "ring_physical_region.xml")
BOUNDARIES = MeshFunction('size_t', MESH, "ring_facet_region.xml")

plt.figure()
plot(MESH, title="GMSH Mesh")
plt.show()

# Setup function spaces and shape functions
P2 = VectorElement("Lagrange", MESH.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", MESH.ufl_cell(), 1)
ME = P2 * P1  # mixed element space
W = FunctionSpace(MESH, ME)

# Boundary conditions # TODO: Insert the right BC
INNER_THETA = 1.0
OUTER_THETA = 0.5
# Inner=3000
BC_INNER_THETA = DirichletBC(W.sub(1), Constant(INNER_THETA), BOUNDARIES, 3000)
BC_INNER_S = DirichletBC(W.sub(0), Constant((0.0, 0.0)), BOUNDARIES, 3000)

# Outer=3100
BC_OUTER_THETA = DirichletBC(W.sub(1), Constant(OUTER_THETA), BOUNDARIES, 3100)
BC_OUTER_S = DirichletBC(W.sub(0), Constant((0.0, 0.0)), BOUNDARIES, 3100)

BCS = [BC_INNER_THETA, BC_INNER_S, BC_OUTER_THETA, BC_OUTER_S]

# %%
# Parameters
TAU = 0.1
XI_TILDE = 0.5  # TODO: What's the value for this?, check against normal xi
THETA_W = 0.5  # TODO: Implement BC

# Define trial and testfunction
(S, THETA) = TrialFunctions(W)
(R, KAPPA) = TestFunctions(W)

# Define source function
F = Expression("2 - pow(x[0],2) - pow(x[1],2)", degree=2)  # f=2-x^2-y^2

print("Setup variational formulation")

# TODO: These maybe have to be "Functions" to allow for nonlinear bilinear form
N = FacetNormal(MESH)
S_N = dot(S, N)
R_N = dot(R, N)
S_T = sqrt(dot(S, S) - S_N)
R_T = sqrt(dot(R, R) - R_N)

A = (inner(dev(grad(S)), grad(R)) + 2/3 * 1/TAU * dot(S, R) - 5/2 * THETA * div(R) -
     div(S)*KAPPA) * dx
    # TODO: Think about how to treat boudary integrals
    #  div(S)*KAPPA) * dx + (5/(4*XI_TILDE) * S_N * R_N + (11*XI_TILDE)/10 * S_T * R_T)*ds

# A += # TODO: Include jump term, see MS thesis from sweden: "jump"

L = (5/2 * THETA_W) * R_N * ds - (F * KAPPA) * dx



# %%
print("Solving system...")
SOL = Function(W)
solve(A == L, SOL, BCS)
(S, THETA) = SOL.split()

print("Plotting solution...")
plt.figure()
plot(S, title="s")
plt.figure()
plot(THETA, title="Theta")
plt.show()

print("Writing output files...")
S_FILE_PVD = dolfin.File("s.pvd")
S_FILE_PVD << S
THETA_FILE_PVD = dolfin.File("theta.pvd")
THETA_FILE_PVD << THETA


#%%
