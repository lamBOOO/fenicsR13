"""Program to solve the heat system"""

# %%
from mshr import Circle, generate_mesh
from dolfin import Point, plot, set_log_level, FiniteElement, VectorElement, \
    FunctionSpace, Mesh, MeshFunction, TrialFunctions, TestFunctions, Expression, \
    dx, dot, grad, dev, div, ds, FacetNormal, inner, DirichletBC, Constant
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

# Boundary conditions
# Inner=3000
BC_INNER_THETA = DirichletBC(W.sub(0), Constant(0.0), BOUNDARIES, 3000)
BC_INNER_S = DirichletBC(W.sub(1), Constant((0.0, 0.0)), BOUNDARIES, 3000)

# Outer=3100
BC_OUTER_THETA = DirichletBC(W.sub(0), Constant(0.0), BOUNDARIES, 3100)
BC_OUTER_S = DirichletBC(W.sub(1), Constant((0.0, 0.0)), BOUNDARIES, 3100)

BCS = [BC_INNER_THETA, BC_INNER_S, BC_OUTER_THETA, BC_OUTER_S]

# %%
# Specifiy BCs

# %%
# Parameters
TAU = 0.1
XI_TILDE = 0.5  # TODO: What's the value for this?, check against normal xi
THETA_W = 0.5  # TODO: Implement BC

# Define trial and testfunction
(THETA, S) = TrialFunctions(W)
(KAPPA, R) = TestFunctions(W)

# Define source function
F = Expression("2 - pow(x[0],2) - pow(x[1],2)", degree=2)  # f=2-x^2-y^2

# Define (combined) variational form
N = FacetNormal(MESH)

# A = (inner(dev(grad(S)), grad(R)) + 2/3 * 1/TAU * dot(S, R) - 5/2 * THETA * div(R) - div(S)*KAPPA) * \
    # dx + (5/(4*XI_TILDE) * dot(S, N) * dot(R, N) +
    #       (11*XI_TILDE)/10 * (S-dot(S, N)) + (R-dot(R, N)))*ds

# L = (5/2 * THETA_W)*dot(R, N)*ds - (F * KAPPA) * dx

# %%
