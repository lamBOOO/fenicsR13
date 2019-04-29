"""Program to solve the heat system"""

# %%
from mshr import Circle, generate_mesh
from dolfin import Point, plot, set_log_level, FiniteElement, VectorElement, \
  FunctionSpace
import matplotlib.pyplot as plt

# Settings
set_log_level(1)  # all logs

# PARAMETERS
# 1) Radii
R1 = 0.5
R2 = 2.0
# 2) Mesh resolution
RES = 10

# Create mesh
DOMAIN = Circle(Point(0.0, 0.0), R2) - Circle(Point(0.0, 0.0), R1)
MESH = generate_mesh(DOMAIN, RES)

# Show mesh
plt.figure()
plot(MESH, title="Mesh")
plt.show()

# Setup function spaces and shape functions
P2 = VectorElement("Lagrange", MESH.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", MESH.ufl_cell(), 1)
ME = P2 * P1 # mixed element space
W = FunctionSpace(MESH, ME)

#%%
