"""Simple mesh test: Obtains same mesh as in westerkamp2019"""

from mshr import Circle, generate_mesh
from dolfin import Point, FunctionSpace, plot
import matplotlib.pyplot as plt

# Parameters
# 1) Radii
R1 = 0.5
R2 = 2.0
# 2) Mesh resolution
RES = 10

DOMAIN = Circle(Point(0.0, 0.0), R2) - Circle(Point(0.0, 0.0), R1)
MESH = generate_mesh(DOMAIN, RES)
V = FunctionSpace(MESH, "Lagrange", 1)

plt.figure()
plot(MESH, title="Mesh")
plt.show()
