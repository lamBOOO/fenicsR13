"""Simple mesh test: Obtains same mesh as in westerkamp2019"""

from mshr import Circle, generate_mesh
import dolfin
import matplotlib.pyplot as plt

# Parameters
# 1) Radii
R1 = 0.5
R2 = 2.0
# 2) Mesh resolution
RES = 10

DOMAIN = Circle(dolfin.Point(0.0, 0.0), R2) - Circle(dolfin.Point(0.0, 0.0), R1)
MESH = generate_mesh(DOMAIN, RES)

plt.figure()
dolfin.plot(MESH, title="Subtraction test")
plt.show()
