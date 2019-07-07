from dolfin import *
import matplotlib.pyplot as plt

mesh = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), 10, 10)

el = FiniteElement("Lagrange", mesh.ufl_cell(), degree=1)
V = FunctionSpace(mesh, el)
x = SpatialCoordinate(mesh)

# Approach 1: Fails as UFL expression with project
phi = interpolate(Expression("2*boost::math::cyl_bessel_ee(0,x[0])", degree=2), V)
p1 = plot(phi, title="project(ufl.atan_2(x[1], x[0]), V)")
plt.colorbar(p1)
file = File("phi1.pvd")
file.write(phi)