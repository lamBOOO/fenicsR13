from dolfin import *
import matplotlib.pyplot as plt
import ufl

mesh = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), 10, 10)

el = FiniteElement("Lagrange", mesh.ufl_cell(), degree=1)
V = FunctionSpace(mesh, el)
x = SpatialCoordinate(mesh)

# Approach 1: Fails as UFL expression with project
phi = project(ufl.atan_2(x[1], x[0]), V)
plt.subplot(211)
p1 = plot(phi, title="project(ufl.atan_2(x[1], x[0]), V)")
plt.colorbar(p1)
file = File("phi1.pvd")
file.write(phi)

# Approach 2: Works as UserExpression with project
class RadialAngle(UserExpression):
    def eval(self, values, x):
        values[0] = ufl.atan_2(x[1], x[0])
phi2 = project(RadialAngle(), V)
plt.subplot(212)
p2 = plot(phi2, title="project(RadialAngle(), V)")
plt.colorbar(p2)
file2 = File("phi2.pvd")
file2.write(phi2)

plt.show()


