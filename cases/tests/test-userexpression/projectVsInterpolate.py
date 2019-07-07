from dolfin import *
import matplotlib.pyplot as plt
import ufl

mesh = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), 10, 10)

el = FiniteElement("Lagrange", mesh.ufl_cell(), degree=1)
V = FunctionSpace(mesh, el)
x = SpatialCoordinate(mesh)

# Approach 1: Fails with project
phi = project(ufl.atan_2(x[1], x[0]), V)
# phi = interpolate(ufl.atan_2(x[1], x[0]), V) # fails, no eval() method..
plt.subplot(221)
p1 = plot(phi, title="phi1")
plt.colorbar(p1)
file = File("phi1.pvd")
file.write(phi)

# Approach 2: Works with interpolate and UserExpression (has eval method)
class RadialAngle(UserExpression):
    def eval(self, values, x):
        values[0] = ufl.atan_2(x[1], x[0])
phi2 = interpolate(RadialAngle(), V)
plt.subplot(222)
p2 = plot(phi2, title="phi2")
plt.colorbar(p2)
file2 = File("phi2.pvd")
file2.write(phi2)

# Approach 3: Works with project
phi3 = project(Expression("std::atan2(x[1],x[0])", element=V.ufl_element()), V)
plt.subplot(223)
p3 = plot(phi3, title="phi3")
plt.colorbar(p3)
file3 = File("phi3.pvd")
file3.write(phi3)

# Approach 4: Works with interpolate
phi4 = interpolate(Expression("std::atan2(x[1],x[0])", element=V.ufl_element()), V)
plt.subplot(224)
p4 = plot(phi4, title="phi4")
plt.colorbar(p4)
file4 = File("phi4.pvd")
file4.write(phi4)

plt.rcParams['figure.figsize'] = [15, 10]
plt.show()

