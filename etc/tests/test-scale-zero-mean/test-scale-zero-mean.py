from dolfin import *
import numpy as np

mesh = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), 10, 10)

el = FiniteElement("Lagrange", mesh.ufl_cell(), degree=1)
V = FunctionSpace(mesh, el)
x = SpatialCoordinate(mesh)

f = Expression("x[0]+1", element=V.ufl_element())
f_i = interpolate(f, V)

f_mean = np.mean(f_i.vector()[:])

f_i_s = f_i.vector()[:] - f_mean


# # Approach 3: Works with project
# phi3 = project(Expression("std::atan2(x[1],x[0])", element=V.ufl_element()), V)
# plt.subplot(223)
# p3 = plot(phi3, title="phi3")
# plt.colorbar(p3)
# file3 = File("phi3.pvd")
# file3.write(phi3)

# # Approach 4: Works with interpolate
# phi4 = interpolate(Expression("std::atan2(x[1],x[0])", element=V.ufl_element()), V)
# plt.subplot(224)
# p4 = plot(phi4, title="phi4")
# plt.colorbar(p4)
# file4 = File("phi4.pvd")
# file4.write(phi4)

# plt.rcParams['figure.figsize'] = [15, 10]
# plt.show()

