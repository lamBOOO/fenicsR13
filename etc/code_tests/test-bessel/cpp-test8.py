from dolfin import *

# parameters["form_compiler"]["external_includes"] = "test:test2"
# parameters["form_compiler"].update({"external_includes": ("header1",)} ) # fails
parameters["form_compiler"].update({"external_includes": 1} )
# parameters["form_compiler"]["external_includes"] = ""
print(parameters["form_compiler"]["optimize"])

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