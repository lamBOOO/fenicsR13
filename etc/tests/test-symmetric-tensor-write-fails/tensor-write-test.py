from dolfin import *

mesh = UnitSquareMesh(16, 16)

# works
u = Function(TensorFunctionSpace(mesh, "Lagrange", 1))
u.vector()[:] = 1.0
File("u.pvd") << u
print("wrote u")

# fails
v = Function(TensorFunctionSpace(mesh, "Lagrange", 1, symmetry=True))
v.vector()[:] = 1.0
File("v.pvd") << v
print("wrote v")