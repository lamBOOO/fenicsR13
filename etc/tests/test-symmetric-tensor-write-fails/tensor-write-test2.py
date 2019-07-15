"This solves the problem with a hack"

from dolfin import *

mesh = UnitSquareMesh(4, 4)

# works
u = Function(TensorFunctionSpace(mesh, "Lagrange", 1))
u.vector()[:] = 1.0
File("u.pvd") << u
print("wrote u")

# fails
v = Function(TensorFunctionSpace(mesh, "Lagrange", 1, symmetry=True))
v.vector()[:] = 1.0
# v.assign(u) # also does not work because v afterwards lives in symmetric space
# File("v.pvd") << v # fails
# print("wrote v")

# projection works up to machine precision
vp = project(v, TensorFunctionSpace(mesh, "Lagrange", 1))
error = vp.vector()[:] - u.vector()[:]
print(error)
File("vp.pvd") << vp
print("wrote vp")
