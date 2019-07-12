from dolfin import *
mesh = UnitIntervalMesh(100)
V = FunctionSpace(mesh,"Lagrange",1)

def lumpedProject(f,V):
    v = TestFunction(V)
    lhs = assemble(inner(Constant(1.0),v)*dx)
    rhs = assemble(inner(f,v)*dx)
    u = Function(V)
    as_backend_type(u.vector())\
        .vec().pointwiseDivide(as_backend_type(rhs).vec(),
                               as_backend_type(lhs).vec())
    return u

# Compare results:
from matplotlib import pyplot as plt
x = SpatialCoordinate(mesh)
f = 0.2*sin(2.0*pi*x[0]) + conditional(gt(x[0],0.5),1.0,Constant(0.0))
plot(project(f,V))
plot(lumpedProject(f,V))
plt.show()