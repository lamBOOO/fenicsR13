import matplotlib.pyplot as plt
import fenics as fe
from ufl import bessel_I

mesh = fe.IntervalMesh(100, 0, 10)
V = fe.FunctionSpace(mesh,'P',1)

x = fe.SpatialCoordinate(mesh)
bs = bessel_I(0,x[0])
expr = fe.Expression('cos(x[0]*b)', b=10, element=V.ufl_element())

expr_bessel = expr*bs
proj_expr = fe.project(expr_bessel, V)

# fe.plot(proj_expr)
# plt.show()



# class MyExpression(fe.UserExpression):
#     def eval(self, values, x):
#         values[0] = x[0] ** 2
#     def value_shape(self):
#         return ()
# expr = MyExpression(degree=1)

# class InitialConditions(fe.UserExpression):
#     def __init__(self, **kwargs):
#         # random.seed(2 + MPI.rank(MPI.comm_world))
#         super().__init__(**kwargs)
#     def eval(self, values, x):
#         values[0] = 0.63
#     def value_shape(self):
#         return (1,)
# expr = InitialConditions(degree=1)

# class K(UserExpression):
#     def set_k_values(self, k_0, k_1):
#         self.k_0, self.k_1 = k_0, k_1
#     def eval(self, value, x):
#         "Set value[0] to value at point x"
#         tol = 1E-14
#         if x[1] <= 0.5 + tol:
#             value[0] = self.k_0
#         else:
#             value[0] = self.k_1
# # Initialize kappa
# kappa = K(degree=0)

# class F0(UserExpression):
#     def eval(self, values, x):
#         values[0] = sin(3.0*x[0])*sin(3.0*x[1])*sin(3.0*x[2])

# f0 = F0(name="f0", label="My expression", degree=2)

class F0(UserExpression):
    def eval(self, values, x):
        values[0] = sin(3.0*x[0])*sin(3.0*x[1])*sin(3.0*x[2])

class F1(UserExpression):
    def __init__(self, mesh, *arg, **kwargs):
        super().__init__(*arg, **kwargs)
        self.mesh = mesh

    def eval_cell(self, values, x, cell):
        c = Cell(self.mesh, cell.index)
        values[0] = sin(3.0*x[0])*sin(3.0*x[1])*sin(3.0*x[2])

e0 = F0(degree=2)
e1 = F1(mesh, degree=2)
e2 = Expression("sin(3.0*x[0])*sin(3.0*x[1])*sin(3.0*x[2])", degree=2)

u0 = interpolate(e0, V)
u1 = interpolate(e1, V)
u2 = interpolate(e2, V)

fp = fe.interpolate(f0, V)
fe.plot(fp)
plt.show()