import dolfin as df

# Create mesh and define function space
mesh = df.UnitSquareMesh(32, 32)
V = df.FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    eps = df.DOLFIN_EPS
    return x[0] < eps or x[0] > 1.0 - eps

# Define boundary condition
u0 = df.Constant(0.0)
bc = df.DirichletBC(V, u0, boundary)

# Define variational problem
u = df.TrialFunction(V)
v = df.TestFunction(V)
f = df.Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)",degree=1)
g = df.Expression("sin(5*x[0])",degree=1)
a = df.inner(df.grad(u), df.grad(v))*df.dx
L = f*v*df.dx + g*v*df.ds

# Compute solution
u = df.Function(V)
df.solve(a == L, u, bc)

# Save solution in VTK format
file = df.File("poisson.pvd")
file << u
