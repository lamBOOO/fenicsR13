import matplotlib.pyplot as plt
from dolfin import *

# Next, some parameters for the form compiler are set::

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# A mesh is created, and a quadratic finite element function space::

# Make mesh ghosted for evaluation of DG terms
parameters["ghost_mode"] = "shared_facet"

# Create mesh and define function space
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "CG", 2)

# A subclass of :py:class:`SubDomain <dolfin.cpp.SubDomain>`,
# ``DirichletBoundary`` is created for later defining the boundary of
# the domain::

# Define Dirichlet boundary
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# A subclass of :py:class:`Expression
# <dolfin.functions.expression.Expression>`, ``Source`` is created for
# the source term :math:`f`::

class Source(UserExpression):
    def eval(self, values, x):
        values[0] = 4.0*pi**4*sin(pi*x[0])*sin(pi*x[1])



# The Dirichlet boundary condition is created::

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, DirichletBoundary())

# On the finite element space ``V``, trial and test functions are
# created::

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)

# A function for the cell size :math:`h` is created, as is a function
# for the average size of cells that share a facet (``h_avg``).  The UFL
# syntax ``('+')`` and ``('-')`` restricts a function to the ``('+')``
# and ``('-')`` sides of a facet, respectively. The unit outward normal
# to cell boundaries (``n``) is created, as is the source term ``f`` and
# the penalty parameter ``alpha``. The penalty parameters is made a
# :py:class:`Constant <dolfin.functions.constant.Constant>` so that it
# can be changed without needing to regenerate code. ::

# Define normal component, mesh size and right-hand side
h = CellDiameter(mesh)
h_avg = (h('+') + h('-'))/2.0
n = FacetNormal(mesh)
f = Source(degree=2)

test = interpolate(f, V)
plot(test)
plt.show()

# works on imac