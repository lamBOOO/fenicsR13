#!/usr/bin/env python3
# pylint: disable=invalid-name
# pylint: disable=unused-import

"""
Calculates analytical solution of heat system
TODO: Add generalization with Christoffel symbols
TODO: Add generalization with contraviriant deriv? see matlab file
"""

from sympy import (Symbol, Function, cos, sin, diff, Matrix,
                   init_printing, eye, simplify, Rational, Array, sqrt)
from sympy.physics.vector import ReferenceFrame, gradient
from sympy.vector import CoordSys3D
from sympy.printing import print_ccode, print_latex
from IPython.display import display
# init_printing()

# # Example
# theta, R, tau, C1, C2 = var("theta, R tau C1 C2")
# expr = (cos(R**2)-1)/(34*C1**C2*theta)

# # C++ code:
# print(cxxcode(expr))
# print_ccode(expr)

# # Nice LaTeX output in Jupyter
# display(expr)

# # LaTeX output
# print_latex(expr)

# # Terminal formatted equation output
# print("1")
# pprint(expr, use_unicode=False)
# print(expr)



# Literature:
# - Mathematical Methods for Physicists, Seventh Edition:
#   A Comprehensive Guide (7th) page 224
# - TODO: Compare Itzkov CM notes
# - http://mathworld.wolfram.com/CylindricalCoordinates.html
# -  CRC Standard Mathematical Tables page 388


R = Symbol("R", positive=True)
phi = Symbol("phi")
z = Symbol("z")
xx = [R, phi, z] # coordinates
dim = len(xx)

P = Matrix([R*cos(phi), R*sin(phi), z]) # parametrization
J = P.jacobian(xx)
J_inv = simplify(J.inv())
g = simplify(J.transpose() * J)
g_inv = simplify(g.inv())
Gamma = Array([
    [[
        Rational(1, 2) * sum([g_inv[n, k] * (diff(g[i, k], xx[j])+diff(g[j, k], xx[i])-diff(g[i, j], xx[k])) for k in range(dim)])
        for j in range(dim)
    ] for i in range(dim)] for n in range(dim)
    ]) # Gamma[(n,i,j)] means looking into n-th Christoffel matrix
display(P)
display(J)
display(J_inv)
display(g)
display(g_inv)
display(Gamma)

tau = Symbol("tau")
s_R = Symbol("s_R")
s_phi = Symbol("s_phi")

alpha_0 = Function('alpha_0')(R)
alpha = Function('alpha')(R)
beta_0 = Function('beta_0')(R)
beta = Function('beta')(R)
c_0 = Function('c_0')(R)
c = Function('c')(R)

theta = c_0 + cos(phi)*c
print("theta:")
display(theta)
grad0_theta = Array([diff(theta, R), diff(theta, phi)])
display(grad0_theta)

def grad0(rank0_):
    "Calculates gradient of rank-0 tensor = scalar"
    return Array([diff(rank0_, x) for x in xx])

def grad1(rank1_):
    "Calculates gradient of rank-0 tensor = scalar"
    return Array([1/sqrt(g[k, k]) * diff(rank1_[i], xx[k]) for k in range(dim)] for i in range(dim))

def tracefreeSym2(rank2_):
    "Calculates deviatoric (tracefree symmateric) parts of rank-2 tensor"
    rank2_ = rank2_.tomatrix()
    return (Rational(1, 2) * (rank2_ + rank2_.transpose())
            - Rational(1, 3) * rank2_.trace() * eye(dim))

def div2(rank2_):
    "Calculates divergence of rank-2 tensor"
    return "nothing"

display(grad0(theta))

print("s:")
s = Array([alpha_0 + cos(phi)*alpha, beta_0 + cos(phi)*beta, 0])
display(diff(s, R))
display(grad1(s))
display(tracefreeSym2(grad1(s)))



# def theta(R_, phi_):
#     "Temperature theta"
#     return c_0(R_) + cos(phi_) * c(R_)

# def s(R_, phi_):
#     "Heatflux vector s"
#     return Matrix([alpha_0, 2])

# display(diff(theta(R,phi)))







# display(diff(c, R, phi))

# R = ReferenceFrame('R')
# v = 3*R.x + 4*R.y

# display(v)

# R = CoordSys3D('Rs')
# v = R.i + 4*R.j + 5*R.k


# test = Matrix([1*phi**2, 2, 3])
# diff(test, R)
# diff(test, phi)
