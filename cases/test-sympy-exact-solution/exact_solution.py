#!/usr/bin/env python3


# ------------------------------------------------------------------------------
# PYLINT SETTINGS
# pylint: disable=invalid-name
# pylint: disable=unused-import
# ------------------------------------------------------------------------------


"""
Exact solution of heat system
=============================================

TODOs
-----
- Add generalization with Christoffel symbols
- Add generalization with contraviriant deriv? see matlab file
- Try to reduce to a first oder system
- Check if bessel ODEs are now implemented in sympy

Informations
------------
- Coupled system not solvaeable, classify_sysode say no type of [Eq,Eq]

Literature
----------
- Mathematical Methods for Physicists, Seventh Edition:
  A Comprehensive Guide (7th) page 224
- http://mathworld.wolfram.com/CylindricalCoordinates.html
- CRC Standard Mathematical Tables page 388
"""


# ------------------------------------------------------------------------------
# IMPORTS
# ------------------------------------------------------------------------------
from sympy import (Symbol, Function, cos, sin, diff, Matrix,
                   init_printing, eye, simplify, Rational, Array, sqrt, dsolve,
                   Derivative, Eq, solve, log, ln, collect, expand,
                   classify_ode, symbols)
from sympy.physics.vector import ReferenceFrame, gradient
from sympy.vector import CoordSys3D
from sympy.solvers.ode import classify_sysode
from sympy.printing import print_ccode, print_latex, cxxcode
from IPython.display import display
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# DEFINE VARIABLES, PARAMETRIZATION, JACOBIAN, METRIC, CHRISTOFFEL
# ------------------------------------------------------------------------------
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
        Rational(1, 2) * sum([
            g_inv[n, k] * (
                diff(g[i, k], xx[j]) + diff(g[j, k], xx[i])
                - diff(g[i, j], xx[k])
            ) for k in range(dim)
        ])
        for j in range(dim)
    ] for i in range(dim)] for n in range(dim)
    ]) # Gamma[(n,i,j)] means looking into n-th Christoffel matrix
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# DEFINE OPERATORS
# ------------------------------------------------------------------------------
def grad0(rank0_):
    "Calculates gradient of rank-0 tensor = scalar"
    return Array([diff(rank0_, x) for x in xx])

def laplace0(rank0_):
    "Copied from Mathematica script"
    temp1 = Array([diff(rank0_, x) for x in xx])
    temp2 = Array([[
        diff(temp1[i], xx[k])
        - sum([Gamma[sigma, i, k] * temp1[sigma] for sigma in range(dim)]) for k in range(dim)
    ] for i in range(dim)])
    return simplify(sum([g_inv[k, k] * temp2[k, k] for k in range(dim)]))

def grad1(rank1_):
    "Calculates gradient of rank-0 tensor = scalar"
    return Array([
        1/sqrt(g[k, k]) * diff(rank1_[i], xx[k]) for k in range(dim)
        ] for i in range(dim))

def tracefreeSym2(rank2_):
    "Calculates deviatoric (tracefree symmateric) parts of rank-2 tensor"
    rank2_ = rank2_.tomatrix()
    return (Rational(1, 2) * (rank2_ + rank2_.transpose())
            - Rational(1, 3) * rank2_.trace() * eye(dim))

def div2():
    "Calculates divergence of rank-2 tensor"
    return "nothing"
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# DEFINE ANSATZ AND EQUATION
# ------------------------------------------------------------------------------
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
s = Array([alpha_0 + cos(phi)*alpha, beta_0 + cos(phi)*beta, 0])
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# SETUP ODES AND SOLVE THEM
# ------------------------------------------------------------------------------
A_0 = Symbol("A0")
A_1 = Symbol("A1")
A_2 = Symbol("A2")
F = A_0 + A_2*R**2 + A_1 * cos(phi)
eqn1 = laplace0(theta) - Rational(2, 3)*laplace0(F) - F
odes = eqn1.coeff(cos(phi), 1), eqn1.coeff(cos(phi), 0)

C_11 = Symbol("C_11")
C_12 = Symbol("C_12")
C_21 = Symbol("C_21")
C_22 = Symbol("C_22")

if odes[1].coeff(c, 1) == 0 and odes[0].coeff(c_0, 1) == 0:
    sol1 = dsolve(odes[0], c).subs([(Symbol("C1"), C_11), (Symbol("C2"), C_12)])
    sol2 = dsolve(odes[1], c_0).subs([(Symbol("C1"), C_21), (symbols("C2"), C_22)])
    theta = expand(theta.subs([(c, sol1.rhs), (c_0, sol2.rhs)]))
    display(theta)
else:
    raise Exception("Theta system is not decoupled but it should be.")
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# CALCUALTE INTEGRATION CONSTANT FROM BOUNDARY CONDITIONS
# Assume theta is fixed, in theorz: theta_wall = const + cos(phi)
# ------------------------------------------------------------------------------
theta_0 = Symbol('theta_0')
theta_1 = Symbol('theta_1')

R_0 = Symbol('R_0')
R_1 = Symbol('R_1')

i_constants = solve([
    Eq(theta.subs(R, R_0).coeff(cos(phi), 0), theta_0),
    Eq(theta.subs(R, R_0).coeff(cos(phi), 1), 0),
    Eq(theta.subs(R, R_1).coeff(cos(phi), 0), theta_1),
    Eq(theta.subs(R, R_1).coeff(cos(phi), 1), 0)
    ], [C_11, C_12, C_21, C_22])
theta = simplify(theta.subs([
    (C_11, i_constants[C_11]),
    (C_12, i_constants[C_12]),
    (C_21, i_constants[C_21]),
    (C_22, i_constants[C_22])
    ]))
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# FINAL SOLUTION FOR THETA WITH ACTUAL VALUES
# ------------------------------------------------------------------------------
theta_0_v = Rational(1)
theta_1_v = Rational(1, 2)
R_0_v = Rational(1, 2)
R_1_v = Rational(2)
A_0_v = Rational(2)
A_1_v = Rational(0)
A_2_v = Rational(-1)
theta_final = simplify(theta.subs([
    (theta_0, theta_0_v),
    (theta_1, theta_1_v),
    (R_0, R_0_v),
    (R_1, R_1_v),
    (A_0, A_0_v),
    (A_1, A_1_v),
    (A_2, A_2_v)
    ]))
display(theta_final)
print(cxxcode(theta_final)) # with namespaces
print_ccode(theta_final)
# ------------------------------------------------------------------------------
