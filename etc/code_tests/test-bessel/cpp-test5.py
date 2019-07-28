from dolfin import *
import ufl
print(ufl.bessel_I(1, x[0])(1))
R_ = Expression("sqrt(pow(x[0],2)+pow(x[1],2))", degree=2) # radius
f = Expression("3.14 * I1", degree=2, I1=ufl.bessel_I(1, R_), R=R_)
print(f(1.23)) # fails