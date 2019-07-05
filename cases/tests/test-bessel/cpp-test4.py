from dolfin import *
f = Expression("std::cyl_bessel_i(1,x[0])", degree=2)
print(f(1.23))