from dolfin import *
import numpy

# System size definition
# system = (x1, x1+Lx) x (y1, y1+Ly)
Lx = 2.0; Ly = 2.0
x1 = 10.0; y1 = 10.0

c0 = 1.0
D0 = 1.0
x0 = 1.5
y0 = 1.5

# Create mesh and define function space
dx0 = dx1 = 0.1
nx = int(Lx/dx0)
ny = int(Ly/dx1)
mesh = RectangleMesh(Point(x1, y1), Point(x1+Lx, y1+Ly), nx, ny,"right/left")
V = FunctionSpace(mesh, 'Lagrange', 1)

code = '''
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
using boost::math::cyl_bessel_i;
using boost::math::cyl_bessel_j;
using boost::math::cyl_bessel_k;
using boost::math::cyl_neumann;

namespace dolfin {
      class MyFun : public Expression
    {
            double c0,D0,x0,y0;
            public:
                      MyFun(): Expression() {};
            void eval(Array<double>& values, const Array<double>& x) const {
                    double f = D0*sqrt( logf(x[0]/x0)*logf(x[0]/x0) + logf(x[1]/y0)*logf(x[1]/y0) )/c0 ;
                    values[0] = cyl_bessel_k(0,f);
            }

            void update(double _c0, double _D0, double _x0, double _y0) {
                  c0 = _c0;
                  D0 = _D0;
                  x0 = _x0;
                  y0 = _y0;
          }
    };
}'''

f=Expression(code)
f.update(c0,D0,x0,y0)

ss = interpolate(f,V)
plot(ss, interactive=True)