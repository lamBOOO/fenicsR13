print("start0")

import dolfin as df

print("start")

mesh = df.RectangleMesh(df.Point(1.0, 1.0), df.Point(2.0, 2.0), 10, 10)
# mesh = UnitIntervalMesh(10)

el = df.FiniteElement("Lagrange", mesh.ufl_cell(), degree=1)
V = df.FunctionSpace(mesh, el)
x = df.SpatialCoordinate(mesh)

p_code = """
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

class Pressure : public dolfin::Expression
{
public:

  // Create expression with 3 components
  Pressure() : dolfin::Expression() {}

  // Function for evaluating expression on each cell
  void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
  {
    double C_0 = -50.80230139855979;
    double C_1 = 0.6015037593984962;
    double C_2 = 0;
    double C_3 = -444.7738727200452;
    double C_4 = -0.12443443849461801;
    double C_5 = 9.38867688999618;
    double C_6 = -0.6917293233082705;
    double C_7 = 0;
    double C_8 = 0;
    double C_9 = 0;
    double C_10 = 0;
    double C_11 = 0;
    double C_12 = 2.255312046238658E-11;
    double C_13 = 7.2248457002586;
    double C_14 = -104.89346597195336;
    double C_15 = 4.870715709115059E-7;

    double tau = 0.1;
    double A_1 = 0.4;

    double lambda_1 = sqrt(5.0/9.0);
    double lambda_2 = sqrt(5.0/6.0);
    double lambda_3 = sqrt(3.0/2.0);

    double R = sqrt(pow(x[0],2)+pow(x[1],2));
    double phi = atan2(x[1],x[0]);

    double d_0 = C_9 + C_2*cyl_bessel_k(0, R*lambda_2/tau) + C_8*cyl_bessel_i(0, R*lambda_2/tau);
    double d = - (10*A_1 * pow(R,2))/(27*tau) + (4*C_4*R)/(tau) - (2*C_5*tau)/R + C_14*cyl_bessel_k(1, R*lambda_2/tau) + C_15* cyl_bessel_i(1, R*lambda_2/tau);

    values[0] = d_0 + cos(phi) * d;
  }

    // The data stored in mesh functions
  double R;

};

PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<Pressure, std::shared_ptr<Pressure>, dolfin::Expression>
    (m, "Pressure")
  .def(py::init<>());
}
"""

c = df.CompiledExpression(df.compile_cpp_code(p_code).Pressure(), degree=1)
print(2*c([1,1]))

p_i = df.interpolate(c, V)
df.plot(p_i)

file = df.File("out.pvd")
file.write(p_i)