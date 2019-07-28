print("start0")

import dolfin as df

print("start")

mesh = df.RectangleMesh(df.Point(-1.0, -1.0), df.Point(1.0, 1.0), 10, 10)
# mesh = UnitIntervalMesh(10)

el = df.FiniteElement("Lagrange", mesh.ufl_cell(), degree=1)
V = df.FunctionSpace(mesh, el)
x = df.SpatialCoordinate(mesh)

conductivity_code = """
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <dolfin/function/Expression.h>

class RadialAngle : public dolfin::Expression
{
public:

  // Create expression with 3 components
  RadialAngle() : dolfin::Expression() {}

  // Function for evaluating expression on each cell
  void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
  {
    values[0] = atan2(x[1],x[0]);
  }

    // The data stored in mesh functions
  double R;

};

PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<RadialAngle, std::shared_ptr<RadialAngle>, dolfin::Expression>
    (m, "RadialAngle")
    .def(py::init<>())
    .def_readwrite("R", &RadialAngle::R);
}
"""

c = df.CompiledExpression(df.compile_cpp_code(conductivity_code).RadialAngle(), degree=1, R=1.23)
print(2*c(1))
