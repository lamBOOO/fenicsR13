from dolfin import *

conductivity_code = """

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <dolfin/function/Expression.h>

class Conductivity : public dolfin::Expression
{
public:

  // Create expression with 3 components
  Conductivity() : dolfin::Expression(1) {}

  // Function for evaluating expression on each cell
  void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
  {
    values[0] = R;
  }

    // The data stored in mesh functions
  double R;

};

PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<Conductivity, std::shared_ptr<Conductivity>, dolfin::Expression>
    (m, "Conductivity")
    .def(py::init<>())
    .def_readwrite("R", &Conductivity::R);
}

"""

c = CompiledExpression(compile_cpp_code(conductivity_code).Conductivity(), degree=0, R=1.23)
print(2*c(1))

f = Expression("c_", degree=2, c_=c)
print(f(1))