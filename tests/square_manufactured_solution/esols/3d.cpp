#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

class Temperature : public dolfin::Expression {
    public:
    Temperature() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

    values[0] = -(std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2));  // temperature
    }
};

class Heatflux : public dolfin::Expression {
    public:
    Heatflux() : dolfin::Expression(3) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        values[0] = std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // heat flux x
        values[1] = std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // heat flux y
        values[2] = std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // heat flux z
    }
};

class Pressure : public dolfin::Expression {
    public:
    Pressure() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        values[0] = 0.000244140625 - std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // pressure
    }
};

class Velocity : public dolfin::Expression {
    public:
    Velocity() : dolfin::Expression(3) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        values[0] = -2*std::pow(-1 + x[0],2)*std::pow(x[0],2)*(-1 + x[1])*x[1]*(-1 + x[2])*x[2]*(x[1] + (-1 + x[2])*x[2] - 2*x[1]*std::pow(x[2],2) + std::pow(x[1],2)*(-1 + 2*x[2]));  // velocity x
        values[1] = 2*(-1 + x[0])*x[0]*std::pow(-1 + x[1],2)*std::pow(x[1],2)*(-1 + x[2])*x[2]*(x[0] + (-1 + x[2])*x[2] - 2*x[0]*std::pow(x[2],2) + std::pow(x[0],2)*(-1 + 2*x[2]));  // velocity y
        values[2] = -2*(-1 + x[0])*x[0]*(-1 + x[1])*x[1]*(x[0] + (-1 + x[1])*x[1] - 2*x[0]*std::pow(x[1],2) + std::pow(x[0],2)*(-1 + 2*x[1]))*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // velocity z
    }
};

class Stress : public dolfin::Expression {
    public:
    Stress() : dolfin::Expression(3,3) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        values[0] = std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor xx
        values[1] = std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor xy
        values[2] = -2*std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor xz
        values[3] = std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor yx
        values[4] = std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor yy
        values[5] = -2*std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor yz
        values[6] = -2*std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor zx
        values[7] = -2*std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor zy
        values[8] = -2*std::pow(-1 + x[0],2)*std::pow(x[0],2)*std::pow(-1 + x[1],2)*std::pow(x[1],2)*std::pow(-1 + x[2],2)*std::pow(x[2],2);  // stress tensor zz
    }
};

PYBIND11_MODULE(SIGNATURE, m) {

    py::class_<Temperature, std::shared_ptr<Temperature>, dolfin::Expression>
        (m, "Temperature")
    .def(py::init<>());

    py::class_<Heatflux, std::shared_ptr<Heatflux>, dolfin::Expression>
        (m, "Heatflux")
    .def(py::init<>());

    py::class_<Pressure, std::shared_ptr<Pressure>, dolfin::Expression>
        (m, "Pressure")
    .def(py::init<>());

    py::class_<Velocity, std::shared_ptr<Velocity>, dolfin::Expression>
        (m, "Velocity")
    .def(py::init<>());

    py::class_<Stress, std::shared_ptr<Stress>, dolfin::Expression>
        (m, "Stress")
    .def(py::init<>());

}
