#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

double tau = 0.1;
double C_1 = -0.40855716127979214;
double C_2 = 2.4471587630476663;

class Temperature : public dolfin::Expression {
    public:
    Temperature() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double theta = (-76*std::pow(R,2))/15. + (5*std::pow(R,4))/8. + (11*(-36817 + 148480*std::log(5) - 24880*std::log(20)))/(1920.*(-23 + 80*std::log(5) - 80*std::log(20))) - (5665*std::log(10*R))/(8.*(-23 + 80*std::log(5) - 80*std::log(20)));

        values[0] = theta;
    }
};

class Heatflux : public dolfin::Expression {
    public:
    Heatflux() : dolfin::Expression(2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double s_R = R - std::pow(R,3)/4. + 1133/(16.*R*(-23 + 80*std::log(5) - 80*std::log(20)));
        double s_phi = 0;

        double s_x = s_R * cos(phi) - s_phi * sin(phi);
        double s_y = s_R * sin(phi) + s_phi * cos(phi);

        values[0] = s_x ;
        values[1] = s_y ;
    }
};


PYBIND11_MODULE(SIGNATURE, m) {

    py::class_<Temperature, std::shared_ptr<Temperature>, dolfin::Expression>
        (m, "Temperature")
    .def(py::init<>());

    py::class_<Heatflux, std::shared_ptr<Heatflux>, dolfin::Expression>
        (m, "Heatflux")
    .def(py::init<>());

}