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

        double c_0 = C_2 + (- (20.0*(C_1)*std::log(R)) + ((5.0/4.0)*std::pow(R, 4)) - (2.0*std::pow(R, 2)*(24.0*std::pow(tau, 2) + 5.0)))/(tau*75.0);
        double c = 0.0;

        double theta = c_0 + cos(phi) * c;

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

        double alpha_0 = (C_1)/R + (pow(R,2) - (pow(R,4)/4))/R;
        double alpha = 0;

        double beta_0 = 0;
        double beta = 0;

        double s_R = alpha_0 + cos(phi) * alpha;
        double s_phi = beta_0 - sin(phi) * beta;

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