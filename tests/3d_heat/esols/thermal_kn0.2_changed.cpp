#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

double BI( int n, double x ) { return boost::math::cyl_bessel_i(n,x); }
double BK( int n, double x ) { return boost::math::cyl_bessel_k(n,x); }

class Temperature : public dolfin::Expression {
    public:
    Temperature() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double r=sqrt(pow(z[0],2)+pow(z[1],2)+pow(z[2],2));
        double teta = atan2(z[1],z[0]);
	double phi1 = acos(z[2]);
	double x = z[0];
	double y = z[1];
	double v = z[2];
	double E = 2.7182818284;

        values[0] = 2.0824823984706486 - 0.18650627127337152/sqrt(std::pow(v,2) + std::pow(x,2) + std::pow(y,2)) ;
    }
};

class Heatflux : public dolfin::Expression {
    public:
    Heatflux() : dolfin::Expression(3) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double r=sqrt(pow(z[0],2)+pow(z[1],2)+pow(z[2],2));
        double teta = atan2(z[1],z[0]);
	double phi1 = acos(z[2]);
	double x = z[0];
	double y = z[1];
	double v = z[2];
	double E = 2.7182818284;

        values[0] = 0. + (x*(0. - 0.13987970345502865/(std::pow(v,2) + std::pow(x,2) + std::pow(y,2)))*sqrt(1 - std::pow(v,2)/(std::pow(v,2) + std::pow(x,2) + std::pow(y,2))))/sqrt(std::pow(x,2) + std::pow(y,2)) ;
        values[1] = 0. + (y*(0. - 0.13987970345502865/(std::pow(v,2) + std::pow(x,2) + std::pow(y,2)))*sqrt(1 - std::pow(v,2)/(std::pow(v,2) + std::pow(x,2) + std::pow(y,2))))/sqrt(std::pow(x,2) + std::pow(y,2)) ;
        values[2] = 0. + (v*(0. - 0.13987970345502865/(std::pow(v,2) + std::pow(x,2) + std::pow(y,2))))/sqrt(std::pow(v,2) + std::pow(x,2) + std::pow(y,2)) ;
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
