#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

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

class Pressure : public dolfin::Expression {
    public:
    Pressure() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double d_0 = C_9 + C_2*cyl_bessel_k(0, R*lambda_2/tau) + C_8*cyl_bessel_i(0, R*lambda_2/tau);
        double d = - (10*A_1 * pow(R,2))/(27*tau) + (4*C_4*R)/(tau) - (2*C_5*tau)/R + C_14*cyl_bessel_k(1, R*lambda_2/tau) + C_15* cyl_bessel_i(1, R*lambda_2/tau);

        values[0] = d_0 + cos(phi) * d;
    }
};

class Velocity : public dolfin::Expression {
    public:
    Velocity() : dolfin::Expression(2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double a_0 = C_7*tau/R;
        double a = -A_1*R*((2*pow(R,2))/(27*pow(tau,2)) - 2.0/3.0) + C_0 - (C_3*pow(tau,2))/(2*pow(R,2)) + (C_4 * pow(R,2))/(2*pow(tau,2)) + C_5 * (std::log(R/tau)-1.0/2.0);

        double b_0 = C_1/(2*R*tau) + C_6*R;
        double b = -A_1*R*((pow(R,2))/(54*pow(tau,2)) - 1.0/3.0) + C_0 + (C_3*pow(tau,2))/(2*pow(R,2)) + (3*C_4 * pow(R,2))/(2*pow(tau,2)) + C_5 * (std::log(R/tau)+1.0/2.0);


        double u_R = a_0 + cos(phi) * a;
        double u_phi = b_0 + cos(phi) * b;

        // https://ocw.mit.edu/courses/aeronautics-and-astronautics/
        // ...16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec05.pdf
        double u_x = u_R * cos(phi) - u_phi * sin(phi);
        double u_y = u_R * sin(phi) + u_phi * cos(phi);

        values[0] = u_x ;
        values[1] = u_y ;
    }
};

PYBIND11_MODULE(SIGNATURE, m) {
    py::class_<Pressure, std::shared_ptr<Pressure>, dolfin::Expression>
        (m, "Pressure")
    .def(py::init<>());
    py::class_<Velocity, std::shared_ptr<Velocity>, dolfin::Expression>
        (m, "Velocity")
    .def(py::init<>());
}