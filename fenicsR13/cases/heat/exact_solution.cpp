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
double C_5 = 39.38867688999618;
double C_6 = -0.6917293233082705;
double C_7 = 0;
double C_8 = 0;
double C_9 = 0;
double C_10 = 0;
double C_11 = 0;
double C_12 = 2.255312046238658E-11;
double C_13 = 407.2248457002586;
double C_14 = -104.89346597195336;
double C_15 = 4.870715709115059E-7;

double tau = 0.1;
double A_1 = 0.4;
// double A_1 = 0.0;

double lambda_1 = sqrt(5.0/9.0);
double lambda_2 = sqrt(5.0/6.0);
double lambda_3 = sqrt(3.0/2.0);

double BI( int n, double x ) { return boost::math::cyl_bessel_i(n,x); }
double BK( int n, double x ) { return boost::math::cyl_bessel_k(n,x); }

double I_0( double x) { return BI(0,x); }
double I_1( double x) { return BI(1,x); }
double I_2( double x) { return BI(2,x); }
double I_3( double x) { return BI(3,x); }

double K_0( double x) { return BK(0,x); }
double K_1( double x) { return BK(1,x); }
double K_2( double x) { return BK(2,x); }
double K_3( double x) { return BK(3,x); }

class Pressure : public dolfin::Expression {
    public:
    Pressure() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double d_0 = C_9 + C_2*K_0( R*lambda_2/tau) + C_8*I_0(R*lambda_2/tau);
        double d = - (10*A_1 * pow(R,2))/(27*tau) + (4*C_4*R)/(tau) - (2*C_5*tau)/R + C_14*K_1( R*lambda_2/tau) + C_15* I_1( R*lambda_2/tau);

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
        double u_phi = b_0 - sin(phi) * b;

        // https://ocw.mit.edu/courses/aeronautics-and-astronautics/
        // ...16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec05.pdf
        double u_x = u_R * cos(phi) - u_phi * sin(phi);
        double u_y = u_R * sin(phi) + u_phi * cos(phi);

        values[0] = u_x ;
        values[1] = u_y ;
    }
};

class Stress : public dolfin::Expression {
    public:
    Stress() : dolfin::Expression(2,2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double gamma_0 = (2*C_7*pow(tau,2))/(pow(R,2)) + (C_10*K_1((R*lambda_3)/tau))/(lambda_3*R) + (C_11*I_1((R*lambda_3)/tau))/(lambda_3*R) + C_2*( (tau*K_1((R*lambda_2/tau)))/(2*lambda_2*R) - K_2(R*lambda_2/tau) ) + C_8*( -(tau*I_1(R*lambda_2/tau))/(2*lambda_2*R) - I_2(R*lambda_2/tau) );
        double gamma = (7*A_1*pow(R,2))/(27*tau) - (2*pow(tau,3)*(C_3-(64*C_5)/15))/pow(R,3) - (2*C_4*R)/tau - (2*C_5*tau)/R + (C_12*tau*I_2(R*lambda_3/tau))/(lambda_3*R) + (C_13*tau*K_2(R*lambda_3/tau))/(lambda_3*R) + C_14*( -K_1(R*lambda_2/tau)-(3*tau*K_2(R*lambda_2/tau))/(2*lambda_2*R) ) + C_15*( (3*tau*I_2(R*lambda_2/tau))/(2*lambda_2*R)-I_1(R*lambda_2/tau) );

        double kappa_0 = C_1/pow(R,2);
        double kappa = -(A_1*pow(R,2))/(9*tau) - (2*pow(tau,3)*(C_3-(64*C_5)/15))/pow(R,3) + (2*C_4*R)/tau + (C_12*tau*I_2(R*lambda_3/tau))/(lambda_3*R) + (C_13*tau*K_2(R*lambda_3/tau))/(lambda_3*R) - (3*C_14*tau*K_2(R*lambda_2/tau))/(2*lambda_2*R) + (3*C_15*tau*I_2(R*lambda_2/tau))/(2*lambda_2*R);

        double omega_0 = (2*C_7*pow(tau,2))/pow(R,2) + C_10*(K_0(R*lambda_3/tau)+(tau*K_1(R*lambda_3/tau)/(lambda_3*R))) + C_11*( (tau*I_1(R*lambda_3/tau))/(lambda_3*R) -I_0(R*lambda_3/tau) ) + C_2*( -1./2.*K_0(R*lambda_2/tau) - (3*K_1(R*lambda_2/tau))/(2*lambda_2*R) ) + C_8*( (3*I_1(R*lambda_2/tau))/(2*lambda_2*R) - 1./2.*I_0(R*lambda_2/tau));
        double omega = (2*A_1*pow(R,2))/(27*tau) - (2*pow(tau,3)*(C_3-(64*C_5)/15.))/(pow(R,3)) - (2*C_4*R)/tau - (2*C_5*tau)/R + C_12*((tau*I_2(R*lambda_3/tau))/(lambda_3*R) - I_1(R*lambda_3/tau)) + C_13*( - K_1(R*lambda_3/tau) + (tau*K_2(R*lambda_3/tau))/(lambda_3*R) ) + C_14*((tau*K_2(R*lambda_2/tau))/(2*lambda_2*R) - 1/2.*K_3(R*lambda_2/tau) ) + C_15*(-(tau*I_2(R*lambda_2/tau))/(2*lambda_2*R) - 1/2.*I_3(R*lambda_2/tau) );

        double sigma_RR = gamma_0 + cos(phi) * gamma;
        double sigma_Rphi = kappa_0 + sin(phi) * kappa;
        double sigma_phiphi = -(omega_0 + cos(phi) * omega);

        // Heinz Schade p377: Assume symmetry of tensor <=> a_Rphi = a_phiR
        double sigma_xx = sigma_RR * pow(cos(phi),2) - (sigma_Rphi) * sin(2*phi) + sigma_phiphi * pow(sin(phi),2);
        double sigma_xy = sigma_Rphi * cos(2*phi) + 0.5*(sigma_RR-sigma_phiphi) * sin(2*phi); // same as in MomentsDG System_xy
        double sigma_yy = sigma_phiphi * pow(cos(phi),2) - (sigma_Rphi) * sin(2*phi) + sigma_RR * pow(sin(phi),2);

        values[0] = sigma_xx ;
        values[1] = sigma_xy ;
        values[2] = sigma_yy ;
    }
};

PYBIND11_MODULE(SIGNATURE, m) {

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