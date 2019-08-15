#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

double tau = 1.0;
double A_1 = 1.0;

// Constants
double C_0 = -1.004419665507595;
double C_1 = 0.4448013383516245;
double C_2 = 0.2112530227632378;
double C_3 = -0.4315096542509863;
double C_4 = 0.08819320821751526;
double C_5 = -0.06836832840389916;
double C_6 = 0.09335102518735922;
double C_7 = 0;
double C_8 = 0.5183777145707696;
double C_9 = -0.0026285115884281396; // FIXME: Check scaling
double C_10 = -0.06963847872220634;
double C_11 = 0.01839832993369312;
double C_12 = -0.03561221007805056;
double C_13 = 0.2707355692561467;
double C_14 = 0.04291262895440975;
double C_15 = -0.02500281007035313;
double C_16 = 0.03610593775576753;
double C_17 = -0.006316859711714699;
double C_18 = -0.01923465390597724;
double C_19 = -0.001999331311951152;
double C_20 = -0.09987535948000502;
double C_21 = -0.002620826571728536;
double C_22 = -0.01585859733640677;
double C_23 = 0.8086861581439194;

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

class Temperature : public dolfin::Expression {
    public:
    Temperature() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double c_0 = C_8 + (2*C_17*I_0((R*lambda_2)/tau))/5. + (2*C_16*K_0((R*lambda_2)/tau))/5. - (4*C_6*std::log(R/tau))/15.;
        double c = (-4*C_2*R)/(15.*tau) + (8*C_4*R)/(5.*tau) - (4*A_1*std::pow(R,2))/(27.*tau) - (2*C_1*tau)/(15.*R) - (4*C_5*tau)/(5.*R) + (2*C_23*I_1((R*lambda_2)/tau))/5. + (2*C_22*K_1((R*lambda_2)/tau))/5.;

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

        double alpha_0 = (C_6*tau)/R;
        double alpha = C_2 - (C_1*std::pow(tau,2))/(2.*std::pow(R,2)) - (5*((2*C_14*tau*I_1((R*lambda_1)/tau))/(R*lambda_1) + (2*C_15*tau*K_1((R*lambda_1)/tau))/(R*lambda_1)))/2.;

        double beta_0 = C_11*I_1((R*lambda_1)/tau) + C_12*K_1((R*lambda_1)/tau);
        double beta = C_2 + (C_1*std::pow(tau,2))/(2.*std::pow(R,2)) - (5*(C_14*(I_0((R*lambda_1)/tau) + I_2((R*lambda_1)/tau)) - C_15*(K_0((R*lambda_1)/tau) + K_2((R*lambda_1)/tau))))/2.;

        double s_R = alpha_0 + cos(phi) * alpha;
        double s_phi = beta_0 - sin(phi) * beta;

        // https://ocw.mit.edu/courses/aeronautics-and-astronautics/
        // ...16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec05.pdf
        double s_x = s_R * cos(phi) - s_phi * sin(phi);
        double s_y = s_R * sin(phi) + s_phi * cos(phi);

        values[0] = s_x ;
        values[1] = s_y ;
    }
};

class Pressure : public dolfin::Expression {
    public:
    Pressure() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R   = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double d_0 = C_9 + C_17*I_0((R*lambda_2)/tau) + C_16*K_0((R*lambda_2)/tau);
        double d = (4*C_4*R)/tau - (10*A_1*std::pow(R,2))/(27.*tau) - (2*C_5*tau)/R + C_23*I_1((R*lambda_2)/tau) + C_22*K_1((R*lambda_2)/tau);

        values[0] = d_0 + cos(phi) * d;
    }
};

class Velocity : public dolfin::Expression {
    public:
    Velocity() : dolfin::Expression(2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R   = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double a_0 = (C_7*tau)/R;
        double a = C_0 - A_1*R*(-2/3. + (2*std::pow(R,2))/(27.*std::pow(tau,2))) + (C_4*std::pow(R,2))/(2.*std::pow(tau,2)) - (C_3*std::pow(tau,2))/(2.*std::pow(R,2)) + (2*C_14*tau*I_1((R*lambda_1)/tau))/(R*lambda_1) + (2*C_15*tau*K_1((R*lambda_1)/tau))/(R*lambda_1) + C_5*(-1/2. + std::log(R/tau));

        double b_0 = C_10*R + C_13/(2.*R*tau) + (C_11*R)/(3.*std::sqrt(5)*tau) + (-6*C_11*I_1((R*lambda_1)/tau) - 6*C_12*K_1((R*lambda_1)/tau))/15.;
        double b = C_0 - A_1*R*(-1/3. + std::pow(R,2)/(54.*std::pow(tau,2))) + (3*C_4*std::pow(R,2))/(2.*std::pow(tau,2)) + (C_3*std::pow(tau,2))/(2.*std::pow(R,2)) + C_14*(I_0((R*lambda_1)/tau) + I_2((R*lambda_1)/tau)) - C_15*(K_0((R*lambda_1)/tau) + K_2((R*lambda_1)/tau)) + C_5*(1/2. + std::log(R/tau));

        double u_R   = a_0 + cos(phi) * a;
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

        double gamma_0 = (4*C_6*std::pow(tau,2))/(5.*std::pow(R,2)) + (2*C_7*std::pow(tau,2))/std::pow(R,2) + (C_19*tau*I_1((R*lambda_3)/tau))/(R*lambda_3) + C_17*(-(tau*I_1((R*lambda_2)/tau))/(2.*R*lambda_2) - I_2((R*lambda_2)/tau)) + (C_18*tau*K_1((R*lambda_3)/tau))/(R*lambda_3) + C_16*((tau*K_1((R*lambda_2)/tau))/(2.*R*lambda_2) - K_2((R*lambda_2)/tau));
        double gamma = (-2*C_4*R)/tau + (7*A_1*std::pow(R,2))/(27.*tau) - (2*C_5*tau)/R - (4*C_1*std::pow(tau,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(tau,3))/std::pow(R,3) + C_23*(-I_1((R*lambda_2)/tau) + (3*tau*I_2((R*lambda_2)/tau))/(2.*R*lambda_2)) + (C_20*tau*I_2((R*lambda_3)/tau))/(R*lambda_3) + C_22*(-K_1((R*lambda_2)/tau) - (3*tau*K_2((R*lambda_2)/tau))/(2.*R*lambda_2)) + (C_21*tau*K_2((R*lambda_3)/tau))/(R*lambda_3);

        double omega_0 = (4*C_6*std::pow(tau,2))/(5.*std::pow(R,2)) + (2*C_7*std::pow(tau,2))/std::pow(R,2) + C_17*(-I_0((R*lambda_2)/tau)/2. + (3*tau*I_1((R*lambda_2)/tau))/(2.*R*lambda_2)) + C_19*(-I_0((R*lambda_3)/tau) + (tau*I_1((R*lambda_3)/tau))/(R*lambda_3)) + C_16*(-K_0((R*lambda_2)/tau)/2. - (3*tau*K_1((R*lambda_2)/tau))/(2.*R*lambda_2)) + C_18*(K_0((R*lambda_3)/tau) + (tau*K_1((R*lambda_3)/tau))/(R*lambda_3));
        double omega = (-2*C_4*R)/tau + (2*A_1*std::pow(R,2))/(27.*tau) - (2*C_5*tau)/R - (4*C_1*std::pow(tau,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(tau,3))/std::pow(R,3) + C_20*(-I_1((R*lambda_3)/tau) + (tau*I_2((R*lambda_3)/tau))/(R*lambda_3)) + C_23*(-(tau*I_2((R*lambda_2)/tau))/(2.*R*lambda_2) - I_3((R*lambda_2)/tau)/2.) + C_21*(K_1((R*lambda_3)/tau) + (tau*K_2((R*lambda_3)/tau))/(R*lambda_3)) + C_22*((tau*K_2((R*lambda_2)/tau))/(2.*R*lambda_2) - K_3((R*lambda_2)/tau)/2.);

        double kappa_0 = C_13/std::pow(R,2);
        double kappa = (2*C_4*R)/tau - (A_1*std::pow(R,2))/(9.*tau) - (4*C_1*std::pow(tau,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(tau,3))/std::pow(R,3) + (3*C_23*tau*I_2((R*lambda_2)/tau))/(2.*R*lambda_2) + (C_20*tau*I_2((R*lambda_3)/tau))/(R*lambda_3) - (3*C_22*tau*K_2((R*lambda_2)/tau))/(2.*R*lambda_2) + (C_21*tau*K_2((R*lambda_3)/tau))/(R*lambda_3);

        double sigma_RR     =   gamma_0 + cos(phi) * gamma;
        double sigma_Rphi   =   kappa_0 + sin(phi) * kappa;
        double sigma_phiphi = -(omega_0 + cos(phi) * omega);

        // Heinz Schade p377: Assume symmetry of tensor <=> a_Rphi = a_phiR
        double sigma_xx =
            + sigma_RR * pow(cos(phi),2)
            - (sigma_Rphi + sigma_Rphi) * cos(phi) * sin(phi)
            + sigma_phiphi * pow(sin(phi),2);
        double sigma_xy =
            + sigma_Rphi * pow(cos(phi),2)
            + (sigma_RR - sigma_phiphi) * cos(phi) * sin(phi)
            - sigma_Rphi * pow(sin(phi),2);
        double sigma_yy =
            + sigma_phiphi * pow(cos(phi),2)
            + (sigma_Rphi + sigma_Rphi) * cos(phi) * sin(phi)
            + sigma_RR * pow(sin(phi),2);

        values[0] = sigma_xx ;
        values[1] = sigma_xy ;
        values[2] = sigma_yy ;
        values[3] = sigma_xy ; // no difference because symmetry
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