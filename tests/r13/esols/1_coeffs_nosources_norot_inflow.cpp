#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

double tau = 1.0;

// // Constants eps=10^-5
// double C_0 = 1.088418574133493;
// double C_1 = -0.4708784837854026;
// double C_2 = -0.1060780078877844;
// double C_3 = 0.6397363738946663;
// double C_4 = -0.04268373030990083;
// double C_5 = 0.1027291801998076;
// double C_6 = -0.1867020667903523;
// double C_7 = 2.864923360570221e-7;
// double C_8 = 1.963244537638472;
// double C_9 = 0.007704644954804194;
// double CK_1 = -0.01841555365364274;
// double CK_2 = 0.02647436311982035;
// double CK_3 = -0.0722116290730638;
// double CK_4 = 0.01263372155827893;
// double CK_5 = 0.03846910453592139;
// double CK_6 = 0.003998659067908465;
// double CK_7 = 0.05469703377538136;
// double CK_8 = 0.04448263813229264;
// double CK_9 = -0.001607567940383826;
// double CK_10 = 0.1781491089766695;

// // Constants eps=10^-4
// double C_0 = 1.088458465830652;
// double C_1 = -0.4708795265589713;
// double C_2 = -0.1060844454754048;
// double C_3 = 0.6397268622237638;
// double C_4 = -0.0426888328939176;
// double C_5 = 0.1027285697835575;
// double C_6 = -0.1867022145254646;
// double C_7 = 2.864825757859002e-6;
// double C_8 = 1.963244238669897;
// double C_9 = 0.007703775736883896;
// double CK_1 = -0.01841695056397639;
// double CK_2 = 0.02647443224729624;
// double CK_3 = -0.07220941121077909;
// double CK_4 = 0.01263374077119675;
// double CK_5 = 0.03846727512087732;
// double CK_6 = 0.003998627065175332;
// double CK_7 = 0.0547022288973836;
// double CK_8 = 0.04447888927384885;
// double CK_9 = -0.001604798505932295;
// double CK_10 = 0.178170030606713;

// // Constants eps=10^-3
double C_0 = 1.088856712167725;
double C_1 = -0.4708896102041703;
double C_2 = -0.1061486462553457;
double C_3 = 0.6396313746480878;
double C_4 = -0.04273974910778836;
double C_5 = 0.102722409367338;
double C_6 = -0.1867036913227729;
double C_7 = 0.00002863849458707334;
double C_8 = 1.963241250104882;
double C_9 = 0.007695099574563038;
double CK_1 = -0.01843088435918791;
double CK_2 = 0.02647510412125701;
double CK_3 = -0.07218724090202373;
double CK_4 = 0.01263393282835227;
double CK_5 = 0.03844898782835787;
double CK_6 = 0.003998307157812554;
double CK_7 = 0.05475405925519971;
double CK_8 = 0.04444138977514516;
double CK_9 = -0.001577106876225382;
double CK_10 = 0.178378791898016;

// // Constants eps=10^-2
// double C_0 = 1.09277278324226;
// double C_1 = -0.4709565408029582;
// double C_2 = -0.1067733521869535;
// double C_3 = 0.6386400415463009;
// double C_4 = -0.04323806992636355;
// double C_5 = 0.1026552684836559;
// double C_6 = -0.1867184037607001;
// double C_7 = 0.0002854059606868768;
// double C_8 = 1.963211476840119;
// double C_9 = 0.007609935703161957;
// double CK_1 = -0.0185667325855968;
// double CK_2 = 0.02647991106816222;
// double CK_3 = -0.07196637153183492;
// double CK_4 = 0.0126358461775725;
// double CK_5 = 0.03826680259919404;
// double CK_6 = 0.003995120114343059;
// double CK_7 = 0.05526041387534725;
// double CK_8 = 0.04406535024733006;
// double CK_9 = -0.001300480668267124;
// double CK_10 = 0.180421409499547;

// // Constants eps=10^-1
// double C_0 = 1.125921604594353;
// double C_1 = -0.4686941217654443;
// double C_2 = -0.1114820681020402;
// double C_3 = 0.6256505450040041;
// double C_4 = -0.04724969749859539;
// double C_5 = 0.1015107112011764;
// double C_6 = -0.1868598401676765;
// double C_7 = 0.002753811809744579;
// double C_8 = 1.962925254827631;
// double C_9 = 0.006913904797024802;
// double CK_1 = -0.01961423092232162;
// double CK_2 = 0.0263626002451694;
// double CK_3 = -0.06984306808715408;
// double CK_4 = 0.01265423995026528;
// double CK_5 = 0.03651538489775871;
// double CK_6 = 0.003964481821802226;
// double CK_7 = 0.05925618682225351;
// double CK_8 = 0.04024172550289172;
// double CK_9 = 0.001420920976733996;
// double CK_10 = 0.1968170158183649;

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

        double c_0 = C_8 + (2*CK_4*I_0((R*lambda_2)/tau))/5. + (2*CK_3*K_0((R*lambda_2)/tau))/5. - (4*C_6*std::log(R/tau))/15.;
        double c = (-4*C_2*R)/(15.*tau) + (8*C_4*R)/(5.*tau) - (2*C_1*tau)/(15.*R) - (4*C_5*tau)/(5.*R) + (2*CK_10*I_1((R*lambda_2)/tau))/5. + (2*CK_9*K_1((R*lambda_2)/tau))/5.;

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
        double alpha = C_2 - (C_1*std::pow(tau,2))/(2.*std::pow(R,2)) - (5*((2*CK_1*tau*I_1((R*lambda_1)/tau))/(R*lambda_1) + (2*CK_2*tau*K_1((R*lambda_1)/tau))/(R*lambda_1)))/2.;

        double beta_0 = 0;
        double beta = C_2 + (C_1*std::pow(tau,2))/(2.*std::pow(R,2)) - (5*(CK_1*(I_0((R*lambda_1)/tau) + I_2((R*lambda_1)/tau)) - CK_2*(K_0((R*lambda_1)/tau) + K_2((R*lambda_1)/tau))))/2.;

        double s_R = alpha_0 + cos(phi) * alpha;
        double s_phi = beta_0 - sin(phi) * beta;

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

        double d_0 = C_9 + CK_4*I_0((R*lambda_2)/tau) + CK_3*K_0((R*lambda_2)/tau);
        double d = (4*C_4*R)/tau - (2*C_5*tau)/R + CK_10*I_1((R*lambda_2)/tau) + CK_9*K_1((R*lambda_2)/tau);

        double p = d_0 + cos(phi) * d;

        values[0] = p;
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
        double a = C_0 + (C_4*std::pow(R,2))/(2.*std::pow(tau,2)) - (C_3*std::pow(tau,2))/(2.*std::pow(R,2)) + (2*CK_1*tau*I_1((R*lambda_1)/tau))/(R*lambda_1) + (2*CK_2*tau*K_1((R*lambda_1)/tau))/(R*lambda_1) + C_5*(-0.5 + std::log(R/tau));

        double b_0 = 0;
        double b = C_0 + (3*C_4*std::pow(R,2))/(2.*std::pow(tau,2)) + (C_3*std::pow(tau,2))/(2.*std::pow(R,2)) + CK_1*(I_0((R*lambda_1)/tau) + I_2((R*lambda_1)/tau)) - CK_2*(K_0((R*lambda_1)/tau) + K_2((R*lambda_1)/tau)) + C_5*(0.5 + std::log(R/tau));

        double u_R   = a_0 + cos(phi) * a;
        double u_phi = b_0 - sin(phi) * b;

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

        double gamma_0 = (4*C_6*std::pow(tau,2))/(5.*std::pow(R,2)) + (2*C_7*std::pow(tau,2))/std::pow(R,2) + (CK_6*tau*I_1((R*lambda_3)/tau))/(R*lambda_3) + CK_4*(-(tau*I_1((R*lambda_2)/tau))/(2.*R*lambda_2) - I_2((R*lambda_2)/tau)) + (CK_5*tau*K_1((R*lambda_3)/tau))/(R*lambda_3) + CK_3*((tau*K_1((R*lambda_2)/tau))/(2.*R*lambda_2) - K_2((R*lambda_2)/tau));
        double gamma = (-2*C_4*R)/tau - (2*C_5*tau)/R - (4*C_1*std::pow(tau,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(tau,3))/std::pow(R,3) + CK_10*(-I_1((R*lambda_2)/tau) + (3*tau*I_2((R*lambda_2)/tau))/(2.*R*lambda_2)) + (CK_7*tau*I_2((R*lambda_3)/tau))/(R*lambda_3) + CK_9*(-K_1((R*lambda_2)/tau) - (3*tau*K_2((R*lambda_2)/tau))/(2.*R*lambda_2)) + (CK_8*tau*K_2((R*lambda_3)/tau))/(R*lambda_3);

        double omega_0 = (4*C_6*std::pow(tau,2))/(5.*std::pow(R,2)) + (2*C_7*std::pow(tau,2))/std::pow(R,2) + CK_4*(-I_0((R*lambda_2)/tau)/2. + (3*tau*I_1((R*lambda_2)/tau))/(2.*R*lambda_2)) + CK_6*(-I_0((R*lambda_3)/tau) + (tau*I_1((R*lambda_3)/tau))/(R*lambda_3)) + CK_3*(-K_0((R*lambda_2)/tau)/2. - (3*tau*K_1((R*lambda_2)/tau))/(2.*R*lambda_2)) + CK_5*(K_0((R*lambda_3)/tau) + (tau*K_1((R*lambda_3)/tau))/(R*lambda_3));
        double omega = (-2*C_4*R)/tau - (2*C_5*tau)/R - (4*C_1*std::pow(tau,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(tau,3))/std::pow(R,3) + CK_7*(-I_1((R*lambda_3)/tau) + (tau*I_2((R*lambda_3)/tau))/(R*lambda_3)) + CK_10*(-(tau*I_2((R*lambda_2)/tau))/(2.*R*lambda_2) - I_3((R*lambda_2)/tau)/2.) + CK_8*(K_1((R*lambda_3)/tau) + (tau*K_2((R*lambda_3)/tau))/(R*lambda_3)) + CK_9*((tau*K_2((R*lambda_2)/tau))/(2.*R*lambda_2) - K_3((R*lambda_2)/tau)/2.);

        double kappa_0 = 0;
        double kappa = (2*C_4*R)/tau - (4*C_1*std::pow(tau,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(tau,3))/std::pow(R,3) + (3*CK_10*tau*I_2((R*lambda_2)/tau))/(2.*R*lambda_2) + (CK_7*tau*I_2((R*lambda_3)/tau))/(R*lambda_3) - (3*CK_9*tau*K_2((R*lambda_2)/tau))/(2.*R*lambda_2) + (CK_8*tau*K_2((R*lambda_3)/tau))/(R*lambda_3);

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