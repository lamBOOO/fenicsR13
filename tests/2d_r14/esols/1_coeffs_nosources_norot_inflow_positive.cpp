#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

double kn = 1.0;

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
double C_0 = 0.;
double C_1 = -1.013912644124032;
double C_2 = -1.2677166806836977;
double C_3 = 0.04992492338715168;
double C_4 = 1.450672783912046;
double C_5 = 0.;
double C_6 = -0.5758748079393637;
double C_7 = 0.10460960523903531;
double C_8 = -0.08439250784359258;
double C_9 = -0.05818511245891174;
double C_10 = 0.044231601637359494;
double C_11 = -0.12041762247811193;
double CK_1 = -0.12330240809678536;
double CK_2 = 2.672949922828857;
double CK_3 = 0.01079637736731986;
double CK_4 = -0.049087424652564535;
double CK_5 = -0.0016650371020496696;
double CK_6 = 0.6120973059811231;
double CK_7 = 0.018332699081602767;
double CK_8 = 0.049452314164954506;
double CK_9 = -0.00002172097726484143;
double CK_10 = 0.27405520113895676;
double CK_11 = 0.0018031130578746815;
double CK_12 = 0.06975633394673138;
double CK_13 = 0.;
double CK_14 = 0.;
double CK_15 = -0.017056297688993596;
double CK_16 = -0.12307639475900282;

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

double lambda_1 = 0.6546536707079771;
double lambda_2 = 0.9128709291752769;
double lambda_3 = 1.224744871391589;
double lambda_4 = 0.7453559924999299;

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

        double c_0 = C_4 - (CK_1*I_0((lambda_1*R)/kn))/15. - (8*CK_5*I_0((lambda_2*R)/kn))/15. - (CK_2*K_0((lambda_1*R)/kn))/15. - (8*CK_6*K_0((lambda_2*R)/kn))/15. - (4*C_2*log(R/kn))/15.;
        double c = (C_9*kn)/R + (C_10*R)/kn - (CK_3*I_1((lambda_1*R)/kn))/15. - (16*CK_7*I_1((lambda_2*R)/kn))/15. - (CK_4*K_1((lambda_1*R)/kn))/15. - (16*CK_8*K_1((lambda_2*R)/kn))/15.;

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

        double alpha_0 = (C_2*kn)/R;
        double alpha = (-15*C_10)/4. - 3*C_8 - (1.5*C_7*pow(kn,2))/pow(R,2) + (3.75*C_9*pow(kn,2))/pow(R,2) + CK_15*(I_0((lambda_4*R)/kn) - I_2((lambda_4*R)/kn)) + CK_16*(K_0((lambda_4*R)/kn) - K_2((lambda_4*R)/kn));

        double beta_0 = CK_13*I_1((lambda_4*R)/kn) + CK_14*K_1((lambda_4*R)/kn);
        double beta = (-15*C_10)/4. - 3*C_8 + (1.5*C_7*pow(kn,2))/pow(R,2) - (3.75*C_9*pow(kn,2))/pow(R,2) + CK_15*(I_0((lambda_4*R)/kn) + I_2((lambda_4*R)/kn)) + CK_16*(K_0((lambda_4*R)/kn) + K_2((lambda_4*R)/kn));

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

        double d_0 = C_9 + CK_4*I_0((R*lambda_2)/kn) + CK_3*K_0((R*lambda_2)/kn);
        double d = (4*C_4*R)/kn - (2*C_5*kn)/R + CK_10*I_1((R*lambda_2)/kn) + CK_9*K_1((R*lambda_2)/kn);

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

        double a_0 = (0.5*C_1*kn)/R - (0.4*C_2*kn)/R;
        double a = C_11 + (0.25*C_6*pow(kn,2))/pow(R,2) - (1.5*C_9*pow(kn,2))/pow(R,2) - (0.25*C_8*pow(R,2))/pow(kn,2) + CK_15*((-2*I_0((lambda_4*R)/kn))/5. + (2*I_2((lambda_4*R)/kn))/5.) + CK_16*((-2*K_0((lambda_4*R)/kn))/5. + (2*K_2((lambda_4*R)/kn))/5.) + C_7*(0.5 + (1.6666666666666667*pow(kn,2))/pow(R,2) - log(R/kn)/2.);

        double b_0 = (0.5*C_0*kn)/R + (C_5*R)/kn - (2*CK_13*I_1((lambda_4*R)/kn))/5. - (2*CK_14*K_1((lambda_4*R)/kn))/5.;
        double b = C_11 - (0.25*C_6*pow(kn,2))/pow(R,2) + (1.5*C_9*pow(kn,2))/pow(R,2) - (0.75*C_8*pow(R,2))/pow(kn,2) + CK_15*((-2*I_0((lambda_4*R)/kn))/5. - (2*I_2((lambda_4*R)/kn))/5.) + CK_16*((-2*K_0((lambda_4*R)/kn))/5. - (2*K_2((lambda_4*R)/kn))/5.) + C_7*((-1.6666666666666667*pow(kn,2))/pow(R,2) - log(R/kn)/2.);

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

        double gamma_0 = (C_1*pow(kn,2))/pow(R,2) + CK_5*(I_0((lambda_2*R)/kn)/3. + I_2((lambda_2*R)/kn)) + CK_9*(-I_0((lambda_3*R)/kn) + I_2((lambda_3*R)/kn)) + CK_6*(K_0((lambda_2*R)/kn)/3. + K_2((lambda_2*R)/kn)) + CK_10*(-K_0((lambda_3*R)/kn) + K_2((lambda_3*R)/kn));
        double gamma = (C_6*pow(kn,3))/pow(R,3) + (C_7*kn)/R + (C_8*R)/kn + CK_7*((5*I_1((lambda_2*R)/kn))/3. + I_3((lambda_2*R)/kn)) + CK_11*(-I_1((lambda_3*R)/kn) + I_3((lambda_3*R)/kn)) + CK_8*((5*K_1((lambda_2*R)/kn))/3. + K_3((lambda_2*R)/kn)) + CK_12*(-K_1((lambda_3*R)/kn) + K_3((lambda_3*R)/kn));

        double kappa_0 = (C_0*pow(kn,2))/pow(R,2);
        double kappa = (C_6*pow(kn,3))/pow(R,3) - (C_8*R)/kn + CK_7*(-I_1((lambda_2*R)/kn) + I_3((lambda_2*R)/kn)) + CK_11*(-I_1((lambda_3*R)/kn) + I_3((lambda_3*R)/kn)) + CK_8*(-K_1((lambda_2*R)/kn) + K_3((lambda_2*R)/kn)) + CK_12*(-K_1((lambda_3*R)/kn) + K_3((lambda_3*R)/kn));

        double omega_0 = (C_1*pow(kn,2))/pow(R,2) + CK_5*(-0.3333333333333333*I_0((lambda_2*R)/kn) + I_2((lambda_2*R)/kn)) + CK_9*(I_0((lambda_3*R)/kn) + I_2((lambda_3*R)/kn)) + CK_6*(-0.3333333333333333*K_0((lambda_2*R)/kn) + K_2((lambda_2*R)/kn)) + CK_10*(K_0((lambda_3*R)/kn) + K_2((lambda_3*R)/kn));
        double omega = (C_6*pow(kn,3))/pow(R,3) + (C_7*kn)/R + (C_8*R)/kn + CK_7*(I_1((lambda_2*R)/kn)/3. + I_3((lambda_2*R)/kn)) + CK_11*(3*I_1((lambda_3*R)/kn) + I_3((lambda_3*R)/kn)) + CK_8*(K_1((lambda_2*R)/kn)/3. + K_3((lambda_2*R)/kn)) + CK_12*(3*K_1((lambda_3*R)/kn) + K_3((lambda_3*R)/kn));

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
        values[2] = sigma_xy ;
        values[3] = sigma_yy ;
    }
};

class ScriptR : public dolfin::Expression {
    public:
    ScriptR() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override {

        double R = sqrt(pow(x[0],2)+pow(x[1],2));
        double phi = atan2(x[1],x[0]);

        double e_0 = CK_1*I_0((lambda_1*R)/kn) + CK_2*K_0((lambda_1*R)/kn);
        double e = CK_3*I_1((lambda_1*R)/kn) + CK_4*K_1((lambda_1*R)/kn);

        double scriptR = e_0 + cos(phi) * e;

        values[0] = scriptR;
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

    py::class_<ScriptR, std::shared_ptr<ScriptR>, dolfin::Expression>
        (m, "ScriptR")
    .def(py::init<>());
}
