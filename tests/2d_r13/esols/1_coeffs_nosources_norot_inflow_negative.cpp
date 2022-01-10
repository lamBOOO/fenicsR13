#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

double kn = 1.0;

// // Constants eps=10^-5
// double C_0 = 1.337436195119152;
// double C_1 = -0.6064772869609508;
// double C_2 = -0.1726980022374072;
// double C_3 = 0.7901754852718373;
// double C_4 = -0.08396168002205009;
// double C_5 = 0.1252124495285327;
// double C_6 = -0.1867020765534401;
// double C_7 = 4.568817251907919e-7;
// double C_8 = 1.963244517881108;
// double C_9 = -0.03843441110894395;
// double CK_1 = -0.03180963566002944;
// double CK_2 = 0.03411564204278963;
// double CK_3 = -0.07221148250544122;
// double CK_4 = 0.01263372282796629;
// double CK_5 = 0.03846898363886858;
// double CK_6 = 0.003998656953005047;
// double CK_7 = 0.1003792356748235;
// double CK_8 = 0.05030863773067603;
// double CK_9 = -0.001651068357981648;
// double CK_10 = 0.3495787712905514;

// // Constants eps=10^-4
// double C_0 = 1.337429881309671;
// double C_1 = -0.6064620302518795;
// double C_2 = -0.17269389293132;
// double C_3 = 0.7901524616369147;
// double C_4 = -0.0839597704511531;
// double C_5 = 0.1252094721825813;
// double C_6 = -0.1867023121536032;
// double C_7 = 4.568671844035206e-6;
// double C_8 = 1.963244041101796;
// double C_9 = -0.03843414811236764;
// double CK_1 = -0.03180886812453335;
// double CK_2 = 0.03411478821870104;
// double CK_3 = -0.0722079455756748;
// double CK_4 = 0.01263375346771417;
// double CK_5 = 0.03846606618426868;
// double CK_6 = 0.003998605916734596;
// double CK_7 = 0.1003768640953838;
// double CK_8 = 0.05030581411268845;
// double CK_9 = -0.001649480514012068;
// double CK_10 = 0.3495706866039064;

// // Constants eps=10^-3
double C_0 = 1.337366717236668;
double C_1 = -0.6063094894924456;
double C_2 = -0.1726528012092394;
double C_3 = 0.7899222753828902;
double C_4 = -0.08394067339153539;
double C_5 = 0.1251797046222099;
double C_6 = -0.1867046673301068;
double C_7 = 0.00004567217257213145;
double C_8 = 1.963239274978468;
double C_9 = -0.03843149872072692;
double CK_1 = -0.0318011928378749;
double CK_2 = 0.03410625144135048;
double CK_3 = -0.07217258866517415;
double CK_4 = 0.012634059757885;
double CK_5 = 0.03843690185588499;
double CK_6 = 0.003998095732771025;
double CK_7 = 0.1003531476907441;
double CK_8 = 0.05027758648967144;
double CK_9 = -0.00163360750417993;
double CK_10 = 0.3494898346560165;

// // Constants eps=10^-2
// double C_0 = 1.336732533609719;
// double C_1 = -0.6047867380273104;
// double C_2 = -0.1722420311465265;
// double C_3 = 0.7876254443155443;
// double C_4 = -0.08374957633036394;
// double C_5 = 0.1248826224847917;
// double C_6 = -0.1867281362946452;
// double C_7 = 0.0004552621117990246;
// double C_8 = 1.963191781306897;
// double C_9 = -0.03840306975572652;
// double CK_1 = -0.03172444950554837;
// double CK_2 = 0.03402103132100606;
// double CK_3 = -0.07182026259642857;
// double CK_4 = 0.01263711189140766;
// double CK_5 = 0.03814628389686272;
// double CK_6 = 0.003993011829577151;
// double CK_7 = 0.1001159319865724;
// double CK_8 = 0.04999616566544083;
// double CK_9 = -0.001475419473718654;
// double CK_10 = 0.3486808426586021;

// // Constants eps=10^-1
// double C_0 = 1.330185932037747;
// double C_1 = -0.5898443496048875;
// double C_2 = -0.1681607691585836;
// double C_3 = 0.7651780840964286;
// double C_4 = -0.08183362422804108;
// double C_5 = 0.1219741008236003;
// double C_6 = -0.1869542885544411;
// double C_7 = 0.004402163569437981;
// double C_8 = 1.962734121526709;
// double C_9 = -0.0379330049303497;
// double CK_1 = -0.03096036349825703;
// double CK_2 = 0.03318469878613935;
// double CK_3 = -0.06842516876936609;
// double CK_4 = 0.01266652294157217;
// double CK_5 = 0.03534582339667038;
// double CK_6 = 0.003944022186510511;
// double CK_7 = 0.09774692911694641;
// double CK_8 = 0.04726696889230132;
// double CK_9 = 0.0000532462109035967;
// double CK_10 = 0.3405754068324816;

// // Constants eps=10^-0
// double C_0 = 1.26608703305865;
// double C_1 = -0.4723447945692413;
// double C_2 = -0.1340987755805701;
// double C_3 = 0.5922083059862845;
// double C_4 = -0.06517741879346212;
// double C_5 = 0.09936322021509666;
// double C_6 = -0.1883345440537925;
// double C_7 = 0.02849094538051043;
// double C_8 = 1.959940926229848;
// double C_9 = -0.02236820885608115;
// double CK_1 = -0.02452200335350193;
// double CK_2 = 0.02660474893156432;
// double CK_3 = -0.04770418667594679;
// double CK_4 = 0.01284602485664462;
// double CK_5 = 0.01825401538380585;
// double CK_6 = 0.003645027941151212;
// double CK_7 = 0.07750866083895219;
// double CK_8 = 0.02706654817009965;
// double CK_9 = 0.01115352331657436;
// double CK_10 = 0.2703240873448836;

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

        double c_0 = C_8 + (2*CK_4*I_0((R*lambda_2)/kn))/5. + (2*CK_3*K_0((R*lambda_2)/kn))/5. - (4*C_6*std::log(R/kn))/15.;
        double c = (-4*C_2*R)/(15.*kn) + (8*C_4*R)/(5.*kn) - (2*C_1*kn)/(15.*R) - (4*C_5*kn)/(5.*R) + (2*CK_10*I_1((R*lambda_2)/kn))/5. + (2*CK_9*K_1((R*lambda_2)/kn))/5.;

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

        double alpha_0 = (C_6*kn)/R;
        double alpha = C_2 - (C_1*std::pow(kn,2))/(2.*std::pow(R,2)) - (5*((2*CK_1*kn*I_1((R*lambda_1)/kn))/(R*lambda_1) + (2*CK_2*kn*K_1((R*lambda_1)/kn))/(R*lambda_1)))/2.;

        double beta_0 = 0;
        double beta = C_2 + (C_1*std::pow(kn,2))/(2.*std::pow(R,2)) - (5*(CK_1*(I_0((R*lambda_1)/kn) + I_2((R*lambda_1)/kn)) - CK_2*(K_0((R*lambda_1)/kn) + K_2((R*lambda_1)/kn))))/2.;

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

        double a_0 = (C_7*kn)/R;
        double a = C_0 + (C_4*std::pow(R,2))/(2.*std::pow(kn,2)) - (C_3*std::pow(kn,2))/(2.*std::pow(R,2)) + (2*CK_1*kn*I_1((R*lambda_1)/kn))/(R*lambda_1) + (2*CK_2*kn*K_1((R*lambda_1)/kn))/(R*lambda_1) + C_5*(-0.5 + std::log(R/kn));

        double b_0 = 0;
        double b = C_0 + (3*C_4*std::pow(R,2))/(2.*std::pow(kn,2)) + (C_3*std::pow(kn,2))/(2.*std::pow(R,2)) + CK_1*(I_0((R*lambda_1)/kn) + I_2((R*lambda_1)/kn)) - CK_2*(K_0((R*lambda_1)/kn) + K_2((R*lambda_1)/kn)) + C_5*(0.5 + std::log(R/kn));

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

        double gamma_0 = (4*C_6*std::pow(kn,2))/(5.*std::pow(R,2)) + (2*C_7*std::pow(kn,2))/std::pow(R,2) + (CK_6*kn*I_1((R*lambda_3)/kn))/(R*lambda_3) + CK_4*(-(kn*I_1((R*lambda_2)/kn))/(2.*R*lambda_2) - I_2((R*lambda_2)/kn)) + (CK_5*kn*K_1((R*lambda_3)/kn))/(R*lambda_3) + CK_3*((kn*K_1((R*lambda_2)/kn))/(2.*R*lambda_2) - K_2((R*lambda_2)/kn));
        double gamma = (-2*C_4*R)/kn - (2*C_5*kn)/R - (4*C_1*std::pow(kn,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(kn,3))/std::pow(R,3) + CK_10*(-I_1((R*lambda_2)/kn) + (3*kn*I_2((R*lambda_2)/kn))/(2.*R*lambda_2)) + (CK_7*kn*I_2((R*lambda_3)/kn))/(R*lambda_3) + CK_9*(-K_1((R*lambda_2)/kn) - (3*kn*K_2((R*lambda_2)/kn))/(2.*R*lambda_2)) + (CK_8*kn*K_2((R*lambda_3)/kn))/(R*lambda_3);

        double omega_0 = (4*C_6*std::pow(kn,2))/(5.*std::pow(R,2)) + (2*C_7*std::pow(kn,2))/std::pow(R,2) + CK_4*(-I_0((R*lambda_2)/kn)/2. + (3*kn*I_1((R*lambda_2)/kn))/(2.*R*lambda_2)) + CK_6*(-I_0((R*lambda_3)/kn) + (kn*I_1((R*lambda_3)/kn))/(R*lambda_3)) + CK_3*(-K_0((R*lambda_2)/kn)/2. - (3*kn*K_1((R*lambda_2)/kn))/(2.*R*lambda_2)) + CK_5*(K_0((R*lambda_3)/kn) + (kn*K_1((R*lambda_3)/kn))/(R*lambda_3));
        double omega = (-2*C_4*R)/kn - (2*C_5*kn)/R - (4*C_1*std::pow(kn,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(kn,3))/std::pow(R,3) + CK_7*(-I_1((R*lambda_3)/kn) + (kn*I_2((R*lambda_3)/kn))/(R*lambda_3)) + CK_10*(-(kn*I_2((R*lambda_2)/kn))/(2.*R*lambda_2) - I_3((R*lambda_2)/kn)/2.) + CK_8*(K_1((R*lambda_3)/kn) + (kn*K_2((R*lambda_3)/kn))/(R*lambda_3)) + CK_9*((kn*K_2((R*lambda_2)/kn))/(2.*R*lambda_2) - K_3((R*lambda_2)/kn)/2.);

        double kappa_0 = 0;
        double kappa = (2*C_4*R)/kn - (4*C_1*std::pow(kn,3))/(5.*std::pow(R,3)) - (2*(C_3 - (64*C_5)/15.)*std::pow(kn,3))/std::pow(R,3) + (3*CK_10*kn*I_2((R*lambda_2)/kn))/(2.*R*lambda_2) + (CK_7*kn*I_2((R*lambda_3)/kn))/(R*lambda_3) - (3*CK_9*kn*K_2((R*lambda_2)/kn))/(2.*R*lambda_2) + (CK_8*kn*K_2((R*lambda_3)/kn))/(R*lambda_3);

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
