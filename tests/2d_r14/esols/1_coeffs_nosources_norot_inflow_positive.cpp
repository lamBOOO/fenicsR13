#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

namespace py = pybind11;

#include <dolfin/function/Expression.h>

double kn = 1.0;

// Lambdas
double lambda_1 = 0.5345224838248488;
double lambda_2 = 0.9128709291752769;
double lambda_3 = 1.224744871391589;
double lambda_4 = 0.7453559924999299;

// Constants
double C_0 = 0.;
double C_1 = -0.1477752433373134;
double C_2 = -0.18478991652292415;
double C_3 = 0.007616222250050322;
double C_4 = 1.9613523915545736;
double C_5 = 0.;
double C_6 = -0.025990132046856542;
double C_7 = -0.20544333345973814;
double C_8 = 0.08546266837583295;
double C_9 = -0.01965608692371418;
double C_10 = -0.040261109257103515;
double C_11 = 1.1399409455629295;
double CK_1 = -0.12432069916127567;
double CK_2 = 0.24749516983247868;
double CK_3 = -0.02593994988782829;
double CK_4 = -0.006941070313859298;
double CK_5 = -0.009378323117905304;
double CK_6 = 0.05358547329379231;
double CK_7 = -0.06686833734375043;
double CK_8 = 0.0005915756434353135;
double CK_9 = -0.0019786615912372857;
double CK_10 = 0.019027436471274847;
double CK_11 = -0.013685132661250488;
double CK_12 = 0.011110746687159737;
double CK_13 = 0.;
double CK_14 = 0.;
double CK_15 = 0.045818412638178437;
double CK_16 = 0.06590993098932302;

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

        // Temperature
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

        double d_0 = C_3 - (4*CK_5*I_0((lambda_2*R)/kn))/3. - (4*CK_6*K_0((lambda_2*R)/kn))/3.;
        double d = (C_7*kn)/R - (2.*C_8*R)/kn - (8*CK_7*I_1((lambda_2*R)/kn))/3. - (8*CK_8*K_1((lambda_2*R)/kn))/3.;

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
