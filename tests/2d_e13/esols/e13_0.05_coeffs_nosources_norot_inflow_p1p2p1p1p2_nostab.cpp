// generated with:
// https://github.com/lamBOOO/Analytical-Solution-Generating-Script-for-Lin-R13

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

        double x = z[0];
        double y = z[1];
        using boost::math::cyl_bessel_i;
        using boost::math::cyl_bessel_k;
        const double r   = std::sqrt(x*x + y*y);

        values[0] = 1.1475416669755338 - 1.4136499460087132e-19/std::pow(r,5) - \
                    1.9855291553960966e-18/std::pow(r,4) - 8.503184547098384e-18/std::pow(r,3) - \
                    7.51454581320324e-16/std::pow(r,2) - 0.006824199418861071/r + \
                    0.009607734992218002*r + 1.1347378340243468e-9*std::pow(r,2) - \
                    2.861003252615381e-7*std::pow(r,3) + 7.835105240179381e-8*std::pow(r,4) - \
                    3.2050302992908996e-6*std::pow(r,5) + \
                    2.1335502466342462e-17*cyl_bessel_i(0,18.250101513540354*r) - \
                    2.3955884082321502e-17*cyl_bessel_i(1,18.250101513540354*r) - \
                    7.515152743849166e6*cyl_bessel_k(0,18.250101513540354*r) - \
                    1.7297773670560922e6*cyl_bessel_k(1,18.250101513540354*r) + \
                    1.1304068783278136*std::log(r);
    }
};

class Heatflux : public dolfin::Expression {
    public:
    Heatflux() : dolfin::Expression(2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];

        values[0] = -0.18517279447922175*x/(pow(x, 2) + pow(y, 2)) - y*(-0.0046391317838279369*BI(1, 3.333333333333333*sqrt(5)*sqrt(pow(x, 2) + pow(y, 2)))/BI(1, 6.6666666666666661*sqrt(5)) + 0.060193212257875592*BK(1, 3.333333333333333*sqrt(5)*sqrt(pow(x, 2) + pow(y, 2)))/BK(1, 1.6666666666666665*sqrt(5)))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[1] = x*(-0.0046391317838279369*BI(1, 3.333333333333333*sqrt(5)*sqrt(pow(x, 2) + pow(y, 2)))/BI(1, 6.6666666666666661*sqrt(5)) + 0.060193212257875592*BK(1, 3.333333333333333*sqrt(5)*sqrt(pow(x, 2) + pow(y, 2)))/BK(1, 1.6666666666666665*sqrt(5)))/sqrt(pow(x, 2) + pow(y, 2)) - 0.18517279447922175*y/(pow(x, 2) + pow(y, 2)) ;
    }
};

class Pressure : public dolfin::Expression {
    public:
    Pressure() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];
        using boost::math::cyl_bessel_i;
        using boost::math::cyl_bessel_k;
        const double r   = std::sqrt(x*x + y*y);

        values[0] = 0.0018993199772926886 + 2.196493140489761e-19/std::pow(r,5) - \
                    1.820618508617261e-17/std::pow(r,4) + 1.7806835163503867e-16/std::pow(r,3) + \
                    1.4497382343428474e-15/std::pow(r,2) - 0.6945985178625124/r - \
                    0.2277552405043594*r + 2.2637072826466894e-9*std::pow(r,2) - \
                    2.033833348503882e-7*std::pow(r,3) + 8.411870914721436e-7*std::pow(r,4) - \
                    6.155955617129015e-6*std::pow(r,5) + \
                    1.7810877458902648e-17*cyl_bessel_i(0,18.250101513540354*r) - \
                    1.999837203192199e-17*cyl_bessel_i(1,18.250101513540354*r) - \
                    6.273649510565284e6*cyl_bessel_k(0,18.250101513540354*r) - \
                    1.4440181460184255e6*cyl_bessel_k(1,18.250101513540354*r) - \
                    2.0537662284439946e-12*std::log(r);
        }
};

class Velocity : public dolfin::Expression {
    public:
    Velocity() : dolfin::Expression(2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];
        using boost::math::cyl_bessel_i;
        using boost::math::cyl_bessel_k;
        const double r   = std::sqrt(x*x + y*y);

        values[0] = -2.0447636553038926 + 1.0359201653139577e-19/std::pow(r,5) + \
                    2.0767426137201643e-18/std::pow(r,4) + 3.1555921783409397e-16/std::pow(r,3) \
                    + 2.608427543187341/std::pow(r,2) - 2.270488285671194e-6/r + \
                    3.4283192024707477e-12*r - 0.5693881120084758*std::pow(r,2) - \
                    2.135570368091424e-9*std::pow(r,3) - 1.0672365177466179e-7*std::pow(r,4) + \
                    4.480857814008834e-6*std::pow(r,5) - \
                    2.8359656347824803e-14*cyl_bessel_i(0,14.920428347746274*r) + \
                    2.8359656347824803e-14*cyl_bessel_i(2,14.920428347746274*r) - \
                    388304.5478035767*cyl_bessel_k(0,14.920428347746274*r) + \
                    388304.5478035762*cyl_bessel_k(2,14.920428347746274*r) + \
                    6.945985178621474*std::log(r);
        values[1] = -9.029638580084187e-12 - 4.7257426603090694e-20/std::pow(r,5) + \
                    1.2549535982543759e-18/std::pow(r,4) - 1.3310760046952658e-16/std::pow(r,3) \
                    + 4.908584471553888e-15/std::pow(r,2) + 1.1114114849250902/r - \
                    0.2592858292136669*r + 8.46678012321097e-11*std::pow(r,2) + \
                    6.313585336147986e-9*std::pow(r,3) - 3.827994226840259e-8*std::pow(r,4) + \
                    2.8097142147287375e-7*std::pow(r,5) - \
                    7.033822365876907e-15*cyl_bessel_i(1,14.920428347746274*r) + \
                    155945.4280877491*cyl_bessel_k(1,14.920428347746274*r) - \
                    1.7470036479589143e-12*std::log(r);
    }
};

class Stress : public dolfin::Expression {
    public:
    Stress() : dolfin::Expression(2,2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];

        values[0] = x*(x*(-1.8458676902269367e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 5.5376030706808106e-10*BI(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 6.9689048822689771e-14*BI(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 2.0458491401654371*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 6.1375474204963112*BK(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 4.6692102079126236*BK(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 0.014813819184464139/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + 0.074449891379848129*y/pow(pow(x, 2) + pow(y, 2), 3.0/2.0))/sqrt(pow(x, 2) + pow(y, 2)) - y*(-0.074449891379848129*x/pow(pow(x, 2) + pow(y, 2), 3.0/2.0) - y*(-1.8458676902269367e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 5.5376030706808106e-10*BI(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 2.0458491401654371*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 6.1375474204963112*BK(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 0.014813819184464139/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[1] = x*(-0.074449891379848129*x/pow(pow(x, 2) + pow(y, 2), 3.0/2.0) - y*(-1.8458676902269367e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 5.5376030706808106e-10*BI(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 2.0458491401654371*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 6.1375474204963112*BK(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 0.014813819184464139/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + y*(x*(-1.8458676902269367e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 5.5376030706808106e-10*BI(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 6.9689048822689771e-14*BI(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 2.0458491401654371*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 6.1375474204963112*BK(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 4.6692102079126236*BK(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 0.014813819184464139/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + 0.074449891379848129*y/pow(pow(x, 2) + pow(y, 2), 3.0/2.0))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[2] = x*(-0.074449891379848129*x/pow(pow(x, 2) + pow(y, 2), 3.0/2.0) - y*(-1.8458676902269367e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 5.5376030706808106e-10*BI(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 2.0458491401654371*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 6.1375474204963112*BK(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 0.014813819184464139/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + y*(x*(-1.8458676902269367e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 5.5376030706808106e-10*BI(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 6.9689048822689771e-14*BI(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 2.0458491401654371*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 6.1375474204963112*BK(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 4.6692102079126236*BK(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 0.014813819184464139/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + 0.074449891379848129*y/pow(pow(x, 2) + pow(y, 2), 3.0/2.0))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[3] = x*(x*(-1.8458676902269367e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 5.5376030706808106e-10*BI(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 2.0458491401654371*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 6.1375474204963112*BK(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 0.014813819184464139/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) - 0.074449891379848129*y/pow(pow(x, 2) + pow(y, 2), 3.0/2.0))/sqrt(pow(x, 2) + pow(y, 2)) + y*(-0.074449891379848129*x/pow(pow(x, 2) + pow(y, 2), 3.0/2.0) + y*(-1.8458676902269367e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 6.9689048822689771e-14*BI(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 5.5376030706808106e-10*BI(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 6.9689048822689771e-14*BI(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 2.0458491401654371*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 4.6692102079126236*BK(0, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) + 6.1375474204963112*BK(2, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 4.6692102079126236*BK(2, 12.24744871391589*sqrt(pow(x, 2) + pow(y, 2))) - 0.014813819184464139/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) ;
    }
};

PYBIND11_MODULE(SIGNATURE, m) {

    py::class_<Temperature, std::shared_ptr<Temperature>, dolfin::Expression>
        (m,"Temperature")
    .def(py::init<>());

    py::class_<Heatflux, std::shared_ptr<Heatflux>, dolfin::Expression>
        (m,"Heatflux")
    .def(py::init<>());

    py::class_<Pressure, std::shared_ptr<Pressure>, dolfin::Expression>
        (m,"Pressure")
    .def(py::init<>());

    py::class_<Velocity, std::shared_ptr<Velocity>, dolfin::Expression>
        (m,"Velocity")
    .def(py::init<>());

    py::class_<Stress, std::shared_ptr<Stress>, dolfin::Expression>
        (m,"Stress")
    .def(py::init<>());

}
