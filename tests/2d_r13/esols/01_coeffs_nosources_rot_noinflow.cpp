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

        values[0] = 0.49379411861125794*std::log(10.0*sqrt(pow(x, 2) + pow(y, 2))) + 2.953388304363099e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 3.2733586242646995*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 0.47880449201690833 ;
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

        values[0] = 7.3834707609077468e-10*BI(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) - 8.1833965606617483*BK(0, 9.1287092917527684*sqrt(pow(x, 2) + pow(y, 2))) + 0.00090394958226724855 ;
    }
};

class Velocity : public dolfin::Expression {
    public:
    Velocity() : dolfin::Expression(2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];

        values[0] = 2.1869368024840519e-8*x/(pow(x, 2) + pow(y, 2)) - y*(0.58515350123209375*sqrt(pow(x, 2) + pow(y, 2)) + 0.0018556527135311749*BI(1, 3.333333333333333*sqrt(5)*sqrt(pow(x, 2) + pow(y, 2)))/BI(1, 6.6666666666666661*sqrt(5)) - 0.02407728490315024*BK(1, 3.333333333333333*sqrt(5)*sqrt(pow(x, 2) + pow(y, 2)))/BK(1, 1.6666666666666665*sqrt(5)) - 0.3722494568992406/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[1] = x*(0.58515350123209375*sqrt(pow(x, 2) + pow(y, 2)) + 0.0018556527135311749*BI(1, 3.333333333333333*sqrt(5)*sqrt(pow(x, 2) + pow(y, 2)))/BI(1, 6.6666666666666661*sqrt(5)) - 0.02407728490315024*BK(1, 3.333333333333333*sqrt(5)*sqrt(pow(x, 2) + pow(y, 2)))/BK(1, 1.6666666666666665*sqrt(5)) - 0.3722494568992406/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + 2.1869368024840519e-8*y/(pow(x, 2) + pow(y, 2)) ;
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
