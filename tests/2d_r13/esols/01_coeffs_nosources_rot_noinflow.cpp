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

        values[0] = /* Not supported in C: *//* BI *//* BK */0.49379411861125799*std::log(10.0*sqrt(pow(x, 2) + pow(y, 2))) + 2.9533883043631001e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 3.2733586242646999*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 0.478804492016908 ;
    }
};

class Heatflux : public dolfin::Expression {
    public:
    Heatflux() : dolfin::Expression(2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];

        values[0] = /* Not supported in C: *//* BI *//* BI *//* BK *//* BK */-0.185172794479222*x/(pow(x, 2) + pow(y, 2)) - y*(-0.0046391317838279404*BI(1, 3.33333333333333*sqrt(5)*sqrt(pow(x,2) + pow(y,2)))/BI(1, 6.66666666666667*sqrt(5)) + 0.060193212257875599*BK(1, 3.33333333333333*sqrt(5)*sqrt(pow(x,2) + pow(y,2)))/BK(1, 1.66666666666667*sqrt(5)))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[1] = /* Not supported in C: *//* BI *//* BI *//* BK *//* BK */x*(-0.0046391317838279404*BI(1, 3.33333333333333*sqrt(5)*sqrt(pow(x,2) + pow(y,2)))/BI(1, 6.66666666666667*sqrt(5)) + 0.060193212257875599*BK(1, 3.33333333333333*sqrt(5)*sqrt(pow(x,2) + pow(y,2)))/BK(1, 1.66666666666667*sqrt(5)))/sqrt(pow(x, 2) + pow(y, 2)) - 0.185172794479222*y/(pow(x, 2) + pow(y, 2)) ;
    }
};

class Pressure : public dolfin::Expression {
    public:
    Pressure() : dolfin::Expression() {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];

        values[0] = /* Not supported in C: *//* BI *//* BK */7.3834707609077499e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 8.18339656066175*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 0.00503803940928663 ;
    }
};

class Velocity : public dolfin::Expression {
    public:
    Velocity() : dolfin::Expression(2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];

        values[0] = /* Not supported in C: *//* BI *//* BI *//* BK *//* BK */2.1869368024840499e-8*x/(pow(x, 2) + pow(y, 2)) - y*(0.58515350123209398*sqrt(pow(x, 2) + pow(y, 2)) + 0.0018556527135311701*BI(1, 3.33333333333333*sqrt(5)*sqrt(pow(x,2) + pow(y,2)))/BI(1, 6.66666666666667*sqrt(5)) - 0.024077284903150201*BK(1, 3.33333333333333*sqrt(5)*sqrt(pow(x,2) + pow(y,2)))/BK(1, 1.66666666666667*sqrt(5)) - 0.37224945689924099/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[1] = /* Not supported in C: *//* BI *//* BI *//* BK *//* BK */x*(0.58515350123209398*sqrt(pow(x, 2) + pow(y, 2)) + 0.0018556527135311701*BI(1, 3.33333333333333*sqrt(5)*sqrt(pow(x,2) + pow(y,2)))/BI(1, 6.66666666666667*sqrt(5)) - 0.024077284903150201*BK(1, 3.33333333333333*sqrt(5)*sqrt(pow(x,2) + pow(y,2)))/BK(1, 1.66666666666667*sqrt(5)) - 0.37224945689924099/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + 2.1869368024840499e-8*y/(pow(x, 2) + pow(y, 2)) ;
    }
};

class Stress : public dolfin::Expression {
    public:
    Stress() : dolfin::Expression(2,2) {}
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> z) const override {

        double x = z[0];
        double y = z[1];

        values[0] = /* Not supported in C: *//* BI *//* BI *//* BI *//* BI *//* BK *//* BK *//* BK *//* BK */x*(x*(-1.8458676902269401e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 5.5376030706808096e-10*BI(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 6.9689048822689796e-14*BI(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 2.0458491401654402*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 6.1375474204963103*BK(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 4.6692102079126201*BK(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 0.014813819184464099/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + 0.074449891379848102*y/pow(pow(x, 2) + pow(y, 2), 3.0/2.0))/sqrt(pow(x, 2) + pow(y, 2)) - y*(-0.074449891379848102*x/pow(pow(x, 2) + pow(y, 2), 3.0/2.0) - y*(-1.8458676902269401e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 5.5376030706808096e-10*BI(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 2.0458491401654402*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 6.1375474204963103*BK(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 0.014813819184464099/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[1] = /* Not supported in C: *//* BI *//* BI *//* BI *//* BI *//* BK *//* BK *//* BK *//* BK */x*(-0.074449891379848102*x/pow(pow(x, 2) + pow(y, 2), 3.0/2.0) - y*(-1.8458676902269401e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 5.5376030706808096e-10*BI(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 2.0458491401654402*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 6.1375474204963103*BK(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 0.014813819184464099/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + y*(x*(-1.8458676902269401e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 5.5376030706808096e-10*BI(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 6.9689048822689796e-14*BI(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 2.0458491401654402*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 6.1375474204963103*BK(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 4.6692102079126201*BK(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 0.014813819184464099/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + 0.074449891379848102*y/pow(pow(x, 2) + pow(y, 2), 3.0/2.0))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[2] = /* Not supported in C: *//* BI *//* BI *//* BI *//* BI *//* BK *//* BK *//* BK *//* BK */x*(-0.074449891379848102*x/pow(pow(x, 2) + pow(y, 2), 3.0/2.0) - y*(-1.8458676902269401e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 5.5376030706808096e-10*BI(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 2.0458491401654402*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 6.1375474204963103*BK(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 0.014813819184464099/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + y*(x*(-1.8458676902269401e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 5.5376030706808096e-10*BI(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 6.9689048822689796e-14*BI(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 2.0458491401654402*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 6.1375474204963103*BK(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 4.6692102079126201*BK(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 0.014813819184464099/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) + 0.074449891379848102*y/pow(pow(x, 2) + pow(y, 2), 3.0/2.0))/sqrt(pow(x, 2) + pow(y, 2)) ;
        values[3] = /* Not supported in C: *//* BI *//* BI *//* BI *//* BI *//* BK *//* BK *//* BK *//* BK */x*(x*(-1.8458676902269401e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 5.5376030706808096e-10*BI(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 2.0458491401654402*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 6.1375474204963103*BK(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 0.014813819184464099/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) - 0.074449891379848102*y/pow(pow(x, 2) + pow(y, 2), 3.0/2.0))/sqrt(pow(x, 2) + pow(y, 2)) + y*(-0.074449891379848102*x/pow(pow(x, 2) + pow(y, 2), 3.0/2.0) + y*(-1.8458676902269401e-10*BI(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 6.9689048822689796e-14*BI(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 5.5376030706808096e-10*BI(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 6.9689048822689796e-14*BI(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 2.0458491401654402*BK(0, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) - 4.6692102079126201*BK(0, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) + 6.1375474204963103*BK(2, 9.12870929175277*sqrt(pow(x,2) + pow(y,2))) + 4.6692102079126201*BK(2, 12.2474487139159*sqrt(pow(x,2) + pow(y,2))) - 0.014813819184464099/(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)))/sqrt(pow(x, 2) + pow(y, 2)) ;
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
