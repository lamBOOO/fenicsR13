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

        const double x2 = x * x,   y2 = y * y;
        const double x3 = x2 * x,  y3 = y2 * y;
        const double x4 = x3 * x,  y4 = y3 * y;
        const double x5 = x4 * x,  y5 = y4 * y;
        const double x6 = x5 * x,  y6 = y5 * y;
        const double x7 = x6 * x,  y7 = y6 * y;
        const double x8 = x7 * x,  y8 = y7 * y;
        const double x9 = x8 * x,  y9 = y8 * y;
        const double x10 = x9 * x, y10 = y9 * y;
        const double x11 = x10 * x,y11 = y10 * y;
        const double x12 = x11 * x,y12 = y11 * y;
        const double x13 = x12 * x,y13 = y12 * y;
        const double x14 = x13 * x,y14 = y13 * y;
        const double x15 = x14 * x;

        const double r2  = x2 + y2;
        const double r   = std::sqrt(r2);
        const double lr  = std::log(r2);
        const double arg = 18.2501 * r;

        const double bi0 = cyl_bessel_i(0, arg);
        const double bi1 = cyl_bessel_i(1, arg);
        const double bk0 = cyl_bessel_k(0, arg);
        const double bk1 = cyl_bessel_k(1, arg);

        // Plain polynomial terms (no sqrt, no log, no bessel)
        double poly =
            1.86265e-20 * x5    + 7.10234e-20 * x6    + 4.13774e-17 * x7
            + 9.34785e-16 * x8    - 0.00680351  * x9    + 1.14756      * x10
            + 0.00958524  * x11   + 3.53952e-8  * x12   - 3.61242e-7   * x13
            + 8.29357e-7  * x14   - 5.02999e-6  * x15
            + 3.72529e-20 * x3*y2 + 2.1307e-19  * x4*y2 + 1.24132e-16  * x5*y2
            + 3.73914e-15 * x6*y2 - 0.0272141   * x7*y2 + 5.73781       * x8*y2
            + 0.0479262   * x9*y2 + 2.12371e-7  * x10*y2- 2.16745e-6   * x11*y2
            + 5.8055e-6   * x12*y2- 3.52099e-5  * x13*y2
            + 1.86265e-20 * x*y4  + 2.1307e-19  * x2*y4 + 1.24132e-16  * x3*y4
            + 5.60871e-15 * x4*y4 - 0.0408211   * x5*y4 + 11.4756       * x6*y4
            + 0.0958524   * x7*y4 + 5.30928e-7  * x8*y4 - 5.41863e-6   * x9*y4
            + 1.74165e-5  * x10*y4- 1.0563e-4   * x11*y4
            + 7.10234e-20 * y6    + 4.13774e-17 * x*y6  + 3.73914e-15  * x2*y6
            - 0.0272141   * x3*y6 + 11.4756     * x4*y6 + 0.0958524     * x5*y6
            + 7.07904e-7  * x6*y6 - 7.22484e-6  * x7*y6 + 2.90275e-5   * x8*y6
            - 1.7605e-4   * x9*y6
            + 9.34785e-16 * y8    - 0.00680351  * x*y8  + 5.73781       * x2*y8
            + 0.0479262   * x3*y8 + 5.30928e-7  * x4*y8 - 5.41863e-6   * x5*y8
            + 2.90275e-5  * x6*y8 - 1.7605e-4   * x7*y8
            + 1.14756     * y10   + 0.00958524  * x*y10 + 2.12371e-7    * x2*y10
            - 2.16745e-6  * x3*y10+ 1.74165e-5  * x4*y10- 1.0563e-4    * x5*y10
            + 3.53952e-8  * y12   - 3.61242e-7  * x*y12 + 5.8055e-6     * x2*y12
            - 3.52099e-5  * x3*y12
            + 8.29357e-7  * y14   - 5.02999e-6  * x*y14;

        // Terms with sqrt(x^2+y^2)
        double sqrtTerms = r * (
            - 5.44199e-20 * x4    + 1.43567e-18 * x5    - 6.71292e-17 * x6
            + 4.61042e-16 * x7    + 2.36759e-13 * x8    + 4.27934e-13 * x9
            + 1.54024e-11 * x10   - 4.5631e-10  * x11   - 3.11377e-9  * x12
            - 9.10296e-8  * x13   + 4.24103e-7  * x14
            - 1.0884e-19  * x2*y2 + 2.87133e-18 * x3*y2 - 2.01388e-16 * x4*y2
            + 1.38313e-15 * x5*y2 + 9.47034e-13 * x6*y2 + 1.71174e-12 * x7*y2
            + 7.70122e-11 * x8*y2 - 2.28155e-9  * x9*y2 - 1.86826e-8  * x10*y2
            - 5.46177e-7  * x11*y2+ 2.96872e-6  * x12*y2
            - 5.44199e-20 * y4    + 1.43567e-18 * x*y4  - 2.01388e-16 * x2*y4
            + 1.38313e-15 * x3*y4 + 1.42055e-12 * x4*y4 + 2.5676e-12  * x5*y4
            + 1.54024e-10 * x6*y4 - 4.5631e-9   * x7*y4 - 4.67066e-8  * x8*y4
            - 1.36544e-6  * x9*y4 + 8.90617e-6  * x10*y4
            - 6.71292e-17 * y6    + 4.61042e-16 * x*y6  + 9.47034e-13 * x2*y6
            + 1.71174e-12 * x3*y6 + 1.54024e-10 * x4*y6 - 4.5631e-9   * x5*y6
            - 6.22754e-8  * x6*y6 - 1.82059e-6  * x7*y6 + 1.48436e-5  * x8*y6
            + 2.36759e-13 * y8    + 4.27934e-13 * x*y8  + 7.70122e-11 * x2*y8
            - 2.28155e-9  * x3*y8 - 4.67066e-8  * x4*y8 - 1.36544e-6  * x5*y8
            + 1.48436e-5  * x6*y8
            + 1.54024e-11 * y10   - 4.5631e-10  * x*y10 - 1.86826e-8  * x2*y10
            - 5.46177e-7  * x3*y10+ 8.90617e-6  * x4*y10
            - 3.11377e-9  * y12   - 9.10296e-8  * x*y12 + 2.96872e-6  * x2*y12
            + 4.24103e-7  * y14
        );

        // Log terms (no sqrt)
        double logTerms = lr * (
            0.565296  * x10   + 2.82648  * x8*y2  + 5.65296  * x6*y4
            + 5.65296   * x4*y6 + 2.82648  * x2*y8  + 0.565296 * y10
        );

        // Log * sqrt terms
        double logSqrtTerms = lr * r * (
            - 3.40893e-13 * x9    - 1.36357e-12 * x7*y2 - 2.04536e-12 * x5*y4
            - 1.36357e-12 * x3*y6 - 3.40893e-13 * x*y8
        );

        // Bessel I0 terms
        double besselI0Terms = bi0 * (
            2.11916e-17 * x10   + 1.05958e-16 * x8*y2  + 2.11916e-16 * x6*y4
            + 2.11916e-16 * x4*y6 + 1.05958e-16 * x2*y8  + 2.11916e-17 * y10
        );

        // Bessel I1 terms (has extra x*r factor)
        double besselI1Terms = bi1 * x * r * (
            - 2.38918e-17 * x8    - 9.55671e-17 * x6*y2  - 1.43351e-16 * x4*y4
            - 9.55671e-17 * x2*y6 - 2.38918e-17 * y8
        );

        // Bessel K0 terms
        double besselK0Terms = bk0 * (
            - 7.51623e6  * x10   - 3.75811e7   * x8*y2  - 7.51623e7   * x6*y4
            - 7.51623e7  * x4*y6 - 3.75811e7   * x2*y8  - 7.51623e6   * y10
        );

        // Bessel K1 terms (has extra r factor)
        double besselK1Terms = bk1 * r * (
            - 1.72973e6  * x9    - 6.9189e6    * x7*y2  - 1.03784e7   * x5*y4
            - 6.9189e6   * x3*y6 - 1.72973e6   * x*y8
        );

        double numerator = poly + sqrtTerms + logTerms + logSqrtTerms
                        + besselI0Terms + besselI1Terms
                        + besselK0Terms + besselK1Terms;

        values[0] = 1.1475416669755338 - 1.4136499460087132e-19/Power(r,5) - \
                    1.9855291553960966e-18/Power(r,4) - 8.503184547098384e-18/Power(r,3) - \
                    7.51454581320324e-16/Power(r,2) - 0.006824199418861071/r + \
                    0.009607734992218002*r + 1.1347378340243468e-9*Power(r,2) - \
                    2.861003252615381e-7*Power(r,3) + 7.835105240179381e-8*Power(r,4) - \
                    3.2050302992908996e-6*Power(r,5) + \
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
