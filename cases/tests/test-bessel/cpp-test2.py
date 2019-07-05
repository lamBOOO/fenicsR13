from dolfin import *
code = '''
#include <boost/math/special_functions/bessel.hpp>
using boost::math::cyl_bessel_i;
using boost::math::cyl_bessel_j;
using boost::math::cyl_bessel_k;
using boost::math::cyl_neumann;

namespace dolfin {
    class MyFun : public Expression
    {
        public:
            MyFun(): Expression() {};
        void eval(Eigen::Ref<Eigen::VectorXd> values,
                      Eigen::Ref<const Eigen::VectorXd> x) const {
            values[0] = cyl_bessel_j(0,x[0]);
        }
    };
}'''

f=CompiledExpression(compile_cpp_code(code).MyFun(),degree=2)
print(f(1.0))