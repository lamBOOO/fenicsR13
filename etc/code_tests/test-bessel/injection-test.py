import dolfin as df

testInjection = df.Expression("0;\n#include <boost/math/special_functions/bessel.hpp>\n;values[0]=x[0]", degree=2)
df.plot(df.interpolate(testInjection, V))