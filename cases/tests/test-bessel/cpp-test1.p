from dolfin import *

cppcode = """
class K : public Expression
{
public:
void eval(Array<double>& values,
const Array<double>& x,
const ufc::cell& cell) const
{
if ((*materials)[cell.index] == 0)
values[0] = k_0;
else
values[0] = k_1;
}
std::shared_ptr<MeshFunction<std::size_t>> materials;
double k_0;
double k_1;
};
"""
kappa = Expression(cppcode=cppcode, degree=0)