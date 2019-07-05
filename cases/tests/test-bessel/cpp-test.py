from dolfin import *
import numpy as np

mesh = UnitSquareMesh(5,5)
V = FunctionSpace(mesh,'CG',1)
V2 = VectorFunctionSpace(mesh,'CG',1)
g = interpolate(Expression(('x[0]+x[1]','-x[0]')),V2)
f = interpolate(Expression('1.+x[0]*x[0]'),V)

code = '''
class MyFunc : public Expression
{
public:
  boost::shared_ptr<Function> f;
  boost::shared_ptr<Function> g;

  MyFunc() : Expression() { }

  void eval(Array<double>& values, const Array<double>& x,
            const ufc::cell& c) const
  {
      Array<double> val(2);
      g->eval(val, x, c);
      f->eval(values, val, c);
  }
};
'''

fg = Expression(code)
fg.f = f
fg.g = g

# Find cell that intersects Point(0.1, 0.2)
p = Point(0.1, 0.2)
bt = mesh.bounding_box_tree()
cell_id = bt.compute_first_entity_collision(p)
cell = Cell(mesh, cell_id)
values = np.zeros(1)
x = np.array([p.x(), p.y()])
fg.eval_cell(values, x, cell)
print(values)