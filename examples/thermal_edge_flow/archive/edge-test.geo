

res = 1000;
nnodes1 = 10;
nnodes2 = 11;
prog = 1.2;

Point(1001) = {0.0, 0.0, 0, res};
Point(1002) = {1.0, 0.0, 0, res};
Point(1003) = {0.0, 1.0, 0, res};

Point(1004) = {0.0, -1.0, 0, res};
Point(1005) = {1.0, -1.0, 0, res};

Line(2001) = {1001,1002}; Transfinite Line{2001} = nnodes1 Using Progression prog;
Line(2002) = {1001,1003}; Transfinite Line{2002} = nnodes1 Using Progression prog;
Circle(2003) = {1002,1001,1003}; Transfinite Line{2003} = nnodes2;

Curve Loop(1) = {2001, 2003, -2002};
Plane Surface(1) = {1}; Transfinite Surface{1} = {} Alternate;
//+
Line(2004) = {1001, 1004}; Transfinite Line{2004} = nnodes2;
//+
Line(2005) = {1004, 1005}; Transfinite Line{2005} = nnodes1 Using Progression prog;
//+
Line(2006) = {1005, 1002}; Transfinite Line{2006} = nnodes2;
//+
Curve Loop(2) = {2001, -2006, -2005, -2004};
//+
Plane Surface(2) = {2}; Transfinite Surface{2} = {} Alternate;
