//+
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(2) = {-0, 0, 0, 2, -Pi/2, Pi/2, 2*Pi};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
//+
Physical Surface("inner1", 10000) = {1};
//+
Physical Surface("outer", 10001) = {2};
