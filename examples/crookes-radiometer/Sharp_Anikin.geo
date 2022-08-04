// Gmsh project created on Thu Jan 27 11:23:26 2022
SetFactory("OpenCASCADE");
//+
L = 1;
//+
Circle(1) = {0, 0, 0, 5.66*L, 0, 2*Pi};
//+
Point(2) = {0.1*L, 1.1*L, 0, 0.01};
//+
Point(3) = {-0.1*L, 1.1*L, 0, 0.01};
//+
Point(4) = {-0.1*L, -1.1*L, 0, 0.01};
//+
Point(5) = {0.1*L, -1.1*L, 0, 0.01};
//+
Point(6) = {1.1*L, 0.1*L, 0, 0.01};
//+
Point(7) = {-1.1*L, 0.1*L, 0, 0.01};
//+
Point(8) = {-1.1*L, -0.1*L, 0, 0.01};
//+
Point(9) = {1.1*L, -0.1*L, 0, 0.01};
//+
Point(10) = {0.1*L, 0.1*L, 0, 0.01};
//+
Point(11) = {-0.1*L, 0.1*L, 0, 0.01};
//+
Point(12) = {-0.1*L, -0.1*L, 0, 0.01};
//+
Point(13) = {0.1*L, -0.1*L, 0, 0.01};
//+
Point(14) = {0.1*L, 0.1*L, 0, 0.01};
//+
Point(15) = {-0.1*L, 0.1*L, 0, 0.01};
//+
Point(16) = {-0.1*L, -0.1*L, 0, 0.01};
//+
Point(17) = {0.1*L, -0.1*L, 0, 0.01};
//+
Line(2) = {7, 11};
//+
Line(3) = {11, 3};
//+
Line(4) = {3, 2};
//+
Line(5) = {2, 10};
//+
Line(6) = {10, 6};
//+
Line(7) = {6, 9};
//+
Line(8) = {9, 13};
//+
Line(9) = {13, 5};
//+
Line(10) = {5, 4};
//+
Line(11) = {4, 12};
//+
Line(12) = {12, 8};
//+
Line(13) = {8, 7};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("Cold_Side", 1000) = {5, 8, 11, 2};
//+
Physical Curve("Hot_Side", 10001) = {6, 9, 12, 3};
//+
Physical Curve("Outer", 10002) = {1};
//+
Physical Curve("Trans X-Top", 10003) = {4};
//+
Physical Curve("Trans X-Bottom", 10004) = {10};
//+
Physical Curve("Trans Y-Left", 10005) = {13};
//+
Physical Curve("Trans Y-Right", 10006) = {7};
//+
Physical Surface("Inner", 10007) = {1};
