If(!Exists(p))
  p = 0;
EndIf

Mesh.MshFileVersion = 2.0;

// Define the corners of the box
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0, 0, 1, 1.0};
Point(6) = {1, 0, 1, 1.0};
Point(7) = {1, 1, 1, 1.0};
Point(8) = {0, 1, 1, 1.0};

// Define the lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Define the surfaces
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};Physical Surface(4000) = 14;
Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};Physical Surface(4001) = 16;
Line Loop(17) = {1, 10, -5, -9};
Plane Surface(18) = {17};Physical Surface(4002) = 18;
Line Loop(19) = {2, 11, -6, -10};
Plane Surface(20) = {19};Physical Surface(4003) = 20;
Line Loop(21) = {3, 12, -7, -11};
Plane Surface(22) = {21};Physical Surface(4004) = 22;
Line Loop(23) = {4, 9, -8, -12};
Plane Surface(24) = {23}; Physical Surface(4005) = 24;

// Define the volume
Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25}; Physical Volume(6000) = 26;

// Mesh refinements
Mesh.CharacteristicLengthMax = 1/(2^p);
