// Command line Parameters
If(!Exists(p))
  p = 0;
EndIf

// Settings
res = 100;
Mesh.MshFileVersion = 2.0;

// Parameters
R1 = 0.5;
R2 = 2.0;

// Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.CharacteristicLengthMax = 2 * Pi * R1 * ( 1 / (2^(p+1)+2) );

Point(1) = {0, 0, 0, res};

Point(1000) = { R1, 0, 0, res};
Point(1001) = {-R1, 0, 0, res};

Point(1100) = { R2, 0, 0, res};
Point(1101) = {-R2, 0, 0, res};

Circle(2000) = {1000,1,1001};
Circle(2001) = {1001,1,1000};

Circle(2100) = {1100,1,1101};
Circle(2101) = {1101,1,1100};

Line Loop(3000) = {2000,2001}; Physical Curve("inner",3000) = {2001, 2000};
Line Loop(3100) = {2100,2101}; Physical Curve("outer",3100) = {2101, 2100};

Plane Surface(4000) = {3100,3000}; Physical Surface("mesh",4000) = {4000};
