// Command line Parameters
If(!Exists(p))
  p = 5;
EndIf

// Settings
res1 = 0.005;
res2 = 0.5;
// Mesh.CharacteristicLengthMax = 1.0;
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

Point(1001) = {0, 0, 0, res1};
Point(1002) = {0.5, 0, 0, res1};
Point(1003) = {1.5, 0, 0, res1};
Point(1004) = {2, 0, 0, res1};
Point(1005) = {2, 10, 0, res2};
Point(1006) = {0, 10, 0, res2};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1005};
Line(2005) = {1005,1006};
Line(2006) = {1006,1001};

Line Loop(3000) = {2006}; Physical Curve("left",3000) = {2006};
Line Loop(3001) = {2005}; Physical Curve("upper",3001) = {2005};
Line Loop(3002) = {2004}; Physical Curve("right",3002) = {2004};
Line Loop(3003) = {2001}; Physical Curve("lower_left",3003) = {2001};
Line Loop(3004) = {2002}; Physical Curve("lower_middle",3004) = {2002};
Line Loop(3005) = {2003}; Physical Curve("lower_right",3005) = {2003};

Plane Surface(4000) = {3003,3004,3005,3002,3001,3000}; Physical Surface("mesh",4000) = {4000};
