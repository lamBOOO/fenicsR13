// Command line Parameters
If(!Exists(p))
  p = 2;
EndIf

// Settings
res = 100;
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

Point(1001) = {0, 0, 0, res};
Point(1002) = {8, 0, 0, res};
Point(1003) = {8, 8, 0, res};
Point(1004) = {0, 8, 0, res};

Point(1011) = {1, 1, 0, res};
Point(1012) = {3, 1, 0, res};
Point(1013) = {3, 3, 0, res};
Point(1014) = {1, 3, 0, res};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1001};

Line(2011) = {1011,1012};
Line(2012) = {1012,1013};
Line(2013) = {1013,1014};
Line(2014) = {1014,1011};

Line Loop(3000) = {2001,2002,2003,2004}; Physical Curve("outer",3000) = {2001,2002,2003,2004};
Line Loop(3100) = {2011,2012,2013,2014}; Physical Curve("inner",3100) = {2011,2012,2013,2014};

Plane Surface(4000) = {3000,3100}; Physical Surface("mesh",4000) = {4000};
