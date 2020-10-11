// Command line Parameters
If(!Exists(p))
  p = 2;
EndIf

// Settings
res1 = 0.0005;
res2 = 0.04;
res3 = 0.1;
// Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

Point(1001) = {0, 0, 0, res2};
Point(1002) = {4, 0, 0, res2};
Point(1003) = {8, 0, 0, res3};
Point(1004) = {8, 8, 0, res3};
Point(1005) = {0, 8, 0, res3};
Point(1006) = {0, 4, 0, res2};

Point(1011) = {1, 1, 0, res1};
Point(1012) = {3, 1, 0, res1};
Point(1013) = {3, 3, 0, res1};
Point(1014) = {1, 3, 0, res1};

Point(1021) = {4, 4, 0, res2};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1005};
Line(2005) = {1005,1006};
Line(2006) = {1006,1001};

Line(2011) = {1011,1012};
Line(2012) = {1012,1013};
Line(2013) = {1013,1014};
Line(2014) = {1014,1011};

Line(2021) = {1002,1021};
Line(2022) = {1021,1006};

Line Loop(3000) = {2001:2006}; Physical Curve("outer",3000) = {2001:2006};
Line Loop(3100) = {2011,2012,2013,2014}; Physical Curve("inner",3100) = {2011,2012,2013,2014};

Line Loop(3200) = {2001,2021,2022,2006};
Line Loop(3300) = {2002:2005,-2022,-2021};

Plane Surface(4001) = {3200,3100};
Plane Surface(4002) = {3300};
Physical Surface("mesh",4000) = {4001,4002};
