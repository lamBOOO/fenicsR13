// Command line Parameters
If(!Exists(p1))
  p1 = -5;
EndIf
If(!Exists(p2))
  p2 = -6;
EndIf
If(!Exists(p3))
  p3 = -11;
EndIf
Printf("p1=%g", p1);
Printf("p2=%g", p2);
Printf("p3=%g", p3);

// Settings
res1 = 2^p1;
res2 = 2^p2;
res3 = 2^p3;

// Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

Point(1001) = {0, 0, 0, res3};
Point(1002) = {4, 0, 0, res3};
Point(1003) = {8, 0, 0, res2};
Point(1004) = {8, 8, 0, res1};
Point(1005) = {0, 8, 0, res2};
Point(1006) = {0, 4, 0, res3};

Point(1011) = {1, 1, 0, res3};
Point(1012) = {3, 1, 0, res3};
Point(1013) = {3, 3, 0, res3};
Point(1014) = {1, 3, 0, res3};

Point(1021) = {4, 4, 0, res2};

Point(1031) = {0.5, 0.5, 0, res2};
Point(1032) = {3.5, 0.5, 0, res2};
Point(1033) = {3.5, 3.5, 0, res2};
Point(1034) = {0.5, 3.5, 0, res2};

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

Line(2031) = {1031,1032};
Line(2032) = {1032,1033};
Line(2033) = {1033,1034};
Line(2034) = {1034,1031};

Line Loop(3000) = {2001:2006}; Physical Curve("outer",3000) = {2001:2006};
Line Loop(3100) = {2011,2012,2013,2014}; Physical Curve("inner",3100) = {2011,2012,2013,2014};

Line Loop(3200) = {2001,2021,2022,2006};
Line Loop(3300) = {2002:2005,-2022,-2021};
Line Loop(3400) = {2031,2032,2033,2034};

Plane Surface(4001) = {3200,3400};
Plane Surface(4002) = {3300};
Plane Surface(4003) = {3400,3100};
Physical Surface("mesh",4000) = {4001,4002,4003};
