// Command line Parameters
If(!Exists(p))
  p = 5;
EndIf

// Settings
res = 100;
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

Point(1001) = {0, 0, 0, res};
Point(1002) = {1, 0, 0, res};
Point(1003) = {1, 0.5, 0, res};
Point(1004) = {1, 1, 0, res};
Point(1005) = {0, 1, 0, res};
Point(1006) = {0, 0.5, 0, res};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1005};
Line(2005) = {1005,1006};
Line(2006) = {1006,1001};

Line(2007) = {1003,1006};

Line Loop(3000) = {2004}; Physical Curve("upper",3000) = {2004};
Line Loop(3100) = {2005,2006,2001,2002,2003}; Physical Curve("lower",3100) = {2005,2006,2001,2002,2003};
// Line Loop(3900) = {2007}; Physical Curve("mid",3900) = {2007};

Line Loop(3200) = {2005,-2007,2003,2004};
Line Loop(3300) = {2006,2001,2002,2007};
Plane Surface(4000) = {3200}; Physical Surface("top",4000) = {4000};
Plane Surface(4100) = {3300}; Physical Surface("bot",4100) = {4100};

