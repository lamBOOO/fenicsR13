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
Point(1003) = {1, 1, 0, res};
Point(1004) = {0, 1, 0, res};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1001};

Line Loop(3000) = {2003}; Physical Curve("top",3000) = {2003};
Line Loop(3100) = {2002}; Physical Curve("right",3100) = {2002};
Line Loop(3200) = {2001}; Physical Curve("bot",3200) = {2001};
Line Loop(3300) = {2004}; Physical Curve("left",3300) = {2004};

Plane Surface(4000) = {3000,3100,3200,3300}; Physical Surface("mesh",4000) = {4000};
