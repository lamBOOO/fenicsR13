// Command line Parameters
If(!Exists(p))
  p = 0;
EndIf

// Settings
res = 100;
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

// Parameters
R1 = 0.5;
R2 = 2.0;

l = 1.0;
L = 2*l;

Point(1) = {0,  0, 0, res};
Point(2) = {l,  0, 0, res};
Point(3) = {-l, 0, 0, res};

Point(1010) = {l,   R1, 0, res};
Point(1011) = {-l,  R1, 0, res};
Point(1012) = {-l, -R1, 0, res};
Point(1013) = {l,  -R1, 0, res};

Point(1020) = {l,   R2, 0, res};
Point(1021) = {-l,  R2, 0, res};
Point(1022) = {-l, -R2, 0, res};
Point(1023) = {l,  -R2, 0, res};

Line(1010) = {1010,1011};
Circle(1011) = {1011,3,1012};
Line(1012) = {1012,1013};
Circle(1013) = {1013,2,1010};

Line(1020) = {1020,1021};
Circle(1021) = {1021,3,1022};
Line(1022) = {1022,1023};
Circle(1023) = {1023,2,1020};

Line Loop(3010) = {1010}; Physical Curve("it",3010) = {1010};
Line Loop(3011) = {1011}; Physical Curve("il",3011) = {1011};
Line Loop(3012) = {1012}; Physical Curve("ib",3012) = {1012};
Line Loop(3013) = {1013}; Physical Curve("ir",3013) = {1013};

Line Loop(3020) = {1020}; Physical Curve("ot",3020) = {1020};
Line Loop(3021) = {1021}; Physical Curve("ol",3021) = {1021};
Line Loop(3022) = {1022}; Physical Curve("ob",3022) = {1022};
Line Loop(3023) = {1023}; Physical Curve("or",3023) = {1023};

Line Loop(4010) = {1010,1011,1012,1013};
Line Loop(4020) = {1020,1021,1022,1023};
Plane Surface(4000) = {4020,4010}; Physical Surface("mesh",4000) = {4000};
