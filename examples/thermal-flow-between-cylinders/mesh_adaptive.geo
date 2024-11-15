// Command line Parameters
If(!Exists(p))
  p = 5;
EndIf
If(!Exists(d))
  d = 0.5;
EndIf

// Settings
extra_i = 3;
extra_o = 1;
extra_m = 1;
extra_m2 = -1;
extra_d = 0.2;
res = 1.0 * 2^(-(p+extra_i));
res3 = 1.0 * 2^(-(p+extra_o));
res4 = 1.0 * 2^(-(p+extra_m2));
res2 = 1.0 * 2^(-(p+extra_m));
// Mesh.CharacteristicLengthMax = 100;
Mesh.MshFileVersion = 2.0;

// Parameters
R1 = 1.0;
R2 = 2.0;

Point(1) = {0, 0, 0, res};
Point(2) = {0, -d, 0, res};

Point(1000) = { R1, -d, 0, res};
Point(1001) = {-R1, -d, 0, res};

Point(1200) = { (R1+extra_d), -d, 0, res2};
Point(1201) = {-(R1+extra_d), -d, 0, res2};

Point(1300) = { (R2-extra_d), 0, 0, res4};
Point(1301) = {-(R2-extra_d), 0, 0, res4};

Point(1100) = { R2, 0, 0, res3};
Point(1101) = {-R2, 0, 0, res3};

Circle(2000) = {1000,2,1001};
Circle(2001) = {1001,2,1000};

Circle(2200) = {1200,2,1201};
Circle(2201) = {1201,2,1200};

Circle(2300) = {1300,1,1301};
Circle(2301) = {1301,1,1300};

Circle(2100) = {1100,1,1101};
Circle(2101) = {1101,1,1100};

Line Loop(3000) = {2000,2001}; Physical Curve("inner",3000) = {2001, 2000};
Line Loop(3100) = {2100,2101}; Physical Curve("outer",3100) = {2101, 2100};
Line Loop(3200) = {2200,2201};
Line Loop(3300) = {2300,2301};

Plane Surface(4000) = {3100,3300};
Plane Surface(4200) = {3300,3200};
Plane Surface(4100) = {3200,3000};
Physical Surface("mesh",4000) = {4000,4100,4200};
