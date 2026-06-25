// Structured triangular mesh for the thermal-edge-flow simple geometry.
// This variant alternates left/right diagonals to avoid a uniform bias.
//
// Domain: (0, 8)^2 \ [1, 3]^2
// Physical IDs match simple_mesh.geo:
// - 3000: outer wall
// - 3100: inner wall
// - 4000: gas domain

If(!Exists(p))
  p = 5;
EndIf

Mesh.MshFileVersion = 2.0;

n = 2^p;

Point(1001) = {0, 0, 0};
Point(1002) = {1, 0, 0};
Point(1003) = {3, 0, 0};
Point(1004) = {8, 0, 0};
Point(1011) = {0, 1, 0};
Point(1012) = {1, 1, 0};
Point(1013) = {3, 1, 0};
Point(1014) = {8, 1, 0};
Point(1021) = {0, 3, 0};
Point(1022) = {1, 3, 0};
Point(1023) = {3, 3, 0};
Point(1024) = {8, 3, 0};
Point(1031) = {0, 8, 0};
Point(1032) = {1, 8, 0};
Point(1033) = {3, 8, 0};
Point(1034) = {8, 8, 0};

Line(2001) = {1001, 1002};
Line(2002) = {1002, 1003};
Line(2003) = {1003, 1004};
Line(2011) = {1011, 1012};
Line(2012) = {1012, 1013};
Line(2013) = {1013, 1014};
Line(2021) = {1021, 1022};
Line(2022) = {1022, 1023};
Line(2023) = {1023, 1024};
Line(2031) = {1031, 1032};
Line(2032) = {1032, 1033};
Line(2033) = {1033, 1034};

Line(2101) = {1001, 1011};
Line(2102) = {1011, 1021};
Line(2103) = {1021, 1031};
Line(2111) = {1002, 1012};
Line(2112) = {1012, 1022};
Line(2113) = {1022, 1032};
Line(2121) = {1003, 1013};
Line(2122) = {1013, 1023};
Line(2123) = {1023, 1033};
Line(2131) = {1004, 1014};
Line(2132) = {1014, 1024};
Line(2133) = {1024, 1034};

// One interval per physical length unit, refined by 2^p.
Transfinite Curve{2001, 2011, 2021, 2031} = 1 * n + 1;
Transfinite Curve{2002, 2012, 2022, 2032} = 2 * n + 1;
Transfinite Curve{2003, 2013, 2023, 2033} = 5 * n + 1;
Transfinite Curve{2101, 2111, 2121, 2131} = 1 * n + 1;
Transfinite Curve{2102, 2112, 2122, 2132} = 2 * n + 1;
Transfinite Curve{2103, 2113, 2123, 2133} = 5 * n + 1;

Curve Loop(3001) = {2001, 2111, -2011, -2101};
Curve Loop(3002) = {2002, 2121, -2012, -2111};
Curve Loop(3003) = {2003, 2131, -2013, -2121};
Curve Loop(3004) = {2011, 2112, -2021, -2102};
Curve Loop(3005) = {2013, 2132, -2023, -2122};
Curve Loop(3006) = {2021, 2113, -2031, -2103};
Curve Loop(3007) = {2022, 2123, -2032, -2113};
Curve Loop(3008) = {2023, 2133, -2033, -2123};

Plane Surface(4001) = {3001};
Plane Surface(4002) = {3002};
Plane Surface(4003) = {3003};
Plane Surface(4004) = {3004};
Plane Surface(4005) = {3005};
Plane Surface(4006) = {3006};
Plane Surface(4007) = {3007};
Plane Surface(4008) = {3008};

Transfinite Surface{4001} AlternateRight;
Transfinite Surface{4002} AlternateRight;
Transfinite Surface{4003} AlternateRight;
Transfinite Surface{4004} AlternateRight;
Transfinite Surface{4005} AlternateRight;
Transfinite Surface{4006} AlternateRight;
Transfinite Surface{4007} AlternateRight;
Transfinite Surface{4008} AlternateRight;

Physical Curve("outer", 3000) = {
  2001, 2002, 2003,
  2131, 2132, 2133,
  2031, 2032, 2033,
  2101, 2102, 2103
};

Physical Curve("inner", 3100) = {
  2012, 2122, 2022, 2112
};

Physical Surface("mesh", 4000) = {
  4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008
};
