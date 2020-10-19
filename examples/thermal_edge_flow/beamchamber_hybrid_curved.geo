Mesh.MshFileVersion = 2.0;

// Command line Parameters
If(!Exists(exp))
  exp = 2;
EndIf
If(!Exists(nnodes))
  nnodes = 4 + (7-exp); // for exp=7, n=10 works well
EndIf
Printf("exp=%g", exp);
Printf("nnodes=%g", nnodes);

If (exp<2)
  Printf("exp>=2 needed!");
  Abort;
EndIf

// TODO REMOVE 1/prog!!!! maybe with -LINEID = prog?
// TODO Exlcude middle

dist1 = 0.15;
dist2 = 0.25;

rho1 = 10;
rho2 = 30;

nn1 = 10;
nn2 = 30;
// nn3 = 20;

prog1 = 2.0;
prog2 = 1.1;
prog3 = 1.5;

res1 = 10000;
res2 = 1000;

Point(1001) = {0.0, 0.0, 0, res1};
Point(1002) = {0.0+dist1, 0.0, 0, res1};
Point(1003) = {8.0-dist1, 0.0, 0, res1};
Point(1004) = {8.0, 0.0, 0, res1};
Point(1005) = {8.0, 0.0+dist1, 0, res1};
Point(1006) = {8.0, 8.0-dist1, 0, res1};
Point(1007) = {8.0, 8.0, 0, res1};
Point(1008) = {8.0-dist1, 8.0, 0, res1};
Point(1009) = {0.0+dist1, 8.0, 0, res1};
Point(1010) = {0.0, 8.0, 0, res1};
Point(1011) = {0.0, 8.0-dist1, 0, res1};
Point(1012) = {0.0, 0.0+dist1, 0, res1};

Point(1101) = {0.0+dist1, 0.0+dist1, 0, res1};
Point(1102) = {8.0-dist1, 0.0+dist1, 0, res1};
Point(1103) = {8.0-dist1, 8.0-dist1, 0, res1};
Point(1104) = {0.0+dist1, 8.0-dist1, 0, res1};

// Point(1201) = {1.0+0.0, 1.0+0.0, 0, res2};
Point(1202) = {1.0+0.0+dist2, 1.0+0.0, 0, res2};
Point(1203) = {1.0+2.0-dist2, 1.0+0.0, 0, res2};
// Point(1204) = {1.0+2.0, 1.0+0.0, 0, res2};
Point(1205) = {1.0+2.0, 1.0+0.0+dist2, 0, res2};
Point(1206) = {1.0+2.0, 1.0+2.0-dist2, 0, res2};
// Point(1207) = {1.0+2.0, 1.0+2.0, 0, res2};
Point(1208) = {1.0+2.0-dist2, 1.0+2.0, 0, res2};
Point(1209) = {1.0+0.0+dist2, 1.0+2.0, 0, res2};
// Point(1210) = {1.0+0.0, 1.0+2.0, 0, res2};
Point(1211) = {1.0+0.0, 1.0+2.0-dist2, 0, res2};
Point(1212) = {1.0+0.0, 1.0+0.0+dist2, 0, res2};

Point(1301) = {1.0+0.0+dist2, 1.0+0.0+dist2, 0, res2};
Point(1302) = {1.0+2.0-dist2, 1.0+0.0+dist2, 0, res2};
Point(1303) = {1.0+2.0-dist2, 1.0+2.0-dist2, 0, res2};
Point(1304) = {1.0+0.0+dist2, 1.0+2.0-dist2, 0, res2};


Line(2001) = {1001, 1002}; Transfinite Line {2001} = nn1 Using Progression prog1;
Line(2002) = {1002, 1003}; Transfinite Line {2002} = 8.0*rho1;
Line(2003) = {1003, 1004}; Transfinite Line {-2003} = nn1 Using Progression prog1;
Line(2004) = {1004, 1005}; Transfinite Line {2004} = nn1 Using Progression prog1;
Line(2005) = {1005, 1006}; Transfinite Line {2005} = 8.0*rho1;
Line(2006) = {1006, 1007}; Transfinite Line {-2006} = nn1 Using Progression prog1;
Line(2007) = {1007, 1008}; Transfinite Line {2007} = nn1 Using Progression prog1;
Line(2008) = {1008, 1009}; Transfinite Line {2008} = 8.0*rho1;
Line(2009) = {1009, 1010}; Transfinite Line {-2009} = nn1 Using Progression prog1;
Line(2010) = {1010, 1011}; Transfinite Line {2010} = nn1 Using Progression prog1;
Line(2011) = {1011, 1012}; Transfinite Line {2011} = 8.0*rho1;
Line(2012) = {1012, 1001}; Transfinite Line {-2012} = nn1 Using Progression prog1;


Line(2101) = {1101, 1102}; Transfinite Line {2101} = 8.0*rho1;
Line(2102) = {1102, 1103}; Transfinite Line {2102} = 8.0*rho1;
Line(2103) = {1103, 1104}; Transfinite Line {2103} = 8.0*rho1;
Line(2104) = {1104, 1101}; Transfinite Line {2104} = 8.0*rho1;


// Line(2201) = {1201, 1202}; Transfinite Line {-2201} = nn2 Using Progression prog2;
Line(2202) = {1202, 1203}; Transfinite Line {2202} = 2.0*rho2;
// Line(2203) = {1203, 1204}; Transfinite Line {2203} = nn2 Using Progression prog2;
// Line(2204) = {1204, 1205}; Transfinite Line {-2204} = nn2 Using Progression prog2;
Line(2205) = {1205, 1206}; Transfinite Line {2205} = 2.0*rho2;
// Line(2206) = {1206, 1207}; Transfinite Line {2206} = nn2 Using Progression prog2;
// Line(2207) = {1207, 1208}; Transfinite Line {-2207} = nn2 Using Progression prog2;
Line(2208) = {1208, 1209}; Transfinite Line {2208} = 2.0*rho2;
// Line(2209) = {1209, 1210}; Transfinite Line {2209} = nn2 Using Progression prog2;
// Line(2210) = {1210, 1211}; Transfinite Line {-2210} = nn2 Using Progression prog2;
Line(2211) = {1211, 1212}; Transfinite Line {2211} = 2.0*rho2;
// Line(2212) = {1212, 1201}; Transfinite Line {2212} = nn2 Using Progression prog2;

Line(2301) = {1301, 1302}; Transfinite Line {2301} = 2.0*rho2;
Line(2302) = {1302, 1303}; Transfinite Line {2302} = 2.0*rho2;
Line(2303) = {1303, 1304}; Transfinite Line {2303} = 2.0*rho2;
Line(2304) = {1304, 1301}; Transfinite Line {2304} = 2.0*rho2;

Line(2401) = {1002, 1101}; Transfinite Line {2401} = nn1 Using Progression prog1;
Line(2402) = {1003, 1102}; Transfinite Line {2402} = nn1 Using Progression prog1;
Line(2403) = {1005, 1102}; Transfinite Line {2403} = nn1 Using Progression prog1;
Line(2404) = {1006, 1103}; Transfinite Line {2404} = nn1 Using Progression prog1;
Line(2405) = {1008, 1103}; Transfinite Line {2405} = nn1 Using Progression prog1;
Line(2406) = {1009, 1104}; Transfinite Line {2406} = nn1 Using Progression prog1;
Line(2407) = {1011, 1104}; Transfinite Line {2407} = nn1 Using Progression prog1;
Line(2408) = {1012, 1101}; Transfinite Line {2408} = nn1 Using Progression prog1;

Line(2501) = {1202, 1301}; Transfinite Line {-2501} = nn2 Using Progression prog2;
Line(2502) = {1203, 1302}; Transfinite Line {-2502} = nn2 Using Progression prog2;
Line(2503) = {1205, 1302}; Transfinite Line {-2503} = nn2 Using Progression prog2;
Line(2504) = {1206, 1303}; Transfinite Line {-2504} = nn2 Using Progression prog2;
Line(2505) = {1208, 1303}; Transfinite Line {-2505} = nn2 Using Progression prog2;
Line(2506) = {1209, 1304}; Transfinite Line {-2506} = nn2 Using Progression prog2;
Line(2507) = {1211, 1304}; Transfinite Line {-2507} = nn2 Using Progression prog2;
Line(2508) = {1212, 1301}; Transfinite Line {-2508} = nn2 Using Progression prog2;

Circle(2601) = {1212, 1301, 1202}; Transfinite Line {2601} = rho2 * Pi * 2 * dist2 / 4 + 3;
Circle(2602) = {1203, 1302, 1205}; Transfinite Line {2602} = rho2 * Pi * 2 * dist2 / 4 + 3;
Circle(2603) = {1206, 1303, 1208}; Transfinite Line {2603} = rho2 * Pi * 2 * dist2 / 4 + 3;
Circle(2604) = {1209, 1304, 1211}; Transfinite Line {2604} = rho2 * Pi * 2 * dist2 / 4 + 3;


Curve Loop(3001) = {2001, 2401, -2408, 2012}; Plane Surface(4001) = {3001};
Curve Loop(3002) = {2002, 2402, -2101, -2401}; Plane Surface(4002) = {3002};
Curve Loop(3003) = {2003, 2004, 2403, -2402}; Plane Surface(4003) = {3003};
Curve Loop(3004) = {2403, 2102, -2404, -2005}; Plane Surface(4004) = {3004};
Curve Loop(3005) = {2404, -2405, -2007, -2006}; Plane Surface(4005) = {3005};
Curve Loop(3006) = {2103, -2406, -2008, 2405}; Plane Surface(4006) = {3006};
Curve Loop(3007) = {2407, -2406, 2009, 2010}; Plane Surface(4007) = {3007};
Curve Loop(3008) = {2408, -2104, -2407, 2011}; Plane Surface(4008) = {3008};

Transfinite Surface {4001} = {} AlternateRight;
Transfinite Surface {4002} = {} AlternateRight;
Transfinite Surface {4003} = {} AlternateRight;
Transfinite Surface {4004} = {} AlternateRight;
Transfinite Surface {4005} = {} AlternateRight;
Transfinite Surface {4006} = {} AlternateRight;
Transfinite Surface {4007} = {} AlternateRight;
Transfinite Surface {4008} = {} AlternateRight;

// Curve Loop(3101) = {2201, 2501, -2508, 2212}; Plane Surface(4101) = {3101};
Curve Loop(3102) = {2202, 2502, -2301, -2501}; Plane Surface(4102) = {3102};
// Curve Loop(3103) = {2203, 2204, 2503, -2502}; Plane Surface(4103) = {3103};
Curve Loop(3104) = {2503, 2302, -2504, -2205}; Plane Surface(4104) = {3104};
// Curve Loop(3105) = {2504, -2505, -2207, -2206}; Plane Surface(4105) = {3105};
Curve Loop(3106) = {2303, -2506, -2208, 2505}; Plane Surface(4106) = {3106};
// Curve Loop(3107) = {2507, -2506, 2209, 2210}; Plane Surface(4107) = {3107};
Curve Loop(3108) = {2508, -2304, -2507, 2211}; Plane Surface(4108) = {3108};

// Transfinite Surface {4101} = {} AlternateRight;
Transfinite Surface {4102} = {} AlternateRight;
// Transfinite Surface {4103} = {} AlternateRight;
Transfinite Surface {4104} = {} AlternateRight;
// Transfinite Surface {4105} = {} AlternateRight;
Transfinite Surface {4106} = {} AlternateRight;
// Transfinite Surface {4107} = {} AlternateRight;
Transfinite Surface {4108} = {} AlternateRight;


Curve Loop(3201) = {2501, -2508, 2601};
Curve Loop(3202) = {2503, -2502, 2602};
Curve Loop(3203) = {2505, -2504, 2603};
Curve Loop(3204) = {2507, -2506, 2604};

Plane Surface(4201) = {3201};
Plane Surface(4202) = {3202};
Plane Surface(4203) = {3203};
Plane Surface(4204) = {3204};

Transfinite Surface {4201} = {} Alternate;
Transfinite Surface {4202} = {} Alternate;
Transfinite Surface {4203} = {} Alternate;
Transfinite Surface {4204} = {} AlternateLeft;



Curve Loop(3205) = {2104, 2101, 2102, 2103};
Curve Loop(3206) = {2202, 2602, 2205, 2603, 2208, 2604, 2211, 2601};
Plane Surface(4205) = {3205, 3206};

// Physical Curve("outer",3000) = {
//   210101:210110,
//   220111:221011:100,
//   -211110:-211101,
//   -221001:-220101:100
// };

// Physical Curve("inner",3100) = {
//   210404:210407,
//   220408:220708:100,
//   -210807:-210804,
//   -220704:-220404:100
// };

// Physical Surface("mesh",4000) = {
//   310101:310110,
//   310201:310210,
//   310301:310310,
//   310401:310403,310408:310410,
//   310501:310503,310508:310510,
//   310601:310603,310608:310610,
//   310701:310703,310708:310710,
//   310801:310810,
//   310901:310910,
//   311001:311010
// };
