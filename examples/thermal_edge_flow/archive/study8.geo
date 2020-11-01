Mesh.MshFileVersion = 2.0;

// Command line Parameters
If(!Exists(exp1))
  exp1 = 6;
EndIf
Printf("exp1=%g", exp1);
If(!Exists(exp2))
  exp2 = 6;
EndIf
Printf("exp2=%g", exp2);
If(!Exists(exp3))
  exp3 = 13;
EndIf
Printf("exp3=%g", exp3);
If(!Exists(split))
  split = 0;
EndIf
Printf("split=%g", split);

res1 = 2^-exp1;
res2 = 2^-exp2;
res3 = 2^-exp3;

dist1 = 0.10 * 2^-(exp1-4);
dist2 = 0.30;
dist3 = 0.45;
rho1 = 1 / res1;
nn1 = 10;
prog1 = 1.5;


Point(1001) = {0.0, 0.0, 0, 10000};
Point(1002) = {0.0+dist1, 0.0, 0, 10000};
Point(1003) = {8.0-dist1, 0.0, 0, 10000};
Point(1004) = {8.0, 0.0, 0, 10000};
Point(1005) = {8.0, 0.0+dist1, 0, 10000};
Point(1006) = {8.0, 8.0-dist1, 0, 10000};
Point(1007) = {8.0, 8.0, 0, 10000};
Point(1008) = {8.0-dist1, 8.0, 0, 10000};
Point(1009) = {0.0+dist1, 8.0, 0, 10000};
Point(1010) = {0.0, 8.0, 0, 10000};
Point(1011) = {0.0, 8.0-dist1, 0, 10000};
Point(1012) = {0.0, 0.0+dist1, 0, 10000};

Point(1101) = {0.0+dist1, 0.0+dist1, 0, 10000};
Point(1102) = {8.0-dist1, 0.0+dist1, 0, 10000};
Point(1103) = {8.0-dist1, 8.0-dist1, 0, 10000};
Point(1104) = {0.0+dist1, 8.0-dist1, 0, 10000};

Point(1202) = {1.0, 1.0-dist2, 0, res2};
Point(1203) = {3.0, 1.0-dist2, 0, res2};
Point(1205) = {3.0+dist2, 1.0, 0, res2};
Point(1206) = {3.0+dist2, 3.0, 0, res2};
Point(1208) = {3.0, 3.0+dist2, 0, res2};
Point(1209) = {1.0, 3.0+dist2, 0, res2};
Point(1211) = {1.0-dist2, 3.0, 0, res2};
Point(1212) = {1.0-dist2, 1.0, 0, res2};

Point(1301) = {1.0, 1.0, 0, res3};
Point(1302) = {3.0, 1.0, 0, res3};
Point(1303) = {3.0, 3.0, 0, res3};
Point(1304) = {1.0, 3.0, 0, res3};

Point(1402) = {1.0, 1.0-(dist2+dist3), 0, res1};
Point(1403) = {3.0, 1.0-(dist2+dist3), 0, res1};
Point(1405) = {3.0+(dist2+dist3), 1.0, 0, res1};
Point(1406) = {3.0+(dist2+dist3), 3.0, 0, res1};
Point(1408) = {3.0, 3.0+(dist2+dist3), 0, res1};
Point(1409) = {1.0, 3.0+(dist2+dist3), 0, res1};
Point(1411) = {1.0-(dist2+dist3), 3.0, 0, res1};
Point(1412) = {1.0-(dist2+dist3), 1.0, 0, res1};


Line(2001) = {1001, 1002}; Transfinite Line {2001} = nn1 Using Progression prog1;
Line(2002) = {1002, 1003}; Transfinite Line {2002} = 8.0*rho1*2^+0;
Line(2003) = {1003, 1004}; Transfinite Line {-2003} = nn1 Using Progression prog1;
Line(2004) = {1004, 1005}; Transfinite Line {2004} = nn1 Using Progression prog1;
Line(2005) = {1005, 1006}; Transfinite Line {2005} = 8.0*rho1*2^-1;
Line(2006) = {1006, 1007}; Transfinite Line {-2006} = nn1 Using Progression prog1;
Line(2007) = {1007, 1008}; Transfinite Line {2007} = nn1 Using Progression prog1;
Line(2008) = {1008, 1009}; Transfinite Line {2008} = 8.0*rho1*2^-1;
Line(2009) = {1009, 1010}; Transfinite Line {-2009} = nn1 Using Progression prog1;
Line(2010) = {1010, 1011}; Transfinite Line {2010} = nn1 Using Progression prog1;
Line(2011) = {1011, 1012}; Transfinite Line {2011} = 8.0*rho1*2^+0;
Line(2012) = {1012, 1001}; Transfinite Line {-2012} = nn1 Using Progression prog1;


Line(2101) = {1101, 1102}; Transfinite Line {2101} = 8.0*rho1*2^+0;
Line(2102) = {1102, 1103}; Transfinite Line {2102} = 8.0*rho1*2^-1;
Line(2103) = {1103, 1104}; Transfinite Line {2103} = 8.0*rho1*2^-1;
Line(2104) = {1104, 1101}; Transfinite Line {2104} = 8.0*rho1*2^+0;


Line(2202) = {1202, 1203};
Line(2205) = {1205, 1206};
Line(2208) = {1208, 1209};
Line(2211) = {1211, 1212};

Line(2301) = {1301, 1302};
Line(2302) = {1302, 1303};
Line(2303) = {1303, 1304};
Line(2304) = {1304, 1301};

Line(2401) = {1002, 1101}; Transfinite Line {2401} = nn1 Using Progression prog1;
Line(2402) = {1003, 1102}; Transfinite Line {2402} = nn1 Using Progression prog1;
Line(2403) = {1005, 1102}; Transfinite Line {2403} = nn1 Using Progression prog1;
Line(2404) = {1006, 1103}; Transfinite Line {2404} = nn1 Using Progression prog1;
Line(2405) = {1008, 1103}; Transfinite Line {2405} = nn1 Using Progression prog1;
Line(2406) = {1009, 1104}; Transfinite Line {2406} = nn1 Using Progression prog1;
Line(2407) = {1011, 1104}; Transfinite Line {2407} = nn1 Using Progression prog1;
Line(2408) = {1012, 1101}; Transfinite Line {2408} = nn1 Using Progression prog1;

Line(2501) = {1202, 1301};
Line(2502) = {1203, 1302};
Line(2503) = {1205, 1302};
Line(2504) = {1206, 1303};
Line(2505) = {1208, 1303};
Line(2506) = {1209, 1304};
Line(2507) = {1211, 1304};
Line(2508) = {1212, 1301};

Circle(2601) = {1212, 1301, 1202};
Circle(2602) = {1203, 1302, 1205};
Circle(2603) = {1206, 1303, 1208};
Circle(2604) = {1209, 1304, 1211};

Circle(2701) = {1412, 1301, 1402};
Circle(2702) = {1403, 1302, 1405};
Circle(2703) = {1406, 1303, 1408};
Circle(2704) = {1409, 1304, 1411};

Line(2802) = {1402, 1403};
Line(2805) = {1405, 1406};
Line(2808) = {1408, 1409};
Line(2811) = {1411, 1412};


Line(2901) = {1402, 1202};
Line(2902) = {1403, 1203};
Line(2903) = {1405, 1205};
Line(2904) = {1406, 1206};
Line(2905) = {1408, 1208};
Line(2906) = {1409, 1209};
Line(2907) = {1411, 1211};
Line(2908) = {1412, 1212};


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

Curve Loop(3102) = {2202, 2502, -2301, -2501}; Plane Surface(4102) = {3102};
Curve Loop(3104) = {2503, 2302, -2504, -2205}; Plane Surface(4104) = {3104};
Curve Loop(3106) = {2303, -2506, -2208, 2505}; Plane Surface(4106) = {3106};
Curve Loop(3108) = {2508, -2304, -2507, 2211}; Plane Surface(4108) = {3108};

Curve Loop(3201) = {2501, -2508, 2601};
Curve Loop(3202) = {2503, -2502, 2602};
Curve Loop(3203) = {2505, -2504, 2603};
Curve Loop(3204) = {2507, -2506, 2604};

Plane Surface(4201) = {3201};
Plane Surface(4202) = {3202};
Plane Surface(4203) = {3203};
Plane Surface(4204) = {3204};

Curve Loop(3307) = {2802, 2902, -2202, -2901}; Plane Surface(4306) = {3307};
Curve Loop(3308) = {2903, 2205, -2904, -2805}; Plane Surface(4307) = {3308};
Curve Loop(3309) = {2208, -2906, -2808, 2905}; Plane Surface(4308) = {3309};
Curve Loop(3310) = {2908, -2211, -2907, 2811}; Plane Surface(4309) = {3310};

Curve Loop(3411) = {2901, -2601, -2908, 2701}; Plane Surface(4410) = {3411};
Curve Loop(3412) = {2903, -2602, -2902, 2702}; Plane Surface(4411) = {3412};
Curve Loop(3413) = {2905, -2603, -2904, 2703}; Plane Surface(4412) = {3413};
Curve Loop(3414) = {2907, -2604, -2906, 2704}; Plane Surface(4413) = {3414};


Curve Loop(3205) = {2104, 2101, 2102, 2103};
Curve Loop(3206) = {2802, 2702, 2805, 2703, 2808, 2704, 2811, 2701};
Plane Surface(4205) = {3205, 3206};

Physical Curve("outer",3000) = {
  2001:2012
};

Physical Curve("inner",3100) = {
  2301:2304
};

Physical Surface("mesh",4000) = {
  4001:4008,
  4102, 4104, 4106, 4108,
  4201:4204,
  4205,
  4306:4309,
  4410:4413
};

If (split != 0)
  Mesh 2;
EndIf
For i In{1:split}
  Printf("Split=%g", i);
  RefineMesh;
EndFor
