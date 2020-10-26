Mesh.MshFileVersion = 2.0;

// Command line Parameters

If(!Exists(split))
  split = 0;
EndIf
Printf("split=%g", split);


dist1 = 0.12;
dist2 = 0.4;

res_tmp1 = 2^-split * 2^-6; // leftbot edge
res_tmp2 = 2^-split * 2^-3; // bulk
res_tmp3 = 2^-split * 2^-4; // topright edge
res_tmp4 = 2^-split * 2^-7; // leftbot focus
res_tmp5 = 2^-split * 2^-12; // inner edge focus
res_tmp6 = 2^-split * 2^-7; // inner edge
res_tmp7 = 2^-split * 2^-4; // focus bulk


Point(1001) = {0.0, 0.0, 0, res_tmp1};
Point(1002) = {2.0, 0.0, 0, res_tmp4};
Point(1003) = {4.0, 0.0, 0, res_tmp1};
Point(1004) = {8.0-dist1, 0.0, 0, res_tmp3};
Point(1005) = {8.0, 0.0, 0, res_tmp3};
Point(1006) = {8.0, 0.0+dist1, 0, res_tmp3};
Point(1007) = {8.0, 8.0-dist1, 0, res_tmp3};
Point(1008) = {8.0, 8.0, 0, res_tmp3};
Point(1009) = {8.0-dist1, 8.0, 0, res_tmp3};
Point(1010) = {0.0+dist1, 8.0, 0, res_tmp3};
Point(1011) = {0.0, 8.0, 0, res_tmp3};
Point(1012) = {0.0, 8.0-dist1, 0, res_tmp3};
Point(1013) = {0.0, 4.0, 0, res_tmp1};
Point(1014) = {0.0, 2.0, 0, res_tmp4};

Point(1101) = {0.0+dist1, 0.0+dist1, 0, res_tmp7};
Point(1102) = {4.0, 0.0+dist1, 0, res_tmp7};
Point(1103) = {8.0-dist1, 0.0+dist1, 0, res_tmp2};
Point(1104) = {8.0-dist1, 8.0-dist1, 0, res_tmp2};
Point(1105) = {0.0+dist1, 8.0-dist1, 0, res_tmp2};
Point(1106) = {0.0+dist1, 4.0, 0, res_tmp7};

Point(1202) = {1.0, 1.0-dist2, 0, res_tmp7};
Point(1203) = {3.0, 1.0-dist2, 0, res_tmp7};
Point(1205) = {3.0+dist2, 1.0, 0, res_tmp7};
Point(1206) = {3.0+dist2, 3.0, 0, res_tmp7};
Point(1208) = {3.0, 3.0+dist2, 0, res_tmp7};
Point(1209) = {1.0, 3.0+dist2, 0, res_tmp7};
Point(1211) = {1.0-dist2, 3.0, 0, res_tmp7};
Point(1212) = {1.0-dist2, 1.0, 0, res_tmp7};

Point(1301) = {1.0, 1.0, 0, res_tmp5};
Point(1302) = {1.0+dist2, 1.0, 0, res_tmp6};
Point(1303) = {3.0-dist2, 1.0, 0, res_tmp6};
Point(1304) = {3.0, 1.0, 0, res_tmp5};
Point(1305) = {3.0, 1.0+dist2, 0, res_tmp6};
Point(1306) = {3.0, 3.0-dist2, 0, res_tmp6};
Point(1307) = {3.0, 3.0, 0, res_tmp5};
Point(1308) = {3.0-dist2, 3.0, 0, res_tmp6};
Point(1309) = {1.0+dist2, 3.0, 0, res_tmp6};
Point(1310) = {1.0, 3.0, 0, res_tmp5};
Point(1311) = {1.0, 3.0-dist2, 0, res_tmp6};
Point(1312) = {1.0, 1.0+dist2, 0, res_tmp6};


Point(1401) = {4.0, 4.0, 0, res_tmp7};


Line(2001) = {1001, 1002};
Line(2002) = {1002, 1003};
Line(2003) = {1003, 1004};
Line(2004) = {1004, 1005};
Line(2005) = {1005, 1006};
Line(2006) = {1006, 1007};
Line(2007) = {1007, 1008};
Line(2008) = {1008, 1009};
Line(2009) = {1009, 1010};
Line(2010) = {1010, 1011};
Line(2011) = {1011, 1012};
Line(2012) = {1012, 1013};
Line(2013) = {1013, 1014};
Line(2014) = {1014, 1001};


Line(2101) = {1101, 1102};
Line(2102) = {1102, 1103};
Line(2103) = {1103, 1104};
Line(2104) = {1104, 1105};
Line(2105) = {1105, 1106};
Line(2106) = {1106, 1101};


Line(2202) = {1202, 1203};
Line(2205) = {1205, 1206};
Line(2208) = {1208, 1209};
Line(2211) = {1211, 1212};

Line(2301) = {1301, 1302};
Line(2302) = {1302, 1303};
Line(2303) = {1303, 1304};
Line(2304) = {1304, 1305};
Line(2305) = {1305, 1306};
Line(2306) = {1306, 1307};
Line(2307) = {1307, 1308};
Line(2308) = {1308, 1309};
Line(2309) = {1309, 1310};
Line(2310) = {1310, 1311};
Line(2311) = {1311, 1312};
Line(2312) = {1312, 1301};

Line(2401) = {1106, 1401};
Line(2402) = {1401, 1102};

Circle(2601) = {1212, 1301, 1202};
Circle(2602) = {1203, 1304, 1205};
Circle(2603) = {1206, 1307, 1208};
Circle(2604) = {1209, 1310, 1211};


Curve Loop(3207) = {2202, 2602, 2205, 2603, 2208, 2604, 2211, 2601};
Curve Loop(3208) = {2301:2312};
Plane Surface(4207) = {3207, 3208};
Curve Loop(3209) = {2105, 2106, 2101, 2102, 2103, 2104};
Curve Loop(3210) = {2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2001};
Plane Surface(4209) = {3210, 3209};
Curve Loop(3211) = {2102, 2103, 2104, 2105, 2401, 2402};
Plane Surface(4210) = {3211};
Curve Loop(3212) = {2101, -2402, -2401, 2106};
Plane Surface(4211) = {3212, 3207};


Physical Curve("outer",3000) = {
  2001:2014
};

Physical Curve("inner",3100) = {
  2301:2312
};

Physical Surface("mesh",4000) = {
  4207,
  4209,
  4210,
  4211
};
