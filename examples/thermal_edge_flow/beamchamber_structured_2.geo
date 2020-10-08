Mesh.MshFileVersion = 2.0;

res = 100;
nnodes = 5;
progval = 2;
bumpval = 0.05;

LB = 0.5;
NB = 3/LB + 1;

// TODO make loop

Point(1111) = {0.0, 0.0, 0, res};
Point(1112) = {0.5, 0.0, 0, res};
Point(1113) = {1.0, 0.0, 0, res};
Point(1114) = {2.0, 0.0, 0, res};
Point(1115) = {3.0, 0.0, 0, res};
Point(1116) = {4.0, 0.0, 0, res};
Point(1117) = {7.0, 0.0, 0, res};
Point(1118) = {8.0, 0.0, 0, res};

Point(1121) = {0.0, 0.5, 0, res};
Point(1122) = {0.5, 0.5, 0, res};
Point(1123) = {1.0, 0.5, 0, res};
Point(1124) = {2.0, 0.5, 0, res};
Point(1125) = {3.0, 0.5, 0, res};
Point(1126) = {4.0, 0.5, 0, res};
Point(1127) = {7.0, 0.5, 0, res};
Point(1128) = {8.0, 0.5, 0, res};

Point(1131) = {0.0, 1.0, 0, res};
Point(1132) = {0.5, 1.0, 0, res};
Point(1133) = {1.0, 1.0, 0, res};
Point(1134) = {2.0, 1.0, 0, res};
Point(1135) = {3.0, 1.0, 0, res};
Point(1136) = {4.0, 1.0, 0, res};
Point(1137) = {7.0, 1.0, 0, res};
Point(1138) = {8.0, 1.0, 0, res};

Point(1141) = {0.0, 2.0, 0, res};
Point(1142) = {0.5, 2.0, 0, res};
Point(1143) = {1.0, 2.0, 0, res};
Point(1144) = {2.0, 2.0, 0, res};
Point(1145) = {3.0, 2.0, 0, res};
Point(1146) = {4.0, 2.0, 0, res};
Point(1147) = {7.0, 2.0, 0, res};
Point(1148) = {8.0, 2.0, 0, res};

Point(1151) = {0.0, 3.0, 0, res};
Point(1152) = {0.5, 3.0, 0, res};
Point(1153) = {1.0, 3.0, 0, res};
Point(1154) = {2.0, 3.0, 0, res};
Point(1155) = {3.0, 3.0, 0, res};
Point(1156) = {4.0, 3.0, 0, res};
Point(1157) = {7.0, 3.0, 0, res};
Point(1158) = {8.0, 3.0, 0, res};

Point(1161) = {0.0, 4.0, 0, res};
Point(1162) = {0.5, 4.0, 0, res};
Point(1163) = {1.0, 4.0, 0, res};
Point(1164) = {2.0, 4.0, 0, res};
Point(1165) = {3.0, 4.0, 0, res};
Point(1166) = {4.0, 4.0, 0, res};
Point(1167) = {7.0, 4.0, 0, res};
Point(1168) = {8.0, 4.0, 0, res};

Point(1171) = {0.0, 7.0, 0, res};
Point(1172) = {0.5, 7.0, 0, res};
Point(1173) = {1.0, 7.0, 0, res};
Point(1174) = {2.0, 7.0, 0, res};
Point(1175) = {3.0, 7.0, 0, res};
Point(1176) = {4.0, 7.0, 0, res};
Point(1177) = {7.0, 7.0, 0, res};
Point(1178) = {8.0, 7.0, 0, res};

Point(1181) = {0.0, 8.0, 0, res};
Point(1182) = {0.5, 8.0, 0, res};
Point(1183) = {1.0, 8.0, 0, res};
Point(1184) = {2.0, 8.0, 0, res};
Point(1185) = {3.0, 8.0, 0, res};
Point(1186) = {4.0, 8.0, 0, res};
Point(1187) = {7.0, 8.0, 0, res};
Point(1188) = {8.0, 8.0, 0, res};





For x In {1:8}
  For y In {1:7}
    Line(2100 + 10*x + y) = {1100 + 10*x + y, 1100 + 10*x + y + 1};
  EndFor
EndFor

For y In {1:7}
  For x In {1:8}
    Line(2200 + 10*y + x) = {1100 + 10*y + x, 1110 + 10*y + x};
  EndFor
EndFor




For y In {1:8}
  Printf("t=%g", y);
  Transfinite Line {2100 + 10*y + 1} = nnodes Using Progression progval;
  Transfinite Line {2100 + 10*y + 2} = nnodes Using Progression 1/progval;
  If (y!=4)
    Transfinite Line {2100 + 10*y + 3} = nnodes Using Progression progval;
  EndIf
  If (y!=4)
    Transfinite Line {2100 + 10*y + 4} = nnodes Using Progression 1/progval;
  EndIf
  Transfinite Line {2100 + 10*y + 5} = nnodes Using Progression progval;
  Transfinite Line {2100 + 10*y + 6} = NB;
  Transfinite Line {2100 + 10*y + 7} = nnodes Using Progression 1/progval;
EndFor

For x In {1:8}
  Printf("t=%g", x);
  Transfinite Line {2200 + 10 + x} = nnodes Using Progression progval;
  Transfinite Line {2200 + 20 + x} = nnodes Using Progression 1/progval;
  If (x!=4)
    Transfinite Line {2200 + 30 + x} = nnodes Using Progression progval;
  EndIf
  If (x!=4)
    Transfinite Line {2200 + 40 + x} = nnodes Using Progression 1/progval;
  EndIf
  Transfinite Line {2200 + 50 + x} = nnodes Using Progression progval;
  Transfinite Line {2200 + 60 + x} = NB;
  Transfinite Line {2200 + 70 + x} = nnodes Using Progression 1/progval;
EndFor





For y In {1:7}
  Printf("t=%g", y);
  For x In {1:7}
    If (!((x==3 || x==4) && (y==3 || y==4)))
      Printf("u=%g", x);
      Line Loop(3100 + 10*y + x) = {2100 + 10*y + x, 2201 + 10*y + x, -(2110 + 10*y + x), -(2200 + 10*y + x)};
      Plane Surface(3100 + 10*y + x) = {3100 + 10*y + x};
      Transfinite Surface {3100 + 10*y + x};
    EndIf
  EndFor
EndFor










// // outer
// p01 = newp; Point(p01) = {0.0, 0.0, 0, res};
// p02 = newp; Point(p02) = {0.5, 0.0, 0, res};
// p03 = newp; Point(p03) = {1.0, 0.0, 0, res};
// p04 = newp; Point(p04) = {2.0, 0.0, 0, res};
// p05 = newp; Point(p05) = {3.0, 0.0, 0, res};
// p06 = newp; Point(p06) = {4.0, 0.0, 0, res};
// p07 = newp; Point(p07) = {7.0, 0.0, 0, res};

// p11 = newp; Point(p11) = {8.0, 0.0, 0, res};
// p12 = newp; Point(p12) = {8.0, 0.5, 0, res};
// p13 = newp; Point(p13) = {8.0, 1.0, 0, res};
// p14 = newp; Point(p14) = {8.0, 2.0, 0, res};
// p15 = newp; Point(p15) = {8.0, 3.0, 0, res};
// p16 = newp; Point(p16) = {8.0, 4.0, 0, res};
// p17 = newp; Point(p17) = {8.0, 7.0, 0, res};

// p21 = newp; Point(p21) = {8.0, 8.0, 0, res};
// p22 = newp; Point(p22) = {7.0, 8.0, 0, res};
// p23 = newp; Point(p23) = {4.0, 8.0, 0, res};
// p24 = newp; Point(p24) = {3.0, 8.0, 0, res};
// p25 = newp; Point(p25) = {2.0, 8.0, 0, res};
// p26 = newp; Point(p26) = {1.0, 8.0, 0, res};
// p27 = newp; Point(p27) = {0.5, 8.0, 0, res};

// p31 = newp; Point(p31) = {0.0, 8.0, 0, res};
// p32 = newp; Point(p32) = {0.0, 7.0, 0, res};
// p33 = newp; Point(p33) = {0.0, 3.0, 0, res};
// p34 = newp; Point(p34) = {0.0, 4.0, 0, res};
// p35 = newp; Point(p35) = {0.0, 2.0, 0, res};
// p36 = newp; Point(p36) = {0.0, 1.0, 0, res};
// p37 = newp; Point(p37) = {0.0, 0.5, 0, res};

// // // inner
// // p41 = newp; Point(p41) = {1.0, 1.0, 0, res};
// // p42 = newp; Point(p42) = {2.0, 1.0, 0, res};
// // p43 = newp; Point(p43) = {3.0, 1.0, 0, res};
// // p44 = newp; Point(p44) = {3.0, 2.0, 0, res};
// // p45 = newp; Point(p45) = {3.0, 3.0, 0, res};
// // p46 = newp; Point(p46) = {2.0, 3.0, 0, res};
// // p47 = newp; Point(p47) = {1.0, 3.0, 0, res};
// // p48 = newp; Point(p48) = {1.0, 2.0, 0, res};

// // // connection points
// // p51 = newp; Point(p51) = {1.0, 0.5, 0, res};
// // p51 = newp; Point(p51) = {2.0, 0.5, 0, res};
// // p51 = newp; Point(p51) = {3.0, 0.5, 0, res};
// // p51 = newp; Point(p51) = {4.0, 1.0, 0, res};

// l01 = newl; Line(l01) = {p01, p02};
// l02 = newl; Line(l02) = {p02, p03};

// l04 = newl; Line(l04) = {p03, p04};
// l05 = newl; Line(l05) = {p04, p05};
// l06 = newl; Line(l06) = {p05, p06};
// l07 = newl; Line(l07) = {p06, p07};
// l08 = newl; Line(l08) = {p07, p11};

// l11 = newl; Line(l11) = {p11, p12};
// l12 = newl; Line(l12) = {p12, p13};

// l14 = newl; Line(l14) = {p13, p14};
// l15 = newl; Line(l15) = {p14, p15};
// l16 = newl; Line(l16) = {p15, p16};
// l17 = newl; Line(l17) = {p16, p17};
// l18 = newl; Line(l18) = {p17, p21};

// l21 = newl; Line(l21) = {p21, p22};
// l22 = newl; Line(l22) = {p22, p23};

// l24 = newl; Line(l24) = {p23, p24};
// l25 = newl; Line(l25) = {p24, p25};
// l26 = newl; Line(l26) = {p25, p26};
// l27 = newl; Line(l27) = {p26, p27};
// l28 = newl; Line(l28) = {p27, p31};

// l31 = newl; Line(l31) = {p31, p32};
// l32 = newl; Line(l32) = {p32, p33};

// l34 = newl; Line(l34) = {p33, p34};
// l35 = newl; Line(l35) = {p34, p35};
// l36 = newl; Line(l36) = {p35, p36};
// l37 = newl; Line(l37) = {p36, p37};
// l38 = newl; Line(l38) = {p37, p01};

// Transfinite Line {l01} = nnodes Using Progression progval;
// Transfinite Line {l02} = nnodes Using Progression progval;
// // Transfinite Line {l03} = nnodes Using Progression progval;
// Transfinite Line {l04} = nnodes Using Progression progval;
// Transfinite Line {l05} = nnodes Using Progression progval;
// Transfinite Line {l06} = nnodes Using Progression progval;
// Transfinite Line {l07} = nnodes Using Progression progval;
// Transfinite Line {l08} = nnodes Using Progression progval;
// Transfinite Line {l11} = nnodes Using Progression progval;
// Transfinite Line {l12} = nnodes Using Progression progval;
// // Transfinite Line {l13} = nnodes Using Progression progval;
// Transfinite Line {l14} = nnodes Using Progression progval;
// Transfinite Line {l15} = nnodes Using Progression progval;
// Transfinite Line {l16} = nnodes Using Progression progval;
// Transfinite Line {l17} = nnodes Using Progression progval;
// Transfinite Line {l18} = nnodes Using Progression progval;
// Transfinite Line {l21} = nnodes Using Progression progval;
// Transfinite Line {l22} = nnodes Using Progression progval;
// // Transfinite Line {l23} = nnodes Using Progression progval;
// Transfinite Line {l24} = nnodes Using Progression progval;
// Transfinite Line {l25} = nnodes Using Progression progval;
// Transfinite Line {l26} = nnodes Using Progression progval;
// Transfinite Line {l27} = nnodes Using Progression progval;
// Transfinite Line {l28} = nnodes Using Progression progval;
// Transfinite Line {l31} = nnodes Using Progression progval;
// Transfinite Line {l32} = nnodes Using Progression progval;
// // Transfinite Line {l33} = nnodes Using Progression progval;
// Transfinite Line {l34} = nnodes Using Progression progval;
// Transfinite Line {l35} = nnodes Using Progression progval;
// Transfinite Line {l36} = nnodes Using Progression progval;
// Transfinite Line {l37} = nnodes Using Progression progval;
// Transfinite Line {l38} = nnodes Using Progression progval;






// Point(1002) = {1, 0, 0, res};
// Point(1003) = {3, 0, 0, res};
// Point(1004) = {8, 0, 0, res};
// Point(1005) = {8, 1, 0, res};
// Point(1006) = {8, 3, 0, res};
// Point(1007) = {8, 8, 0, res};
// Point(1008) = {3, 8, 0, res};
// Point(1009) = {1, 8, 0, res};
// Point(1010) = {0, 8, 0, res};
// Point(1011) = {0, 3, 0, res};
// Point(1012) = {0, 1, 0, res};

// Point(1101) = {1, 1, 0, res};
// Point(1102) = {3, 1, 0, res};
// Point(1103) = {3, 3, 0, res};
// Point(1104) = {1, 3, 0, res};

// Line(2001) = {1001,1002};
// Line(2002) = {1002,1003};
// Line(2003) = {1003,1004};
// Line(2004) = {1004,1005};
// Line(2005) = {1005,1006};
// Line(2006) = {1006,1007};
// Line(2007) = {1007,1008};
// Line(2008) = {1008,1009};
// Line(2009) = {1009,1010};
// Line(2010) = {1010,1011};
// Line(2011) = {1011,1012};
// Line(2012) = {1012,1001};
// Transfinite Line {2001} = nnodes Using Bump bumpval;
// Transfinite Line {2002} = nnodes Using Bump bumpval;
// Transfinite Line {2003} = nnodes Using Bump bumpval;
// Transfinite Line {2004} = nnodes Using Bump bumpval;
// Transfinite Line {2005} = nnodes Using Bump bumpval;
// Transfinite Line {2006} = nnodes Using Bump bumpval;
// Transfinite Line {2007} = nnodes Using Bump bumpval;
// Transfinite Line {2008} = nnodes Using Bump bumpval;
// Transfinite Line {2009} = nnodes Using Bump bumpval;
// Transfinite Line {2010} = nnodes Using Bump bumpval;
// Transfinite Line {2011} = nnodes Using Bump bumpval;
// Transfinite Line {2012} = nnodes Using Bump bumpval;

// Line(2101) = {1101,1102};
// Line(2102) = {1102,1103};
// Line(2103) = {1103,1104};
// Line(2104) = {1104,1101};
// Transfinite Line {2101} = nnodes Using Bump bumpval;
// Transfinite Line {2102} = nnodes Using Bump bumpval;
// Transfinite Line {2103} = nnodes Using Bump bumpval;
// Transfinite Line {2104} = nnodes Using Bump bumpval;

// Line(2201) = {1002,1101};
// Line(2202) = {1003,1102};
// Line(2203) = {1005,1102};
// Line(2204) = {1006,1103};
// Line(2205) = {1008,1103};
// Line(2206) = {1009,1104};
// Line(2207) = {1011,1104};
// Line(2208) = {1012,1101};
// Transfinite Line {2201} = nnodes Using Bump bumpval;
// Transfinite Line {2202} = nnodes Using Bump bumpval;
// Transfinite Line {2203} = nnodes Using Bump bumpval;
// Transfinite Line {2204} = nnodes Using Bump bumpval;
// Transfinite Line {2205} = nnodes Using Bump bumpval;
// Transfinite Line {2206} = nnodes Using Bump bumpval;
// Transfinite Line {2207} = nnodes Using Bump bumpval;
// Transfinite Line {2208} = nnodes Using Bump bumpval;

// Line Loop(3201) = {2001,2201,-2208,2012};
// Plane Surface(4201) = {3201};
// Transfinite Surface {4201} = {} Right;
// //Recombine Surface {4200};

// Line Loop(3202) = {2002,2202,-2101,-2201};
// Plane Surface(4202) = {3202};
// Transfinite Surface {4202};

// Line Loop(3203) = {2003,2004,2203,-2202};
// Plane Surface(4203) = {3203};
// Transfinite Surface {4203};

// Line Loop(3204) = {-2203,2005,2204,-2102};
// Plane Surface(4204) = {3204};
// Transfinite Surface {4204};

// Line Loop(3205) = {-2204,2006,2007,2205};
// Plane Surface(4205) = {3205};
// Transfinite Surface {4205};

// Line Loop(3206) = {-2103,-2205,2008,2206};
// Plane Surface(4206) = {3206};
// Transfinite Surface {4206};

// Line Loop(3207) = {2207,-2206,2009,2010};
// Plane Surface(4207) = {3207};
// Transfinite Surface {4207};

// Line Loop(3208) = {2208,-2104,-2207,2011};
// Plane Surface(4208) = {3208};
// Transfinite Surface {4208};

// Line Loop(3000) = {2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012}; Physical Curve("outer",3000) = {2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012};

// Line Loop(3100) = {2101,2102,2103,2104}; Physical Curve("inner",3100) = {2101,2102,2103,2104};

// //Plane Surface(4000) = {3201,3202,3203,3204,3205,3206,3207,3208};
// Physical Surface("mesh",4000) = {4201,4202,4203,4204,4205,4206,4207,4208};

// //Recomine Surface {4000};
//+

