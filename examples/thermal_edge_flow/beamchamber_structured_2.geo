Mesh.MshFileVersion = 2.0;

res = 1000;
nnodes = 5;
progval = 2;
bumpval = 0.05;

exp = 1;
shift = 1 - 2^(-exp);
Printf("shift=%g", shift);

LB = 0.5^exp;
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
  Transfinite Line {2100 + 10*y + 1} = nnodes-1 Using Progression progval;
  Transfinite Line {2100 + 10*y + 2} = nnodes-1 Using Progression 1/progval;
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
  Transfinite Line {2200 + 10 + x} = nnodes-1 Using Progression progval;
  Transfinite Line {2200 + 20 + x} = nnodes-1 Using Progression 1/progval;
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




Line Loop(3000) = {
  2111,
  2112,
  2113,
  2114,
  2115,
  2116,
  2117,
  2218,
  2228,
  2238,
  2248,
  2258,
  2268,
  2278,
  -2178,
  -2177,
  -2176,
  -2175,
  -2174,
  -2173,
  -2172,
  -2171,
  -2271,
  -2261,
  -2251,
  -2241,
  -2231,
  -2221,
  -2211
};
Physical Curve("outer",3000) = {
  2111,
  2112,
  2113,
  2114,
  2115,
  2116,
  2117,
  2218,
  2228,
  2238,
  2248,
  2258,
  2268,
  2278,
  -2178,
  -2177,
  -2176,
  -2175,
  -2174,
  -2173,
  -2172,
  -2171,
  -2271,
  -2261,
  -2251,
  -2241,
  -2231,
  -2221,
  -2211
};



Line Loop(3100) = {
  2133,
  2134,
  2235,
  2245,
  -2154,
  -2153,
  -2243,
  -2233
};
Physical Curve("inner",3100) = {
  2133,
  2134,
  2235,
  2245,
  -2154,
  -2153,
  -2243,
  -2233
};

// Plane Surface(4000) = {
//   3111:3117,
//   3121:3127,
//   3131:3132,3135:3137,
//   3141:3142,3145:3147,
//   3151:3157,
//   3161:3167,
//   3171:3177
// };
Physical Surface("mesh",4000) = {
  3111:3117,
  3121:3127,
  3131:3132,3135:3137,
  3141:3142,3145:3147,
  3151:3157,
  3161:3167,
  3171:3177
};

