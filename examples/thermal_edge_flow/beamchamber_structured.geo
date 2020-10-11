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

res = 1000;
progval = 2;
bumpval = 0.05;

shift = 1 - 2^(-exp);
shift05 = shift - 0.5;
Printf("shift=%g", shift);
Printf("shift05=%g", shift05);

Lbulk = 0.5-shift05;

N1 = (3+2*shift)/Lbulk + 1;
N2 = (0+2*shift05)/Lbulk + 1;
N3 = (0.5+shift05)/Lbulk + 1;
Printf("N1=%g", N1);
Printf("N2=%g", N2);
Printf("N3=%g", N3);

// TODO make loop

Point(110101) = {0.0, 0.0, 0, res};
Point(110102) = {0.5-shift05, 0.0, 0, res};
Point(110103) = {0.5+shift05, 0.0, 0, res};
Point(110104) = {1.0, 0.0, 0, res};
Point(110105) = {2.0-shift, 0.0, 0, res};
Point(110106) = {2.0, 0.0, 0, res};
Point(110107) = {2.0+shift, 0.0, 0, res};
Point(110108) = {3.0, 0.0, 0, res};
Point(110109) = {4.0-shift, 0.0, 0, res};
Point(110110) = {7.0+shift, 0.0, 0, res};
Point(110111) = {8.0, 0.0, 0, res};

Point(110201) = {0.0, 0.5-shift05, 0, res};
Point(110202) = {0.5-shift05, 0.5-shift05, 0, res};
Point(110203) = {0.5+shift05, 0.5-shift05, 0, res};
Point(110204) = {1.0, 0.5-shift05, 0, res};
Point(110205) = {2.0-shift, 0.5-shift05, 0, res};
Point(110206) = {2.0, 0.5-shift05, 0, res};
Point(110207) = {2.0+shift, 0.5-shift05, 0, res};
Point(110208) = {3.0, 0.5-shift05, 0, res};
Point(110209) = {4.0-shift, 0.5-shift05, 0, res};
Point(110210) = {7.0+shift, 0.5-shift05, 0, res};
Point(110211) = {8.0, 0.5-shift05, 0, res};

Point(110301) = {0.0, 0.5+shift05, 0, res};
Point(110302) = {0.5-shift05, 0.5+shift05, 0, res};
Point(110303) = {0.5+shift05, 0.5+shift05, 0, res};
Point(110304) = {1.0, 0.5+shift05, 0, res};
Point(110305) = {2.0-shift, 0.5+shift05, 0, res};
Point(110306) = {2.0, 0.5+shift05, 0, res};
Point(110307) = {2.0+shift, 0.5+shift05, 0, res};
Point(110308) = {3.0, 0.5+shift05, 0, res};
Point(110309) = {4.0-shift, 0.5+shift05, 0, res};
Point(110310) = {7.0+shift, 0.5+shift05, 0, res};
Point(110311) = {8.0, 0.5+shift05, 0, res};

Point(110401) = {0.0, 1.0, 0, res};
Point(110402) = {0.5-shift05, 1.0, 0, res};
Point(110403) = {0.5+shift05, 1.0, 0, res};
Point(110404) = {1.0, 1.0, 0, res};
Point(110405) = {2.0-shift, 1.0, 0, res};
Point(110406) = {2.0, 1.0, 0, res};
Point(110407) = {2.0+shift, 1.0, 0, res};
Point(110408) = {3.0, 1.0, 0, res};
Point(110409) = {4.0-shift, 1.0, 0, res};
Point(110410) = {7.0+shift, 1.0, 0, res};
Point(110411) = {8.0, 1.0, 0, res};

Point(110501) = {0.0, 2.0-shift, 0, res};
Point(110502) = {0.5-shift05, 2.0-shift, 0, res};
Point(110503) = {0.5+shift05, 2.0-shift, 0, res};
Point(110504) = {1.0, 2.0-shift, 0, res};
Point(110505) = {2.0-shift, 2.0-shift, 0, res};
Point(110506) = {2.0, 2.0-shift, 0, res};
Point(110507) = {2.0+shift, 2.0-shift, 0, res};
Point(110508) = {3.0, 2.0-shift, 0, res};
Point(110509) = {4.0-shift, 2.0-shift, 0, res};
Point(110510) = {7.0+shift, 2.0-shift, 0, res};
Point(110511) = {8.0, 2.0-shift, 0, res};

Point(110601) = {0.0, 2.0, 0, res};
Point(110602) = {0.5-shift05, 2.0, 0, res};
Point(110603) = {0.5+shift05, 2.0, 0, res};
Point(110604) = {1.0, 2.0, 0, res};
Point(110605) = {2.0-shift, 2.0, 0, res};
Point(110606) = {2.0, 2.0, 0, res};
Point(110607) = {2.0+shift, 2.0, 0, res};
Point(110608) = {3.0, 2.0, 0, res};
Point(110609) = {4.0-shift, 2.0, 0, res};
Point(110610) = {7.0+shift, 2.0, 0, res};
Point(110611) = {8.0, 2.0, 0, res};

Point(110701) = {0.0, 2.0+shift, 0, res};
Point(110702) = {0.5-shift05, 2.0+shift, 0, res};
Point(110703) = {0.5+shift05, 2.0+shift, 0, res};
Point(110704) = {1.0, 2.0+shift, 0, res};
Point(110705) = {2.0-shift, 2.0+shift, 0, res};
Point(110706) = {2.0, 2.0+shift, 0, res};
Point(110707) = {2.0+shift, 2.0+shift, 0, res};
Point(110708) = {3.0, 2.0+shift, 0, res};
Point(110709) = {4.0-shift, 2.0+shift, 0, res};
Point(110710) = {7.0+shift, 2.0+shift, 0, res};
Point(110711) = {8.0, 2.0+shift, 0, res};

Point(110801) = {0.0, 3.0, 0, res};
Point(110802) = {0.5-shift05, 3.0, 0, res};
Point(110803) = {0.5+shift05, 3.0, 0, res};
Point(110804) = {1.0, 3.0, 0, res};
Point(110805) = {2.0-shift, 3.0, 0, res};
Point(110806) = {2.0, 3.0, 0, res};
Point(110807) = {2.0+shift, 3.0, 0, res};
Point(110808) = {3.0, 3.0, 0, res};
Point(110809) = {4.0-shift, 3.0, 0, res};
Point(110810) = {7.0+shift, 3.0, 0, res};
Point(110811) = {8.0, 3.0, 0, res};

Point(110901) = {0.0, 4.0-shift, 0, res};
Point(110902) = {0.5-shift05, 4.0-shift, 0, res};
Point(110903) = {0.5+shift05, 4.0-shift, 0, res};
Point(110904) = {1.0, 4.0-shift, 0, res};
Point(110905) = {2.0-shift, 4.0-shift, 0, res};
Point(110906) = {2.0, 4.0-shift, 0, res};
Point(110907) = {2.0+shift, 4.0-shift, 0, res};
Point(110908) = {3.0, 4.0-shift, 0, res};
Point(110909) = {4.0-shift, 4.0-shift, 0, res};
Point(110910) = {7.0+shift, 4.0-shift, 0, res};
Point(110911) = {8.0, 4.0-shift, 0, res};

Point(111001) = {0.0, 7.0+shift, 0, res};
Point(111002) = {0.5-shift05, 7.0+shift, 0, res};
Point(111003) = {0.5+shift05, 7.0+shift, 0, res};
Point(111004) = {1.0, 7.0+shift, 0, res};
Point(111005) = {2.0-shift, 7.0+shift, 0, res};
Point(111006) = {2.0, 7.0+shift, 0, res};
Point(111007) = {2.0+shift, 7.0+shift, 0, res};
Point(111008) = {3.0, 7.0+shift, 0, res};
Point(111009) = {4.0-shift, 7.0+shift, 0, res};
Point(111010) = {7.0+shift, 7.0+shift, 0, res};
Point(111011) = {8.0, 7.0+shift, 0, res};

Point(111101) = {0.0, 8.0, 0, res};
Point(111102) = {0.5-shift05, 8.0, 0, res};
Point(111103) = {0.5+shift05, 8.0, 0, res};
Point(111104) = {1.0, 8.0, 0, res};
Point(111105) = {2.0-shift, 8.0, 0, res};
Point(111106) = {2.0, 8.0, 0, res};
Point(111107) = {2.0+shift, 8.0, 0, res};
Point(111108) = {3.0, 8.0, 0, res};
Point(111109) = {4.0-shift, 8.0, 0, res};
Point(111110) = {7.0+shift, 8.0, 0, res};
Point(111111) = {8.0, 8.0, 0, res};





For y In {1:11}
  For x In {1:10}
    Line(210000 + 100*y + x) = {110000 + 100*y + x, 110000 + 100*y + x + 1};
  EndFor
EndFor

For y In {1:10}
  For x In {1:11}
    Line(220000 + 100*y + x) = {110000 + 100*y + x, 110100 + 100*y + x};
  EndFor
EndFor




For y In {1:11}
  Printf("t=%g", y);
  Transfinite Line {210000 + 100*y + 1} = nnodes Using Progression progval;
  Transfinite Line {210000 + 100*y + 2} = N2;
  Transfinite Line {210000 + 100*y + 3} = nnodes Using Progression 1/progval;
  Transfinite Line {210000 + 100*y + 4} = nnodes Using Progression progval;
  Transfinite Line {210000 + 100*y + 5} = N3;
  Transfinite Line {210000 + 100*y + 6} = N3;
  Transfinite Line {210000 + 100*y + 7} = nnodes Using Progression 1/progval;
  Transfinite Line {210000 + 100*y + 8} = nnodes Using Progression progval;
  Transfinite Line {210000 + 100*y + 9} = N1;
  Transfinite Line {210000 + 100*y + 10} = nnodes Using Progression 1/progval;
EndFor

For x In {1:11}
  Printf("t=%g", x);
  Transfinite Line {220000 + 100 + x} = nnodes Using Progression progval;
  Transfinite Line {220000 + 200 + x} = N2;
  Transfinite Line {220000 + 300 + x} = nnodes Using Progression 1/progval;
  Transfinite Line {220000 + 400 + x} = nnodes Using Progression progval;
  Transfinite Line {220000 + 500 + x} = N3;
  Transfinite Line {220000 + 600 + x} = N3;
  Transfinite Line {220000 + 700 + x} = nnodes Using Progression 1/progval;
  Transfinite Line {220000 + 800 + x} = nnodes Using Progression progval;
  Transfinite Line {220000 + 900 + x} = N1;
  Transfinite Line {220000 + 1000 + x} = nnodes Using Progression 1/progval;
EndFor





For y In {1:10}
  Printf("t=%g", y);
  For x In {1:10}
    // If (!((x==3 || x==4) && (y==3 || y==4)))
    Printf("u=%g", x);
    Line Loop(310000 + 100*y + x) = {210000 + 100*y + x, 220001 + 100*y + x, -(210100 + 100*y + x), -(220000 + 100*y + x)};
    Plane Surface(310000 + 100*y + x) = {310000 + 100*y + x};
    Transfinite Surface {310000 + 100*y + x} = {} AlternateLeft;
    // EndIf
  EndFor
EndFor




Line Loop(3000) = {
  210101:210110,
  220111:221011:100,
  -211110:-211101,
  -221001:-220101:100
};
Physical Curve("outer",3000) = {
  210101:210110,
  220111:221011:100,
  -211110:-211101,
  -221001:-220101:100
};

Line Loop(3100) = {
  210404:210407,
  220408:220708:100,
  -210807:-210804,
  -220704:-220404:100
};
Physical Curve("inner",3100) = {
  210404:210407,
  220408:220708:100,
  -210807:-210804,
  -220704:-220404:100
};

Physical Surface("mesh",4000) = {
  310101:310110,
  310201:310210,
  310301:310310,
  310401:310403,310408:310410,
  310501:310503,310508:310510,
  310601:310603,310608:310610,
  310701:310703,310708:310710,
  310801:310810,
  310901:310910,
  311001:311010
};

