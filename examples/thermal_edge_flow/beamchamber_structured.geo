// Command line Parameters
If(!Exists(p))
  p = 3;
EndIf

// Settings
res = 100;
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;
nnodes = 10;
bumpval = 0.25;

Point(1001) = {0, 0, 0, res};
Point(1002) = {1, 0, 0, res};
Point(1003) = {3, 0, 0, res};
Point(1004) = {8, 0, 0, res};
Point(1005) = {8, 1, 0, res};
Point(1006) = {8, 3, 0, res};
Point(1007) = {8, 8, 0, res};
Point(1008) = {3, 8, 0, res};
Point(1009) = {1, 8, 0, res};
Point(1010) = {0, 8, 0, res};
Point(1011) = {0, 3, 0, res};
Point(1012) = {0, 1, 0, res};

Point(1101) = {1, 1, 0, res};
Point(1102) = {3, 1, 0, res};
Point(1103) = {3, 3, 0, res};
Point(1104) = {1, 3, 0, res};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1005};
Line(2005) = {1005,1006};
Line(2006) = {1006,1007};
Line(2007) = {1007,1008};
Line(2008) = {1008,1009};
Line(2009) = {1009,1010};
Line(2010) = {1010,1011};
Line(2011) = {1011,1012};
Line(2012) = {1012,1001};
Transfinite Line {2001} = nnodes Using Bump bumpval;
Transfinite Line {2002} = nnodes Using Bump bumpval;
Transfinite Line {2003} = nnodes Using Bump bumpval;
Transfinite Line {2004} = nnodes Using Bump bumpval;
Transfinite Line {2005} = nnodes Using Bump bumpval;
Transfinite Line {2006} = nnodes Using Bump bumpval;
Transfinite Line {2007} = nnodes Using Bump bumpval;
Transfinite Line {2008} = nnodes Using Bump bumpval;
Transfinite Line {2009} = nnodes Using Bump bumpval;
Transfinite Line {2010} = nnodes Using Bump bumpval;
Transfinite Line {2011} = nnodes Using Bump bumpval;
Transfinite Line {2012} = nnodes Using Bump bumpval;

Line(2101) = {1101,1102};
Line(2102) = {1102,1103};
Line(2103) = {1103,1104};
Line(2104) = {1104,1101};
Transfinite Line {2101} = nnodes Using Bump bumpval;
Transfinite Line {2102} = nnodes Using Bump bumpval;
Transfinite Line {2103} = nnodes Using Bump bumpval;
Transfinite Line {2104} = nnodes Using Bump bumpval;

Line(2201) = {1002,1101};
Line(2202) = {1003,1102};
Line(2203) = {1005,1102};
Line(2204) = {1006,1103};
Line(2205) = {1008,1103};
Line(2206) = {1009,1104};
Line(2207) = {1011,1104};
Line(2208) = {1012,1101};
Transfinite Line {2201} = nnodes Using Bump bumpval;
Transfinite Line {2202} = nnodes Using Bump bumpval;
Transfinite Line {2203} = nnodes Using Bump bumpval;
Transfinite Line {2204} = nnodes Using Bump bumpval;
Transfinite Line {2205} = nnodes Using Bump bumpval;
Transfinite Line {2206} = nnodes Using Bump bumpval;
Transfinite Line {2207} = nnodes Using Bump bumpval;
Transfinite Line {2208} = nnodes Using Bump bumpval;

Line Loop(3201) = {2001,2201,-2208,2012};
Plane Surface(4201) = {3201};
Transfinite Surface {4201} = {} Right;
//Recombine Surface {4200};

Line Loop(3202) = {2002,2202,-2101,-2201};
Plane Surface(4202) = {3202};
Transfinite Surface {4202};

Line Loop(3203) = {2003,2004,2203,-2202};
Plane Surface(4203) = {3203};
Transfinite Surface {4203};

Line Loop(3204) = {-2203,2005,2204,-2102};
Plane Surface(4204) = {3204};
Transfinite Surface {4204};

Line Loop(3205) = {-2204,2006,2007,2205};
Plane Surface(4205) = {3205};
Transfinite Surface {4205};

Line Loop(3206) = {-2103,-2205,2008,2206};
Plane Surface(4206) = {3206};
Transfinite Surface {4206};

Line Loop(3207) = {2207,-2206,2009,2010};
Plane Surface(4207) = {3207};
Transfinite Surface {4207};

Line Loop(3208) = {2208,-2104,-2207,2011};
Plane Surface(4208) = {3208};
Transfinite Surface {4208};

Line Loop(3000) = {2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012}; Physical Curve("outer",3000) = {2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012};

Line Loop(3100) = {2101,2102,2103,2104}; Physical Curve("inner",3100) = {2101,2102,2103,2104};

//Plane Surface(4000) = {3201,3202,3203,3204,3205,3206,3207,3208};
Physical Surface("mesh",4000) = {4201,4202,4203,4204,4205,4206,4207,4208};

//Recomine Surface {4000};
