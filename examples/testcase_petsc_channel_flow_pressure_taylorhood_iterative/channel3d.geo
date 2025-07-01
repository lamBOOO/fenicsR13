// ---------- Parameters & global mesh settings ----------
If(!Exists(p))
  p = 8;                     // mesh refinement exponent
EndIf

res  = 100;
Mesh.CharacteristicLengthMax = 1.6 * (2^(1/3))^(-p);
Mesh.MshFileVersion         = 2.0;

// Try a robust unstructured tetra algorithm (HXT, fall-back is Delaunay)
Mesh.Algorithm3D = 10;

// ---------- Bar dimensions ----------
length = 4;   // x-dimension
height = 1;   // y-dimension
depth  = 1;   // z-dimension  <-- NEW

// ---------- Points ----------
Point(1001) = {0,      0,      0,     res};
Point(1002) = {length, 0,      0,     res};
Point(1003) = {length, height, 0,     res};
Point(1004) = {0,      height, 0,     res};

Point(1005) = {0,      0,      depth, res};
Point(1006) = {length, 0,      depth, res};
Point(1007) = {length, height, depth, res};
Point(1008) = {0,      height, depth, res};

// ---------- Edges ----------
Line(2001) = {1001,1002};   // bottom face
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1001};

Line(2005) = {1005,1006};   // top face
Line(2006) = {1006,1007};
Line(2007) = {1007,1008};
Line(2008) = {1008,1005};

Line(2009) = {1001,1005};   // vertical edges
Line(2010) = {1002,1006};
Line(2011) = {1003,1007};
Line(2012) = {1004,1008};

// ---------- Line loops & surfaces ----------
Line Loop(3001) = {2001,2002,2003,2004};         Plane Surface(4001) = {3001}; // bottom
Line Loop(3002) = {2005,2006,2007,2008};         Plane Surface(4002) = {3002}; // top
Line Loop(3003) = {2001,2010,-2005,-2009};       Plane Surface(4003) = {3003}; // front
Line Loop(3004) = {2002,2011,-2006,-2010};       Plane Surface(4004) = {3004}; // right
Line Loop(3005) = {2003,2012,-2007,-2011};       Plane Surface(4005) = {3005}; // back
Line Loop(3006) = {2004,2009,-2008,-2012};       Plane Surface(4006) = {3006}; // left

// ---------- Volume ----------
Surface Loop(5000) = {4001,4002,4003,4004,4005,4006};
Volume(6000)       = {5000};

// ---------- Physical groups ----------
Physical Surface(4001) = {4001};
Physical Surface(4002) = {4002};
Physical Surface(4003) = {4003};
Physical Surface(4004) = {4004};
Physical Surface(4005) = {4005};
Physical Surface(4006) = {4006};
Physical Volume(6000)  = {6000};
