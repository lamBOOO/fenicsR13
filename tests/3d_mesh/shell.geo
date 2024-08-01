// Command line Parameters
If(!Exists(p))
  p = 1;
EndIf

SetFactory("OpenCASCADE");

// Settings
res = 100; // Mesh resolution
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

// Define Parameters for the Spherical Surfaces
R1 = 0.5; // Inner radius
R2 = 2.0; // Outer radius

// Define Points at the center
Point(1) = {0, 0, 0, res};

// Define spherical surfaces
Sphere(1001) = {0, 0, 0, R1};  // Inner sphere
Sphere(1002) = {0, 0, 0, R2};  // Outer sphere


// Boolean operation to create the spherical shell
BooleanDifference(1234) = { Volume{1002}; Delete; }{ Volume{1001}; Delete; };
Physical Volume(10002) = {1234};

// Define Physical Surfaces
// Attention: BooleanDifference changes id's of the surfaces! check in GUI!
Physical Surface(10000) = {2};  // inner
Physical Surface(10001) = {1};  // outer
