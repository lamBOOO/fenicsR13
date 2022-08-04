// Command line Parameters
If(!Exists(p))
  p = 7;
EndIf

// Settings
res = 100;
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

length_top=0.1;
increase_bot=1/Sqrt(3);
height=1;

Point(1001) = {-length_top-increase_bot,      0,      0, res};
Point(1002) = {length_top+increase_bot, 0,      0, res};
Point(1003) = {length_top, height, 0, res};
Point(1004) = {-length_top,      height, 0, res};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1001};

Line Loop(3001) = {2001}; Physical Curve("bot",3001) = {2001};
Line Loop(3002) = {2002}; Physical Curve("right",3002) = {2002};
Line Loop(3003) = {2003}; Physical Curve("top",3003) = {2003};
Line Loop(3004) = {2004}; Physical Curve("left",3004) = {2004};

numvert = 22;
numhori = 5;
distvert = height / (numvert-1);
disthori = 2*length_top / (numhori-1);
rad = 0.01;

holes = {};
circs = {};

For i In {1:numvert}
For t In {1:numhori+(i-1)}
  Printf("%f,%f", i, t);
  Point(101 + 10*t + 1000*i) = {-length_top     - disthori/2*(i-1) + (t-1)*disthori, height-distvert-(i-1)*disthori*Sqrt(3)/2, 0, res};
  Point(102 + 10*t + 1000*i) = {-length_top+rad - disthori/2*(i-1) + (t-1)*disthori, height-distvert-(i-1)*disthori*Sqrt(3)/2, 0, res};
  Point(103 + 10*t + 1000*i) = {-length_top-rad - disthori/2*(i-1) + (t-1)*disthori, height-distvert-(i-1)*disthori*Sqrt(3)/2, 0, res};

  Circle(301 + 10*t + 1000*i) = {102 + 10*t + 1000*i,101+10*t + 1000*i,103+10*t + 1000*i};
  Circle(302 + 10*t + 1000*i) = {103 + 10*t + 1000*i,101+10*t + 1000*i,102+10*t + 1000*i};
  Line Loop(401 + 10*t + 1000*i) = {301+10*t + 1000*i,302+10*t + 1000*i};

  holes += 401 + 10*t + 1000*i;
  circs += 301 + 10*t + 1000*i;
  circs += 302 + 10*t + 1000*i;
EndFor
EndFor

Printf("%s", holes);

Line Loop(9999) = {1411,1421};

Physical Curve("holes",3005) = {circs[]};

test = {1411,1421};

Printf("%f", test[0]);

Plane Surface(4000) = {
  3001:3004,
  holes[]
}; Physical Surface("mesh",4000) = {4000};
