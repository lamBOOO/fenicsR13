#!/usr/bin/env wolframscript
(* ::Package:: *)

(* this requires a specific problem setup for R13 *)


Print["*******"];
Print["Analytical solution of " <> theory <> " equations."];

Get[theory <> "AnalyticalSolution.wls"];


Print["Solving for integration constants."];

Kn = \[Tau];


solbc = Solve[(BCeqns) == 0, coefs][[1]];


(*Export solution to Cpp*)

(* cppRules = {
   "Log" -> "log",
    "Cos" -> "cos",
    "Sin" -> "sin",
    "Power" -> "pow",
    "Sqrt" -> "sqrt"
};

temperaturec[x_, y_] = temperature[Sqrt[x^2 + y^2], ArcTan[x, y]] /. scaling /. solbc;

temperaturecCode =
  ToString[CForm[temperaturec[x, y]]];

temperaturecCode = StringReplace[temperaturecCode, cppRules];

Print["temperaturecCode: ", temperaturecCode];

{uxc[x_, y_], uyc[x_, y_]} = {ur[Sqrt[x^2 + y^2], ArcTan[x, y]], 
                            uTheta[Sqrt[x^2 + y^2], ArcTan[x, y]]}.basis[ArcTan[x, y]] /. scaling /. solbc;
uxcCode =
  ToString[CForm[uxc[x, y]]];
  
uycCode =
  ToString[CForm[uyc[x, y]]];
  
uxcCode = StringReplace[uxcCode, cppRules];
  
uycCode = StringReplace[uycCode, cppRules];
  
Print["uxcCode: ", uxcCode];
Print["uycCode: ", uycCode];

pressurec[x_, y_] = p[Sqrt[x^2 + y^2], ArcTan[x, y]] /. scaling /. solbc;

pressurecCode =
  ToString[CForm[pressurec[x, y]]];

pressurecCode = StringReplace[pressurecCode, cppRules];

Print["pressurecCode: ", pressurecCode];

rxy = Sqrt[x^2 + y^2];
thxy = ArcTan[x, y];

SigmaC[x_, y_] =
  Simplify[
    Transpose[basis[thxy]] .
    Sigma[rxy, thxy][[1 ;; 2, 1 ;; 2]] .
    basis[thxy]
    /. scaling /. solbc
  ];
  
sigmaxxCode = ToString[CForm[SigmaC[x, y][[1, 1]]]];
sigmaxyCode = ToString[CForm[SigmaC[x, y][[1, 2]]]];
sigmayxCode = ToString[CForm[SigmaC[x, y][[2, 1]]]];
sigmayyCode = ToString[CForm[SigmaC[x, y][[2, 2]]]];

sigmaxxCode = StringReplace[sigmaxxCode, cppRules];
sigmaxyCode = StringReplace[sigmaxyCode, cppRules];
sigmayxCode = StringReplace[sigmayxCode, cppRules];
sigmayyCode = StringReplace[sigmayyCode, cppRules];

Print["sigmaxxCode: ", sigmaxxCode];
Print["sigmaxyCode: ", sigmaxyCode];
Print["sigmayxCode: ", sigmayxCode];
Print["sigmayyCode: ", sigmayyCode]; *)

(* cppRules = {
   "C(0)" -> C0,
   "C(1)" -> C2,
   C[3] -> CC2,
   C[4] -> CC3,
   C[5] -> CC4,
   C[6] -> CC5,
   C[7] -> CC6,
   C[8] -> CC7,
   C[9] -> CC8,
   C[10] -> CC9,

   C[1, "I", 1] -> CKK1,
   C[1, "K", 1] -> CKK2,
   C[2, "I", 1] -> CKK3,
   C[2, "K", 1] -> CKK4
   "CC" -> "C_",
   "CKK" -> "CK_",

   "BI(0," -> "I_0(",
   "BI(1," -> "I_1(",
   "BI(2," -> "I_2(",

   "BK(0," -> "K_0(",
   "BK(1," -> "K_1(",
   "BK(2," -> "K_2(",

   "Log" -> "std::log"
}; *)

Print["// Lambdas"];

Do[
  Print[
    "double lambda_", i, " = ",
    ToString[CForm[N[lambdaList[[i]], 16]]],
    ";"
  ],
  {i, 1, Length[lambdaList]}
];

Print[""];
Print["// Scaled constants"];

(* 1. Build coefficient lists *)
bulkCoefs = C /@ Range[nC];

knCoefs = Flatten[
  Table[
    {C[m, "I", i], C[m, "K", i]},
    {m, 1, Length[lambdaList]},
    {i, 1, nspcL[[m]]}
  ]
];

(* Bulk constants with scaling *)
scaledBulkRules = {
  C[1] -> C[1]*Kn^-1,
  C[4] -> C[4]*Kn^2,
  C[9] -> C[9]*Kn^2,
  C[11] -> C[11]*Kn^2
};

Do[
  Print[
    "double C_", i - 1, " = ",
    ToString[
      CForm[
        N[
          bulkCoefs[[i]] /. scaledBulkRules /. solbc,
          16
        ]
      ]
    ],
    ";"
  ],
  {i, 1, Length[bulkCoefs]}
];

(* Knudsen constants with Bessel normalization absorbed *)
knScaledValues =
  Flatten[
    Table[
      {
        C[m, "I", i]*
          BI[0, lambdaList[[m]]*R1/Kn]^-1,

        C[m, "K", i]*
          BK[0, lambdaList[[m]]*R0/Kn]^-1
      },
      {m, 1, Length[lambdaList]},
      {i, 1, nspcL[[m]]}
    ]
  ] /. ruleBessel;

Do[
  Print[
    "double CK_", i, " = ",
    ToString[
      CForm[
        N[
          knScaledValues[[i]] /. solbc,
          16
        ]
      ]
    ],
    ";"
  ],
  {i, 1, Length[knScaledValues]}
];

(* 2. Rename Mathematica constants to temporary valid Mathematica symbols *)
nameRules =
  Join[
   Table[
    bulkCoefs[[i]] -> Symbol["CC" <> ToString[i - 1]],
    {i, 1, Length[bulkCoefs]}
   ],
   Table[
    knCoefs[[i]] -> Symbol["CKK" <> ToString[i]],
    {i, 1, Length[knCoefs]}
   ]
  ];

(* 3. Replace numeric lambda/Kn combinations by symbolic lambda_i/kn *)
lambdaKnRules =
  Flatten[
   Table[
    With[{ii = i},
     {
      x_ /; NumericQ[x] &&
          Abs[N[x] - N[lambdaList[[ii]]/Kn]] < 10^-10
        :> Symbol["lambda" <> ToString[ii]]/Kn,

      x_ /; NumericQ[x] &&
          Abs[N[x] - N[lambdaList[[ii]]]] < 10^-10
        :> Symbol["lambda" <> ToString[ii]]
     }
    ],
    {i, 1, Length[lambdaList]}
   ]
  ];
  
(* Write the r/Kn as symbolic *)

rScaleRule = R -> Kn*R/kn;

(* 4. String replacements after CForm *)
cppStringRules = {
   "CC" -> "C_",
   "CKK" -> "CK_",

   "lambda1" -> "lambda_1",
   "lambda2" -> "lambda_2",
   "lambda3" -> "lambda_3",
   "lambda4" -> "lambda_4",
   "lambda5" -> "lambda_5",
   "lambda6" -> "lambda_6",
   "lambda7" -> "lambda_7",

   "BI(0," -> "I_0(",
   "BI(1," -> "I_1(",
   "BI(2," -> "I_2(",
   "BI(3," -> "I_3(",
   "BI(4," -> "I_4(",
   "BI(5," -> "I_5(",
   "BI(6," -> "I_6(",

   "BK(0," -> "K_0(",
   "BK(1," -> "K_1(",
   "BK(2," -> "K_2(",
   "BK(3," -> "K_3(",
   "BK(4," -> "K_4(",
   "BK(5," -> "K_5(",
   "BK(6," -> "K_6",

   "Log" -> "log",
   "Power" -> "pow",
   "Sqrt" -> "sqrt",
   
   "1.*" -> "",
   "((R)/kn)" -> "(R/kn)",
   "(R)" -> "R"
};

(* 5. Export helper *)
toCpp[expr_] :=
 Module[{e, s},
  e = expr /. nameRules /. lambdaKnRules /. rScaleRule;
  s = ToString[CForm[e]];
  s = StringReplace[s, cppStringRules];
  s = StringReplace[s, {
     "((R)/kn)" -> "(R/kn)",
     "(R)" -> "R"
  }];
  s
 ];
 
 toCppExpl[expr_] :=
 Module[{e, s},
  e = expr /. scaling /. solbc /. ruleBessel;
  s = ToString[CForm[e]];
  s = StringReplace[s, cppStringRules];
  s = StringReplace[s, {
     "((R)/kn)" -> "(R/kn)",
     "(R)" -> "R"
  }];
  s
 ];
 
Print["// Temperature"];
Print["double c_0 = ", toCpp[c0[R]], ";"];
Print["double c = ", toCpp[c[R]], ";"];

Print[""];
Print["// Pressure / density"];
Print["double d_0 = ", toCpp[d0[R]], ";"];
Print["double d = ", toCpp[d[R]], ";"];

Print[""];
Print["// Velocity"];
Print["double a_0 = ", toCpp[a0[R]], ";"];
Print["double a = ", toCpp[a[R]], ";"];
Print["double b_0 = ", toCpp[b0[R]], ";"];
Print["double b = ", toCpp[b[R]], ";"];

Print[""];
Print["// Heat flux"];
Print["double alpha_0 = ", toCpp[alpha0[R]], ";"];
Print["double alpha = ", toCpp[alpha[R]], ";"];
Print["double beta_0 = ", toCpp[beta0[R]], ";"];
Print["double beta = ", toCpp[beta[R]], ";"];

Print[""];
Print["// Scalar R moment"];
Print["double e_0 = ", toCpp[e0[R]], ";"];
Print["double e = ", toCpp[e[R]], ";"];

Print[""];
Print["// Stress"];
Print["double gamma_0 = ", toCpp[gamma0[R]], ";"];
Print["double gamma = ", toCpp[gamma[R]], ";"];
Print["double kappa_0 = ", toCpp[kappa0[R]], ";"];
Print["double kappa = ", toCpp[kappa[R]], ";"];
Print["double omega_0 = ", toCpp[omega0[R]], ";"];
Print["double omega = ", toCpp[omega[R]], ";"];

(* C++ constants/lambdas *)
lambdaLines =
  Table[
    "double lambda_" <> ToString[i] <> " = " <>
      ToString[CForm[N[lambdaList[[i]], 16]]] <> ";",
    {i, 1, Length[lambdaList]}
  ];

scaledBulkCoefs =
  bulkCoefs /. scaledBulkRules;
bulkLines =
  Table[
    "double C_" <> ToString[i - 1] <> " = " <>
      ToString[CForm[N[scaledBulkCoefs[[i]] /. solbc, 16]]] <> ";",
    {i, 1, Length[bulkCoefs]}
  ];

knLines =
  Table[
    "double CK_" <> ToString[i] <> " = " <>
      ToString[CForm[N[knScaledValues[[i]] /. solbc, 16]]] <> ";",
    {i, 1, Length[knCoefs]}
  ];
  
fieldLines = {
  "// Temperature",
  "double c_0 = " <> toCpp[c0[R]] <> ";",
  "double c = " <> toCpp[c[R]] <> ";",
  "",
  "// Pressure",
  "double d_0 = " <> toCpp[d0[R]] <> ";",
  "double d = " <> toCpp[d[R]] <> ";",
  "",
  "// Velocity",
  "double a_0 = " <> toCpp[a0[R]] <> ";",
  "double a = " <> toCpp[a[R]] <> ";",
  "double b_0 = " <> toCpp[b0[R]] <> ";",
  "double b = " <> toCpp[b[R]] <> ";",
  "",
  "// Heat flux",
  "double alpha_0 = " <> toCpp[alpha0[R]] <> ";",
  "double alpha = " <> toCpp[alpha[R]] <> ";",
  "double beta_0 = " <> toCpp[beta0[R]] <> ";",
  "double beta = " <> toCpp[beta[R]] <> ";",
  "",
  "// Scalar R moment",
  "double e_0 = " <> toCpp[e0[R]] <> ";",
  "double e = " <> toCpp[e[R]] <> ";",
  "",
  "// Stress",
  "double gamma_0 = " <> toCpp[gamma0[R]] <> ";",
  "double gamma = " <> toCpp[gamma[R]] <> ";",
  "double kappa_0 = " <> toCpp[kappa0[R]] <> ";",
  "double kappa = " <> toCpp[kappa[R]] <> ";",
  "double omega_0 = " <> toCpp[omega0[R]] <> ";",
  "double omega = " <> toCpp[omega[R]] <> ";",
  
  "// Explicit expressions:",
  "// Temperature",
  "double c_0 = " <> toCppExpl[c0[R]] <> ";",
  "double c = " <> toCppExpl[c[R]] <> ";",
  "",
  "// Pressure",
  "double d_0 = " <> toCppExpl[d0[R]] <> ";",
  "double d = " <> toCppExpl[d[R]] <> ";",
  "",
  "// Velocity",
  "double a_0 = " <> toCppExpl[a0[R]] <> ";",
  "double a = " <> toCppExpl[a[R]] <> ";",
  "double b_0 = " <> toCppExpl[b0[R]] <> ";",
  "double b = " <> toCppExpl[b[R]] <> ";",
  "",
  "// Heat flux",
  "double alpha_0 = " <> toCppExpl[alpha0[R]] <> ";",
  "double alpha = " <> toCppExpl[alpha[R]] <> ";",
  "double beta_0 = " <> toCppExpl[beta0[R]] <> ";",
  "double beta = " <> toCppExpl[beta[R]] <> ";",
  "",
  "// Scalar R moment",
  "double e_0 = " <> toCppExpl[e0[R]] <> ";",
  "double e = " <> toCppExpl[e[R]] <> ";",
  "",
  "// Stress",
  "double gamma_0 = " <> toCppExpl[gamma0[R]] <> ";",
  "double gamma = " <> toCppExpl[gamma[R]] <> ";",
  "double kappa_0 = " <> toCppExpl[kappa0[R]] <> ";",
  "double kappa = " <> toCppExpl[kappa[R]] <> ";",
  "double omega_0 = " <> toCppExpl[omega0[R]] <> ";",
  "double omega = " <> toCppExpl[omega[R]] <> ";"
};

cppLines =
  Join[
    {"// Lambdas"},
    lambdaLines,
    {""},
    {"// Constants"},
    bulkLines,
    knLines,
    {""},
    fieldLines
  ];

Export[
  "C:\\Users\\joshu\\Bachelorarbeit\\r14_exact_export.cpp",
  StringRiffle[cppLines, "\n"],
  "Text"
]



Print["Converting to cartesian coordinates."];

basis[cphi_] = {{Cos[cphi], Sin[cphi]},{-Sin[cphi], Cos[cphi]}};

temperatureC[x_, y_] = temperature[Sqrt[x^2 + y^2], ArcTan[x, y]] /. scaling /. solbc /. ruleBessel;

rhoC[x_, y_] = rho[Sqrt[x^2 + y^2], ArcTan[x, y]] /. scaling /. solbc /. ruleBessel;

{ux[x_, y_], uy[x_, y_]} = {ur[Sqrt[x^2 + y^2], ArcTan[x, y]], 
                            uTheta[Sqrt[x^2 + y^2], ArcTan[x, y]]} . basis[ArcTan[x, y]] /. scaling /. solbc /. ruleBessel;

plotLegends = Placed[
   SwatchLegend[{Blue, Red, Green}, 
            (Style[#1, FontSize -> 12] & ) /@ {"\[Theta]=0", "\[Theta]=\[Pi]/4", "\[Theta]=\[Pi]/2"}], Scaled[{0.85, 0.16}]]

PlotExactTemperature = Plot[
    {temperatureC[r, 0], temperatureC[r*Cos[Pi/4], r*Sin[Pi/4]], temperatureC[r*Cos[Pi/2], r*Sin[Pi/2]]}, {r, R0, R1}, 
    Frame -> True, FrameLabel -> {"radius r", "temperature"}, PlotStyle -> {{Blue, style}, {Red, style}, {Green, style}},
    PlotLegends -> plotLegends ]; 

PlotExactVelocity = Plot[
    {Norm[{ux[r, 0], uy[r, 0]}], Norm[{ux[r*Cos[Pi/4], r*Sin[Pi/4]], uy[r*Cos[Pi/4], r*Sin[Pi/4]]}], 
     Norm[{ux[r*Cos[Pi/2], r*Sin[Pi/2]], uy[r*Cos[Pi/2], r*Sin[Pi/2]]}]}, {r, R0, R1}, 
     Frame -> True, FrameLabel -> {"radius r", "speed"}, PlotStyle -> {{Blue, style}, {Red, style}, {Green, style}},
     PlotLegends -> plotLegends ]; 


(* if necessary the solution can be saved as 3D surface plot here *)
Get["AnalyticalPlot3D.wls"]; 

