#!/usr/bin/env wolframscript
(* ::Package:: *)

(* this requires a specific problem setup for R13 *)


Print["*******"];
Print["Analytical solution of " <> theory <> " equations."];

Get[theory <> "AnalyticalSolution.wls"];


Print["Solving for integration constants."];

Kn = \[Tau];


solbc = Solve[BCeqns == 0, coefs][[1]];


Print["Converting to cartesian coordinates."];

basis[cphi_] = {{Cos[cphi], Sin[cphi]},{-Sin[cphi], Cos[cphi]}};

temperatureC[x_, y_] = temperature[Sqrt[x^2 + y^2], ArcTan[x, y]] /. scaling /. solbc /. ruleBessel;
rhoC[x_, y_] = rho[Sqrt[x^2 + y^2], ArcTan[x, y]] /. scaling /. solbc /. ruleBessel;

{ux[x_, y_], uy[x_, y_]} = {ur[Sqrt[x^2 + y^2], ArcTan[x, y]], 
                            uTheta[Sqrt[x^2 + y^2], ArcTan[x, y]]}.basis[ArcTan[x, y]] /. scaling /. solbc /. ruleBessel;

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
