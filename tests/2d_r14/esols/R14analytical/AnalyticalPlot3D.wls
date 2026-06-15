#!/usr/bin/env wolframscript
(* ::Package:: *)

(* this requires a specific problem setup and analytical solution for Rxy *)
(* 3D plot is saved to AnalyticalPlot3D.png *)


Print["Saving the analytical solution as 3D surface plot."];

output=Export[theory <> "AnalyticalPlot3D.png",
 GraphicsRow[{
 Plot3D[temperatureC[x, y], {x, -2.1, 2.1}, {y, -2.1, 2.1}, 
 PlotPoints -> 30, 
 ColorFunction -> Function[{x, y, z}, ColorData["SouthwestColors"][z]], 
 RegionFunction -> Function[{x, y, z}, R0 < Sqrt[x^2 + y^2] < R1], 
 PlotRange -> All,PlotLabel->theory <> " temperature"],
 Plot3D[ux[x, y], {x, -2.1, 2.1}, {y, -2.1, 2.1}, 
 PlotPoints -> 30, 
 ColorFunction -> Function[{x, y, z}, ColorData["SouthwestColors"][z]], 
 RegionFunction -> Function[{x, y, z}, R0 < Sqrt[x^2 + y^2] < R1], 
 PlotRange -> All,PlotLabel->theory <> " x-velocity"],
 Plot3D[uy[x, y], {x, -2.1, 2.1}, {y, -2.1, 2.1}, 
 PlotPoints -> 30, 
 ColorFunction -> Function[{x, y, z}, ColorData["SouthwestColors"][z]], 
 RegionFunction -> Function[{x, y, z}, R0 < Sqrt[x^2 + y^2] < R1], 
 PlotRange -> All,PlotLabel->theory <> " y-velocity"]}]];

Print["Done: "<>output];
