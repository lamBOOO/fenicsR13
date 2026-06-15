#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Set the working directory to the location of this script/notebook *)
Check[
  SetDirectory[
    DirectoryName[
      If[$InputFileName =!= "", $InputFileName, NotebookFileName[]]
    ]
  ],
  Print["Could not set directory. Ensure the file is saved."]
];


Print["Hello World!"];

ClearAll["Global`*"];


Print["Setup for a specific solution."];

R0 = 1; 
R1 = 2; 
\[Tau] = 0.5; 
facBC = 1; 
chi2value = 1; 
v0 = 2.; 
vR = -1.; 
p0 = 0.; 
eps = 10^(-5); 
temperature0 = 1.; 
temperature1 = 2.; 

chi2BC[r_, theta_] = Piecewise[{
                            {chi2value, r == R0}, {chi2value, r == R1}}]; 
temperatureBC[r_, theta_] = Piecewise[{
                            {temperature0, r == R0}, {temperature1, r == R1}}]; 
vtBC[r_, theta_] = Piecewise[{
                            {vR, r == R0}, {(-v0)*Sin[theta], r == R1}}]; 
vnBC[r_, theta_] = Piecewise[{
                            {0, r == R0}, {v0*Cos[theta], r == R1}}]; 
pBC[r_, theta_] = Piecewise[{
                            {0, r == R0}, {p0*Cos[theta], r == R1}}]; 
epspBC[r_, theta_] = Piecewise[{
                            {1, r == R0}, {1, r == R1}}]; 
epsvBC[r_, theta_] = Piecewise[{
                            {eps, r == R0}, {1, r == R1}}]; 

nn[R0] = -1; 
nn[R1] = 1; 
tt[R0] = -1; 
tt[R1] = 1; 
setup[r_, theta_] := {
  chi2 -> chi2BC[r, theta],
  temperatureW -> temperatureBC[r, theta], 
  vnW -> vnBC[r, theta], 
  vtW -> vtBC[r, theta], 
  pW -> pBC[r, theta], 
  epsp -> epspBC[r, theta], 
  epsv -> epsvBC[r, theta]}

(* default initialization *)
PlotR13Temperature = Nothing;
PlotR13Velocity = Nothing;
PlotR14Temperature = Nothing;
PlotR14Velocity = Nothing;


(* analytical solution R14 *)
style = Dashing[0.01];
theory = theory1 = "R13";
Get["AnalyticalExample.wls"];  
PlotR13Temperature = PlotExactTemperature;
PlotR13Velocity = PlotExactVelocity;


(* analytical solution R14 *)
style = Dashing[{}];
theory = theory2 = "R14";
Get["AnalyticalExample.wls"]; 
PlotR14Temperature = PlotExactTemperature;
PlotR14Velocity = PlotExactVelocity;

Print["Saving a comparison as plot."];
GraphicsRow[{
   Show[{PlotR14Temperature, PlotR13Temperature}, PlotLabel -> theory1 <> " and " <> theory2 <> " analytical solutions", 
    ImageSize -> 900], 
   Show[{PlotR14Velocity, PlotR13Velocity},  PlotLabel -> theory1 <> " and " <> theory2 <> " analytical solutions", 
   ImageSize -> 900]}]
output=Export["speed_temperature_plots.png",
GraphicsRow[{
   Show[{PlotR14Temperature, PlotR13Temperature}, PlotLabel -> theory1 <> " and " <> theory2 <> " analytical solutions", 
    ImageSize -> 900], 
   Show[{PlotR14Velocity, PlotR13Velocity},  PlotLabel -> theory1 <> " and " <> theory2 <> " analytical solutions", 
   ImageSize -> 900]}]];

Print["Done: "<>output];




