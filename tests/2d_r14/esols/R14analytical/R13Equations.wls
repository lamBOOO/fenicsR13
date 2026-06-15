#!/usr/bin/env wolframscript
(* ::Package:: *)

Print["Definition of 2-dimensional R13 equations and boundary conditions."];

(* MFS-specifics: we remove delta[] here and the diag-Trafo of the variable vector,
 Note: no ansatz necessary for delta[]
 ultimately we change the order of the variables in BC and U *)


(* R13 closure *)
m[r_, theta_] = -2*Kn*TraceFreeSym3[Grad2[Sigma[r, theta]]]; 
R[r_, theta_] = (-(24/5))*Kn*TraceFreeSym2[Grad1[s[r, theta]]]; 

(* R13 equations *)
eqns0 = {
   Div1[u[r, theta]],
   Grad0[p[r, theta]] + Div2[Sigma[r, theta]],
   Div1[s[r, theta]] + Div1[u[r, theta]],
   2*TraceFreeSym2[Grad1[u[r, theta]]] + (4/5)*TraceFreeSym2[Grad1[s[r, theta]]] 
      + (1/Kn)*Sigma[r, theta] + Div3[m[r, theta]],
   (5/2)*Grad0[temperature[r, theta]] + Div2[Sigma[r, theta]] + (1/2)*Div2[R[r, theta]] 
      + (2/(3*Kn))*s[r, theta]
   };
