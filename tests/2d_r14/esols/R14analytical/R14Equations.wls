#!/usr/bin/env wolframscript
(* ::Package:: *)

Print["Definition of 2-dimensional R14 equations and boundary conditions."];

(* we add delta[]=scriptR[] here and the diag-Trafo of the variable vector,
 Note: ansatz necessary for delta[]=scriptR[]
 be careful with the order of the variables in BC and U *)


(* R14 closure *)
m[r_, theta_] = -2*Kn*TraceFreeSym3[Grad2[Sigma[r, theta]]]; 
R[r_, theta_] = (-(24/5))*Kn*TraceFreeSym2[Grad1[s[r, theta]]];
Omega[r_, theta_] = (-(7/3))*Kn*Grad0[scriptR[r, theta]]; 

(* R14 equations *)
eqns0 = {
   Div1[u[r, theta]],
   Grad0[p[r, theta]] + Div2[Sigma[r, theta]],
   Div1[s[r, theta]] + Div1[u[r, theta]],
   2*TraceFreeSym2[Grad1[u[r, theta]]] + (4/5)*TraceFreeSym2[Grad1[s[r, theta]]] 
      + (1/Kn)*Sigma[r, theta] + Div3[m[r, theta]],
   (5/2)*Grad0[temperature[r, theta]] + Div2[Sigma[r, theta]] + (1/2)*Div2[R[r, theta]] 
      + (1/6)*Grad0[scriptR[r, theta]] + (2/(3*Kn))*s[r, theta],
    8*Div1[s[r, theta]] + Div1[Omega[r, theta]] + (2/3*Kn)*scriptR[r, theta]  
   };
