#!/usr/bin/env wolframscript
(* ::Package:: *)

(* MFS-specifics: we remove delta in the BC below and also in the BC-variables of CoefficientArrays[] and the diag-Trafo,
 we also switch from p to rho+temperature in BCs
 ultimately we change the order of the variables in BC and U *)


(* R13 boundary conditions *)
BCeqn0 = { 
   (-epsp)*(vn - vnW) + epsv*(rho + temperature + facBC*Sigmann + pW), 
   -Sigmant + chi2*((vt - vtW) + qt/5 + mnnt), 
   -qn + chi2*(2*(temperature - temperatureW) + Sigmann/2 + (2*Rnn)/5 ), 
   -mnnn + chi2*((-(2/5))*(temperature - temperatureW) + (7*Sigmann)/5 - (2*Rnn)/25 ), 
   -((1/2)*chi2*((-(2/5))*(temperature - temperatureW) + (7*Sigmann)/5 - (2*Rnn)/25) + mntt) + chi2*(Sigmann/2 + Sigmatt), 
   -Rnt + chi2*(-(vt - vtW) + (11*qt)/5 - mnnt)
   }; 

Print["Definition of R13 boundary conditions."];

{BCrhs, BC} = Normal[CoefficientArrays[
    BCeqn0, {rho, vn, vt, temperature, Sigmann, Sigmant, Sigmatt, 
     qn, qt, Rnn, Rnt, Rtt, mnnn, mnnt, mntt, mttt}]];


(* variable transformation into general coordinate system 
Trafo = {
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
    {0, nx, ny, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
    {0, tx, ty, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 0, nx^2, 2*nx*ny, ny^2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 0, nx*tx, ny*tx + nx*ty, ny*ty, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 0, tx^2, 2*tx*ty, ty^2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 0, 0, nx, ny, 0, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 0, 0, tx, ty, 0, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, nx^2, 2*nx*ny, ny^2, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, nx*tx, ny*tx + nx*ty, ny*ty, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, tx^2, 2*tx*ty, ty^2, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, nx^3, 3*nx^2*ny, 3*nx*ny^2, ny^3}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, nx^2*tx, 2*nx*ny*tx + nx^2*ty, ny^2*tx + 2*nx*ny*ty, ny^2*ty}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, nx*tx^2, ny*tx^2 + 2*nx*tx*ty, 2*ny*tx*ty + nx*ty^2, ny*ty^2}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, tx^3, 3*tx^2*ty, 3*tx*ty^2, ty^3}
    }; 
*)
