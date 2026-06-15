#!/usr/bin/env wolframscript
(* ::Package:: *)

(*  *)


(* R14 boundary conditions *)
BCeqn0 = { 
   (-epsp)*(vn - vnW) + epsv*(rho + temperature + facBC*Sigmann + pW), 
   -Sigmant + chi2*((vt - vtW) + qt/5 + mnnt), 
   -qn + chi2*(2*(temperature - temperatureW) + Sigmann/2 + (2*Rnn)/5 + RR/15 ), 
   -mnnn    +    chi2*((-(2/5))*(temperature - temperatureW) + (7*Sigmann)/5 - (2*Rnn)/25 - RR/75), 
   -mntt - (1/2)*chi2*((-(2/5))*(temperature - temperatureW) + (7*Sigmann)/5 - (2*Rnn)/25 - RR/75) + chi2*(Sigmann/2 + Sigmatt), 
   -Rnt + chi2*(-(vt - vtW) + (11*qt)/5 - mnnt),
   -Omegan + chi2*(-8*(temperature - temperatureW) -2*Sigmann - (8*Rnn)/5 + (4*RR)/3 )
   }; 

Print["Definition of R14 boundary conditions."];

{BCrhs, BC} = Normal[CoefficientArrays[
    BCeqn0, {rho, vn, vt, temperature, Sigmann, Sigmant, Sigmatt, 
     qn, qt, RR, Rnn, Rnt, Rtt, mnnn, mnnt, mntt, mttt, Omegan, Omegat}]];



