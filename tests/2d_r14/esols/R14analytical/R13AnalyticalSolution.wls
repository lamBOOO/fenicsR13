#!/usr/bin/env wolframscript
(* ::Package:: *)

Get["DifferentialOperators.wls"];

Print["Constructing analytical solution for R13 between two circles of radii (R0,R1)."];
Print["(radial symmetric part and one angular Fourier mode)"];

coef[expr_, var_, idx_] := Flatten[Normal[CoefficientArrays[Extract[expr, idx], 
      var]]]; 
TraceFreeSym2[mat_] := Block[{perm, sym1, trace}, perm = Permutations[{1, 2}]; 
    sym1 = Expand[(1/Length[perm])*Sum[Transpose[mat, idx], {idx, perm}]]; 
    trace = Sum[sym1[[jj,jj]], {jj, 1, 3}]; 
    Expand[sym1 - (1/3)*trace*Table[KroneckerDelta[ii, jj], {ii, 1, 3}, {jj, 1, 3}]]]
TraceFreeSym3[mat_] := Block[{perm, sym1, trace}, perm = Permutations[{1, 2, 3}]; 
    sym1 = Expand[(1/Length[perm])*Sum[Transpose[mat, idx], {idx, perm}]]; 
    trace = Sum[sym1[[All,jj,jj]], {jj, 1, 3}]; 
    Expand[sym1 - (1/5)*Table[trace[[ii]]*KroneckerDelta[jj, kk] + 
         trace[[jj]]*KroneckerDelta[ii, kk] + trace[[kk]]*KroneckerDelta[ii, jj], 
        {ii, 1, 3}, {jj, 1, 3}, {kk, 1, 3}]]]
scalarIdx = {{}}; 
vectorIdx = {{1}, {2}}; 
tensor2Idx = {{1, 1}, {1, 2}, {2, 2}}; 
tensor3Idx = {{1, 1, 1}, {1, 1, 2}, {1, 2, 2}, {2, 2, 2}}; 

Clear[Kn, a, b, alpha, beta, c, d, gamma, kappa, omega, a0, b0, c0, d0, alpha0, beta0, gamma0, omega0, kappa0, 
  BK, BI, lambda, lambdaList]
nVar = 18; 
ansatz = {a, b, alpha, beta, c, d, gamma, kappa, omega, a0, b0, c0, d0, alpha0, beta0, gamma0, omega0, kappa0}; 
temperature[r_, theta_] = c0[r] + c[r]*Cos[theta]; 
sr[r_, theta_] = alpha0[r] + Cos[theta]*alpha[r]; 
sTheta[r_, theta_] = beta0[r] - Sin[theta]*beta[r]; 
s[r_, theta_] = {sr[r, theta], sTheta[r, theta], 0}; 
ur[r_, theta_] = a0[r] + Cos[theta]*a[r]; 
uTheta[r_, theta_] = b0[r] - Sin[theta]*b[r]; 
u[r_, theta_] = {ur[r, theta], uTheta[r, theta], 0}; 
p[r_, theta_] = d0[r] + d[r]*Cos[theta]; 
rho[r_, theta_] = p[r, theta] - temperature[r, theta]; 
Sigmarr[r_, theta_] = Cos[theta]*gamma[r] + gamma0[r]; 
SigmaThetaTheta[r_, theta_] = -(Cos[theta]*omega[r] + omega0[r]); 
SigmarTheta[r_, theta_] = kappa0[r] + Sin[theta]*kappa[r]; 
Sigmazz[r_, theta_] = -Sigmarr[r, theta] - SigmaThetaTheta[r, theta]; 
Sigma[r_, theta_] = {{Sigmarr[r, theta], SigmarTheta[r, theta], 0}, {SigmarTheta[r, theta], SigmaThetaTheta[r, theta], 0}, 
    {0, 0, Sigmazz[r, theta]}}; 

(* equations and BCs *)
Get["R13Equations.wls"];
Get["R13BoundaryConditions.wls"];

eqns = {coef[eqns0[[1]], {Cos[theta], Sin[theta]}, scalarIdx], 
   coef[eqns0[[2]], {Cos[theta], Sin[theta]}, vectorIdx],
   coef[eqns0[[3]], {Cos[theta], Sin[theta]}, scalarIdx ], 
   coef[eqns0[[4]], {Cos[theta], Sin[theta]}, tensor2Idx],
   coef[eqns0[[5]], {Cos[theta], Sin[theta]}, vectorIdx]
   }; 

{Rrr[r_, theta_], RrTheta[r_, theta_], RThetaTheta[r_, theta_]} = Extract[R[r, theta], tensor2Idx]; 
{mrrr[r_, theta_], mrrTheta[r_, theta_], mrThetaTheta[r_, theta_], mThetaThetaTheta[r_, theta_]} = 
   Extract[m[r, theta], tensor3Idx]; 
   
Print["Computing eigenvalues of solution modes."]
Clear[lambda]
eqnsL = eqns /. Power[r, n_?Negative] -> 0; 
{a[r_], b[r_], alpha[r_], beta[r_], c[r_], d[r_], gamma[r_], kappa[r_], omega[r_], a0[r_], b0[r_], 
    c0[r_], d0[r_], alpha0[r_], beta0[r_], gamma0[r_], omega0[r_], kappa0[r_]} = 
   Table[v[j]*Exp[lambda*r], {j, 1, nVar}]; 
coefs = Table[v[j], {j, 1, nVar}]; 
eqnsL1 = Select[Flatten[Expand[Exp[(-lambda)*r]*eqnsL /. Kn -> 1]], Length[#1] > 0 & ]; 
vars2 = Sort[DeleteDuplicates[(Extract[eqnsL1, Append[#1, -1]] & ) /@ 
      Position[eqnsL1, (val_)*lambda^2*v[_]]]]; 
eqnsL2 = Join[eqnsL1 /. Thread[Expand[lambda^2*vars2] -> lambda*w /@ Range[Length[vars2]]], 
    w /@ Range[Length[vars2]] - Expand[lambda*vars2]]; 
lambdaList = Select[DeleteDuplicates[
     lambda /. Solve[Det[CoefficientArrays[eqnsL2, Join[coefs, 
            w /@ Range[Length[vars2]]]][[2]]] == 0, lambda]], N[#1] > 0 & ]; 
            
Print[lambdaList];

Print["Computing bulk solution."]
kmax = 5; 
phi = Flatten[Join[{1, Log[r]}, Table[{r^k, r^(-k)}, {k, 1, kmax}]]]; 
{a[r_], b[r_], alpha[r_], beta[r_], c[r_], d[r_], gamma[r_], kappa[r_], omega[r_], a0[r_], b0[r_], 
    c0[r_], d0[r_], alpha0[r_], beta0[r_], gamma0[r_], omega0[r_], kappa0[r_]} = 
   Table[Table[v[j, k], {k, 1, Length[phi]}] . phi, {j, 1, nVar}]; 
coefs = Flatten[Table[Table[v[j, k], {k, 1, Length[phi]}], {j, 1, nVar}]]; 
{bb, AA} = Normal[CoefficientArrays[Expand[Flatten[eqns /. {Kn -> 1}]] /. 
        {Log[r] -> Logr[0]} /. {r^(j_.) -> rr[j]} /. {Logr[0]*rr[k_] -> Logr[k]}, 
     Join[Table[rr[k], {k, -2*kmax, 2*kmax}], Table[Logr[k], 
       {k, -2*kmax, 2*kmax}]]]]; 
{bb, AA} = Normal[CoefficientArrays[Select[Flatten[{bb, AA}], Length[#1] > 0 & ], 
     coefs]]; 
nspc = NullSpace[AA]; 
nC = Length[nspc]; 
sol0 = Sum[C[i]*Table[Table[v[j, k], {k, 1, Length[phi]}] . phi, {j, 1, nVar}] /. 
     Thread[coefs -> nspc[[i]]], {i, 1, nC}]; 
     
(* construct individual rules for Bessel functions *)
Clear[lam$, qq]
BesselRule0["I"] = {BI[-1, z_] -> BI[1, z], BI[-2, z_] -> BI[2, z]}; 
BesselRule1a["I"] = {Derivative[0, 1][BI][i_, (lam$_)*r] -> 
     (1/2)*(BI[i - 1, lam$*r] + BI[i + 1, lam$*r])}; 
BesselRule1b["I"] = {Derivative[0, 2][BI][i_, (lam$_)*r] -> 
     (1/2)*(Derivative[0, 1][BI][i - 1, lam$*r] + Derivative[0, 1][BI][i + 1, lam$*r])}; 
BesselRule2["I"] = {BI[(i_)?(#1 > 0 & ), (lam$_)*r]*r^(qq_.) -> 
     (lam$/(2*i))*(BI[i - 1, lam$*r] - BI[i + 1, lam$*r])*r^(qq + 1)}; 
BesselRule3["I"] = {BI[0, (lam$_)*r]*r^(qq_.) -> BI[0, lam$*r^(1 + Abs[qq])]}; 
BesselRule0["K"] = {BK[-1, z_] -> BK[1, z], BK[-2, z_] -> BK[2, z]}; 
BesselRule1a["K"] = {Derivative[0, 1][BK][i_, (lam$_)*r] -> 
     (-(1/2))*(BK[i - 1, lam$*r] + BK[i + 1, lam$*r])}; 
BesselRule1b["K"] = {Derivative[0, 2][BK][i_, (lam$_)*r] -> 
     (-(1/2))*(Derivative[0, 1][BK][i - 1, lam$*r] + Derivative[0, 1][BK][i + 1, lam$*r])}; 
BesselRule2["K"] = {BK[(i_)?(#1 > 0 & ), (lam$_)*r]*r^(qq_.) -> 
     (-(lam$/(2*i)))*(BK[i - 1, lam$*r] - BK[i + 1, lam$*r])*r^(qq + 1)}; 
BesselRule3["K"] = {BK[0, (lam$_)*r]*r^(qq_.) -> BK[0, lam$*r^(1 + Abs[qq])]}; 
ruleBessel = {BK -> BesselK, BI -> BesselI}; 

Print["Computing Knudsen layers."];
Clear[l, case]
kmax = 7; 
Do[nspcL = {}; 
   Do[lambda = lambdaList[[l]]; Bessel[k_, z_] := 
      ToExpression[StringJoin["B", case]][k, z]; 
     {a[r_], b[r_], alpha[r_], beta[r_], c[r_], d[r_], gamma[r_], kappa[r_], omega[r_], a0[r_], 
       b0[r_], c0[r_], d0[r_], alpha0[r_], beta0[r_], gamma0[r_], omega0[r_], kappa0[r_]} = 
      Table[Sum[v[l, j, k]*Bessel[k, lambda*r], {k, 0, kmax}], {j, 1, nVar}]; 
     coefs = Flatten[Table[Table[v[l, j, k], {k, 0, kmax}], {j, 1, nVar}]]; 
     eqns3 = Expand[Expand[Expand[Expand[Flatten[eqns /. {Kn -> 1}]] /. BesselRule1b[
                case] /. BesselRule1a[case] /. BesselRule0[case]] /. 
           BesselRule2[case]] /. BesselRule2[case]] /. BesselRule3[case]; 
     {bb, AA} = Normal[CoefficientArrays[
        Select[Flatten[Normal[CoefficientArrays[eqns3, DeleteDuplicates[
              Cases[eqns3, Bessel[_, _], 3]]]][[2]]], Length[#1] > 0 & ], coefs]]; 
     nspc = NullSpace[AA]; nspcL = Append[nspcL, Length[nspc]]; 
     Do[sol[l, case, i] = Chop[Table[Sum[v[l, j, k]*Bessel[k, lambda*r], {k, 0, kmax}], 
          {j, 1, nVar}] /. Thread[coefs -> nspc[[i]]]], {i, 1, Length[nspc]}]; 
   , {l, 1, Length[lambdaList]}]; 
, {case, {"I", "K"}}]
    
{a[r_], b[r_], alpha[r_], beta[r_], c[r_], 
   d[r_], gamma[r_], kappa[r_], omega[r_], a0[r_], b0[r_], 
   c0[r_], d0[r_], alpha0[r_], beta0[r_], gamma0[
    r_], omega0[r_], kappa0[r_]} = 
  sol0 + Sum[
     Sum[C[m, "I", i] sol[m, "I", i] + 
       C[m, "K", i] sol[m, "K", i], {i, 1, nspcL[[m]]}], {m, 1, 
      Length[lambdaList]}] /. {r -> r/Kn};
      
Print["Validating analytical solution:"];
Print[Table[Chop[Expand[eqns /. {Kn -> 0.1, r -> rr}] /. ruleBessel, 10^-7], {rr, {0.5, 0.3, 0.25, 1.0}}]];

coefs = Join[C /@ Range[nC], 
  Flatten[Table[
    Table[{C[m, "I", i], C[m, "K", i]}, {i, 1, nspcL[[m]]}], {m, 1, 
     Length[lambdaList]}]]];
           
U[r_, theta_] := {rho[r, theta], ur[r, theta], uTheta[r, theta], 
    temperature[r, theta], Sigmarr[r, theta], SigmarTheta[r, theta], SigmaThetaTheta[r, theta], 
    sr[r, theta], sTheta[r, theta], 
    Rrr[r, theta], RrTheta[r, theta], 
    RThetaTheta[r, theta], mrrr[r, theta], 
    mrrTheta[r, theta], mrThetaTheta[r, theta], 
    mThetaThetaTheta[r, theta]} /. ruleBessel;
          
(* including 1D transformation into "general" coordinates *)
Un[r_, theta_] := 
  DiagonalMatrix[{1, nn[r], tt[r], 1, 1, tt[r] nn[r], 1, nn[r], tt[r], 
      1, tt[r] nn[r], 1, nn[r], tt[r], nn[r], tt[r]}] . U[r, theta];

scaling := 
  Join[{C[1] -> C[1] Kn^-1, C[4] -> C[4] Kn^2, 
    C[9] -> C[9] Kn^2, C[11] -> C[11] Kn^2}, 
   Flatten[Table[{C[m, "I", i] -> 
        C[m, "I", i] BI[0, lambdaList[[m]] R1/Kn]^-1, 
       C[m, "K", i] -> 
        C[m, "K", i] BK[0, lambdaList[[m]] R0/Kn]^-1} /. 
      ruleBessel, {m, 1, Length[lambdaList]}, {i, 1, nspcL[[m]]}]]];

BCeqns := (Join[
     coef[
      Expand[BC . Un[R0, theta] + BCrhs /. 
        setup[R0, theta]], {Cos[theta], Sin[theta]}, {{}}],
     coef[
      Expand[BC . Un[R1, theta] + BCrhs /. 
        setup[R1, theta]], {Cos[theta], Sin[theta]}, {{}}]
     ] /. scaling);

