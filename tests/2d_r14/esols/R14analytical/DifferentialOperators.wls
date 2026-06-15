#!/usr/bin/env wolframscript
(* ::Package:: *)

Print["Defining differential operators for tensors in cylindrical coordinates."]

$Assumptions = {r > 0, 0 < theta < Pi, Element[{theta, r, z}, Reals]}; 
xx = {r, theta, z}; 
P = {r*Cos[theta], r*Sin[theta], z}; 
J = Simplify[Table[D[P[[nu]], xx[[mu]]], {nu, 1, 3}, {mu, 1, 3}]]; 
Ji = Simplify[Inverse[J]]; 
(* metric coefficients *)
gg = Simplify[Transpose[J] . J]; 
gi = Simplify[Inverse[gg]]; 
(* Christoffel symbols*)
Gam = Simplify[Table[(1/2)*Sum[gi[[lambda,sigma]]*(D[gg[[sigma,nu]], xx[[mu]]] + 
         D[gg[[mu,sigma]], xx[[nu]]] - D[gg[[mu,nu]], xx[[sigma]]]), 
       {sigma, 1, 3}], {lambda, 1, 3}, {mu, 1, 3}, {nu, 1, 3}]]; 
       
(* gradient, Laplace of scalars *)
Grad0[S_] := (temp = Simplify[Table[D[S, xx[[i]]], {i, 1, 3}]]; 
   Simplify[Table[Sqrt[gi[[i,i]]]*temp[[i]], {i, 1, 3}]])
   
Laplace0[S_] := (temp1 = Simplify[Table[D[S, xx[[i]]], {i, 1, 3}]]; 
   temp2 = Simplify[Table[D[temp1[[i]], xx[[k]]] - Sum[Gam[[sigma,i,k]]*temp1[[sigma]], 
        {sigma, 1, 3}], {i, 1, 3}, {k, 1, 3}]]; 
   Simplify[Sum[gi[[k,k]]*temp2[[k,k]], {k, 1, 3}]])
   
(* gradient, divergence, Laplace of vectors *)
Grad1[V_] := (temp = Simplify[Table[V[[i]]/Sqrt[gg[[i,i]]], {i, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i]], xx[[k]]] + Sum[Gam[[i,sigma,k]]*temp[[sigma]], 
        {sigma, 1, 3}], {i, 1, 3}, {k, 1, 3}]]; 
   Simplify[Table[Sqrt[gg[[i,i]]]*Sqrt[gi[[k,k]]]*temp1[[i,k]], {i, 1, 3}, 
     {k, 1, 3}]])

Div1[V_] := (temp = Simplify[Table[V[[i]]/Sqrt[gg[[i,i]]], {i, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i]], xx[[k]]] + Sum[Gam[[i,sigma,k]]*temp[[sigma]], 
        {sigma, 1, 3}], {i, 1, 3}, {k, 1, 3}]]; Simplify[Sum[temp1[[k,k]], {k, 1, 3}]])

Laplace1[V_] := (temp = Simplify[Table[V[[i]]/Sqrt[gg[[i,i]]], {i, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i]], xx[[k]]] + Sum[Gam[[i,sigma,k]]*temp[[sigma]], 
        {sigma, 1, 3}], {i, 1, 3}, {k, 1, 3}]]; 
   temp2 = Simplify[Table[D[temp1[[i,j]], xx[[k]]] + Sum[Gam[[i,sigma,k]]*temp1[[sigma,j]], 
        {sigma, 1, 3}] - Sum[Gam[[sigma,j,k]]*temp1[[i,sigma]], {sigma, 1, 3}], {i, 1, 3}, 
      {j, 1, 3}, {k, 1, 3}]]; 
   Simplify[Table[Sqrt[gg[[i,i]]]*Sum[gi[[k,k]]*temp2[[i,k,k]], {k, 1, 3}], 
     {i, 1, 3}]])
     
(* gradient, divergence, Laplace of matrices *)
Div2[V_] := (temp = Simplify[Table[V[[i,j]]/Sqrt[gg[[i,i]]*gg[[j,j]]], {i, 1, 3}, 
      {j, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i,j]], xx[[k]]] + 
       Sum[Gam[[i,sigma,k]]*temp[[sigma,j]], {sigma, 1, 3}] + Sum[Gam[[j,sigma,k]]*temp[[i,sigma]], 
        {sigma, 1, 3}], {i, 1, 3}, {j, 1, 3}, {k, 1, 3}]]; 
   Simplify[Table[Sqrt[gg[[i,i]]]*Sum[temp1[[i,k,k]], {k, 1, 3}], {i, 1, 3}]])

Grad2[V_] := (temp = Simplify[Table[V[[i,j]]/Sqrt[gg[[i,i]]*gg[[j,j]]], {i, 1, 3}, 
      {j, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i,j]], xx[[k]]] + 
       Sum[Gam[[i,sigma,k]]*temp[[sigma,j]], {sigma, 1, 3}] + Sum[Gam[[j,sigma,k]]*temp[[i,sigma]], 
        {sigma, 1, 3}], {i, 1, 3}, {j, 1, 3}, {k, 1, 3}]]; 
   Simplify[Table[Sqrt[gg[[i,i]]]*Sqrt[gg[[j,j]]]*Sqrt[gi[[k,k]]]*temp1[[i,j,k]], 
     {i, 1, 3}, {j, 1, 3}, {k, 1, 3}]])

Laplace2[V_] := (temp = Simplify[Table[V[[i,j]]/Sqrt[gg[[i,i]]*gg[[j,j]]], 
      {i, 1, 3}, {j, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i,j]], xx[[k]]] + Sum[Gam[[i,sigma,k]]*temp[[sigma,j]], 
        {sigma, 1, 3}] + Sum[Gam[[j,sigma,k]]*temp[[i,sigma]], {sigma, 1, 3}], {i, 1, 3}, 
      {j, 1, 3}, {k, 1, 3}]]; 
   temp2 = Simplify[Table[D[temp1[[i,j,k]], xx[[l]]] + 
       Sum[Gam[[i,sigma,l]]*temp1[[sigma,j,k]], {sigma, 1, 3}] + 
       Sum[Gam[[j,sigma,l]]*temp1[[i,sigma,k]], {sigma, 1, 3}] - 
       Sum[Gam[[sigma,k,l]]*temp1[[i,j,sigma]], {sigma, 1, 3}], {i, 1, 3}, {j, 1, 3}, 
      {k, 1, 3}, {l, 1, 3}]]; 
   Simplify[Table[Sqrt[gg[[i,i]]]*Sqrt[gg[[j,j]]]*
      Sum[gi[[k,k]]*temp2[[i,j,k,k]], {k, 1, 3}], {i, 1, 3}, {j, 1, 3}]])

(* gradient, divergence, Laplace of 3-tensors *) 
Div3[V_] := (temp = Simplify[Table[V[[i,j,k]]/Sqrt[gg[[i,i]]*gg[[j,j]]*gg[[k,k]]], 
      {i, 1, 3}, {j, 1, 3}, {k, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i,j,k]], xx[[l]]] + 
       Sum[Gam[[i,sigma,l]]*temp[[sigma,j,k]], {sigma, 1, 3}] + 
       Sum[Gam[[j,sigma,l]]*temp[[i,sigma,k]], {sigma, 1, 3}] + 
       Sum[Gam[[k,sigma,l]]*temp[[i,j,sigma]], {sigma, 1, 3}], {i, 1, 3}, {j, 1, 3}, 
      {k, 1, 3}, {l, 1, 3}]]; Simplify[Table[Sqrt[gg[[i,i]]*gg[[j,j]]]*
      Sum[temp1[[i,j,k,k]], {k, 1, 3}], {i, 1, 3}, {j, 1, 3}]])
      
Grad3[V_] := (temp = Simplify[Table[V[[i,j,k]]/Sqrt[gg[[i,i]]*gg[[j,j]]*gg[[k,k]]], 
      {i, 1, 3}, {j, 1, 3}, {k, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i,j,k]], xx[[l]]] + 
       Sum[Gam[[i,sigma,l]]*temp[[sigma,j,k]], {sigma, 1, 3}] + 
       Sum[Gam[[j,sigma,l]]*temp[[i,sigma,k]], {sigma, 1, 3}] + 
       Sum[Gam[[k,sigma,l]]*temp[[i,j,sigma]], {sigma, 1, 3}], {i, 1, 3}, {j, 1, 3}, 
      {k, 1, 3}, {l, 1, 3}]]; Simplify[Table[Sqrt[gg[[i,i]]]*Sqrt[gg[[j,j]]]*
      Sqrt[gg[[k,k]]]*Sqrt[gi[[l,l]]]*temp1[[i,j,k,l]], {i, 1, 3}, {j, 1, 3}, 
     {k, 1, 3}, {l, 1, 3}]])
     
Laplace3[V_] := 
  (temp = Simplify[Table[V[[i,j,k]]/Sqrt[gg[[i,i]]*gg[[j,j]]*gg[[k,k]]], {i, 1, 3}, 
      {j, 1, 3}, {k, 1, 3}]]; 
   temp1 = Simplify[Table[D[temp[[i,j,k]], xx[[l]]] + 
       Sum[Gam[[i,sigma,l]]*temp[[sigma,j,k]], {sigma, 1, 3}] + 
       Sum[Gam[[j,sigma,l]]*temp[[i,sigma,k]], {sigma, 1, 3}] + 
       Sum[Gam[[k,sigma,l]]*temp[[i,j,sigma]], {sigma, 1, 3}], {i, 1, 3}, {j, 1, 3}, 
      {k, 1, 3}, {l, 1, 3}]]; 
   temp2 = Simplify[Table[D[temp1[[i,j,k,l]], xx[[m]]] + 
       Sum[Gam[[i,sigma,m]]*temp1[[sigma,j,k,l]], {sigma, 1, 3}] + 
       Sum[Gam[[j,sigma,m]]*temp1[[i,sigma,k,l]], {sigma, 1, 3}] + 
       Sum[Gam[[k,sigma,m]]*temp1[[i,j,sigma,l]], {sigma, 1, 3}] - 
       Sum[Gam[[sigma,l,m]]*temp1[[i,j,k,sigma]], {sigma, 1, 3}], {i, 1, 3}, {j, 1, 3}, 
      {k, 1, 3}, {l, 1, 3}, {m, 1, 3}]]; 
   Simplify[Table[Sqrt[gg[[i,i]]]*Sqrt[gg[[j,j]]]*Sqrt[gg[[k,k]]]*
      Sum[gi[[l,l]]*temp2[[i,j,k,l,l]], {l, 1, 3}], {i, 1, 3}, {j, 1, 3}, 
     {k, 1, 3}]])
     
Print["Done: grad/div/Laplace up to Rank 3"];