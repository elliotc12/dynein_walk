(* ::Package:: *)

ClearAll["Global`*"]
mlx[t_] := Ls*Cos[bla[t]]+blx;
mly[t_] := Ls*Sin[bla[t]]+bly;

tx[t_] := Lt*Cos[mla[t]]+mlx[t];
ty[t_] := Lt*Sin[mla[t]]+mly[t];

mrx[t_] := Lt*Cos[-mra[t]]+tx[t];
mry[t_] := Lt*Sin[-mra[t]]+ty[t];

brx[t_] := Ls*Cos[-bra[t]]+mrx[t];
bry[t_] := Ls*Sin[-bra[t]]+mry[t];

T[t_] := 
(1/2)*mm*(mlx'[t]^2+mly'[t]^2) +
(1/2)*mt*(tx'[t]^2+ty'[t]^2) +
(1/2)*mm*(mrx'[t]^2+mry'[t]^2) + 
(1/2)*mb*(brx'[t]^2+bry'[t]^2);

U[t_] := 
  (1/2)*kbl*(bla[t]-bA)^2+
(1/2)*kml*((Pi-bla[t])+mla[t]-mA)^2+
(1/2)*kt*(Pi - mla[t]-mra[t]-tA)^2+
(1/2)*kmr*((Pi-bra[t])+mra[t]-mA)^2 +
	Fmlx*mlx[t]+Fmly*mly[t] +
	Ftx*tx[t] + Fty*ty[t] +
	Fmrx*mrx[t]+Fmry*mry[t] +
	Fbrx*brx[t]+Fbry*bry[t];

L[t_] := T[t] - U[t];

out =Solve[D[D[L[t],bla'[t]],t]==D[L[t],bla[t]] && D[D[L[t],mla'[t]],t]==D[L[t],mla[t]] && D[D[L[t],mra'[t]],t]==D[L[t],mra[t]] &&D[D[L[t],bra'[t]],t]==D[L[t],bra[t]],{bla''[t], mla''[t],mra''[t], bra''[t]}];

outfile = OpenWrite["Motion_Equations/DyneinBrownianLeftboundSolutionsUnsimplified.txt"];

Put[out,outfile];
Close[outfile];
