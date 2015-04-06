(* ::Package:: *)

ClearAll["Global`*"]
mrx[t_] := Ls*Cos[bra[t]]+brx;
mry[t_] := Ls*Sin[bra[t]]+bry;

tx[t_] := Lt*Cos[mra[t]]+mrx[t];
ty[t_]:= Lt*Sin[mra[t]]+mry[t];

mlx[t_] := Lt*Cos[-mla[t]]+tx[t];
mly[t_] := Lt*Sin[-mla[t]]+ty[t];

blx[t_] := Ls*Cos[-bla[t]]+mlx[t];
bly[t_] := Ls*Sin[-bla[t]]+mly[t];

T[t_] := 
(1/2)*mm*(mrx'[t]^2+mry'[t]^2) +
(1/2)*mt*(tx'[t]^2+ty'[t]^2) +
(1/2)*mm*(mlx'[t]^2+mly'[t]^2) + 
(1/2)*mb*(blx'[t]^2+bly'[t]^2);

U[t_] := 
  (1/2)*kbr*(bra[t]-bA)^2+
(1/2)*kmr*((Pi-bra[t])+mra[t]-mA)^2+
(1/2)*kt*(Pi - mra[t]-mla[t]-tA)^2+
(1/2)*kml*((Pi-bla[t])+mla[t]-mA)^2 +
	Fmrx*mrx[t]+Fmry*mry[t] +
	Ftx*tx[t] + Fty*ty[t] +
	Fmlx*mlx[t]+Fmly*mly[t] +
	Fblx*blx[t]+Fbly*bly[t];

L[t_] := T[t] - U[t];

out =Solve[D[D[L[t],bla'[t]],t]==D[L[t],bla[t]] && D[D[L[t],mla'[t]],t]==D[L[t],mla[t]] && D[D[L[t],mra'[t]],t]==D[L[t],mra[t]] &&D[D[L[t],bra'[t]],t]==D[L[t],bra[t]],{bla''[t], mla''[t],mra''[t], bra''[t]}];

outfile = OpenWrite["Motion_Equations/DyneinBrownianRightboundSolutionsUnsimplified.txt"];

Put[out,outfile];
Close[outfile];
