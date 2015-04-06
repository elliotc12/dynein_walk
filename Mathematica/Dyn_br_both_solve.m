(* ::Package:: *)

ClearAll["Global`*"]

mrx[t_] := Ls*Cos[bra[t]]+brx[t];
mry[t_] := Ls*Sin[bra[t]]+bry[t];

mlx[t_] := Ls*Cos[bla[t]]+blx[t];
mly[t_] := Ls*Sin[bla[t]]+bly[t];

tx[t_] :=(1/2)*(Lt*Cos[mra[t]]+mrx[t] + Lt*Cos[mla[t]]+mlx[t])
ty[t_] :=(1/2)*(Lt*Sin[mra[t]]+mry[t] + Lt*Sin[mla[t]]+mly[t])

T[t_] := 
(1/2)*mm*(mlx'[t]^2+mly'[t]^2) +
(1/2)*mt*(tx'[t]^2+ty'[t]^2) +
(1/2)*mm*(mrx'[t]^2+mry'[t]^2);

U[t_] := 
  (1/2)*kbl*(bla[t]-bA)^2+
(1/2)*kml*((Pi-bla[t])+mla[t]-mA)^2+
(1/2)*kt*(Pi - mla[t]-mra[t]-tA)^2+
(1/2)*kmr*((Pi-bra[t])+mra[t]-mA)^2 +
(1/2)*kbr*(bra[t]-bA)^2+
	Fmlx*mlx[t]+Fmly*mly[t] +
	Ftx*tx[t] + Fty*ty[t] +
	Fmrx*mrx[t]+Fmry*mry[t];

L[t_] := T[t] - U[t];

out =Solve[D[D[L[t],bla'[t]],t]==D[L[t],bla[t]] && D[D[L[t],mla'[t]],t]==D[L[t],mla[t]] && D[D[L[t],mra'[t]],t]==D[L[t],mra[t]] &&D[D[L[t],bra'[t]],t]==D[L[t],bra[t]],{bla''[t], mla''[t],mra''[t], bra''[t]}];

outfile = OpenWrite["Motion_Equations/DyneinBrownianBothboundSolutionsUnsimplified.txt"];

Put[out,outfile];
Close[outfile];
