(* ::Package:: *)

mlx[t_] := Ls*Cos[bla[t]] + blx;
mly[t_] := Ls*Sin[bla[t]] + bly;

tx[t_] := Lt*Cos[mla[t]] + mlx[t];
ty[t_] := Lt*Sin[mla[t]] + mly[t];

mrx[t_] := Lt*Cos[-mra[t]] + tx[t];
mry[t_] := Lt*Sin[-mra[t]] + ty[t];

brx[t_] := Ls*Cos[-bra[t]] + mrx[t];
bry[t_] := Ls*Sin[-bra[t]] + mry[t];

T[t_] := 
  (1/2)*mm*(mlx'[t]^2 + mly'[t]^2) +
   (1/2)*mt*(tx'[t]^2 + ty'[t]^2) +
   (1/2)*mm*(mrx'[t]^2 + mry'[t]^2) + 
   (1/2)*mb*(brx'[t]^2 + bry'[t]^2);

U[t_] := 
    (1/2)*kbl*(bla[t] - bA)^2 +
   (1/2)*kml*((Pi - bla[t]) + mla[t] - mA)^2 +
   (1/2)*kt*(Pi - mla[t] - mra[t] - tA)^2 +
   (1/2)*kmr*((Pi - bra[t]) + mra[t] - mA)^2 +
   	Fmlx*mlx[t] + Fmly*mly[t] +
   	Ftx*tx[t] + Fty*ty[t] +
   	Fmrx*mrx[t] + Fmry*mry[t] +
   	Fbrx*brx[t] + Fbry*bry[t];

L[t_] := T[t] - U[t];

out = Solve[
   D[D[L[t], bla'[t]], t] == D[L[t], bla[t]] && 
    D[D[L[t], mla'[t]], t] == D[L[t], mla[t]] && 
    D[D[L[t], mra'[t]], t] == D[L[t], mra[t]] && 
    D[D[L[t], bra'[t]], t] == D[L[t], bra[t]], {bla''[t], mla''[t], 
    mra''[t], bra''[t]}];

ddbla[t] = bla''[t]  /. out[[1]];
ddmla[t] = mla''[t]  /. out[[1]];
ddmra[t] = mra''[t]  /. out[[1]];
ddbra[t] = bra''[t]  /. out[[1]];

ddbla[t] = ddbla[t] /. bla[t] -> blaVal;	
ddmla[t] = ddmla[t] /. bla[t] -> blaVal;	
ddmra[t] = ddmra[t] /. bla[t] -> blaVal;	
ddbra[t] = ddbra[t] /. bla[t] -> blaVal;
ddbla[t] = ddbla[t] /. mla[t] -> mlaVal;	
ddmla[t] = ddmla[t] /. mla[t] -> mlaVal;	
ddmra[t] = ddmra[t] /. mla[t] -> blaVal;	
ddbra[t] = ddbra[t] /. mla[t] -> blaVal;
ddbla[t] = ddbla[t] /. mra[t] -> mraVal;	
ddmla[t] = ddmla[t] /. mra[t] -> mraVal;	
ddmra[t] = ddmra[t] /. mra[t] -> blaVal;	
ddbra[t] = ddbra[t] /. mra[t] -> blaVal;
ddbla[t] = ddbla[t] /. bra[t] -> braVal;	
ddmla[t] = ddmla[t] /. bra[t] -> braVal;	
ddmra[t] = ddmra[t] /. bra[t] -> blaVal;	
ddbra[t] = ddbra[t] /. bra[t] -> blaVal;

ddbla[t] = ddbla[t] /. bla'[t] -> dblaVal;	
ddmla[t] = ddmla[t] /. bla'[t] -> dblaVal;	
ddmra[t] = ddmra[t] /. bla'[t] -> dblaVal;	
ddbra[t] = ddbra[t] /. bla'[t] -> dblaVal;
ddbla[t] = ddbla[t] /. mla'[t] -> dmlaVal;	
ddmla[t] = ddmla[t] /. mla'[t] -> dmlaVal;	
ddmra[t] = ddmra[t] /. mla'[t] -> dblaVal;	
ddbra[t] = ddbra[t] /. mla'[t] -> dblaVal;
ddbla[t] = ddbla[t] /. mra'[t] -> dmraVal;	
ddmla[t] = ddmla[t] /. mra'[t] -> dmraVal;	
ddmra[t] = ddmra[t] /. mra'[t] -> dblaVal;	
ddbra[t] = ddbra[t] /. mra'[t] -> dblaVal;
ddbla[t] = ddbla[t] /. bra'[t] -> dbraVal;	
ddmla[t] = ddmla[t] /. bra'[t] -> dbraVal;	
ddmra[t] = ddmra[t] /. bra'[t] -> dblaVal;	
ddbra[t] = ddbra[t] /. bra'[t] -> dblaVal;

Lt = 10.0;
Ls = 10.0;

kt = 1.0;
kml = 1.0;
kmr = 1.0;
kbl = 1.0;
kbr = 1.0;

mb = 1.0;
mm = 1.0;
mt = 1.0;

Fblx = 0;
Fbly = 0;
Fmlx = 0;
Fmly = 0;
Ftx = 0;
Fty = 0;
Fmrx = 0;
Fmry = 0;
Fbrx = 0;
Fbry = 0;

mA = (108.0/180)*Pi;
bA = (108.0/180)*Pi;
tA = (108.0/180)*Pi;

braVal = (180.0/180)*Pi;
mlaVal = (36.0/180)*Pi;
mraVal = (36.0/180)*Pi;
blaVal = (180.0/180)*Pi;

dbraVal = 0;
dmlaVal = 0;
dmraVal = 0;
dblaVal = 0;

WriteString[$Output, "\n\nMathematica LeftBound Accelerations:\n\n"];

WriteString[$Output, "ddbla: "];
Print[ddbla[t]];

WriteString[$Output, "ddmla: "];
Print[ddmla[t]];

WriteString[$Output, "ddmra: "];
Print[ddmra[t]];

WriteString[$Output, "ddbra: "];
Print[ddbra[t]];
