(* ::Package:: *)

ClearAll["Global`*"]

Xbl= 0;
Xml= Ls*Cos[Tbl];
Xt = Xml + Lt*Cos[Tml];
Xmr = Xt + Lt*Cos[Tmr-Pi];
Xbr = Xmr + Ls*Cos[Tbr-Pi];

Ybl = 0;
Yml = Ls*Sin[Tbl];
Yt = Yml + Lt*Sin[Tml];
Ymr = Yt+Lt*Sin[Tmr-Pi];
Ybr = Ymr + Ls*Sin[Tbr-Pi];

Xvbl = 0;
Xvml = -Ls*Sin[Tbl]*Tvbl;
Xvt = Xvml + -Lt*Sin[Tml]*Tvml;
Xvmr=Xvt+-Lt*Sin[Tmr-Pi]*Tvmr;
Xvbr=Xvmr+-Ls*Sin[Tbr-Pi]*Tvbr;

Yvbl = 0;
Yvml=Ls*Cos[Tbl]*Tvbl;
Yvt=Yvml+Lt*Cos[Tml]*Tvml;
Yvmr=Yvt+Lt*Cos[Tmr-Pi]*Tvmr;
Yvbr=Yvmr+Ls*Cos[Tbr-Pi]*Tvbr;

Xbml=g*(Fxml+-Lls*(Xml-Xbl)+Llt*(Xt-Xml))+Rxml;
Xbt=g*(Ft+-Llt*(Xt-Xbl)+Lrt*(Xmr-Xt))+Rxt;
Xbmr=g*(Fmr+-Lrt*(Xmr-Xt)+Lrs*(Xbr-Xmr))+Rxmr;
Xbbr=g*(Fbr+-Lrs*(Xbr-Xmr))+Rxbr;

Ybml=g*(Fyml+-Lls*(Yml-Ybl)+Llt*(Yt-Yml))+Ryml;
Ybt=g*(Fyt+-Llt*(Yt-Yml)+Lrt*(Ymr-Yt))+Ryt;
Ybmr=g*(Fymr+-Lrt*(Ymr-Yt)+Lrs*(Ybr-Ymr))+Rmr;
Ybbr=g*(Fbr+-Lrs*(Ybr-Ymr))+Rbr;

A={{Dxbl, 0, 0, 0, g(Xml-Xbl), -g(Xt-Xml), 0, 0},
  {Dxbl, Dxml, 0, 0, 0, g(Xt-Xml), -g(Xmr-Xt), 0},
  {Dxbl, Dxml, Dxmr, 0, 0, 0, g(Xmr-Xt), -g(Xbr-Xmr)},
  {Dxbl, Dxml, Dxmr, Dxbr, 0, 0, 0, g(Xbr-Xmr)},
  {Dybl, 0, 0, 0, g(Yml-Ybl), -g(Yt-Ybl), 0, 0},
  {Dybl, Dyml, 0, 0, 0, g(Yt-Ybl), -g(Yml-Yt), 0},
  {Dybl, Dyml, Dymr, 0, 0, 0, g(Ymr-Yt), -g(Ybr-Ymr)},
  {Dybl, Dyml, Dymr, Dybr, 0, 0, 0, g(Ybr-Yml)}};

B={{g*Fxml+Rxml},
    {g*Fxt+Rxt},
    {g*Fxmr+Rxmr},
    {g*Fxbr+Rxbr},
    {g*Fyml+Ryml},
    {g*Fyt+Ryt},
    {g*Fymr+Rymr},
    {g*Fybr+Rybr}};

out = Inverse[A] . B

print 
