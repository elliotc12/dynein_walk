#! /usr/bin/env python

import sympy as s

blx, bly = s.symbols("blx, bly")
mm, mt, mb = s.symbols("mm, mt, mb")
kbl, kml, kt, kmr = s.symbols("kbl, kml, kt, kmr")
ba, ma, ta = s.symbols("ba, ma, ta")
ls, lt = s.symbols("ls, lt")

t = s.symbols('t')
bla = s.Function("bla")(t)
mla = s.Function("mla")(t)
mra = s.Function("mra")(t)
bra = s.Function("bra")(t)

mlx = ls * s.cos(bla) +  blx
mly = ls * s.sin(bla) +  bly
tx  = lt * s.cos(mla) +  mlx
ty  = lt * s.sin(mla) +  mly
mrx = lt * s.cos(-mra) + tx
mry = lt * s.sin(-mra) + ty
brx = ls * s.cos(-bra) + mrx
bry = ls * s.sin(-bra) + mry

d_mlx = s.diff(mlx, t)
d_tx  = s.diff(tx,  t)
d_mrx = s.diff(mrx, t)
d_brx = s.diff(brx, t)
d_mly = s.diff(mly, t)
d_ty  = s.diff(ty,  t)
d_mry = s.diff(mry, t)
d_bry = s.diff(bry, t)

dd_bla = s.diff(bla, t, t)
dd_mla = s.diff(mla, t, t)
dd_mra = s.diff(mra, t, t)
dd_bra = s.diff(bra, t, t)


T = 0.5 * mm * (d_mlx **2 + d_mly **2) + \
	0.5 * mt * (d_tx  **2 + d_ty  **2) + \
	0.5 * mm * (d_mrx **2 + d_mry **2) + \
	0.5 * mb * (d_brx **2 + d_mry **2)
	
U = 0.5 * kbl * (bla - ba) **2 + \
	0.5 * kml * ((s.pi - bla) + mla - ma) **2 + \
	0.5 * kt  * (s.pi - mla - mra - ta) **2 + \
	0.5 * kmr * ((s.pi - bra) + mra - ma) **2 \
	
	
L = T - U

diff1 = s.diff( s.diff(L, s.diff(bla, t) ), t) - s.diff(L, bla)
diff2 = s.diff( s.diff(L, s.diff(mla, t) ), t) - s.diff(L, mla)
diff3 = s.diff( s.diff(L, s.diff(mra, t) ), t) - s.diff(L, mra)
diff4 = s.diff( s.diff(L, s.diff(bra, t) ), t) - s.diff(L, bra)


sol = s.solve([diff1, diff2, diff3, diff4], dd_bla, dd_mla, dd_mra, dd_bra)

print sol

f = fopen('LagrangianSolutions.txt', 'w')
f.write(sol)
f.close()
