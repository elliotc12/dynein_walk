#! /usr/bin/env python

import sympy as s


# Symbol Declaration
mm, mt, mb = s.symbols("mm, mt, mb")
kbl, kml, kt, kmr, kbr = s.symbols("kbl, kml, kt, kmr, kbr")
b_blx, b_mlx, b_tx, b_mrx, b_brx = s.symbols("b_blx, b_mlx, b_tx, b_mrx, b_brx")
b_bly, b_mly, b_ty, b_mry, b_bry = s.symbols("b_bly, b_mly, b_ty, b_mry, b_bry")
ba, ma, ta = s.symbols("ba, ma, ta")
ls, lt = s.symbols("ls, lt")

t = s.symbols('t')
bla = s.Function("bla")(t)
mla = s.Function("mla")(t)
mra = s.Function("mra")(t)
bra = s.Function("bra")(t)

blx = s.Function("blx")(t)
bly = s.Function("bly")(t)
mlx = ls * s.cos(bla) +  blx
mly = ls * s.sin(bla) +  bly
tx  = lt * s.cos(mla) +  mlx
ty  = lt * s.sin(mla) +  mly
mrx = lt * s.cos(-mra) + tx
mry = lt * s.sin(-mra) + ty
brx = ls * s.cos(-bra) + mrx
bry = ls * s.sin(-bra) + mry

d_blx = s.diff(blx, t)
d_mlx = s.diff(mlx, t)
d_tx  = s.diff(tx,  t)
d_mrx = s.diff(mrx, t)
d_brx = s.diff(brx, t)
d_bly = s.diff(blx, t)
d_mly = s.diff(mly, t)
d_ty  = s.diff(ty,  t)
d_mry = s.diff(mry, t)
d_bry = s.diff(bry, t)

dd_bla = s.diff(bla, t, t)
dd_mla = s.diff(mla, t, t)
dd_mra = s.diff(mra, t, t)
dd_bra = s.diff(bra, t, t)

# Leftbound Case ################
lb_T = \
	0.5 * mm * (d_mlx **2 + d_mly **2) + \
	0.5 * mt * (d_tx  **2 + d_ty  **2) + \
	0.5 * mm * (d_mrx **2 + d_mry **2) + \
	0.5 * mb * (d_brx **2 + d_bry **2)
	
lb_U = \
	0.5 * kbl * (bla - ba) **2 + \
	0.5 * kml * ((s.pi - bla) + mla - ma) **2 + \
	0.5 * kt  * (s.pi - mla - mra - ta) **2 + \
	0.5 * kmr * ((s.pi - bra) + mra - ma) **2 + \
	b_mlx * mlx + b_tx * tx + b_mrx * mrx + b_brx * brx +\
	b_mly * mly + b_ty * ty + b_mry * mry + b_bry * bry
	
lb_L = lb_T - lb_U
lb_diff1 = s.diff( s.diff(lb_L, s.diff(bla, t) ), t) - s.diff(lb_L, bla)
lb_diff2 = s.diff( s.diff(lb_L, s.diff(mla, t) ), t) - s.diff(lb_L, mla)
lb_diff3 = s.diff( s.diff(lb_L, s.diff(mra, t) ), t) - s.diff(lb_L, mra)
lb_diff4 = s.diff( s.diff(lb_L, s.diff(bra, t) ), t) - s.diff(lb_L, bra)

# Rightbound Case ################


rb_T = \
	0.5 * mb * (d_blx **2 + d_bly **2) + \
	0.5 * mm * (d_mlx **2 + d_mly **2) + \
	0.5 * mt * (d_tx  **2 + d_ty  **2) + \
	0.5 * mm * (d_mrx **2 + d_mry **2)
	
rb_U = \
	0.5 * kml * ((s.pi - bla) + mla - ma) **2 + \
	0.5 * kt  * (s.pi - mla - mra - ta) **2 + \
	0.5 * kmr * ((s.pi - bra) + mra - ma) **2 + \
	0.5 * kbr * (bra - ba) **2 + \
	b_blx * blx + b_mlx * mlx + b_tx * tx + b_mrx * mrx +\
	b_bly * bly + b_mly * mly + b_ty * ty + b_mry * mry
	
rb_L = rb_T - rb_U
rb_diff1 = s.diff( s.diff(rb_L, s.diff(bla, t) ), t) - s.diff(rb_L, bla)
rb_diff2 = s.diff( s.diff(rb_L, s.diff(mla, t) ), t) - s.diff(rb_L, mla)
rb_diff3 = s.diff( s.diff(rb_L, s.diff(mra, t) ), t) - s.diff(rb_L, mra)
rb_diff4 = s.diff( s.diff(rb_L, s.diff(bra, t) ), t) - s.diff(rb_L, bra)


# Bothbound Case ################

bb_T = \
	0.5 * mm * (d_mlx **2 + d_mly **2) + \
	0.5 * mt * (d_tx  **2 + d_ty  **2) + \
	0.5 * mm * (d_mrx **2 + d_mry **2)
	
bb_U = \
	0.5 * kbl * (bla - ba) **2 + \
	0.5 * kml * ((s.pi - bla) + mla - ma) **2 + \
	0.5 * kt  * (s.pi - mla - mra - ta) **2 + \
	0.5 * kmr * ((s.pi - bra) + mra - ma) **2 + \
	0.5 * kbr * (bra - ba) **2 + \
	b_mlx * mlx + b_tx * tx + b_mrx * mrx +\
	b_mly * mly + b_ty * ty + b_mry * mry
	
bb_L = bb_T - bb_U
bb_diff1 = s.diff( s.diff(bb_L, s.diff(bla, t) ), t) - s.diff(bb_L, bla)
bb_diff2 = s.diff( s.diff(bb_L, s.diff(mla, t) ), t) - s.diff(bb_L, mla)
bb_diff3 = s.diff( s.diff(bb_L, s.diff(mra, t) ), t) - s.diff(bb_L, mra)
bb_diff4 = s.diff( s.diff(bb_L, s.diff(bra, t) ), t) - s.diff(bb_L, bra)


#Solutions

lb_sol = s.solve([lb_diff1, lb_diff2, lb_diff3, lb_diff4], dd_bla, dd_mla, dd_mra, dd_bra)
rb_sol = s.solve([rb_diff1, rb_diff2, rb_diff3, rb_diff4], dd_bla, dd_mla, dd_mra, dd_bra)
bb_sol = s.solve([bb_diff1, bb_diff2, bb_diff3, bb_diff4], dd_bla, dd_mla, dd_mra, dd_bra)

print sol

f = open('LagrangianSolutions.txt', 'w')

f.write(lb_sol)
f.write('\n')
f.write(rb_sol)
f.write('\n')
f.write(bb_sol)

f.close()
