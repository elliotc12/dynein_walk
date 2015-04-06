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

lb_U = (0.5 * kbl * (bla - ba) **2 +
	0.5 * kml * ((s.pi - bla) + mla - ma) **2 +
	0.5 * kt  * (s.pi - mla - mra - ta) **2 +
	0.5 * kmr * ((s.pi - bra) + mra - ma) **2 +
	b_mlx*mlx + b_tx*tx + b_mrx*mrx + b_brx*brx +
	b_mly*mly + b_ty*ty + b_mry*mry + b_bry*bry)

lb_L = lb_T - lb_U

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


def euler_lagrange_equation(L, qfunc):		# Does derivative substitution for compatibility w/ old Scipy versions
    qdot = s.symbols('qdot')
    q = s.symbols('q')
    L = L.subs(s.diff(qfunc,t), qdot).subs(qfunc, q)
    dL_dq = s.diff(L, q)
    dL_dqdot = s.diff(L, qdot)
    return (s.diff(dL_dqdot, t) - dL_dq).subs(q, qfunc).subs(qdot, s.diff(qfunc,t))
    
    
# Solutions

lb_bla_deq = euler_lagrange_equation(lb_L, bla)
lb_mla_deq = euler_lagrange_equation(lb_L, mla)
lb_bra_deq = euler_lagrange_equation(lb_L, bra)
lb_mra_deq = euler_lagrange_equation(lb_L, mra)

rb_bla_deq = euler_lagrange_equation(rb_L, bla)
rb_mla_deq = euler_lagrange_equation(rb_L, mla)
rb_bra_deq = euler_lagrange_equation(rb_L, bra)
rb_mra_deq = euler_lagrange_equation(rb_L, mra)

bb_bla_deq = euler_lagrange_equation(bb_L, bla)
bb_mla_deq = euler_lagrange_equation(bb_L, mla)
bb_bra_deq = euler_lagrange_equation(bb_L, bra)
bb_mra_deq = euler_lagrange_equation(bb_L, mra)

f = open('LagrangianSolutions.txt', 'w')

print "Solving leftbound case." 
lb_sol = s.solve([lb_bla_deq, lb_mla_deq, lb_mra_deq, lb_bra_deq], dd_bla, dd_mla, dd_mra, dd_bra)
print lb_sol
f.write(lb_sol)
f.write('\n')

print "Solving rightbound case."
rb_sol = s.solve([rb_bla_deq, rb_mla_deq, rb_mra_deq, rb_bra_deq], dd_bla, dd_mla, dd_mra, dd_bra)
print rb_sol
f.write(rb_sol)
f.write('\n')

print "Solving bothbound case."
bb_sol = s.solve([bb_bla_deq, bb_mla_deq, bb_mra_deq, bb_bra_deq], dd_bla, dd_mla, dd_mra, dd_bra)
print bb_sol
f.write(bb_sol)
f.write('\n')

f.close()
