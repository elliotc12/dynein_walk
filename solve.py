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

def euler_lagrange_equation(L, qfunc):
    qdot = s.symbols('qdot')
    q = s.symbols('q')
    L = L.subs(s.diff(qfunc,t), qdot).subs(qfunc, q)
    dL_dq = s.diff(L, q)
    dL_dqdot = s.diff(L, qdot)
    return (s.diff(dL_dqdot, t) - dL_dq).subs(q, qfunc).subs(qdot, s.diff(qfunc,t))

el_bla = euler_lagrange_equation(L, bla)
el_mla = euler_lagrange_equation(L, mla)
el_bra = euler_lagrange_equation(L, bra)
el_mra = euler_lagrange_equation(L, mra)

def make_prettier(e):
    e = e.subs(s.diff(bla,t), s.symbols('bladot'))
    e = e.subs(s.diff(mla,t), s.symbols('mladot'))
    e = e.subs(s.diff(bra,t), s.symbols('bradot'))
    e = e.subs(s.diff(mra,t), s.symbols('mradot'))

    e = e.subs(s.diff(bla,t,t), s.symbols('blaDDot'))
    e = e.subs(s.diff(mla,t,t), s.symbols('mlaDDot'))
    e = e.subs(s.diff(bra,t,t), s.symbols('braDDot'))
    e = e.subs(s.diff(mra,t,t), s.symbols('mraDDot'))
    return e

#print make_prettier(el_bla)

print 'starting to solve these silly things...'
sol = s.solve([el_bla, el_mla, el_bra, el_mra], dd_bla, dd_mla, dd_mra, dd_bra)

print sol

f = fopen('LagrangianSolutions.txt', 'w')
f.write(sol)
f.close()
