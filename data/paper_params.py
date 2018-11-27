#!/usr/bin/python3

k_b = 3e6 # s^-1
k_ub = 30 # s^-1

trajectory_k_b  = 9.5e9 # s^-1
trajectory_k_ub = 1e20 # s^-1
trajectory_seed = 14 # 1 and 10 and 14 and 17 20 are okay

stroke_angle_k_b  = k_b # s^-1
stroke_angle_k_ub = k_ub # s^-1
stroke_angle_seed = 1

cexp = -0.35

cb = 0.08 # kbT
cm = 1.05 # kbT
ct = 0.36 # kbT

# other good combos of ct/cm/cb: 0.03/1.01/0.2 and 0.23/1.36/0.09

eqb = 120
eqmpre = 197
eqmpost = 242
eqt = 0

ls = 20.75 # nm
lt = 23 # nm

radius_t = 12*.5 # nm
radius_m = 11*.5 # nm
radius_b = 3.5*.5 # nm
