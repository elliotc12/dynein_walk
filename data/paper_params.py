#!/usr/bin/python3

k_b = 1.4e8 # s^-1
k_ub = 120 # s^-1

trajectory_k_b  = 9.5e9 # s^-1
trajectory_k_ub = 1e20 # s^-1
trajectory_seed = 14 # 1 and 10 and 14 and 17 20 are okay

cexp = -0.35

cb = 0.08 # kbT
cm = 1.05 # kbT
ct = 0.36 # kbT

# other good combos of ct/cm/cb: 0.03/1.01/0.2 and 0.23/1.36/0.09

eqb = 120
eqmpre = 197 # 204.5 from burgess
eqmpost = 242 # 222.5 from burgess, but as much as 234.5 -- they say 15nm difference for pre/post ub domain, so maybe latter is better
eqt = 0

ls = 20.75 # nm
lt = 23 # nm

radius_t = 12 # nm
radius_m = 11 # nm
radius_b = 3.5 # nm
