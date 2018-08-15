#!/usr/bin/python3

k_b = 9e7 # s^-1
k_ub = 75 # s^-1

trajectory_k_b  = 9.5e9 # s^-1
trajectory_k_ub = 1e20 # s^-1
trajectory_seed = 14 # 1 and 10 and 14 and 17 20 are okay

cexp = -0.35

cb = 0.1 # kbT
cm = 0.4 # kbT
ct = 0.2 # kbT

eqb = 120
eqmpre = 197 # 204.5 from burgess
eqmpost = 242 # 222.5 from burgess, but as much as 234.5 -- they say 15nm difference for pre/post ub domain, so maybe latter is better
eqt = 0

ls = 20.75 # nm
lt = 23 # nm

radius_t = 8 # nm
radius_m = 11 # nm
radius_b = 3.5 # nm
