#!/usr/bin/python2.7

import numpy as np
import subprocess

ls_min = 12.0 # nm
ls_max = 22.1 # nm
ls_num = 3

lt_min = 7.0 # nm
lt_max = 11.15 # nm
lt_num = 3

k_b_min = 180 # s^-1
k_b_max = 1000 # s^-1
k_b_num = 3

ls_range = np.linspace(ls_min, ls_max, num=ls_num)
lt_range = np.linspace(lt_min, lt_max, num=lt_num)
k_b_range = np.linspace(k_b_min, k_b_max, num=k_b_num)

for permutation in [{"ls": ls,"lt": lt,"k_b": k_b} for ls in ls_range for lt in lt_range for k_b in k_b_range]:
    cmd = [
        "./generate_stepping_data",
        "--Ls", str(permutation["ls"]),
        "--Lt", str(permutation["lt"]),
        "--k_b", str(permutation["k_b"]),
    ]
    print "Running: ", cmd
    subprocess.Popen(cmd)
