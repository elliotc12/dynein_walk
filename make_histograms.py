#!/usr/bin/python2.7

import numpy as np
import subprocess

assert(subprocess.call("make histogram-stuff", shell=True) == 0)

have_slurm = True

try:
    subprocess.check_call("slurm", shell=True)
except (OSError, subprocess.CalledProcessError):
    print "Not using slurm..."
    have_slurm = False

ls_min = 12.0 # nm
ls_max = 12.0 # nm
ls_num = 1

lt_min = 9.075 # nm
lt_max = 9.075 # nm
lt_num = 1

k_b_min = 180 # s^-1
k_b_max = 180 # s^-1
k_b_num = 1

T_min = 310.15 # K
T_max = 310.15 # K
T_num = 1

ls_range = np.linspace(ls_min, ls_max, num=ls_num)
lt_range = np.linspace(lt_min, lt_max, num=lt_num)
k_b_range = np.linspace(k_b_min, k_b_max, num=k_b_num)
T_range = np.linspace(T_min, T_max, num=T_num)

for permutation in [{"ls": ls,"lt": lt,"k_b": k_b, "T": T} for ls in ls_range for lt in lt_range for k_b in k_b_range for T in T_range]:
    cmd = ["srun"] if have_slurm else []
    cmd.extend([
        "./generate_stepping_data",
        "--Ls", str(permutation["ls"]),
        "--Lt", str(permutation["lt"]),
        "--k_b", str(permutation["k_b"]),
        "--T", str(permutation["T"]),
        "--movie"
    ])
    print "Running: ", cmd

    basename = 'ls-%.3g,lt-%.3g,k_b-%s,cb-%s,cm-%s,ct-%s,T-%s' % (permutation['ls'],
                                                                  permutation['lt'],
                                                                  permutation["k_b"],
                                                                  'cb', 'cm', 'ct', permutation['T'])
    out = open('runlogs/' + basename + '.out', 'w')
    subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
