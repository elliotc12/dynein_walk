#!/usr/bin/python2.7
from __future__ import division, print_function

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import sys

datafiles = sys.argv[1:]
dts = np.empty(len(datafiles))
PE_avgs = np.empty(len(datafiles))
PE_stdvs = np.empty(len(datafiles))

for i, f in enumerate(datafiles):
    f_idx_start = f.index(",dt-") + 4
    f_idx_end = f[f_idx_start:].index(".") + f_idx_start
    dt = float(f[f_idx_start:f_idx_end])
    data = np.loadtxt(f)
    PE = data[:,2]
    PE_avgs[i] = np.mean(PE)
    PE_stdvs[i] = np.std(PE)
    dts[i] = dt

benchmark_PE = PE_avgs[np.argmin(dts)]

plt.scatter(dts, PE_avgs/benchmark_PE)
plt.errorbar(dts, PE_avgs/benchmark_PE, yerr=PE_stdvs/benchmark_PE, linestyle="None")
plt.gca().set_xlim([-1e-11,1.3e-10])
plt.title("PE ratio of higher dts")
plt.show()
