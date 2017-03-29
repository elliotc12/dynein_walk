#!/usr/bin/python2.7
from __future__ import division, print_function

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import sys

label = sys.argv[1]
datafiles = os.listdir("data/")
datafiles = [s for s in datafiles if label in s]
datafiles = [s for s in datafiles if "stepping_data" in s]

t_step = []
t_ob = []
t_bb = []
t_ob_uncertainty = []
t_bb_uncertainty = []
t_proc = []
dts = []
empty_dts = []

for i in range(len(datafiles)):
    data = np.loadtxt("data/" + datafiles[i])

    start_dt_idx = datafiles[i].index('dt-') + 3
    if ',' in datafiles[i][start_dt_idx:]:
        end_dt_idx = datafiles[i][start_dt_idx:].index(',') + start_dt_idx
    else:
        end_dt_idx = datafiles[i][start_dt_idx:].index('.') + start_dt_idx

    if len(data) < 3 or str(type(data[0])) == "<type 'numpy.float64'>":
        print("Not enough steps in data file.")
        empty_dts.append(float(datafiles[i][start_dt_idx:end_dt_idx]))
        continue
    else:
        dts.append(float(datafiles[i][start_dt_idx:end_dt_idx]))

    bind_times = np.array(data[:,1])
    unbind_times = np.array(data[:,0])
    near_positions = np.array(data[:,2])
    far_positions = np.array(data[:,3])
    near_step_idxs = near_positions[1:] != near_positions[:-1]
    far_step_idxs = far_positions[1:] != far_positions[:-1]
    near_step_lens = (near_positions[1:] - near_positions[:-1])[near_step_idxs]
    far_step_lens = (far_positions[1:] - far_positions[:-1])[far_step_idxs]

    step_times = np.array(bind_times[1:] - bind_times[:-1])
    onebound_times = bind_times - unbind_times
    bothbound_times = bind_times[1:] - unbind_times[:-1]
    step_lengths = np.concatenate((near_step_lens, far_step_lens))

    num_steps = len(step_lengths)

    t_step.append(np.mean(step_times))
    t_ob.append(np.mean(onebound_times))
    t_bb.append(np.mean(bothbound_times))
    t_proc.append(t_step[-1]*t_bb[-1]/t_ob[-1])

    t_ob_uncertainty.append(np.std(onebound_times)/np.sqrt(num_steps)*1.645) # 95% chance of true average being 1.645 stdevs from the sample average
    t_bb_uncertainty.append(np.std(bothbound_times)/np.sqrt(num_steps)*1.645)

print("dts: ", dts)
print("tobs: ", t_ob)
print("tbbs: ", t_bb)

if len(empty_dts) == 0:
    minx_dt = min(dts)
    maxx_dt = max(dts)
else:
    minx_dt = min(min(dts), min(empty_dts))
    maxx_dt = max(max(dts), max(empty_dts))

maxy_tob = max(max(t_ob), 4.5*10**-4)
maxy_tbb = max(max(t_bb), 0.0595)

plt.figure()
plt.gca().set_ylabel("$<t_{ob}> s$")
plt.scatter(dts, t_ob, color='b', s=50, label="$<t_{ob}>$")
plt.plot([minx_dt*0.7, maxx_dt*1.3], [4.52*10**-4, 4.52*10**-4], color='g', label="$<t_{ob}>$ exp: $4.5*10^{-4}$ s")
plt.scatter(empty_dts, [0]*len(empty_dts), color='k', marker="x", label="no steps", s=50)
plt.title("timestep vs unbound time")
plt.gca().set_xlabel("$k_b$" + " $s^{-1}$")
plt.gca().set_xscale('log')
plt.gca().set_xlim([minx_dt*0.8, maxx_dt*1.2])
plt.gca().set_ylim([-0.1*maxy_tob, 1.5*maxy_tob])
plt.legend(shadow=True)
plt.errorbar(dts, t_ob, yerr=t_ob_uncertainty, linestyle="None")
plt.savefig("plots/" + sys.argv[1] + "-tob-vs-dt-fit.pdf", bbox_inches='tight')

plt.figure()
plt.gca().set_ylabel("$<t_{bb}> s$")
plt.scatter(dts, t_bb, color='r', label="$<t_{bb}>$", s=50)
plt.plot([minx_dt*0.7, maxx_dt*1.3], [0.0595, 0.0595], color='g', label="$<t_{bb}>$ exp: $0.0595$ s")
plt.scatter(empty_dts, [0]*len(empty_dts), color='k', marker="x", label="no steps", s=50)
plt.title("Binding constants vs bound time")
plt.gca().set_xlabel("$k_b$" + " $s^{-1}$")
plt.gca().set_xscale('log')
plt.gca().set_xlim([minx_dt*0.8, maxx_dt*1.2])
plt.gca().set_ylim([-0.1*maxy_tbb, 1.5*maxy_tbb])
plt.legend(shadow=True)
plt.errorbar(dts, t_bb, yerr=t_bb_uncertainty, linestyle="None")
plt.savefig("plots/" + sys.argv[1] + "-tbb-vs-dt-fit.pdf", bbox_inches='tight')

# dts = np.empty(len(datafiles))
# PE_avgs = np.empty(len(datafiles))
# PE_stdvs = np.empty(len(datafiles))

# for i, f in enumerate(datafiles):
#     f_idx_start = f.index(",dt-") + 4
#     f_idx_end = f[f_idx_start:].index(".") + f_idx_start
#     dt = float(f[f_idx_start:f_idx_end])
#     data = np.loadtxt(f)
#     PE = data[:,2]
#     PE_avgs[i] = np.mean(PE)
#     PE_stdvs[i] = np.std(PE)
#     dts[i] = dt

# benchmark_PE = PE_avgs[np.argmin(dts)]

# plt.scatter(dts, PE_avgs/benchmark_PE)
# plt.errorbar(dts, PE_avgs/benchmark_PE, yerr=PE_stdvs/benchmark_PE, linestyle="None")
# plt.gca().set_xlim([-1e-11,1.3e-10])
# plt.title("PE ratio of higher dts")
# plt.show()
