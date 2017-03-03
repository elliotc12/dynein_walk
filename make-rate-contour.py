#!/usr/bin/python2.7
from __future__ import division, print_function

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import sys

datafiles = [s for s in os.listdir("data/contour/") if "stepping_data" in s]

t_step = []
t_ob = []
t_bb = []
t_proc = []
kbs = []
kubs = []

for i in range(len(datafiles)):
    fname = "data/contour/" + datafiles[i]
    if os.stat(fname).st_size == 0:
        continue
    data = np.loadtxt(fname)

    start_kb_idx = datafiles[i].index('k_b-') + 4
    end_kb_idx = datafiles[i][start_kb_idx:].index(',') + start_kb_idx
    kb = float(datafiles[i][start_kb_idx:end_kb_idx])

    start_kub_idx = datafiles[i].index('k_ub-') + 5
    end_kub_idx = datafiles[i][start_kub_idx:].index(',') + start_kub_idx
    kub = float(datafiles[i][start_kub_idx:end_kub_idx])

    if len(data) < 3 or str(type(data[0])) == "<type 'numpy.float64'>":
        print("Not enough steps in data file.")
        kbs.append(kb)
        kubs.append(kub)
        # t_step.append(float("inf"))
        # t_ob.append(float("inf"))
        # t_bb.append(float("inf"))
        t_step.append(0)
        t_ob.append(10) # hokey code to say this is very large!
        t_bb.append(0)
        t_proc.append(0)
    else:
        kbs.append(kb)
        kubs.append(kub)
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

        t_step.append(np.mean(step_times))
        t_ob.append(np.mean(onebound_times))
        t_bb.append(np.mean(bothbound_times))
        t_proc.append(t_step[-1]*t_bb[-1]/t_ob[-1])

print("kbs: ", kbs)
print("kubs: ", kubs)
print("tobs: ", t_ob)
print("tbbs: ", t_bb)

min_kb = min(kbs)
max_kb = max(kbs)
min_kub = min(kubs)
max_kub = max(kubs)

maxy_tob = max(t_ob, 4.5*10**-4)
maxy_tbb = max(t_bb, 0.0595)

KB, KUB = np.meshgrid(kbs, kubs)

TBB = np.zeros((len(kbs), len(kubs)))

for i in range(len(t_ob)):
    TBB[i,i] = t_ob[i]

plt.figure()
ax = plt.gca()

ratio = np.array(t_ob) / (4.5*10**-4)
m = cm.ScalarMappable(cmap=cm.jet)
ratiomax = 10
m.set_array(np.arange(0, ratiomax, 0.01))
for i in range(len(ratio)):
    mycolor = m.cmap(ratio[i]/ratiomax)
    plot = plt.plot(kbs[i], kubs[i], '.',
                    color=mycolor, markeredgecolor=mycolor)
CB = plt.colorbar(m)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$k_b$")
ax.set_ylabel("$k_{ub}$")
ax.set_xlim([0.01*min_kb, 100*max_kb])
ax.set_ylim([0.01*min_kub, 100*max_kub])
plt.title('(Onebound time / experimental) vs un/binding rates')

plt.savefig("plots/kb-kub-contour.pdf")

# print("ratio: ", ratio)

# plt.figure()
# plt.gca().set_ylabel("$<t_{ob}> s$")
# plt.scatter(kbs, t_ob, color='b', s=50, label="$<t_{ob}>$")
# plt.plot([minx_kb*0.7, maxx_kb*1.3], [4.52*10**-4, 4.52*10**-4], color='g', label="$<t_{ob}>$ exp: $4.5*10^{-4}$ s")
# plt.scatter(empty_kbs, [0]*len(empty_kbs), color='k', marker="x", label="no steps", s=50)
# plt.title("Binding constants vs unbound time")
# plt.gca().set_xlabel("$k_b$" + " $s^{-1}$")
# plt.gca().set_xscale('log')
# plt.gca().set_xlim([minx_kb*0.8, maxx_kb*1.2])
# plt.gca().set_ylim([-0.1*maxy_tob, 1.5*maxy_tob])
# plt.legend(shadow=True)
# plt.savefig("plots/" + sys.argv[1] + "-tob-vs-kb-fit.pdf", bbox_inches='tight')

# plt.figure()
# plt.gca().set_ylabel("$<t_{bb}> s$")
# plt.scatter(kbs, t_bb, color='r', label="$<t_{bb}>$", s=50)
# plt.plot([minx_kb*0.7, maxx_kb*1.3], [0.0595, 0.0595], color='g', label="$<t_{bb}>$ exp: $0.0595$ s")
# plt.scatter(empty_kbs, [0]*len(empty_kbs), color='k', marker="x", label="no steps", s=50)
# plt.title("Binding constants vs bound time")
# plt.gca().set_xlabel("$k_b$" + " $s^{-1}$")
# plt.gca().set_xscale('log')
# plt.gca().set_xlim([minx_kb*0.8, maxx_kb*1.2])
# plt.gca().set_ylim([-0.1*maxy_tbb, 1.5*maxy_tbb])
# plt.legend(shadow=True)
# plt.savefig("plots/" + sys.argv[1] + "-tbb-vs-kb-fit.pdf", bbox_inches='tight')
