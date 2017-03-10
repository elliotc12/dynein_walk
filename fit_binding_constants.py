#!/usr/bin/python2.7
from __future__ import division, print_function

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

# def add_config_info(config_txt):
#     ax = plt.gca()
#     box = ax.get_position()
#     ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
#     box = TextArea(config_txt, textprops=dict(size=10))
#     anchored_box = AnchoredOffsetbox(loc=2, child=box, bbox_to_anchor=(1,1),
#                                 bbox_transform=ax.transAxes, frameon=False)
#     ax.add_artist(anchored_box)

label = sys.argv[1]
datafiles = os.listdir("data/")
datafiles = [s for s in datafiles if label in s]
datafiles = [s for s in datafiles if "stepping_data" in s]

t_step = []
t_ob = []
t_bb = []
t_proc = []
kbs = []
kubs = []
empty_kbs = []
empty_kubs = []

for i in range(len(datafiles)):
    data = np.loadtxt("data/" + datafiles[i])

    start_kb_idx = datafiles[i].index('k_b-') + 4
    end_kb_idx = datafiles[i][start_kb_idx:].index(',') + start_kb_idx

    start_kub_idx = datafiles[i].index('k_ub-') + 5
    end_kub_idx = datafiles[i][start_kub_idx:].index(',') + start_kub_idx

    if len(data) < 3 or str(type(data[0])) == "<type 'numpy.float64'>":
        print("Not enough steps in data file.")
        empty_kbs.append(float(datafiles[i][start_kb_idx:end_kb_idx]))
        empty_kubs.append(float(datafiles[i][start_kub_idx:end_kub_idx]))
        continue
    else:
        kbs.append(float(datafiles[i][start_kb_idx:end_kb_idx]))
        kubs.append(float(datafiles[i][start_kub_idx:end_kub_idx]))

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

print("kbs: ", kbs)
print("kubs: ", kubs)
print("tobs: ", t_ob)
print("tbbs: ", t_bb)

if len(empty_kbs) == 0:
    minx_kb = min(kbs)
    maxx_kb = max(kbs)
else:
    minx_kb = min(min(kbs), min(empty_kbs))
    maxx_kb = max(max(kbs), max(empty_kbs))

maxy_tob = max(max(t_ob), 4.5*10**-4)
maxy_tbb = max(max(t_bb), 0.0595)

plt.figure()
plt.gca().set_ylabel("$<t_{ob}> s$")
plt.scatter(kbs, t_ob, color='b', s=50, label="$<t_{ob}>$")
plt.plot([minx_kb*0.7, maxx_kb*1.3], [4.52*10**-4, 4.52*10**-4], color='g', label="$<t_{ob}>$ exp: $4.5*10^{-4}$ s")
plt.scatter(empty_kbs, [0]*len(empty_kbs), color='k', marker="x", label="no steps", s=50)
plt.title("Binding constants vs unbound time")
plt.gca().set_xlabel("$k_b$" + " $s^{-1}$")
plt.gca().set_xscale('log')
plt.gca().set_xlim([minx_kb*0.8, maxx_kb*1.2])
plt.gca().set_ylim([-0.1*maxy_tob, 1.5*maxy_tob])
plt.legend(shadow=True)
plt.savefig("plots/" + sys.argv[1] + "-tob-vs-kb-fit.pdf", bbox_inches='tight')

plt.figure()
plt.gca().set_ylabel("$<t_{bb}> s$")
plt.scatter(kbs, t_bb, color='r', label="$<t_{bb}>$", s=50)
plt.plot([minx_kb*0.7, maxx_kb*1.3], [0.0595, 0.0595], color='g', label="$<t_{bb}>$ exp: $0.0595$ s")
plt.scatter(empty_kbs, [0]*len(empty_kbs), color='k', marker="x", label="no steps", s=50)
plt.title("Binding constants vs bound time")
plt.gca().set_xlabel("$k_b$" + " $s^{-1}$")
plt.gca().set_xscale('log')
plt.gca().set_xlim([minx_kb*0.8, maxx_kb*1.2])
plt.gca().set_ylim([-0.1*maxy_tbb, 1.5*maxy_tbb])
plt.legend(shadow=True)
plt.savefig("plots/" + sys.argv[1] + "-tbb-vs-kb-fit.pdf", bbox_inches='tight')

if len(empty_kubs) == 0:
    minx_kub = min(kubs)
    maxx_kub = max(kubs)
else:
    minx_kub = min(min(kubs), min(empty_kubs))
    maxx_kub = max(max(kubs), max(empty_kubs))

plt.figure()
plt.gca().set_ylabel("$<t_{ob}> s$")
plt.scatter(kubs, t_ob, color='b', s=50, label="$<t_{ob}>$")
plt.plot([minx_kub*0.7, maxx_kub*1.3], [4.52*10**-4, 4.52*10**-4], color='g', label="$<t_{ob}>$ exp: $4.5*10^{-4}$ s")
plt.scatter(empty_kubs, [0]*len(empty_kubs), color='k', marker="x", label="no steps", s=50)
plt.title("Unbinding constants vs unbound time")
plt.gca().set_xlabel("$k_{ub}$" + " $s^{-1}$")
plt.gca().set_xscale('log')
plt.gca().set_xlim([minx_kub*0.8, maxx_kub*1.2])
plt.gca().set_ylim([-0.1*maxy_tob, 1.5*maxy_tob])
plt.legend(shadow=True)
plt.savefig("plots/" + sys.argv[1] + "-tob-vs-kub-fit.pdf", bbox_inches='tight')

plt.figure()
plt.gca().set_ylabel("$<t_{bb}> s$")
plt.scatter(kubs, t_bb, color='r', label="$<t_{bb}>$", s=50)
plt.plot([minx_kub*0.7, maxx_kub*1.3], [0.0595, 0.0595], color='g', label="$<t_{bb}>$ exp: $0.0595$ s")
plt.scatter(empty_kubs, [0]*len(empty_kubs), color='k', marker="x", label="no steps", s=50)
plt.title("Unbinding constants vs bound time")
plt.gca().set_xlabel("$k_{ub}$" + " $s^{-1}$")
plt.gca().set_xscale('log')
plt.gca().set_xlim([minx_kub*0.8, maxx_kub*1.2])
plt.gca().set_ylim([-0.1*maxy_tbb, 1.5*maxy_tbb])
plt.legend(shadow=True)
plt.savefig("plots/" + sys.argv[1] + "-tbb-vs-kub-fit.pdf", bbox_inches='tight')
