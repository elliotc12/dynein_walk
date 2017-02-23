#!/usr/bin/python2.7
from __future__ import division, print_function

import matplotlib
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

start_idx = datafiles[0].index('k_ub-') + 5
end_idx = datafiles[0][start_idx:].index(',') + start_idx
kub = float(datafiles[0][start_idx:end_idx])

for i in range(len(datafiles)):
    data = np.loadtxt("data/" + datafiles[i])

    if len(data) < 3 or str(type(data[0])) == "<type 'numpy.float64'>":
        print("Not enough steps in data file.")
        continue

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

    start_idx = datafiles[i].index('k_b-') + 4
    end_idx = datafiles[i][start_idx:].index(',') + start_idx
    kbs.append(float(datafiles[i][start_idx:end_idx]))

print(t_bb)
print(t_ob)
print(kbs)

fig, axarr = plt.subplots(2, sharex=True)

axarr[0].set_ylabel("$<t_{ob}> s$")
axarr[0].scatter(kbs, t_ob, color='b', label="$<t_{ob}>$")
axarr[0].plot([min(kbs), max(kbs)], [4.52*10**-4, 4.52*10**-4], color='g', label="$<t_{bb}>$ ideal = $4.5*10^{-4}$ s")

axarr[1].set_ylabel("$<t_{bb}> s$")
axarr[1].scatter(kbs, t_bb, color='r', label="$<t_{bb}>$")
axarr[1].plot([min(kbs), max(kbs)], [0.0595, 0.0595], color='g', label="$<t_{bb}>$ ideal = $0.0595$ s")

legend = axarr[0].legend(shadow=True)
legend = axarr[1].legend(shadow=True)

plt.title("Binding constants vs un/bound time, kub = " + str(kub))
axarr[0].set_xlabel("$k_b$" + " $s^{-1}$")
axarr[0].set_xscale('log')
plt.show()
