#!/usr/bin/python2.7
from __future__ import division, print_function

import matplotlib
matplotlib.use("PDF") # use PDF generating backend bc it can handle bigger datasets
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea


def add_config_info(config_txt):
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
    box = TextArea(config_txt, textprops=dict(size=10))
    anchored_box = AnchoredOffsetbox(loc=2, child=box, bbox_to_anchor=(1,1),
                                bbox_transform=ax.transAxes, frameon=False)
    ax.add_artist(anchored_box)

data_name = sys.argv[1]

data = np.loadtxt("data/stepping_data_" + data_name + ".txt")

config_txt = open("data/stepping_config_" + data_name + ".txt", "r").read()

if len(data) < 3 or str(type(data[0])) == "<type 'numpy.float64'>":
    print("Not enough steps in data file.")
    exit(0)

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

#step length histogram
fig = plt.figure()
plt.hist(step_lengths, bins=12)
plt.title("Stepping length histogram " + data_name)
plt.xlabel("Step length (nm)")
plt.ylabel("Frequency")
add_config_info(config_txt + "\n\n<step length>: %.2e" % np.mean(step_lengths))
plt.savefig("plots/stepping_length_histogram_" + data_name + ".pdf", format="pdf")
plt.close(fig)

#step time histogram
fig = plt.figure()
plt.hist(step_times, bins=12)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title("Stepping time histogram " + data_name)
plt.xlabel("Step time (s)")
plt.ylabel("Frequency")

t_step = np.mean(step_times)
t_ob = np.mean(onebound_times)
t_bb = np.mean(bothbound_times)
t_proc = t_step*t_bb/t_ob

config_txt += r"""

$\langle$step times$\rangle$: %.2g s""" % np.mean(step_times)
config_txt += r"""
$\langle$# steps processivity$\rangle$: %.2g""" % (t_proc/t_step)
config_txt += r"""
$\left<t_{ob}\right>$: %.2g s""" % t_ob
config_txt += r"""
$\left<t_{bb}\right>$: %.2g s""" % t_bb
config_txt += r"""
$\left<t_{proc}\right>$: %.2g s""" % t_proc
config_txt += r"""
$f_{ob}$ = %.2g s""" % (t_ob/t_step)
add_config_info(config_txt)
plt.savefig("plots/stepping_time_histogram_" + data_name + ".pdf", format="pdf")
plt.close(fig)

print('t_step = ', t_step, 's')
print('t_ob = ', t_ob, 's')
print('t_bb = ', t_bb, 's')
print('t_proc = ', t_proc, 's')
