#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess, re

if 'show' not in sys.argv:
    matplotlib.use('Agg')

#import draw.balls as cartoon
import draw.cartoon as cartoon
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle

import io

tail = 'tail' in sys.argv

usage = '''
Usage: python2 %s MOVIE_FILENAME STEPS_FILENAME [show]"
       show: show plot in a window
''' % (sys.argv[0])

if len(sys.argv) < 3:
  print(usage)
  sys.exit(1)

data_filename = sys.argv[1]
title = data_filename[data_filename.index("data/")+5:-4]

raw_lines = open(data_filename, 'r').readlines()
lines = sum(1 for l in raw_lines)
start_line = 0
end_line = int(min(1e6, lines))
plot_length = end_line - start_line - 1

raw_data = "".join(raw_lines[start_line:end_line])

run_conditions = open("data/exploration_stepping_parameters.tex").read()
raw_run_conditions = run_conditions.replace("\n", " ").replace("\\\\", "\\")

data = np.genfromtxt(io.BytesIO(raw_data.encode()), delimiter="\t", invalid_raise=False)

avging_window_width = 300
if len(data) < avging_window_width:
       print("Need at least ", avging_window_width, " data points!")
       exit(1)

nbxs =  np.zeros(plot_length)
fbxs =  np.zeros(plot_length)
nbys =  np.zeros(plot_length)
fbys =  np.zeros(plot_length)

times = np.empty(plot_length)

for i in range(plot_length):
    if int(data[i,0]) == 0 or int(data[i,0]) == 2:
        nbxs[i] = data[i,7]
        fbxs[i] = data[i,15]
        fbys[i] = data[i,16]
    elif int(data[i,0]) == 1:
        nbxs[i] = data[i,15]
        fbxs[i] = data[i,7]
        nbys[i] = data[i,16]
    times[i] = data[i,1]

num_windows = plot_length // avging_window_width # floor division

avg_nbxs = np.array([np.mean(nbxs[avging_window_width*i:avging_window_width*(i+1)]) for i in range(num_windows)])
avg_fbxs = np.array([np.mean(fbxs[avging_window_width*i:avging_window_width*(i+1)]) for i in range(num_windows)])
avg_nbys = np.array([np.mean(nbys[avging_window_width*i:avging_window_width*(i+1)]) for i in range(num_windows)])
avg_fbys = np.array([np.mean(fbys[avging_window_width*i:avging_window_width*(i+1)]) for i in range(num_windows)])

avg_times = np.array([times[int(np.floor((i+0.5)*avging_window_width))] for i in range(num_windows)])

y_min_xproj = np.min([np.min(avg_nbxs), np.min(avg_fbxs)])
y_max_xproj = np.max([np.max(avg_nbxs), 20])

y_min_yproj = np.min([np.min(avg_nbys), np.min(avg_fbys)])
y_max_yproj = np.max([np.max(avg_nbys), np.max(avg_fbys)])

fig = plt.figure()
plt.rc('text', usetex=True)

gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
plt.setp([ax0.get_xticklabels()], visible=False)

### Trajectory plots

# x projection
ax0.set_ylabel("x-projection (nm)")
ax0.set_ylim(y_min_xproj-1,y_max_xproj+1)
plt.setp(ax0.get_xticklabels(), visible=False)

ax0.plot(avg_times, avg_nbxs, label="near foot", c='b')
ax0.plot(avg_times, avg_fbxs, label="far foot", c='r')

ax0.legend(loc="upper right")

# y projection
ax1.set_xlabel("time (s)")
ax1.set_ylabel("y-projection (nm)")

ax1.set_ylim(y_min_yproj-1,y_max_yproj+1)

ax1.plot(avg_times, avg_nbys, label="near foot", c='b')
ax1.plot(avg_times, avg_fbys, label="far foot", c='r')

gs.tight_layout(fig, pad=2)

plt.gcf().suptitle(
    raw_run_conditions +
    r' $k_{b}: \kb, k_{ub}: \kub$',
    fontsize=14)

os.system('mkdir -p plots')
plt.savefig("plots/exploration-plot.pdf")
plt.show()

### Histogram plots
stepdata = np.loadtxt(sys.argv[2])

bind_times = np.array(stepdata[:,1])
unbind_times = np.array(stepdata[:,0])
near_positions = np.array(stepdata[:,2])
far_positions = np.array(stepdata[:,3])
near_step_idxs = near_positions[1:] != near_positions[:-1]
far_step_idxs = far_positions[1:] != far_positions[:-1]
near_step_lens = (near_positions[1:] - near_positions[:-1])[near_step_idxs]
far_step_lens = (far_positions[1:] - far_positions[:-1])[far_step_idxs]

step_times = np.array(bind_times[1:] - bind_times[:-1])
onebound_times = bind_times - unbind_times
bothbound_times = unbind_times[1:] - bind_times[:-1]
step_lengths = np.concatenate((near_step_lens, far_step_lens))

for i in range(len(onebound_times)):
    print("lifted off at {}, bound at {}, ob_time: {}".format(unbind_times[i], bind_times[i], onebound_times[i]))

for j in range(len(bothbound_times)):
    print("bound at {}, lifted off at {}, bb_time: {}".format(
        bind_times[:-1][j], unbind_times[1:][j], bothbound_times[j]))

num_steps = len(step_lengths)

t_step = []
t_ob = []
t_bb = []
t_ob_uncertainty = []
t_bb_uncertainty = []
t_proc = []

t_step.append(np.mean(step_times))
t_ob.append(np.mean(onebound_times))
t_bb.append(np.mean(bothbound_times))
t_proc.append(t_step[-1]*t_bb[-1]/t_ob[-1])

t_ob_uncertainty.append(np.std(onebound_times)/np.sqrt(num_steps)*1.645) # 95% chance of true average being 1.645 stdevs from the sample average
t_bb_uncertainty.append(np.std(bothbound_times)/np.sqrt(num_steps)*1.645)

#step length histogram

fig = plt.figure()

plt.hist(step_lengths, bins=50)
plt.xlabel("Step length (nm)")
plt.ylabel("Frequency")
plt.savefig("plots/stepping_length_histogram.pdf", format="pdf")
plt.close(fig)

#step time histogram
fig = plt.figure()
plt.rc('text', usetex=True)

gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])

ax0.hist(step_times, bins=50)
ax0.set_title("Step times")
ax0.set_ylabel("Frequency")
ax0.set_yscale("log")
ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax1.hist(onebound_times, bins=50)
ax1.set_title("onebound times (theory 0.0595s)")
ax1.set_ylabel("Frequency")
ax1.set_yscale("log")
ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax2.hist(bothbound_times, bins=50)
ax2.set_title("bothbound times (theory: 4.52e-4s)")
ax2.set_xlabel("Step time (s)")
ax2.set_ylabel("Frequency")
ax2.set_yscale("log")
ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

print(onebound_times)
print(bothbound_times)

plt.gcf().suptitle(
    raw_run_conditions +
    r' $k_{b}: \kb, k_{ub}: \kub, runtime: \runtime$',
    fontsize=14)

plt.subplots_adjust(hspace=0.6)

plt.show()
plt.savefig("plots/stepping_time_histogram.pdf", format="pdf")
plt.close(fig)
