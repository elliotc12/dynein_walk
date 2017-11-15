#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess, re

if 'show' not in sys.argv:
    matplotlib.use('Agg')

#import draw.balls as cartoon
import dynein.draw.cartoon as cartoon
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle

import io

tail = 'tail' in sys.argv

usage = '''
Usage: python2 %s BASE_FILENAME [show]"
       show: show plot in a window
       note: BASE_FILENAME is typically something like data/paper or data/thesis to
             which _movie_data.txt or similar is appended.
''' % (sys.argv[0])

if len(sys.argv) < 2:
  print(usage)
  sys.exit(1)

base_filename = sys.argv[1]
data_filename = base_filename+'_movie_data.txt'
step_filename = base_filename+'_stepping_data.txt'
title = data_filename[data_filename.index("data/")+5:-4]

run_conditions = open(parameters_filename).read()
raw_run_conditions = run_conditions.replace("\n", " ").replace("\\\\", "\\")

data = np.loadtxt(data_filename)
plot_length = min(int(5e6), len(data[:,0]))
data = data[:plot_length, :]

nbxs =  np.zeros(plot_length)
fbxs =  np.zeros(plot_length)
nbys =  np.zeros(plot_length)
fbys =  np.zeros(plot_length)

times = np.empty(plot_length)

for i in range(plot_length-1):
    if int(data[i,0]) == 0 or int(data[i,0]) == 2:
        nbxs[i] = data[i,7]
        fbxs[i] = data[i,15]
        fbys[i] = data[i,16]
    elif int(data[i,0]) == 1:
        nbxs[i] = data[i,15]
        fbxs[i] = data[i,7]
        nbys[i] = data[i,16]
    times[i] = data[i,1]

num_points = 100
if (plot_length < num_points):
    print("Error, need more data points to make trajectory plot.")
    exit(0)
sample_points = np.floor(np.linspace(0, plot_length-1, num_points))

avg_nbxs = np.array([nbxs[int(t)] for t in sample_points])
avg_fbxs = np.array([fbxs[int(t)] for t in sample_points])
avg_nbys = np.array([nbys[int(t)] for t in sample_points])
avg_fbys = np.array([fbys[int(t)] for t in sample_points])
avg_times = np.array([times[int(t)] for t in sample_points])


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
plt.savefig("plots/paper-trajectory-plot.pdf")
plt.show()
