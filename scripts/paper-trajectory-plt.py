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
from matplotlib.ticker import ScalarFormatter

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
parameters_filename = base_filename+'_stepping_parameters.tex'
title = data_filename[data_filename.index("data/")+5:-4]

run_conditions = open(parameters_filename).read()
raw_run_conditions = run_conditions.replace("\n", " ").replace("\\\\", "\\")

data = np.loadtxt(data_filename)
plot_length = min(int(5e6), len(data[:,0]))
data = data[:plot_length, :]

num_points = 300
if (plot_length < num_points):
    print("Error, need more data points to make trajectory plot.")
    exit(0)

sample_points = map(int, np.floor(np.linspace(0, plot_length-1, num_points)))

nbxs =  np.zeros(num_points)
fbxs =  np.zeros(num_points)
nbys =  np.zeros(num_points)
fbys =  np.zeros(num_points)
times = np.empty(num_points)

i = 0
for p in sample_points:
    if int(data[p,0]) == 0 or int(data[p,0]) == 2:
        nbxs[i] = data[p,7]
        fbxs[i] = data[p,15]
        fbys[i] = data[p,16]
    elif int(data[p,0]) == 1:
        nbxs[i] = data[p,15]
        fbxs[i] = data[p,7]
        nbys[i] = data[p,16]
    times[i] = data[p,1]
    i += 1

y_min_xproj = np.min([np.min(nbxs), np.min(fbxs)])
y_max_xproj = np.max([np.max(nbxs), 20])

y_min_yproj = np.min([np.min(nbys), np.min(fbys)])
y_max_yproj = np.max([np.max(nbys), np.max(fbys)])

fig = plt.figure()
plt.rc('text', usetex=True)

gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
plt.setp([ax0.get_xticklabels()], visible=False)

### Trajectory plots

# x projection
ax0.set_ylabel("x-projection (nm)")
# ax0.set_ylim(y_min_xproj-1,y_max_xproj+1)
ax0.set_ylim(-150,150)
plt.setp(ax0.get_xticklabels(), visible=False)

ax0.plot(times*1e6, nbxs, label="near foot", c='b')
ax0.plot(times*1e6, fbxs, label="far foot", c='r')

ax0.legend(loc="lower right")
ax0.axes.get_xaxis().set_visible(False)

# ax = plt.gca().xaxis
# ax.set_major_formatter(ScalarFormatter()) 
# ax0.set_xticklabels([])

# y projection
ax1.set_xlabel("time ($\mu s$)")
ax1.set_ylabel("y-projection (nm)")

ax1.set_ylim(y_min_yproj-1,y_max_yproj+1)

ax1.plot(times*1e6, nbys, label="near foot", c='b')
ax1.plot(times*1e6, fbys, label="far foot", c='r')

ax1.set_ylim(-20,65)

# ax = plt.gca().xaxis
# ax.set_major_formatter(ScalarFormatter())
# ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

gs.tight_layout(fig, pad=2)

# plt.gcf().suptitle(
#     raw_run_conditions +
#     r' $k_{b}: \kb, k_{ub}: \kub$',
#     fontsize=14)

os.system('mkdir -p plots')
plt.savefig("plots/paper_trajectory_plot.pdf")
plt.show()
