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

sample_points = list(map(int, np.floor(np.linspace(0, plot_length-1, num_points))))

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

gs = gridspec.GridSpec(3, 1, height_ratios=[6, 2, 5])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
ax2 = fig.add_subplot(gs[2], sharex=ax0)
plt.setp([ax0.get_xticklabels(), ax1.get_xticklabels()], visible=False)

### Trajectory plots

# cartoons
# ax1.axis('off')
ax1.set_aspect('equal', 'datalim')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)

x_axes_size = ax1.get_xlim()[1] - ax1.get_xlim()[0]
y_axes_size = ax1.get_ylim()[1] - ax1.get_ylim()[0]

x_scaling = 2e-2 #
y_scaling = 2e-2


# cartoon_draw_times_x_proj = np.array([9.241e-07, 3.0*1e-6, 4.899*1e-6, 7.0155e-06, 9.0*1e-6])
cartoon_draw_times_idxs = [sample_points[num_points*1//10], sample_points[num_points*3//10], sample_points[num_points*5//10], sample_points[num_points*7//10], sample_points[num_points*9//10]]
cartoon_draw_times_x_proj = data[:,1][cartoon_draw_times_idxs]

plt.sca(ax1)
for t in cartoon_draw_times_x_proj:
    idx = np.where(data[:,1] >= t)[0][0]
    Xs = np.abs(x_scaling*data[idx,7:16:2])
    Ys = y_scaling*data[idx,8:17:2]
    state = int(data[idx, 0])
    if state == 1:
        Xs = Xs[::-1]
        Ys = Ys[::-1]
        Xs -= Xs[0]
        Xs = Xs + t*1e6
        cartoon.draw(ax1, Xs, Ys)
        ax1.plot([Xs[4]-0.5, Xs[4]+0.5], [Ys[4], Ys[4]], linewidth=1, c="teal")
    else:
        Xs -= Xs[0]
        Xs = Xs + t*1e6
        cartoon.draw(ax1, Xs, Ys)
        ax1.plot([Xs[0]-0.5, Xs[0]+0.5], [Ys[0], Ys[0]], linewidth=1, c="teal")

ax1.axis('off')

# x projection
ax0.set_ylabel("x-projection (nm)")
ax0.set_ylim(-30,110)
plt.setp(ax0.get_xticklabels(), visible=False)

ax0.plot(times*1e6, nbxs, label="far foot", c='b')
ax0.plot(times*1e6, fbxs, label="near foot", c='r')

ax0.legend(loc="lower right")
ax0.axes.get_xaxis().set_visible(False)

# y projection
ax2.set_xlabel("time ($\mu s$)")
ax2.set_ylabel("y-projection (nm)")

ax2.set_ylim(y_min_yproj-1,y_max_yproj+1)

ax2.plot(times*1e6, nbys, label="near foot", c='b')
ax2.plot(times*1e6, fbys, label="far foot", c='r')

ax2.set_ylim(-10,75)

for t in cartoon_draw_times_x_proj:
    idx = np.where(data[:,1] >= t)[0][0]
    X = np.min([data[idx,7], data[idx,15]])
    Y = np.min([data[idx,8], data[idx,16]])
    ax0.annotate('', xy=(t*1e6, X-7), xytext=(t*1e6, X-7.01), arrowprops=dict(facecolor='black', width=1, headwidth=3, headlength=3, shrink=0))
    ax2.annotate('', xy=(t*1e6, Y-3), xytext=(t*1e6, Y-3.01), arrowprops=dict(facecolor='black', width=1, headwidth=3, headlength=3, shrink=0))

gs.tight_layout(fig, pad=0.01)

os.system('mkdir -p plots')
plt.savefig("plots/paper_trajectory_plot.pdf")
plt.show()
