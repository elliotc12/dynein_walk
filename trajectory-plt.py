#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess

if 'show' not in sys.argv:
    matplotlib.use('Agg')

import draw_cartoon
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle

import io

tail = 'tail' in sys.argv

usage = '''
Usage: python2 TITLE %s [show] [tail]"
       show: show animation in a window while generating movie
	     omitting show makes %s faster but less exciting to watch
''' % (sys.argv[0], sys.argv[0])

if len(sys.argv) < 2:
  print(usage)
  sys.exit(1)

data_filename = sys.argv[1]
title = data_filename[data_filename.index("data/")+5:-4]

raw_lines = open(data_filename, 'r').readlines()
lines = sum(1 for l in raw_lines)
start_line = 0
end_line = int(min(1e5,lines))
plot_length = end_line - start_line - 1

raw_data = "".join(raw_lines[start_line:end_line])

data = np.genfromtxt(io.BytesIO(raw_data.encode()), delimiter="\t", invalid_raise=False)

if len(data) == 0:
       print("Very short run!")
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
    times[i] = data[i,1]*1e6

avging_window_width = 300
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
gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 3])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
ax2 = fig.add_subplot(gs[2], sharex=ax0)

plt.setp([ax0.get_xticklabels(), ax1.get_xticklabels()], visible=False)

# x projection
ax0.set_ylabel("x-projection (nm)")
ax0.set_ylim(y_min_xproj-1,y_max_xproj+1)
plt.setp(ax0.get_xticklabels(), visible=False)

ax0.plot(avg_times, avg_nbxs, label="near foot", c='b')
ax0.plot(avg_times, avg_fbxs, label="far foot", c='r')

ax0.legend(loc="upper right")

# cartoons
# ax1.axis('off')
ax1.set_aspect('equal', 'datalim')
plt.setp(ax1.get_xticklabels(), visible=False)

x_axes_size = ax1.get_xlim()[1] - ax1.get_xlim()[0]
y_axes_size = ax1.get_ylim()[1] - ax1.get_ylim()[0]

x_scaling = 0.04
y_scaling = 0.04

cartoon_draw_times_x_proj = np.array([1.0*1e-6, 3.0*1e-6, 4.899*1e-6, 7.0*1e-6, 9.0*1e-6])

plt.sca(ax1)
for t in cartoon_draw_times_x_proj:
    idx = np.where(data[:,1] == t)[0][0]
    Xs = data[idx,7:16:2]
    Ys = data[idx,8:17:2]
    draw_cartoon.draw_cartoon([t*1e6, 0], Xs, Ys, x_scaling, y_scaling)

ax1.set_ylim([0,2])
ax1.add_patch(Rectangle((0, 0), data[-1,1]*1e6, 0.4))

# y projection
ax2.set_xlabel("time ($\mu$s)")
ax2.set_ylabel("y-projection (nm)")

ax2.set_ylim(y_min_yproj-1,y_max_yproj+1)

ax2.plot(avg_times, avg_nbys, label="near foot", c='b')
ax2.plot(avg_times, avg_fbys, label="far foot", c='r')

# x_axes_size = ax2.get_xlim()[1] - ax2.get_xlim()[0]
# y_axes_size = ax2.get_ylim()[1] - ax2.get_ylim()[0]

# x_scaling = 0.005*x_axes_size
# y_scaling = 0.005*y_axes_size

# cartoon_draw_times_y_proj = np.array([0.3*1e-6, 0.6*1e-6])

# plt.sca(ax2)
# for t in cartoon_draw_times_y_proj:
#     idx = np.where(data[:,1] == t)[0][0]
#     Xs = data[idx,7:16:2]
#     Ys = data[idx,8:17:2]
#     draw_cartoon.draw_cartoon([t*1e6,-10], Xs, Ys, x_scaling, y_scaling)

gs.tight_layout(fig, h_pad=0)

os.system('mkdir -p plots')
plt.savefig("plots/trajectory-plot.pdf")
