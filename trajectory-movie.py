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
Usage: python3 TITLE %s [show] [tail]"
       show: show animation in a window while generating movie
	     omitting show makes %s faster but less exciting to watch
''' % (sys.argv[0], sys.argv[0])

if len(sys.argv) < 2:
  print(usage)
  sys.exit(1)

os.system("rm -rf PNGs") # ensure the PNGs directory is empty.
os.system("mkdir -p PNGs") # ensure the PNGs directory exists.

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
gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax0)
ax2 = fig.add_subplot(gs[2])
gs.update(wspace=0.025, hspace=0.05)

# x projection
plt.setp([ax0.get_xticklabels()], visible=False)

ax0.set_ylabel("x-projection (nm)")
ax0.set_ylim(y_min_xproj-1,y_max_xproj+1)
plt.setp(ax0.get_xticklabels(), visible=False)

ax0.plot(avg_times, avg_nbxs, label="near foot", c='b')
ax0.plot(avg_times, avg_fbxs, label="far foot", c='r')

ax0.legend(loc="upper right")

# y projection
ax1.set_xlabel("time ($\mu$s)", labelpad=-3)
ax1.set_ylabel("y-projection (nm)")

ax1.set_ylim(y_min_yproj-1,y_max_yproj+1)

ax1.plot(avg_times, avg_nbys, label="near foot", c='b')
ax1.plot(avg_times, avg_fbys, label="far foot", c='r')

# cartoons
# ax2.axis('off')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)

x_axes_size = ax2.get_xlim()[1] - ax2.get_xlim()[0]
y_axes_size = ax2.get_ylim()[1] - ax2.get_ylim()[0]

x_scaling = 0.5
y_scaling = 1

cartoon_draw_times_x_proj = np.array([9.241e-07, 3.0*1e-6, 4.899*1e-6, 7.0155e-06, 9.0*1e-6])

# gs.tight_layout(fig, h_pad=0)
plt.sca(ax2)

i = 0
savefigframe = 0
dt = data[1,1] - data[0,1]

while i*dt < 9.0*1e-6:
    ax2.cla()

    ax2.add_patch(Rectangle((30, -5), -90, 5, alpha=0.8, zorder=-1))
    # for silhouette_time in cartoon_draw_times_x_proj:
    #     if (silhouette_time > i*dt):
    #         continue
    #     idx = np.where(data[:,1] == silhouette_time)[0][0]
    #     Xs = data[idx,7:16:2]
    #     Ys = data[idx,8:17:2]
    #     if int(data[idx, 0]) == 1:
    #         Xs = Xs[::-1]
    #         Ys = Ys[::-1]
    #     alpha = 0.6
    #     draw_cartoon.draw_cartoon([silhouette_time*1e6, 0], Xs, Ys, x_scaling, y_scaling, alpha)

    Xs = data[i,7:16:2]
    Ys = data[i,8:17:2]
    if int(data[i, 0]) == 1:
        Xs = Xs[::-1]
        Ys = Ys[::-1]
    alpha = 1.0
    draw_cartoon.draw_cartoon_movie_edition([0, 0], int(data[i, 0]), Xs, Ys, x_scaling, y_scaling, alpha)

    ax2.set_xlim([-60,30])
    ax2.set_ylim([-5,40])
    ax2.axis('off')
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)

    fname = 'PNGs/movie-%06d.png' % savefigframe
    savefigframe += 1
    plt.savefig(fname)
    sys.stdout.write("video progress: %.1f%%\r" % (i*dt/(9.0*1e-6)*100))
    sys.stdout.flush()
    i += 100

os.system('mkdir -p plots')

framerate = 24
have_avconv = True

try:
    subprocess.check_call("avconv --help > /dev/null", shell=True)
except (OSError, subprocess.CalledProcessError):
    print("Not using avconv...")
    have_avconv = False

if have_avconv:
    #movie_cmd = "avconv -loglevel quiet -y -framerate %g -i PNGs/movie-%%06d.png -b 1000k movies/movie.mp4" % framerate
    movie_cmd = "avconv -y -framerate %g -i PNGs/movie-%%06d.png -b 1000k movies/movie.mp4" % framerate
else:
    movie_cmd = "ffmpeg -loglevel quiet -y -framerate %g -i PNGs/movie-%%06d.png -b 1000k movies/movie.mp4" % framerate

print(movie_cmd)
os.system(movie_cmd) # make the movie
