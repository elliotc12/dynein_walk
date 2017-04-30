#! /usr/bin/env python2.7

import numpy as np
import time, signal, sys, os, matplotlib, subprocess

if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

if tail:
    data = np.genfromtxt(data_filename, delimiter="\t", invalid_raise=False, skip_header=skiplen)
else:
    data = np.genfromtxt(data_filename, delimiter="\t", invalid_raise=False)

timesteps = len(data)

if tail and sys.stdin.isatty():
    skiplen = sum(1 for line in open(data_filename)) - 100
    if skiplen < 0:
        skiplen = 1

if len(data) == 0:
       print("Very short run!")
       exit(1)

avging_window_width = 5000
num_windows = len(data) // avging_window_width # floor division

nbys =  np.zeros(timesteps)
fbys =  np.zeros(timesteps)
times = np.empty(timesteps)

for i in range(timesteps):
    if int(data[i,0]) == 0:
        fbys[i] = data[i,16]
    elif int(data[i,0]) == 1:
        nbys[i] = data[i,16]
    times[i] = data[i,1]*1000

avg_nbys = np.array([np.mean(nbys[avging_window_width*i:avging_window_width*(i+1)]) for i in range(num_windows)])
avg_fbys = np.array([np.mean(fbys[avging_window_width*i:avging_window_width*(i+1)]) for i in range(num_windows)])
avg_times = np.array([times[int(np.floor((i+0.5)*avging_window_width))] for i in range(num_windows)])

y_min = np.min([np.min(avg_nbys), np.min(avg_fbys)])
y_max = np.max([np.max(avg_nbys), np.max(avg_fbys)])

plt.xlabel("time (ms)")
plt.ylabel("Binding domain y-projection (nm)")

plt.gca().set_ylim(y_min-1,y_max+1)

plt.plot(avg_times, avg_nbys, label="nby")
plt.plot(avg_times, avg_fbys, label="fby")

plt.legend()
plt.tight_layout()

os.system('mkdir -p plots')
plt.savefig("plots/y-trajectory-" + title + ".pdf")
