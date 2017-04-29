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

print("argv: ", sys.argv)

if len(sys.argv) < 2:
  print usage
  sys.exit(1)

data_filename = sys.argv[1]
title = data_filename[data_filename.index("data/")+5:-4]

if tail:
    data = np.genfromtxt(data_filename, delimiter="\t", invalid_raise=False, skip_header=skiplen)
else:
    data = np.genfromtxt(data_filename, delimiter="\t", invalid_raise=False)

if tail and sys.stdin.isatty():
    skiplen = sum(1 for line in open(data_filename)) - 100
    if skiplen < 0:
	skiplen = 1

if len(data) == 0:
       print "Very short run!"
       exit(1)

nbxs = data[:,7]
fbxs = data[:,15]

x_min = np.min([np.min(nbxs), np.min(fbxs)])
x_max = np.max([np.max(nbxs), np.max(fbxs)])

plt.title("Trajectories of two feet")
plt.xlabel("time (s)")
plt.ylabel("foot position")

plt.gca().set_ylim(x_min-1,x_max+1)

plt.plot(nbxs, label="nbx")
plt.plot(fbxs, label="fbx")

plt.legend()
plt.tight_layout()

os.system('mkdir -p plots')
plt.savefig("plots/x-trajectory-" + title + ".pdf")
