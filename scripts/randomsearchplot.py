#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess, re

# matplotlib.use('TkAgg')
# matplotlib.use('pdf')

import argparse
import datetime
import dynein.draw.cartoon as cartoon
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import gridspec
from matplotlib.patches import Rectangle

import dynein.data as data
import io

plt.rc('text', usetex=True)

def get_step_lengths(fname):
    return data.SteppingData(fname).step_lengths

def get_springs_from_fname(fname):
    cbidx = fname.index("cb")
    cmidx = fname.index("cm")
    ctidx = fname.index("ct")
    cb = fname[cbidx+3:cbidx+3+fname[cbidx+3:].index(",")]
    cm = fname[cmidx+3:cmidx+3+fname[cmidx+3:].index(",")]
    ct = fname[ctidx+3:ctidx+3+fname[ctidx+3:].index(",")]
    return cb, cm, ct

def plot_steps(ct, cm, cb):
    yildiz_step_lengths = np.array([])
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-37]*1)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-35]*1)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-34]*1)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-33]*2)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-31]*2)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-30]*3)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-29]*1)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-28]*1)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-27]*4)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-26]*4)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-25]*2)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-24]*3)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-23]*4)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-21]*4)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-20]*3)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-19]*3)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-18]*5)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-17]*3)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-16]*3)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-15]*7)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-14]*5)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-13]*7)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-12]*12)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-11]*16)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-10]*14)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-9]*20)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-8]*14)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-7]*10)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-6]*9)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-5]*11)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-4]*8)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [-3]*2)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [4]*6)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [5]*7)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [6]*12)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [7]*20)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [8]*19)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [9]*22)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [10]*30)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [11]*34)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [12]*26)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [13]*21)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [14]*23)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [15]*22)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [16]*30)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [17]*29)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [18]*23)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [19]*22)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [20]*26)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [21]*12)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [22]*21)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [23]*16)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [24]*7)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [25]*8)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [26]*7)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [27]*8)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [28]*5)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [29]*9)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [30]*7)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [31]*8)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [32]*6)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [33]*2)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [34]*2)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [35]*9)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [36]*4)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [37]*9)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [38]*5)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [39]*1)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [40]*2)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [41]*1)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [42]*3)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [43]*2)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [44]*4)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [45]*4)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [46]*1)
    yildiz_step_lengths = np.append(yildiz_step_lengths, [47]*1)

    data_file = "data/stepping_data_randomsearch__k_b-1e+08,k_ub-1e+09,c--0.35,cb-{},cm-{},ct-{},ls-20.75,lt-23,seed-1,dt-1e-10.txt".format(cb, cm, ct)
    step_lengths = get_step_lengths(data_file)

    bins = np.histogram(np.hstack((yildiz_step_lengths)), bins=20)[1]

    plt.figure()
    plt.hist(yildiz_step_lengths, bins, alpha=0.5, label="Yildiz 2012", normed=True, stacked=True)
    plt.hist(step_lengths, bins, alpha=0.5, label="Model", normed=True, stacked=True)

    plt.legend(loc="upper right")
    plt.xlabel("Step length (nm)")
    plt.ylabel("Frequency")
    
    plt.savefig("plots/rs-{},{},{}.pdf".format(ct, cm, cb), format="pdf")


## get Yildiz stepping histogram

yildiz_step_lengths = np.array([])
yildiz_step_lengths = np.append(yildiz_step_lengths, [-37]*1)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-35]*1)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-34]*1)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-33]*2)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-31]*2)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-30]*3)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-29]*1)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-28]*1)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-27]*4)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-26]*4)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-25]*2)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-24]*3)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-23]*4)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-21]*4)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-20]*3)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-19]*3)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-18]*5)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-17]*3)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-16]*3)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-15]*7)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-14]*5)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-13]*7)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-12]*12)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-11]*16)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-10]*14)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-9]*20)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-8]*14)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-7]*10)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-6]*9)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-5]*11)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-4]*8)
yildiz_step_lengths = np.append(yildiz_step_lengths, [-3]*2)
yildiz_step_lengths = np.append(yildiz_step_lengths, [4]*6)
yildiz_step_lengths = np.append(yildiz_step_lengths, [5]*7)
yildiz_step_lengths = np.append(yildiz_step_lengths, [6]*12)
yildiz_step_lengths = np.append(yildiz_step_lengths, [7]*20)
yildiz_step_lengths = np.append(yildiz_step_lengths, [8]*19)
yildiz_step_lengths = np.append(yildiz_step_lengths, [9]*22)
yildiz_step_lengths = np.append(yildiz_step_lengths, [10]*30)
yildiz_step_lengths = np.append(yildiz_step_lengths, [11]*34)
yildiz_step_lengths = np.append(yildiz_step_lengths, [12]*26)
yildiz_step_lengths = np.append(yildiz_step_lengths, [13]*21)
yildiz_step_lengths = np.append(yildiz_step_lengths, [14]*23)
yildiz_step_lengths = np.append(yildiz_step_lengths, [15]*22)
yildiz_step_lengths = np.append(yildiz_step_lengths, [16]*30)
yildiz_step_lengths = np.append(yildiz_step_lengths, [17]*29)
yildiz_step_lengths = np.append(yildiz_step_lengths, [18]*23)
yildiz_step_lengths = np.append(yildiz_step_lengths, [19]*22)
yildiz_step_lengths = np.append(yildiz_step_lengths, [20]*26)
yildiz_step_lengths = np.append(yildiz_step_lengths, [21]*12)
yildiz_step_lengths = np.append(yildiz_step_lengths, [22]*21)
yildiz_step_lengths = np.append(yildiz_step_lengths, [23]*16)
yildiz_step_lengths = np.append(yildiz_step_lengths, [24]*7)
yildiz_step_lengths = np.append(yildiz_step_lengths, [25]*8)
yildiz_step_lengths = np.append(yildiz_step_lengths, [26]*7)
yildiz_step_lengths = np.append(yildiz_step_lengths, [27]*8)
yildiz_step_lengths = np.append(yildiz_step_lengths, [28]*5)
yildiz_step_lengths = np.append(yildiz_step_lengths, [29]*9)
yildiz_step_lengths = np.append(yildiz_step_lengths, [30]*7)
yildiz_step_lengths = np.append(yildiz_step_lengths, [31]*8)
yildiz_step_lengths = np.append(yildiz_step_lengths, [32]*6)
yildiz_step_lengths = np.append(yildiz_step_lengths, [33]*2)
yildiz_step_lengths = np.append(yildiz_step_lengths, [34]*2)
yildiz_step_lengths = np.append(yildiz_step_lengths, [35]*9)
yildiz_step_lengths = np.append(yildiz_step_lengths, [36]*4)
yildiz_step_lengths = np.append(yildiz_step_lengths, [37]*9)
yildiz_step_lengths = np.append(yildiz_step_lengths, [38]*5)
yildiz_step_lengths = np.append(yildiz_step_lengths, [39]*1)
yildiz_step_lengths = np.append(yildiz_step_lengths, [40]*2)
yildiz_step_lengths = np.append(yildiz_step_lengths, [41]*1)
yildiz_step_lengths = np.append(yildiz_step_lengths, [42]*3)
yildiz_step_lengths = np.append(yildiz_step_lengths, [43]*2)
yildiz_step_lengths = np.append(yildiz_step_lengths, [44]*4)
yildiz_step_lengths = np.append(yildiz_step_lengths, [45]*4)
yildiz_step_lengths = np.append(yildiz_step_lengths, [46]*1)
yildiz_step_lengths = np.append(yildiz_step_lengths, [47]*1)

bins = np.histogram(yildiz_step_lengths, bins=20)[1]

yildiz_fraction_per_bin = np.zeros(len(bins))

for b in range(len(bins)):
    counts = np.sum([1 for s in yildiz_step_lengths if s > bins[b] and s <= bins[b+1]])
    fraction = counts / len(yildiz_step_lengths)
    yildiz_fraction_per_bin[b] = fraction

data_files = []
for fname in os.listdir("data/"):
    if os.path.isfile("data/" + fname):
        if ("stepping_data_randomsearch__" in fname and ".txt" in fname):
            if ("~" not in fname and "movie" not in fname and "config" not in fname):
                data_files.append("data/" + fname)

cts = np.zeros(len(data_files))
cms = np.zeros(len(data_files))
cbs = np.zeros(len(data_files))
errors = np.zeros(len(data_files))

for df in range(len(data_files)):
    cb, cm, ct = get_springs_from_fname(data_files[df])
    step_lengths = get_step_lengths(data_files[df])
    error = 0
    cts[df] = float(ct)
    cbs[df] = float(cb)
    cms[df] = float(cm)
    for b in range(len(bins)-1):
        counts = np.sum([1 for s in step_lengths if s > bins[b] and s <= bins[b+1]])
        fraction = counts / len(step_lengths)
        error += (fraction-yildiz_fraction_per_bin[b])**2
    errors[df] = error

fig = plt.figure()
ax = plt.axes(projection='3d')
plot = ax.scatter(cbs, cms, cts, c=errors, cmap='viridis')
ax.set_xlabel("cbs")
ax.set_ylabel("cms")
ax.set_zlabel("cts")
plt.colorbar(plot)
plt.show()
plt.savefig("plots/randomsearchplot.png")

best_error_idxs = np.where(errors < np.sort(errors)[15])
for n, i in enumerate(best_error_idxs[0]):
    print("i: {}, ct: {}, cm: {}, cb: {}, error: {}".format(i, cts[i], cms[i], cbs[i], errors[i]))
    plot_steps(cts[i], cms[i], cbs[i])
