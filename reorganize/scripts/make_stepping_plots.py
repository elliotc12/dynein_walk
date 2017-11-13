#!/usr/bin/python3
from __future__ import division, print_function

import matplotlib
matplotlib.use("PDF") # use PDF generating backend bc it can handle bigger datasets
import matplotlib.pyplot as plt
import os
import numpy as np
import sys

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

file_string = sys.argv[1]
data_files = []

for fname in os.listdir("data/"):
    if os.path.isfile("data/" + fname):
        if (file_string in "data/" + fname) and ("stepping_data" in "data/" + fname):
            data_files.append("data/" + fname)

step_times = []
step_lengths = []
onebound_times = []
bothbound_times = []

for data_file in data_files:
    data = np.loadtxt(data_file)
    if len(data) < 3 or str(type(data[0])) == "<type 'numpy.float64'>":
        continue

    bind_times = np.array(data[:,1])
    unbind_times = np.array(data[:,0])
    near_positions = -1*np.array(data[:,2])
    far_positions = -1*np.array(data[:,3])
    near_step_idxs = near_positions[1:] != near_positions[:-1]
    far_step_idxs = far_positions[1:] != far_positions[:-1]
    near_step_lens = (near_positions[1:] - near_positions[:-1])[near_step_idxs]
    far_step_lens = (far_positions[1:] - far_positions[:-1])[far_step_idxs]

    step_times = np.concatenate((step_times, np.array(bind_times[1:] - bind_times[:-1])))
    onebound_times = np.concatenate((onebound_times, bind_times - unbind_times))
    bothbound_times = np.concatenate((bothbound_times, bind_times[1:] - unbind_times[:-1]))
    step_lengths = np.concatenate((step_lengths, near_step_lens, far_step_lens))

#step length histogram
fig = plt.figure()
plt.hist(step_lengths, bins=50)
plt.xlabel("Step length (nm)")
plt.ylabel("Frequency")
plt.savefig("plots/stepping_length_histogram.pdf", format="pdf")
plt.close(fig)

#step time histogram
fig = plt.figure()
plt.hist(step_times, bins=50)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel("Step time (s)")
plt.ylabel("Frequency")
plt.savefig("plots/stepping_time_histogram.pdf", format="pdf")
plt.close(fig)
