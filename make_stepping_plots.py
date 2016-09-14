#!/usr/bin/python2.7

import matplotlib
matplotlib.use("PDF") # use PDF generating backend bc it can handle bigger datasets
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

def add_config_info(config_txt):
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
    box = TextArea(config_txt, textprops=dict(size=10))
    anchored_box = AnchoredOffsetbox(loc=2, child=box, bbox_to_anchor=(1,1),
                                bbox_transform=ax.transAxes, frameon=False)
    ax.add_artist(anchored_box)

data_name = sys.argv[1]

data_arr = np.loadtxt("data/stepping_data_" + data_name + ".txt", delimiter=",")
config_txt = open("data/stepping_config_" + data_name + ".txt", "r").read()

num_steps = data_arr[0]
step_lengths = data_arr[1:num_steps+1]
step_times = data_arr[num_steps+1:]

#step length histogram
fig = plt.figure()
plt.hist(step_lengths,
   bins=[-64, -56, -48, -40, -32, -24, -16, -8, 0, 8, 16, 24, 32, 40, 48, 56, 64])
plt.title("Stepping length histogram " + data_name)
plt.xlabel("Step length (nm)")
plt.ylabel("Frequency")
add_config_info(config_txt)
plt.savefig("plots/stepping_length_histogram_" + data_name + ".pdf", format="pdf")
plt.close(fig)

#step time histogram
fig = plt.figure()
plt.hist(step_times, bins=12)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title("Stepping time histogram " + data_name)
plt.xlabel("Step time (s)")
plt.ylabel("Frequency")
add_config_info(config_txt)
plt.savefig("plots/stepping_time_histogram_" + data_name + ".pdf", format="pdf")
plt.close(fig)
