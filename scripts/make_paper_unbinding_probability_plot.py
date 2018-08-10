#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess, re

if 'show' not in sys.argv:
    matplotlib.use('Agg')

import argparse
import datetime
import dynein.draw.cartoon as cartoon
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle

import io

EPSILON = 1e-7

def equal(f1, f2):
    return abs(f1-f2) < EPSILON

data_files = []
for fname in os.listdir("data"):
    if os.path.isfile("data/" + fname):
        if ("paper_unbinding_probability__" in fname and ".txt" in fname):
            if ("~" not in fname):
                data_files.append("data/" + fname)

if len(data_files) == 0:
    print("No files of form data/*.txt found. Exiting.")
    exit(1)

Ls = []
mean_lagging_probability_per_L = []
mean_leading_probability_per_L = []
leading_probabilities_vs_time = []
lagging_probabilities_vs_time = []

for data_file in data_files:
    data = np.loadtxt(data_file, dtype = np.float64)
    if len(data) == 15:
        continue
    start_L_idx = data_file.find('L-')+2
    end_L_idx = data_file[start_L_idx:].find(',') + start_L_idx
    L = float(data_file[start_L_idx:end_L_idx])

    times = data[:,0]
    near_unbinding_probabilities = data[:,1]
    far_unbinding_probabilities = data[:,2]

    near_binding_xs = data[:,5]
    far_binding_xs = data[:,9]

    lagging_unbinding_probabilities = np.array([far_unbinding_probabilities[s] if near_binding_xs[s] > far_binding_xs[s] else near_unbinding_probabilities[s] for s in range(len(times))])
    leading_unbinding_probabilities = np.array([near_unbinding_probabilities[s] if near_binding_xs[s] > far_binding_xs[s] else far_unbinding_probabilities[s] for s in range(len(times))])
    
    Ls.append(L)
    mean_lagging_probability_per_L.append(np.mean(lagging_unbinding_probabilities))
    mean_leading_probability_per_L.append(np.mean(leading_unbinding_probabilities))
    leading_probabilities_vs_time.append(leading_unbinding_probabilities)
    lagging_probabilities_vs_time.append(lagging_unbinding_probabilities)

mean_lagging_probability_per_L = np.array(mean_lagging_probability_per_L)
mean_leading_probability_per_L = np.array(mean_leading_probability_per_L)
leading_probabilities_vs_time = np.array(leading_probabilities_vs_time)
lagging_probabilities_vs_time = np.array(lagging_probabilities_vs_time)

### lagging fraction vs L
fig = plt.figure()

yildiz_displacements = [10, 20, 30, 40, 50]
yildiz_lagging_fractions = [0.525, 0.545, 0.61, 0.59, 0.67]
yildiz_lagging_uncertainty = [0.06, 0.04, 0.035, 0.045, 0.075]

plt.errorbar(yildiz_displacements, yildiz_lagging_fractions, yerr=yildiz_lagging_uncertainty, label="Experimental (Yildiz 2012)", fmt='o-',)

plt.scatter(Ls, mean_lagging_probability_per_L / (mean_lagging_probability_per_L + mean_leading_probability_per_L), label="simulation")
plt.xlabel("|L| (nm)")
plt.ylabel("Lagging $k_{ub}$ fraction")
plt.legend()

plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("plots/paper_unbinding_probability_vs_L.pdf", format="pdf")
plt.close(fig)
