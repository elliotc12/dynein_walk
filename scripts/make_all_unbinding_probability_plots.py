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

parser = argparse.ArgumentParser(description = 'script to generate various histograms from stepping data.')

parser.add_argument('-d', '--data-directory', dest = 'data_directory', action='store', type = str,
                    default="", help='data file directory', required = True)
parser.add_argument('-b', '--data-basename', dest = 'data_basename', action='store', type = str,
                    default="", help='data file basename', required = True)
parser.add_argument('-p', '--param-file', dest = 'parameters_filename', action='store', type = str,
                    default="", help='parameter filename (.tex)')

args = parser.parse_args()

if args.parameters_filename != "":
    run_conditions = open(args.parameters_filename).read().replace("\n", " ").replace("\\\\", "\\")

data_files = []
for fname in os.listdir(args.data_directory):
    if os.path.isfile(args.data_directory + "/" + fname):
        if (args.data_basename + "__" in fname and ".txt" in fname):
            if ("~" not in fname):
                data_files.append(args.data_directory + "/" + fname)

if len(data_files) == 0:
    print("No files of form " + args.data_directory + "/*" + args.data_basename + "*.txt found. Exiting.")
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
yildiz_fractions = [0.525, 0.545, 0.61, 0.59, 0.67]
yildiz_uncertainty = [0.06, 0.04, 0.035, 0.045, 0.075]

plt.errorbar(yildiz_displacements, yildiz_fractions, yerr=yildiz_uncertainty, label="Yildiz 2012", fmt='o-',)

plt.scatter(Ls, mean_lagging_probability_per_L / (mean_lagging_probability_per_L + mean_leading_probability_per_L), label="simulation")
plt.xlabel("|L| (nm)")
plt.ylabel("Lagging $k_{ub}$ fraction")
plt.legend()

plt.savefig("plots/unbinding_probability/plots_for_latex/lagging_fraction_vs_L.pdf", format="pdf")
plt.close(fig)

### <Probability> vs L
fig = plt.figure()

plt.scatter(Ls, mean_lagging_probability_per_L, label="lagging")
plt.scatter(Ls, mean_leading_probability_per_L, label="leading")
plt.xlabel("|L| (nm)")
plt.ylabel("<$k_{ub}$>")
plt.legend()

plt.savefig("plots/unbinding_probability/plots_for_latex/unbinding_rate_vs_L.pdf", format="pdf")
plt.close(fig)

### Probability vs t
fig = plt.figure()

# for s in range(len(Ls)):
#     plt.plot(times, leading_probabilities_vs_time[s], label = "leading, L = " + str(Ls[s]), alpha=0.1)
#     plt.plot(times, lagging_probabilities_vs_time[s], label = "lagging, L = " + str(Ls[s]), alpha=0.1)

plt.plot(times, leading_probabilities_vs_time[0], label = "leading, L = " + str(Ls[0]))

plt.legend()
plt.xlabel("t (s)")
plt.ylabel("$k_{ub}$")

plt.savefig("plots/unbinding_probability/plots_for_latex/rate_vs_t.pdf", format="pdf")
plt.close(fig)
