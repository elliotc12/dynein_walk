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

plt.rc('text', usetex=True)

EPSILON = 1e-7

def equal(f1, f2):
    return abs(f1-f2) < EPSILON

plt.rcParams.update({'font.size': 8})

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
        if (args.data_basename in fname and ".txt" in fname):
            if ("~" not in fname):
                data_files.append(args.data_directory + "/" + fname)

if len(data_files) == 0:
    print("No files of form " + args.data_directory + "/*" + args.data_basename + "*.txt found. Exiting.")
    exit(1)

step_times = []
onebound_times = []
bothbound_times = []
step_lengths = []
initial_displacements = []
step_times = []
run_velocities = []

alternating_passing = 0
alternating_not_passing = 0
not_alternating_passing = 0
not_alternating_not_passing = 0

leading_foot_steps = 0
trailing_foot_steps = 0

for data_file in data_files:
    data = np.loadtxt(data_file, dtype = np.float64)
    if len(data) < 3 or str(type(data[0])) == "<type 'numpy.float64'>":
        continue

    bind_times = np.array(data[:,1])
    unbind_times = np.array(data[:,0])
    near_foot_positions = np.around(np.array(data[:,2]), decimals=12)  #need to figure out why fixing number of decimals is necessary
    far_foot_positions = np.around(np.array(data[:,3]), decimals=12)
    near_step_lens = near_foot_positions[1:] - near_foot_positions[:-1]  #reduces total length by one. Will include 0 step lengths
    far_step_lens = far_foot_positions[1:] - far_foot_positions[:-1]

    for s in range(1,len(near_foot_positions)):
        assert(equal(near_foot_positions[s-1],near_foot_positions[s]) or equal(far_foot_positions[s-1],far_foot_positions[s]))
        if equal(near_foot_positions[s-1],near_foot_positions[s]) and equal(far_foot_positions[s-1],far_foot_positions[s]):
            continue

        if not equal(near_foot_positions[s-1],near_foot_positions[s]):  # i.e. if near foot stepped 
            step_lengths.append(data[s,2]-data[s-1,2])
            initial_displacements.append(data[s-1,2]-data[s-1,3])  #near_foot_positions[s-1]-far_foot_positions[s-1]
            #final_displacements.append(near_foot_positions[s]-far_foot_positions[s])
        elif not equal(far_foot_positions[s-1],far_foot_positions[s]):  # i.e. if far foot stepped 
            step_lengths.append(data[s,3]-data[s-1,3])
            initial_displacements.append(data[s-1,3]-data[s-1,2])  # far_foot_positions[s-1]-near_foot_positions[s-1]
            #final_displacements.append(far_foot_positions[s]-near_foot_positions[s]) 
        onebound_times = np.concatenate((onebound_times, [bind_times[s]-unbind_times[s]]))
        bothbound_times = np.concatenate((bothbound_times, [unbind_times[s]-bind_times[s-1]]))
        step_times = np.concatenate((step_times, onebound_times + bothbound_times))

    run_velocities.append((data[-1,2] + data[-1,2]) / 2 / data[-1,1])

    assert(len(near_foot_positions) > 5)
    for s in range(2, len(near_foot_positions)):
        if near_foot_positions[s-1] < far_foot_positions[s-1]:
            trailing_foot = near_foot_positions
            leading_foot = far_foot_positions
        else:
            trailing_foot = far_foot_positions
            leading_foot = near_foot_positions
        if not equal(near_foot_positions[s], near_foot_positions[s-1]) and not equal(far_foot_positions[s], far_foot_positions[s-1]):
            print("Error, both feet moved in a step.")
            exit(1)
        if not equal(trailing_foot[s], trailing_foot[s-1]): #must've been a leading foot step
            leading_foot_steps += 1
            if not equal(trailing_foot[s-1], trailing_foot[s-2]): # not alternating, the last was the same foot
                # the leading foot moved twice in a row, that makes this "not alternating"
                not_alternating_not_passing += 1
            else:
                # It is alternating because leading foot moved, but before that the other foot moved.
                alternating_not_passing += 1
        else: # must've been a trailing foot step
            trailing_foot_steps += 1
            if equal(trailing_foot[s-1], trailing_foot[s-2]): # alternating, other foot moved last time
                if trailing_foot[s] > leading_foot[s]:
                    not_alternating_passing += 1
                else:
                    not_alternating_not_passing += 1
            else:
                if trailing_foot[s] > leading_foot[s]:
                    alternating_passing += 1
                else:
                    alternating_not_passing += 1

num_steps = len(step_lengths)
initial_displacements = np.array(initial_displacements)

t_step = []
t_ob = []
t_bb = []
t_ob_uncertainty = []
t_bb_uncertainty = []
t_proc = []

t_step.append(np.mean(step_times))
t_ob.append(np.mean(onebound_times))
t_bb.append(np.mean(bothbound_times))
t_proc.append(t_step[-1]*t_bb[-1]/t_ob[-1])

t_ob_uncertainty.append(np.std(onebound_times)/np.sqrt(num_steps)*1.645) # 95% chance of true average being 1.645 stdevs from the sample average
t_bb_uncertainty.append(np.std(bothbound_times)/np.sqrt(num_steps)*1.645)

#####
## Model behavior plot
#####
fig = plt.figure(figsize=(4,4))
gs = gridspec.GridSpec(3, 2)
ax0 = fig.add_subplot(gs[0:2, 0:2])
ax1 = fig.add_subplot(gs[2, 0])
ax2 = fig.add_subplot(gs[2, 1], sharey=ax1)

#step length histogram

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

if len(step_lengths) == 0:
    print("No steps to put in histogram!")

bins = np.histogram(np.hstack((yildiz_step_lengths, step_lengths)), bins=20)[1]

ax0.hist(yildiz_step_lengths, bins, alpha=0.5, label="Yildiz 2012", normed=True, stacked=True)
ax0.hist(step_lengths, bins, alpha=0.5, label="Model", normed=True, stacked=True)

ax0.scatter([np.mean(step_lengths)], [0], label=r'$\overline{\Delta x} = ' + str(np.around(np.mean(step_lengths), decimals=2)) + r'$ \textit{nm}')

ax0.legend(loc="upper right")
ax0.set_xlabel("Step length (nm)")
ax0.set_ylabel("Frequency")

if args.parameters_filename != "":
    plt.gcf().suptitle(run_conditions + r' $k_{b}: \kb, k_{ub}: \kub, cb: \cb, cm: \cm, ct: \ct, runtime: \runtime$')

ax0.spines["top"].set_visible(False)
ax0.spines["right"].set_visible(False)
# plt.savefig("plots/stepping_length_histogram.pdf", format="pdf")

# OB time histogram
if len(step_times) > 0:
    ax1.hist(onebound_times, bins=np.logspace(np.log10(1e-9),np.log10(1e-3), 50))
else:
    print("Error, no step_times")
    exit(1)

ax1.set_ylabel("Frequency")
ax1.set_xlabel("Onebound time (s)")
ax1.set_xscale('log')
# ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

if args.parameters_filename != "":
    plt.gcf().suptitle(run_conditions + r' $k_{b}: \kb, k_{ub}: \kub, cb: \cb, cm: \cm, ct: \ct, runtime: \runtime$')

ob_theory_avg = 6e-6
ax1.axvline(ob_theory_avg, color='red', linestyle='dashed', linewidth=1)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
# ax1.legend()
# ax1.tight_layout()
# plt.savefig("plots/onebound_time_histogram.pdf", format="pdf")

# BB (dwell) time histogram
if len(step_times) > 0:
    ax2.hist(bothbound_times, bins=np.logspace(np.log10(1e-6),np.log10(1e-0), 50))

bb_theory_avg = 2.5e-3
ax2.axvline(bb_theory_avg, color='red', linestyle='dashed', linewidth=1)
# ax2.spines["left"].set_visible(False)
plt.setp([ax2.get_yticklabels()], visible=False)
ax2.tick_params(axis='y', which="both", left=False, right=False, labelbottom=False)
ax2.set_xlabel("Bothbound time (s)")
ax2.set_xscale('log')
# ax2.legend()
# ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

if args.parameters_filename != "":
    plt.gcf().suptitle(run_conditions + r' $k_{b}: \kb, k_{ub}: \kub, cb: \cb, cm: \cm, ct: \ct, runtime: \runtime$')

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.tight_layout(w_pad=2, h_pad=0.5)
plt.savefig("plots/model_behavior.pdf", format="pdf")

#step time histogram
gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, 1, 1])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])
ax3 = fig.add_subplot(gs[3])

if len(step_times) > 0:
    ax0.hist(step_times, bins=np.logspace(np.log10(1e-10),np.log10(1e-2), 50))
    ax1.hist(onebound_times, bins=np.logspace(np.log10(1e-10),np.log10(1e-2), 50))
    ax2.hist(bothbound_times, bins=np.logspace(np.log10(1e-10),np.log10(1e-2), 50))
    ax3.hist(run_velocities, bins=50)

ax0.set_ylabel("Frequency")
ax0.set_xscale('log')
# ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax1.set_ylabel("Frequency")
ax1.set_xscale('log')
# ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax2.set_xlabel("Step time (s)")
ax2.set_ylabel("Frequency")
ax2.set_xscale('log')
# ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax3.set_xlabel("velocity (nm/s)")
ax3.set_ylabel("Frequency")

if args.parameters_filename != "":
    plt.gcf().suptitle(run_conditions + r' $k_{b}: \kb, k_{ub}: \kub, cb: \cb, cm: \cm, ct: \ct, runtime: \runtime$')

plt.subplots_adjust(hspace=0.6)

plt.show()
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
# plt.tight_layout()
plt.savefig("plots/stepping_time_histogram.pdf", format="pdf")
plt.close(fig)

# OB_time vs step_length scatter
assert len(onebound_times) == len(step_lengths)
fig = plt.figure()
plt.scatter(onebound_times, step_lengths, alpha=0.1)
plt.gca().set_xscale('log')
plt.xlabel("Onebound time (s)")
plt.ylabel("Step length (nm)")

plt.gca().set_xlim((1e-7, 1e-2))

if args.parameters_filename != "":
    plt.gcf().suptitle(run_conditions + r' $k_{b}: \kb, k_{ub}: \kub, cb: \cb, cm: \cm, ct: \ct, runtime: \runtime$')

plt.tight_layout()
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.savefig("plots/ob-vs-length-scatter.pdf", format="pdf")
plt.close(fig)

# BB_time vs step_length scatter
assert len(bothbound_times) == len(step_lengths)
fig = plt.figure()
plt.scatter(bothbound_times, step_lengths, alpha=0.1)
plt.gca().set_xscale('log')
plt.xlabel("Bothbound time (s), theory= 5e-3s")
plt.ylabel("Step length (nm)")

plt.gca().set_xlim((1e-5, 1))

if args.parameters_filename != "":
    plt.gcf().suptitle(run_conditions + r' $k_{b}: \kb, k_{ub}: \kub, cb: \cb, cm: \cm, ct: \ct, runtime: \runtime$')

plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("plots/bb-vs-length-scatter.pdf", format="pdf")
plt.close(fig)

# initial displacement vs motor step length scatter
fig = plt.figure()
plt.plot(initial_displacements, initial_displacements+step_lengths, '.', alpha=0.3)
plt.xlabel("Initial displacement (nm)")
plt.ylabel("Step length (nm)")
plt.ylabel("Final displacement (nm)")
plt.axes().set_aspect('equal')

if args.parameters_filename != "":
    plt.gcf().suptitle(run_conditions + r' $k_{b}: \kb, k_{ub}: \kub, cb: \cb, cm: \cm, ct: \ct, runtime: \runtime$')

plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("plots/displacement_vs_step_length.pdf", format="pdf")
plt.close(fig)

fig = plt.figure(figsize=(8*.6, 6*.6), dpi=300)
plt.rc('text', usetex=True)

gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])

width = 0.5

ax1.bar([0, 1], [alternating_passing + alternating_not_passing, not_alternating_passing + not_alternating_not_passing,], width)
ax1.set_ylabel('frequency')
ax1.set_xticks([0, 1])
ax1.set_xticklabels(('alternating', 'not\nalternating'), rotation=45)

ax2.bar([0, 1], [alternating_passing + not_alternating_passing, alternating_not_passing + not_alternating_not_passing,], width)
ax2.set_ylabel('frequency')
ax2.set_xticks([0, 1])
ax2.set_xticklabels(('passing', 'not passing'), rotation=45)

ax3.bar([0, 1, 2, 3], [alternating_passing, alternating_not_passing, not_alternating_passing, not_alternating_not_passing,], width)
ax3.set_ylabel('frequency')
ax3.set_xticks([0, 1, 2, 3])
ax3.set_xticklabels(('alternating,\n passing', 'alternating,\n not passing', 'not alternating,\n passing', 'not alternating,\n not passing'), rotation=45)

plt.tight_layout()
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.savefig("plots/stepping_analysis.pdf", format="pdf")
plt.close(fig)


# fig = plt.figure()

# initial_displacements = np.array(initial_displacements)
# indices = np.argsort(np.abs(initial_displacements))
# sorted_displacements = initial_displacements[indices]

# Nbin = 50
# L = np.zeros(int((len(sorted_displacements)-1)/Nbin)+1)
# ntrailing = np.zeros_like(L)
# nleading = np.zeros_like(L)
# for i in range(len(L)):
#     bunch = sorted_displacements[i*Nbin:(i+1)*Nbin]
#     ntrailing[i] = (bunch < 0).sum()
#     nleading[i] = (bunch > 0).sum()
#     L[i] = np.abs(bunch).mean()

# fraction_trailing = ntrailing / (ntrailing + nleading)

# yildiz_displacements = [10, 20, 30, 40, 50]
# yildiz_fractions = [0.525, 0.545, 0.61, 0.59, 0.67]
# yildiz_uncertainty = [0.06, 0.04, 0.035, 0.045, 0.075]

# plt.plot(L, fraction_trailing, 'o-', label="Model")
# plt.errorbar(yildiz_displacements, yildiz_fractions, yerr=yildiz_uncertainty, label="Experimental (Yildiz 2012)", fmt='o-',)
# plt.ylim(0,1)

# plt.xlabel("FIXME Initial foot x-displacement (unstepping - stepping) (nm)")
# plt.ylabel("Fraction trailing")

# plt.ylim([0,1.1])

# plt.legend()

# if args.parameters_filename != "":
#     plt.gcf().suptitle(run_conditions + r' $k_{b}: \kb, k_{ub}: \kub, cb: \cb, cm: \cm, ct: \ct, runtime: \runtime$')

# plt.gca().spines["top"].set_visible(False)
# plt.gca().spines["right"].set_visible(False)
# plt.tight_layout()
# plt.savefig("plots/displacement_histogram.pdf", format="pdf")
# plt.close(fig)

plt.figure()

data = np.loadtxt(data_files[0], dtype = np.float64)
bind_times = np.array(data[:,1])
near_foot_positions = np.around(np.array(data[:,2]), decimals=12)
far_foot_positions = np.around(np.array(data[:,3]), decimals=12)

bind_times_duplicated = np.zeros(len(bind_times)*2-1)
near_positions_duplicated = np.zeros(len(bind_times)*2-1)
far_positions_duplicated = np.zeros(len(bind_times)*2-1)

for t in range(1, len(bind_times)):
    bind_times_duplicated[2*t-1] = bind_times[t]
    bind_times_duplicated[2*t] = bind_times[t]
    near_positions_duplicated[2*t-1] = near_foot_positions[t-1]
    near_positions_duplicated[2*t] = near_foot_positions[t]
    far_positions_duplicated[2*t-1] = far_foot_positions[t-1]
    far_positions_duplicated[2*t] = far_foot_positions[t]

bind_times_duplicated[0] = bind_times[0]
near_positions_duplicated[0] = near_foot_positions[0]
far_positions_duplicated[0] = far_foot_positions[0]

plt.plot(bind_times_duplicated, near_positions_duplicated, 'o-', label="Near BD", markersize=.01, linewidth=0.5)
plt.plot(bind_times_duplicated, far_positions_duplicated, 'o-', label="Far BD", markersize=.01, linewidth=0.5)

plt.xlabel("time (s)")
plt.ylabel("Position (nm)")

plt.legend()

plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("plots/stepping_trajectory.pdf", format="pdf")
plt.close(fig)
