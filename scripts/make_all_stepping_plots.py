#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess, re
from io import BytesIO

if 'show' not in sys.argv:
    matplotlib.use('Agg')

import argparse
import datetime
import dynein.draw.cartoon as cartoon
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle

from basic_units import radians

import dynein.data as datalib

import io

plt.rc('text', usetex=True)

EPSILON = 1e-7

def equal(f1, f2):
    return abs(f1-f2) < EPSILON

def get_stepping_data(args):
    data_files = []
    for fname in os.listdir(args.data_directory):
        if os.path.isfile(args.data_directory + "/" + fname):
            if (args.data_basename in fname and ".txt" in fname):
                if ("~" not in fname and "movie" not in fname and "config" not in fname):
                    data_files.append(args.data_directory + "/" + fname)

    if len(data_files) == 0:
        print("No files of form " + args.data_directory + "/*" + args.data_basename + "*.txt found. Exiting.")
        exit(1)

    data_objects = []

    for data_file in data_files:
        data_objects.append(datalib.SteppingData(data_file))

    stepping_data = {}

    stepping_data["onebound_times"] = np.concatenate([d.onebound_times for d in data_objects])
    stepping_data["bothbound_times"] = np.concatenate([d.bothbound_times for d in data_objects])
    stepping_data["step_times"] = stepping_data["onebound_times"] + stepping_data["bothbound_times"]
    stepping_data["step_lengths"] = np.concatenate([d.step_lengths for d in data_objects])
    stepping_data["initial_displacements"] = np.concatenate([d.initial_displacements for d in data_objects])

    stepping_data["alternating_passing"] = sum([d.alternating_passing for d in data_objects])
    stepping_data["alternating_not_passing"] = sum([d.alternating_not_passing for d in data_objects])
    stepping_data["not_alternating_passing"] = sum([d.not_alternating_passing for d in data_objects])
    stepping_data["not_alternating_not_passing"] = sum([d.not_alternating_not_passing for d in data_objects])

    stepping_data["leading_foot_steps"] = np.sum([d.leading_foot_steps for d in data_objects])
    stepping_data["trailing_foot_steps"] = np.sum([d.trailing_foot_steps for d in data_objects])

    stepping_data["num_steps"] = len(stepping_data["step_lengths"])
    stepping_data["initial_displacements"] = np.array(stepping_data["initial_displacements"])
    return stepping_data

def get_unbinding_probability_data(args):
    data_files = []
    for fname in os.listdir("data"):
        if os.path.isfile("data/" + fname):
            if ("paper_unbinding_probability__" in fname and ".txt" in fname):
                if ("~" not in fname):
                    data_files.append("data/" + fname)

    if len(data_files) == 0:
        print("No files of form data/*.txt found. Exiting.")
        exit(1)

    up_data = {}
    up_data["Ls"] = []
    up_data["mean_lagging_probability_per_L"] = []
    up_data["mean_leading_probability_per_L"] = []

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

        up_data["Ls"].append(L)
        up_data["mean_lagging_probability_per_L"].append(np.mean(lagging_unbinding_probabilities))
        up_data["mean_leading_probability_per_L"].append(np.mean(leading_unbinding_probabilities))

    up_data["mean_lagging_probability_per_L"] = np.array(up_data["mean_lagging_probability_per_L"])
    up_data["mean_leading_probability_per_L"] = np.array(up_data["mean_leading_probability_per_L"])
    return up_data

def get_force_data(args):
    data_files = []
    for fname in os.listdir("data"):
        if os.path.isfile("data/" + fname):
            if ("paper_force_stepping_data" in fname and ".txt" in fname):
                if ("~" not in fname):
                    data_files.append("data/" + fname)

    if len(data_files) == 0:
        print("No files of form data/stepping_data_paper_force__*.txt found. Exiting.")
        exit(1)

    force_data = {}
    force_data["forces"] = []
    force_data["velocities"] = []
    force_data["seeds"] = []

    for data_file in data_files:
        data = np.loadtxt(data_file, dtype = np.float64)
        if len(data) == 0 or len(data[1, :]) <= 5:
            continue
        start_F_idx = data_file.find('F-')+2
        end_F_idx = data_file[start_F_idx:].find(',') + start_F_idx
        F = float(data_file[start_F_idx:end_F_idx])

        start_s_idx = data_file.find(',s-')+3
        end_s_idx = data_file[start_s_idx:].find('.') + start_s_idx
        s = float(data_file[start_s_idx:end_s_idx])

        final_nbx = float(data[-1,2])
        final_fbx = float(data[-1,3])
        final_t   = float(data[-1,1])
        velocity = (final_nbx + final_fbx) / (2 * final_t)

        force_data["forces"].append([F])
        force_data["seeds"].append([s])
        force_data["velocities"].append([velocity])
    return force_data

def get_stroke_angles_data(args, longrun=False):
    if longrun:
        basename = "paper_long_stroke_angle_data"
    else:
        basename = "paper_stroke_angle_data"

    data_files = []
    for fname in os.listdir("data"):
        if os.path.isfile("data/" + fname):
            if (basename in fname and ".txt" in fname):
                if ("~" not in fname):
                    data_files.append("data/" + fname)

    if len(data_files) == 0:
        print("Please run 'python3 scripts/generate-angle-data.py' to generate angle data. Replacing angle figures with dummies.")
        return "NODATAFILE"

    data_txt = ""
    for data_file in data_files:
        with open(data_file, 'r') as fd:
            data_txt = data_txt + fd.read()

    onebound_angles_data = {}
    onebound_angles_data["state"] = []
    onebound_angles_data["times"] = []
    onebound_angles_data["dtailx"] = []
    onebound_angles_data["taily"] = []
    onebound_angles_data["ma"] = []
    onebound_angles_data["d_ma"] = []

    bothbound_angles_data = {}
    bothbound_angles_data["state"] = []
    bothbound_angles_data["times"] = []
    bothbound_angles_data["dtailx"] = []
    bothbound_angles_data["taily"] = []
    bothbound_angles_data["ma"] = []
    bothbound_angles_data["d_ma"] = []
    bothbound_angles_data["ba"] = []

    onebound_longest_run_length = 0
    onebound_longest_run_times = []

    bothbound_longest_run_length = 0
    bothbound_longest_run_times = []

    onebound_data_txt = ""
    bothbound_data_txt = ""

    for line in data_txt.split('\n'):
        if "NEARBOUN" in line or "FARBOUND" in line or "NEWUNBINDING" in line:
            onebound_data_txt += line + '\n'
        elif '#' not in line:
            bothbound_data_txt += line + '\n'

    while (onebound_data_txt.find("NEWUNBINDING") != -1):
        step_idx = onebound_data_txt.find("NEWUNBINDING")
        if (onebound_data_txt[step_idx+13:].find("NEWUNBINDING") != -1):
            next_step_idx = onebound_data_txt[step_idx+13:].find("NEWUNBINDING")
            stroke_onebound_data_txt = onebound_data_txt[step_idx+13:step_idx+13+next_step_idx]
            onebound_data_txt = onebound_data_txt[step_idx+13+next_step_idx:]
        else:
            stroke_onebound_data_txt = onebound_data_txt[step_idx+13:]
            onebound_data_txt = ""

        if (len(stroke_onebound_data_txt) == 0):
            continue

        state = stroke_onebound_data_txt[0:9]
        stroke_data = np.genfromtxt(BytesIO(stroke_onebound_data_txt.encode()))

        if (stroke_data.size < 7):
            continue

        onebound_angles_data["state"].append(state)
        onebound_angles_data["times"].append(stroke_data[:,1] - stroke_data[0,1])
        onebound_angles_data["dtailx"].append(stroke_data[:,2] - stroke_data[0,2])
        onebound_angles_data["taily"].append(stroke_data[:,3])
        onebound_angles_data["ma"].append(stroke_data[:,4])
        onebound_angles_data["d_ma"].append(stroke_data[:,4] - stroke_data[0,4])
        if len(stroke_data[:,1]) > onebound_longest_run_length:
            onebound_longest_run_length = len(stroke_data[:,1])
            onebound_longest_run_times = stroke_data[:,1] - stroke_data[0,1]

    onebound_angles_data["times"] = np.array([np.pad(s, (0,onebound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in onebound_angles_data["times"]]) / 1e-9
    onebound_angles_data["dtailx"] = np.array([np.pad(s, (0,onebound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in onebound_angles_data["dtailx"]])
    onebound_angles_data["taily"] = np.array([np.pad(s, (0,onebound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in onebound_angles_data["taily"]])
    onebound_angles_data["ma"] = np.array([np.pad(s, (0,onebound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in onebound_angles_data["ma"]])
    onebound_angles_data["d_ma"] = np.array([np.pad(s, (0,onebound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in onebound_angles_data["d_ma"]])
    onebound_angles_data["longest_times"] = onebound_longest_run_times / 1e-9

    while (bothbound_data_txt.find("NEWBINDING") != -1):
        step_idx = bothbound_data_txt.find("NEWBINDING")
        if (bothbound_data_txt[step_idx+10:].find("NEWBINDING") != -1):
            next_step_idx = bothbound_data_txt[step_idx+10:].find("NEWBINDING")
            stroke_bothbound_data_txt = bothbound_data_txt[step_idx+10:step_idx+10+next_step_idx]
            bothbound_data_txt = bothbound_data_txt[step_idx+10+next_step_idx:]
        else:
            stroke_bothbound_data_txt = bothbound_data_txt[step_idx+10:]
            bothbound_data_txt = ""

        if (len(stroke_bothbound_data_txt) < 5):
            continue

        state = stroke_bothbound_data_txt[0:9]
        stroke_data = np.genfromtxt(BytesIO(stroke_bothbound_data_txt.encode()))

        if (stroke_data.size < 7):
            continue

        bothbound_angles_data["state"].append(state)
        bothbound_angles_data["times"].append(stroke_data[:,1] - stroke_data[0,1])
        bothbound_angles_data["dtailx"].append(stroke_data[:,2] - stroke_data[0,2])
        bothbound_angles_data["taily"].append(stroke_data[:,3])
        bothbound_angles_data["ma"].append(stroke_data[:,4])
        bothbound_angles_data["d_ma"].append(stroke_data[:,4] - stroke_data[0,4])
        bothbound_angles_data["ba"].append(stroke_data[:,5])
        if len(stroke_data[:,1]) > bothbound_longest_run_length:
            bothbound_longest_run_length = len(stroke_data[:,1])
            bothbound_longest_run_times = stroke_data[:,1] - stroke_data[0,1]

    bothbound_angles_data["times"] = np.array([np.pad(s, (0,bothbound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in bothbound_angles_data["times"]]) / 1e-9
    bothbound_angles_data["dtailx"] = np.array([np.pad(s, (0,bothbound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in bothbound_angles_data["dtailx"]])
    bothbound_angles_data["taily"] = np.array([np.pad(s, (0,bothbound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in bothbound_angles_data["taily"]])
    bothbound_angles_data["ma"] = np.array([np.pad(s, (0,bothbound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in bothbound_angles_data["ma"]])
    bothbound_angles_data["d_ma"] = np.array([np.pad(s, (0,bothbound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in bothbound_angles_data["d_ma"]])
    bothbound_angles_data["ba"] = np.array([np.pad(s, (0,bothbound_longest_run_length-len(s)), mode="constant", constant_values=np.nan) for s in bothbound_angles_data["ba"]])
    bothbound_angles_data["longest_times"] = bothbound_longest_run_times / 1e-9

    return (onebound_angles_data, bothbound_angles_data)

def get_cli_arguments():
    parser = argparse.ArgumentParser(description = 'script to generate various histograms from stepping data.')

    parser.add_argument('-d', '--data-directory', dest = 'data_directory', action='store', type = str,
                        default="data/", help='data file directory', required = False)
    parser.add_argument('-b', '--data-basename', dest = 'data_basename', action='store', type = str,
                        default="", help='data file basename', required = True)
    parser.add_argument('-p', '--param-file', dest = 'parameters_filename', action='store', type = str,
                        default="", help='parameter filename (.tex)')
    parser.add_argument('-q', '--quick', dest='quick', action='store_true', default=False,
                        help='Only make the stepping length histogram')

    return parser.parse_args()

def make_behavior_plot(args, stepping_data, up_data):
    fig = plt.figure(figsize=(4,6))
    gs = gridspec.GridSpec(5, 2)
    ax0 = fig.add_subplot(gs[0:2, 0:2])
    ax1 = fig.add_subplot(gs[2:4, 0:2])
    ax2 = fig.add_subplot(gs[4, 0])
    ax3 = fig.add_subplot(gs[4, 1], sharey=ax2)

    data_files = []
    for fname in os.listdir(args.data_directory):
        if os.path.isfile(args.data_directory + "/" + fname):
            if (args.data_basename in fname and ".txt" in fname):
                if ("~" not in fname and "movie" not in fname and "config" not in fname):
                    data_files.append(args.data_directory + "/" + fname)

    if len(data_files) == 0:
        print("No files of form " + args.data_directory + "/*" + args.data_basename + "*.txt found. Exiting.")
        exit(1)

    yildiz_data = np.loadtxt("data/yildiz_tracking_data.txt", dtype = np.float64)
    ax0.plot(yildiz_data[:,0]/1e3, yildiz_data[:,1], 'o-', label="Experiment", markersize=0, linewidth=0.5, color="C0", alpha=0.8)

    # model position trace
    for i in [1, 2, 3, 4, 5, 6, 7, 8]:
        df = "data/paper_main_stepping_data-{}.txt".format(i)
        if df not in data_files:
            continue
        data = np.loadtxt(df, dtype = np.float64)
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

        if (i == 1):
            ax0.plot(bind_times, near_foot_positions, 'o-', label="Model", markersize=0, linewidth=0.5, c="C1", alpha=0.8)
        else:
            ax0.plot(bind_times, near_foot_positions, 'o-', markersize=0, linewidth=0.5, c="C1", alpha=0.8)

    ax0.set_xlabel("time (s)")
    ax0.set_ylabel("Foot $\hat{x}$ position (nm)")

    ax0.legend()

    ax0.spines["top"].set_visible(False)
    ax0.spines["right"].set_visible(False)

    #step length histogram
    yildiz_step_lengths = np.concatenate(([-37]*1, [-35]*1, [-34]*1, [-33]*2, [-31]*2, [-30]*3, [-29]*1, [-28]*1, [-27]*4, [-26]*4, [-25]*2, [-24]*3, [-23]*4, [-21]*4, [-20]*3,
                                          [-19]*3, [-18]*5, [-17]*3, [-16]*3, [-15]*7, [-14]*5, [-13]*7, [-12]*12, [-11]*16, [-10]*14, [-9]*20, [-8]*14, [-7]*10, [-6]*9, [-5]*11,
                                          [-4]*8, [-3]*2, [4]*6, [5]*7, [6]*12, [7]*20, [8]*19, [9]*22, [10]*30, [11]*34, [12]*26, [13]*21, [14]*23, [15]*22, [16]*30, [17]*29,
                                          [18]*23, [19]*22, [20]*26, [21]*12, [22]*21, [23]*16, [24]*7, [25]*8, [26]*7, [27]*8, [28]*5, [29]*9, [30]*7, [31]*8, [32]*6, [33]*2,
                                          [34]*2, [35]*9, [36]*4, [37]*9, [38]*5, [39]*1, [40]*2, [41]*1, [42]*3, [43]*2, [44]*4, [45]*4, [46]*1, [47]*1))

    if len(stepping_data["step_lengths"]) == 0:
        print("No steps to put in histogram!")

    bins = np.histogram(np.hstack((yildiz_step_lengths, stepping_data["step_lengths"])), bins=20)[1]

    ax1.hist(yildiz_step_lengths, bins, alpha=0.5, label="Experiment", normed=True, stacked=True, color="C0")
    ax1.hist(stepping_data["step_lengths"], bins, alpha=0.5, label="Model", normed=True, stacked=True, color="C1")

    # ax1.scatter([np.mean(stepping_data["step_lengths"])], [0], label=r'$\overline{\Delta x} = ' + str(np.around(np.mean(stepping_data["step_lengths"]), decimals=2)) + r'$ \textit{nm}')

    ax1.legend(loc="upper right")
    ax1.set_xlabel("Step length (nm)")
    ax1.set_ylabel("Frequency")

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # OB time histogram
    if len(stepping_data["step_times"]) > 0:
        ax2.hist(stepping_data["onebound_times"], bins=np.logspace(np.log10(1e-9),np.log10(1e-3), 50))
    else:
        print("Error, no step_times")
        exit(1)

    ax2.set_ylabel("Frequency")
    ax2.set_xlabel("Onebound time (s)")
    ax2.set_xscale('log')

    ob_theory_avg = 5.2e-4
    ax2.axvline(ob_theory_avg, color='red', linestyle='dashed', linewidth=1)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # BB (dwell) time histogram
    if len(stepping_data["step_times"]) > 0:
        ax3.hist(stepping_data["bothbound_times"], bins=np.logspace(np.log10(1e-6),np.log10(1e-0), 50))

    bb_theory_avg = 6.4e-2
    ax3.axvline(bb_theory_avg, color='red', linestyle='dashed', linewidth=1)
    plt.setp([ax3.get_yticklabels()], visible=False)
    ax3.tick_params(axis='y', which="both", left=False, right=False, labelbottom=False)
    ax3.set_xlabel("Bothbound time (s)")
    ax3.set_xscale('log')

    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)

    plt.tight_layout(w_pad=2, h_pad=0.5)
    plt.savefig("plots/model_behavior.pdf", format="pdf")

    # unbinding probability vs L plot
    fig = plt.figure(figsize=(4,2))
    yildiz_displacements = [10, 20, 30, 40, 50]
    yildiz_lagging_fractions = [0.525, 0.545, 0.61, 0.59, 0.67]
    yildiz_lagging_uncertainty = [0.06, 0.04, 0.035, 0.045, 0.075]

    plt.errorbar(yildiz_displacements, yildiz_lagging_fractions, yerr=yildiz_lagging_uncertainty, label="Experiment", fmt='o-', c='C0', markersize=4, linestyle='', capsize=1, elinewidth=0.3, markeredgewidth=0.3)

    plt.scatter(up_data["Ls"], up_data["mean_lagging_probability_per_L"] / (up_data["mean_lagging_probability_per_L"] + up_data["mean_leading_probability_per_L"]), c='C1', label="Model", zorder=2, s=4**2)
    plt.gca().set_xlabel("Binding domain separation (nm)")
    fig.gca().set_ylabel("P(lagging step)")
    plt.legend()

    fig.gca().spines["top"].set_visible(False)
    fig.gca().spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig("plots/unbinding_probability_vs_displacement.pdf")

def make_analysis_plot(args, stepping_data):
    fig = plt.figure(figsize=(8*.6, 6*.6), dpi=300)
    plt.rc('text', usetex=True)

    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    width = 0.5

    ax1.bar([0, 1], [stepping_data["alternating_passing"] + stepping_data["alternating_not_passing"], stepping_data["not_alternating_passing"] + stepping_data["not_alternating_not_passing"],], width)
    ax1.set_ylabel('frequency')
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(('alternating', 'not\nalternating'), rotation=45)

    ax2.bar([0, 1], [stepping_data["alternating_passing"] + stepping_data["not_alternating_passing"], stepping_data["alternating_not_passing"] + stepping_data["not_alternating_not_passing"],], width)
    ax2.set_ylabel('frequency')
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(('passing', 'not passing'), rotation=45)

    ax3.bar([0, 1, 2, 3], [stepping_data["alternating_passing"], stepping_data["alternating_not_passing"], stepping_data["not_alternating_passing"], stepping_data["not_alternating_not_passing"],], width)
    ax3.set_ylabel('frequency')
    ax3.set_xticks([0, 1, 2, 3])
    ax3.set_xticklabels(('alternating,\n passing', 'alternating,\n not passing', 'not alternating,\n passing', 'not alternating,\n not passing'), rotation=45)

    plt.tight_layout()
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.savefig("plots/stepping_analysis.pdf", format="pdf")
    plt.close(fig)

def make_force_plot(args, force_data):
    fig = plt.figure(figsize=(4, 2.2), dpi=300)
    gennerich_forces = [4, 0, -1, -2, -3, -4, -5, -6, -7, -8, -10]
    gennerich_velocities = [50.76, 45.38, 43.07, 36.92, 22.3, 8.46, 4.6, 3.07, -.6, -5.38, -15.38]
    gennerich_errors = [2.3, 1.8, 3, 3, 1.6, 1.4, 1.5, 1.5, 1.4, 1.4, 2.3]
    plt.errorbar(gennerich_forces, gennerich_velocities, yerr=gennerich_errors, label="Gennerich 2007", fmt='o-', c='b', markersize=4, linestyle='', capsize=1, elinewidth=0.3, markeredgewidth=0.3)
    plt.scatter(force_data["forces"], force_data["velocities"], c='r', label="Model", zorder=2, s=4**2)
    plt.xlabel("$\hat{x}$ Force (pN)")
    plt.ylabel("Velocity (nm/s)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/force_vs_velocity.pdf", format="pdf")

def window_avg(times, array, width):
    if len(times) < width*3:
        return [], []
    indices = np.arange(0, len(times), width)
    w_times = np.zeros(len(indices)-2)
    w_array = np.zeros(len(indices)-2)
    for i, n in enumerate(indices[1:-1]):
        w_times[i] = times[int(n-width//2)]
        w_array[i] = np.mean(array[int(n-width):int(n)])
    return w_times, w_array

def make_stroke_plots(args, angles_data, longrun=False):
    plt.figure()

    if (angles_data == "NODATAFILE"):
        plt.figure(figsize=(4,4))
        plt.text(0.5, 0.5, "No data file \ndata/paper\_stroke\_angle\_data.txt.\n Run \nscripts/generate-angle-data.py \nto generate.",
                 horizontalalignment="center", verticalalignment="center", fontsize=17)

        plt.axis("off")
        plt.savefig("plots/onebound_stroke_taily_positions.pdf", format="pdf")
        plt.savefig("plots/onebound_stroke_tailx_positions.pdf", format="pdf")
        plt.savefig("plots/onebound_stroke_angles.pdf", format="pdf")
        plt.savefig("plots/bothbound_stroke_angles_bd.pdf", format="pdf")
        plt.savefig("plots/bothbound_stroke_taily_positions.pdf", format="pdf")
        plt.savefig("plots/bothbound_stroke_tailx_positions.pdf", format="pdf")
        plt.savefig("plots/bothbound_stroke_angles_md.pdf", format="pdf")
        return

    (onebound_angles_data, bothbound_angles_data) = angles_data

    if (len(onebound_angles_data["ma"]) <= 1 or len(bothbound_angles_data["ma"]) <= 1):
        print("paper_stroke_data* is too short, run for longer.")
        exit(1)

    plt.rcParams.update({'font.size': 21})

    dtailx_avg = np.nanmean(bothbound_angles_data["dtailx"], axis=0)
    taily_avg = np.nanmean(bothbound_angles_data["taily"], axis=0)
    ma_avg = np.nanmean(bothbound_angles_data["ma"], axis=0)
    ba_avg = np.nanmean(bothbound_angles_data["ba"], axis=0)

    if longrun:
        plt.figure()
        stds = np.zeros_like(bothbound_angles_data["longest_times"])
        for i, time in enumerate(bothbound_angles_data["longest_times"]):
            t_idx = np.where(bothbound_angles_data["longest_times"] == time)[0]
            stds[i] = np.nanstd(np.array([arr[int(t_idx)] for arr in bothbound_angles_data["ba"]]))
        plt.fill_between(bothbound_angles_data["longest_times"], (ba_avg+stds)*radians, (ba_avg-stds)*radians, zorder=2, yunits=radians, alpha=0.5)
        plt.plot(bothbound_angles_data["longest_times"], ba_avg*radians, linewidth=3, zorder=2, yunits=radians)
        plt.xlabel("time (ns)")
        plt.ylabel(r"$\theta_{ub}$ poststroke")
        plt.gca().set_xlim(0, bothbound_angles_data["longest_times"][-1])
        plt.gca().axhline(120 / 180 * np.pi * radians, color='blue', linestyle='dashed', linewidth=1)
        plt.tight_layout()
        plt.savefig("plots/bothbound_long_stroke_angles_bd.pdf", format="pdf")
        return

    # for num, ma in enumerate(bothbound_angles_data["ma"]):
    #     w_times, w_ma = window_avg(bothbound_angles_data["times"][num], ma, 1e2)
    #     plt.plot(w_times, w_ma*radians, zorder=1, yunits=radians)
    stds = np.zeros_like(bothbound_angles_data["longest_times"])
    for i, time in enumerate(bothbound_angles_data["longest_times"]):
        t_idx = np.where(bothbound_angles_data["longest_times"] == time)[0]
        stds[i] = np.nanstd(np.array([arr[int(t_idx)] for arr in bothbound_angles_data["ma"]]))# / (len(bothbound_angles_data["ma"]) - 1)
    plt.fill_between(bothbound_angles_data["longest_times"], (ma_avg+stds)*radians, (ma_avg-stds)*radians, zorder=2, yunits=radians, alpha=0.5)
    plt.plot(bothbound_angles_data["longest_times"], ma_avg*radians, linewidth=3, zorder=2, yunits=radians)
    plt.xlabel("time (ns)")
    plt.ylabel(r"$\theta_{um}$ poststroke")
    plt.gca().axhline(197 / 180 * np.pi * radians, color='red', linestyle='dashed', linewidth=1, yunits=radians)
    plt.gca().axhline(242 / 180 * np.pi * radians, color='blue', linestyle='dashed', linewidth=1, yunits=radians)
    plt.gca().set_xlim(0, bothbound_angles_data["longest_times"][-1])
    plt.tight_layout()
    plt.savefig("plots/bothbound_stroke_angles_md.pdf", format="pdf")

    plt.figure()
    # for num, dtailx in enumerate(bothbound_angles_data["dtailx"]):
    #     w_times, w_dtailx = window_avg(bothbound_angles_data["times"][num], dtailx, 1)
    #     # plt.scatter(0, dtailx[0], zorder=3, color="k")
    #     plt.plot(w_times, w_dtailx, zorder=1)
    stds = np.zeros_like(bothbound_angles_data["longest_times"])
    for i, time in enumerate(bothbound_angles_data["longest_times"]):
        t_idx = np.where(bothbound_angles_data["longest_times"] == time)[0]
        stds[i] = np.nanstd(np.array([arr[int(t_idx)] for arr in bothbound_angles_data["dtailx"]]))
    plt.fill_between(bothbound_angles_data["longest_times"], dtailx_avg+stds, dtailx_avg-stds, zorder=2, alpha=0.5)
    plt.plot(bothbound_angles_data["longest_times"], dtailx_avg, linewidth=3, zorder=2)
    plt.xlabel("time (ns)")
    plt.ylabel("tail x poststroke (nm)")
    plt.gca().set_xlim(0, bothbound_angles_data["longest_times"][-1])
    plt.tight_layout()
    plt.savefig("plots/bothbound_stroke_tailx_positions.pdf", format="pdf")

    plt.figure()
    # for num, taily in enumerate(bothbound_angles_data["taily"]):
    #     w_times, w_taily = window_avg(bothbound_angles_data["times"][num], taily, 1e2)
    #     # plt.scatter(0, taily[0], zorder=3, color="k")
    #     plt.plot(w_times, w_taily, zorder=1)
    stds = np.zeros_like(bothbound_angles_data["longest_times"])
    for i, time in enumerate(bothbound_angles_data["longest_times"]):
        t_idx = np.where(bothbound_angles_data["longest_times"] == time)[0]
        stds[i] = np.nanstd(np.array([arr[int(t_idx)] for arr in bothbound_angles_data["taily"]]))# / (len(bothbound_angles_data["ma"]) - 1)
    plt.fill_between(bothbound_angles_data["longest_times"], taily_avg+stds, taily_avg-stds, zorder=2, alpha=0.5)
    plt.plot(bothbound_angles_data["longest_times"], taily_avg, linewidth=3, zorder=2)
    plt.xlabel("time (ns)")
    plt.ylabel("tail y poststroke (nm)")
    plt.gca().set_xlim(0, bothbound_angles_data["longest_times"][-1])
    plt.tight_layout()
    plt.savefig("plots/bothbound_stroke_taily_positions.pdf", format="pdf")

    plt.figure()
    # for num, ba in enumerate(bothbound_angles_data["ba"]):
    #     w_times, w_ba = window_avg(bothbound_angles_data["times"][num], ba, 1e2)
    #     # plt.scatter(0, ba[0]*radians, zorder=3, color="k", yunits=radians)
    #     plt.plot(w_times, w_ba*radians, zorder=1, yunits=radians)
    stds = np.zeros_like(bothbound_angles_data["longest_times"])
    for i, time in enumerate(bothbound_angles_data["longest_times"]):
        t_idx = np.where(bothbound_angles_data["longest_times"] == time)[0]
        stds[i] = np.nanstd(np.array([arr[int(t_idx)] for arr in bothbound_angles_data["ba"]]))
    plt.fill_between(bothbound_angles_data["longest_times"], (ba_avg+stds)*radians, (ba_avg-stds)*radians, zorder=2, yunits=radians, alpha=0.5)
    plt.plot(bothbound_angles_data["longest_times"], ba_avg*radians, linewidth=3, zorder=2, yunits=radians)
    plt.xlabel("time (ns)")
    plt.ylabel(r"$\theta_{ub}$ poststroke")
    plt.gca().set_xlim(0, bothbound_angles_data["longest_times"][-1])
    plt.gca().axhline(120 / 180 * np.pi * radians, color='blue', linestyle='dashed', linewidth=1)
    plt.tight_layout()
    plt.savefig("plots/bothbound_stroke_angles_bd.pdf", format="pdf")

    dtailx_avg = np.nanmean(onebound_angles_data["dtailx"], axis=0)
    taily_avg = np.nanmean(onebound_angles_data["taily"], axis=0)
    ma_avg = np.nanmean(onebound_angles_data["ma"], axis=0)

    plt.figure()
    # for num, ma in enumerate(onebound_angles_data["ma"]):
    #     w_times, w_ma = window_avg(onebound_angles_data["times"][num], ma, 1e2)
    #     # plt.scatter(0, ma[0]*radians, zorder=3, color="k", yunits=radians)
    #     plt.plot(w_times, w_ma*radians, zorder=1, yunits=radians)
    stds = np.zeros_like(onebound_angles_data["longest_times"])
    for i, time in enumerate(onebound_angles_data["longest_times"]):
        t_idx = np.where(onebound_angles_data["longest_times"] == time)[0]
        stds[i] = np.nanstd(np.array([arr[int(t_idx)] for arr in onebound_angles_data["ma"]]))
    plt.fill_between(onebound_angles_data["longest_times"], (ma_avg+stds)*radians, (ma_avg-stds)*radians, zorder=2, yunits=radians, alpha=0.5)
    plt.plot(onebound_angles_data["longest_times"], ma_avg*radians, linewidth=3, zorder=2, yunits=radians)
    plt.xlabel("time (ns)")
    plt.ylabel(r"$\theta_{um}$ prestroke")
    plt.gca().axhline(197 / 180 * np.pi * radians, color='red', linestyle='dashed', linewidth=1)
    plt.gca().axhline(242 / 180 * np.pi * radians, color='blue', linestyle='dashed', linewidth=1)
    plt.gca().set_xlim(0, onebound_angles_data["longest_times"][-1])
    plt.tight_layout()
    plt.savefig("plots/onebound_stroke_angles.pdf", format="pdf")

    plt.figure()
    # for num, dtailx in enumerate(onebound_angles_data["dtailx"]):
    #     w_times, w_dtailx = window_avg(onebound_angles_data["times"][num], dtailx, 1e2)
    #     # plt.scatter(0, dtailx[0], zorder=3, color="k")
    #     plt.plot(w_times, w_dtailx, zorder=1)
    stds = np.zeros_like(onebound_angles_data["longest_times"])
    for i, time in enumerate(onebound_angles_data["longest_times"]):
        t_idx = np.where(onebound_angles_data["longest_times"] == time)[0]
        stds[i] = np.nanstd(np.array([arr[int(t_idx)] for arr in onebound_angles_data["dtailx"]])) # nanstd fixes this but why is it diff than others?
    plt.fill_between(onebound_angles_data["longest_times"], dtailx_avg+stds, dtailx_avg-stds, zorder=2, alpha=0.5)
    plt.plot(onebound_angles_data["longest_times"], dtailx_avg, linewidth=3, zorder=2)
    plt.xlabel("time (ns)")
    plt.ylabel("tail x prestroke (nm)")
    plt.gca().set_xlim(0, onebound_angles_data["longest_times"][-1])
    plt.tight_layout()
    plt.savefig("plots/onebound_stroke_tailx_positions.pdf", format="pdf")

    plt.figure()
    # for num, taily in enumerate(onebound_angles_data["taily"]):
    #     w_times, w_taily = window_avg(onebound_angles_data["times"][num], taily, 1e2)
    #     plt.plot(w_times, w_taily, zorder=1)
    stds = np.zeros_like(onebound_angles_data["longest_times"])
    for i, time in enumerate(onebound_angles_data["longest_times"]):
        t_idx = np.where(onebound_angles_data["longest_times"] == time)[0]
        stds[i] = np.nanstd(np.array([arr[int(t_idx)] for arr in onebound_angles_data["taily"]]))
    plt.fill_between(onebound_angles_data["longest_times"], taily_avg+stds, taily_avg-stds, zorder=2, alpha=0.5)
    plt.plot(onebound_angles_data["longest_times"], taily_avg, linewidth=3, zorder=2)
    plt.xlabel("time (ns)")
    plt.ylabel("tail y prestroke (nm)")
    plt.gca().set_xlim(0, onebound_angles_data["longest_times"][-1])
    plt.tight_layout()
    plt.savefig("plots/onebound_stroke_taily_positions.pdf", format="pdf")

# is it possible that our current results are true, and that dynein moves mostly on its prestroke, not poststroke?
# the best way to figure out this would be to get the average x-displacement during onebound and bothbound
# perhaps just measure the tx-displacement histogram for onebound and bothbound and see if there are interesting trends
# make ba plot
# augment the stepping data files to also show the ob vs bb tail displacement

# Maybe have make automatically make the data file since it's pretty small; have it in the repository
# Make an alternate trajectory movie which is smaller so you can see steps better

# at the current sim conditions (1e5 logging iterations, kb 1e8 kub 1e4), tailx displaces about 17.6nm in about .5us per step during the powerstroke, suggesting a max speed of 35 mm/s
# the binding domain during powerstroke reliably changes by about pi/6 or 30 degrees, suggesting a partial powerstroke mechanism (as opposed to winch)
# the tail actually increases its y coordinate by about 10nm, in opposition to the winch mechanism

def initialize():
    plt.rcParams.update({'font.size': 8})

def main():
    initialize()
    args = get_cli_arguments()
    stepping_data = get_stepping_data(args)
    unbinding_probability_data = get_unbinding_probability_data(args)
    force_data = get_force_data(args)
    angles_data = get_stroke_angles_data(args)
    long_angles_data = get_stroke_angles_data(args, longrun=True)
    make_behavior_plot(args, stepping_data, unbinding_probability_data)
    if args.quick:
        exit()
    make_analysis_plot(args, stepping_data)
    force_data = get_force_data(args)
    make_force_plot(args, force_data)
    make_stroke_plots(args, angles_data)
    make_stroke_plots(args, long_angles_data, longrun=True)

if __name__ == "__main__":
    main()
