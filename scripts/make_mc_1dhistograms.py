import os
from os import path
import numpy as np
from numpy.linalg import matrix_power
import matplotlib.pyplot as plt
from matplotlib import gridspec
from statistics import mean
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants
import sys, re
sys.path.append("../data")
import importlib
import argparse
import subprocess
import bb_energy_distribution
from glob import glob

params = importlib.import_module("params")

plotpath = '../plots/mc_plots/'
datapath = '../data/compressed_mc_data/'
alldatafiles = glob('{}/mc_data_file_*.npz'.format(datapath))

kb_arr = []
kstk_arr = []
all_ob_time_dict = {}
i = 0
for data in alldatafiles:
    mc_data = np.load(data, allow_pickle=True)
    data = data[len(datapath):]
    first_int = re.search(r"\d", data).start()
    first_ = data.find('_', first_int)

    k_b = float(data[first_int:first_])
    kb_arr.append(k_b)
    k_stk = float(data[first_+1:-4])
    kstk_arr.append(k_stk)
    # Putting all ob times in 1 dict
    ob_time_data = mc_data['ob_time_dict'].item()
    list_ob_time_data = list(ob_time_data.values())
    all_ob_time_dict[i] = np.concatenate(list_ob_time_data)
    i += 1

def make_combined_ob_time_1dhist():
    max_time = max(i for m in all_ob_time_dict.values() for i in m)
    # 1D Time Histogram of all simulations to compare rate parameters
    plt.figure('One Bound Times')
    for j in range(len(kb_arr)):
        plt.hist(all_ob_time_dict[j], label='kb: {0:.2e}, kstk: {1:.2e}'.format(kb_arr[j], kstk_arr[j]), density=True,
                stacked=True, histtype='step')
        plt.xticks(fontsize=10)
        # plt.xscale('log')
    plt.xlabel('time (s)')
    plt.ylabel('Probability')
    plt.legend()
    plt.savefig(plotpath+'onebound_time_1dhist_all_rates.pdf')

def make_individual_ob_time_1dhist():
    for j in range(len(kb_arr)):
        k_b = kb_arr[j]
        k_stk = kstk_arr[j]
        datafile = '../data/compressed_mc_data/mc_data_file_{0:.2e}_{1:.2e}.npz'.format(k_b, k_stk)
        mc_data = np.load(datafile, allow_pickle=True)
        ob_time_dict = mc_data['ob_time_dict'].item()
        # print(ob_time_dict)
        initial_disp = ob_time_dict.keys()
        max_time = max(i for m in ob_time_dict.values() for i in m)

        # 1D Time Histogram
        plt.figure()
        for i_disp in initial_disp:
            plt.hist(ob_time_dict[i_disp], label='init disp: {}'.format(i_disp), bins=np.geomspace(1e-13, max_time, num=50), density=True,
                    stacked=True, histtype='step')
            plt.xticks(fontsize=8)
        plt.xlabel('time (s)')
        plt.xscale('log')
        plt.ylabel('Probability')
        plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(k_b, k_stk))
        plt.legend()
        plt.savefig(plotpath+'onebound_time_1dhist_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))

def make_step_length_1dhist():
    #step length histogram from Yildiz 2012 figure 1B top head-labelled panel
    yildiz_step_lengths = np.concatenate(([-37]*1, [-35]*1, [-34]*1, [-33]*2, [-31]*2, [-30]*3, [-29]*1, [-28]*1, [-27]*4, [-26]*4, [-25]*2, [-24]*3, [-23]*4, [-21]*4, [-20]*3,
                                          [-19]*3, [-18]*5, [-17]*3, [-16]*3, [-15]*7, [-14]*5, [-13]*7, [-12]*12, [-11]*16, [-10]*14, [-9]*20, [-8]*14, [-7]*10, [-6]*9, [-5]*11,
                                          [-4]*8, [-3]*2, [4]*6, [5]*7, [6]*12, [7]*20, [8]*19, [9]*22, [10]*30, [11]*34, [12]*26, [13]*21, [14]*23, [15]*22, [16]*30, [17]*29,
                                          [18]*23, [19]*22, [20]*26, [21]*12, [22]*21, [23]*16, [24]*7, [25]*8, [26]*7, [27]*8, [28]*5, [29]*9, [30]*7, [31]*8, [32]*6, [33]*2,
                                          [34]*2, [35]*9, [36]*4, [37]*9, [38]*5, [39]*1, [40]*2, [41]*1, [42]*3, [43]*2, [44]*4, [45]*4, [46]*1, [47]*1))

    for j in range(len(kb_arr)):
        k_b = kb_arr[j]
        k_stk = kstk_arr[j]
        datafile = '../data/compressed_mc_data/mc_data_file_{0:.2e}_{1:.2e}.npz'.format(k_b, k_stk)
        mc_data = np.load(datafile, allow_pickle=True)
        final_disp_dict = mc_data['final_disp_dict'].item()
        initial_disp = final_disp_dict.keys()

        step_length = []
        for i_disp in initial_disp:
            step_length.append(final_disp_dict[i_disp] - i_disp)
        step_length = np.concatenate(np.asarray(step_length), axis=None)

        plt.figure()
        plt.hist(yildiz_step_lengths, bins=np.linspace(min(step_length), max(step_length), num=50),
                alpha=0.5, label='Experiment', density=True, stacked=True, color="C0")
        plt.hist(step_length, label='Model', bins=np.linspace(min(step_length), max(step_length), num=50),
                 alpha=0.5, density=True, stacked=True, color="C1")
        plt.xticks(fontsize=8)
        plt.xlabel('Step length (nm)')
        plt.ylabel('Frequency')
        plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(k_b, k_stk))
        plt.legend()
        plt.savefig(plotpath+'step_length_1hist_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))


# make_combined_ob_time_1dhist()
make_individual_ob_time_1dhist()
# make_step_length_1dhist()
plt.show()
