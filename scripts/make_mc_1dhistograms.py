import os
from os import path
import numpy as np
from numpy.linalg import matrix_power
import matplotlib.pyplot as plt
from matplotlib import gridspec
from statistics import mean
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants
import sys
sys.path.append("../data")
import importlib
import argparse
import subprocess
import bb_energy_distribution
from glob import glob
import re

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
        # plt.legend()
        plt.savefig(plotpath+'onebound_time_1dhist_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))

def make_step_length_1dhist():
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
        plt.hist(step_length, label='Experiment', bins=np.linspace(min(step_length), max(step_length), num=50),
                 alpha=0.5, density=True, stacked=True, color="C1")
        plt.xticks(fontsize=8)
        plt.xlabel('Step length (nm)')
        plt.ylabel('Frequency')
        plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(k_b, k_stk))
        # plt.legend()
        plt.savefig(plotpath+'step_length_1hist_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))


# make_combined_ob_time_1dhist()
# make_individual_ob_time_1dhist()
make_step_length_1dhist()
plt.show()
