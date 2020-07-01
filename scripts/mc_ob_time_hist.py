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
datafiles = glob('{}/mc_data_file_*.npz'.format(datapath))

all_ob_time_dict = {}
kb_arr = []
kstk_arr = []
i = 0
for data in datafiles:
    mc_data = np.load(data, allow_pickle=True)
    data = data[len(datapath):]
    first_int = re.search(r"\d", data).start()
    first_ = data.find('_', first_int)
    k_b = float(data[first_int:first_])
    kb_arr.append(k_b)
    k_stk = float(data[first_+1:-4])
    kstk_arr.append(k_stk)
    ob_time_data = mc_data['ob_time_dict'].item()
    list_ob_time_data = list(ob_time_data.values())
    all_ob_time_dict[i] = np.concatenate(list_ob_time_data)
    i += 1

max_time = max(i for m in all_ob_time_dict.values() for i in m)
# 1D Time Histogram of all simulations to compare rate parameters
plt.figure('One Bound Times')
for j in range(len(kb_arr)):
    plt.hist(all_ob_time_dict[j], label='kb: {0:.2e}, kstk: {1:.2e}'.format(kb_arr[j], kstk_arr[j]), density=False,
            weights=np.ones(len(all_ob_time_dict[j]))/len(all_ob_time_dict[j]), stacked=True, histtype='step')
    plt.xticks(fontsize=10)
    plt.xscale('log')
plt.xlabel('time (s)')
plt.ylabel('Probability')
plt.legend()
plt.savefig(plotpath+'onebound_time_hist_all_rates.pdf'.format(float(k_b), float(k_stk)))
plt.show()
