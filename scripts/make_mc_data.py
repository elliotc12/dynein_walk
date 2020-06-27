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


params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-k", "--kb", type=float, help="Binding const", default=params.for_simulation['k_b'])
parser.add_argument("-s", "--ks", type=float, help="Sticky const", default=params.for_simulation['k_stk'])
args = parser.parse_args()

k_b = float(args.kb)        # Binding Rate Constant
k_stk = float(args.ks)      # Sticky Rate Constant


basepath = '../data/mc_data_{0:.2e}_{1:.2e}/'.format(k_b, k_stk)
plotpath = '../plots/mc_plots/'
datapath = '../data/compressed_mc_data/mc_data_file_{0:.2e}_{1:.2e}'.format(k_b, k_stk)
leading_files = glob('{}/l_*.txt'.format(basepath))

initial_disp = []
final_disp_dict = {}
ob_time_dict = {} 

# probability of being a leading step
P_leading = []

for leading in leading_files:
    leading = leading[len(basepath):]
    leading_data = {'L': [], 't': []}
    trailing_data = {'L': [], 't': []}
    trailing = leading.replace('l_', 't_')
    first_ = leading.find('_',2)
    second_ = leading.find('_', first_+1)
    third_ = leading.find('_', second_+1)
    fourth_ = leading.find('_', third_+1)
    fifth_ = leading.find('_', fourth_+1)
    sixth_ = leading.find('_', fifth_+1)
    seventh_ = leading.find('_', sixth_+1)
    eigth_ = leading.find('_', seventh_+1)
    iL = float(leading[2:first_])

    if iL in initial_disp:
        print('woopsies, we have two files with the same L', iL, 'one of them is', leading)
        exit(1)
    N = leading[second_+1:third_]
    k_b = leading[third_+1:fourth_]
    dt = leading[fourth_+1:fifth_]
    cb = leading[fifth_+1:sixth_]
    cm = leading[sixth_+1:seventh_]
    ct = leading[seventh_+1:eigth_]
    C = leading[eigth_+1:leading.rfind('.')]
    leading_data['L'] = np.loadtxt(basepath+leading)[0]
    leading_data['t'] = np.loadtxt(basepath+leading)[1]
    try:
        trailing_data['L'] = np.loadtxt(basepath+trailing)[0]
        trailing_data['t'] = np.loadtxt(basepath+trailing)[1]
        initial_disp.append(iL)
        initial_disp.append(-iL)
        final_disp_dict[iL] = leading_data['L']
        final_disp_dict[-iL] = trailing_data['L']

        ob_time_dict[iL] = leading_data['t']
        ob_time_dict[-iL] = trailing_data['t']

        # calculate probability for step to be leading or trailing
        leading_data_length = len(leading_data['L'])
        trailing_data_length = len(trailing_data['L'])
        Prob_ld = leading_data_length / (leading_data_length + trailing_data_length)
        P_leading.append(Prob_ld)
        if path.exists(plotpath+'hist_final_L_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.pdf'.format(iL, N, k_b, dt, cb, cm, ct, C)):
            print('about to plot_hist', leading)
            plot_hist(iL, N, k_b, dt, cb, cm, ct, C)
    except:
        if path.exists(plotpath+'hist_final_L_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.pdf'.format(iL, N, k_b, dt, cb, cm, ct, C)):
            print('unable to load trailing data for', leading)

# Probability of Leading and Trailing Steps Based on Data
P_leading = np.array(P_leading)
P_trailing = 1-P_leading
P_unbinding = {'leading': P_leading, 'trailing': P_trailing}

if len(P_leading) == 0:
    print("we have no data!!! :(")
    exit(1)

print(final_disp_dict)
print(ob_time_dict)
print(P_unbinding)
np.savez_compressed(datapath, final_disp_dict=final_disp_dict, ob_time_dict=ob_time_dict, P_unbinding=P_unbinding)
