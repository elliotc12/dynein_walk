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
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=1e20)
parser.add_argument("-b", "--kb", type=float, help="Binding const", default=params.for_simulation['k_b'])
parser.add_argument("-s", "--ks", type=float, help="Sticky const", default=params.for_simulation['k_stk'])
parser.add_argument("-cb", "--cb", type=float, help="Spring constant binding domain", default=params.for_simulation['cb'])
parser.add_argument("-cm", "--cm", type=float, help="Spring constant motor domain", default=params.for_simulation['cm'])
parser.add_argument("-ct", "--ct", type=float, help="Spring constant tail domain", default=params.for_simulation['ct'])
parser.add_argument("--eqb", type=float, help="Binding equilibrium angle", default=params.for_simulation['eqb'])
parser.add_argument("--eqmpre", type=float, help="Motor pre equilibrium angle", default=params.for_simulation['eqmpre'])
parser.add_argument("--eqmpost", type=float, help="Motor post equilibrium angle", default=params.for_simulation['eqmpost'])
parser.add_argument("-t", "--dt", type=float, help="Time step dt", default=params.for_simulation['dt'])
parser.add_argument("-C", "--C", type=float, help="exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
parser.add_argument("--underMT", action="store_false", help="Plot sims where binding domain can go under MT", default=True)
args = parser.parse_args()

k_b = float(args.kb)        # Binding Rate Constant
k_stk = float(args.ks)      # Sticky Rate Constant

u = ''
if args.underMT == False:
    u = 'u_'
params_string =  "{0:.2e}_{1:.2e}_{2}_{3}_{4}_{5}_{6}_{7}_{8}".format(k_b, k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C)


basepath = '../data/mc_data_' + params_string
datapath = '../data/compressed_mc_data/mc_data_file_' + u + params_string
plottingdatapath = '../data/mc_plotting_data/mc_plotting_data_' + u + params_string
leading_files = glob('{}/{}l_*.txt'.format(basepath, u))

initial_disp = []
final_disp_dict = {}
ob_time_dict = {}
affinity_time_dict = {}
max_final_disp = []

print('Making data for data file: \n ', basepath)

for leading in leading_files:
    leading = leading[len(basepath):]
    leading_data = {'L': [], 't': [], 't_affinity': []}
    trailing_data = {'L': [], 't': [], 't_affinity': []}
    trailing = leading.replace('l_', 't_')
    first_ = 2
    if args.underMT == False:
        first_ = 4
    second_ = leading.find('_', first_+1)
    third_ = leading.find('_', second_+1)
    fourth_ = leading.find('_', third_+1)
    fifth_ = leading.find('_', fourth_+1)
    sixth_ = leading.find('_', fifth_+1)
    seventh_ = leading.find('_', sixth_+1)
    eigth_ = leading.find('_', seventh_+1)
    iL = float(leading[first_+1:second_])
    if iL in initial_disp:
        print('woopsies, we have two files with the same L', iL, 'one of them is', leading)
        exit(1)
    # N = leading[second_+1:third_]
    # k_b = leading[third_+1:fourth_]
    # dt = leading[fourth_+1:fifth_]
    # cb = leading[fifth_+1:sixth_]
    # cm = leading[sixth_+1:seventh_]
    # ct = leading[seventh_+1:eigth_]
    # C = leading[eigth_+1:leading.rfind('.')]
    if len(np.loadtxt(basepath+leading)) > 0:
        leading_data['L'] = np.loadtxt(basepath+leading)[0]
        leading_data['t'] = np.loadtxt(basepath+leading)[1]
        # leading_data['t_affinity'] = np.loadtxt(basepath+leading)[2]
    try:
        if len(np.loadtxt(basepath+trailing)) > 0:
            trailing_data['L'] = np.loadtxt(basepath+trailing)[0]
            trailing_data['t'] = np.loadtxt(basepath+trailing)[1]
            # trailing_data['t_affinity'] = np.loadtxt(basepath+trailing)[2]
        initial_disp.append(iL)
        initial_disp.append(-iL)
        final_disp_dict[iL] = leading_data['L']
        final_disp_dict[-iL] = trailing_data['L']
        max_final_disp.append(np.max(leading_data['L']))
        max_final_disp.append(np.max(trailing_data['L']))

        if iL == 8.0 or iL == 16.0:
            ob_time_dict[iL] = leading_data['t']
            ob_time_dict[-iL] = trailing_data['t']
            leading_data['t_affinity'] = np.loadtxt(basepath+leading)[2]
            trailing_data['t_affinity'] = np.loadtxt(basepath+leading)[2]
            affinity_time_dict[iL] = leading_data['t_affinity']
            affinity_time_dict[-iL] = trailing_data['t_affinity']
    except:
        if not path.exists(leading.replace('l', 't', 1)):
            print('unable to load trailing data for ', leading)


if len(initial_disp) == 0:
    print("we have no data!!! :(")
    exit(1)

# make bin center the data point (for pcolor)
initial_disp = np.array(sorted(initial_disp)) # -50 to 50 array
final_bin_center = initial_disp*1.0 # set final_bin_center to initial_disp
final_bin_edges = np.zeros(len(final_bin_center)+1)
for i in range(1,len(final_bin_edges)-1):
    final_bin_edges[i] = (final_bin_center[i-1] + final_bin_center[i])*0.5
final_bin_edges[0] = 2*final_bin_center[0] - final_bin_edges[1]
final_bin_edges[-1] = 2*final_bin_center[-1] - final_bin_edges[-2]

# obtain meshgrid for pcolor
initial_disp_center, final_disp_center = np.meshgrid(initial_disp, final_bin_center)

final_disp_bin_width = final_bin_edges[1:] - final_bin_edges[:-1]    # a 1D array giving final displacement bin width

# Make 2d histogram plotting data for final displacements
hist = np.zeros_like(initial_disp_center)
normalized_hist = np.zeros_like(initial_disp_center)    # Dimensions: 1/distance

for i_disp in initial_disp:
    i_disp_index = 0
    for i in range(1,len(final_bin_edges)):
        if i_disp < final_bin_edges[i]:
            i_disp_index = i-1
            break
    total_counts = len(final_disp_dict[i_disp])
    for f_disp in final_disp_dict[i_disp]:
        f_disp_index = None
        for i in range(1,len(final_bin_edges)):
            if f_disp < final_bin_edges[i]:
                f_disp_index = i-1
                break
        if f_disp_index is None or f_disp < final_bin_edges[0] or f_disp > final_bin_edges[-1]:
            # Data that is outside of range goes into the bin edges.
            if f_disp < final_bin_edges[0]:
                hist[0, i_disp_index] += 1/total_counts
                normalized_hist[0, i_disp_index] += 1/total_counts/final_disp_bin_width[0]
            if f_disp > final_bin_edges[-1]:
                hist[-1, i_disp_index] += 1/total_counts
                normalized_hist[-1, i_disp_index] += 1/total_counts/final_disp_bin_width[-1]
            continue
            # print("crazasges", f_disp, 'vs', final_bin_edges[0], 'and', final_bin_edges[-1])
            # Possibly think about making a infinite bin for final_L that goes outside plot

        else:
            # Collect unnormalized counts in hist for Transition Matrix
            hist[f_disp_index, i_disp_index] += 1/total_counts  # Dimensionless

            # Normalized by bin width (length)
            normalized_hist[f_disp_index, i_disp_index] += 1/total_counts/final_disp_bin_width[f_disp_index]     # Dimensions: 1/distance

# Make 1d histogram plotting data for ob times
max_time = 1e-6 # 10 micro s
increment = 5e-9 #5 ns
time_bin_edges = np.arange(0.0, max_time+increment, increment, dtype=float)
time_bin_center = np.arange(increment/2, max_time, increment, dtype=float)
time_hists = {}
time_hists['max_time'] = max_time
time_hists['increment'] = increment
avg_affinity_time = {}

if len(ob_time_dict.keys()) == 0:
    print('Do not have data for initial L of 1nm or 8nm')
else:
    for i_disp in ob_time_dict.keys():
        avg_affinity_time[i_disp] = np.mean(affinity_time_dict[i_disp])
        time_hists[i_disp] = np.zeros_like(time_bin_center)
        total_time_counts = len(ob_time_dict[i_disp])

        for time in ob_time_dict[i_disp]:
            t_index = None
            for i in range(1, len(time_bin_edges)):
                if time < time_bin_edges[i]:
                    t_index = i-1
                    break
            if t_index is None:
                continue
            else:
                time_hists[i_disp][t_index] += 1/total_time_counts


# Saving 'plotting data' into ../data/mc_plotting_data
np.savez_compressed(plottingdatapath,
                    initial_disp=initial_disp,
                    hist=hist,
                    normalized_hist=normalized_hist,
                    time_hists=time_hists,
                    avg_affinity_time=avg_affinity_time)

# Saving 'Compressed data' into ../data/compressed_mc_data
np.savez_compressed(datapath,
                    final_disp_dict=final_disp_dict,
                    ob_time_dict=ob_time_dict)
