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

def best_fit(x,y):
    m = (((mean(x)*mean(y)) - mean(x*y))/
        ((mean(x)*mean(x)) - mean(x*x)))
    b = mean(y) - m*mean(x)
    return m, b

def rad_to_deg(angle):
    # array = angle*180/np.pi
    return angle

def L_to_initial_displacement(P_l, P_t):
    num_col = 2*len(P_l)
    num_rows = len(P_l)
    prob_step = np.zeros((num_col,num_rows))
    for i in range(num_col):
        if i < num_col/2:
            prob_step[num_rows-1-i,i] = P_t[num_rows-1-i]              # P_t array starts at 1 to 50
            prob_step[num_rows+i, i] = P_l[num_rows-1-i]        # P_l array starts at 1 to 50
    return prob_step

def L_to_L(T, P_l, P_t):
    num_col = 2*len(P_l)
    num_rows = len(P_l)
    abs = np.zeros((num_rows,num_col))
    for i in range(num_col):
        if i < num_col/2:
            abs[num_rows-1-i,i] = 1
        else:
            abs[i-num_rows,i] = 1
    prob_step = L_to_initial_displacement(P_l, P_t)
    T_L = abs*T*prob_step
    # plt.figure('abs')
    # plt.pcolor(abs);
    # plt.colorbar();
    # plt.figure('prob_step')
    # plt.pcolor(prob_step);
    # plt.colorbar();
    # plt.show()
    return T_L


params = importlib.import_module("params")

parser = argparse.ArgumentParser()
# parser.add_argument("-L", "--L", type=float, help="displacement in nm", required=True)
# parser.add_argument("-k", "--kb", type=float, help="binding rate", required=True)
# parser.add_argument("-t", "--dt", type=float, help="dt", required=True)
args = parser.parse_args()


basepath = '../data/mc_data/'
plotpath = '../plots/mc_plots/'
leading_files = glob('{}/l_*.txt'.format(basepath))

initial_L = []
final_L_lists = {}

# probability of being a leading step
P_lt = []
P_lt_1 = []

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

    if iL in initial_L:
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
        initial_L.append(iL)
        initial_L.append(-iL)
        final_L_lists[iL] = leading_data['L']
        final_L_lists[-iL] = trailing_data['L']

        # calculate probability for step to be leading or trailing
        leading_data_length = len(leading_data['L'])
        trailing_data_length = len(trailing_data['L'])
        Prob_lt = leading_data_length / (leading_data_length + trailing_data_length)
        P_lt.append(Prob_lt)
        P_lt_1.append(Prob_lt)
        if path.exists(plotpath+'hist_final_L_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.pdf'.format(iL, N, k_b, dt, cb, cm, ct, C)):
            print('about to plot_hist', leading)
            plot_hist(iL, N, k_b, dt, cb, cm, ct, C)
    except:
        if path.exists(plotpath+'hist_final_L_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.pdf'.format(iL, N, k_b, dt, cb, cm, ct, C)):
            print('unable to load trailing data for', leading)

# combine to make leading/trailing probability
P_lt = list(reversed(P_lt))
P_lt.extend(P_lt_1)
P_l = np.array(P_lt_1)
P_t = int(1)-P_l


# make bin center the data point (for pcolor)
initial_L = np.array(sorted(initial_L)) # -50 to 50 array
final_L_center = initial_L*1.0 # set final_L_center to initial_L
final_L_edges = np.zeros(len(final_L_center)+1)
for i in range(1,len(final_L_edges)-1):
    final_L_edges[i] = (final_L_center[i-1] + final_L_center[i])*0.5
final_L_edges[0] = 2*final_L_center[0] - final_L_edges[1]
final_L_edges[-1] = 2*final_L_center[-1] - final_L_edges[-2]
# obtain meshgrid for pcolor
i_LLcenter, f_LLcenter = np.meshgrid(initial_L, final_L_center)

bin_width = final_L_edges[1:] - final_L_edges[:-1]

hist = np.zeros_like(i_LLcenter)
filtered_hist = np.zeros_like(i_LLcenter)
hist2 = np.zeros_like(i_LLcenter)
normalized_hist = np.zeros_like(i_LLcenter)
filtered_norm_hist = np.zeros_like(i_LLcenter)
epsilon = 4

for iL in initial_L:
    iL_index = 0
    for i in range(1,len(final_L_edges)):
        if iL < final_L_edges[i]:
            iL_index = i-1
            break
    total_counts = len(final_L_lists[iL])
    for fL in final_L_lists[iL]:
        fL_index = None
        for i in range(1,len(final_L_edges)):
            if fL < final_L_edges[i]:
                fL_index = i-1
                break
        if fL_index is None or fL < final_L_edges[0] or fL > final_L_edges[-1]:
            # This seems fishy...
            if fL < final_L_edges[0]:
                normalized_hist[0, iL_index] += 1/total_counts/bin_width[0]
                hist[0, iL_index] += 1/total_counts
                if abs(fL-iL) > epsilon:
                    filtered_norm_hist[0, iL_index] += 1/total_counts/bin_width[0]
                    filtered_hist[0, iL_index] += 1/total_counts
            if fL > final_L_edges[-1]:
                normalized_hist[-1, iL_index] += 1/total_counts/bin_width[-1]
                hist[-1, iL_index] += 1/total_counts
                if abs(fL-iL) > epsilon:
                    filtered_norm_hist[-1, iL_index] += 1/total_counts/bin_width[-1]
                    filtered_hist[-1, iL_index] += 1/total_counts
            continue
            # print("crazasges", fL, 'vs', final_L_edges[0], 'and', final_L_edges[-1])
            # Possibly think about making a infinite bin for final_L that goes outside plot
            # Will have some normalization issues

        else:
            hist[fL_index, iL_index] += 1/total_counts
            hist2[fL_index, iL_index] += 1

            # Normalized by area
            normalized_hist[fL_index, iL_index] += 1/total_counts/bin_width[fL_index]
            if abs(fL-iL) > epsilon:
                filtered_norm_hist[fL_index, iL_index] += 1/total_counts/bin_width[fL_index]
                filtered_hist[fL_index, iL_index] += 1/total_counts

# for i in range(len(normalized_hist)):
    # print('norm', i, (normalized_hist[:,i]*bin_width).sum())
    # print('hist', i, hist[:,i].sum())

i_LLedge, f_LLedge = np.meshgrid(final_L_edges, final_L_edges)


plt.close('all')
plt.figure('From Data')
plt.pcolor(i_LLedge, f_LLedge, normalized_hist)
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()
plt.savefig(plotpath+'2dhist_initL_vs_finalL.pdf')

plt.figure('Filtered Instant Rebinds Data')
plt.pcolor(i_LLedge, f_LLedge, filtered_norm_hist)
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()

T = np.matrix(hist)
filtered_T = np.matrix(filtered_hist)
P = np.matrix(np.zeros((int(len(T)/2),1)))
f_L = np.array(final_L_center[len(final_L_center)//2:])
f_L_bin_width = np.zeros_like(f_L)
f_L_bin_width[1:-1] = (f_L[2:] - f_L[:-2])/2
f_L_bin_width[0] = f_L[0] + (f_L[1] - f_L[0])/2
f_L_bin_width[-1] = f_L[-1] - f_L[-2]

# Get bin widths for initial displacement histogram
i_L = np.array(initial_L[len(initial_L)//2:])
i_L_bin_width = np.zeros_like(i_L)
i_L_bin_width[1:-1] = (i_L[2:] - i_L[:-2])/2
i_L_bin_width[0] = i_L[0] + (i_L[1] - i_L[0])/2
i_L_bin_width[-1] = i_L[-1] - i_L[-2]

new_hist = []

num_steps = 16

plt.figure('prob density')
# Plot L to L probability density
plt.legend(loc='best')
plt.xlabel('L')
plt.ylabel('Prob per L')

# Obtain L to L probability density
for i in range(len(P)):
    P[:,:]= 0
    P[i] = 1
    # prob = (T**num_steps)*P
    prob = (L_to_L(T, P_l, P_t)**num_steps)*P
    # WORKING ON IT, FIXME:
    # prob = prob_steps(T, P, num_steps, P_lt)
    prob_flat = np.array(prob)[:,0] # convert to a 1D array from a column vector
    print('prob_flat.sum()', prob_flat.sum())
    prob_flat_norm = prob_flat.sum()*f_L_bin_width # sum of prob flat * bin width of both axis
    # print('prob_flat_norm', prob_flat_norm)
    prob_den = prob_flat/prob_flat_norm
    print('prob_den SUM:', prob_den.sum())
    plt.plot(f_L, prob_den, label=f'i is {i}')


# Filtered Prob Density Function
plt.figure('filtered prob density')
plt.legend(loc='best')
plt.xlabel('L')
plt.ylabel('Prob per L')

# Obtain L to L probability density
for i in range(len(P)):
    P[:,:]= 0
    P[i] = 1
    # prob = (T**num_steps)*P
    filtered_prob = (L_to_L(filtered_T, P_l, P_t)**num_steps)*P
    filtered_prob_flat = np.array(filtered_prob)[:,0] # convert to a 1D array from a column vector
    # print('filtered_prob_flat.sum()', prob_flat.sum())
    filtered_prob_flat_norm = filtered_prob_flat.sum()*f_L_bin_width # sum of prob flat * bin width of both axis
    filtered_prob_den = filtered_prob_flat/filtered_prob_flat_norm

    plt.plot(f_L, filtered_prob_den, label=f'i is {i}')

# print(prob_den)
prob_dx = L_to_initial_displacement(P_l, P_t).dot(prob_den)
print(prob_dx.shape)
plt.figure('prob_dx')
plt.plot(initial_L, prob_dx)
plt.xlabel('displacement')
plt.ylabel('probability density')

filtered_prob_dx = list(reversed(filtered_prob_den))
filtered_prob_dx.extend(filtered_prob_den)
filtered_prob_dx = np.array(filtered_prob_dx)

bin_width = list(reversed(i_L_bin_width*1.0))
bin_width.extend(i_L_bin_width)


# plot the normalized histogram multiplied by the probability
final_normalized_hist = np.zeros_like(normalized_hist)
for i in range(final_normalized_hist.shape[0]):
    final_normalized_hist[i,:] = normalized_hist[i,:]*prob_dx

filtered_final_norm_hist = np.zeros_like(filtered_norm_hist)
for i in range(filtered_final_norm_hist.shape[0]):
    filtered_final_norm_hist[i,:] = filtered_norm_hist[i,:]*filtered_prob_dx

# print(np.sum(final_normalized_hist))

plt.figure('Probability Distribution to Match Yildiz')
plt.pcolor(i_LLedge, f_LLedge, final_normalized_hist)
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()

plt.figure('Filtered Probability Distribution to Match Yildiz')
plt.pcolor(i_LLedge, f_LLedge, filtered_final_norm_hist)
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()

print('SUM:', np.sum(normalized_hist))
print('FINAL SUM:', np.sum(final_normalized_hist))

dx = 1
v = 0
for i in range(len(final_normalized_hist)):
    for j in range(len(final_normalized_hist[i])):
        v += final_normalized_hist[i][j]*bin_width[i]*bin_width[j]

print('Integral of P(x,y)', v)  # FIXME!

# print(final_normalized_hist)
# print(np.shape(final_normalized_hist))

# print(i_LLcenter)
# print(np.shape(i_LLcenter))

# line_best = np.polyfit(i_LLcenter[0], final_normalized_hist, 1)
# print(line_best)


plt.show()
