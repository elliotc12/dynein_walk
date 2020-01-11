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

def make_hist(ax, stacked_hist, data, data0, bin, Label, Label0, tof, Color, Color0, Title, xlabel):
    ax.hist(data, bins=bin, alpha=0.5, label=Label, normed=tof, stacked=True)
    if stacked_hist == True:
        ax.hist(data0, bins=bin, alpha=0.5, label=Label0, normed=tof, stacked=True)
    ax.legend(loc="upper right")
    ax.set_title(Title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")

def plot_hist(L, N, k_ub, k_b, dt, cb, cm, ct, C):
    fig = plt.figure()
    gs = gridspec.GridSpec(1,1)
    make_hist(fig.add_subplot(gs[:,:]), True, trailing_data['L'], leading_data['L'], 50,
                "Trailing Step", "Leading Step", False, "C0", "C1",
                "Initial Displacement: {}nm\nBinding Rate: {}{}\t dt: {}s".format(L,
                k_b, r'$s^{-1}$', dt), "Final Displacement (nm)")
    plt.savefig(plotpath+'hist_final_L_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.pdf'.format(L,
                N, k_b, dt, cb, cm, ct, C))

def best_fit(x,y):
    m = (((mean(x)*mean(y)) - mean(x*y))/
        ((mean(x)*mean(x)) - mean(x*x)))
    b = mean(y) - m*mean(x)
    return m, b

# DRAFT OF PROB FUNCTION, FIXME
#
# # Determining final location after a set number of steps
# def prob_steps(T, P, num_steps, P_lt):
#     P = np.matrix(np.zeros((len(T),1)))
#     # T_invs = np.asarray(list(reversed(T))
#     for i in range(len(P)): # for each initial displacement ...
#         P = np.matrix(np.zeros((len(T), 1)))
#         P[i] = 1
#         # set prob to P
#         prob = P
#         for j in range(num_steps): # step a number of times...
#             prob = T*P
#
#             if np.random.random_sample() > P_lt[i]:
#                 # then it is a trailing step.
#                 # if it's already a trailing step, leave same
#                 # otherwise, invert the position on probability matrix.
#
#
#                 # T *= T
#             # else:
#                 # T *= T_invs
#                 a = 1

def L_to_L(T, P_l, P_t):
    abs = np.zeros((50,100))
    prob_step = np.zeros((100,50))
    num_col = np.shape(abs)[1]
    num_rows = np.shape(abs)[0]
    for i in range(num_col):
        if i < num_col/2:
            abs[num_rows-1-i,i] = 1
            prob_step[num_rows-1-i,i] = P_t[num_rows-1-i]              # P_t array starts at 1 to 50
            prob_step[num_rows+i, i] = P_l[num_rows-1-i]        # P_l array starts at 1 to 50
        else:
            abs[i-num_rows,i] = 1
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

def rad_to_deg(angle):
    # array = angle*180/np.pi
    return angle

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
P_l = P_lt_1
P_t = int(1)-np.asarray(P_l)


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

hist = np.zeros_like(i_LLcenter)
hist2 = np.zeros_like(i_LLcenter)
normalized_hist = np.zeros_like(i_LLcenter)

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
        if fL_index is None or fL < final_L_edges[0]:
            continue
            # print("crazasges", fL, 'vs', final_L_edges[0], 'and', final_L_edges[-1])
            # Possibly think about making a infinite bin for final_L that goes outside plot
            # Will have some normalization issues
        else:
            hist[fL_index, iL_index] += 1/total_counts
            hist2[fL_index, iL_index] += 1

            # Normalized by area
            normalized_hist[fL_index, iL_index] += 1/total_counts/(final_L_edges[fL_index+1] - final_L_edges[fL_index])

i_LLedge, f_LLedge = np.meshgrid(final_L_edges, final_L_edges)

# slope, intercept = best_fit(i_LLedge, f_LLedge)
# rgrsn_line = [(slope*x)+intercept for x in np.asarray(i_LLedge)]
plt.close('all')
# plt.plot(i_LLedge, rgrsn_line, label='Model: y = ({:.3}) + ({:.3})x'.format(intercept,slope), linestyle=":")
plt.figure('From Data')
plt.pcolor(i_LLedge, f_LLedge, normalized_hist)
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()
plt.savefig(plotpath+'2dhist_initL_vs_finalL.pdf')

T = np.matrix(hist)
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

fig = plt.figure('prob plot')
prob_plot = fig.add_subplot(111)
new_hist = []
prob_den = []

num_steps = 8

# Obtain L to L probability density
for i in range(len(P)):
    P = np.matrix(np.zeros((int(len(T)/2), 1)))
    P[i] = 1
    # prob = (T**num_steps)*P
    prob = (L_to_L(T, P_l, P_t)**num_steps)*P
    # WORKING ON IT, FIXME:
    # prob = prob_steps(T, P, num_steps, P_lt)
    prob_flat = np.array(prob).flatten()
    prob_flat_norm = prob_flat.sum()*f_L_bin_width # sum of prob flat * bin width of both axis
    prob_den.append(np.array(prob_flat/prob_flat_norm))
    prob_plot.plot(f_L, prob_flat/prob_flat.sum(), label=f'i is {i}')

    plt.figure('prob density')
    plt.plot(f_L, prob_flat/prob_flat_norm, label=f'i is {i}')

# Plot L to L probability density
plt.figure('prob density')
plt.legend(loc='best')

prob_plot.legend()

prob_den_1 = list(reversed(prob_den))
prob_den_1.extend(prob_den)
prob_dx = []

for i in prob_den_1:
    reverse = list(reversed(i))
    reverse.extend(list(i))
    prob_dx.append(reverse)

bin_width = list(reversed(i_L_bin_width*1.0))
bin_width.extend(i_L_bin_width)


# plot the normalized histogram multiplied by the probability
final_normalized_hist = np.multiply(normalized_hist, prob_dx)
plt.figure('Match Yildiz divided by area of box')
plt.pcolor(i_LLedge, f_LLedge, final_normalized_hist/(np.power(bin_width,2)))
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()

plt.figure('Match Yildiz Not divided')
plt.pcolor(i_LLedge, f_LLedge, final_normalized_hist)
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()

plt.show()










# fig = plt.figure(figsize=(10,15))
#
# # make contourf graph
# ax1 = fig.add_subplot(111, projection='3d')
# ax1.scatter(init_ang[0], init_ang[1], final_L)
# # contour = ax1.contour(rad_to_deg(init_ang[0]), rad_to_deg(init_ang[1]), final_L, np.linspace(0, 1, 5), colors='w', linewidth=10)
# ax1.set_xlabel(r'$\theta_{nm}$', fontsize=60)
# ax1.set_ylabel(r'$\theta_{fm}$', fontsize=60)
# ax1.set_zlabel(r'Final L')
# # ax1.set_xticks(np.linspace(0, 360, 13))
# # ax1.set_yticks(np.linspace(0, 360, 13))
# # cb = plt.colorbar(FinalLPlot)
# # cb.set_label(r"Final L", fontsize=40)
# # cb.set_ticks(np.linspace(0, 1, 5))
# # cb.add_lines(contour)
#
# # FIXME PLEASE
#
#
# fig0 = plt.figure(0, figsize=(12,8))
# gs0 = gridspec.GridSpec(2, 21)
# ax0 = fig0.add_subplot(gs0[0, 0:10])
# ax1 = fig0.add_subplot(gs0[1, 0:10])
# ax2 = fig0.add_subplot(gs0[0, 11:21])
# ax3 = fig0.add_subplot(gs0[1, 11:21])
#
# separate_step_hist = make_hist(ax0, True, trailing_data['L'], leading_data['L'], 30,
#                     "Trailing Step", "Leading Step", True, "C0", "C1",
#                     "Initial Displacement {0}nm".format(int(L)), "Final Displacement (nm)")
# step_hist = make_hist(ax1, False, final_data['step_length'], None, 50,
#                     None, None, True, "C3", None,
#                     "", "Step Length (nm)")
# separate_time_hist = make_hist(ax2, True, trailing_data['t'], leading_data['t'], 30,
#                     "Trailing time", "Leading time", False, "C0", "C1",
#                     "k_b: {0:e}".format(k_b), "time (s)")
# time_hist = make_hist(ax3, False, final_data['t'], None, 50,
#                     None, None, False, "C3", None,
#                     "", "time (s)")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_onebound_length_time.pdf'.format(int(L), k_b, dt, N), transparent=False)
#
# # ax1.hist(final_L_arr, bins=50, alpha=0.5, normed=True, stacked=True, color="C2")
# # ax1.legend(loc="upper right")
# # ax1.set_xlabel("Final Displacement (nm)")
# # ax1.set_ylabel("Frequency")
#
# # ax2.hist(trailing_data[1], bins=50, alpha=0.5, label="Trailing time", normed=False, stacked=True, color="C0")
# # ax2.hist(leading_data[1], bins=50, alpha=0.5, label="Leading time", normed=False, stacked=True, color="C1")
# # ax2.legend(loc="upper right")
# # ax2.set_xlabel("time (s)")
# # ax2.set_ylabel("Frequency")
# #
# # ax3.hist(ob_t_arr, bins=50, alpha=0.5, normed=False, stacked=True, color="C2")
# # ax3.legend(loc="upper right")
# # ax3.set_xlabel("time (s)")
# # ax3.set_ylabel("Frequency")
#
# fig1 = plt.figure(1)
# gs1 = gridspec.GridSpec(1,1)
# ax4 = fig1.add_subplot(gs1[:,:])
#
# initial_angle_hist = make_hist(ax4, True, angles['nma'], angles['fma'], 30,
#                     "nma", "fma", True, "C0", "C1",
#                     "Initial Both Bound Angles", "Initial Angles (rad)")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_init_ang.pdf'.format(int(L), k_b, dt, N), transparent=False)
#
#
# # ax4.hist(angles[0], bins=50, alpha=0.5, label="nma", normed=True, stacked=True, color="C0")
# # ax4.hist(angles[1], bins=50, alpha=0.5, label="fma", normed=True, stacked=True, color="C1")
# # ax4.legend(loc="upper right")
# # ax4.set_title("Initial Displacement 8nm")
# # ax4.set_xlabel("Initial Angles")
# # ax4.set_ylabel("Frequency")
#
# fig2 = plt.figure(2, figsize=(6,8))
# gs2 = gridspec.GridSpec(2,1)
# ax5 = fig2.add_subplot(gs2[0,:])
# ax6 = fig2.add_subplot(gs2[1,:])
#
# tx_position_hist = make_hist(ax5, False, r_t['x'], None, 30,
#                     "tx", None, True, "C0", None,
#                     "Initial Both Bound Tail Position", "Tail x Positions")
#
# ty_position_hist = make_hist(ax6, False, r_t['y'], None, 30,
#                     "ty", None, True, "C1", None,
#                     "", "Tail y Positions")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_tail_position.pdf'.format(int(L), k_b, dt, N), transparent=False)
#
#
# fig3 = plt.figure(3, figsize=(6,8))
# ax7 = fig3.add_subplot(gs2[0,:])
# ax8 = fig3.add_subplot(gs2[1,:])
#
# nmx_position_hist = make_hist(ax7, False, r_nm['x'], None, 30,
#                     "nmx", None, True, "C0", None,
#                     "Initial Both Bound Near Motor Position", "Near Motor x Positions")
#
# nmy_position_hist = make_hist(ax8, False, r_nm['y'], None, 30,
#                     "nmy", None, True, "C1", None,
#                     "", "Near Motor y Positions")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_nm_position.pdf'.format(int(L), k_b, dt, N), transparent=False)
#
#
# fig4 = plt.figure(4, figsize=(6,8))
# ax9 = fig4.add_subplot(gs2[0,:])
# ax10 = fig4.add_subplot(gs2[1,:])
#
# fmx_position_hist = make_hist(ax9, False, r_fm['x'], None, 30,
#                     "fmx", None, True, "C0", None,
#                     "Initial Both Bound Far Motor Position", "Far Motor x Positions")
#
# fmy_position_hist = make_hist(ax10, False, r_fm['y'], None, 30,
#                     "fmy", None, True, "C1", None,
#                     "", "Far Motor y Positions")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_fm_position.pdf'.format(int(L), k_b, dt, N), transparent=False)
#
#
# fig5 = plt.figure(5)
# ax11 = fig5.add_subplot(gs1[:,:])
#
# Energy_hist = make_hist(ax11, False, E['bb'], None, 30,
#                     "Energies", None, True, "C0", None,
#                     "Initial Both Bound Energy", "Energies")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_energy.pdf'.format(int(L), k_b, dt, N), transparent=False)
#
# fig6 = plt.figure(6, figsize=(6,8))
# ax12 = fig6.add_subplot(gs2[0,:])
# ax13 = fig6.add_subplot (gs2[1,:])
# prob_separate_hist = make_hist(ax12, True, prob_unbinding['trailing'], prob_unbinding['leading'], 30,
#                     "Trailing", "Leading", True, "C0", "C1",
#                     "Unbinding Probabilities", "Probability")
# prob_hist = make_hist(ax13, False, prob_unbinding['unbinding'], None, 30,
#                     "Probabilities", None, True, "C0", None,
#                     "Cumulative Unbinding Probabilities", "Probability")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_unbinding_prob.pdf'.format(int(L), k_b, dt, N), transparent=False)
#
# plt.show()
