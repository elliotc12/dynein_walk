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

#####################################################
####### WARNING!! THIS CODE DOES NOT WORK! ##########
#####################################################

# This script is a temporary place holder for old code


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


fig = plt.figure(figsize=(10,15))

# make contourf graph
ax1 = fig.add_subplot(111, projection='3d')
ax1.scatter(init_ang[0], init_ang[1], final_L)
# contour = ax1.contour(rad_to_deg(init_ang[0]), rad_to_deg(init_ang[1]), final_L, np.linspace(0, 1, 5), colors='w', linewidth=10)
ax1.set_xlabel(r'$\theta_{nm}$', fontsize=60)
ax1.set_ylabel(r'$\theta_{fm}$', fontsize=60)
ax1.set_zlabel(r'Final L')
# ax1.set_xticks(np.linspace(0, 360, 13))
# ax1.set_yticks(np.linspace(0, 360, 13))
# cb = plt.colorbar(FinalLPlot)
# cb.set_label(r"Final L", fontsize=40)
# cb.set_ticks(np.linspace(0, 1, 5))
# cb.add_lines(contour)

# FIXME PLEASE


fig0 = plt.figure(0, figsize=(12,8))
gs0 = gridspec.GridSpec(2, 21)
ax0 = fig0.add_subplot(gs0[0, 0:10])
ax1 = fig0.add_subplot(gs0[1, 0:10])
ax2 = fig0.add_subplot(gs0[0, 11:21])
ax3 = fig0.add_subplot(gs0[1, 11:21])

separate_step_hist = make_hist(ax0, True, trailing_data['L'], leading_data['L'], 30,
                    "Trailing Step", "Leading Step", True, "C0", "C1",
                    "Initial Displacement {0}nm".format(int(L)), "Final Displacement (nm)")
step_hist = make_hist(ax1, False, final_data['step_length'], None, 50,
                    None, None, True, "C3", None,
                    "", "Step Length (nm)")
separate_time_hist = make_hist(ax2, True, trailing_data['t'], leading_data['t'], 30,
                    "Trailing time", "Leading time", False, "C0", "C1",
                    "k_b: {0:e}".format(k_b), "time (s)")
time_hist = make_hist(ax3, False, final_data['t'], None, 50,
                    None, None, False, "C3", None,
                    "", "time (s)")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_onebound_length_time.pdf'.format(int(L), k_b, dt, N), transparent=False)

# ax1.hist(final_L_arr, bins=50, alpha=0.5, normed=True, stacked=True, color="C2")
# ax1.legend(loc="upper right")
# ax1.set_xlabel("Final Displacement (nm)")
# ax1.set_ylabel("Frequency")

# ax2.hist(trailing_data[1], bins=50, alpha=0.5, label="Trailing time", normed=False, stacked=True, color="C0")
# ax2.hist(leading_data[1], bins=50, alpha=0.5, label="Leading time", normed=False, stacked=True, color="C1")
# ax2.legend(loc="upper right")
# ax2.set_xlabel("time (s)")
# ax2.set_ylabel("Frequency")
#
# ax3.hist(ob_t_arr, bins=50, alpha=0.5, normed=False, stacked=True, color="C2")
# ax3.legend(loc="upper right")
# ax3.set_xlabel("time (s)")
# ax3.set_ylabel("Frequency")

fig1 = plt.figure(1)
gs1 = gridspec.GridSpec(1,1)
ax4 = fig1.add_subplot(gs1[:,:])

initial_angle_hist = make_hist(ax4, True, angles['nma'], angles['fma'], 30,
                    "nma", "fma", True, "C0", "C1",
                    "Initial Both Bound Angles", "Initial Angles (rad)")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_init_ang.pdf'.format(int(L), k_b, dt, N), transparent=False)


# ax4.hist(angles[0], bins=50, alpha=0.5, label="nma", normed=True, stacked=True, color="C0")
# ax4.hist(angles[1], bins=50, alpha=0.5, label="fma", normed=True, stacked=True, color="C1")
# ax4.legend(loc="upper right")
# ax4.set_title("Initial Displacement 8nm")
# ax4.set_xlabel("Initial Angles")
# ax4.set_ylabel("Frequency")

fig2 = plt.figure(2, figsize=(6,8))
gs2 = gridspec.GridSpec(2,1)
ax5 = fig2.add_subplot(gs2[0,:])
ax6 = fig2.add_subplot(gs2[1,:])

tx_position_hist = make_hist(ax5, False, r_t['x'], None, 30,
                    "tx", None, True, "C0", None,
                    "Initial Both Bound Tail Position", "Tail x Positions")

ty_position_hist = make_hist(ax6, False, r_t['y'], None, 30,
                    "ty", None, True, "C1", None,
                    "", "Tail y Positions")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_tail_position.pdf'.format(int(L), k_b, dt, N), transparent=False)


fig3 = plt.figure(3, figsize=(6,8))
ax7 = fig3.add_subplot(gs2[0,:])
ax8 = fig3.add_subplot(gs2[1,:])

nmx_position_hist = make_hist(ax7, False, r_nm['x'], None, 30,
                    "nmx", None, True, "C0", None,
                    "Initial Both Bound Near Motor Position", "Near Motor x Positions")

nmy_position_hist = make_hist(ax8, False, r_nm['y'], None, 30,
                    "nmy", None, True, "C1", None,
                    "", "Near Motor y Positions")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_nm_position.pdf'.format(int(L), k_b, dt, N), transparent=False)


fig4 = plt.figure(4, figsize=(6,8))
ax9 = fig4.add_subplot(gs2[0,:])
ax10 = fig4.add_subplot(gs2[1,:])

fmx_position_hist = make_hist(ax9, False, r_fm['x'], None, 30,
                    "fmx", None, True, "C0", None,
                    "Initial Both Bound Far Motor Position", "Far Motor x Positions")

fmy_position_hist = make_hist(ax10, False, r_fm['y'], None, 30,
                    "fmy", None, True, "C1", None,
                    "", "Far Motor y Positions")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_fm_position.pdf'.format(int(L), k_b, dt, N), transparent=False)


fig5 = plt.figure(5)
ax11 = fig5.add_subplot(gs1[:,:])

Energy_hist = make_hist(ax11, False, E['bb'], None, 30,
                    "Energies", None, True, "C0", None,
                    "Initial Both Bound Energy", "Energies")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_energy.pdf'.format(int(L), k_b, dt, N), transparent=False)

fig6 = plt.figure(6, figsize=(6,8))
ax12 = fig6.add_subplot(gs2[0,:])
ax13 = fig6.add_subplot (gs2[1,:])
prob_separate_hist = make_hist(ax12, True, prob_unbinding['trailing'], prob_unbinding['leading'], 30,
                    "Trailing", "Leading", True, "C0", "C1",
                    "Unbinding Probabilities", "Probability")
prob_hist = make_hist(ax13, False, prob_unbinding['unbinding'], None, 30,
                    "Probabilities", None, True, "C0", None,
                    "Cumulative Unbinding Probabilities", "Probability")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_unbinding_prob.pdf'.format(int(L), k_b, dt, N), transparent=False)

plt.show()

### FOUND AFTER MC_SIM.py
#
# def make_hist(ax, stacked_hist, data, data0, bin, Label, Label0, tof, Color, Color0, Title, xlabel):
#     ax.hist(data, bins=bin, alpha=0.5, label=Label, normed=tof, stacked=True, color=Color)
#     if stacked_hist == True:
#         ax.hist(data0, bins=bin, alpha=0.5, label=Label0, normed=tof, stacked=True, color=Color0)
#     ax.legend(loc="upper right")
#     ax.set_title(Title)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel("Frequency")

# def plot_hist(L, k_b, dt, N):
# fig0 = plt.figure(0, figsize=(12,8))
# gs0 = gridspec.GridSpec(2, 22)
# gs0 = gridspec.GridSpec(1,1)
# ax0 = fig0.add_subplot(gs0[0, 0:10])
# ax1 = fig0.add_subplot(gs0[1, 0:10])
# ax2 = fig0.add_subplot(gs0[0, 12:22])
# ax3 = fig0.add_subplot(gs0[1, 12:22])
#
# fig0 = plt.figure(0)
# ax0 = fig0.add_subplot(gs0[:,:])
# separate_step_hist = make_hist(ax0, True, trailing_data['L'], leading_data['L'], 50,
#                     "Trailing Step", "Leading Step", False, "C0", "C1",
#                     "Initial Displacement: {}nm\nBinding Rate: {:.1e}{}\t dt: {}s".format(int(L), k_b, r'$s^{-1}$', dt), "Final Displacement (nm)")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_hist_final_L.png'.format(int(L), k_b, dt, N), transparent=True)
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_hist_final_L.svg'.format(int(L), k_b, dt, N), transparent=True)
#
# fig1 = plt.figure(1)
# ax1 = fig1.add_subplot(gs0[:,:])
# step_hist = make_hist(ax1, True, trailing_data['step_length'], leading_data['step_length'], 50,
#                     "Trailing Step", "Leading Step", False, "C0", "C1",
#                     "Initial Displacement: {}nm\nBinding Rate: {:.1e}{}\t dt: {}s".format(int(L), k_b, r'$s^{-1}$', dt), "Step Length (nm)")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_hist_step_length.png'.format(int(L), k_b, dt, N), transparent=True)
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_hist_step_length.svg'.format(int(L), k_b, dt, N), transparent=True)
#
# fig2 = plt.figure(2)
# ax2 = fig2.add_subplot(gs0[:,:])
# separate_time_hist = make_hist(ax2, True, np.array(trailing_data['t'])*ts, np.array(leading_data['t'])*ts, 50,
#                     "Trailing time", "Leading time", False, "C0", "C1",
#                     "Initial Displacement: {}nm\nBinding Rate: {:.1e}{}\t dt: {}s".format(int(L), k_b, r'$s^{-1}$', dt), r'time ($\mu$s)')
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_hist_time.png'.format(int(L), k_b, dt, N), transparent=True)
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_hist_time.svg'.format(int(L), k_b, dt, N), transparent=True)
#
# fig3 = plt.figure(3)
# ax3 = fig3.add_subplot(gs0[:,:])
# time_hist = make_hist(ax3, False, np.array(final_data['t'])*ts, None, 50,
#                     None, None, False, "C3", None,
#                     "Initial Displacement: {}nm\nBinding Rate: {:.1e}{}\t dt: {}s".format(int(L), k_b, r'$s^{-1}$', dt), r'time ($\mu$s)')
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_hist_time_all.png'.format(int(L), k_b, dt, N), transparent=True)
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_hist_time_all.svg'.format(int(L), k_b, dt, N), transparent=True)
#
# plt.show()

# plot_hist(L, k_b, dt, N)
