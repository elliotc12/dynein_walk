import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from statistics import mean
import scipy.constants
import sys
sys.path.append("../data")
import importlib
import argparse
import subprocess
import bb_energy_distribution

"""
Monte Carlo simulation for dynein taking a step
"""

def run_onebound(bba, bma, uma, uba, state, k):
        """
        Runs onebound.cpp with bb configuration and params.py
        """
        global seed
        seed += 1 # use a different seed every time.  ugh, global variables!
        # print('running with inputs', bba, bma, uma, uba, state, k)
        process = subprocess.Popen(['../onebound',
                                    str(k_b),
                                    str(params.for_simulation['cb']),
                                    str(params.for_simulation['cm']),
                                    str(params.for_simulation['ct']),
                                    str(params.for_simulation['ls']),
                                    str(params.for_simulation['lt']),
                                    str(params.for_simulation['rt']),
                                    str(params.for_simulation['rm']),
                                    str(params.for_simulation['rb']),
                                    str(seed),
                                    str(dt),
                                    str(params.for_simulation['eqb']),
                                    str(params.for_simulation['eqmpre']),
                                    str(params.for_simulation['eqmpost']),
                                    str(params.for_simulation['eqt']),
                                    str(params.for_simulation['force']),
                                    str(params.for_simulation['exp-unbinding-constant']),
                                    str(bba), str(bma), str(uma), str(uba), str(state), str(k),
        ], stdout=subprocess.PIPE)
        (output, err) = process.communicate()
        exit_code = process.wait()
        output_data = eval(output.decode('utf-8'))
        # print('onebound for a step (which kind?) gives:\n%s\nEND OUTPUT' % output.decode('utf-8'))
        # print('data', output_data)
        assert(exit_code == 0);
        return output_data


def collect_bothbound_data(k, self, P, state, nma, fma, prob, bb_data_file):
        """
        Collect bothbound statistics
        """
        P_arr.append(P)
        # Storing bb angles
        angles['nma'].append(nma)
        angles['fma'].append(fma)

        # Sum calculation for averages
        r_t['x_avg'] += self.r_t[0]*P
        r_t['y_avg'] += self.r_t[1]*P
        r_nm['x_avg'] += self.r_nm[0]*P
        r_nm['y_avg'] += self.r_nm[1]*P
        r_fm['x_avg'] += self.r_fm[0]*P
        r_fm['y_avg'] += self.r_fm[1]*P
        E['avg'] += self.E_total*P
        prob_unbinding['avg'] +=  prob*P

        # Storing data for histograms
        r_t['x'].append(self.r_t[0])
        r_t['y'].append(self.r_t[1])
        r_nm['x'].append(self.r_nm[0])
        r_nm['y'].append(self.r_nm[1])
        r_fm['x'].append(self.r_fm[0])
        r_fm['y'].append(self.r_fm[1])
        E['bb'].append(self.E_total)
        # FIXME These are wrong I'm pretty sure
        if state == 0:
            # NEARBOUND State - Leading step
            prob_unbinding['leading'].append(prob)
            if bb_data_file == True:
                data_file_bb_leading.write("{0:f}\t{1:f}\n".format(prob_unbinding['leading'][k[0]-len(prob_unbinding['trailing'])-1],
                prob_unbinding['unbinding'][k[0]]))
        else:
            # FARBOUND State - Trailing step
            prob_unbinding['trailing'].append(prob)
            if bb_data_file == True:
                data_file_bb_trailing.write("{0:f}\t{1:f}\n".format(prob_unbinding['trailing'][k[0]-len(prob_unbinding['leading'])-1],
                prob_unbinding['unbinding'][k[0]]))

def collect_onebound_data(k, state, bba, bma, uma, uba, L, step_data, ob_data_file):
        """
        Call run_onebound function and collect onebound statistics
        """
        # print('\n\nbothbound angles ',bba, bma, uma, uba, state)
        step = run_onebound(bba, bma, uma, uba, state, k[0])

        final_data['L'].append(step['L'])          # Final L array
        final_data['t'].append(step['t'])          # Onebound time array

        if state == 0:
            # NEARBOUND State - Leading step ddata
            print('leading stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            step_data['L'].append(step['L'])
            step_data['t'].append(step['t'])
            step_data['step_length'].append(step['L']-L)
            final_data['step_length'].append(step['L']-L)         # Contains both steps data

            # Storing final motor angles
            # if step['L'] < 0:   #FIXME uma != nma !!! need to calculate new bb nma
            #     final_data['nma'].append(step['uma'])
            #     final_data['fma'].append(step['bma'])
            # else:
            #     final_data['nma'].append(step['bma'])
            #     final_data['fma'].append(step['uma'])

            # if ob_data_file == True:
            #     data_file_ob_leading.write("{0:f}\t{1:e}\n".format(final_data['L'][k[0]],final_data['t'][k[0]]))
        else:
            # FARBOUND State - Trailing step data
            print('trailing stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            step_data['L'].append(step['L'])
            step_data['t'].append(step['t'])
            step_data['step_length'].append(step['L']+L)
            final_data['step_length'].append(step['L']+L)         # Contains both steps data

            # Storing final motor angles
            # if step['L'] > 0:
            #     final_data['nma'].append(step['bma'])
            #     final_data['fma'].append(step['uma'])
            # else:
            #     final_data['nma'].append(step['uma'])
            #     final_data['fma'].append(step['bma'])

            # if ob_data_file == True:
            #     data_file_ob_trailing.write("{0:f}\t{1:e}\n".format(final_data['L'][k[0]],final_data['t'][k[0]]))
        k[0]+=1
        sys.stdout.flush()


def plot_bb_before_step(self, dynein_color_nm, dynein_color_fm):
        """
        Plot just the figure of dynein for the both bound configuration given an
        array of motor angles and an initial displacement before the step.
        """
        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        x_coords_nm = [self.r_nb[0],
                        self.r_nm[0],
                        self.r_t[0]]

        y_coords_nm = [self.r_nb[1],
                        self.r_nm[1],
                        self.r_t[1]]

        x_coords_fm = [self.r_t[0],
                        self.r_fm[0],
                        self.r_fb[0]]

        y_coords_fm = [self.r_t[1],
                        self.r_fm[1],
                        self.r_fb[1]]


        ax.plot(x_coords_nm, y_coords_nm, color= dynein_color_nm, linewidth=3)
        ax.plot(x_coords_fm, y_coords_fm, color= dynein_color_fm, linewidth=3)
        ax.plot([-6, 30], [0, 0], color = 'black', linestyle='-', linewidth=3)
        ax.axis('off')
        ax.axis('equal')
        ax.legend()


def plot_bb_after_step(nbx, nby, nmx, nmy, tx, ty, fmx, fmy, fbx, fby, dynein_color_nm, dynein_color_fm):
        """
        Plot just the figure of dynein for the both bound configuration given an
        array of motor angles and an initial displacement after the step.
        """
        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        x_coords_nm = [nbx,
                        nmx,
                        tx]

        y_coords_nm = [nby,
                        nmy,
                        ty]

        x_coords_fm = [tx,
                        fmx,
                        fbx]

        y_coords_fm = [ty,
                        fmy,
                        fby]


        ax.plot(x_coords_nm, y_coords_nm, color= dynein_color_nm, linewidth=3)
        ax.plot(x_coords_fm, y_coords_fm, color= dynein_color_fm, linewidth=3)
        ax.plot([-6, 30], [0, 0], color = 'black', linestyle='-', linewidth=3)
        ax.axis('off')
        ax.axis('equal')
        ax.legend()


def best_fit(x,y):
    m = (((mean(x)*mean(y)) - mean(x*y))/
        ((mean(x)*mean(x)) - mean(x*x)))
    b = mean(y) - m*mean(x)
    return m, b


params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=float, help="displacement in nm", default=32)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=1e18)
parser.add_argument("-u", "--kub", type=float, help="Manually set the unbinding const", default=params.for_simulation['k_ub'])
parser.add_argument("-k", "--kb", type=float, help="Manually set the binding const", default=params.for_simulation['k_b'])
parser.add_argument("-cb", "--cb", type=float, help="Spring constant binding domain", default=params.for_simulation['cb'])
parser.add_argument("-cm", "--cm", type=float, help="Spring constant motor domain", default=params.for_simulation['cm'])
parser.add_argument("-ct", "--ct", type=float, help="Spring constant tail domain", default=params.for_simulation['ct'])
parser.add_argument("--eqb", type=float, help="Binding equilibrium angle", default=params.for_simulation['eqb'])
parser.add_argument("--eqmpre", type=float, help="Motor pre equilibrium angle", default=params.for_simulation['eqmpre'])
parser.add_argument("--eqmpost", type=float, help="Motor post equilibrium angle", default=params.for_simulation['eqmpost'])
parser.add_argument("-t", "--dt", type=float, help="Manually set the dt", default=params.for_simulation['dt'])
parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
parser.add_argument("-b", "--bb", type=bool, help="Collect Bothbound data", default=False)
parser.add_argument("-o", "--ob", type=bool, help="Colelct Onebound data", default=False)
args = parser.parse_args()

params.for_simulation['k_ub'] = args.kub
k_b = args.kb        # Binding Rate Constant
params.for_simulation['cb'] = args.cb
params.for_simulation['cm'] = args.cm
params.for_simulation['ct'] = args.ct
params.for_simulation['eqb'] = args.eqb
params.for_simulation['eqmpre'] = args.eqmpre
params.for_simulation['eqmpost'] = args.eqmpost
dt = args.dt          # Time Step
bb_data_file = args.bb
ob_data_file = args.ob

L = args.L           # Initial Length
N = args.N           # Count
Z = 0                # Partition Function
k = [0]              # Dynein Count
ts = 1e4
if dt < 1e-10:
    ts = 1e7

max_unbinding = 1
b = 11.82733524          # thermodynamic beta from default_parameters.h
eqb_angle = params.for_simulation['eqb']
if bb_energy_distribution.eq_in_degrees:
        eqb_angle = eqb_angle*np.pi/180

rate_unbinding_leading = []                 # Leading (Far) Unbinding Rates
rate_unbinding_trailing = []                # Trailing (Near) Unbinding Rates

max_rate_trailing = 0
max_rate_leading = 0

# Bothbound Data
P_arr = []
angles = {'nma': [],'fma': []}          # Pair of Angles
r_t = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}     # Tail data
r_nm = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}    # Near motor data
r_fm = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}    # Far motor data
E = {'bb': [], 'avg': 0}               # Energy data
rate_unbinding = {      # Unbinding rate
        'trailing': [],
        'leading': []
}
prob_unbinding = {      # Unbinding probability
        'trailing': [],
        'leading': [],
        'unbinding': [],
        'avg': 0
}

# Onebound Data
trailing_data = {   # Trailing step data
        'L': [],
        't': [],
        'step_length': [],
}
leading_data = {    # Leading step data
        'L': [],
        't': [],
        'step_length': [],
}
final_data = {
        'nma': [],
        'fma': [],
        'L': [],
        't': [],
        'step_length': [],
        'L_avg': 0,
        't_avg': 0,
        'step_length_avg': 0,
}

seed = 0
np.random.seed(0)

# Creating Data File for Specific L
if bb_data_file == True:
    data_file_bb_trailing = open("../data/mc_data_kb_{0:e}/mc_bb_trailing_data_{1}_{2}_{3}_{4}.txt".format(k_b, int(L), k_b, dt, N), "w")
    data_file_bb_trailing.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}, N-{3}, C-{4}********\n\n\n".format(L,
                    k_b, dt, N, C))
    data_file_bb_trailing.write("unbinding prob\t cumulative unbinding prob\n")
    data_file_bb_leading = open("../data/mc_data_kb_{0:e}/mc_bb_trailing_data_{1}_{2}_{3}_{4}.txt".format(k_b, int(L), k_b, dt, N), "w")
    data_file_bb_leading.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}, N-{3}, C-{4}********\n\n\n".format(L,
                    k_b, dt, N, C))
    data_file_bb_leading.write("unbinding prob\t cumulative unbinding prob\n")
# if ob_data_file == True:
#     data_file_ob_trailing = open("../data/mc_data_kb_{0:e}/mc_ob_trailing_data_{1}_{2}_{3}_{4}.txt".format(k_b, int(L), k_b, dt, N), "w")
#     data_file_ob_trailing.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}, N-{3}, C-{4}********\n\n\n".format(L,
#                     k_b, dt, N, C))
#     data_file_ob_trailing.write("final L\t t\n")
#     data_file_ob_leading = open("../data/mc_data_kb_{0:e}/mc_ob_leading_data_{1}_{2}_{3}_{4}.txt".format(k_b, int(L), k_b, dt, N), "w")
#     data_file_ob_leading.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}, N-{3}, C-{4}********\n\n\n".format(L,
#                     k_b, dt, N, C))
#     data_file_ob_leading.write("final L\t t\n")

while Z < N:
        # Making random motor angles
        nma = np.random.uniform(0, 2*np.pi)
        fma = np.random.uniform(0, 2*np.pi)

        dynein = bb_energy_distribution.DyneinBothBound(nma, fma, params, L)

        # Checking if energy is nan
        if np.isnan(dynein.E_total) == True:
                continue
        else:
                # Calculating partition function
                P = np.exp(-b*dynein.E_total)
                Z += P

                rate_trailing = np.exp(args.C*(dynein.nba - eqb_angle))
                rate_leading = np.exp(args.C*(dynein.fba - eqb_angle))
                rate_unbinding['trailing'].append(rate_trailing)
                rate_unbinding['leading'].append(rate_leading)
                max_rate_leading = max(rate_leading, max_rate_leading)
                max_rate_trailing = max(rate_trailing, max_rate_trailing)

                prob_trailing = P*rate_trailing     #   Unnormalized
                prob_leading = P*rate_leading       #   Unnormalized


                new_nma = nma-(np.pi-dynein.nba)
                new_fma = fma-(np.pi-dynein.fba)

                if np.random.random() < prob_trailing: # FIXME need to normalize this a tad so it is never > 1.
                        # FARBOUND State
                        state = 1
                        # collect_bothbound_data(k, dynein, P, state, nma, fma, prob_trailing, bb_data_file)


                        collect_onebound_data(k, state, dynein.fba, new_fma, new_nma, dynein.nba,
                                                L, trailing_data, ob_data_file)

                        # plot_bb_before_step(dynein, 'red', 'blue')
                        # plt.savefig('../plots/mc_plots/trailing_{}a_before_step.png'.format(k), transparent=False)

                        # plot_bb_after_step(step['ubx'], step['uby'], step['umx'], step['umy'],
                        #                 step['tx'], step['ty'], step['bmx'], step['bmy'],
                        #                 step['bbx'], step['bby'], 'red', 'blue')
                        # plt.savefig('../plots/mc_plots/trailing_{}b_after_step.png'.format(k), transparent=False)
                        # plt.show()
                        if k[0] % 10000 == 0:
                            np.savetxt('../data/mc_data/t_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.txt'.format(int(L),
                                        N, args.kub, k_b, dt, args.cb, args.cm, args.ct, args.C),
                                        (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')


                if np.random.random() < prob_leading:
                        # NEARBOUND State
                        state = 0
                        # collect_bothbound_data(k, dynein, P, state, nma, fma, prob_leading, bb_data_file)


                        collect_onebound_data(k, state, dynein.nba, new_nma, new_fma, dynein.fba,
                                                L, leading_data, ob_data_file)

                        # plot_bb_before_step(dynein, 'red', 'blue')
                        # plt.savefig('../plots/mc_plots/leading_{}a_before_step.png'.format(k), transparent=False)

                        # plot_bb_after_step(step['bbx'], step['bby'], step['bmx'], step['bmy'],
                        #                 step['tx'], step['ty'], step['umx'], step['umy'],
                        #                 step['ubx'], step['uby'], 'red', 'blue')
                        # plt.savefig('../plots/mc_plots/leading_{}b_after_step.png'.format(k), transparent=False)
                        # plt.show()
                        if k[0] % 10000 == 0:
                            np.savetxt('../data/mc_data/l_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.txt'.format(int(L),
                                        N, args.kub, k_b, dt, args.cb, args.cm, args.ct, args.C),
                                        (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')



print("FINAL DISPLACEMENTS: {0} \n".format(final_data['L']))
# for i in range(len(final_data['L'])):
#     final_data['L_avg'] += final_data['L'][i]*P_arr[i]
#     final_data['t_avg'] += final_data['t'][i]*P_arr[i]
#     final_data['step_length_avg'] += final_data['step_length'][i]*P_arr[i]

# What to collect and output or visualize?

### Bothbound data
# Mean angles while bothbound? (no stepping required)
# Mean motor/tail domain locations


### Stepping data (separately for leading/trailing)
# Final displacement (mean/histogram/list)
# Onebound time (mean/histogram/list)
# Rate of stepping

# print("rate_unbinding_leading: ", rate_unbinding_leading)
# print("rate_unbinding_trailing: ", rate_unbinding_trailing)
# print('max_rate_trailing', max_rate_trailing)
# print('max_rate_leading', max_rate_leading)

### What to export, and in what format?
# Histograms of final displacements?
#
# print('Averages')
# tx_avg = r_t['x_avg']/Z          # Tail x
# ty_avg = r_t['y_avg']/Z          # Tail y
# nmx_avg = r_nm['x_avg']/Z        # Near motor x
# nmy_avg = r_nm['y_avg']/Z        # Near Motor y
# fmx_avg = r_fm['x_avg']/Z        # Far motor x
# fmy_avg = r_fm['y_avg']/Z        # Far Motor y
# E_avg = E['avg']/Z          # Average energy
# prob_unbinding_avg = prob_unbinding['avg']/Z
# final_L_avg = final_data['L_avg']/Z
# step_length_avg = final_data['step_length_avg']/Z
# obt_avg = final_data['t_avg']/Z
#
# print("BOTHBOUND AVERAGES")
# print("Avg Tail x:", tx_avg)
# print("Avg Tail y:", ty_avg)
# print("Avg nmx:", nmx_avg)
# print("Avg nmy:", nmy_avg)
# print("Avg fmx:", fmx_avg)
# print("Avg fmy:", fmy_avg)
# print("Avg E:", E_avg)
# print("Avg prob_unbinding:", prob_unbinding_avg)
# print("Avg Final Displacement:", final_L_avg)
# print("Avg Step Length:", step_length_avg)
# print("Avg ob time:", obt_avg)


if bb_data_file == True:
    data_file_bb_trailing.close()
    data_file_bb_leading.close()
# if ob_data_file == True:
#     data_file_ob_trailing.close()
#     data_file_ob_leading.close()

np.savetxt('../data/mc_data/t_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.txt'.format(int(L),
            N, args.kub, k_b, dt, args.cb, args.cm, args.ct, args.C),
            (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
np.savetxt('../data/mc_data/l_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.txt'.format(int(L),
            N, args.kub, k_b, dt, args.cb, args.cm, args.ct, args.C),
            (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')

# END OF SIM


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
