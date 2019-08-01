import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
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
        print('running with inputs', bba, bma, uma, uba, state, k)
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


def collect_bothbound_data(self, P, nma, fma):
        """
        Collect bothbound statistics
        """
        P_arr.append(P)
        # Storing bb angles
        angles[0].append(nma)
        angles[1].append(fma)

        # Sum calculation for averages
        r_tx[0] += self.r_t[0]*P
        r_ty[0] += self.r_t[1]*P
        r_nmx[0] += self.r_nm[0]*P
        r_nmy[0] += self.r_nm[1]*P
        r_fmx[0] += self.r_fm[0]*P
        r_fmy[0] += self.r_fm[1]*P
        E[0] += self.E_total*P

        # Storing data for histograms
        r_t_arr[0].append(self.r_t[0])
        r_t_arr[1].append(self.r_t[1])
        r_nm_arr[0].append(self.r_nm[0])
        r_nm_arr[1].append(self.r_nm[1])
        r_fm_arr[0].append(self.r_fm[0])
        r_fm_arr[1].append(self.r_fm[1])
        E_arr.append(self.E_total)


def collect_onebound_data(k, state, bba, bma, uma, uba, L, step_data):
        """
        Call run_onebound function and collect onebound statistics
        """
        print('\n\nbothbound angles ',bba, bma, uma, uba, state)
        step = run_onebound(bba, bma, uma, uba, state, k[0])

        if state == 0:      # NEARBOUND State
            print('leading stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            step_data['L'].append(step['L'])        # Just leading step data
            step_data['step_length'].append(step['L']-L)        # Just leading step data
            step_length.append(step['L']-L)         # Contains both steps data
            # Storing final motor angles
            if step['L'] < 0:
                final_ang_arr[0].append(step['uma'])
                final_ang_arr[1].append(step['bma'])
            else:
                final_ang_arr[0].append(step['bma'])
                final_ang_arr[1].append(step['uma'])
        else:               # FARBOUND State
            print('trailing stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            step_data['L'].append(step['L'])        # Just trailing step data
            step_data['step_length'].append(step['L']+L)        # Just trailing step data
            step_length.append(step['L']+L)         # Contains both steps data
            # Storing final motor angles
            if step['L'] > 0:
                final_ang_arr[0].append(step['bma'])
                final_ang_arr[1].append(step['uma'])
            else:
                final_ang_arr[0].append(step['uma'])
                final_ang_arr[1].append(step['bma'])

        final_L_arr.append(step['L'])       # Final L array
        step_data['t'].append(step['t'])      # Specific type of step time array
        ob_t_arr.append(step['t'])          # Onebound time array
        k[0]+=1


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



params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=float, help="displacement in nm", default=8)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=100)
parser.add_argument("-k", "--kb", type=float, help="Manually set the binding rate", default=params.for_simulation['k_b'])
parser.add_argument("-t", "--dt", type=float, help="Manually set the dt", default=params.for_simulation['dt'])
args = parser.parse_args()

C =  params.for_simulation['exp-unbinding-constant']         # exponential binding constant from paper_params.py April 12

k_b = args.kb        # Binding Rate Constant
dt = args.dt         # Time Step

# Creating Data File for All L
data_file = open("../data/mc_data_{0}_{1}.txt".format(k_b, dt), "w")
data_file.write("#********mc_data: k_b-{0}, dt-{1}********\n\n\n".format(k_b, params.for_simulation['dt']))
data_file.write("#init L\t\t init nma\t init fma\t state\t\t final L\t final nma\t final fma\t step length\t t\n")


# for L in range(1, 33):
Z = 0                # Partition Function
L = args.L           # Initial Length
k = [0]              # Dynein Count & RNG Seed

# Sums for averages in Bothbound
r_tx = [0]                # Tail x position
r_ty = [0]                # Tail y position
r_nmx = [0]               # Near motor x position
r_nmy = [0]               # Near motor y position
r_fmx = [0]               # Far motor x position
r_fmy = [0]               # Far motor y position
E = [0]                   # Bothbound Energy
final_L = 0             # Final Displacements
step_L = 0              # Step Lengths
ob_t = 0                # One bound times

b = 11.82733524          # thermodynamic beta from default_parameters.h

rate_unbinding_leading = []                 # Leading (Far) Unbinding Rates
rate_unbinding_trailing = []                # Trailing (Near) Unbinding Rates

max_rate_trailing = 0
max_rate_leading = 0

P_arr = []
angles = [[] for i in range(2)]                # Pair of angles
trailing_data = { # Just trailing data
        'L': [],
        't': [],
        'step_length': [],
}
leading_data = { # Just leading data
        'L': [],
        't': [],
        'step_length': [],
}

r_t_arr = [[] for i in range(2)]                # Tail position array
r_nm_arr = [[] for i in range(2)]               # Near motor position array
r_fm_arr = [[] for i in range(2)]               # Far motor position array
E_arr = []                                      # Bothbound Energy array
final_L_arr = []                                # Final L array
final_ang_arr = [[] for i in range(2)]          # Final motor angles array
step_length = []                                # Step length array
ob_t_arr = []                                   # Onebound time array

seed = 0
np.random.seed(0)

eqb_angle = params.for_simulation['eqb']
if bb_energy_distribution.eq_in_degrees:
        eqb_angle = eqb_angle*np.pi/180

# Creating Data File for Specific L
# data_file = open("../data/mc_data_{0}_{1}_{2}.txt".format(int(L), k_b, dt), "w")
# data_file.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}********\n\n\n".format(L,
#                 k_b, params.for_simulation['dt']))
# data_file.write("#init L\t\t init nma\t init fma\t state\t\t final L\t final nma\t final fma\t step length\t t\n")

while Z < args.N:
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

                rate_trailing = np.exp(C*(dynein.nba - eqb_angle))
                rate_leading = np.exp(C*(dynein.fba - eqb_angle))
                rate_unbinding_trailing.append(rate_trailing)
                rate_unbinding_leading.append(rate_leading)
                max_rate_leading = max(rate_leading, max_rate_leading)
                max_rate_trailing = max(rate_trailing, max_rate_trailing)

                prob_trailing = P*rate_trailing
                prob_leading = P*rate_leading

                # print("prob_leading: ", prob_leading)
                # print("prob_trailing: ", prob_trailing)

                new_nma = nma-(np.pi-dynein.nba)
                new_fma = fma-(np.pi-dynein.fba)

                if np.random.random() < prob_trailing: # FIXME need to normalize this a tad so it is never > 1.
                        # FARBOUND State
                        state = 1

                        # while True:
                        collect_bothbound_data(dynein, P, nma, fma)


                        collect_onebound_data(k, state, dynein.fba, new_fma, new_nma, dynein.nba,
                                                L, trailing_data)
                            # if final_L_arr[k[0]-1] < -8 or -8 < final_L_arr[k[0]-1] < 8 or 8 < final_L_arr[k[0]-1]:
                            #     break
                        # plot_bb_before_step(dynein, 'red', 'blue')
                        # plt.savefig('../plots/mc_plots/trailing_{}a_before_step.png'.format(k), transparent=False)

                        # plot_bb_after_step(step['ubx'], step['uby'], step['umx'], step['umy'],
                        #                 step['tx'], step['ty'], step['bmx'], step['bmy'],
                        #                 step['bbx'], step['bby'], 'red', 'blue')
                        # plt.savefig('../plots/mc_plots/trailing_{}b_after_step.png'.format(k), transparent=False)
                        # plt.show()

                        data_file.write("{0:f}\t{1:f}\t{2:f}\t{3:s}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:f}\n".format(int(L),
                                        nma, fma, "FARBOUND", final_L_arr[k[0]-1], final_ang_arr[0][k[0]-1], final_ang_arr[1][k[0]-1],
                                        step_length[k[0]-1], ob_t_arr[k[0]-1]))

                if np.random.random() < prob_leading:
                        # NEARBOUND State
                        state = 0

                        # while True:
                        collect_bothbound_data(dynein, P, nma, fma)


                        collect_onebound_data(k, state, dynein.nba, new_nma, new_fma, dynein.fba,
                                                L, leading_data)
                            # if final_L_arr[k[0]-1] < -8 or -8 < final_L_arr[k[0]-1] < 8 or 8 < final_L_arr[k[0]-1]:
                            #     break
                        # plot_bb_before_step(dynein, 'red', 'blue')
                        # plt.savefig('../plots/mc_plots/leading_{}a_before_step.png'.format(k), transparent=False)

                        # plot_bb_after_step(step['bbx'], step['bby'], step['bmx'], step['bmy'],
                        #                 step['tx'], step['ty'], step['umx'], step['umy'],
                        #                 step['ubx'], step['uby'], 'red', 'blue')
                        # plt.savefig('../plots/mc_plots/leading_{}b_after_step.png'.format(k), transparent=False)
                        # plt.show()

                        data_file.write("{0:f}\t{1:f}\t{2:f}\t{3:s}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:f}\n".format(int(L),
                                        nma, fma, "NEARBOUND", final_L_arr[k[0]-1], final_ang_arr[0][k[0]-1], final_ang_arr[1][k[0]-1],
                                        step_length[k[0]-1], ob_t_arr[k[0]-1]))


# data_file.close()

print("FINAL DISPLACEMENTS: {0} \n".format(final_L_arr))
for i in range(len(final_L_arr)):
    final_L += final_L_arr[i]*P_arr[i]
    step_L += step_length[i]*P_arr[i]
    ob_t += ob_t_arr[i]*P_arr[i]

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

tx = r_tx[0]/Z          # Tail x
ty = r_ty[0]/Z          # Tail y
nmx = r_nmx[0]/Z        # Near motor x
nmy = r_nmy[0]/Z        # Near Motor y
fmx = r_fmx[0]/Z        # Far motor x
fmy = r_fmy[0]/Z        # Far Motor y
E_avg = E[0]/Z          # Average energy
mean_final_L = final_L/Z
mean_step_L = step_L/Z
mean_obt = ob_t/Z

print("BOTHBOUND AVERAGES")
print("Avg Tail x:", tx)
print("Avg Tail y:", ty)
print("Avg nmx:", nmx)
print("Avg nmy:", nmy)
print("Avg fmx:", fmx)
print("Avg fmy:", fmy)
print("Avg E:", E_avg)
print("Avg Final Displacement:", mean_final_L)
print("Avg Step Length:", mean_step_L)
print("Avg ob time:", mean_obt)

def make_hist(ax, stacked_hist, data, data0, bin, Label, Label0, tof, Color, Color0, Title, xlabel):
    ax.hist(data, bins=bin, alpha=0.5, label=Label, normed=tof, stacked=True, color=Color)
    if stacked_hist == True:
        ax.hist(data0, bins=bin, alpha=0.5, label=Label0, normed=tof, stacked=True, color=Color0)
    ax.legend(loc="upper right")
    ax.set_title(Title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")


fig0 = plt.figure(0, figsize=(12,8))
gs0 = gridspec.GridSpec(2, 21)
ax0 = fig0.add_subplot(gs0[0, 0:10])
ax1 = fig0.add_subplot(gs0[1, 0:10])
ax2 = fig0.add_subplot(gs0[0, 11:21])
ax3 = fig0.add_subplot(gs0[1, 11:21])

separate_step_hist = make_hist(ax0, True, trailing_data['L'], leading_data['L'], 30,
                    "Trailing Step", "Leading Step", True, "C0", "C1",
                    "Initial Displacement {0}nm".format(int(L)), "Final Displacement (nm)")
step_hist = make_hist(ax1, False, step_length, None, 50,
                    None, None, True, "C3", None,
                    "", "Step Length (nm)")
separate_time_hist = make_hist(ax2, True, trailing_data['t'], leading_data['t'], 30,
                    "Trailing time", "Leading time", False, "C0", "C1",
                    "k_b: {0:e}".format(k_b), "time (s)")
time_hist = make_hist(ax3, False, ob_t_arr, None, 50,
                    None, None, False, "C3", None,
                    "", "time (s)")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_onebound_length_time.pdf'.format(int(L), k_b), transparent=False)

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

initial_angle_hist = make_hist(ax4, True, angles[0], angles[1], 30,
                    "nma", "fma", True, "C0", "C1",
                    "Initial Both Bound Angles", "Initial Angles (rad)")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_bothbound_init_ang.pdf'.format(int(L), k_b), transparent=False)


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

tx_position_hist = make_hist(ax5, False, r_t_arr[0], None, 30,
                    "tx", None, True, "C0", None,
                    "Initial Both Bound Tail Position", "Tail x Positions")

ty_position_hist = make_hist(ax6, False, r_t_arr[1], None, 30,
                    "ty", None, True, "C1", None,
                    "", "Tail y Positions")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_bothbound_tail_position.pdf'.format(int(L), k_b), transparent=False)


fig3 = plt.figure(3, figsize=(6,8))
ax7 = fig3.add_subplot(gs2[0,:])
ax8 = fig3.add_subplot(gs2[1,:])

nmx_position_hist = make_hist(ax7, False, r_nm_arr[0], None, 30,
                    "nmx", None, True, "C0", None,
                    "Initial Both Bound Near Motor Position", "Near Motor x Positions")

nmy_position_hist = make_hist(ax8, False, r_nm_arr[1], None, 30,
                    "nmy", None, True, "C1", None,
                    "", "Near Motor y Positions")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_bothbound_nm_position.pdf'.format(int(L), k_b), transparent=False)


fig4 = plt.figure(4, figsize=(6,8))
ax9 = fig4.add_subplot(gs2[0,:])
ax10 = fig4.add_subplot(gs2[1,:])

fmx_position_hist = make_hist(ax9, False, r_fm_arr[0], None, 30,
                    "fmx", None, True, "C0", None,
                    "Initial Both Bound Far Motor Position", "Far Motor x Positions")

fmy_position_hist = make_hist(ax10, False, r_fm_arr[1], None, 30,
                    "fmy", None, True, "C1", None,
                    "", "Far Motor y Positions")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_bothbound_fm_position.pdf'.format(int(L), k_b), transparent=False)


fig5 = plt.figure(5)
ax11 = fig5.add_subplot(gs1[:,:])

Energy_hist = make_hist(ax11, False, E_arr, None, 30,
                    "Energies", None, True, "C0", None,
                    "Initial Both Bound Energy", "Energies")
plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_bothbound_energy.pdf'.format(int(L), k_b), transparent=False)

data_file.close()

plt.show()
