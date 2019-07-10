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

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", type=float, help="displacement in nm", required=True)
args = parser.parse_args()

C =  params.for_simulation['exp-unbinding-constant']         # exponential binding constant from paper_params.py April 12
# k_b = 3e7

Z = 0                # Partition Function
N = 100              # Count
L = args.L           # Length

# Sums for averages
r_tx = 0                # Tail x position
r_ty = 0                # Tail y position
r_nmx = 0               # Near motor x position
r_nmy = 0               # Near motor y position
r_fmx = 0               # Far motor x position
r_fmy = 0               # Far motor y position
E = 0                   # Energy Average
final_L = 0             # Final Displacements
step_L = 0              # Step Lengths
ob_t = 0                # One bound times

b = 11.82733524          # thermodynamic beta from default_parameters.h

rate_unbinding_leading = []                        # Leading (Far) Unbinding Rates
rate_unbinding_trailing = []                # Trailing (Near) Unbinding Rates

max_rate_trailing = 0
max_rate_leading = 0

P_arr = []
angles = [[] for i in range(2)]                # Pair of angles
trailing_data = [[] for i in range(2)]
leading_data = [[] for i in range(2)]

r_t_arr = [[] for i in range(2)]                # Tail position array
r_nm_arr = [[] for i in range(2)]               # Near motor position array
r_fm_arr = [[] for i in range(2)]               # Far motor position array
E_arr = []                                      # Energy Average array
final_L_arr = []
step_length = []
ob_t_arr = []


eqb_angle = params.for_simulation['eqb']
if bb_energy_distribution.eq_in_degrees:
        eqb_angle = eqb_angle*np.pi/180

seed = 0 # The initial seed for the C++ onebound code.
np.random.seed(0)

def run_onebound(bba, bma, uma, uba, state):
        global seed
        seed += 1 # use a different seed every time.  ugh, global variables!
        print('running with inputs', bba, bma, uma, uba, state)
        process = subprocess.Popen(['../onebound',
                                    str(params.for_simulation['k_b']),
                                    str(params.for_simulation['cb']),
                                    str(params.for_simulation['cm']),
                                    str(params.for_simulation['ct']),
                                    str(params.for_simulation['ls']),
                                    str(params.for_simulation['lt']),
                                    str(params.for_simulation['rt']),
                                    str(params.for_simulation['rm']),
                                    str(params.for_simulation['rb']),
                                    str(seed), # FIXME
                                    str(params.for_simulation['dt']),
                                    str(params.for_simulation['eqb']),
                                    str(params.for_simulation['eqmpre']),
                                    str(params.for_simulation['eqmpost']),
                                    str(params.for_simulation['eqt']),
                                    str(params.for_simulation['force']),
                                    str(params.for_simulation['exp-unbinding-constant']),
                                    str(bba), str(bma), str(uma), str(uba), str(state),
        ], stdout=subprocess.PIPE)
        (output, err) = process.communicate()
        exit_code = process.wait()
        output_data = eval(output.decode('utf-8'))
        # print('onebound for a step (which kind?) gives:\n%s\nEND OUTPUT' % output.decode('utf-8'))
        # print('data', output_data)
        assert(exit_code == 0);
        return output_data

# bool = True
# while bool == True:
#     Z = 0
#     final_L_arr = []
#     print("OLD K_B:", k_b)
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
                        P_arr.append(P)
                        angles[0].append(nma)
                        angles[1].append(fma)

                        # Calculation for averages
                        r_tx += dynein.r_t[0]*P
                        r_ty += dynein.r_t[1]*P
                        r_nmx += dynein.r_nm[0]*P
                        r_nmy += dynein.r_nm[1]*P
                        r_fmx += dynein.r_fm[0]*P
                        r_fmy += dynein.r_fm[1]*P
                        E += dynein.E_total*P

                        r_t_arr[0].append(dynein.r_t[0])
                        r_t_arr[1].append(dynein.r_t[1])
                        r_nm_arr[0].append(dynein.r_nm[0])
                        r_nm_arr[1].append(dynein.r_nm[1])
                        r_fm_arr[0].append(dynein.r_fm[0])
                        r_fm_arr[1].append(dynein.r_fm[1])
                        E_arr.append(dynein.E_total)

                        print('\n\ntrailing',dynein.fba,
                                            new_fma,
                                            new_nma,
                                            dynein.nba,
                                            state)
                        step = run_onebound(dynein.fba,
                                            new_fma,
                                            new_nma,
                                            dynein.nba,
                                            state)
                        print('trailing stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
                        trailing_data[0].append(step['L']+L)
                        trailing_data[1].append(step['t'])
                        final_L_arr.append(step['L'])
                        step_length.append(step['L']+L)
                        ob_t_arr.append(step['t'])

                if np.random.random() < prob_leading:
                        # NEARBOUND State
                        state = 0
                        P_arr.append(P)
                        angles[0].append(nma)
                        angles[1].append(fma)

                        # Calculation for averages
                        r_tx += dynein.r_t[0]*P
                        r_ty += dynein.r_t[1]*P
                        r_nmx += dynein.r_nm[0]*P
                        r_fmx += dynein.r_fm[0]*P
                        E += dynein.E_total*P

                        r_t_arr[0].append(dynein.r_t[0])
                        r_t_arr[1].append(dynein.r_t[1])
                        r_nm_arr[0].append(dynein.r_nm[0])
                        r_nm_arr[1].append(dynein.r_nm[1])
                        r_fm_arr[0].append(dynein.r_fm[0])
                        r_fm_arr[1].append(dynein.r_fm[1])
                        E_arr.append(dynein.E_total)

                        print('\n\nleading', dynein.nba,
                                            new_nma,
                                            new_fma,
                                            dynein.fba,
                                            state)
                        step = run_onebound(dynein.nba,
                                            new_nma,
                                            new_fma,
                                            dynein.fba,
                                            state)
                        print('leading stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
                        leading_data[0].append(step['L']-L)
                        leading_data[1].append(step['t'])
                        final_L_arr.append(step['L'])
                        step_length.append(step['L']-L)
                        ob_t_arr.append(step['t'])

                    #FIXME! Does not work
                    # print("Old K_B:", k_b)
                    # print("FINAL DISPLACEMENTS: {0} \n".format(final_L_arr))
                    # for i in range(0, len(final_L_arr)):
                    #     print("i:", i)
                    #     if i == len(final_L_arr) - 1:
                    #         bool = False
                    #         break
                    #     for j in range(i+1, len(final_L_arr)):
                    #         print("j:", j)
                    #         if final_L_arr[i] == final_L_arr[j]:
                    #             k_b = k_b * 1.3
                    #             Z = 30
                    #             break
                    #     break


# print("FINAL K_B:", k_b)
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

tx = r_tx/Z          # Tail x
ty = r_ty/Z          # Tail y
nmx = r_nmx/Z        # Near motor x
nmy = r_nmy/Z        # Near Motor y
fmx = r_fmx/Z        # Far motor x
fmy = r_fmy/Z        # Far Motor y
E_avg = E/Z          # Average energy
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
    plt.savefig('../plots/mc_plots/mc_{0}_{1}.pdf'.format(L, xlabel), transparent=False)


fig0 = plt.figure(0, figsize=(12,8))
gs0 = gridspec.GridSpec(2, 21)
ax0 = fig0.add_subplot(gs0[0, 0:10])
ax1 = fig0.add_subplot(gs0[1, 0:10])
ax2 = fig0.add_subplot(gs0[0, 11:21])
ax3 = fig0.add_subplot(gs0[1, 11:21])

separate_step_hist = make_hist(ax0, True, trailing_data[0], leading_data[0], 50,
                    "Trailing Step", "Leading Step", True, "C0", "C1",
                    "Initial Displacement {0}nm".format(L), "Step Length (nm)")
step_hist = make_hist(ax1, False, step_length, None, 50,
                    None, None, True, "C3", None,
                    "", "Step Length (nm)")
separate_time_hist = make_hist(ax2, True, trailing_data[1], leading_data[1], 50,
                    "Trailing time", "Leading time", False, "C0", "C1",
                    "Initial Displacement {0}nm".format(L), "time (s)")
time_hist = make_hist(ax3, False, ob_t_arr, None, 50,
                    None, None, False, "C3", None,
                    "", "time (s)")

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

initial_angle_hist = make_hist(ax4, True, angles[0], angles[1], 50,
                    "nma", "fma", True, "C0", "C1",
                    "Initial Both Bound Angles", "Initial Angles (rad)")

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

tx_position_hist = make_hist(ax5, False, r_t_arr[0], None, 50,
                    "tx", None, True, "C0", None,
                    "Initial Both Bound Tail Position", "Tail x Positions")

ty_position_hist = make_hist(ax6, False, r_t_arr[1], None, 50,
                    "ty", None, True, "C1", None,
                    "", "Tail y Positions")

fig3 = plt.figure(3, figsize=(6,8))
ax7 = fig3.add_subplot(gs2[0,:])
ax8 = fig3.add_subplot(gs2[1,:])

nmx_position_hist = make_hist(ax7, False, r_nm_arr[0], None, 50,
                    "nmx", None, True, "C0", None,
                    "Initial Both Bound Near Motor Position", "Near Motor x Positions")

nmy_position_hist = make_hist(ax8, False, r_nm_arr[1], None, 50,
                    "nmy", None, True, "C1", None,
                    "", "Near Motor y Positions")

fig4 = plt.figure(4, figsize=(6,8))
ax9 = fig4.add_subplot(gs2[0,:])
ax10 = fig4.add_subplot(gs2[1,:])

fmx_position_hist = make_hist(ax9, False, r_fm_arr[0], None, 50,
                    "fmx", None, True, "C0", None,
                    "Initial Both Bound Far Motor Position", "Far Motor x Positions")

fmy_position_hist = make_hist(ax10, False, r_fm_arr[1], None, 50,
                    "fmy", None, True, "C1", None,
                    "", "Far Motor y Positions")

fig5 = plt.figure(5)
ax11 = fig5.add_subplot(gs1[:,:])

Energy_hist = make_hist(ax11, False, E_arr, None, 50,
                    "Energies", None, True, "C0", None,
                    "Initial Both Bound Energy", "Energies")


plt.show()
