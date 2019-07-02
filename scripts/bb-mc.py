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

Z = 0                # Partition Function
N = 500              # Count
L = args.L           # Length

# Sums for averages
r_tx = 0                # Tail x position
r_ty = 0                # Tail y position
r_nmx = 0               # Near motor x position
r_nmy = 0               # Near motor y position
r_fmx = 0               # Far motor x position
r_fmy = 0               # Far motor y position
E_avg = 0               # Energy Average
final_L = 0             # Final Displacements
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

r_tx_arr = []                # Tail x position array
r_ty_arr = []                # Tail y position array
r_nmx_arr = []               # Near motor x position array
r_nmy_arr = []               # Near motor y position array
r_fmx_arr = []               # Far motor x position array
r_fmy_arr = []               # Far motor y position array
E_avg_arr = []               # Energy Average array
final_L_arr = []
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
                        E_avg += dynein.E_total*P

                        r_tx_arr.append(dynein.r_t[0])
                        r_ty_arr.append(dynein.r_t[1])
                        r_nmx_arr.append(dynein.r_nm[0])
                        r_nmy_arr.append(dynein.r_nm[1])
                        r_fmx_arr.append(dynein.r_fm[0])
                        r_fmy_arr.append(dynein.r_fm[1])
                        E_avg_arr.append(dynein.E_total)

                        print('\n\ntrailing',dynein.nba,
                                            new_nma,
                                            new_fma,
                                            dynein.fba,
                                            state)
                        step = run_onebound(dynein.nba,
                                            new_nma,
                                            new_fma,
                                            dynein.fba,
                                            state)
                        print('trailing stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
                        trailing_data[0].append(step['L'])
                        trailing_data[1].append(step['t'])
                        final_L_arr.append(step['L'])
                        ob_t_arr.append(step['t'])

                if np.random.random() < prob_leading:
                        # NEARBOUND State
                        state = 2
                        P_arr.append(P)
                        angles[0].append(nma)
                        angles[1].append(fma)

                        # Calculation for averages
                        r_tx += dynein.r_t[0]*P
                        r_ty += dynein.r_t[1]*P
                        r_nmx += dynein.r_nm[0]*P
                        r_fmx += dynein.r_fm[0]*P
                        E_avg += dynein.E_total*P

                        r_tx_arr.append(dynein.r_t[0])
                        r_ty_arr.append(dynein.r_t[1])
                        r_nmx_arr.append(dynein.r_nm[0])
                        r_nmy_arr.append(dynein.r_nm[1])
                        r_fmx_arr.append(dynein.r_fm[0])
                        r_fmy_arr.append(dynein.r_fm[1])
                        E_avg_arr.append(dynein.E_total)

                        print('\n\nleading', dynein.fba,
                                            new_fma,
                                            new_nma,
                                            dynein.nba,
                                            state)
                        step = run_onebound(dynein.fba,
                                            new_fma,
                                            new_nma,
                                            dynein.nba,
                                            state)
                        print('leading stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
                        leading_data[0].append(step['L'])
                        leading_data[1].append(step['t'])
                        final_L_arr.append(step['L'])
                        ob_t_arr.append(step['t'])


print("FINAL DISPLACEMENTS:", final_L_arr)
for i in range(len(final_L_arr)):
    final_L += final_L_arr[i]*P_arr[i]
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

tx = r_tx/Z          # Tail x array
ty = r_ty/Z          # Tail y array
nmx = r_nmx/Z        # Near motor array
fmx = r_fmx/Z        # Far motor array
E_avg_arr = E_avg/Z  # Average energy array
mean_L = final_L/Z
mean_obt = ob_t/Z

print("Avg Tail x:", tx)
print("Avg Tail y:", ty)
print("Avg nmx:", nmx)
print("Avg fmx:", fmx)
print("Avg E:", E_avg_arr)
print("Avg Final Displacement:", mean_L)
print("Avg ob t:", mean_obt)

fig = plt.figure(1)
gs = gridspec.GridSpec(2, 4)
ax0 = fig.add_subplot(gs[0, 0:2])
ax1 = fig.add_subplot(gs[1, 0:2])
ax2 = fig.add_subplot(gs[0, 2:4])
ax3 = fig.add_subplot(gs[1, 2:4])

ax0.hist(trailing_data[0], bins=50, alpha=0.5, label="Trailing step", normed=True, stacked=True, color="C0")
ax0.hist(leading_data[0], bins=50, alpha=0.5, label="Leading step", normed=True, stacked=True, color="C1")
ax0.legend(loc="upper right")
ax0.set_title("Initial Displacement 8nm")
ax0.set_xlabel("Final Displacement (nm)")
ax0.set_ylabel("Frequency")

ax1.hist(final_L_arr, bins=50, alpha=0.5, normed=True, stacked=True, color="C2")
ax1.legend(loc="upper right")
ax1.set_xlabel("Final Displacement (nm)")
ax1.set_ylabel("Frequency")

ax2.hist(trailing_data[1], bins=50, alpha=0.5, label="Trailing time", normed=False, stacked=True, color="C0")
ax2.hist(leading_data[1], bins=50, alpha=0.5, label="Leading time", normed=False, stacked=True, color="C1")
ax2.legend(loc="upper right")
ax2.set_xlabel("time (s)")
ax2.set_ylabel("Frequency")

ax3.hist(ob_t_arr, bins=50, alpha=0.5, normed=False, stacked=True, color="C2")
ax3.legend(loc="upper right")
ax3.set_xlabel("time (s)")
ax3.set_ylabel("Frequency")

fig1 = plt.figure(2)
gs1 = gridspec.GridSpec(1,1)
ax4 = fig1.add_subplot(gs1[:,:])

ax4.hist(angles[0], bins=50, alpha=0.5, label="nma", normed=True, stacked=True, color="C0")
ax4.hist(angles[1], bins=50, alpha=0.5, label="fma", normed=True, stacked=True, color="C1")
ax4.legend(loc="upper right")
ax4.set_title("Initial Displacement 8nm")
ax4.set_xlabel("Initial Angles")
ax4.set_ylabel("Frequency")

plt.show()
