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

def run_onebound(bba, bma, uma, uba, state, k, count):
        """
        Runs onebound.cpp with bb configuration and params.py
        """
        global seed
        seed += 1 # use a different seed every time.  ugh, global variables!
        print('running with inputs', bba, bma, uma, uba, state, k)
        process = subprocess.Popen(['../onebound',
                                    str(k_b[count]),
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
        assert(exit_code == 0);
        return output_data


def collect_onebound_data(k, state, bba, bma, uma, uba, L, count):
        """
        Call run_onebound function and collect onebound statistics
        """
        print('\n\nbothbound angles ',bba, bma, uma, uba, state)
        step = run_onebound(bba, bma, uma, uba, state, k[0], count)

        parent_data[count]['t'].append(step['t'])

        if state == 1:
            # FARBOUND State - Trailing step data
            parent_data[count]['init_L'].append(-L)
            print('trailing stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            parent_data[count]['final_L'].append(step['L'])
            parent_data[count]['step_length'].append(step['L']+L)
        else:
            # NEARBOUND State - Leading step ddata
            parent_data[count]['init_L'].append(L)
            print('leading stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            parent_data[count]['final_L'].append(step['L'])
            parent_data[count]['step_length'].append(step['L']-L)

        k[0]+=1


def best_fit(x,y):
    m = (((mean(x)*mean(y)) - mean(x*y))/
        ((mean(x)*mean(x)) - mean(x*x)))
    b = mean(y) - m*mean(x)
    return m, b


params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=int, help="max initial displacement in nm", default=34)
parser.add_argument("-l", "--dL", type=int, help="intervals in L", default=8)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=100)
parser.add_argument("-cb", "--cb", type=float, help="Spring constant binding domain", default=params.for_simulation['cb'])
parser.add_argument("-cm", "--cm", type=float, help="Spring constant motor domain", default=params.for_simulation['cm'])
parser.add_argument("-ct", "--ct", type=float, help="Spring constant tail domain", default=params.for_simulation['ct'])
parser.add_argument("--eqb", type=float, help="Binding equilibrium angle", default=params.for_simulation['eqb'])
parser.add_argument("--eqmpre", type=float, help="Motor pre equilibrium angle", default=params.for_simulation['eqmpre'])
parser.add_argument("--eqmpost", type=float, help="Motor post equilibrium angle", default=params.for_simulation['eqmpost'])
parser.add_argument("-t", "--dt", type=float, help="Manually set the dt", default=params.for_simulation['dt'])
parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
args = parser.parse_args()

k_b = [3e6, 5.5e6, 8e6, 3e7, 5.5e7, 8e7, 3e8, 5.5e8, 8e8, 3e9, 5.5e9, 8e9]       # Binding Rate Constant
dt = args.dt          # Time Step
params.for_simulation['cb'] = args.cb
params.for_simulation['cm'] = args.cm
params.for_simulation['ct'] = args.ct
params.for_simulation['eqb'] = args.eqb
params.for_simulation['eqmpre'] = args.eqmpre
params.for_simulation['eqmpost'] = args.eqmpost
parent_data = {0:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                1:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                2:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                3:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                4:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                5:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                6:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                7:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                8:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                9:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                10:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                11:{'init_L': [], 'final_L': [], 'step_length': [], 't': []},
                12:{'init_L': [], 'final_L': [], 'step_length': [], 't': []}}
count = 0

for x in k_b:
    if k_b[count] > 1e9:
        dt = 1e-13

    for L in range(1, args.L, args.Ls):
        N = args.N           # Count
        Z = 0                # Partition Function
        k = [0]              # Dynein Count & RNG Seed

        max_unbinding = 1
        b = 11.82733524          # thermodynamic beta from default_parameters.h
        eqb_angle = params.for_simulation['eqb']
        if bb_energy_distribution.eq_in_degrees:
                eqb_angle = eqb_angle*np.pi/180

        seed = 0
        np.random.seed(0)


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
                        prob_trailing = P*rate_trailing     #   Unnormalized
                        prob_leading = P*rate_leading       #   Unnormalized

                        new_nma = nma-(np.pi-dynein.nba)
                        new_fma = fma-(np.pi-dynein.fba)

                        if np.random.random() < prob_trailing: # FIXME need to normalize this a tad so it is never > 1.
                                # FARBOUND State
                                state = 1
                                collect_onebound_data(k, state, dynein.fba, new_fma, new_nma, dynein.nba,
                                                        L, count)

                        if np.random.random() < prob_leading:
                                # NEARBOUND State
                                state = 0
                                collect_onebound_data(k, state, dynein.nba, new_nma, new_fma, dynein.fba,
                                                        L, count)
    count += 1


def make_hist2d(tof, ax, x_data, y_data, k_b, dt, label):
    slope, intercept = best_fit(np.asarray(x_data), np.asarray(y_data))
    rgrsn_line = [(slope*x)+intercept for x in np.asarray(x_data)]
    ax.hist2d(x_data, y_data, bins=(range(-34,34), 60), cmap=plt.cm.jet)
    ax.plot(x_data, rgrsn_line, label='Model: y = ({:.3}) + ({:.3})x'.format(intercept,slope), linestyle=":")
    if tof == True:
        yildiz_line = [(0.6*x)+8.7 for x in np.asarray(x_data)]
        ax.plot(x_data, yildiz_line, label='Experiment: y = ({:.3}) + ({:.3})x'.format(8.7, 0.6), linestyle=":")
    if tof == False:
        yildiz_line = [(-0.4*x)+9.1 for x in np.asarray(x_data)]
        ax.plot(x_data, yildiz_line, label='Experiment: y = ({:.3}) + ({:.3})x'.format(9.1, -0.4), linestyle=":")
    ax.set_ylabel(label)
    ax.set_title('Binding Rate: {:.0e}    dt: {}'.format(k_b, dt))
    ax.legend()

fig7, ax = plt.subplots(4,3, figsize=(12,14))
fig7.tight_layout()

for i in range(12):
    if i < 3:
        make_hist2d(True, ax[0, i], parent_data[i]['init_L'], parent_data[i]['final_L'], k_b[i], 1e-10, "Final L")
    if 3 <= i < 6:
        make_hist2d(True, ax[1, i-3], parent_data[i]['init_L'], parent_data[i]['final_L'], k_b[i], 1e-10, "Final L")
    if 6 <= i < 9:
        make_hist2d(True, ax[2, i-6], parent_data[i]['init_L'], parent_data[i]['final_L'], k_b[i], 1e-10, "Final L")
    if 9 <= i < 12:
        make_hist2d(True, ax[3, i-9], parent_data[i]['init_L'], parent_data[i]['final_L'], k_b[i], 1e-13, "Final L")
plt.savefig('../plots/mc_plots/mc_{}_{}_{}_{}_{}_{}_{}_{}_init_vs_final.pdf'.format(N, dt, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost), transparent=False)

fig8, ax = plt.subplots(4,3, figsize=(12,14))
fig8.tight_layout()
for i in range(12):
    if i < 3:
        make_hist2d(False, ax[0, i], parent_data[i]['init_L'], parent_data[i]['step_length'], k_b[i], 1e-10, "Step Length")
    if 3 <= i < 6:
        make_hist2d(False, ax[1, i-3], parent_data[i]['init_L'], parent_data[i]['step_length'], k_b[i], 1e-10, "Step Length")
    if 6 <= i < 9:
        make_hist2d(False, ax[2, i-6], parent_data[i]['init_L'], parent_data[i]['step_length'], k_b[i], 1e-10, "Step Length")
    if 9 <= i < 12:
        make_hist2d(False, ax[3, i-9], parent_data[i]['init_L'], parent_data[i]['step_length'], k_b[i], 1e-13, "Step Length")
plt.savefig('../plots/mc_plots/mc_{}_{}_{}_{}_{}_{}_{}_{}_init_vs_step_length.pdf'.format(N, dt, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost), transparent=False)
