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
        assert(exit_code == 0);
        return output_data


def collect_onebound_data(k, state, bba, bma, uma, uba, L):
        """
        Call run_onebound function and collect onebound statistics
        """
        print('\n\nbothbound angles ',bba, bma, uma, uba, state)
        step = run_onebound(bba, bma, uma, uba, state, k[0])

        parent_data['t'].append(step['t'])

        if state == 1:
            # FARBOUND State - Trailing step data
            parent_data['init_L'].append(-L)
            print('trailing stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            parent_data['final_L'].append(step['L'])
            parent_data['step_length'].append(step['L']+L)
        else:
            # NEARBOUND State - Leading step ddata
            parent_data['init_L'].append(L)
            print('leading stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            parent_data['final_L'].append(step['L'])
            parent_data['step_length'].append(step['L']-L)

        k[0]+=1


def best_fit(x,y):
    m = (((mean(x)*mean(y)) - mean(x*y))/
        ((mean(x)*mean(x)) - mean(x*x)))
    b = mean(y) - m*mean(x)
    return m, b


params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=int, help="max initial displacement in nm", default=34)
parser.add_argument("-l", "--Ls", type=int, help="intervals in L", default=8)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=100)
parser.add_argument("-k", "--kb", type=float, help="Manually set the binding const", default=params.for_simulation['k_b'])
parser.add_argument("-b", "--cb", type=float, help="Spring constant binding domain", default=params.for_simulation['cb'])
parser.add_argument("-m", "--cm", type=float, help="Spring constant motor domain", default=params.for_simulation['cm'])
parser.add_argument("-s", "--ct", type=float, help="Spring constant tail domain", default=params.for_simulation['ct'])
parser.add_argument("--eqb", type=float, help="Binding equilibrium angle", default=params.for_simulation['eqb'])
parser.add_argument("--eqmpre", type=float, help="Motor pre equilibrium angle", default=params.for_simulation['eqmpre'])
parser.add_argument("--eqmpost", type=float, help="Motor post equilibrium angle", default=params.for_simulation['eqmpost'])
parser.add_argument("-t", "--dt", type=float, help="Manually set the dt", default=params.for_simulation['dt'])
parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
args = parser.parse_args()

k_b = args.kb        # Binding Rate Constant
dt = args.dt          # Time Step
params.for_simulation['cb'] = args.cb
params.for_simulation['cm'] = args.cm
params.for_simulation['ct'] = args.ct
params.for_simulation['eqb'] = args.eqb
params.for_simulation['eqmpre'] = args.eqmpre
params.for_simulation['eqmpost'] = args.eqmpost
parent_data = {'init_L': [], 'final_L': [], 'step_length': [], 't': []}



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
                                                    L)

                    if np.random.random() < prob_leading:
                            # NEARBOUND State
                            state = 0
                            collect_onebound_data(k, state, dynein.nba, new_nma, new_fma, dynein.fba,
                                                    L)



def make_hist(ax, stacked_hist, data, data0, bin, Label, Label0, tof, Color, Color0, Title, xlabel):
    ax.hist(data, bins=bin, alpha=0.5, label=Label, normed=tof, stacked=True, color=Color)
    if stacked_hist == True:
        ax.hist(data0, bins=bin, alpha=0.5, label=Label0, normed=tof, stacked=True, color=Color0)
    ax.legend(loc="upper right")
    ax.set_title(Title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")

def make_hist2d(tof, ax, x_data, y_data, k_b, label):
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
    ax.set_xlabel("Initial Displacement (nm)")
    ax.legend()

fig7, ax = plt.subplots(2,1, figsize=(5,7))
# fig7.tight_layout()
make_hist2d(True, ax[0], parent_data['init_L'], parent_data['final_L'], k_b, "Final L")
ax[0].set_title('Binding Rate: {:.0e}    dt: {}'.format(k_b, dt))
make_hist2d(False, ax[1], parent_data['init_L'], parent_data['step_length'], k_b, "Step Length")
# plt.savefig('../plots/mc_plots/mc_{}_{}_{}_{}_{}_{}_fitting_kb.png'.format(N, dt, k_b, args.cb, args.cm, args.ct), transparent=True)
# plt.savefig('../plots/mc_plots/mc_{}_{}_{}_{}_{}_{}_fitting_kb.svg'.format(N, dt, k_b, args.cb, args.cm, args.ct), transparent=True)
plt.savefig('../plots/mc_plots/mc_{}_{}_{}_{}_{}_{}_{}_{}_{}_fitting_kb.pdf'.format(N, dt, k_b, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost), transparent=True)

# fig8, ax3 = plt.subplots(1,1, figsize=(5,8))
# time_hist = make_hist(ax3, False, parent_data['t'], None, 50,
#                     None, None, False, "C3", None,
#                     "k_b: {0:e}".format(k_b), "time (s)")
# plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_onebound_length_time.pdf'.format(int(L), k_b, dt, N), transparent=False)
