import os
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
                                    str(k_stk),
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

def collect_onebound_data(k, state, bba, bma, uma, uba, L, step_data):
        """
        Call run_onebound function and collect onebound statistics
        """
        # print('\n\nbothbound angles ',bba, bma, uma, uba, state)
        step = run_onebound(bba, bma, uma, uba, state, k[0])

        final_data['L'].append(step['L'])          # Final L array
        final_data['t'].append(step['t'])          # Onebound time array

        if state == 0:
            # NEARBOUND State - Leading step ddata
            # print('leading stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
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


        else:
            # FARBOUND State - Trailing step data
            # print('trailing stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
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
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=1e12)
parser.add_argument("-u", "--kub", type=float, help="Manually set the unbinding const", default=params.for_simulation['k_ub'])
parser.add_argument("-k", "--kb", type=float, help="Manually set the binding const", default=params.for_simulation['k_b'])
parser.add_argument("-s", "--ks", type=float, help="Manually set the sticky const", default=params.for_simulation['k_stk'])
parser.add_argument("-cb", "--cb", type=float, help="Spring constant binding domain", default=params.for_simulation['cb'])
parser.add_argument("-cm", "--cm", type=float, help="Spring constant motor domain", default=params.for_simulation['cm'])
parser.add_argument("-ct", "--ct", type=float, help="Spring constant tail domain", default=params.for_simulation['ct'])
parser.add_argument("--eqb", type=float, help="Binding equilibrium angle", default=params.for_simulation['eqb'])
parser.add_argument("--eqmpre", type=float, help="Motor pre equilibrium angle", default=params.for_simulation['eqmpre'])
parser.add_argument("--eqmpost", type=float, help="Motor post equilibrium angle", default=params.for_simulation['eqmpost'])
parser.add_argument("-t", "--dt", type=float, help="Manually set the dt", default=params.for_simulation['dt'])
parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
args = parser.parse_args()

params.for_simulation['k_ub'] = args.kub
k_b = args.kb        # Binding Rate Constant
k_stk = args.ks      # Sticky Rate Constant
params.for_simulation['cb'] = args.cb
params.for_simulation['cm'] = args.cm
params.for_simulation['ct'] = args.ct
params.for_simulation['eqb'] = args.eqb
params.for_simulation['eqmpre'] = args.eqmpre
params.for_simulation['eqmpost'] = args.eqmpost
dt = args.dt          # Time Step

# Create MC Data Directory if don't exist
mc_data_dir = '../data/mc_data_{0:.2e}_{1:.2e}'.format(k_b, k_stk)
if not os.path.exists(mc_data_dir):
    os.mkdir(mc_data_dir)

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

while Z < N:
        # Making random motor angles
        dynein = bb_energy_distribution.generate_random_bb(L-0.5, L+0.5, params)

        # Checking if energy is nan
        if np.isnan(dynein.E_total) == True:
            continue
        else:
            # Calculating partition function
            P = np.exp(-b*dynein.E_total)
            Z += P

            rate_trailing = np.exp(args.C*(dynein.nba - eqb_angle))
            rate_leading = np.exp(args.C*(dynein.fba - eqb_angle)) # rate_max = np.exp(args.C*(-np.pi (or 0 or np.pi???) - eqb_angle)) if C is negative (check and figure this out!!!!!!!!!)
            max_rate_leading = max(rate_leading, max_rate_leading)
            max_rate_trailing = max(rate_trailing, max_rate_trailing)

            prob_trailing = P*rate_trailing     #   Unnormalized
            prob_leading = P*rate_leading       #   Unnormalized
            print('P: ', P)
            print('Prob trailing: ', prob_trailing)
            print('Prob leading: ', prob_leading)


            new_nma = nma-(np.pi-dynein.nba)
            new_fma = fma-(np.pi-dynein.fba)

            assert(prob_trailing <= 1) # if this crashes, we could add a factor to reduce the prob_ to be always less than 1
            assert(prob_leading <= 1)
            if np.random.random() < prob_trailing: # Maybe should adjust this a tad so it is never > 1.
                    # FARBOUND State
                    state = 1

                    collect_onebound_data(k, state, dynein.fba, new_fma, new_nma, dynein.nba,
                                            L, trailing_data)

                    # plot_bb_before_step(dynein, 'red', 'blue')
                    # plt.savefig('../plots/mc_plots/trailing_{}a_before_step.png'.format(k), transparent=False)

                    # plot_bb_after_step(step['ubx'], step['uby'], step['umx'], step['umy'],
                    #                 step['tx'], step['ty'], step['bmx'], step['bmy'],
                    #                 step['bbx'], step['bby'], 'red', 'blue')
                    # plt.savefig('../plots/mc_plots/trailing_{}b_after_step.png'.format(k), transparent=False)
                    # plt.show()

                    if k[0] % 100 == 0:
                        np.savetxt('../data/mc_data_{0:.2e}_{1:.2e}/t_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk, int(L),
                                    N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C),
                                    (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
                        np.savetxt('../data/mc_data_{0:.2e}_{1:.2e}/l_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk, int(L),
                                    N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C),
                                    (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
                        if os.path.getsize('../data/mc_data_{0:.2e}_{1:.2e}/t_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk, int(L), N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C)) > 700000:
                            break
                        if os.path.getsize('../data/mc_data_{0:.2e}_{1:.2e}/l_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk, int(L), N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C)) > 700000:
                            break
                    print('shtep t')



            if np.random.random() < prob_leading:
                    # NEARBOUND State
                    state = 0

                    collect_onebound_data(k, state, dynein.nba, new_nma, new_fma, dynein.fba,
                                            L, leading_data)

                    # plot_bb_before_step(dynein, 'red', 'blue')
                    # plt.savefig('../plots/mc_plots/leading_{}a_before_step.png'.format(k), transparent=False)

                    # plot_bb_after_step(step['bbx'], step['bby'], step['bmx'], step['bmy'],
                    #                 step['tx'], step['ty'], step['umx'], step['umy'],
                    #                 step['ubx'], step['uby'], 'red', 'blue')
                    # plt.savefig('../plots/mc_plots/leading_{}b_after_step.png'.format(k), transparent=False)
                    # plt.show()
                    if k[0] % 100 == 0:
                        np.savetxt('../data/mc_data_{0:.2e}_{1:.2e}/t_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk, int(L),
                                    N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C),
                                    (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
                        np.savetxt('../data/mc_data_{0:.2e}_{1:.2e}/l_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk, int(L),
                                    N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C),
                                    (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
                        if os.path.getsize('../data/mc_data_{0:.2e}_{1:.2e}/t_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk,int(L), N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C)) > 700000:
                            break
                        if os.path.getsize('../data/mc_data_{0:.2e}_{1:.2e}/l_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk,int(L), N, args.kub, k_b, k_stk,  dt, args.cb, args.cm, args.ct, args.C)) > 700000:
                            break
                    print('shtep l')




# print("FINAL DISPLACEMENTS: {0} \n".format(final_data['L']))
# for i in range(len(final_data['L'])):
#     final_data['L_avg'] += final_data['L'][i]*P_arr[i]
#     final_data['t_avg'] += final_data['t'][i]*P_arr[i]
#     final_data['step_length_avg'] += final_data['step_length'][i]*P_arr[i]


# print("rate_unbinding_leading: ", rate_unbinding_leading)
# print("rate_unbinding_trailing: ", rate_unbinding_trailing)
# print('max_rate_trailing', max_rate_trailing)
# print('max_rate_leading', max_rate_leading)

# Calculating Averages
# final_L_avg = final_data['L_avg']/Z
# step_length_avg = final_data['step_length_avg']/Z
# obt_avg = final_data['t_avg']/Z
#
# print("ONEBOUND AVERAGES")
# print("Avg prob_unbinding:", prob_unbinding_avg)
# print("Avg Final Displacement:", final_L_avg)
# print("Avg Step Length:", step_length_avg)
# print("Avg ob time:", obt_avg)

np.savetxt('../data/mc_data_{0:.2e}_{1:.2e}/t_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk,int(L),
            N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C),
            (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
np.savetxt('../data/mc_data_{0:.2e}_{1:.2e}/l_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk,int(L),
            N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C),
            (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')

# END OF SIM
