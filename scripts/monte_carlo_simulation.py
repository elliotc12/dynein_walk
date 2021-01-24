import os
import numpy as np
import sys
sys.path.append("../data")
import importlib
import argparse
import subprocess
import bb_energy_distribution
import time

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
        assert(exit_code == 0);
        return output_data

def collect_onebound_data(k, state, bba, bma, uma, uba, L, step_data):
        """
        Call run_onebound function and collect onebound statistics
        """
        step = run_onebound(bba, bma, uma, uba, state, k[0])
        step_data['L'].append(step['L'])
        step_data['t'].append(step['t'])

        if k[0] % 50 == 0:
            pictures['bb_final'].append(np.array([step['bbx'], step['bby']]))
            pictures['bm_final'].append(np.array([step['bmx'], step['bmy']]))
            pictures['t_final'].append(np.array([step['tx'], step['ty']]))
            pictures['um_final'].append(np.array([step['umx'], step['umy']]))
            pictures['ub_final'].append(np.array([step['ubx'], step['uby']]))


        k[0]+=1
        sys.stdout.flush()

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=float, help="displacement in nm", default=32)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=1e20)
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
parser.add_argument("--underMT", action="store_false", help="Plot sims where binding domain can go under MT", default=True)

args = parser.parse_args()

k_ub = args.kub
k_b = args.kb        # Binding Rate Constant
k_stk = args.ks      # Sticky Rate Constant
params.for_simulation['cb'] = args.cb
params.for_simulation['cm'] = args.cm
params.for_simulation['ct'] = args.ct
params.for_simulation['eqb'] = args.eqb
params.for_simulation['eqmpre'] = args.eqmpre
params.for_simulation['eqmpost'] = args.eqmpost
dt = args.dt          # Time Step


eqb_angle = params.for_simulation['eqb']
if bb_energy_distribution.eq_in_degrees:
        eqb_angle = eqb_angle*np.pi/180

# Create MC Data Directory if don't exist
mc_data_dir = '../data/mc_data_{0}_{1:.2e}_{2:.2e}_{3}_{4}_{5}_{6}_{7}_{8}_{9}/'.format(k_ub, k_b, k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C)
if not os.path.exists(mc_data_dir):
    os.mkdir(mc_data_dir)

L = args.L           # Initial Length
N = args.N           # Count
Z = 0                # Partition Function
k = [0]              # Dynein Count
L_err = 0.5          # Margin of error for Length
P_factor = 1.0e100

b = 1/(params.for_simulation['boltzmann-constant']*params.for_simulation['T'])       # thermodynamic beta from default_parameters.h

# Strings for data file name
u = ''
if args.underMT == False:
    u = 'u_'

t_data_file = mc_data_dir + u + 't_{0}_{1}_{2}_{3:.2e}_{4:.2e}_{5}_{6}_{7}_{8}_{9}_{10}_{11}_{12}.txt'.format(int(L), N, k_ub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C)
l_data_file = mc_data_dir + u + 'l_{0}_{1}_{2}_{3:.2e}_{4:.2e}_{5}_{6}_{7}_{8}_{9}_{10}_{11}_{12}.txt'.format(int(L), N, k_ub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C)
pictures_data_file = mc_data_dir + u + 'pictures_{0}_{1}_{2}_{3:.2e}_{4:.2e}_{5}_{6}_{7}_{8}_{9}_{10}_{11}_{12}'.format(int(L), N, k_ub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C)

seed = 0
np.random.seed(0)

while Z < N:
        if k[0] == 0:
            # Onebound Data
            trailing_data = {'L': [], 't': []} # Trailing step data
            leading_data = {'L': [], 't': []}  # Leading step data
            pictures = {'bb_init': [], 'bm_init': [], 't_init': [], 'um_init': [], 'ub_init': [],
                        'bb_final': [], 'bm_final': [], 't_final': [], 'um_final': [], 'ub_final': [],}

        # Making random motor angles
        dynein = bb_energy_distribution.generate_random_bb(L-L_err, L+L_err, params)

        # Checking if energy is nan
        if np.isnan(dynein.E_total) == True:
            continue
        else:
            # Calculating partition function
            P = np.exp(-b*dynein.E_total)
            Z += P
            rate_trailing = np.exp(args.C*(dynein.nba - eqb_angle))
            rate_leading = np.exp(args.C*(dynein.fba - eqb_angle))

            prob_trailing = P*rate_trailing*P_factor
            prob_leading = P*rate_leading*P_factor

            # if (len(trailing_data['L'])) == 7:
            #     print('nba: {}, nma: {}, ta: {}, fma: {}, fba: {}'.format(dynein.nba*57.3, dynein.ob_nma*57.3, dynein.ta*57.3, dynein.ob_fma*57.3, dynein.fba*57.3))
            #     print('nba: {}, bb_nma: {}, ta: {}, bb_fma: {}, fba: {}'.format(dynein.nba*57.3, dynein.bb_nma*57.3, dynein.ta*57.3, dynein.bb_fma*57.3, dynein.fba*57.3))
            #     print('nb: {}, nm: {}, t: {}, fm: {}, fb: {}'.format(dynein.r_nb, dynein.r_nm, dynein.r_t, dynein.r_fm, dynein.r_fb))

            assert(prob_trailing > 0 and prob_leading > 0), "Underflow of probability! Need to fix P_factor"
            if prob_trailing > 1 or prob_leading > 1:
                P_factor = P_factor*0.5
                k[0] = 0
                Z = 0
                continue
            if np.random.random() < prob_trailing:
                    # FARBOUND State
                    state = 1
                    if k[0] % 50 == 0:
                        pictures['bb_init'].append(np.array([dynein.r_fb[0]-L, dynein.r_fb[1]]))
                        pictures['bm_init'].append(np.array([dynein.r_fm[0]-L, dynein.r_fm[1]]))
                        pictures['t_init'].append(np.array([dynein.r_t[0]-L, dynein.r_t[1]]))
                        pictures['um_init'].append(np.array([dynein.r_nm[0]-L, dynein.r_nm[1]]))
                        pictures['ub_init'].append(np.array([dynein.r_nb[0]-L, dynein.r_nb[1]]))
                    collect_onebound_data(k, state, dynein.fba, dynein.ob_fma, dynein.ob_nma, dynein.nba,
                                            L, trailing_data)

            if np.random.random() < prob_leading:
                    # NEARBOUND State
                    state = 0
                    if k[0] % 50 == 0:
                        pictures['bb_init'].append(dynein.r_nb)
                        pictures['bm_init'].append(dynein.r_nm)
                        pictures['t_init'].append(dynein.r_t)
                        pictures['um_init'].append(dynein.r_fm)
                        pictures['ub_init'].append(dynein.r_fb)
                    collect_onebound_data(k, state, dynein.nba, dynein.ob_nma, dynein.ob_fma, dynein.fba,
                                            L, leading_data)

            if k[0] % 50 == 0 and k[0]>0:
                    np.savetxt(t_data_file, (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
                    np.savetxt(l_data_file, (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
                    np.savez_compressed(pictures_data_file, pictures=pictures)
                    if os.path.getsize(t_data_file) > 200000:
                        break
                    if os.path.getsize(l_data_file) > 200000:
                        break


np.savetxt(t_data_file, (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
np.savetxt(l_data_file, (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
np.savez_compressed(pictures_data_file, pictures=pictures)


# END OF SIM
