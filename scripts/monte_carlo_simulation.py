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

        if state == 0:
            # NEARBOUND State - Leading step ddata
            step_data['L'].append(step['L'])
            step_data['t'].append(step['t'])

        else:
            # FARBOUND State - Trailing step data
            step_data['L'].append(step['L'])
            step_data['t'].append(step['t'])

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
mc_data_dir = '../data/mc_data_{0:.2e}_{1:.2e}_{2}_{3}_{4}'.format(k_b, k_stk, args.cb, args.cm, args.ct)
if not os.path.exists(mc_data_dir):
    os.mkdir(mc_data_dir)

L = args.L           # Initial Length
N = args.N           # Count
Z = 0                # Partition Function
k = [0]              # Dynein Count
L_err = 0.5          # Margin of error for Length

max_unbinding = 1
b = 1/(params.for_simulation['boltzmann-constant']*params.for_simulation['T'])       # thermodynamic beta from default_parameters.h
eqb_angle = params.for_simulation['eqb']
if bb_energy_distribution.eq_in_degrees:
        eqb_angle = eqb_angle*np.pi/180

# Onebound Data
trailing_data = {   # Trailing step data
        'L': [],
        't': []
}
leading_data = {    # Leading step data
        'L': [],
        't': []
}

# Strings for data file name
t_data_file = '../data/mc_data_{0:.2e}_{1:.2e}_{2}_{3}_{4}/t_{5}_{6}_{7}_{8:.2e}_{9:.2e}_{10}_{11}_{12}_{13}_{14}.txt'.format(k_b, k_stk, args.cb, args.cm, args.ct, int(L), N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C)
l_data_file = '../data/mc_data_{0:.2e}_{1:.2e}_{2}_{3}_{4}/l_{5}_{6}_{7}_{8:.2e}_{9:.2e}_{10}_{11}_{12}_{13}_{14}.txt'.format(k_b, k_stk, args.cb, args.cm, args.ct, int(L), N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C)

seed = 0
np.random.seed(0)

while Z < N:
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
            rate_leading = np.exp(args.C*(dynein.fba - eqb_angle)) # rate_max = np.exp(args.C*(-np.pi (or 0 or np.pi???) - eqb_angle)) if C is negative (check and figure this out!!!!!!!!!)

            prob_trailing = P*rate_trailing     #   Unnormalized (*0.4 in order for prob < 1)
            prob_leading = P*rate_leading       #   Unnormalized (*0.4 in order for prob < 1)
            # print('nba: {}, nma: {}, ta: {}, fma: {}, fba: {}'.format(dynein.nba*57.3, dynein.ob_nma*57.3, dynein.ta*57.3, dynein.ob_fma*57.3, dynein.fba*57.3))
            # print('nba: {}, bb_nma: {}, ta: {}, bb_fma: {}, fba: {}'.format(dynein.nba*57.3, dynein.bb_nma*57.3, dynein.ta*57.3, dynein.bb_fma*57.3, dynein.fba*57.3))
            # print('nb: {}, nm: {}, t: {}, fm: {}, fb: {}'.format(dynein.r_nb, dynein.r_nm, dynein.r_t, dynein.r_fm, dynein.r_fb))
            assert(prob_trailing <= 1),"prob trailing > 1" # if this crashes, we could add a factor to reduce the prob_ to be always less than 1
            assert(prob_leading <= 1),"prob leading > 1"
            if np.random.random() < prob_trailing:
                    # FARBOUND State
                    state = 1

                    collect_onebound_data(k, state, dynein.fba, dynein.ob_fma, dynein.ob_nma, dynein.nba,
                                            L, trailing_data)

            if np.random.random() < prob_leading:
                    # NEARBOUND State
                    state = 0

                    collect_onebound_data(k, state, dynein.nba, dynein.ob_nma, dynein.ob_fma, dynein.fba,
                                            L, leading_data)
                    
            if k[0] % 10 == 0 and k[0]>0:
                    print('Saving data!')
                    np.savetxt(t_data_file, (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
                    np.savetxt(l_data_file, (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
                    if os.path.getsize(t_data_file) > 200000:
                        break
                    if os.path.getsize(l_data_file) > 200000:
                        break


np.savetxt(t_data_file, (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
np.savetxt(l_data_file, (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')

# END OF SIM
