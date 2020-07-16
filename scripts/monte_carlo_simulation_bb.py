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
        k[0]+=1



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

max_rate_trailing = 0
max_rate_leading = 0

# Bothbound Data
P_arr = []
angles = {'nma': [],'fma': []}                       # Pair of Angles
r_t = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}     # Tail data
r_nm = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}    # Near motor data
r_fm = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}    # Far motor data
E = {'bb': [], 'avg': 0}                             # Energy data
rate_unbinding = {                                   # Unbinding rates
        'trailing': [],
        'leading': []
}
prob_unbinding = {                                   # Unbinding probabilities
        'trailing': [],
        'leading': [],
        'unbinding': [],
        'avg': 0
}


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
                        collect_bothbound_data(k, dynein, P, state, nma, fma, prob_trailing)

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


                if np.random.random() < prob_leading:
                        # NEARBOUND State
                        state = 0
                        collect_bothbound_data(k, dynein, P, state, nma, fma, prob_leading)

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



# What to collect and output or visualize?

### Bothbound data
# Mean angles while bothbound? (no stepping required)
# Mean motor/tail domain locations

# print("rate_unbinding_leading: ", rate_unbinding_leading)
# print("rate_unbinding_trailing: ", rate_unbinding_trailing)
# print('max_rate_trailing', max_rate_trailing)
# print('max_rate_leading', max_rate_leading)

# BOTHBOUND AVERAGES
# tx_avg = r_t['x_avg']/Z          # Tail x
# ty_avg = r_t['y_avg']/Z          # Tail y
# nmx_avg = r_nm['x_avg']/Z        # Near motor x
# nmy_avg = r_nm['y_avg']/Z        # Near Motor y
# fmx_avg = r_fm['x_avg']/Z        # Far motor x
# fmy_avg = r_fm['y_avg']/Z        # Far Motor y
# E_avg = E['avg']/Z          # Average energy
# prob_unbinding_avg = prob_unbinding['avg']/Z


# print("BOTHBOUND AVERAGES")
# print("Avg Tail x:", tx_avg)
# print("Avg Tail y:", ty_avg)
# print("Avg nmx:", nmx_avg)
# print("Avg nmy:", nmy_avg)
# print("Avg fmx:", fmx_avg)
# print("Avg fmy:", fmy_avg)
# print("Avg E:", E_avg)
# print("Avg prob_unbinding:", prob_unbinding_avg)


np.savetxt('../data/mc_data_{0:.2e}_{1:.2e}/t_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk,int(L),
            N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C),
            (trailing_data['L'], trailing_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')
np.savetxt('../data/mc_data_{0:.2e}_{1:.2e}/l_{2}_{3}_{4}_{5:.2e}_{6:.2e}_{7}_{8}_{9}_{10}_{11}.txt'.format(k_b, k_stk,int(L),
            N, args.kub, k_b, k_stk, dt, args.cb, args.cm, args.ct, args.C),
            (leading_data['L'], leading_data['t']), fmt='%.6e', delimiter=' ', newline='\n\n')

# END OF SIM
