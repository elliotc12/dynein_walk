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
import time

"""
Monte Carlo simulation for dynein taking a step
"""


def collect_bothbound_data(k, self, P, nma, fma, prob):
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

        # Storing data for histograms
        r_t['x'].append(self.r_t[0])
        r_t['y'].append(self.r_t[1])
        r_nm['x'].append(self.r_nm[0])
        r_nm['y'].append(self.r_nm[1])
        r_fm['x'].append(self.r_fm[0])
        r_fm['y'].append(self.r_fm[1])
        E['bb'].append(self.E_total)

        k[0]+=1



params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=float, help="displacement in nm", default=np.arange(1,51,1))
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=1e9)
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
mc_bb_data_dir = '../data/mc_bb_data/'
if not os.path.exists(mc_bb_data_dir):
    os.mkdir(mc_bb_data_dir)
bbdatapath = mc_bb_data_dir + 'bb_{0:.2e}_{1:.2e}'.format(k_b, k_stk)

L_arr = np.array(args.L)               # All initial lengths
rate_leading = np.zeros(len(L_arr))
rate_trailing = np.zeros_like(rate_leading)

for i in range(len(L_arr)):
    start_time = time.time()
    L = L_arr[i]         # Initial Length
    print(L)
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

                    this_rate_trailing = np.exp(args.C*(dynein.nba - eqb_angle))
                    this_rate_leading = np.exp(args.C*(dynein.fba - eqb_angle))

                    max_rate_leading = max(this_rate_leading, max_rate_leading)
                    max_rate_trailing = max(this_rate_trailing, max_rate_trailing)

                    rate_trailing[i] += P*this_rate_trailing     #   Not yet normalized
                    rate_leading[i] += P*this_rate_leading       #   Not yet normalized

                    new_nma = nma-(np.pi-dynein.nba)
                    new_fma = fma-(np.pi-dynein.fba)

                #     if np.random.random() < prob_trailing: # Should normalize this a tad so it is never > 1.
                #             # FARBOUND State
                #             state = 1
                #             collect_bothbound_data(k, dynein, P, nma, fma, prob_trailing)



                #     if np.random.random() < prob_leading:
                #             # NEARBOUND State
                #             state = 0
                #             collect_bothbound_data(k, dynein, P, nma, fma, prob_leading)

    rate_leading[i] /= Z # Normalize our average, but we're still missing the unbinding rate factor
    rate_trailing[i] /= Z
    print('saving to', bbdatapath)
    np.savez_compressed(bbdatapath, L=L_arr, rate_leading=rate_leading, rate_trailing=rate_trailing)
    print('TIME: {}s for N = {}'.format(time.time()-start_time), N)

np.savez_compressed(bbdatapath, L=L_arr, rate_leading=rate_leading, rate_trailing=rate_trailing)

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


# END OF SIM
