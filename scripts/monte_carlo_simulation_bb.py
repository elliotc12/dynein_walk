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

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=float, help="displacement in nm", default=50)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=1e4)
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

dL = 1.0 # 1 nm resolution ???
L_arr = np.arange(dL, args.L + dL/2, dL)               # All initial lengths
rate_leading = np.zeros_like(L_arr)
rate_trailing = np.zeros_like(rate_leading)

start_time = time.time()
Z = np.zeros_like(L_arr) # Partition Function
Ndata = np.zeros_like(Z)

b = 11.82733524          # thermodynamic beta from default_parameters.h
eqb_angle = params.for_simulation['eqb']
if bb_energy_distribution.eq_in_degrees:
        eqb_angle = eqb_angle*np.pi/180

seed = 0
np.random.seed(0)

while Z.min() < args.N:
        # Making random motor angles
        dynein = bb_energy_distribution.generate_random_bb_any_L(params)

        # Checking if energy is nan
        if np.isnan(dynein.E_total) == True:
            continue
        else:
            # Calculating partition function
            L = dynein.L # this is the actual length
            i = int(L/dL+ 0.5) # FIXME check that this goes in the right places!!!
            if i >= len(L_arr):
                continue
            P = np.exp(-b*dynein.E_total)
            Z[i] += P
            Ndata[i] += 1

            this_rate_trailing = np.exp(args.C*(dynein.nba - eqb_angle))
            this_rate_leading = np.exp(args.C*(dynein.fba - eqb_angle))

            rate_trailing[i] += P*this_rate_trailing     #   Not yet normalized
            rate_leading[i] += P*this_rate_leading       #   Not yet normalized
            print('at L={} we are {} done'.format(L, Z.min()/args.N))

rate_leading /= Z # Normalize our average, but we're still missing the unbinding rate factor
rate_trailing /= Z
print('saving to', bbdatapath)
np.savez_compressed(bbdatapath, L=L_arr, rate_leading=rate_leading, rate_trailing=rate_trailing)
print('TIME: {}s for N = {}'.format(time.time()-start_time, args.N))

np.savez_compressed(bbdatapath, L=L_arr, rate_leading=rate_leading, rate_trailing=rate_trailing)
