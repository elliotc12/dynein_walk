import os 
import numpy as np
import argparse
import sys
sys.path.append("../data")
import importlib
import bb_energy_distribution
import time

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=float, help="Array of initial L", default=50)
parser.add_argument("-i", "--i", type=float, help="Increment", default=5)
parser.add_argument("-o", "--o", type=float, help="Starting value", default=1)
parser.add_argument("-k", "--kb", type=float, help="binding const", default=params.for_simulation['k_b'])
parser.add_argument("-s", "--ks", type=float, help="sticky const", default=params.for_simulation['k_stk'])
parser.add_argument("-c", "-cancel", type=bool, help="cancel all jobs", default=False)
args = parser.parse_args()

L = args.L
i = args.i
o = args.o
k_b = args.kb
k_stk = args.ks

init_L = np.arange(o, L, i, dtype=int)

if args.c == True:
    for i in range(len(init_L)):
        os.system('rq cancel -J dynein-mc-{0:.1e}-{1:.1e}-{2}'.format(k_b, k_stk, init_L[i]))
else: 
    for i in range(len(init_L)):
        os.system('rq run -J dynein-mc-{0:.1e}-{1:.1e}-{2} -o L-{3:.1e}-{4:.1e}-{5} python3 monte_carlo_simulation.py -k {6} -s {7} -L {8}'.format(k_b, k_stk, init_L[i], k_b, k_stk, init_L[i], k_b, k_stk, init_L[i]))
        time.sleep(0.1)
