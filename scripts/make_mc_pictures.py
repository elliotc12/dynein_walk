import os
from os import path
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
from glob import glob

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=float, help="Initial L", default=8)
parser.add_argument("-u", "--kub", type=float, help="Unbinding const", default=params.for_simulation['k_ub'])
parser.add_argument("-k", "--kb", type=float, help="Binding const", default=params.for_simulation['k_b'])
parser.add_argument("-s", "--ks", type=float, help="Sticky const", default=params.for_simulation['k_stk'])
parser.add_argument("-cb", "--cb", type=float, help="Spring constant binding domain", default=params.for_simulation['cb'])
parser.add_argument("-cm", "--cm", type=float, help="Spring constant motor domain", default=params.for_simulation['cm'])
parser.add_argument("-ct", "--ct", type=float, help="Spring constant tail domain", default=params.for_simulation['ct'])
parser.add_argument("--eqb", type=float, help="Binding equilibrium angle", default=params.for_simulation['eqb'])
parser.add_argument("--eqmpre", type=float, help="Motor pre equilibrium angle", default=params.for_simulation['eqmpre'])
parser.add_argument("--eqmpost", type=float, help="Motor post equilibrium angle", default=params.for_simulation['eqmpost'])
parser.add_argument("-t", "--dt", type=float, help="Time step dt", default=params.for_simulation['dt'])
args = parser.parse_args()

basepath = '../data/mc_data_{0}_{1:.2e}_{2:.2e}_{3}_{4}_{5}_{6}_{7}_{8}/'.format(args.k_ub, args.k_b, args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost)

pictures_files = glob('{}/pictures_*.txt'.format(basepath))
