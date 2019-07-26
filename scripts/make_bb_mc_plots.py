import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy.constants
import sys
sys.path.append("../data")
import importlib
import argparse
import subprocess
import bb_energy_distribution

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=float, help="displacement in nm", required=True)
parser.add_argument("-k", "--kb", type=float, help="binding rate", required=True)
parser.add_argument("-t", "--dt", type=float, help="dt", required=True)
args = parser.parse_args()

data_file = open("../data/mc_data_{0}_{1}_{2}.txt".format(int(args.L), args.kb, args.dt), "r")

print(data_file.readline(6000))

data_file.close()
