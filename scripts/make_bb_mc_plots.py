import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
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

def rad_to_deg(angle):
    # array = angle*180/np.pi
    return angle

data_file = open("../data/mc_data_{0}_{1}_{2}.txt".format(int(args.L), args.kb, args.dt), "r")

data = data_file.readlines()
init_ang = [[] for i in range(2)]
final_L = []
final_ang = [[] for i in range(2)]
step_length = []

for i in range(3, 100):
    line = data[i].split("\t")
    init_ang[0].append(float(line[1]))
    init_ang[1].append(float(line[2]))
    final_L.append(float(line[3]))
    final_ang[0].append(float(line[4]))
    final_ang[1].append(float(line[5]))
    step_length.append(float(line[6]))

print("nma angles: ", init_ang[0])
print("fma angles: ", init_ang[1])
print("final L: ", final_L)


fig = plt.figure(figsize=(10,15))

# make contourf graph
ax1 = fig.add_subplot(111, projection='3d')
ax1.scatter(init_ang[0], init_ang[1], final_L)
# contour = ax1.contour(rad_to_deg(init_ang[0]), rad_to_deg(init_ang[1]), final_L, np.linspace(0, 1, 5), colors='w', linewidth=10)
ax1.set_xlabel(r'$\theta_{nm}$', fontsize=60)
ax1.set_ylabel(r'$\theta_{fm}$', fontsize=60)
ax1.set_zlabel(r'Final L')
# ax1.set_xticks(np.linspace(0, 360, 13))
# ax1.set_yticks(np.linspace(0, 360, 13))
# cb = plt.colorbar(FinalLPlot)
# cb.set_label(r"Final L", fontsize=40)
# cb.set_ticks(np.linspace(0, 1, 5))
# cb.add_lines(contour)

plt.show()

data_file.close()
