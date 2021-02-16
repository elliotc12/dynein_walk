import os
from os import path
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import sys
sys.path.append("../data")
import importlib
import argparse
import subprocess
from glob import glob

import bb_energy_distribution
import dynein.draw.cartoon as cartoon

from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter

def animate(i):
    ax.cla()
    ax.set(xlim = (-50,50), ylim = (0,80))
    cartoon.dyneinCircles(bb[0][i], bb[1][i], params.for_simulation['rb'],
                            bm[0][i], bm[1][i], params.for_simulation['rm'],
                            t[0][i], t[1][i], params.for_simulation['rt'],
                            'red', 0.5, fig.gca())

    cartoon.dyneinCircles(ub[0][i], ub[1][i], params.for_simulation['rb'],
                            um[0][i], um[1][i], params.for_simulation['rm'],
                            t[0][i], t[1][i], params.for_simulation['rt'],
                            'blue', 0.3, fig.gca())

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-k", "--k_b", type=float, help="Binding const", default=params.for_simulation['k_b'])
parser.add_argument("-s", "--k_stk", type=float, help="Sticky const", default=params.for_simulation['k_stk'])
parser.add_argument("-cb", "--cb", type=float, help="Spring constant binding domain", default=params.for_simulation['cb'])
parser.add_argument("-cm", "--cm", type=float, help="Spring constant motor domain", default=params.for_simulation['cm'])
parser.add_argument("-ct", "--ct", type=float, help="Spring constant tail domain", default=params.for_simulation['ct'])
parser.add_argument("--eqb", type=float, help="Binding equilibrium angle", default=params.for_simulation['eqb'])
parser.add_argument("--eqmpre", type=float, help="Motor pre equilibrium angle", default=params.for_simulation['eqmpre'])
parser.add_argument("--eqmpost", type=float, help="Motor post equilibrium angle", default=params.for_simulation['eqmpost'])
parser.add_argument("--underMT", action="store_false", help="Plot sims where binding domain can go under MT", default=True)
args = parser.parse_args()

u = ''
if args.underMT == False:
    u = 'u_'


movie_data = '../data/mc_movie_data_5.5e+09_1e+08_0_1_1_120_197_242.txt'
# movie_file = '../data/mc_movie_data_{0:.1e}_{1:.0e}_{2:.0n}_{3:.0n}_{4:.0n}_{5:.0n}_{6:.0n}_{7:.0n}.txt'.format(args.k_b,
#                 args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost)

time = np.loadtxt(movie_data)[:,1]
bb = np.array([np.loadtxt(movie_data)[:,7], np.loadtxt(movie_data)[:,8]])
bm = np.array([np.loadtxt(movie_data)[:,9], np.loadtxt(movie_data)[:,10]])
t = np.array([np.loadtxt(movie_data)[:,11], np.loadtxt(movie_data)[:,12]])
um = np.array([np.loadtxt(movie_data)[:,13], np.loadtxt(movie_data)[:,14]])
ub = np.array([np.loadtxt(movie_data)[:,15], np.loadtxt(movie_data)[:,16]])

fig, ax = plt.subplots(figsize = (10,8))
ax.set(xlim = (-50,50), ylim = (0,80))
cartoon.dyneinCircles(bb[0][0], bb[1][0], params.for_simulation['rb'],
                            bm[0][0], bm[1][0], params.for_simulation['rm'],
                            t[0][0], t[1][0], params.for_simulation['rt'],
                            'red', 0.5, fig.gca())
cartoon.dyneinCircles(ub[0][0], ub[1][0], params.for_simulation['rb'],
                            um[0][0], um[1][0], params.for_simulation['rm'],
                            t[0][0], t[1][0], params.for_simulation['rt'],
                            'blue', 0.3, fig.gca())

movie = FuncAnimation(fig, animate, frames = len(time), interval = 50)


# plt.show()

movie_save_fname =r"../plots/mc_plots/mc_movie.mp4"
writervideo = FFMpegWriter(fps=10)

movie.save(movie_save_fname, writer = writervideo)
