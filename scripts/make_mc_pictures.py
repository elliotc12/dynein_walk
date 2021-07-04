import os
from os import path
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../data")
import importlib
import argparse

import bb_energy_distribution
import dynein.draw.cartoon as cartoon

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", "--L", type=int, help="Initial L", default=8)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=1e20)
parser.add_argument("-u", "--k_ub", type=float, help="Unbinding const", default=params.for_simulation['k_ub'])
parser.add_argument("-b", "--k_b", type=float, help="Binding const", default=params.for_simulation['k_b'])
parser.add_argument("-s", "--k_stk", type=float, help="Sticky const", default=params.for_simulation['k_stk'])
parser.add_argument("-cb", "--cb", type=float, help="Spring constant binding domain", default=params.for_simulation['cb'])
parser.add_argument("-cm", "--cm", type=float, help="Spring constant motor domain", default=params.for_simulation['cm'])
parser.add_argument("-ct", "--ct", type=float, help="Spring constant tail domain", default=params.for_simulation['ct'])
parser.add_argument("--eqb", type=float, help="Binding equilibrium angle", default=params.for_simulation['eqb'])
parser.add_argument("--eqmpre", type=float, help="Motor pre equilibrium angle", default=params.for_simulation['eqmpre'])
parser.add_argument("--eqmpost", type=float, help="Motor post equilibrium angle", default=params.for_simulation['eqmpost'])
parser.add_argument("-t", "--dt", type=float, help="Time step dt", default=params.for_simulation['dt'])
parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
parser.add_argument("--underMT", action="store_false", help="Plot sims where binding domain can go under MT", default=True)
args = parser.parse_args()

u = ''
if args.underMT == False:
    u = 'u_'

data_dir = '../data/mc_data_{0:.2e}_{1:.2e}_{2}_{3}_{4}_{5}_{6}_{7}_{8}/'.format(args.k_b,
            args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C)

params_extenstion = 'pictures_{0}_{1}_{2:.2e}_{3:.2e}_{4}_{5}_{6}_{7}_{8}_{9}_{10}_{11}.npz'.format(args.L,
                    args.N, args.k_b, args.k_stk, args.dt, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C)

pictures_file = data_dir + u + params_extenstion

pictures_data = np.load(pictures_file, allow_pickle=True)

pictures = pictures_data['pictures'].item()

print(len(pictures['bb_init']))

for i in range(len(pictures['bb_init'])):
    fig = plt.figure()
    if pictures['bb_init'][i][0] < pictures['ub_init'][i][0]:
        # it is a leading step
        leading = 'u'
        trailing = 'b'
    else:
        # it is a trailing step
        leading = 'b'
        trailing = 'u'
    cartoon.dyneinCircles(pictures[trailing+'b_init'][i][0], pictures[trailing+'b_init'][i][1], params.for_simulation['rb'],
                        pictures[trailing+'m_init'][i][0], pictures[trailing+'m_init'][i][1], params.for_simulation['rm'],
                        pictures['t_init'][i][0], pictures['t_init'][i][1], params.for_simulation['rt'],
                        'red', 0.5, fig.gca())
    cartoon.dyneinCircles(pictures[leading+'b_init'][i][0], pictures[leading+'b_init'][i][1], params.for_simulation['rb'],
                        pictures[leading+'m_init'][i][0], pictures[leading+'m_init'][i][1], params.for_simulation['rm'],
                        pictures['t_init'][i][0], pictures['t_init'][i][1], params.for_simulation['rt'],
                        'blue', 0.3, fig.gca())

    cartoon.dyneinCircles(pictures[trailing+'b_final'][i][0], pictures[trailing+'b_final'][i][1], params.for_simulation['rb'],
                        pictures[trailing+'m_final'][i][0], pictures[trailing+'m_final'][i][1], params.for_simulation['rm'],
                        pictures['t_final'][i][0], pictures['t_final'][i][1], params.for_simulation['rt'],
                        'orange', 0.5, fig.gca())
    cartoon.dyneinCircles(pictures[leading+'b_final'][i][0], pictures[leading+'b_final'][i][1], params.for_simulation['rb'],
                        pictures[leading+'m_final'][i][0], pictures[leading+'m_final'][i][1], params.for_simulation['rm'],
                        pictures['t_final'][i][0], pictures['t_final'][i][1], params.for_simulation['rt'],
                        'purple', 0.2, fig.gca())
    plt.arrow(pictures['ub_init'][i][0], 0, pictures['ub_final'][i][0] - pictures['ub_init'][i][0], 0,
                head_width=2, length_includes_head=True)
    plt.gca().set_aspect('equal')

    plt.show()
