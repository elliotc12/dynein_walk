#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess, re

if 'show' not in sys.argv:
    matplotlib.use('Agg')

import dynein.draw.cartoon as cartoon

sys.path.insert(0, os.getcwd() + "/data/")
import paper_params as params

import argparse
import datetime
import dynein.draw.cartoon as cartoon
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy import ndimage


def plot_dynein_equilibrium_onebound(fig, start_x_px, start_y_px, px_per_nm, start_bb_angle):
    print(px_per_nm)
    Xs = [0, 0, 0, 0, 0]
    Ys = [0, 0, 0, 0, 0]
    bba_abs = start_bb_angle + params.eqb*np.pi/180.0
    bma_abs = bba_abs - np.pi + params.eqmpost*np.pi/180.0
    ta_abs = bma_abs - np.pi
    uma_abs = bma_abs - params.eqmpre*np.pi/180.0

    Xs[0] = start_x_px
    Xs[1] = Xs[0] + params.ls*np.cos(bba_abs)*px_per_nm
    Xs[2] = Xs[1] + params.lt*np.cos(bma_abs)*px_per_nm
    Xs[3] = Xs[2] + params.lt*np.cos(ta_abs) *px_per_nm
    Xs[4] = Xs[3] + params.ls*np.cos(uma_abs)*px_per_nm

    Ys[0] = start_y_px
    Ys[1] = Ys[0] + params.ls*np.sin(bba_abs)*px_per_nm
    Ys[2] = Ys[1] + params.lt*np.sin(bma_abs)*px_per_nm
    Ys[3] = Ys[2] + params.lt*np.sin(ta_abs) *px_per_nm
    Ys[4] = Ys[3] + params.ls*np.sin(uma_abs)*px_per_nm

    Rs = np.array([params.radius_b, params.radius_m, params.radius_t, params.radius_m, params.radius_b])*px_per_nm

    plt.figure(fig.number)
    plt.plot(Xs, Ys, c="black", zorder=1)
    plt.scatter(Xs, Ys, s=Rs*Rs, zorder=2)

merged_burgess_img = mpimg.imread('papers/paper/figures/model-raw-images/burgess-fig-4-cropped.png')
merged_chowdhury_img = mpimg.imread('papers/paper/figures/model-raw-images/chowdhury-fig-1-cropped.png')
merged_redwine_img = mpimg.imread('papers/paper/figures/model-raw-images/redwine-supplemental-cropped.png')

px_per_nm = 4.13

fig = plt.figure()
plt.imshow(merged_burgess_img, origin="lower") # angles rotate ccw
plot_dynein_equilibrium_onebound(fig, 57, 29, px_per_nm, 60*np.pi/180.0-params.eqb*np.pi/180.0)
plt.plot([57, 57+15*px_per_nm], [10, 10])

plt.axis('off')
plt.savefig("plots/burgess-model-figure.pdf", bbox_inches='tight', format="pdf")




px_per_nm = 4.13

fig = plt.figure()
plt.imshow(merged_chowdhury_img) # angles rotate cw
plot_dynein_equilibrium_onebound(fig, 68, 122, px_per_nm, np.pi)

plt.axis('off')
plt.savefig("plots/chowdhury-model-figure.pdf", bbox_inches='tight', format="pdf")




# px_per_nm = 4.13

# fig = plt.figure()
# plt.imshow(merged_redwine_img)
# plot_dynein_equilibrium_onebound(fig, 68, 122, px_per_nm, 60*np.pi/180.0)

# plt.axis('off')
# plt.savefig("plots/redwine-model-figure.pdf", bbox_inches='tight', format="pdf")
