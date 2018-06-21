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
    bba_abs = params.eqb
    bma_abs = params.eqb + params.eqmpost - np.pi
    ta_abs = bma_abs + np.pi
    uma_abs = params.eqb + params.eqmpre - np.pi

    Xs[0] = start_x_px
    Xs[1] = Xs[0] + params.ls*np.cos(start_bb_angle+bba_abs)*px_per_nm
    Xs[2] = Xs[1] + params.lt*np.cos(start_bb_angle+bma_abs)*px_per_nm
    Xs[3] = Xs[2] + params.lt*np.cos(start_bb_angle+ta_abs) *px_per_nm
    Xs[4] = Xs[3] + params.ls*np.cos(start_bb_angle+uma_abs)*px_per_nm

    Ys[0] = start_y_px
    Ys[1] = Ys[0] + params.ls*np.sin(start_bb_angle+bba_abs)*px_per_nm
    Ys[2] = Ys[1] + params.lt*np.sin(start_bb_angle+bma_abs)*px_per_nm
    Ys[3] = Ys[2] + params.lt*np.sin(start_bb_angle+ta_abs) *px_per_nm
    Ys[4] = Ys[3] + params.ls*np.sin(start_bb_angle+uma_abs)*px_per_nm

    Rs = np.array([params.radius_b, params.radius_m, params.radius_t, params.radius_m, params.radius_b])*px_per_nm

    plt.figure(fig.number)
    plt.plot(Xs, Ys)
    plt.scatter(Xs, Ys, s=Rs*Rs)

merged_burgess_img = mpimg.imread('papers/paper/figures/model-raw-images/burgess-fig-4-cropped.jpg')

px_per_nm = 4.13

fig = plt.figure()
plt.imshow(merged_burgess_img, origin="lower")
plot_dynein_equilibrium_onebound(fig, 50, 25, px_per_nm, np.pi*1.23)
plt.plot([57, 57+15*px_per_nm], [10, 10])

plt.axis('off')
plt.savefig("plots/burgess-model-figure.pdf", bbox_inches='tight', format="pdf")
