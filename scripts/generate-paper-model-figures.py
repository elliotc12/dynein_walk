#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess, re

# if 'show' not in sys.argv:
#     matplotlib.use('Agg')

import dynein.draw.cartoon as cartoon

from mpl_toolkits.axes_grid1 import Divider, LocatableAxes, Size

sys.path.insert(0, os.getcwd() + "/data/")
import paper_params as params

import argparse
import datetime
import dynein.draw.cartoon as cartoon
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy import ndimage


def plot_dynein_equilibrium_onebound(fig, start_x_units, start_y_units, units_per_nm, mt_angle):
    Xs = [0, 0, 0, 0, 0]
    Ys = [0, 0, 0, 0, 0]
    bba_abs = mt_angle + params.eqb*np.pi/180.0
    bma_abs = bba_abs - np.pi + params.eqmpost*np.pi/180.0
    ta_abs = bma_abs - np.pi
    uma_abs = bma_abs - params.eqmpre*np.pi/180.0

    Xs[0] = start_x_units
    Xs[1] = Xs[0] + params.ls*np.cos(bba_abs)*units_per_nm
    Xs[2] = Xs[1] + params.lt*np.cos(bma_abs)*units_per_nm
    Xs[3] = Xs[2] + params.lt*np.cos(ta_abs) *units_per_nm
    Xs[4] = Xs[3] + params.ls*np.cos(uma_abs)*units_per_nm

    Ys[0] = start_y_units
    Ys[1] = Ys[0] + params.ls*np.sin(bba_abs)*units_per_nm
    Ys[2] = Ys[1] + params.lt*np.sin(bma_abs)*units_per_nm
    Ys[3] = Ys[2] + params.lt*np.sin(ta_abs) *units_per_nm
    Ys[4] = Ys[3] + params.ls*np.sin(uma_abs)*units_per_nm

    bbox = plt.gca().get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axes_width_inches = bbox.width

    axes_width_units = fig.gca().get_xlim()[1] - fig.gca().get_xlim()[0]
    points_per_axes_units = axes_width_inches / axes_width_units * 72.0

    Rs = np.array([params.radius_b, params.radius_m, params.radius_t, params.radius_m, params.radius_b])*units_per_nm*points_per_axes_units

    plt.figure(fig.number)
    plt.plot(Xs, Ys, c="white", zorder=1)
    plt.plot([Xs[0], Xs[0] + 5*units_per_nm*np.cos(mt_angle)], [Ys[0], Ys[0] + 5*units_per_nm*np.sin(mt_angle)], c="red", linewidth=1, zorder=5)
    plt.scatter(Xs, Ys, c="#aeaae5", s=Rs*Rs, zorder=2, edgecolor='white', alpha=0.3)
    plt.scatter(Xs, Ys, s=Rs*Rs, zorder=3, edgecolor='white', facecolors="none", linewidth=1)

def plot_image(img, org, dpi):
    fig = plt.figure(figsize = (5,5), dpi=dpi)

    imwidth = np.shape(img)[1]
    imheight = np.shape(img)[0]
    if (imwidth > imheight):
        ax = plt.Axes(fig, [0., 0., 1., imheight / imwidth])
    else:
        ax = plt.Axes(fig, [0., 0., imwidth / imheight, 1.])

    ax.set_axis_off()
    fig.add_axes(ax)

    ax.imshow(img, origin=org, aspect='auto') # angles rotate ccw
    return fig

merged_burgess_img = mpimg.imread('papers/paper/figures/model-raw-images/burgess-fig-4-cropped.png')
merged_chowdhury_img = mpimg.imread('papers/paper/figures/model-raw-images/chowdhury-fig-1-cropped.png')
merged_redwine_img = mpimg.imread('papers/paper/figures/model-raw-images/redwine-supplemental-cropped.png')
merged_crystalstruct_img = mpimg.imread('papers/paper/figures/model-raw-images/pymol-cytoplasmic-superimpose.png')

# burgess fig
units_per_nm = 4.13
scalebar_nm = 15
fig = plot_image(merged_burgess_img, "lower", dpi=100)
plt.plot([57, 57+scalebar_nm*units_per_nm], [10, 10])
plt.axis('off')
plot_dynein_equilibrium_onebound(fig, 57, 29, units_per_nm, -np.pi*0.35)
plt.savefig("plots/burgess-model-figure.pdf", format="pdf", interpolation='none', dpi=100, bbox_inches='tight')

# chowdhury fig
units_per_nm = 2.78
scalebar_nm = 25
fig = plot_image(merged_chowdhury_img, "upper", dpi=100)
plt.imshow(merged_chowdhury_img) # angles rotate cw
plot_dynein_equilibrium_onebound(fig, 68, 122, units_per_nm, np.pi)
plt.plot([142, 142+scalebar_nm*units_per_nm], [175, 175])
plt.axis('off')
plt.savefig("plots/chowdhury-model-figure.pdf", bbox_inches='tight', format="pdf", interpolation='none', dpi=100)

# crystal struct fig
units_per_nm = 35
scalebar_nm = 28.8
fig = plot_image(merged_crystalstruct_img, "upper", dpi=600)
plot_dynein_equilibrium_onebound(fig, 1720, 1542, units_per_nm, np.pi*0.64)
plt.plot([692, 692+scalebar_nm*units_per_nm], [1480, 1480])
plt.axis('off')
plt.savefig("plots/crystal-model-figure.pdf", bbox_inches='tight', format="pdf", interpolation='none', dpi=600)
