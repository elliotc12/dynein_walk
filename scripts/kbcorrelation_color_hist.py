#!/usr/bin/python3


import numpy as np
import matplotlib
import sys

if 'show' not in sys.argv:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import argparse
import os
import glob

import dynein.data as data

parser = argparse.ArgumentParser(description='Script to generate 2 dimensional histogram from dynein stepping data')

parser.add_argument('-d', '--datawildcard', dest='data_wc', action='store', default='paper_main',
                    help='data file wildcard', type=str)
parser.add_argument('-r', '--datadir', dest='data_dir', action='store', default='data/',
                    help='data file src directory', type=str)
parser.add_argument('-b', '--bins', dest='bins', action='store', default=50,
                    help='number of bins for x and y axes', type=int)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
                    help='see prints in console')
parser.add_argument('-s', '--show', dest='show', action='store_true', default=False,
                    help='show graphs in matplotib windows')
# parser.add_argument('-a', '--all', dest='All', action='store_true', default=False,
#                     help='generate plots for all paper_exponential_stepping_data files')
parser.add_argument('-c', '--colormap', dest='cmap', action='store', type=str,
                    default=None, help='set color map for plots')
parser.add_argument('-k', '--kbcor', dest='kbcor', action='store', type=bool,
                    default=False, help='kbcor mode')


args = parser.parse_args()

VERBOSE = args.verbose
DATAWC = args.data_wc
NUM_BINS = args.bins
SHOW = args.show
#ALL = args.All
CMAP = args.cmap

if os.path.exists('color_hist2.py'):
    if VERBOSE: print("navigating to root directory")
    os.chdir('../')

def getBinIndex(p, bins):
    if p <= bins[0]:
        return 0
    elif p >= bins[len(bins)-1]:
        return len(bins)-1
    for i in range(0, len(bins)-1):
        if p > bins[i] and p <= bins[i+1]:  # deal with case for p=bins[0], p=bins[-1]
            return i
    print("Couldn't put the following value in a bin: ", p)
    print("Bins: ", bins)
    assert(False)  # throw exception if we can't find a bin for a value


def getCounts(X, Y, xIsTimeValue, yIsTimeValue):
    # insure that times start at zero and not the min time
    X = np.asarray(X)
    Y = np.asarray(Y)
    if VERBOSE: print(xIsTimeValue, yIsTimeValue) 

    if VERBOSE: print(X.shape, Y.shape)
    assert X.shape == Y.shape

    if VERBOSE: print("binning data")

    if xIsTimeValue is True:
        xbins = np.logspace(-10, np.log10(X.max()), NUM_BINS+1)

    else:
        xbins = np.linspace(X.min(), X.max(), NUM_BINS+1)

    if yIsTimeValue is True:
        ybins = np.linspace(0, Y.max(), NUM_BINS+1)
    else:
        ybins = np.linspace(Y.min(), Y.max(), NUM_BINS+1)

    if VERBOSE: print("counting")
    counts = np.zeros((len(xbins),len(ybins)))

    for k in range(0, len(X)):
        i = getBinIndex(X[k], xbins)
        j = getBinIndex(Y[k], ybins)
        counts[j, i] += 1  # rows are y values, columns are x values
    return xbins, ybins, counts

def plotkbcor():
    plt.rcParams.update({'font.size': 6})
    kb_datafiles = []
    for fname in os.listdir('data/'):
        if os.path.isfile('data/' + fname):
            if ("kbcor" in fname and ".txt" in fname and "~" not in fname and "stepping_data" in fname):
                kb_datafiles.append('data/' + fname)

    if len(kb_datafiles) == 0:
        print("No files matching wildcard kbcor")
        exit(1)

    datasets = {}
    for df in kb_datafiles:
        kb = df[34:34+df[34:].index('.')]

        if kb not in datasets.keys():
            datasets[kb] = []

        datasets[kb].append(data.SteppingData(df))

    plt.close('all')
    fig, ax = plt.subplots(5, 3, figsize=(10, 15))
    ax = ax.flatten()

    for (i, kb) in enumerate([str(s) for s in sorted([int(k) for k in list(datasets.keys())])]):
        initial_displacements = np.concatenate([s.initial_displacements for s in datasets[kb]])
        step_lengths = np.concatenate([s.step_lengths for s in datasets[kb]])

        x_bins, y_bins, counts = getCounts(initial_displacements, step_lengths, False, False)

        # ax[i].set_xlabel("Initial displacement (nm)")
        # ax[i].set_ylabel("Final displacement (nm)")
        ax[i].set_xlim(x_bins[0], x_bins[-1])
        ax[i].set_ylim(y_bins[0], y_bins[-1])
        ax[i].set_title(str(kb))

        ax[i].pcolor(x_bins, y_bins, counts, cmap=CMAP)
        # cb = ax[i].colorbar()
        # cb.set_label('counts')

        A = np.vstack([initial_displacements, np.ones(len(initial_displacements))]).T
        m, c = np.linalg.lstsq(A, step_lengths)[0]
        eq = "Model: y = {:.2} + {:.2}x".format(c, m)
        if m < 0:
            eq = "Model: y = {:.2} - {:.2}x".format(c, -m)
        ax[i].plot([x_bins[0], x_bins[-1]], [x_bins[0]*m, x_bins[-1]*m]+c, label=eq, linestyle=":", color='C1')

        m = -0.4 # yildiz 2012
        c = 9.05 # yildiz 2012
        eq = "Experiment: y = {:.2} - {:.2}x".format(c, -m)
        ax[i].plot(np.array([x_bins[0], x_bins[-1]]), np.array([x_bins[0]*m, x_bins[-1]*m])+c, label=eq, linestyle=":", color='C0')
        ax[i].legend()

    plt.tight_layout()
    filename = "plots/kbcorplot.pdf"
    plt.savefig(filename)

plotkbcor()
