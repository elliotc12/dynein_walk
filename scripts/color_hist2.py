#!/usr/bin/python3

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os
import glob

import dynein.data as data

plt.rcParams.update({'font.size': 14})

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

data_files = []
for fname in os.listdir(args.data_dir):
    if os.path.isfile(args.data_dir + fname):
        if (DATAWC in fname and ".txt" in fname and "~" not in fname and ("stepping_data" in fname or args.data_dir != "data/")):
            data_files.append(args.data_dir + fname)

if len(data_files) == 0:
    print("No files matching wildcard " + args.data_wc)
    exit(1)

# load in data
if VERBOSE: print("Data file found- loading data...")
datasets = []
for df in data_files:
    datasets.append(data.SteppingData(df))


#--------------------------------------------------------------------------------------#

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


def plotCounts(x, y, graph_label, x_label, y_label,
               xIsTimeValue, yIsTimeValue,filename=None, drawline=False):

    x_bins, y_bins, counts = getCounts(x, y, xIsTimeValue, yIsTimeValue)
    if VERBOSE: print("graphing")
    plt.figure()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(x_bins[0], x_bins[-1])
    plt.ylim(y_bins[0], y_bins[-1])
    plt.title(graph_label)
    # plt.axes().set_aspect('equal')

    plt.pcolor(x_bins, y_bins, counts, cmap=CMAP)
    cb = plt.colorbar()
    cb.set_label('counts')

    if xIsTimeValue:
        plt.gca().set_xscale('log')

    if drawline:
        A = np.vstack([x, np.ones(len(x))]).T
        m, c = np.linalg.lstsq(A, y)[0]
        eq = "Y = {:.2}x + {:.2}".format(m, c)
        plt.plot([x_bins[0], x_bins[-1]], [x_bins[0]*m, x_bins[-1]*m]+c, label=eq, linestyle=":")
        plt.legend()

    if filename is None:
        filename = 'plots/'+graph_label.replace(' ',  '_')+".pdf"
    plt.savefig(filename)


seed_label = ''
#if ALL:
#    seed_label = '-multiple-seeds'

# plotCounts(Data.step_times,
#            Data.step_lengths,
#            "total step time vs step length",
#            'step time',
#            'step length',
#            xIsTimeValue=True,
#            yIsTimeValue=False,
#            filename='plots/time-vs-length{}.pdf'.format(seed_label))

initial_displacements = np.concatenate([s.initial_displacements for s in datasets])
final_displacements = np.concatenate([s.final_displacements for s in datasets])
onebound_times = np.concatenate([s.onebound_times for s in datasets])
bothbound_times = np.concatenate([s.bothbound_times for s in datasets])
step_lengths = np.concatenate([s.step_lengths for s in datasets])

plotCounts(initial_displacements,
           final_displacements,
           "",
           "Initial displacement (nm)",
           "Final displacement (nm)",
           xIsTimeValue=False,
           yIsTimeValue=False,
           filename='plots/initial-vs-final-displacement{}.pdf'.format(seed_label))

plotCounts(onebound_times,
           step_lengths,
           "",
           "Onebound time (s)",
           "Step length (nm)",
           xIsTimeValue=True,
           yIsTimeValue=False,
           filename='plots/onebound-time-vs-step-length{}.pdf'.format(seed_label))

plotCounts(bothbound_times,
           step_lengths,
           "",
           "Bothbound time (s)",
           "Step length (nm)",
           xIsTimeValue=True,
           yIsTimeValue=False,
           filename='plots/bothbound-time-vs-step-length{}.pdf'.format(seed_label))


plotCounts(initial_displacements,
           step_lengths,
           "",
           "initial displacement",
           "step length",
           xIsTimeValue=False,
           yIsTimeValue=False,
           drawline=True,
           filename='plots/initial-displacement-vs-step-length{}.pdf'.format(seed_label))

if SHOW:
    plt.show()
if VERBOSE: print("finished")
