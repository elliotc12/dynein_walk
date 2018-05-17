#!/usr/bin/python3

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os
import glob

import dynein.data as data


parser = argparse.ArgumentParser(description='Script to generate 2 dimensional histogram from dynein stepping data')

parser.add_argument('-d', '--datafile', dest='data_file', action='store', default='data/paper_static_stepping_data-1.txt',
                    help='path to data file', type=str)
parser.add_argument('-b', '--bins', dest='bins', action='store', default=20,
                    help='number of bins for x and y axes', type=int)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
                    help='see prints in console')
parser.add_argument('-s', '--show', dest='show', action='store_true', default=False,
                    help='show graphs in matplotib windows')
# parser.add_argument('-a', '--all', dest='All', action='store_true', default=False,
#                     help='generate plots for all paper_static_stepping_data files')
parser.add_argument('-c', '--colormap', dest='cmap', action='store', type=str,
                    default=None, help='set color map for plots')


args = parser.parse_args()

VERBOSE = args.verbose
DATAFILE = args.data_file
NUM_BINS = args.bins
SHOW = args.show
#ALL = args.All
CMAP = args.cmap

if os.path.exists('color_hist.py'):
    if VERBOSE: print("navigating to root directory")
    os.chdir('../')


if not os.path.isfile(DATAFILE):
    print("Could not find data file. Please specify path using -d")
    exit(1)
        # load in data
if VERBOSE: print("Data file found- loading data...")
Data = data.SteppingData(DATAFILE)
   


#--------------------------------------------------------------------------------------#

def getBinIndex(p, bins):
    for i in range(0, len(bins)-1):
        if p >= bins[i] and p <= bins[i+1]:  # deal with case for p=bins[0], p=bins[-1]
            return i
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
        xbins = np.linspace(0, X.max(), NUM_BINS+1)
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
               xIsTimeValue, yIsTimeValue,filename=None, ):
 
    x_bins, y_bins, counts = getCounts(x, y, xIsTimeValue, yIsTimeValue)
    # print('counts', np.sum(counts), len(x))
    if VERBOSE: print("graphing")
    plt.figure()
    plt.pcolor(x_bins, y_bins, counts, cmap=CMAP)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(x_bins[0], x_bins[-1])
    plt.ylim(y_bins[0], y_bins[-1])
    # print(x_bins)
    plt.title(graph_label)
    cb = plt.colorbar()
    cb.set_label('counts')
    plt.axes().set_aspect('equal')
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

plotCounts(Data.initial_displacements,
           Data.final_displacements,
           "initial disp vs final disp",
           "initial displacement",
           "final displacement",
           xIsTimeValue=False,
           yIsTimeValue=False,
           filename='plots/initial-vs-final{}.pdf'.format(seed_label))

print(len(Data.initial_displacements), len(Data.step_lengths))
plotCounts(Data.initial_displacements,
           Data.step_lengths,
           "initial disp vs step length",
           "initial displacement",
           "step length",
           xIsTimeValue=False,
           yIsTimeValue=False,
           filename='plots/initial-vs-length{}.pdf'.format(seed_label))

if SHOW:
    plt.show()
if VERBOSE: print("finished")
