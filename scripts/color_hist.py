#!/usr/bin/python3

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os
import glob


parser = argparse.ArgumentParser(description='Script to generate 2 dimensional histogram from dynein stepping data')

parser.add_argument('-d', '--datafile', dest='data_file', action='store', default='data/paper_static_stepping_data-1.txt',
                    help='path to data file', type=str)
parser.add_argument('-b', '--bins', dest='bins', action='store', default=20,
                    help='number of bins for x and y axes', type=int)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
                    help='see prints in console')
parser.add_argument('-s', '--show', dest='show', action='store_true', default=False,
                    help='show graphs in matplotib windows')
parser.add_argument('-a', '--all', dest='All', action='store_true', default=False,
                    help='generate plots for all paper_static_stepping_data files')
parser.add_argument('-c', '--colormap', dest='cmap', action='store', type=str,
                    default=None, help='set color map for plots')


args = parser.parse_args()

VERBOSE = args.verbose
DATAFILE = args.data_file
NUM_BINS = args.bins
SHOW = args.show
ALL = args.All
CMAP = args.cmap

if os.path.exists('color_hist.py'):
    if VERBOSE: print("navigating to root directory")
    os.chdir('../')

step_times = []
onebound_times = []
bothbound_times = []
step_lengths = []
initial_displacements = []
final_displacements = []

if not ALL:
    if not os.path.isfile(DATAFILE):
        print("Could not find data file. Please specify path using -d")
        exit(1)
        # load in data
    if VERBOSE: print("Data file found- loading data...")
    data = np.loadtxt(DATAFILE)
    if VERBOSE: print("Data loaded- formatting...")
    bind_times = np.array(data[:,1])
    unbind_times = np.array(data[:,0])
    near_positions = np.around(np.array(data[:,2]), decimals=7)  #need to figure out why fixing number of decimals is necessary
    far_positions = np.around(np.array(data[:,3]), decimals=7)
    near_step_lens = near_positions[1:] - near_positions[:-1]  #reduces total length by one. Will include 0 step lengths
    far_step_lens = far_positions[1:] - far_positions[:-1]
    
    assert (len(near_positions) == len(far_positions))
    for s in range(1, len(near_positions)):
        assert((near_positions[s-1] == near_positions[s]) or (far_positions[s-1] == far_positions[s]))
        if (near_positions[s-1] == near_positions[s]) and (far_positions[s-1] == far_positions[s]):
            continue
        if not (near_positions[s-1] == near_positions[s]):  # i.e. if near foot stepped 
            initial_displacements.append(near_positions[s-1]-far_positions[s-1]) 
            final_displacements.append(near_positions[s]-far_positions[s])
        elif not (far_positions[s-1] == far_positions[s]):  # i.e. if far foot stepped 
            initial_displacements.append(far_positions[s-1]-near_positions[s-1]) 
            final_displacements.append(far_positions[s]-near_positions[s]) 
        

  
    onebound_times = bind_times[1:]-unbind_times[1:]
    bothbound_times = unbind_times[1:]-bind_times[:-1]
    step_lengths = near_step_lens + far_step_lens
    step_times = onebound_times + bothbound_times
else:
    data_files = []
    for fname in glob.glob("data/paper_static_stepping_data*.txt"):
        data_files.append(fname)
    if len(data_files) == 0:
        print("Error, no files of form data/paper_static_stepping_data*.txt found. Exiting.")
        exit(1)

    for data_file in data_files:
        data = np.loadtxt(data_file)
        if VERBOSE: print("file found: ", data_file)
        bind_times = np.array(data[:,1])
        unbind_times = np.array(data[:,0])
        near_positions = np.around(np.array(data[:,2]), decimals=7)
        far_positions = np.around(np.array(data[:,3]), decimals=7)
        near_step_lens = (near_positions[1:] - near_positions[:-1])  # still want zero step lengths
        far_step_lens = (far_positions[1:] - far_positions[:-1])

        assert (len(near_positions)==len(far_positions))
        for s in range(1, len(near_positions)):
            assert((near_positions[s-1] == near_positions[s]) or (far_positions[s-1] == far_positions[s]))
            if(near_positions[s-1] == near_positions[s] and far_positions[s-1] == far_positions[s]): 
                continue
            if not(near_positions[s-1] == near_positions[s]):
                initial_displacements.append(near_positions[s-1]-far_positions[s-1])
                final_displacements.append(near_positions[s]-far_positions[s])
            elif not(far_positions[s-1] == far_positions[s]):
                initial_displacements.append(far_positions[s-1] - near_positions[s-1])
                final_displacements.append(far_positions[s]-near_positions[s])
 

        onebound_times = np.concatenate((onebound_times, bind_times[1:] - unbind_times[1:]))
        bothbound_times = np.concatenate((bothbound_times, unbind_times[1:] - bind_times[:-1]))
        step_lengths = np.concatenate((step_lengths, near_step_lens + far_step_lens))
    # times are already concentated in order so all we need to do is add them together
    step_times = onebound_times + bothbound_times


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
    if filename is None:
        filename = 'plots/'+graph_label.replace(' ',  '_')+".pdf"
    plt.savefig(filename)


seed_label = ''
if ALL:
    seed_label = '-multiple-seeds'

plotCounts(step_times,
           step_lengths,
           "total step time vs step length",
           'step time',
           'step length',
           xIsTimeValue=True,
           yIsTimeValue=False,
           filename='plots/time-vs-length{}.pdf'.format(seed_label))

plotCounts(initial_displacements,
           final_displacements,
           "initial disp vs final disp",
           "initial displacement",
           "final displacement",
           xIsTimeValue=False,
           yIsTimeValue=False,
           filename='plots/initial-vs-final{}.pdf'.format(seed_label))

if SHOW:
    plt.show()
if VERBOSE: print("finished")
