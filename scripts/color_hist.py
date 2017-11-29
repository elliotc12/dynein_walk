#!/usr/bin/python3
from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt 
import argparse
import os 

parser = argparse.ArgumentParser(description = 'Script to generate 2 dimensional histogram from dynein stepping data')

parser.add_argument('-d', '--datafile', dest = 'data_file', action='store', type= str,
                    default='data/paper_histogram_stepping_data-1.txt', help='path to data file')
parser.add_argument('-b', '--bins', dest = 'bins', action='store', type = int, default=20,
                    help='number of bins for x and y axes')
parser.add_argument('-v', '--verbose', dest = 'verbose', action='store_true', default = False,
                    help = 'see prints in console')
parser.add_argument('-s', '--show', dest = 'show', action='store_true', default = False,
                    help = 'show graphs in matplotib windows') 
parser.add_argument('-a', '--all', dest = 'All', action='store_true', default = False,
                    help = 'generate plots for all paper_histogram_stepping_data files')
parser.add_argument('-c', '--colormap', dest = 'cmap', action='store', type=str,
                    default='plasma', help='set color map for plots') 


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

if not ALL:
    if not os.path.isfile(DATAFILE):
        print("Could not find data file. Please specify path using -d")
        exit(1)
        #load in data
    if VERBOSE: print("Data file found- loading data...")
    data = np.loadtxt(DATAFILE)
    if VERBOSE: print("Data loaded- formatting...")
    bind_times = np.array(data[:,1])
    unbind_times = np.array(data[:,0])
    near_positions = np.around(np.array(data[:,2]), decimals=7)
    far_positions = np.around(np.array(data[:,3]), decimals=7)
    near_step_lens = near_positions[1:] - near_positions[:-1]
    far_step_lens = far_positions[1:] - far_positions[:-1]

    onebound_times = np.concatenate((onebound_times, bind_times[1:] - unbind_times[1:]))
    bothbound_times = np.concatenate((bothbound_times, unbind_times[1:] - bind_times[:-1]))
    step_lengths = near_step_lens + far_step_lens
    step_times = onebound_times + bothbound_times 
else:
    data_files = []
    for fname in os.listdir("data/"):
        if os.path.isfile("data/" + fname):
            if ("data/paper_histogram_stepping_data" in "data/" + fname):
                data_files.append("data/" + fname)

                if len(data_files) == 0:
                    print("Error, no files of form data/paper_histogram_stepping_data*.txt found. Exiting.")
                    exit(1)


    for data_file in data_files:
        data = np.loadtxt(data_file)
        if len(data) < 3 or str(type(data[0])) == "<type 'numpy.float64'>":
            continue

        bind_times = np.array(data[:,1])
        unbind_times = np.array(data[:,0])
        near_positions = np.around(np.array(data[:,2]), decimals=7)
        far_positions = np.around(np.array(data[:,3]), decimals=7)
        near_step_idxs = near_positions[1:] != near_positions[:-1]
        far_step_idxs = far_positions[1:] != far_positions[:-1]
        near_step_lens = (near_positions[1:] - near_positions[:-1])[near_step_idxs]
        far_step_lens = (far_positions[1:] - far_positions[:-1])[far_step_idxs]

        onebound_times = np.concatenate((onebound_times, bind_times[1:] - unbind_times[1:]))
        bothbound_times = np.concatenate((bothbound_times, unbind_times[1:] - bind_times[:-1]))
        step_lengths = np.concatenate((step_lengths, near_step_lens, far_step_lens))

    step_times = onebound_times + bothbound_times   
##....................## 


def getBinIndex(p, bins):
    for i in range(0, len(bins)-1):
        if p >= bins[i] and p <= bins[i+1]: # deal with case for p=bins[0], p=bins[-1] 
            return i
    assert(False) # throw exception if we can't find a bin for a value 

def getCounts(X,Y):
    if VERBOSE: print(X.shape, Y.shape) 
    assert X.shape == Y.shape  
        
    if VERBOSE: print("binning data") 
    xbins = np.linspace(0, X.max(), NUM_BINS+1) # insure that times start at zero and not the min time
    ybins = np.linspace(Y.min(), Y.max(), NUM_BINS+1)
    
    if VERBOSE: print("counting") 
    counts = np.zeros((len(xbins),len(ybins)))
    
    for k in range(0, len(X)):
        i = getBinIndex(X[k], xbins)
        j = getBinIndex(Y[k], ybins)
        counts[j,i] += 1 # rows are y values, columns are x values
    return xbins, ybins, counts

def plotCounts(x,y, graph_label, x_label, y_label):
    x_bins, y_bins, counts = getCounts(x,y)
    print('counts', np.sum(counts), len(x))
    
    if VERBOSE: print("graphing")
    plt.figure()
    print(counts.shape)
    print(x_bins.shape, y_bins.shape)
    plt.pcolor(x_bins, y_bins, counts, cmap=CMAP)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    print(x_bins)
    plt.title(graph_label)
    cb = plt.colorbar()
    cb.set_label('counts') 
    plt.savefig('plots/'+graph_label.replace(' ','_')+".pdf")


 

seed_label = ''
if ALL:
    seed_label = '(multiple seeds)' 
plotCounts(step_times, step_lengths, "total step time vs step length {}".format(seed_label),
           'step time', 'step length')
plt.figure()
plt.hist2d(step_times,step_lengths, NUM_BINS, cmap=CMAP)
cb = plt.colorbar()
cb.set_label('counts')
plt.title("what it should look like")


if SHOW: plt.show() 
if VERBOSE: print("finished")

