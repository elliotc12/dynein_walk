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


args = parser.parse_args()

VERBOSE = args.verbose
DATAFILE = args.data_file
NUM_BINS = args.bins
SHOW = args.show
ALL = args.All 

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
    near_step_idxs = near_positions[1:] != near_positions[:-1]
    far_step_idxs = far_positions[1:] != far_positions[:-1]
    near_step_lens = (near_positions[1:] - near_positions[:-1])[near_step_idxs]
    far_step_lens = (far_positions[1:] - far_positions[:-1])[far_step_idxs]

    onebound_times = np.concatenate((onebound_times, bind_times[1:] - unbind_times[1:]))
    bothbound_times = np.concatenate((bothbound_times, unbind_times[1:] - bind_times[:-1]))
    step_lengths = np.concatenate((step_lengths, near_step_lens, far_step_lens))
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




def getBins(x):
    xmin = np.min(x)
    xmax = np.max(x)
    return np.linspace(xmin, xmax, NUM_BINS)

def getBinIndex(p, bins):
    index = 0
    for i in range(0, len(bins)):
        if p <= bins[i]: 
            index = i
            break 
    return index 

def getCounts(X,Y):
    if (len(X)!=len(Y)):
        print('axis do not have same length. Make sure you are using stepping_data not stepping_movie_data')
        exit(1)
        
    if VERBOSE: print("binning data") 
    xbins = getBins(X)
    ybins = getBins(Y)
    
    if VERBOSE: print("counting") 
    counts = np.zeros((len(xbins),len(ybins)))
    
    for k in range(0, len(X)):
        i = getBinIndex(X[k], xbins)
        j = getBinIndex(Y[k], ybins)
        counts[i,j] += 1
    return xbins, ybins, counts

def plotCounts(x,y, graph_label, x_label, y_label):
    x_bins, y_bins, counts = getCounts(x,y)
    
    if VERBOSE: print("graphing")
    plt.figure()
    print(counts.shape)
    print(len(x_bins), len(y_bins))
    plt.pcolor(x_bins, y_bins, counts) 
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(graph_label)
    cb = plt.colorbar()
    cb.set_label('counts') 
    plt.savefig('plots/'+graph_label.replace(' ','_')+".pdf")



# x = np.linspace(0,2*np.pi, 1000)
# y = np.sin(x)

# plotCounts(x,y,'test graph', 'x', 'y')
# plt.show() 

seed_label = ''
if ALL:
    seed_label = '(multiple seeds)' 
plotCounts(onebound_times, step_lengths, "OB times vs Step Lengths {}".format(seed_label), "t_ob", "step length")
plotCounts(bothbound_times, step_lengths, "BB times vs Step Lengths {}".format(seed_label), "t_bb", "step length")
plotCounts(step_times, step_lengths, "total step time vs step length {}".format(seed_label), 'step time', 'step length') 
plt.figure()
plt.hist2d(step_times,step_lengths, NUM_BINS)
plt.title("what it should look like") 
if SHOW: plt.show() 
if VERBOSE: print("finished")

