#!/usr/bin/python3
from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt 
import argparse
import os 

parser = argparse.ArgumentParser(description = 'Script to generate 2 dimensional histogram from dynein stepping data')

parser.add_argument('-d', '--datafile', dest = 'data_file', action='store', type= str,default='data/paper_histogram_stepping_data-1.txt', help='path to data file')
parser.add_argument('-b', '--bins', dest = 'bins', action='store', type = int, default=20, help='number of bins for x and y axes')
parser.add_argument('-v', '--verbose', dest = 'verbose', action='store_true', default = False, help = 'see prints in console')
parser.add_argument('-s', '--show', dest = 'show', action='store_true', default = False, help = 'show graphs in matplotib windows') 

args = parser.parse_args()

VERBOSE = args.verbose
DATAFILE = args.data_file
NUM_BINS = args.bins
SHOW = args.show


if VERBOSE: print("moving to root") 
# navigate to root directory 
os.chdir('../')
# check to see if data file exists
if not os.path.isfile(DATAFILE):
    print("Could not find data file. Please specify path using -d")
    exit(1) 
#load in data
if VERBOSE: print("Data file found- loading data...")
data = np.loadtxt(DATAFILE, skiprows=1, comments='#', dtype=np.float64)
if VERBOSE: print("Data loaded- formatting...")


t_ub = data[:,0]    #time unbound
t_b = data[:,1] #time bothbound
nx = data[:,2]
fx = data[:,3]
L = nx-fx
L_abs = np.abs(L)

def getBins(x):
    xmin = np.min(x)
    xmax = np.max(x)
    return np.linspace(xmin, xmax, NUM_BINS)

def getBinIndex(p, bins):
    index = 0
    for i in range(0, len(bins)):
        if p < bins[i]: #need to think about logic for element equal to a bin value 
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
    plt.pcolor(x_bins,y_bins,counts)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(graph_label)
    plt.colorbar() 
    plt.savefig('figs/'+graph_label+".pdf")


plotCounts(t_ub, L, "L vs time to unbind", "t_ub", "L")
plotCounts(t_b, L, "L vs time to bind", "t_b", "L")
if SHOW: plt.show() 
if VERBOSE: print("finished")

for i in range(0, len(t_ub)):
    print(t_ub[i], t_b[i]) 
