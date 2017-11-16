#!/usr/bin/python3
from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt 
import argparse
import os 

parser = argparse.ArgumentParser(description = 'Script to generate 2 dimensional histogram from dynein stepping data')

parser.add_argument('-d', '--datafile', dest = 'data_file', action='store', type= str,default='data/paper_trajectory_movie_data.txt', help='path to data file')
parser.add_argument('-b', '--bins', dest = 'bins', action='store', type = int, default=20, help='number of bins for x and y axes')
parser.add_argument('-v', '--verbose', dest = 'verbose', action='store_true', default = False, help = 'see prints in console') 

args = parser.parse_args()

VERBOSE = args.verbose
DATAFILE = args.data_file
NUM_BINS = args.bins 


if VERBOSE: print("moving to root") 
# navigate to root directory 
os.chdir('../')
# check to see if data file exists
if not os.path.isfile(DATAFILE):
    print("Could not find data file. Please specify path using -d")
    exit(1) 
#load in data
if VERBOSE: print("Data file found- loading data...")
data = np.loadtxt(DATAFILE, delimiter='\t', skiprows=1)
if VERBOSE: print("Data loaded- graphing...")
print(data) 

