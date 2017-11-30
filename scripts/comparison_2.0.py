#!/usr/bin/python3
import matplotlib
matplotlib.use('Agg') # to let me do this over ssh 
import subprocess, os
import argparse 
import numpy as np
import matplotlib.pyplot as plt
import dynein.run as run 


#generate autocorrelation function 
def ac(data, Nmax = None):
    if Nmax is not None:
        f_t = data[:Nmax]
    else:
        f_t = data[::]
    mu = np.mean(f_t)
    f_w = np.fft.fft(f_t-mu)
    f_w_conj = np.conjugate(f_w)
    norm2 = f_w*f_w_conj
    ac = np.fft.ifft(norm2)
    ac = ac/ac[0] #normalize by first element
    return ac


if __name__ == '__main__':

    #set up argparse to handle command line flags 
    parser = argparse.ArgumentParser(description = 'A script to test simulation performance at various values of dt')

    #add flags
    parser.add_argument("-L", "--label", dest = 'label', default = 'compare', help = 'label for output data and graphs')
    parser.add_argument('-v', '--verbose', action='store_true', dest = 'verbose', default = False, help = 'print output to console')
    parser.add_argument('-p', '--plot', action = 'store_true', dest = 'plot', default = False, help = 'print output to console')
    parser.add_argument('-s', '--seed', action = 'store', type = int, default = 1, help = 'Set seed for random number generator')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-o', '--onebound', action='store_true', dest = 'onebound', default = False, help = 'check dt behavior for onebound state')
    group.add_argument('-b', '--bothbound', action='store_true', dest = 'bothbound', default = False, help = 'check dt behavior for bothbound state')

    args = parser.parse_args()

    #set variables so that we don't have to keep calling argparse
    LABEL = args.label
    VERBOSE = args.verbose
    ONEBOUND = args.onebound
    BOTHBOUND = args.bothbound
    PLOT = args.plot
    SEED = args.seed

    #make onebound the default while keeping them mutually exclusive
    if ONEBOUND is False and BOTHBOUND is False:
        ONEBOUND = True

    
    #navigate to root and ensure needed directories exist
    if os.path.exists('comparison_2.0.py'):
        print("navigating to root directory") 
        os.chdir("../")
    if VERBOSE:
        print(os.getcwd(), '\nChecking if data data dir exists:', os.path.exists('data'))
        print('Checking if figs dir exists:', os.path.exists('figs'))
    if not os.path.exists('data'): os.makedirs('data')
    if not os.path.exists('figs'): os.makedirs('figs')
    
    #dts = ['5e-10']
    dts = ['5e-10', '1e-10', '1e-11', '1e-12']
   
    if VERBOSE:
        print('dt values:') 
        for dt in dts: print('\t', dt)

    #set binding constants if onebound or bothbound 
    if ONEBOUND:
        k_b = 1e-10
        k_ub = 1e20
        LABEL = 'oneboundAC' 
    if BOTHBOUND:
        k_b = 1e20
        k_ub = 1e-10
        LABEL = 'bothboundAC'
        
    #set lengths and angles 
    cb = 0.1
    cm = 0.5
    ct = 0.2
    ls = 10.49
    lt = 23.8  # from urnavicius 2015
    eqb = 120  # from redwine 2012 supplemental
    eqmpre = 200# from burgess 2002, 360-160
    eqmpost = 224  # from burgess 2002, 360-136
    eqt = 0

    runtime = 5e-4 

    data_files = [] 
    
    for dt in dts:
        print('\nsaving {} in data/'.format(dt)) 
        #set framerate based on dt 
        if float(dt)>1e-10:
            framerate = dt
        else:
            framerate = 1e-10

        #run simulation 
        basename = run.sim(**{ 'k_b':k_b,
                         'k_ub': k_ub, 
                         'cb': cb,
                         'cm':cm,
                         'ct': ct,
                         'ls': ls,
                         'lt': lt,
                         'eqb': eqb,
                         'eqmpre':eqmpre,
                         'eqmpost':eqmpost,
                         'eqt':eqt,
                         'dt': float(dt),
                         'label': LABEL,
                         'constant-write': True, 
                         'seed': SEED,
                         'runtime': runtime,
                         'framerate': framerate,
                         'crash-movie': False,
                         'no-slurm': True})
        data_files.append(basename) 

    if VERBOSE: print('\nSimulations finished. Reading data')
    
    #set up data structure
    if VERBOSE:
        for dt in dts: print('dt:', dt)
        for file in data_files: print('file name:', file) 
    
    usefullData = {}
    for file in data_files:
        if VERBOSE: print('loading:', file)

        dataTable = np.loadtxt('data/stepping_movie_data_'+file+'.txt', delimiter='\t', skiprows=1)
        if VERBOSE: print('{} has been successfully loaded'.format(file))
        times = dataTable[:,1]
        PE_1 = dataTable[:,2]
        PE_2 = dataTable[:,3]
        PE_3 = dataTable[:,4]
        PE_4 = dataTable[:,5]
        PE_5 = dataTable[:,6]

        if VERBOSE: print('Generating autocorrelation functions')

        Nmax = None
        rho1 = ac(PE_1, Nmax=Nmax)
        rho2 = ac(PE_2, Nmax=Nmax)
        rho3 = ac(PE_3, Nmax=Nmax)
        rho4 = ac(PE_4, Nmax=Nmax)
        rho5 = ac(PE_5, Nmax=Nmax) 

        if VERBOSE: print("putting data into dictionary")

       # make dictionary of all important quantities
        dt_dict = {}
        dt_dict['times'] = times
        dt_dict['PE_1'] = PE_1
        dt_dict['PE_2'] = PE_2
        dt_dict['PE_3'] = PE_3
        dt_dict['PE_4'] = PE_4
        dt_dict['PE_5'] = PE_5
        dt_dict['rho1'] = rho1
        dt_dict['rho2'] = rho2
        dt_dict['rho3'] = rho3
        dt_dict['rho4'] = rho4
        dt_dict['rho5'] = rho5

        # now add this dictionary to usefullData for each file 
        
        usefullData[file] = dt_dict
    if VERBOSE: print("\ndata generated... graphing")

    #generate ac graph # 
    fig1 = plt.figure()
    for key in usefullData:
        if VERBOSE: print(key, '\n',type(key))
        dt_loc = key.find("dt-")
        print('dt_loc',dt_loc) 
        dt = key[dt_loc+3:dt_loc+8]
        print('dt:',dt) 
        plt.plot(usefullData[key]['times'], usefullData[key]['rho1'], label="rho1 dt {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho2'], label="rho2 dt {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho3'], label="rho3 dt {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho4'], label="rho4 dt {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho5'], label="rho5 dt {}".format(dt))
    plt.legend(loc = 0)
    plt.xlim(0, 10**-8)
    plt.ylim(-0.5,1.25) 
    plt.xlabel('t [s]')
    plt.ylabel(r'$\rho(\Delta t)$')
    plt.savefig('figs/'+LABEL+'_ac.pdf')    








        
