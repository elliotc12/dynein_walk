#!usr/bin/env python

import subprocess, os
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt 

def simulate(n): # run simulations
    os.chdir("../run_scripts") # can't change dir with subprocess.call() for some reason...
    print os.getcwd() 
    if n == 1:
        subprocess.call('python onebound_simulations.py', shell=True)
    elif n ==2: 
        subprocess.call('python bothbound_simulations.py', shell=True)
    else:
        print "incorrect usage. n = 1 or 2"
        
def ac(data, Nmax = None): # generate autocorrelation function 
    if Nmax is not None:
        f_t = data[:Nmax]
    else:
        f_t = data[::]
    mu = np.mean(f_t)
    f_w = np.fft.fft(f_t-mu)
    f_w_conj = np.conjugate(f_w)
    norm2 = f_w*f_w_conj
    rho = np.fft.ifft(norm2)
    rho = rho/rho[0]   # normalize by first element 
    return rho 
    

def main():
    # set up command line flags.  
    parser = OptionParser()
    parser.add_option("-l", "--label", dest = 'label',
                      help = 'label for output graphs')
    parser.add_option("-v", "--verbose", action="store_true", dest = 'verbose', default = False, 
                      help = 'print extra status messages')
    parser.add_option("-n", dest ="n", default = 1,
                      help = 'Run onebound or bothbound simulation. 1 = onebound, 2 = bothbound. 1 is the default.') 
    (options, args) = parser.parse_args()


    if options.verbose:
        print "n = {}".format(options.n)
        print "running simulation to generate data" 
    simulate(int(options.n))
    os.chdir("../data") # go back to original directory

    
    if options.verbose:
        print "generating autocorrelation funcitons"

    searchString = os.popen("find stepping_movie_data*").read() 
    rawDataFiles = searchString.split('\n')
    print "Data files:\n"
    print rawDataFiles
    print "\n"
    
    count = 10
    
    usefullData = {} 
    for file in rawDataFiles:
        print file
        if file is not '':
            dataTable = np.loadtxt(file, delimiter='\t', skiprows=1)
            if options.verbose:
                print "data table successfully loaded. Reading data"        
            times = dataTable[:,1]
            PE_1 = dataTable[:,2]
            PE_2 = dataTable[:,3]
            PE_3 = dataTable[:,4]
            PE_4 = dataTable[:,5]
            PE_5 = dataTable[:,6]
            if options.verbose:
                print "Generating autocorrelation function" 
            Nmax = None
            rho1 = ac(PE_1, Nmax=Nmax)
            rho2 = ac(PE_2, Nmax=Nmax)
            rho3 = ac(PE_3, Nmax=Nmax)
            rho4 = ac(PE_4, Nmax=Nmax)
            rho5 = ac(PE_5, Nmax=Nmax)

            if options.verbose:
                print "putting into new data structure"
            #make a dictionary of all the important values
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
            # add this dictionary to usefullData

            usefullData[str(count)+"_"+str(options.n)] = dt_dict
        else:
            pass #for some reason there is a '' empty element in the split line.
                 #this deals with it.
        count += 1 
    if options.verbose:
        print "graphing...\n"

    fig1 = plt.figure()
    
    for key in usefullData:
        print key 
        plt.plot(usefullData[key]['times'], usefullData[key]['rho1'], label="rho1 {}".format(key))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho2'], label="rho2 {}".format(key))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho3'], label="rho3 {}".format(key))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho4'], label="rho4 {}".format(key))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho5'], label="rho5 {}".format(key))
    plt.legend(loc = 0)
   # plt.xlim(0, 5*10**-9)
    plt.xlabel('t [s]')
    plt.ylabel(r'$\rho(\Delta t)$')
    plt.savefig(options.label+'_ac', format='pdf')

    fig2 = plt.figure()
    for key in usefullData:
        print key
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_1'], label="PE_1 {}".format(key))
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_2'], label="PE_2 {}".format(key))
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_3'], label="PE_3 {}".format(key))
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_4'], label="PE_4 {}".format(key))
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_5'], label="PE_5 {}".format(key))
    plt.legend(loc=0)
    plt.xlabel('t [s]')
    plt.ylabel('U(t) ')
    plt.savefig(options.label+'_U',format='pdf') 
    plt.show() 
        
if __name__ == '__main__':
    main() 
    
