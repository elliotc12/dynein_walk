import matplotlib
matplotlib.use('Agg') # to let me do this over ssh 
import subprocess, os
import argparse 
import numpy as np
import matplotlib.pyplot as plt
import dynein.run as run 


def ac(data, Nmax = None): #generate autocorrelation function
    if Nmax is not None:
        f_t = data[:Nmax]
    else:
        f_t = data[::]
    mu = np.mean(f_t)
    f_w = np.fft.fft(f_t-mu)
    f_w_conj = np.conjugate(f_w)
    norm2 = f_w*f_w_conj
    rho = np.fft.ifft(norm2)
    rho = rho/rho[0] #normalize by first element
    return rho

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'A script to test simulation performance at various values of dt')

    parser.add_argument("-L", "--label", dest = 'label', default = 'compare', help = 'label for output graphs')
    parser.add_argument("-v", "--verbose", action = 'store_true', dest = 'verbose', default = False, help = 'print extra status messages')
    parser.add_argument("-p", "--plot", action = 'store_true', dest = 'plot', default = False, help = 'plot output graphs in interactive window')  
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-o", "--onebound", action = "store_true", dest = 'onebound', default = False, help = 'check dt behavior for onebound state')
    group.add_argument("-b", "--bothbound", action = "store_true", dest = 'bothbound', default = False, help = 'check dt behavior for bothbound state')

    args = parser.parse_args()

    LABEL = args.label
    VERBOSE = args.verbose
    ONEBOUND = args.onebound
    BOTHBOUND = args.bothbound
    PLOT = args.plot

    if ONEBOUND is False and BOTHBOUND is False:    # make onebound default while keeping them mutually exclusives args 
        ONEBOUND = True 
    
    if VERBOSE:
        print "Running simulations..."

    os.chdir("../")   # change to folder with generate_stepping_data

    data_files = []   # 
    
    if ONEBOUND:
        if VERBOSE:
            print "ONEBOUND - dt 1e-10"
            print "ONEBOUND - dt 1e-11"
            print "ONEBOUND - dt 1e-11"

        l = 'oneboundAC'
        
        basename1 = run.sim(**{"k_b": 10e-10, "k_ub": 10e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-10, "label": l,
                     "constant-write": True, "runtime": 5e-4})
        basename2 = run.sim(**{"k_b": 10e-10, "k_ub": 10e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-11, "label": l,
                     "constant-write": True, "runtime": 5e-4})
        basename3 = run.sim(**{"k_b": 10e-10, "k_ub": 10e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-12, "label": l,
                     "constant-write": True, "runtime": 5e-4})

        if VERBOSE:
            print "Saving {}, {}, {} in ../data".format(basename1, basename2, basename3)

        data_files.append(basename1)
        data_files.append(basename2)
        data_files.append(basename3) 

    elif BOTHBOUND:
        if VERBOSE:
            print "BOTHBOUND - dt 1e-10"
            print "BOTHBOUND - dt 1e-11"
            print "BOTHBOUND - dt 1e-12"

        l = 'bothboundAC'

        basename1 = run.sim(**{"k_b": 10e20, "k_ub": 10e-10, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-10, "label": l,
                     "constant-write": True, "runtime": 5e-4})
        basename2 = run.sim(**{"k_b": 10e20, "k_ub": 10e-10, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-11, "label": l,
                     "constant-write": True, "runtime": 5e-4})
        basename3 = run.sim(**{"k_b": 10e20, "k_ub": 10e-10, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-12, "label": l,
                     "constant-write": True, "runtime": 5e-4})

        if VERBOSE:
            print "Saving {}, {}, {} in ../data".format(basename1, basename2, basename3)


        data_files.append(basename1)
        data_files.append(basename2)
        data_files.append(basename3)
        
    else:  
        print "Not sure if bothbound or onebound was selected. Run again using flag -o or -b"
        exit(1)

    # process data and generate autocorrelation function
    if VERBOSE:
        print "Simulations finished. Reading data"

    usefullData = {}
    for file in data_files:
        print file
        dataTable = np.loadtxt("data/stepping_movie_data_"+file+".txt", delimiter='\t', skiprows=1)
        if VERBOSE:
            print "{} successfully loaded".format(file)
        times = dataTable[:,1]
        PE_1 = dataTable[:,2]
        PE_2 = dataTable[:,3]
        PE_3 = dataTable[:,4]
        PE_4 = dataTable[:,5]
        PE_5 = dataTable[:,6]

        if VERBOSE:
            print "Generating autocorrelation function"

        Nmax = None
        rho1 = ac(PE_1, Nmax = Nmax)
        rho2 = ac(PE_2, Nmax = Nmax)
        rho3 = ac(PE_3, Nmax = Nmax)
        rho4 = ac(PE_4, Nmax = Nmax)
        rho5 = ac(PE_5, Nmax = Nmax)

        if VERBOSE:
            print "putting into data into dictionary"

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
        
    if VERBOSE:
        print "graphing...\n"

    fig1 = plt.figure()

    for key in usefullData:
        if VERBOSE:
            print key
        dt_loc = key.find("dt-1e")
        dt = key[dt_loc:dt_loc+8]

        plt.plot(usefullData[key]['times'], usefullData[key]['rho1'], label="rho1 {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho2'], label="rho2 {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho3'], label="rho3 {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho4'], label="rho4 {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['rho5'], label="rho5 {}".format(dt))
    plt.legend(loc = 0)
    plt.xlim(0, 5*10**-9)
    plt.xlabel('t [s]')
    plt.ylabel(r'$\rho(\Delta t)$')
    plt.savefig(LABEL+'_ac.pdf')

    fig2 = plt.figure()
    for key in usefullData:
        if VERBOSE:
            print key
        dt_loc = key.find("dt-1e")
        dt = key[dt_loc:dt_loc+8]
        
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_1'], label="PE_1 {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_2'], label="PE_2 {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_3'], label="PE_3 {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_4'], label="PE_4 {}".format(dt))
        plt.plot(usefullData[key]['times'], usefullData[key]['PE_5'], label="PE_5 {}".format(dt))
    plt.legend(loc=0)
    plt.xlabel('t [s]')
    plt.ylabel('U(t)')
    #plt.xlim(0,5e-9)
    plt.savefig(LABEL+'_U.pdf')

    fig3 = plt.figure()
    for key in usefullData:
        if VERBOSE:
            print key
        dt_loc = key.find("dt-1e")
        dt = key[dt_loc:dt_loc+8]
        
        plt.semilogx(usefullData[key]['times'], usefullData[key]['PE_1'], label="PE_1 {}".format(dt))
        plt.semilogx(usefullData[key]['times'], usefullData[key]['PE_2'], label="PE_2 {}".format(dt))
        plt.semilogx(usefullData[key]['times'], usefullData[key]['PE_3'], label="PE_3 {}".format(dt))
        plt.semilogx(usefullData[key]['times'], usefullData[key]['PE_4'], label="PE_4 {}".format(dt))
        plt.semilogx(usefullData[key]['times'], usefullData[key]['PE_5'], label="PE_5 {}".format(dt))
    plt.legend(loc=0)
    plt.xlabel('t [s]')
    plt.ylabel('U(t)')
    #plt.xlim(0,5e-9)

    plt.savefig(options.label+'_U_vs_logt.pdf')

    if options.p is not False: 
        plt.show() 


    plt.savefig(LABEL+'_U_vs_logt.pdf')

 
    if PLOT:
        plt.show() 

        

        
