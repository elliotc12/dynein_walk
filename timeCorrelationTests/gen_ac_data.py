from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import autocorrelation as ac

path = sys.argv[1]
label = sys.argv[2]
    
dataTable = np.loadtxt(path, delimiter='\t', skiprows=1)
print "data table successfully loaded. Fetching energies..."

times = dataTable[:,1]
PE_1 = dataTable[:,2]
PE_2 = dataTable[:,3]
PE_3 = dataTable[:,4]
PE_4 = dataTable[:,5]
PE_5 = dataTable[:,6]

Nmax = 1000

print "Energies fetched. Generating autocorrelation function..."
rho1 = ac.autoCorrelate2(PE_1, Nmax = Nmax, skipIndex = 10, verbose=True)
rho2 = ac.autoCorrelate2(PE_2, Nmax = Nmax, skipIndex = 10, verbose=True)
rho3 = ac.autoCorrelate2(PE_3, Nmax = Nmax, skipIndex = 10, verbose=True) 
rho4 = ac.autoCorrelate2(PE_4, Nmax = Nmax, skipIndex = 10, verbose=True)
rho5 = ac.autoCorrelate2(PE_5, Nmax = Nmax, skipIndex = 10, verbose=True)

print "Functions generated. Plotting..."

plt.figure()
plt.plot(times[:Nmax], rho1, skipIndex = 10, label = 'PE_1')
plt.plot(times[:Nmax], rho2, skipIndex = 10, label = 'PE_2')
plt.plot(times[:Nmax], rho3, skipIndex = 10, label = 'PE_3')
plt.plot(times[:Nmax], rho4, skipIndex = 10, label = 'PE_4')
plt.plot(times[:Nmax], rho5, skipIndex = 10, label = 'PE_5')
plt.title("Autocorrelation functions for potential energies by domain")
plt.xlabel(r'$\Delta t$ [s]')
plt.ylabel(r'$\rho (\Delta t) $')
plt.legend(loc = 0)
plt.savefig("{}.pdf".format(label))


if sys.argv[3] == '-s': #save the data to a txt file to play with later 
    with open('{}_data.txt'.format(label), 'r') as f:
        f.write('rho1, rho2, rho3, rho4, rho5')
        for i in range(len(rho1)):
            f.write(','.join(rho1[i], rho2[i], rho3[i], rho4[i], rho5[i]))

else:
    pass

print "Finished. Saving figure..."

