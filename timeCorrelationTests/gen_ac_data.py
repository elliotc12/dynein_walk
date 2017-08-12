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

Nmax = 6000 # total number of points
S_I = 10  # skip index value
dt = times[1]-times[0]
print

print "Energies fetched. Generating autocorrelation function..."
rho1 = ac.autoCorrelate2(PE_1, Nmax = Nmax, skipIndex = S_I, verbose=True)
rho2 = ac.autoCorrelate2(PE_2, Nmax = Nmax, skipIndex = S_I, verbose=True)
rho3 = ac.autoCorrelate2(PE_3, Nmax = Nmax, skipIndex = S_I, verbose=True) 
rho4 = ac.autoCorrelate2(PE_4, Nmax = Nmax, skipIndex = S_I, verbose=True)
rho5 = ac.autoCorrelate2(PE_5, Nmax = Nmax, skipIndex = S_I, verbose=True)

print "Functions generated. Saving..."

# write the data to a separate file so we can graph it 
with open("ac_data_{}.txt".format(label), 'w') as f:
    f.write("time, rho1, rho2, rho3, rho4, rho5")
    if S_I is not 1:
            times = np.arange(0, Nmax*dt, dt)
    for i in range(0, len(rho1)):
        f.write(','.join([str(times[i]), str(rho1[i]), str(rho2[i]), str(rho3[i]), str(rho4[i]), str(rho5[i])]))
        f.write('\n') 
        
            

# plt.figure()
# plt.plot(times[:Nmax], rho1, label = 'PE_1')
# plt.plot(times[:Nmax], rho2, label = 'PE_2')
# plt.plot(times[:Nmax], rho3, label = 'PE_3')
# plt.plot(times[:Nmax], rho4, label = 'PE_4')
# plt.plot(times[:Nmax], rho5, label = 'PE_5')
# plt.title("Autocorrelation functions for potential energies by domain")
# plt.xlabel(r'$\Delta t$ [s]')
# plt.ylabel(r'$\rho (\Delta t) $')
# plt.legend(loc = 0)
# plt.savefig("{}.pdf".format(label))


# if sys.argv[3] == '-s': #save the data to a txt file to play with later 
#     with open('{}_data.txt'.format(label), 'r') as f:
#         f.write('rho1, rho2, rho3, rho4, rho5')
#         for i in range(len(rho1)):
#             f.write(','.join(rho1[i], rho2[i], rho3[i], rho4[i], rho5[i]))

# else:
#     pass

# print "Finished. Saving figure..."

