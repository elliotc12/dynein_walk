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

Nmax = None # total number of points
S_I = 1  # skip index value

print "Energies fetched. Generating autocorrelation function..."
rho1 = ac.autoCorrelateFFT(PE_1, Nmax = Nmax, skipIndex = S_I)
rho2 = ac.autoCorrelateFFT(PE_2, Nmax = Nmax, skipIndex = S_I)
rho3 = ac.autoCorrelateFFT(PE_3, Nmax = Nmax, skipIndex = S_I)
rho4 = ac.autoCorrelateFFT(PE_4, Nmax = Nmax, skipIndex = S_I)
rho5 = ac.autoCorrelateFFT(PE_5, Nmax = Nmax, skipIndex = S_I)

print "Functions generated. Saving..."

if Nmax is not None:
    times = times[:Nmax:S_I]
    saveData = np.zeros((int(Nmax/S_I), 6))
else:
    saveData = np.zeros((len(rho1), 6)) 
saveData[:,0] = times
saveData[:,1] = rho1
saveData[:,2] = rho2
saveData[:,3] = rho3
saveData[:,4] = rho4
saveData[:,5] = rho5

h = "time, rho1, rho2, rho3, rho4, rho5" 
np.savetxt("{}.txt".format(label), saveData, delimiter=',', header = h,  )

print "finished" 
        

        


