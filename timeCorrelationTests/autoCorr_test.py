from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import autocorrelation as ac

path = sys.argv[1]
dataTable = np.loadtxt(path, delimiter='\t', skiprows=1)
print "data table successfully loaded. Fetching energies..."  

times = dataTable[:,1] 
PE_1 = dataTable[:,2]
PE_2 = dataTable[:,3]
PE_3 = dataTable[:,4]
PE_4 = dataTable[:,5]
PE_5 = dataTable[:,6]

print "Energies fetched. Generating autocorrelation function..."
rho1 = ac.autoCorrelate2(PE_1)


print "Function generated. Plotting..."

plt.figure()
plt.plot(times, rho1)

plt.savefig("testRun.pdf")

print "Finished. Saving figure..." 

