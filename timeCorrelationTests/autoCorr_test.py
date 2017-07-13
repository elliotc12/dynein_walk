from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import autocorrelation as ac

path = sys.argv[1]
dataTable = np.loadtxt(path, delimiter='\t', skiprows=1)
print "Successfully loaded!"
#print np.shape(dataTable) --- 

times = dataTable[:,1] 
PE_1 = dataTable[:,2]
PE_2 = dataTable[:,3]
PE_3 = dataTable[:,4]
PE_4 = dataTable[:,5]
PE_5 = dataTable[:,6]

print "Data loaded... running autocorrelation."


rho1 = ac.autoCorrelate2(PE_1)


plt.figure()
plt.plot(times, rho1)

plt.show()
