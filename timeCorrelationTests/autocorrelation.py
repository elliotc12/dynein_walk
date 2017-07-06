from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


# path = "thesis_movie_data.txt"
# dataTable = np.loadtxt(path, skiprows=2, comments='#', delimiter='\t')

A = np.array([1,1,1,1,1,3,1,4,2,1,2,5,6,2,1,1,1,1,1])
dt = 1

def autoCorrelate(data):
    rho = np.zeros(len(data)) # autocorrelation function to be returned
    n = len(data)
    mu = np.mean(data) #mean 
    sig2 = np.var(data) #variance

    #note: need to check range function for final index
    for k in range(0, n): #index the k's up to n-1
        R = 0
        for t in range(0, n-k): #index t's up to n-k 
            R = R + (data[t]-mu)*(data[t+k]-mu)

        rho[k] = (1/((n-k)*sig2**2))*R
    return rho 
                   
RHO = autoCorrelate(A)
t = np.arrange(0, dt*len(A), dt)
plt.plot(t, RHO)
plt.show()
