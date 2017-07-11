from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


# path = "thesis_movie_data.txt"
# dataTable = np.loadtxt(path, skiprows=2, comments='#', delimiter='\t')

A = np.random.rand(1000)
dt = 1 # this will be changed later to match dynein dt 

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
t = np.arange(0, dt*len(A), dt)
plt.plot(t, RHO)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\rho(\tau)$')
plt.show()
