from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys 

## autocorrelation fucntion definitions ##
def autoCorrelate1(data):
    rho = np.zeros(len(data)) # autocorrelation function to be returned
    n = len(data)
    mu = np.mean(data) #mean 
    sig2 = np.var(data) #variance

    #note: need to check range function for final index
    for k in range(0, n): #index the k's up to n-1
        R = 0
        for t in range(0, n-k): #index t's up to n-k 
            R = R + (data[t]-mu)*(data[t+k]-mu)

        rho[k] = (1/((n-k)*sig2))*R
    return rho 

def autoCorrelate2(data, Nmax = None, verbose = False):
    n = len(data)
    if Nmax is not None:
        n = Nmax
    rho = np.zeros(n)
    mu = np.mean(data)

    percent_done = int(n/100 + 1)
    for k in range(0, n):
        R = 0
        N = 0
        for t in range(0, len(data)-k):
            R = R + (data[t]-mu)*(data[t+k]-mu)
            N = N + (data[t]-mu)**2
        if verbose and k % percent_done == 0:
            print '{}% done...'.format(k/n*100)
        rho[k] = R/N

    return rho

if __name__ == "__main__":
    A = np.random.rand(1000)
    rho1 = autoCorrelate1(A)
    rho2 = autoCorrelate2(A)
    l = len(A)

    plt.figure()
    plt.plot(rho1, 'b')
    plt.xlim(-5, l)
    plt.ylim(-20,20)
    plt.xlabel('n')
    plt.ylabel('np.random.rand(n)')
    plt.title('Wikipedia R(k) method')

    plt.figure()
    plt.plot(rho2, 'k')
    plt.xlim(-5, l)
    plt.ylim(-20,20) 
    plt.xlabel('n')
    plt.ylabel('np.random.rand(n)')
    plt.title('Prof Roundy Method')

    plt.show()
