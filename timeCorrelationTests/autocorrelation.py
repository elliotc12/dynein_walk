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


# use this method -- autoCorrelate1 was used to verify the math was correct
def autoCorrelate2(data, Nmax = None, skipIndex = 0, verbose = False):
    '''
    verbose prints the percent done while running
    skipIndex skips a fixed number of indices during calculation--- 
    this should be used to insure simulations of different dt have 
    the same number of points in calculation 
    ''' 
    n = len(data)
    if Nmax is not None:
        n = Nmax     # total number of points desired 
    rho = np.zeros(n)
    mu = np.mean(data)

    percent_done = int(n/100 + 1)


    if skipIndex is not 0:
        for k in range(0, n):
            R = 0
            N = 0
            for t in range(0, len(data)-k, skipIndex):
                R = R + (data[t]-mu)*(data[t+k]-mu)
                N = N + (data[t]-mu)**2
            if verbose and k % percent_done == 0:
                print '{}% done...'.format(k/n*100)
            rho[k] = R/N
    
    else:
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
    rho10 = autoCorrelate2(A, Nmax = 500, skipIndex = 10)
    rho0 = autoCorrelate2(A, Nmax = 500)
    rho5 = autoCorrelate2(A, Nmax = 500, skipIndex = 5) 
    l = len(A)

    plt.figure()
    plt.plot(rho10, 'k')
    plt.plot(rho0, 'b')
    plt.plot(rho5, 'r') 
    plt.xlabel('n')
    plt.ylabel('np.random.rand(n)')
    plt.title('Test graph')
    plt.xlabel(r'$\Delta t$ [s]')
    plt.ylabel(r'$\rho (\Delta t) $')

    plt.show()
