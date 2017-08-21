from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys 


def autoCorrelate1(data, Nmax = None, skipIndex = 0, verbose = False):
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

# use this method -- autoCorrelate1 was used to verify the math was correct
def autoCorrelate2(data, Nmax = None, skipIndex = 1, verbose = False):
    '''
    verbose prints the percent done while running
    skipIndex skips a fixed number of indices during calculation---
    this should be used to insure simulations of different dt have
    the same number of points in calculation
    '''
    n = len(data)
    if Nmax is not None:
        n = Nmax     # maximum index 
    rho = np.zeros(n)
    mu = np.mean(data)

    percent_done = int(n/100 + 1)
    data_minus_mu = data - mu
    
    for k in range(0, n):
        R = sum(data_minus_mu[0:len(data)-k:skipIndex]*data_minus_mu[k:len(data):skipIndex])
        N = sum((data_minus_mu**2)[0:len(data)-k:skipIndex])
        if verbose and k % percent_done == 0:
            print '{}% done...'.format(k/n*100)
        rho[k] = R/N
    return rho

def autoCorrelateFFT(data, Nmax = None, skipIndex = 1, verbose = False):
    ## need to think about units for Tmax ##
    if Nmax is not None:
        f_t = data[:Nmax:skipIndex]
    else:
        f_t = data[::skipIndex]
    mu = np.mean(f_t)
    f_w = np.fft.fft(f_t-mu)
    f_w_conj = np.conjugate(f_w)
    norm2 = f_w*f_w_conj
    rho = np.fft.ifft(norm2)
    return rho

if __name__ == "__main__":
    x = np.random.rand(10000)
    A = np.exp(x) 
    C = autoCorrelateFFT(A, Nmax = 5000, skipIndex=1)
    rho = autoCorrelate2(A, Nmax = 5000, skipIndex=10, verbose=False)
    plt.figure()
    plt.plot(C, label = 'fft')
    plt.plot(rho, '--', label = 'old')
    plt.legend(loc=0)
    plt.show()
