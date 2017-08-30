from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys

label = sys.argv[1]

data_10 = np.loadtxt('./data/bothbound_dt1e-10.txt', delimiter=',', skiprows = 1)
data_11 = np.loadtxt('./data/bothbound_dt1e-11.txt', delimiter=',', skiprows = 1)
data_12 = np.loadtxt('./data/bothbound_dt1e-12.txt', delimiter=',', skiprows = 1) 

d10 = {}
d11 = {}
d12 = {}

# pull data from saved file and normalize by first value since it should be
# the largest

d10['t'] = data_10[:,0]
d10['rho1'] = data_10[:,1]/data_10[0,1]
d10['rho2'] = data_10[:,2]/data_10[0,2]
d10['rho3'] = data_10[:,3]/data_10[0,3]
d10['rho4'] = data_10[:,4]/data_10[0,4]
d10['rho5'] = data_10[:,5]/data_10[0,5]

d11['t'] = data_11[:,0]
d11['rho1'] = data_11[:,1]/data_11[0,1]
d11['rho2'] = data_11[:,2]/data_11[0,2]
d11['rho3'] = data_11[:,3]/data_11[0,3]
d11['rho4'] = data_11[:,4]/data_11[0,4]
d11['rho5'] = data_11[:,5]/data_11[0,5]

d12['t'] = data_12[:,0]
d12['rho1'] = data_12[:,1]/data_12[0,1]
d12['rho2'] = data_12[:,2]/data_12[0,2]
d12['rho3'] = data_12[:,3]/data_12[0,3]
d12['rho4'] = data_12[:,4]/data_12[0,4]
d12['rho5'] = data_12[:,5]/data_12[0,5]

# plot the data 

plt.figure()

plt.plot(d10['t'], d10['rho1'], 'r', label='rho1 dt1e-10')
plt.plot(d10['t'], d10['rho2'], 'g',label='rho2 dt1e-10')
plt.plot(d10['t'], d10['rho3'], 'b',label='rho3 dt1e-10')
plt.plot(d10['t'], d10['rho4'], 'k',label='rho4 dt1e-10')
plt.plot(d10['t'], d10['rho5'], 'm', label='rho5 dt1e-10')

plt.plot(d11['t'], d11['rho1'], 'r--', label = 'rho1 dt1e-11')
plt.plot(d11['t'], d11['rho2'], 'g--', label = 'rho2 dt1e-11')
plt.plot(d11['t'], d11['rho3'], 'b--', label = 'rho3 dt1e-11')
plt.plot(d11['t'], d11['rho4'], 'k--', label = 'rho4 dt1e-11')
plt.plot(d11['t'], d11['rho5'], 'm--', label = 'rho5 dt1e-11')

plt.plot(d12['t'], d12['rho1'], 'r-.', label = 'rho1 dt1e-12')
plt.plot(d12['t'], d12['rho2'], 'g-.', label = 'rho2 dt1e-12')
plt.plot(d12['t'], d12['rho3'], 'b-.', label = 'rho3 dt1e-12')
plt.plot(d12['t'], d12['rho4'], 'k-.', label = 'rho4 dt1e-12')
plt.plot(d12['t'], d12['rho5'], 'm-.', label = 'rho5 dt1e-12')

plt.title("dt comparison - bothbound")
plt.xlabel('t [s]')
plt.ylabel(r'$\rho(\Delta t)$')
plt.legend(loc = 0)
plt.savefig(label) 
plt.show() 


