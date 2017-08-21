from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys

label = sys.argv[1]

data_10 = np.loadtxt('./data/onebound_dt1e-10.txt', delimiter=',', skiprows = 1)
data_11 = np.loadtxt('./data/onebound_dt1e-11.txt', delimiter=',', skiprows = 1)

d10 = {}
d11 = {}

d10['t'] = data_10[:,0]
d10['rho1'] = data_10[:,1]
d10['rho2'] = data_10[:,2]
d10['rho3'] = data_10[:,3]
d10['rho4'] = data_10[:,4]
d10['rho5'] = data_10[:,5]

d11['t'] = data_11[:,0]
d11['rho1'] = data_11[:,1]
d11['rho2'] = data_11[:,2]
d11['rho3'] = data_11[:,3]
d11['rho4'] = data_11[:,4]
d11['rho5'] = data_11[:,5]

plt.figure()
plt.title("dt = 1e-10 vs dt = 1e-11")

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

plt.xlabel('t [s]')
plt.ylabel(r'$\rho(\Delta t)$')
plt.legend(loc = 0)
plt.savefig(label) 
plt.show()


