from __future__ import division 
import numpy as np

## function that creates tail domain for cartoon ## 

m = 2
N = 500
t = np.linspace(-10,10,N) 
tail = np.zeros((N,2))
tail[:,0] = -np.cos(t)+1
tail[:,1] = np.sin(t)*np.sin(0.5*t)**m

with open("tail.py", "w") as f:
    f.write("from numpy import array \n")
    f.write("array= %s" %repr(tail))

print('All done!')
