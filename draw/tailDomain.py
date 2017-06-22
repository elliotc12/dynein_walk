from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt

m = 2
N = 500
t = np.linspace(-10,10,N) 
tail = np.zeros((N,2))
tail[:,0] = -np.cos(t)+1
tail[:,1] = np.sin(t)*np.sin(0.5*t)**m

with open("tail.py", "w") as f:
    f.write("from numpy import array \n")
    f.write("array= %s" %repr(tail))
    
print 'All done!'

# def resize(x1,x2, y1, y2):
#     length = np.sqrt((x2-x1)**2+(y2-y1)**2) #desired length of tail 
#     tail[:,0] = (0.5*length)*tail[:,0]+(0.5*length) #resize at origin before rotation
#     tail[:,1] = (0.3*length)*tail[:,1]

#     theta = np.arctan2((y2-y1),(x2-x1))
#     rot = np.matrix([[np.cos(theta), -np.sin(theta)],[np.sin(theta),np.cos(theta)]])
#     for i in range(0,N):
#         tail[i]= np.dot(rot, tail[i])

#     tail[:,0] = tail[:,0] + x1
#     tail[:,1] = tail[:,1] + y1
    

# plt.plot(tail[:,0], tail[:,1])
# plt.show()


