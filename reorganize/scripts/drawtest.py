import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import dynein.draw.cartoon as draw

fig1 = plt.figure()
ax1 = fig1.add_subplot(111,aspect='equal')
draw.dyneinPolygon(0,5,10,25,30,35,'blue',1.0,ax1)
draw.dyneinPolygon(10,0,11,21,30,35,'red',1.0,ax1)
plt.xlim(-10,50)
plt.ylim(-10,50)
plt.show()
