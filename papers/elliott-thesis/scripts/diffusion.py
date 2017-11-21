#!/usr/bin/env python2.7


import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

x = np.arange(0,100,1)

plt.plot(x, x*x/(2*.2), 'b', label="Diffusion", linewidth=2.0)
plt.plot(x, x/.1, 'g', label="Dynein", linewidth=2.0)

plt.legend(loc='best')

plt.title("Diffusion vs motor transport of a 140kD protein")
plt.xlabel("log(r) $\mu m$")
plt.ylabel("time (s)")

ax = plt.gca()

ax.set_xscale('log')

ax.set_ylim([0,2000])

plt.show()
