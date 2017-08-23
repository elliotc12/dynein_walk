import numpy as np
import matplotlib.pyplot as plt

a = np.linspace(0, 2*np.pi, 1000)
b = np.sin(a)

plt.figure()
plt.plot(a,b)
plt.savefig("testfig.png") 
