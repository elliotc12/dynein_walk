#! /usr/bin/env python2.7

import getopt
import matplotlib.pyplot as plt
import numpy as np
import sys

plot_params, data_files = getopt.getopt(sys.argv[1:], "f:x:y:s", ["figtitle=", "xlabel=", "ylabel=", "showplot"])

showplot = False

for param, value in plot_params:
    if (param == "--figtitle"):
        title = value
    elif (param == "--xlabel"):
        xlabel = value
    elif (param == "--ylabel"):
        ylabel = value
    elif (param == "--showplot"):
        showplot = True

for data_file in data_files:
    f = open(data_file, 'r')
    legend = f.readline()
    f.close()
    
    line = np.loadtxt(data_file, skiprows=1)
    X = line[:,0]
    Y = line[:,1]
    plt.plot(X, Y, label=legend)

plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.legend()

if showplot:
    plt.show()
else:
    plt.savefig("plots/" + title.replace(" ", "_"))
