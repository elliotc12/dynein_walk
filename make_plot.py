#! /usr/bin/env python2.7

import getopt
import matplotlib.pyplot as plt
import numpy as np
import sys

plot_params, data_files = getopt.getopt(sys.argv[1:], "f:x:y:sh:", ["figtitle=", "xlabel=", "ylabel=", "showplot", "hline="])

showplot = False
hline = False

for param, value in plot_params:
    if (param == "--figtitle"):
        title = value
    elif (param == "--xlabel"):
        xlabel = value
    elif (param == "--ylabel"):
        ylabel = value
    elif (param == "--showplot"):
        showplot = True
    elif (param == "--hline"):
        hline = True
        hlineval = value

max_x_value = 0
        
for data_file in data_files:
    f = open(data_file, 'r')
    legend = f.readline()
    f.close()
    
    line = np.loadtxt(data_file, skiprows=1)
    X = line[:,0]
    Y = line[:,1]
    plt.plot(X, Y, label=legend)

    if max(X) > max_x_value:
        max_x_value = max(X)

if hline:
    plt.plot([0, max_x_value], [hlineval, hlineval])

plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.legend(fontsize="small", loc="best", framealpha=0.0, frameon=False)

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

if showplot:
    plt.show()
else:
    plt.savefig("plots/" + title.replace(" ", "_"))
