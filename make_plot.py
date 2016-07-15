#! /usr/bin/env python2.7

import getopt
import matplotlib.pyplot as plt
import numpy as np
import sys

plot_params, data_files = getopt.getopt(sys.argv[1:], "f:x:y:ph:sxy", ["figtitle=", "xlabel=", "ylabel=", "showplot", "hline=", "scatter", "logx", "logy"])

showplot = False
hline = False
scatter = False
logx = False
logy = False

for param, value in plot_params:
    if (param == "--figtitle"):
        title_list = ([' ' if (a == '_') else a for a in value])
        title = ''.join(title_list)
    elif (param == "--xlabel"):
        xlabel = value
    elif (param == "--ylabel"):
        ylabel = value
    elif (param == "--showplot"):
        showplot = True
    elif (param == "--hline"):
        hline = True
        hlineval = value
    elif (param == "--scatter"):
        scatter = True
    elif (param == "--logx"):
        logx = True
    elif (param == "--logy"):
        logy = True

max_x_value = 0

colors = ['b', 'g', 'r', 'm']

for data_file in data_files:
    f = open(data_file, 'r')
    legend = f.readline()
    f.close()

    line = np.loadtxt(data_file, skiprows=1)
    X = line[:,0]
    Y = line[:,1]
    if (not scatter):
        plt.plot(X, Y, label=legend, color=colors.pop())
    else:
        plt.scatter(X,Y, label=legend, color=colors.pop(), s=1)

    if max(X) > max_x_value:
        max_x_value = max(X)

if hline:
    plt.plot([0, max_x_value], [hlineval, hlineval])

ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
    
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.legend(fontsize="small", loc="center left", framealpha=0.0,
            bbox_to_anchor=(1, 0.5))

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

if logx:
    ax.set_xscale('log')

if logy:
    ax.set_yscale('log')

if showplot:
    plt.show()
else:
    plt.savefig("plots/" + title.replace(" ", "_") + ".pdf", format="pdf")
