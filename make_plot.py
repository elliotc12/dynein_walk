#! /usr/bin/env python2.7

import getopt
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

plot_params, data_files = getopt.getopt(sys.argv[1:], "f:x:y:ph:sxyk:", ["figtitle=", "xlabel=", "ylabel=", "showplot", "hline=", "scatter", "logx", "logy", "skiprows="])

showplot = False
hline = False
scatter = False
logx = False
logy = False

skiprows = 1

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
    elif (param == "--skiprows"):
        skiprows = int(value)

ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.7, box.height])

max_x_value = 0
colors = ['b', 'g', 'r', 'm']

for data_file in data_files:
    if (data_file.find("config") == -1):
        f = open(data_file, 'r')
        legend = f.readline()
        f.close()

        line = np.loadtxt(data_file, skiprows=1)
        X = line[:,0]
        Y = line[:,1]
        X = X[::skiprows]
        Y = Y[::skiprows]
        if (not scatter):
            plt.plot(X, Y, label=legend, color=colors.pop())
        else:
            plt.scatter(X,Y, label=legend, color=colors.pop(), s=1)

        if max(X) > max_x_value:
            max_x_value = max(X)
    else:
        f = open(data_file, 'r')
        config_txt = f.read()
        box = TextArea(config_txt, textprops=dict(size=10))
        anchored_box = AnchoredOffsetbox(loc=2, child=box, bbox_to_anchor=(1,1),
                                         bbox_transform=ax.transAxes, frameon=False)

if hline:
    plt.plot([0, max_x_value], [hlineval, hlineval])

plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.legend(fontsize="small", loc="lower left", framealpha=0.0, bbox_to_anchor=(1, 0))
ax.add_artist(anchored_box)

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

if logx:
    ax.set_xscale('log')

if logy:
    ax.set_yscale('log')

if showplot:
    plt.show()
else:
    plt.savefig("plots/" + title.replace(" ", "_") + ".pdf", format="pdf")
