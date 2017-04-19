#!/usr/bin/python2.7
from __future__ import division, print_function

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math
import os
import sys

if len(sys.argv) != 2:
    print("Usage: ", sys.argv[0], " /path/to/stepping_data/dir/")
    exit(1)

datafiles = [s for s in os.listdir(sys.argv[1]) if "stepping_data" in s]

t_step = []
t_ob = []
t_bb = []
t_proc = []
kbs = []
kubs = []
nan_kbs = []
nan_kubs = []

for i in range(len(datafiles)):
    fname = sys.argv[1] + datafiles[i]
    if os.stat(fname).st_size == 0:
        continue
    print("filename: %s" % fname)

    data = np.loadtxt(fname)

    file_txt = open(fname).read()

    start_kb_idx = file_txt.index('k_b-') + 4
    end_kb_idx = file_txt[start_kb_idx:].index(',') + start_kb_idx
    kb = float(file_txt[start_kb_idx:end_kb_idx])

    start_kub_idx = file_txt.index('k_ub-') + 5
    end_kub_idx = file_txt[start_kub_idx:].index(',') + start_kub_idx
    kub = float(file_txt[start_kub_idx:end_kub_idx])

    if "#EXIT SUCCESSFULLY" not in file_txt:
        print("Simulation either did not exit successfully or is not done.")
        nan_kbs.append(kb)
        nan_kubs.append(kub)
        continue

    if "#Last unbinding time" not in file_txt:
        print("Error, this is an old data file that doesn't have a #Last unbinding time.")
        continue

    start_last_unbinding_idx = file_txt.index("#Last unbinding time: ") + 22
    end_last_unbinding_idx = file_txt[start_last_unbinding_idx:].index('\n') + start_last_unbinding_idx
    last_unbinding_time = float(file_txt[start_last_unbinding_idx:end_last_unbinding_idx])

    start_runtime_idx = file_txt.index("#Runtime: ") + 10
    end_runtime_idx = file_txt[start_runtime_idx:].index('\n') + start_runtime_idx
    runtime = float(file_txt[start_runtime_idx:end_runtime_idx])

    if len(data) == 0: #or str(type(data[0])) == "<type 'numpy.float64'>":
        print("No steps in data file...")
        kbs.append(kb)
        kubs.append(kub)

        t_step.append(0)
        t_ob.append(runtime-last_unbinding_time)
        t_bb.append(last_unbinding_time)
        t_proc.append(float('inf'))

    else:
        print("Found steps in data file...")
        kbs.append(kb)
        kubs.append(kub)
        bind_times = np.array(np.concatenate(([0], data[:,1])))
        unbind_times = np.array(np.concatenate((data[:,0], [last_unbinding_time])))

        step_times = np.array(bind_times[1:] - bind_times[:-1])
        onebound_times = bind_times[1:] - unbind_times[:-1]
        bothbound_times = unbind_times - bind_times

        if (last_unbinding_time != runtime):
            onebound_times = np.append(onebound_times, runtime-last_unbinding_time)

        t_step.append(np.mean(step_times))
        t_ob.append(np.mean(onebound_times))
        t_bb.append(np.mean(bothbound_times))
        t_proc.append(t_step[-1]*t_bb[-1]/t_ob[-1])

    print("\n")

if len(kbs) == 0:
    print("Error, no stepping data in these data files!")
    exit(1)

min_kb = min(kbs)
max_kb = max(kbs)
min_kub = min(kubs)
max_kub = max(kubs)

maxy_tob = max(t_ob, 4.5*10**-4)
maxy_tbb = max(t_bb, 0.0595)

KB, KUB = np.meshgrid(kbs, kubs)

TBB = np.zeros((len(kbs), len(kubs)))

print("onebound times: ", t_ob)
print("bothbound times: ", t_bb)

plt.figure()
ax = plt.gca()

### TOB plot ###
raw_ratio = np.array(t_ob) / 4.5 / 1e-4
ratio = np.log10([s if s != 0 else 1e-100 for s in raw_ratio])

m = cm.ScalarMappable(cmap=cm.jet)
ratiomax = 10

m.set_array(np.linspace(-ratiomax, ratiomax, 100))

for i in range(len(ratio)):
    normalized_color = (ratio[i]+10)/20
    mycolor = m.cmap(normalized_color)
    plt.plot(kbs[i], kubs[i], '.', color=mycolor, markeredgecolor=mycolor)
    plt.annotate("%.2f" % ratio[i], xy=(kbs[i], kubs[i]), textcoords="offset points", xytext=(0,2), fontsize='2', arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

print ("onebound raw ratio: ", raw_ratio)
print ("onebound log ratios: ", ratio)

CB = plt.colorbar(m)

plt.plot(nan_kbs, nan_kubs, 'x', label="Incomplete or NaN-generating simulation")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$k_b$")
ax.set_ylabel("$k_{ub}$")
ax.set_xlim([0.01*min_kb, 100*max_kb])
ax.set_ylim([0.01*min_kub, 100*max_kub])
plt.legend(framealpha=0.5)
plt.title('$\log_{10}$(Onebound time / experimental) vs un/binding rates')

plt.savefig("plots/tob-rate-contour.pdf")

plt.figure()
ax = plt.gca()

### TBB plot ###
raw_ratio = np.array(t_bb) / (0.0595)
ratio = np.log10([s if s != 0 else 1e-100 for s in raw_ratio])

m = cm.ScalarMappable(cmap=cm.jet)
ratiomax = 10

m.set_array(np.linspace(-ratiomax, ratiomax, 100))
for i in range(len(ratio)):
    normalized_color = (ratio[i]+10)/20
    mycolor = m.cmap(normalized_color)
    plt.plot(kbs[i], kubs[i], '.', color=mycolor, markeredgecolor=mycolor)
    plt.annotate("{0:.2f}".format(np.log10(t_bb[i])), xy=(kbs[i], kubs[i]), textcoords="offset points", xytext=(0,2), fontsize='2', arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

CB = plt.colorbar(m)

plt.plot(nan_kbs, nan_kubs, 'x', label="Incomplete or NaN-generating simulation")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$k_b$")
ax.set_ylabel("$k_{ub}$")
ax.set_xlim([0.01*min_kb, 100*max_kb])
ax.set_ylim([0.01*min_kub, 100*max_kub])
plt.legend(framealpha=0.5)
plt.title('$\log_{10}$(Bothbound time / experimental) vs un/binding rates')

plt.savefig("plots/tbb-rate-contour.pdf")

# print("ratio: ", ratio)

# plt.figure()
# plt.gca().set_ylabel("$<t_{ob}> s$")
# plt.scatter(kbs, t_ob, color='b', s=50, label="$<t_{ob}>$")
# plt.plot([minx_kb*0.7, maxx_kb*1.3], [4.52*10**-4, 4.52*10**-4], color='g', label="$<t_{ob}>$ exp: $4.5*10^{-4}$ s")
# plt.scatter(empty_kbs, [0]*len(empty_kbs), color='k', marker="x", label="no steps", s=50)
# plt.title("Binding constants vs unbound time")
# plt.gca().set_xlabel("$k_b$" + " $s^{-1}$")
# plt.gca().set_xscale('log')
# plt.gca().set_xlim([minx_kb*0.8, maxx_kb*1.2])
# plt.gca().set_ylim([-0.1*maxy_tob, 1.5*maxy_tob])
# plt.legend(shadow=True)
# plt.savefig("plots/" + sys.argv[1] + "-tob-vs-kb-fit.pdf", bbox_inches='tight')

# plt.figure()
# plt.gca().set_ylabel("$<t_{bb}> s$")
# plt.scatter(kbs, t_bb, color='r', label="$<t_{bb}>$", s=50)
# plt.plot([minx_kb*0.7, maxx_kb*1.3], [0.0595, 0.0595], color='g', label="$<t_{bb}>$ exp: $0.0595$ s")
# plt.scatter(empty_kbs, [0]*len(empty_kbs), color='k', marker="x", label="no steps", s=50)
# plt.title("Binding constants vs bound time")
# plt.gca().set_xlabel("$k_b$" + " $s^{-1}$")
# plt.gca().set_xscale('log')
# plt.gca().set_xlim([minx_kb*0.8, maxx_kb*1.2])
# plt.gca().set_ylim([-0.1*maxy_tbb, 1.5*maxy_tbb])
# plt.legend(shadow=True)
# plt.savefig("plots/" + sys.argv[1] + "-tbb-vs-kb-fit.pdf", bbox_inches='tight')
