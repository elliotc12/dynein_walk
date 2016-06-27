#! /usr/bin/env python2.7

import math
import numpy
import time
import signal
import sys
import os
import matplotlib.pyplot as plt

def close_windows(*_):
  plt.close()
  sys.exit()

if len(sys.argv) < 2:
  print "Usage: ./plot.py speed=n [loop / step / save'/'savename]"
  sys.exit(1)

loop = False
step = False
savefig = False

if len(sys.argv) == 3:
  if sys.argv[2] == "loop":
    loop = True
  elif sys.argv[2] == "step":
    step = True
  elif sys.argv[2][0:4] == "save":
    savefig = True
    savename = sys.argv[2][5:]
    print "Saving to %s, wait until the program exits!" % savename

X = [0, 1, 2, 3, 4]
Y = [0, 1, 2, 3, 4]

config = numpy.loadtxt("config.txt")
data = numpy.loadtxt("data.txt", delimiter="\t", skiprows=1)
plt.ion()

if str(type(data[0])) == "<type 'numpy.float64'>":
       print "Very short animation!"
       close_windows()

gb = float(config[0]) # FIXME: make point radii based on these
gm = float(config[1])
gt = float(config[2])

ax = plt.gca()
ax.set_aspect("equal", adjustable="box")
ax.set_xlim(-40,40)
ax.set_ylim(-40,40)

microtubule = plt.plot([-400, 400], [-2, -2])
plt.xlabel('$x$ (nm)')
plt.setp(microtubule, color='c', alpha=0.8, linewidth=17.0)

stalk1, = plt.plot([ X[0], X[1] ], [ Y[0], Y[1] ], color="black")
tail1,  = plt.plot([ X[1], X[2] ], [ Y[1], Y[2] ], color="black")
tail2,  = plt.plot([ X[2], X[3] ], [ Y[2], Y[3] ], color="black")
stalk2, = plt.plot([ X[3], X[4] ], [ Y[3], Y[4] ], color="black")

force_line = [i for i in range(5)]
for i in range(5):
  force_line[i], = plt.plot([X[i], X[i]], [Y[i],Y[i]], 'r-')

binding1, = plt.plot([X[0]], [Y[0]], marker='o', color=(data[i][2],0,0), markersize=10)
motor1,   = plt.plot([X[1]], [Y[1]], marker='o', color=(data[i][3],0,0), markersize=20)
tail,     = plt.plot([X[2]], [Y[2]], marker='o', color=(data[i][4],0,0), markersize=5)
motor2,   = plt.plot([X[3]], [Y[3]], marker='o', color=(data[i][5],0,0), markersize=20)
binding2, = plt.plot([X[4]], [Y[4]], marker='o', color=(data[i][6],0,0), markersize=10)

title_text = plt.text(-65, 45, 'State:')
pe_text = plt.text(-65, 40, 'PE: ')
ke_text = plt.text(-65, 35, 'KE: ')
t_text = plt.text(-65, -36, 't=:')

i = 0
savefigframe = 0

frames = int(config[4] / config[3])

if len(data) < frames:
  unbound = True
else:
  unbound = False

signal.signal(signal.SIGINT, close_windows)

while i < len(data) or loop:
  if i >= len(data):
    if unbound:
      print "i: " + str(i)
      title_text.set_text('State: Unbound')
      stalk1.set_color('r')
      tail1.set_color('r')
      tail2.set_color('r')
      stalk2.set_color('r')
      plt.draw()
      print "Protein unbound!"
      plt.pause(3)
      stalk1.set_color('black')
      tail1.set_color('black')
      tail2.set_color('black')
      stalk2.set_color('black')
    i = 0

  nba_scaling = min(max(1 - (data[i][2]) / (0.5*config[5]), 0), 1) # scale based on PE/kbT ratio forced between 0-1
  nma_scaling = min(max(1 - (data[i][3]) / (0.5*config[5]), 0), 1)
  ta_scaling  = min(max(1 - (data[i][4]) / (0.5*config[5]), 0), 1)
  fma_scaling = min(max(1 - (data[i][5]) / (0.5*config[5]), 0), 1)
  fba_scaling = min(max(1 - (data[i][6]) / (0.5*config[5]), 0), 1)

  binding1.set_color((1, nba_scaling, nba_scaling))
  motor1.set_color((1, nma_scaling, nma_scaling))
  tail.set_color((1, ta_scaling, ta_scaling))
  motor2.set_color((1, fma_scaling, fma_scaling))
  binding2.set_color((1, fba_scaling, fba_scaling))

  X[:] = data[i][7:16:2]
  Y[:] = data[i][8:17:2]
  Fx = data[i][17:27:2]
  Fy = data[i][18:28:2]
  # print 'Fx', Fx
  # print 'Fy', Fy

  stalk1.set_data([ X[0], X[1] ], [ Y[0], Y[1] ])
  tail1.set_data([ X[1], X[2] ], [ Y[1], Y[2] ])
  tail2.set_data([ X[2], X[3] ], [ Y[2], Y[3] ])
  stalk2.set_data([ X[3], X[4] ], [ Y[3], Y[4] ])

  force_scaling = 1e-9
  for j in range(5):
    force_line[j].set_data([X[j], X[j]+force_scaling*Fx[j]],
                           [Y[j], Y[j]+force_scaling*Fy[j]])

  binding1.set_data(X[0], Y[0])
  motor1.set_data(X[1], Y[1])
  tail.set_data(X[2], Y[2])
  motor2.set_data(X[3], Y[3])
  binding2.set_data(X[4], Y[4])

  if (data[i][0] == 0):
    title_text.set_text('State: Nearbound')
    stalk1.set_linestyle('-')
    tail1.set_linestyle('-')
    tail2.set_linestyle('--')
    stalk2.set_linestyle('--')

  elif (data[i][0] == 1):
    title_text.set_text('State: Farbound')
    stalk1.set_linestyle('--')
    tail1.set_linestyle('--')
    tail2.set_linestyle('-')
    stalk2.set_linestyle('-')

  elif (data[i][0] == 2):
    title_text.set_text('State: Bothbound')
    stalk1.set_linestyle('-')
    tail1.set_linestyle('-')
    tail2.set_linestyle('--')
    stalk2.set_linestyle('--')

  pe_text.set_text('PE: %.2f' % (data[i][2]+data[i][3]+data[i][4]+data[i][5]+data[i][6]))

  t_text.set_text("Progress: {:3.1f}%".format(data[i][1]/config[4]*100))

  if step:
      k = raw_input("Hit enter to step. [b=back,s=small]")
      if k == 'b':
        i -= 10
      elif k == 's':
        i += 1
      else:
        i += 10
  elif len(sys.argv) >= 2 and sys.argv[1][0:6] == "speed=":
      i += float(sys.argv[1][6:])
      plt.pause(0.001)
  else:
    i += 10
    plt.pause(0.001)
  savefigframe += 1

  if savefig:
    fname = 'PNGs/%s-%03d.png' % (savename, savefigframe)
    plt.savefig(fname)

  plt.draw()

if unbound:
  title_text.set_text('State: Unbound')
  plt.draw()
  print "Protein unbound!"
  plt.pause(3)

if savefig:
  os.system("convert -delay 10 PNGs/%s-*.png movies/%s.gif" % (savename, savename))
  os.system("mencoder PNGs/%s-*.png -mf type=png:fps=10 -ovc lavc"
            " -lavcopts vcodec=wmv2 -oac copy -o movies/%s.mpg" % (savename, savename))
  os.system("rm PNGs/*")
