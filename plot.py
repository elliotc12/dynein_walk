#! /usr/bin/env python2.7

import math
import numpy
import time
import signal
import sys
import matplotlib.pyplot as plt

def close_windows(*_):
  plt.close()
  sys.exit()

if len(sys.argv) < 2:
  print "Usage: ./plot.py speed=n [loop]"
  sys.exit(1)

if len(sys.argv) == 3 and sys.argv[2] == "loop":
  loop = True
else:
  loop = False

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

binding1, = plt.plot([X[0]], [Y[0]], 'ro')
motor1,   = plt.plot([X[1]], [Y[1]], 'bo')
tail,     = plt.plot([X[2]], [Y[2]], 'go')
motor2,   = plt.plot([X[3]], [Y[3]], 'bo')
binding2, = plt.plot([X[4]], [Y[4]], 'ro')

title_text = plt.text(-65, 45, 'State:')
pe_text = plt.text(-65, 40, 'PE: ')
ke_text = plt.text(-65, 35, 'KE: ')
t_text = plt.text(-65, -36, 't=:')

i = 0

frames = config[4] / config[3]
if len(data) < frames:
  unbound = True
else:
  unbound = False

signal.signal(signal.SIGINT, close_windows)

while i < len(data) or loop:
  if i >= len(data):
    if unbound:
      title_text.set_text('State: Unbound')
      plt.draw()
      print "Protein unbound!"
      plt.pause(3)
    i = 0

  X[0] = data[i][4]
  X[1] = data[i][6]
  X[2] = data[i][8]
  X[3] = data[i][10]
  X[4] = data[i][12]

  Y[0] = data[i][5]
  Y[1] = data[i][7]
  Y[2] = data[i][9]
  Y[3] = data[i][11]
  Y[4] = data[i][13]

  stalk1.set_data([ X[0], X[1] ], [ Y[0], Y[1] ])
  tail1.set_data([ X[1], X[2] ], [ Y[1], Y[2] ])
  tail2.set_data([ X[2], X[3] ], [ Y[2], Y[3] ])
  stalk2.set_data([ X[3], X[4] ], [ Y[3], Y[4] ])

  binding1.set_data(X[0], Y[0])
  motor1.set_data(X[1], Y[1])
  tail.set_data(X[2], Y[2])
  motor2.set_data(X[3], Y[3])
  binding2.set_data(X[4], Y[4])

  if (data[i][14] == 0):
    title_text.set_text('State: Nearbound')
    stalk1.set_linestyle('-')
    tail1.set_linestyle('-')
    tail2.set_linestyle('--')
    stalk2.set_linestyle('--')

  elif (data[i][14] == 1):
    title_text.set_text('State: Farbound')
    stalk1.set_linestyle('--')
    tail1.set_linestyle('--')
    tail2.set_linestyle('-')
    stalk2.set_linestyle('-')

  elif (data[i][14] == 2):
    title_text.set_text('State: Bothbound')
    stalk1.set_linestyle('-')
    tail1.set_linestyle('-')
    tail2.set_linestyle('--')
    stalk2.set_linestyle('--')

  pe_text.set_text('PE: ' + str(data[i][1]))
  ke_text.set_text('KE: ' + str(data[i][0]))

  t_text.set_text("Progress: {:3.1f}%".format(data[i][3]/config[4]*100))

  if len(sys.argv) >= 2:
    if sys.argv[1] == "step":
      raw_input("Hit enter to step.")
      i += 10
    elif sys.argv[1][0:6] == "speed=":
      i += float(sys.argv[1][6:])
      plt.pause(0.001)
  else:
    i += 10
    plt.pause(0.001)

  plt.draw()

if unbound:
  title_text.set_text('State: Unbound')
  plt.draw()
  print "Protein unbound!"
  plt.pause(3)
