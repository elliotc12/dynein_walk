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
data = numpy.genfromtxt("data.txt", delimiter="\t", skiprows=1)
plt.ion()

gb = float(config[0]) # FIXME: make point radii based on these
gm = float(config[1])
gt = float(config[2])

ax = plt.gca()
ax.set_aspect("equal", adjustable="box")
ax.set_xlim(-40,40)
ax.set_ylim(-40,40)

microtubule = plt.plot([-40, 40], [-2, -2])
plt.setp(microtubule, color='c', alpha=0.8, linewidth=17.0)

line1, = plt.plot(X, Y)
line2, = plt.plot([X[0]], [Y[0]], 'ro')
line3, = plt.plot([X[1]], [Y[1]], 'bo')
line4, = plt.plot([X[2]], [Y[2]], 'go')
line5, = plt.plot([X[3]], [Y[3]], 'bo')
line6, = plt.plot([X[4]], [Y[4]], 'ro')

title_text = plt.text(-65, 45, 'State:')
pe_text = plt.text(-65, 40, 'PE: ')
ke_text = plt.text(-65, 35, 'KE: ')
t_text = plt.text(-65, -36, 't=:')

i = 0

signal.signal(signal.SIGINT, close_windows)

while i < len(data) or loop:
  if i >= len(data):
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
  
  line1.set_data(X,Y)
  line2.set_data(X[0], Y[0])
  line3.set_data(X[1], Y[1])
  line4.set_data(X[2], Y[2])
  line5.set_data(X[3], Y[3])
  line6.set_data(X[4], Y[4])
  
  if (data[i][4] == 0):
    title_text.set_text('State: Leftbound')
  elif (data[i][4] == 1):
    title_text.set_text('State: Rightbound')
  elif (data[i][4] == 2):
    title_text.set_text('State: Bothbound')
  
  pe_text.set_text('PE: ' + str(data[i][1]))
  ke_text.set_text('KE: ' + str(data[i][0]))

  t_text.set_text('t= ' + str(data[i][3]) + '/' + str(config[4]))
  
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
