#! /usr/bin/python

import math
import numpy
import thread
import matplotlib.pyplot as plib

def wait(flag):
	raw_input()
	flag.append(1)

Ls = 10
Lt = 10
X = [0, 1, 2, 3, 4]
Y = [0, 1, 2, 3, 4]
flag = []
thread.start_new_thread(wait, (flag,))
data = numpy.genfromtxt("data.txt", delimiter="\t", skiprows=1)
plib.ion()
figure = plib.figure()
subplot = figure.add_subplot(111)
subplot.set_aspect("equal", adjustable="box")
subplot.set_xlim(-10,50)
subplot.set_ylim(-10,50)
line1, = subplot.plot(X, Y)
line2, = subplot.plot(X, Y, 'ro')

i = 0
while (True):
	X[0] = data[i][0]
	X[1] = Ls*math.cos(data[i][2]) + X[0]
	X[2] = Lt*math.cos(data[i][3]) + X[1]
	X[3] = Lt*math.cos(-1*data[i][4]) + X[2]
	X[4] = Ls*math.cos(-1*data[i][5]) + X[3]

	Y[0] = data[i][1]
	Y[1] = Ls*math.sin(data[i][2]) + Y[0]
	Y[2] = Lt*math.sin(data[i][3]) + Y[1]
	Y[3] = Lt*math.sin(-1*data[i][4]) + Y[2]
	Y[4] = Ls*math.sin(-1*data[i][5]) + Y[3]
	
	line1.set_data(X,Y)
	line2.set_data(X,Y)
	
	plib.draw()
	
	i += 1
	if (i >= len(data)):
		i = 0
	if flag: break
