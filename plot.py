#! /usr/bin/python

import math
import numpy
import thread
import time
import sys
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
config = numpy.loadtxt("config.txt")
data = numpy.genfromtxt("data.txt", delimiter="\t", skiprows=1)
plib.ion()
figure = plib.figure()
subplot = figure.add_subplot(111)
subplot.set_aspect("equal", adjustable="box")
subplot.set_xlim(-40,40)
subplot.set_ylim(-40,40)
line1, = subplot.plot(X, Y)
line2, = subplot.plot([X[0]], [Y[0]], 'ro')
line3, = subplot.plot([X[1]], [Y[1]], 'bo')
line4, = subplot.plot([X[2]], [Y[2]], 'go')
line5, = subplot.plot([X[3]], [Y[3]], 'bo')
line6, = subplot.plot([X[4]], [Y[4]], 'ro')

title_text = plib.text(-65, 45, 'Dynein BL-bound state')
pe_text = plib.text(-65, 40, 'PE: ')
ke_text = plib.text(-65, 35, 'KE: ')
t_text = plib.text(-65, -36, 't=:')

i = 0
print "Press enter to exit animation."
while (True):
	X[0] = data[i][0]
	X[1] = Ls*math.cos(data[i][2]) + X[0]
	X[2] = Lt*math.cos(data[i][3]) + X[1]
	X[3] = Lt*math.cos(-data[i][4]) + X[2]
	X[4] = Ls*math.cos(-data[i][5]) + X[3]

	Y[0] = data[i][1]
	Y[1] = Ls*math.sin(data[i][2]) + Y[0]
	Y[2] = Lt*math.sin(data[i][3]) + Y[1]
	Y[3] = Lt*math.sin(-data[i][4]) + Y[2]
	Y[4] = Ls*math.sin(-data[i][5]) + Y[3]
	
	line1.set_data(X,Y)
	line2.set_data(X[0], Y[0])
	line3.set_data(X[1], Y[1])
	line4.set_data(X[2], Y[2])
	line5.set_data(X[3], Y[3])
	line6.set_data(X[4], Y[4])
	
	pe_text.set_text('PE: ' + str(data[i][7]))
	ke_text.set_text('KE: ' + str(data[i][6]))

	t_text.set_text('t= ' + str(data[i][8]) + '/' + str(config[1]))
	
	plib.draw()
	
	if len(sys.argv) == 2:
		if sys.argv[1] == "step":
			raw_input("Hit enter to step.")
		elif sys.argv[1] == "slow":
			time.sleep(1)
	
	i += 1
	if (i >= len(data)):
		i = 0
		exit(0)
	if (flag and len(sys.argv) != 2) or (flag and sys.argv[1] != "step"): break
	
	#make text work
