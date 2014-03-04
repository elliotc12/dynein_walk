#! /usr/bin/python

import numpy 
import matplotlib.pyplot as p

data = numpy.genfromtxt("data.txt", delimiter="\t")
p1 = p.plot(data[0:5], data[5:10])

X[0] = data[0]
X[1] = X[0] + cos(data[2])
X[2] = X[1] + cos(data[3])
X[3] = X[2] + cos(data[4] - pi/2)
X[4] = X[3] + cos(data[5] - pi/2)

Y[0] = data[0]
Y[1] = Y[0] + sin(data[2])
Y[2] = Y[1] + sin(data[3])
Y[3] = Y[2] + sin(data[4] - pi/2)
Y[4] = Y[3] + sin(data[5] - pi/2)

p.plot(X, Y, 'ro')
p.text(5, 47, "Potential Energy: " + str(data[6]))
p.ylim(-10,50)
p.xlim(-10,50)
p.gca().set_aspect("equal", adjustable="box")
p.show()
