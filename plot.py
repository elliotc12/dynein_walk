#! /usr/bin/python

import numpy 
import matplotlib.pyplot as p

data = numpy.genfromtxt("data.txt", delimiter="\t")
p.plot(data[0:5], data[5:10])
p.plot(data[0:5], data[5:10], 'ro')
p.axis([0,50,0,50])
p.show()
