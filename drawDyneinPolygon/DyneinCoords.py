import numpy as np
import matplotlib.pyplot as plt
import re

coords_path = "outer_coords.txt"
with open(coords_path, 'r') as poly_data:
    #print circle_data.read() 
    data = []

    for piece in poly_data:
        splitLine = re.split(r'\s|,', piece)
        for strng in splitLine:
            try:
                f =float(strng)
                data.append(f)
            except ValueError:
                print "Cannot convert string to int", strng

x = [0]
y = [0]
for i in range( len(data)):
    if i%2==0:
        y.append(y[-1] + data[i])
    else:
        x.append(x[-1] + data[i])



poly = np.zeros((len(x),2))
for i in range(0, len(x)):
    poly[i,0] = x[i]
    poly[i,1] = y[i]

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')

plt.plot(x,y, '.-')
plt.show()
