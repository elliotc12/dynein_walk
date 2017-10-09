from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import re

## outputs file with data to draw motor domain using outer_coords.txt ##


outer_path = "outer_coords.txt"

actual_radius = 7 # nanometers

with open(outer_path, 'r') as outer_data:
    outer_disps = []  #motor domain displacements for outtermost path

    for piece in outer_data:
        splitLine = re.split(r'\s|,', piece)
        for strng in splitLine:
            try:
                f =float(strng)
                outer_disps.append(f)
            except ValueError:
                print("Cannot convert string to int", strng)
mdo_x = [0]     #motor domain outer path x displacements
mdo_y = [0]     #motor domain outer path y displacements
for i in range( len(outer_disps)):
    if i%2==0:
        mdo_y.append(mdo_y[-1] + outer_disps[i])
    else:
        mdo_x.append(mdo_x[-1] + outer_disps[i])

mdo_array = np.zeros((len(mdo_x),2))
for j in range(len(mdo_x)):
    mdo_array[j,0] = mdo_x[j]
    mdo_array[j,1] = mdo_y[j]

center_mass = mdo_array.sum(axis=0)/len(mdo_array)
#print(center_mass)

for i in range(len(mdo_array)):
    mdo_array[i,:] -= center_mass

radius = 0
for i in range(len(mdo_array)):
    radius = max(radius, np.sqrt(mdo_array[i,0]**2 + mdo_array[i,1]**2))
mdo_array *= actual_radius/radius

with open("motor_domain.py", "w") as f:
    f.write("""from numpy import array
""")
    f.write("array = %s" % repr(mdo_array))
