import numpy as np
import matplotlib.pyplot as plt
import re





outer_path = "outer_coords.txt"
tail_path = "tailCoords.txt"
# with open(outer_path, 'r') as outer_data: 
#     outer_disps = []  #motor domain displacements for outtermost path

#     for piece in outer_data:
#         splitLine = re.split(r'\s|,', piece)
#         for strng in splitLine:
#             try:
#                 f =float(strng)
#                 outer_disps.append(f)
#             except ValueError:
#                 print "Cannot convert string to int", strng


 
    
  
# mdo_x = [0]     #motor domain outer path x displacements
# mdo_y = [0]     #motor domain outer path y displacements
# for i in range( len(outer_disps)):
#     if i%2==0:
#         mdo_y.append(mdo_y[-1] + outer_disps[i])
#     else:
#         mdo_x.append(mdo_x[-1] + outer_disps[i])


def motorDomainArray(x,y,s):
    with open(outer_path, 'r') as outer_data:
        outer_disps = []  #motor domain displacements for outtermost path

        for piece in outer_data:
            splitLine = re.split(r'\s|,', piece)
            for strng in splitLine:
                try:
                    f =float(strng)
                    outer_disps.append(s*f)
                except ValueError:
                    print "Cannot convert string to int", strng
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
    return mdo_array

def tail(x1,y1,x2,y2):
    pass 

# inner_path = "inner_coords.txt"
# with open(inner_path, 'r') as inner_data:
