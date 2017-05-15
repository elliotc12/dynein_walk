import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import re



coords_path = "circle_coords.txt"
circle_data = open(coords_path, 'r')

#print circle_data.read() 
data = []
x = []
y = []

for piece in circle_data:
    splitLine = re.split(r'\s|,', piece)
    for strng in splitLine:
        data.append(strng)
        try:
            f =float(strng)
            data.append(f)
        except ValueError:
            print "Cannot convert string to int", strng

circle_data.close() 

for i in range( len(data)):
    if i%2==0:
        x.append(data[i])
    else:
        y.append(data[i])



circle = np.zeros((len(x),2))
for i in range(0, len(x)):
    circle[i,0] = x[i]
    circle[i,1] = y[i]






fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')


plt.plot(x,y)
plt.show()



ax.add_patch(
    patches.Polygon(
            poly,
            color = 'black',
            alpha = 1.0 
        )
    )





# def Rectangle(x, y, width, height, c, a, ax):
#     ax.add_patch(
#         patches.Rectangle(
#             (x,y),
#             width,
#             height,
#             color=c,
#             alpha=a
#         )
#     )


# def Polygon(xy, c, a, ax):
#     ax.add_patch(
#         patches.Polygon(
#             xy,
#             color = c,
#             alpha = a 
#         )
#     )



def Lower(xL, yL, xU, yU, c, a, ax):
    length = np.sqrt((xU-xL)**2+(yU-yL)**2)
    r1 = 0.05*length
    r2 = 0.1*length

    # binding domain 
    ax.add_patch(
        patches.Circle(
            (xL, yL),
            radius = r1,
            color = c,
            alpha = a
        )
    )

    #leg 
    ax.add_patch(
        patches.Polygon(
            [[xL,yL],[xU,yU]],
            color = c ,
            alpha = a,
            lw = 0.65*length
        )
    )

    #motor domain 
    ax.add_patch(
        patches.Circle(
            (xU, yU),
            radius = r2,
            color = c,
            alpha = a
        )
    )
    

# Lower(0,0, 5,5, 'blue', 1, ax)



# plt.xlim(-10,10)
# plt.ylim(-10,10)
# plt.grid(True) 
plt.show()

