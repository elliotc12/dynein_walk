import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import motorDomain as md



fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
M =md.motorDomainArray(0.5,0.5,0.001)
# ax.add_patch(
#     patches.Polygon(
#         M*0.001,
#         color='blue',
#         alpha= 1.0
#         )
#     )
plt.plot(M[:,0],M[:,1])


plt.show()

def dyneinPolygon(xb, yb, xm, ym,xt,yt, c, a, ax):
    length = np.sqrt((xm-xb)**2+(ym-yb)**2)
    r1 = 0.05*length
    r2 = 0.1*length

    # binding domain 
    ax.add_patch(
        patches.Circle(
            (xb, yb),
            radius = r1,
            color = c,
            alpha = a
        )
    )

    #leg 
    ax.add_patch(
        patches.Polygon(
            [[xb,yb],[xm,ym]],
            color = c ,
            alpha = a,
            lw = 0.65*length
        )
    )

    #motor domain
    md_array = md.motorDomainArray(xm,ym,0.001)
    ax.add_patch(
        patches.Polygon(
            md_array,
            color = c,
            alpha = a
        )
    )
    #tail
    ax.add_patch(
        patches.Polygon(
            [[xm,ym],[xt,yt]],
            color = c,
            alpha = a,
            lw = 0.65*np.sqrt((xt-xm)**2+(yt-ym)**2)
            )
        )






