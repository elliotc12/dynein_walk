import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import polygonData as pd



fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')




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
    md_array = pd.motorDomainArray(xm,ym,0.0025)
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


dyneinPolygon(0,0,1,1,2.5,2,'blue',1.0,ax) 
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.show()



