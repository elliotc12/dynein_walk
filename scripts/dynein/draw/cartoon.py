import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import dynein.draw.motor_domain as md
import dynein.draw.tail as tail

physical_stalk_length = 22.1 # nm

def draw(ax, x_coords, y_coords, alpha=1):
    dyneinPolygon(x_coords[0], y_coords[0], x_coords[1], y_coords[1], x_coords[2], y_coords[2], 'blue', alpha, ax)
    dyneinPolygon(x_coords[4], y_coords[4], x_coords[3], y_coords[3], x_coords[2], y_coords[2], 'red', alpha, ax)

def _draw_circle(x,y,R, color, alpha):
    angles = np.linspace(0, np.pi, 1000)
    plt.fill_between(x + R*np.cos(angles), y + R*np.sin(angles), y - R*np.sin(angles), color=color, alpha=alpha)
    angles = np.linspace(0, 2*np.pi, 1000)
    xs = np.linspace(-R, R, 2000)
    y_outer = np.sqrt(R**2 - xs**2)
    y_inner = np.sqrt((0.95*R)**2 - xs**2)
    y_inner[y_inner != y_inner] = 0
    plt.fill_between(x + xs, y + y_outer, y + y_inner, color='black', alpha=alpha)
    plt.fill_between(x + xs, y - y_outer, y - y_inner, color='black', alpha=alpha)

def dyneinPolygon(xb, yb, xm, ym,xt,yt, c, a, ax):
    stalk_length = np.sqrt((xm-xb)**2+(ym-yb)**2)
    r1 = 0.1*stalk_length

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
            lw = 0.4*stalk_length  #scaling not physically significant. Just for visual appeal
        )
    )

    #tail
    t = np.array(tail.array)
    L = np.sqrt((xt-xm)**2+(yt-ym)**2)
    t[:,0] = (0.5*L)*t[:,0]
    t[:,1] = (0.3*L)*t[:,1]
    theta = np.arctan2((yt-ym),(xt-xm))
    rot = np.matrix([[np.cos(theta), -np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    for i in range(0, len(t[:,0])):
        t[i] = np.dot(rot, t[i])
    t[:,0] = t[:,0] + xm
    t[:,1] = t[:,1] + ym
    ax.add_patch(
        patches.Polygon(
            t,
            color = c, # "black",
            alpha = a
            )
        )

    #motor domain
    motor_domain = np.array(md.array)*stalk_length/physical_stalk_length
    motor_domain[:,0] = motor_domain[:,0] + xm
    motor_domain[:,1] = motor_domain[:,1] + ym
    ax.add_patch(
        patches.Polygon(
            motor_domain,
            color = c,
            alpha = a
            )
        )

def dyneinCircles(xb, yb, Rb, xm, ym, Rm, xt, yt, Rt, color, alpha, ax):
    stalk_length = np.sqrt((xm-xb)**2+(ym-yb)**2)

    for (x,y,R) in [(xb,yb,Rb), (xm,ym,Rm), (xt,yt,Rt)]:
        # circle for domain
        _draw_circle(x,y,R, color, alpha)
    plt.plot([xb,xm,xt], [yb,ym,yt], color='black', alpha=alpha)

if __name__ == "__main__":
    dyneinPolygon(0,5,10,25,30,35,'blue',1.0,ax)
    dyneinPolygon(10,0,11,21,30,35,'red',1.0,ax)

    plt.xlim(-10,50)
    plt.ylim(-10,50)

    plt.show()
