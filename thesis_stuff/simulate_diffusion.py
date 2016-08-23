#!/usr/bin/python2.7

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import os

runtime = 50
dt = 1

#radius = 1e-8
radius = 100
cell_viscosity = 7e-4
gamma = radius * 6.0 * np.pi * cell_viscosity

kb = 1.38e-23
T = 310.15

def get_cell_potential_gradient(pos):
    return np.array([0, 0, 0])

def simulate_particle(init_position):
    t = 0
    pos = np.array(init_position)
    trajectory = [pos]
    while t < runtime:
        F = -get_cell_potential_gradient(pos)
        R = np.random.normal(0, 2*kb*T*gamma/dt, 3)
        v = F/gamma + R
        pos += v*dt
        if int(t/dt) % 10000 == 0:
            print "Position: " + str(pos) +  "completion: " + str(t / runtime)
        trajectory.append(pos.copy())
        t += dt
    return trajectory

def make_gif(trajectory):
    os.system("rm PNGs/*")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.grid(False)
    ax.set_axis_off()

    for t in range(len(trajectory)):
        u = np.linspace(0, np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = 1e-18 * np.outer(np.cos(u), np.sin(v))
        y = 1e-18 * np.outer(np.sin(u), np.sin(v))
        z = 1e-18 * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b', alpha=0.2, linewidth=0.1)
        
        xs = [point[0] for (n, point) in enumerate(trajectory) if n < t] # = (trajectory[0:t])[0]
        ys = [point[1] for (n, point) in enumerate(trajectory) if n < t] # = (trajectory[0:t])[1]
        zs = [point[2] for (n, point) in enumerate(trajectory) if n < t] # = (trajectory[0:t])[2]
        ax.plot(xs, ys, zs)

        ax.grid(False)
        ax.set_axis_off()

        # Get rid of the panes
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

        # Get rid of the spines
        ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

        ax.set_xticks([]) 
        ax.set_yticks([]) 
        ax.set_zticks([])

        ax.set_xlim([-1.5e-18, 1.5e-18])
        ax.set_ylim([-1.5e-18, 1.5e-18])
        ax.set_zlim([-1.5e-18, 1.5e-18])

        plt.savefig('PNGs/diffusion-%03d.png' % t)
        plt.cla()
        print "Making fig: %d of %d" % (t, len(trajectory))
    os.system("convert -delay 10 PNGs/diffusion-*.png diffusion.gif")
    

def main():
    init_position = np.array([0.0, 0.0, 0.0])
    trajectory = simulate_particle(init_position)
    make_gif(trajectory)

if __name__ == "__main__":
    main()
