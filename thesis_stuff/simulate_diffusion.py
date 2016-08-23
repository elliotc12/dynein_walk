#!/usr/bin/python2.7

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import os

runtime = 500
dt = 1

num_frames = 25

#radius = 1e-8
radius = 1000
cell_viscosity = 7e-4
gamma = radius * 6.0 * np.pi * cell_viscosity

cell_radius = 1e-17

kb = 1.38e-23
T = 310.15

plot_motor = True

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
        if int(t/dt) % 10 == 0:
            print "Diffusion position: " + str(pos) +  "completion: " + str(t / float(runtime))
        trajectory.append(pos.copy())
        t += dt
    return trajectory

def simulate_motor(init_position):
    t = 0
    pos = np.array(init_position)
    trajectory = [pos]
    while t < runtime:
        pos += (cell_radius*0.01/np.sqrt(3), cell_radius*0.01/np.sqrt(3), cell_radius*0.01/np.sqrt(3))
        # if (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2] >= cell_radius):
        #     pos = (cell_radius/np.sqrt(3), cell_radius/np.sqrt(3), cell_radius/np.sqrt(3))
        if int(t/dt) % 10 == 0:
            print "Motor position: " + str(pos) +  "completion: " + str(t / float(runtime))
        trajectory.append(pos[:])
        t += dt
    return trajectory

def make_gif(trajectory, motor_trajectory):
    os.system("rm -rf PNGs/*")
    os.system("mkdir PNGs")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.grid(False)
    ax.set_axis_off()

    skipiter = len(trajectory)/num_frames
    for t in range(len(trajectory)):
        if (t % skipiter == 0):

            u = np.linspace(0, np.pi, 100)
            v = np.linspace(0, np.pi, 100)

            x = 1e-17 * np.outer(np.cos(u), np.sin(v))
            y = 1e-17 * np.outer(np.sin(u), np.sin(v))
            z = 1e-17 * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b', alpha=0.2, linewidth=0.1)
        
            xs = [point[0] for (n, point) in enumerate(trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[0]
            ys = [point[1] for (n, point) in enumerate(trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[1]
            zs = [point[2] for (n, point) in enumerate(trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[2]
            ax.plot(xs, ys, zs)

            if (plot_motor):
                m_xs = [point[0] for (n, point) in enumerate(motor_trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[0]
                m_ys = [point[1] for (n, point) in enumerate(motor_trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[1]
                m_zs = [point[2] for (n, point) in enumerate(motor_trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[2]
                ax.plot(m_xs, m_ys, m_zs)

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

            ax.set_xlim([-1.5e-17, 1.5e-17])
            ax.set_ylim([-1.5e-17, 1.5e-17])
            ax.set_zlim([-1.5e-17, 1.5e-17])

            plt.savefig('PNGs/diffusion-%03d.png' % t)
            plt.cla()
            print "Making fig: %d of %d" % (t, len(trajectory))
    os.system("convert -delay 10 PNGs/diffusion-*.png diffusion.gif")
    

def main():
    init_position = np.array([0.0, 0.0, 0.0])
    trajectory = simulate_particle(init_position)
    motor_trajectory = simulate_motor(init_position)
    make_gif(motor_trajectory, motor_trajectory)

if __name__ == "__main__":
    main()
