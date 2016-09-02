#!/usr/bin/python2.7

from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
import os

dt = .00001

num_frames = 300
gif_length = 1500 # cs

kb = 1.38e-23 # J / K
T = 310.15 # K

radius = 1e-8
cell_viscosity = 1.5e-3 # J*s/m^3
gamma = radius * 6.0 * np.pi * cell_viscosity # J*s/m^2 = kg/s

D = kb * T / gamma # m / s

dFdr = 1e-12

kinesin_velocity = 2e-6 # m/s

cell_radius = 1e-5
rod_len_ratio = 9

plot_motor = True
mode = "rod"

def get_cell_potential_gradient(pos):
    if (mode == "sphere"):
        return np.array([0.0, 0.0, 0.0])
    elif (mode == "rod"):
        grad = np.array([0.0, 0.0, 0.0])
        if (pos[0] < 0):
            grad[0] += pos[0]*dFdr/np.abs(pos[0])
        if (pos[1]*pos[1] + pos[2]*pos[2] > cell_radius*cell_radius):
            grad[1] = grad[1] + pos[1]*dFdr/np.abs(pos[1])
            grad[2] = grad[2] + pos[2]*dFdr/np.abs(pos[2])
        return grad

def magnitude(pos):
    return np.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2])

def simulation_finished(pos):
    if (mode == "sphere"):
        return (magnitude(pos) >= cell_radius)
    elif (mode == "rod"):
        return (pos[0] >= rod_len_ratio*cell_radius)

def simulate_particle(init_position):
    t = 0
    pos = np.array(init_position)
    trajectory = [pos]
    np.seterr(all='raise')
    while not simulation_finished(pos):
        F = -get_cell_potential_gradient(pos)
        R = np.random.normal(0, np.sqrt(2*D/dt), 3)
        try:
            v = F/gamma + R
        except FloatingPointError:
            print F, gamma, R
        v = F/gamma + R
        pos += v*dt
        if int(t/dt) % 10000 == 0:
            print "pos: " + str(pos)
        trajectory.append(pos.copy())
        t += dt
    print "final pos: " + str(pos)
    return trajectory

def simulate_motor(init_position):
    t = 0
    pos = np.array(init_position)
    trajectory = [pos]
    while magnitude(pos) < rod_len_ratio*cell_radius:
        pos += np.array([1, 0, 0])*kinesin_velocity*dt
        if (magnitude(pos) >= rod_len_ratio*cell_radius):
            pos = np.array([cell_radius*rod_len_ratio, 0, 0])
        if int(t/dt) % 10000 == 0:
            print "Motor position: " + str(pos)
        trajectory.append(pos.copy())
        t += dt
    return trajectory

def make_gif(trajectory, motor_trajectory):
    os.system("rm -rf PNGs/*")
    os.system("mkdir -p PNGs")
    fig = plt.figure(figsize=(12, 12), dpi=80)
    ax = fig.add_subplot(111, projection='3d')
    ax.grid(False)
    ax.set_axis_off()

    idx = 0
    skipiter = max(len(motor_trajectory)/num_frames, 1)
    for t in range(len(motor_trajectory)):
        if (t % skipiter == 0):

            phi = np.linspace(0, np.pi, 100)
            theta = np.linspace(0, np.pi, 100)

            if (mode == "sphere"):
                x = cell_radius * np.outer(np.cos(theta), np.sin(phi))
                y = cell_radius * np.outer(np.sin(theta), np.sin(phi))
                z = cell_radius * np.outer(np.ones(np.size(theta)), np.cos(phi))
                ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b', alpha=0.1, linewidth=0.01)
                ax.plot_surface(x, -y, z, rstride=4, cstride=4, color='b', alpha=0.1, linewidth=0.01)

            elif (mode == "rod"):
                x = np.linspace(0, rod_len_ratio*cell_radius, 100)
                z = np.linspace(-cell_radius, cell_radius, 100)
                x, z = np.meshgrid(x, z)
                y = np.sqrt(cell_radius*cell_radius - z*z)
                ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b', alpha=0.1, linewidth=0.01)
                ax.plot_surface(x, -y, z, rstride=4, cstride=4, color='b', alpha=0.1, linewidth=0.01)

                cap = patches.Circle((0, 0), cell_radius, color='b', alpha=0.1)
                ax.add_patch(cap)
                art3d.pathpatch_2d_to_3d(cap, z=0, zdir="x")

                cap = patches.Circle((0, 0), cell_radius, color='b', alpha=0.1)
                ax.add_patch(cap)
                art3d.pathpatch_2d_to_3d(cap, z=rod_len_ratio*cell_radius, zdir="x")

            xs = [point[0] for (n, point) in enumerate(trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[0]
            ys = [point[1] for (n, point) in enumerate(trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[1]
            zs = [point[2] for (n, point) in enumerate(trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[2]
            ax.plot(xs, ys, zs, linewidth=0.5)
            if len(xs) > 0:
                ax.scatter([xs[-1]], [ys[-1]], [zs[-1]], s=10, color='r')

            if (plot_motor):
                m_xs = [point[0] for (n, point) in enumerate(motor_trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[0]
                m_ys = [point[1] for (n, point) in enumerate(motor_trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[1]
                m_zs = [point[2] for (n, point) in enumerate(motor_trajectory) if (n < t and n > 0)] # = (trajectory[0:t])[2]
                ax.plot(m_xs, m_ys, m_zs, linewidth=0.5)
                if len(m_xs) > 0:
                    ax.scatter([m_xs[-1]], [m_ys[-1]], [m_zs[-1]], s=10, color='r')

            if (mode == "sphere"):
                ax.text(cell_radius, -2*cell_radius, 0, "seconds: %f" % (t*dt))
            elif (mode == "rod"):
                ax.text(4*cell_radius, -4*cell_radius, 0, "seconds: %f" % (t*dt))

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

            if (mode == "sphere"):
                ax.set_xlim([-cell_radius, cell_radius])
                ax.set_ylim([-cell_radius, cell_radius])
                ax.set_zlim([-cell_radius, cell_radius])
            elif (mode == "rod"):
                ax.set_xlim([0, rod_len_ratio*cell_radius])
                ax.set_ylim([-rod_len_ratio*cell_radius/2, rod_len_ratio*cell_radius/2])
                ax.set_zlim([-rod_len_ratio*cell_radius/2, rod_len_ratio*cell_radius/2])

            plt.savefig('PNGs/diffusion-%03d.png' % idx)
            idx += 1
            plt.cla()
            print "Making fig: %d of %d" % (t, len(motor_trajectory))
    os.system("convert -delay " + str(gif_length/num_frames) + " PNGs/diffusion-*.png diffusion.gif")

def main():
    init_position = np.array([0.0, 0.0, 0.0])
    trajectories = [0]*1
    for n in range(1):
        trajectories[n] = simulate_particle(init_position)
    median_len = np.median([len(t) for t in trajectories])
    median_trajectory = [t for t in trajectories if len(t) == median_len][0]
    motor_trajectory = simulate_motor(init_position)
    make_gif(median_trajectory, motor_trajectory)

if __name__ == "__main__":
    print "D", D
    main()
