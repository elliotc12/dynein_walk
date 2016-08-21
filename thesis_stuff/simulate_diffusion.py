#!/usr/bin/python2.7

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np



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
        if int(t/dt) % 100 == 0:
            print "Position: " + str(pos) +  "completion: " + str(t / runtime)
        trajectory.append(pos.copy())
        t += dt
    return trajectory

def make_gif(trajectories):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for trajectory in trajectories:
        xs = [point[0] for point in trajectory]
        ys = [point[1] for point in trajectory]
        zs = [point[2] for point in trajectory]
        ax.plot(xs, ys, zs)

    # ax.set_xlim([-radius, radius])
    # ax.set_ylim([-radius, radius])
    # ax.set_zlim([-radius, radius])
    plt.title("Diffusion trajectories")
    plt.show()

def main():
    trajectories = []
    for n in range(N):
        init_position = np.array([0.0, 0.0, 0.0])
        trajectory = simulate_particle(init_position)
        trajectories.append(trajectory)
    make_gif(trajectories)

if __name__ == "__main__":
    main()
