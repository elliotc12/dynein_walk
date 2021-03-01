import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../data")
import importlib
params = importlib.import_module("params")

import dynein.draw.cartoon as cartoon

from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter

def animate(i):
    movie_data = '../../../data/mc_movie_data.txt'
    time = np.loadtxt(movie_data)[:,1]
    PE_t = np.loadtxt(movie_data)[:,4]
    ta = np.loadtxt(movie_data)[:,7]

    bb = np.array([np.loadtxt(movie_data)[:,8], np.loadtxt(movie_data)[:,9]])
    bm = np.array([np.loadtxt(movie_data)[:,10], np.loadtxt(movie_data)[:,11]])
    t = np.array([np.loadtxt(movie_data)[:,12], np.loadtxt(movie_data)[:,13]])
    um = np.array([np.loadtxt(movie_data)[:,14], np.loadtxt(movie_data)[:,15]])
    ub = np.array([np.loadtxt(movie_data)[:,16], np.loadtxt(movie_data)[:,17]])

    plt.clf()
    plt.xlim(-50,50)
    plt.ylim(0,80)
    plt.text(-45, 75, 'At time = {}: PE tail = {}'.format(time[i], PE_t[i]))
    plt.text(-45, 72, 'ta = {}'.format(ta[i]*180/np.pi))

    cartoon.dyneinCircles(bb[0][i], bb[1][i], params.for_simulation['rb'],
                            bm[0][i], bm[1][i], params.for_simulation['rm'],
                            t[0][i], t[1][i], params.for_simulation['rt'],
                            'red', 0.5)

    cartoon.dyneinCircles(ub[0][i], ub[1][i], params.for_simulation['rb'],
                            um[0][i], um[1][i], params.for_simulation['rm'],
                            t[0][i], t[1][i], params.for_simulation['rt'],
                            'blue', 0.3)


def main():
    movie_data = '../../../data/mc_movie_data.txt'

    time = np.loadtxt(movie_data)[:,1]
    PE_t = np.loadtxt(movie_data)[:,4]
    ta = np.loadtxt(movie_data)[:,7]

    bb = np.array([np.loadtxt(movie_data)[:,8], np.loadtxt(movie_data)[:,9]])
    bm = np.array([np.loadtxt(movie_data)[:,10], np.loadtxt(movie_data)[:,11]])
    t = np.array([np.loadtxt(movie_data)[:,12], np.loadtxt(movie_data)[:,13]])
    um = np.array([np.loadtxt(movie_data)[:,14], np.loadtxt(movie_data)[:,15]])
    ub = np.array([np.loadtxt(movie_data)[:,16], np.loadtxt(movie_data)[:,17]])

    fig, ax = plt.subplots(figsize = (10,8))
    ax.set(xlim = (-50,50), ylim = (0,80))
    # ax.text(-40, 70, 'At time = {}: PE tail = {}'.format(time[0], PE_t[0]))
    cartoon.dyneinCircles(bb[0][0], bb[1][0], params.for_simulation['rb'],
                                bm[0][0], bm[1][0], params.for_simulation['rm'],
                                t[0][0], t[1][0], params.for_simulation['rt'],
                                'red', 0.5)
    cartoon.dyneinCircles(ub[0][0], ub[1][0], params.for_simulation['rb'],
                                um[0][0], um[1][0], params.for_simulation['rm'],
                                t[0][0], t[1][0], params.for_simulation['rt'],
                                'blue', 0.3)

    movie = FuncAnimation(fig, animate, frames = len(time))



    # print(len((time)))

    movie_save_fname =r"../../../plots/mc_plots/mc_movie.mp4"
    writervideo = FFMpegWriter(fps=10)

    movie.save(movie_save_fname, writer = writervideo)

    plt.show()

if __name__ == '__main__':
    main()
