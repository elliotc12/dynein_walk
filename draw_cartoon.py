#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def draw_cartoon(tail_figure_coords, x_coords, y_coords, x_scaling, y_scaling):
    figure_xs = tail_figure_coords[0] + x_scaling*np.array(x_coords)
    figure_ys = tail_figure_coords[1] + y_scaling*np.array(y_coords)

    plt.plot([figure_xs[0], figure_xs[1]], [figure_ys[0], figure_ys[1]], color="black")
    plt.plot([figure_xs[1], figure_xs[2]], [figure_ys[1], figure_ys[2]], color="black")
    plt.plot([figure_xs[2], figure_xs[3]], [figure_ys[2], figure_ys[3]], color="black")
    plt.plot([figure_xs[3], figure_xs[4]], [figure_ys[3], figure_ys[4]], color="black")

    plt.plot([figure_xs[0]], [figure_ys[0]], marker='o', markeredgecolor='k', color="white", markersize=1)
    plt.plot([figure_xs[1]], [figure_ys[1]], marker='o', markeredgecolor='k', color="white", markersize=12)
    plt.plot([figure_xs[2]], [figure_ys[2]], marker='o', markeredgecolor='k', color="red",   markersize=8)
    plt.plot([figure_xs[3]], [figure_ys[3]], marker='o', markeredgecolor='k', color="white", markersize=12)
    plt.plot([figure_xs[4]], [figure_ys[4]], marker='o', markeredgecolor='k', color="white", markersize=1)
