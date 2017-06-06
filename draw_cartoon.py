#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def draw_cartoon(near_foot_figure_coords, state, x_coords, y_coords, x_scaling, y_scaling, alpha):
    if (state == 0):
        figure_xs = near_foot_figure_coords[0] + x_scaling*(np.array(x_coords) - x_coords[0]) + x_coords[0]
    else:
        figure_xs = near_foot_figure_coords[0] + x_scaling*(np.array(x_coords) - x_coords[4]) + x_coords[4]

    figure_ys = near_foot_figure_coords[1] + y_scaling*np.array(y_coords)

    # print("figure_xs", figure_xs)
    # print("x_coords", x_coords)

    plt.plot([figure_xs[0], figure_xs[1]], [figure_ys[0], figure_ys[1]], color="black", zorder=4, linewidth=0.5, alpha=alpha)
    plt.plot([figure_xs[1], figure_xs[2]], [figure_ys[1], figure_ys[2]], color="black", zorder=4, linewidth=0.5, alpha=alpha)
    plt.plot([figure_xs[2], figure_xs[3]], [figure_ys[2], figure_ys[3]], color="black", zorder=1, linewidth=0.5, linestyle=":", alpha=alpha)
    plt.plot([figure_xs[3], figure_xs[4]], [figure_ys[3], figure_ys[4]], color="black", zorder=1, linewidth=0.5, linestyle=":", alpha=alpha)

    plt.plot([figure_xs[4]], [figure_ys[4]], marker='o', markeredgecolor='k', markeredgewidth=0.4, zorder=2, color="red",   markersize=6*x_scaling, alpha=alpha)
    plt.plot([figure_xs[3]], [figure_ys[3]], marker='o', markeredgecolor='k', markeredgewidth=0.4, zorder=3, color="white", markersize=18*x_scaling, alpha=alpha)
    plt.plot([figure_xs[2]], [figure_ys[2]], marker='o', markeredgecolor='k', markeredgewidth=0.4, zorder=4, color="white", markersize=6*x_scaling, alpha=alpha)
    plt.plot([figure_xs[1]], [figure_ys[1]], marker='o', markeredgecolor='k', markeredgewidth=0.4, zorder=5, color="white", markersize=18*x_scaling, alpha=alpha)
    plt.plot([figure_xs[0]], [figure_ys[0]], marker='o', markeredgecolor='k', markeredgewidth=0.4, zorder=6, color="blue",  markersize=6*x_scaling, alpha=alpha)
