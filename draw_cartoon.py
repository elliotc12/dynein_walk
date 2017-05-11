#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def draw_cartoon(near_foot_figure_coords, x_coords, y_coords, x_scaling, y_scaling):
    figure_xs = near_foot_figure_coords[0] + x_scaling*(np.array(x_coords) - x_coords[0])
    figure_ys = near_foot_figure_coords[1] + y_scaling*np.array(y_coords)

    plt.plot([figure_xs[0], figure_xs[1]], [figure_ys[0], figure_ys[1]], color="black", zorder=1, linewidth=0.5)
    plt.plot([figure_xs[1], figure_xs[2]], [figure_ys[1], figure_ys[2]], color="black", zorder=1, linewidth=0.5)
    plt.plot([figure_xs[2], figure_xs[3]], [figure_ys[2], figure_ys[3]], color="black", zorder=1, linewidth=0.5, linestyle=":")
    plt.plot([figure_xs[3], figure_xs[4]], [figure_ys[3], figure_ys[4]], color="black", zorder=1, linewidth=0.5, linestyle=":")

    plt.plot([figure_xs[4]], [figure_ys[4]], marker='o', markeredgecolor='k', markeredgewidth=0.1, zorder=2, color="red",   markersize=2*x_scaling*50)
    plt.plot([figure_xs[3]], [figure_ys[3]], marker='o', markeredgecolor='k', markeredgewidth=0.1, zorder=3, color="white", markersize=6*x_scaling*50)
    plt.plot([figure_xs[2]], [figure_ys[2]], marker='o', markeredgecolor='k', markeredgewidth=0.1, zorder=4, color="white", markersize=2*x_scaling*50)
    plt.plot([figure_xs[1]], [figure_ys[1]], marker='o', markeredgecolor='k', markeredgewidth=0.1, zorder=5, color="white", markersize=6*x_scaling*50)
    plt.plot([figure_xs[0]], [figure_ys[0]], marker='o', markeredgecolor='k', markeredgewidth=0.1, zorder=6, color="blue",  markersize=2*x_scaling*50)

    plt.gca().add_patch(Rectangle((figure_xs[0], -0.15), 0.003, 0.3, facecolor='k', alpha=0.5, zorder=-1))
