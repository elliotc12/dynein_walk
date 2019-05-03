import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import sys
sys.path.append("../data")
import importlib
import argparse
import subprocess
import bb_energy_distribution

params = importlib.import_module("params")

rate_unbinding_leading = []                        # Leading (Far) Unbinding Rates
rate_unbinding_trailing = []                # Trailing (Near) Unbinding Rates

angles = [[] for i in range(2)]

C =  params.for_simulation['exp-unbinding-constant']         # exponential binding constant from paper_params.py April 12


L = 24         # Length

x = 0
while x < 9:
# Making random motor angles
        nma = np.random.uniform(0, 2*np.pi)
        fma = np.random.uniform(0, 2*np.pi)

        dynein = bb_energy_distribution.DyneinBothBound(nma, fma, params, L)

        # Checking if energy is nan
        if np.isnan(dynein.E_total) == True:
                continue
        else:
                angles[0].append(nma)
                angles[1].append(fma)
                x = x + 1

def plot_bb_energy_distribution(self, d1, d2, d3, d4, d5, d6, d7, d8, d9):
    """Plot the energy distribution for the both bound configuration given an
    array of motor angles and an initial displacement.
    """

    fig = plt.figure()

    # make contourf graph
    ax1 = fig.add_subplot(1, 2, 1)
    energyPlot = ax1.contourf(self.nma, self.fma, self.E_total, 100)
    contour = ax1.contour(self.nma, self.fma, self.E_total, np.arange(1, 30, 1), colors='w', linewidth=10)
    ax1.set_xlabel(r'$\theta_{nm}$')
    ax1.set_ylabel(r'$\theta_{fm}$')
    ax1.set_title('Total Energy distribution for L={0}'.format(self.L))
    ax1.set_xlim(0-0.1, 2*np.pi+0.1)
    ax1.set_ylim(0-0.1, 2*np.pi+0.1)
    cb = plt.colorbar(energyPlot)
    cb.set_label(r"Energy [$k_BT$]")
    cb.set_ticks(np.arange(0, 31, 5))
    cb.add_lines(contour)

    # find the extrema
    j_max, i_max, j_min, i_min = self.find_energy_extrema()
    ax1.scatter(self.nma[i_min, j_min], self.fma[i_min, j_min], color='black')
    ax1.scatter(d1.nma, d1.fma, color='red')
    ax1.scatter(d2.nma, d2.fma, color='orange')
    ax1.scatter(d3.nma, d3.fma, color='lime')
    ax1.scatter(d4.nma, d4.fma, color='blue')
    ax1.scatter(d5.nma, d5.fma, color='gray')
    ax1.scatter(d6.nma, d6.fma, color='magenta')
    ax1.scatter(d7.nma, d7.fma, color='green')
    ax1.scatter(d8.nma, d8.fma, color='purple')
    ax1.scatter(d9.nma, d9.fma, color='sienna')

    ax2 = fig.add_subplot(1, 2, 2)
    x_coords_min = [self.r_nb[0,i_min, j_min],
                    self.r_nm[0,i_min, j_min],
                    self.r_t[0,i_min, j_min],
                    self.r_fm[0,i_min, j_min],
                    self.r_fb[0,i_min, j_min]]
    y_coords_min = [self.r_nb[1,i_min, j_min],
                    self.r_nm[1,i_min, j_min],
                    self.r_t[1,i_min, j_min],
                    self.r_fm[1,i_min, j_min],
                    self.r_fb[1,i_min, j_min]]


    ax2.plot(x_coords_min, y_coords_min, color='black', label="min")
    ax2.axis('off')
    ax2.axis('equal')
    ax2.legend()


def plot_bb_figures(self, dynein_color):
    """Plot just the figure of dynein for the both bound configuration given an
    array of motor angles and an initial displacement.
    """
    fig = plt.figure()

    ax = fig.add_subplot(1, 1, 1)
    x_coords = [self.r_nb[0],
                    self.r_nm[0],
                    self.r_t[0],
                    self.r_fm[0],
                    self.r_fb[0]]
    y_coords = [self.r_nb[1],
                    self.r_nm[1],
                    self.r_t[1],
                    self.r_fm[1],
                    self.r_fb[1]]

    ax.plot(x_coords, y_coords, color= dynein_color, label='dynein')
    ax.axis('off')
    ax.axis('equal')
    ax.legend()

dynein1 = bb_energy_distribution.DyneinBothBound(angles[0][0], angles[1][0], params, L)
dynein2 = bb_energy_distribution.DyneinBothBound(angles[0][1], angles[1][1], params, L)
dynein3 = bb_energy_distribution.DyneinBothBound(angles[0][2], angles[1][2], params, L)
dynein4 = bb_energy_distribution.DyneinBothBound(angles[0][3], angles[1][3], params, L)
dynein5 = bb_energy_distribution.DyneinBothBound(angles[0][4], angles[1][4], params, L)
dynein6 = bb_energy_distribution.DyneinBothBound(angles[0][5], angles[1][5], params, L)
dynein7 = bb_energy_distribution.DyneinBothBound(angles[0][6], angles[1][6], params, L)
dynein8 = bb_energy_distribution.DyneinBothBound(angles[0][7], angles[1][7], params, L)
dynein9 = bb_energy_distribution.DyneinBothBound(angles[0][8], angles[1][8], params, L)

plot_bb_figures(dynein1, 'red')
plot_bb_figures(dynein2, 'orange')
plot_bb_figures(dynein3, 'lime')
plot_bb_figures(dynein4, 'blue')
plot_bb_figures(dynein5, 'gray')
plot_bb_figures(dynein6, 'magenta')
plot_bb_figures(dynein7, 'green')
plot_bb_figures(dynein8, 'purple')
plot_bb_figures(dynein9, 'sienna')

num_points = 500
nma1 = np.linspace(0, 2*np.pi, num_points)
fma1 = np.linspace(0, 2*np.pi, num_points)
NMA, FMA = np.meshgrid(nma1, fma1)
dynein_24 = bb_energy_distribution.DyneinBothBound(NMA, FMA, params, L=24)
plot_bb_energy_distribution(dynein_24, dynein1, dynein2, dynein3, dynein4, dynein5, dynein6, dynein7, dynein8, dynein9)

plt.show()
