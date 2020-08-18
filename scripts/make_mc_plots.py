import os
from os import path
import numpy as np
from numpy.linalg import matrix_power
import matplotlib.pyplot as plt
from matplotlib import gridspec
from statistics import mean
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants
import sys
sys.path.append("../data")
import importlib
import argparse
import subprocess
import bb_energy_distribution
from glob import glob

def integrate_1d(arr, dx):
    sum = 0
    for i in range(len(arr)):
        sum += arr[i]*dx[i]
    return sum

def integrate_2d(arr, dx, dy):
    sum = 0
    for i in range(len(arr)):
        for j in range(len(arr[i])):
            sum += arr[i,j]*dx[j]*dy[i]
    return sum

def least_squares(arr, x, y, dx, dy):
    # X,Y = np.meshgrid(x,y)
    # x_mean = integrate_2d(arr*X, dx, dy)
    # y_mean = integrate_2d(arr*Y, dx, dy)
    # x2_mean = integrate_2d(arr*X**2, dx, dy)
    # xy_mean = integrate_2d(arr*X*Y, dx, dy)
    x_mean = 0
    y_mean = 0
    x2_mean = 0
    xy_mean = 0

    for i in range(len(arr)):
        for j in range(len(arr[i])):
            x_mean += arr[i,j]*x[j]*dx[i]*dy[j]
            y_mean += arr[i,j]*x[i]*dx[i]*dy[j]
            x2_mean += arr[i,j]*y[j]**2*dx[i]*dy[j]
            xy_mean += arr[i,j]*y[j]*x[i]*dx[i]*dy[j]

    # Intercept & Slope for Best-fit
    b = (y_mean*x2_mean-xy_mean*x_mean)/(x2_mean-x_mean**2)
    m = (-b*x_mean+xy_mean)/(x2_mean)
    lin_fit = (m*x)+b

    return b, m, lin_fit


def L_to_initial_displacement(P_leading, P_trailing):
    """
    Returns a matrix (technically array) which when multiplied by a
    probability density of L will result in a probability density of
    displacement.  Thus this is a dimensionless 2D array.

    """
    num_col = 2*len(P_leading)
    num_rows = len(P_leading)
    P_step = np.zeros((num_col,num_rows))
    for i in range(num_col):
        if i < num_col/2:
            P_step[num_rows-1-i,i] = P_trailing[num_rows-1-i]              # P_trailing array starts at 1 to 50
            P_step[num_rows+i, i] = P_leading[num_rows-1-i]        # P_leading array starts at 1 to 50
    return P_step

def L_to_L(T, P_leading, P_trailing):
    num_col = 2*len(P_leading)
    num_rows = len(P_leading)
    abs = np.zeros((num_rows,num_col))
    for i in range(num_col):
        if i < num_col/2:
            abs[num_rows-1-i,i] = 1
        else:
            abs[i-num_rows,i] = 1
    prob_step = L_to_initial_displacement(P_leading, P_trailing)
    T_L = abs*T*prob_step
    # plt.figure('abs')
    # plt.pcolor(abs);
    # plt.colorbar();
    # plt.figure('prob_step')
    # plt.pcolor(prob_step);
    # plt.colorbar();
    # plt.show()
    return T_L

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-k", "--kb", type=float, help="Binding const", default=params.for_simulation['k_b'])
parser.add_argument("-s", "--ks", type=float, help="Sticky const", default=params.for_simulation['k_stk'])
args = parser.parse_args()

k_b = float(args.kb)        # Binding Rate Constant
k_stk = float(args.ks)      # Sticky Rate Constant

plotpath = '../plots/mc_plots/'
if not os.path.exists(plotpath):
    os.mkdir(plotpath)
plottingdatapath = '../data/mc_plotting_data/mc_plotting_data_{0:.2e}_{1:.2e}.npz'.format(k_b, k_stk)

mc_data = np.load(plottingdatapath, allow_pickle=True)
initial_disp = mc_data['initial_disp']
hist = mc_data['hist']
normalized_hist = mc_data['normalized_hist']
P_unbinding = mc_data['P_unbinding'].item()

# make bin center the data point (for pcolor)
initial_disp = np.array(sorted(initial_disp)) # -50 to 50 array
final_bin_center = initial_disp*1.0 # set final_bin_center to initial_disp
final_bin_edges = np.zeros(len(final_bin_center)+1)
for i in range(1,len(final_bin_edges)-1):
    final_bin_edges[i] = (final_bin_center[i-1] + final_bin_center[i])*0.5
final_bin_edges[0] = 2*final_bin_center[0] - final_bin_edges[1]
final_bin_edges[-1] = 2*final_bin_center[-1] - final_bin_edges[-2]

# obtain meshgrid for pcolor
initial_disp_center, final_disp_center = np.meshgrid(initial_disp, final_bin_center)

final_disp_bin_width = final_bin_edges[1:] - final_bin_edges[:-1]    # a 1D array giving final displacement bin width

initial_disp_edge, final_disp_edge = np.meshgrid(final_bin_edges, final_bin_edges)


plt.close('all')
plt.figure('From Data')
plt.pcolor(initial_disp_edge, final_disp_edge, normalized_hist)
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()
plt.savefig(plotpath+'2dhist_initL_vs_finalL_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))

# Transition Matrix
T = np.matrix(hist)     # Dimensionless

# Probability Column Vector for L
P = np.matrix(np.zeros((int(len(T)/2),1)))      # Dimensionless column vector

# Get bin widths for final L
final_L = np.array(final_bin_center[len(final_bin_center)//2:])
final_L_bin_width = np.zeros_like(final_L)
final_L_bin_width[1:-1] = (final_L[2:] - final_L[:-2])/2
final_L_bin_width[0] = final_L[0] + (final_L[1] - final_L[0])/2
final_L_bin_width[-1] = final_L[-1] - final_L[-2]

# Get bin widths for initial L
initial_L = final_L
initial_L_bin_width = final_L_bin_width

# Number of steps for Equilibrium
num_steps = 16

# Plot L to L probability density
plt.figure('prob density')
plt.legend(loc='best')
plt.xlabel('L')
plt.ylabel('probability per L')

# Obtain L to L probability density
for i in range(len(P)):
    # this loop is just to compute the L probability density many
    # different ways so we can tell that num_steps is sufficiently
    # high.
    P[:,:]= 0
    P[i] = 1

    P_L_to_L = (L_to_L(T, P_unbinding['leading'], P_unbinding['trailing'])**num_steps)*P      # Dimensionless 2D Array that only has 1 column vector
    P_L_to_L = np.array(P_L_to_L)[:,0]     # convert to a dimensionless 1D array from a column vector
    # print('P_L_to_L.sum()', P_L_to_L.sum())
    norm_const = 1/((P_L_to_L*final_L_bin_width).sum())     # Dimensions: 1/distance, sum of (P_L_to_L flat * bin width of both axis)
    # print('norm_const', norm_const)
    p_den_L = P_L_to_L*norm_const     # dimensions 1/distance, a probability density
    # print('p_den_L SUM:', (p_den_L*final_L_bin_width).sum())
    plt.plot(final_L, p_den_L, label=f'i is {i}')

plt.savefig(plotpath+'L_to_L_prob_density_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))

# print(p_den_L)
p_den_disp = L_to_initial_displacement(P_unbinding['leading'], P_unbinding['trailing']).dot(p_den_L)   # Dimensions: 1/distance
# print(p_den_disp.shape)
plt.figure('p_den_disp')
plt.plot(initial_disp, p_den_disp)
plt.xlabel('displacement')
plt.ylabel('probability density')
plt.savefig(plotpath+'Probability_density_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))


# Probability Distribution is the normalized histogram multiplied by the probability density
probability_distribution = np.zeros_like(normalized_hist)   # Dimensions: 1/distance**2
for i in range(probability_distribution.shape[0]):
    probability_distribution[i,:] = normalized_hist[i,:]*p_den_disp

filtered_probability_distribution = probability_distribution*1.0    # Dimensions: 1/distance**2

print('prob_den:', np.sum(p_den_disp*final_disp_bin_width))

# Sum of probability distribution
probability_distirbution_sum = integrate_2d(probability_distribution, final_disp_bin_width, final_disp_bin_width)


# Filter steps where final and initial displacements are within 4 nm of each other
# FIXME use actual distances here!!!
for i in range(len(probability_distribution)):
    for j in range(len(probability_distribution[i])):
        if i == j:
            filtered_probability_distribution[i, j] = 0.0
            if i >=4 and i <= 95:
                filtered_probability_distribution[i-4:i+4,j] = 0.0

# Sum of filtered probability distribution
probability_distribution_sum_filt = integrate_2d(filtered_probability_distribution, final_disp_bin_width, final_disp_bin_width)

# Normalize the filtered historgram
filtered_probability_distribution /= probability_distribution_sum_filt

# Intercept & Slope for Best-fit
b, m, lin_fit = least_squares(probability_distribution, initial_disp, initial_disp, final_disp_bin_width, final_disp_bin_width)


# Intercept & Slope for Best-fit Filtered
b_filt, m_filt, lin_fit_filt = least_squares(filtered_probability_distribution, initial_disp, initial_disp, final_disp_bin_width, final_disp_bin_width)



plt.figure('Probability Distribution to Match Yildiz')
plt.pcolor(initial_disp_edge, final_disp_edge, probability_distribution)
plt.plot(initial_disp, lin_fit, label='Model: y = ({:.3}) + ({:.3})x'.format(b,m), linestyle=":", color='r')
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()
plt.legend()
plt.savefig(plotpath+'Match_Yildiz_probability_distribution_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))



plt.figure('Filtered Probability Distribution to Match Yildiz')
plt.pcolor(initial_disp_edge, final_disp_edge, filtered_probability_distribution)
plt.plot(initial_disp, lin_fit_filt, label='Model: y = ({:.3}) + ({:.3})x'.format(b_filt,m_filt), linestyle=":", color='r')
plt.xlabel('initial displacement (nm)')
plt.ylabel('final displacement (nm)')
plt.colorbar()
plt.legend()
plt.savefig(plotpath+'filtered_Match_Yildiz_probability_distribution_{0:.2e}_{1:.2e}.pdf'.format(float(k_b), float(k_stk)))

print('FINAL SUM: ', integrate_2d(probability_distribution, final_disp_bin_width, final_disp_bin_width))


print("""
TO DO ITEMS:

a) Clean code and make less bug-prone:

  - DONE: Move to consistent names (displacement = disp, L only means L,
    etc, maybe P is dimensionless prob and p is probability density,
    bin width is delta_disp, delta_L?)

  - Possibly all caps for 2D arrays? Or some other naming.

  - DONE: Perhaps introduce functions to integrate? 1D & 2D

  - Possibly introduce "random" units.  Define nm = random #, and then
    multiply by nm when you read, divide by nm if you want to print a
    distance in nm.  (low priority)

  - Mainly DONE: Document dimensions and meaning of quantities.

  - Introduce 2D arrays for 1D quantities.

  - Document each functions/more

  - Think of using meshgrid for almost everything.

  - HIGH PRIORITY:  Fix filtering prob dist

  - Add Yildiz fit to the match plot.
  """)

plt.show()
