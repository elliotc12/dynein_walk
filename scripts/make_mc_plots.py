from glob import glob
import bb_energy_distribution
import subprocess
import argparse
import importlib
import os
from os import path, mkdir
import numpy as np
from numpy.linalg import matrix_power
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import gridspec
from statistics import mean
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
import scipy.constants
import sys
sys.path.append("../data")


def integrate_1d(arr, dx):
    sum = 0
    for i in range(len(arr)):
        sum += arr[i]*dx[i]
    return sum


def integrate_2d(arr, dx, dy):
    sum = 0
    for i in range(len(arr)):
        for j in range(len(arr[i])):
            sum += arr[i, j]*dx[j]*dy[i]
    return sum


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


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

    for j in range(len(arr)):
        for i in range(len(arr[j])):
            x_mean += arr[j, i]*x[i]*dx[i]*dy[j]
            y_mean += arr[j, i]*y[j]*dx[i]*dy[j]
            x2_mean += arr[j, i]*x[i]**2*dx[i]*dy[j]
            xy_mean += arr[j, i]*y[j]*x[i]*dx[i]*dy[j]

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
    P_step = np.zeros((num_col, num_rows))
    for i in range(num_col):
        if i < num_col/2:
            # P_step[num_rows-1-i,i] = P_trailing[num_rows-1-i]              # P_trailing array starts at 1 to 50
            # P_step[num_rows+i, i] = P_leading[num_rows-1-i]        # P_leading array starts at 1 to 50
            # P_trailing array starts at 1 to 50
            P_step[num_rows-1-i, i] = P_trailing[i]
            # P_leading array starts at 1 to 50
            P_step[num_rows+i, i] = P_leading[i]
    return P_step


def L_to_L(T, P_leading, P_trailing):
    num_col = 2*len(P_leading)
    num_rows = len(P_leading)
    abs = np.zeros((num_rows, num_col))
    for i in range(num_col):
        if i < num_col/2:
            abs[num_rows-1-i, i] = 1
        else:
            abs[i-num_rows, i] = 1
    prob_step = L_to_initial_displacement(P_leading, P_trailing)
    T_L = abs*T*prob_step
    # print(abs)
    # print(np.shape(abs))
    # print(prob_step)
    # print(np.shape(prob_step))
    # plt.figure('abs')
    # plt.pcolor(abs);
    # plt.colorbar();
    # plt.figure('prob_step')
    # plt.pcolor(prob_step);
    # plt.colorbar();
    # plt.show()
    # exit()
    return T_L

# FIXME Add exp unbinding const C to naming scheme when saving


def get_cli_arguments():
    parser = argparse.ArgumentParser(
        description='Script to generate various plots from Monte Carlo data.')
    parser.add_argument("-u", "--k_ub", type=float,
                        help="Unbinding rate", default=params.for_simulation['k_ub'])
    parser.add_argument("-b", "--k_b", type=float,
                        help="Binding rate", default=params.for_simulation['k_b'])
    parser.add_argument("-s", "--k_stk", type=float,
                        help="Sticky rate", default=params.for_simulation['k_stk'])
    parser.add_argument("-cb", "--cb", type=float,
                        help="Spring const binding domain", default=params.for_simulation['cb'])
    parser.add_argument("-cm", "--cm", type=float,
                        help="Spring const motor domain", default=params.for_simulation['cm'])
    parser.add_argument("-ct", "--ct", type=float,
                        help="Spring const tail domain", default=params.for_simulation['ct'])
    parser.add_argument("--eqb", type=float, help="Binding Equilibrium angle",
                        default=params.for_simulation['eqb'])
    parser.add_argument("--eqmpre", type=float, help="Motor pre Equilibrium angle",
                        default=params.for_simulation['eqmpre'])
    parser.add_argument("--eqmpost", type=float, help="Motor post Equilibrium angle",
                        default=params.for_simulation['eqmpost'])
    parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant",
                        default=params.for_simulation['exp-unbinding-constant'])
    parser.add_argument("-p", "--plot", action="store_false",
                        help="Do not show plots", default=True)
    parser.add_argument("-e", "--extra", action="store_false",
                        help="Show extra plots", default=True)
    parser.add_argument("--underMT", action="store_false",
                        help="Plot sims where binding domain can go under MT", default=True)
    return parser.parse_args()


def get_onebound_data(args, plotting_data_file):
    # Onebound Data
    mc_data = np.load(plotting_data_file, allow_pickle=True)
    return mc_data['initial_disp'], mc_data['hist'], mc_data['normalized_hist'], mc_data['time_hists'].item(), mc_data['avg_affinity_time'].item()


def get_bothbound_data(args, bothbound_data_file):
    # Bothbound Data
    mc_bb_data = np.load(bothbound_data_file, allow_pickle=True)
    bb_L = mc_bb_data['L']
    bb_rate_leading = mc_bb_data['rate_leading']*args.k_ub
    # bb_rate_leading = mc_bb_data['rate_leading']*params.for_simulation['k_ub']
    bb_rate_trailing = mc_bb_data['rate_trailing']*args.k_ub
    # bb_rate_trailing = mc_bb_data['rate_trailing']*params.for_simulation['k_ub']
    bb_P_leading = bb_rate_leading/(bb_rate_leading+bb_rate_trailing)
    bb_P_trailing = bb_rate_trailing/(bb_rate_leading+bb_rate_trailing)
    bb_lt = 1/bb_rate_leading
    bb_tt = 1/bb_rate_trailing
    bb_avg_t = 1/(bb_rate_leading+bb_rate_trailing)

    return bb_L, bb_P_leading, bb_P_trailing, bb_rate_trailing+bb_rate_leading, bb_avg_t, bb_lt, bb_tt


def make_bins_and_edges(initial_disp):
    # make bin center the data point (for pcolor)
    initial_disp = np.array(sorted(initial_disp))  # -50 to 50 array
    final_bin_center = initial_disp*1.0  # set final_bin_center to initial_disp
    final_bin_edges = np.zeros(len(final_bin_center)+1)
    for i in range(1, len(final_bin_edges)-1):
        final_bin_edges[i] = (final_bin_center[i-1] + final_bin_center[i])*0.5
    final_bin_edges[0] = 2*final_bin_center[0] - final_bin_edges[1]
    final_bin_edges[-1] = 2*final_bin_center[-1] - final_bin_edges[-2]

    # obtain meshgrid for pcolor
    initial_disp_center, final_disp_center = np.meshgrid(
        initial_disp, final_bin_center)
    # a 1D array giving final displacement bin width
    final_disp_bin_width = final_bin_edges[1:] - final_bin_edges[:-1]
    initial_disp_edge, final_disp_edge = np.meshgrid(
        final_bin_edges, final_bin_edges)

    # Get bin widths for final L
    final_L = np.array(final_bin_center[len(final_bin_center)//2:])
    final_L_bin_width = np.zeros_like(final_L)
    final_L_bin_width[1:-1] = (final_L[2:] - final_L[:-2])/2
    final_L_bin_width[0] = final_L[0] + (final_L[1] - final_L[0])/2
    final_L_bin_width[-1] = final_L[-1] - final_L[-2]

    # Get bin widths for initial L
    initial_L = final_L
    initial_L_bin_width = final_L_bin_width
    return initial_disp_edge, final_disp_edge, final_L_bin_width, final_disp_bin_width, initial_L


def make_probability_distribution(args, hist, normalized_hist, bb_P_leading, bb_P_trailing, initial_L, final_L_bin_width, initial_disp, **_):
    # Transition Matrix
    T = np.matrix(hist)     # Dimensionless

    # Probability Column Vector for L
    # Dimensionless column vector
    P = np.matrix(np.zeros((int(len(T)/2), 1)))

    # Unbinding probability based on bb simulation
    P_ub_leading = bb_P_leading[initial_L.astype(int)-1]
    P_ub_trailing = bb_P_trailing[initial_L.astype(int)-1]

    # Number of steps for Equilibrium
    num_steps = 16

    # Obtain L to L probability density
    for i in range(len(P)):
        # this loop is just to compute the L probability density many
        # different ways so we can tell that num_steps is sufficiently
        # high.
        P[:, :] = 0
        P[i] = 1

        # Dimensionless 2D Array that only has 1 column vector
        P_L_to_L = (L_to_L(T, P_ub_leading, P_ub_trailing)**num_steps)*P
        # convert to a dimensionless 1D array from a column vector
        P_L_to_L = np.array(P_L_to_L)[:, 0]
        # First normalize the dimensionless probability
        P_L_to_L /= P_L_to_L.sum()
        # dimensions 1/distance, a probability density, divide by size of each bin
        p_den_L = P_L_to_L/final_L_bin_width

    # plt.figure('p_den_L')
    # plt.plot(initial_L, p_den_L)
    # plt.xlabel('L (nm)')
    # plt.ylabel('Probability density')
    # plt.title('kb = {0:.2e}, kstk = {1:.2e}, cb = {2}, cm = {3}, ct = {4}, eqb = {5}, eqmpre = {6}, eqmpost = {7}, C = {8}'.format(args.k_b,
        # args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C), fontsize=7)

    p_den_disp = L_to_initial_displacement(
        P_ub_leading, P_ub_trailing).dot(p_den_L)   # Dimensions: 1/distance
    # plt.figure('p_den_disp')
    # plt.plot(initial_disp, p_den_disp)
    # plt.xlabel('L (nm)')
    # plt.ylabel('Probability density')
    # plt.title('kb = {0:.2e}, kstk = {1:.2e}, cb = {2}, cm = {3}, ct = {4}, eqb = {5}, eqmpre = {6}, eqmpost = {7}, C = {8}'.format(args.k_b,
    # args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C), fontsize=7)

    # Probability Distribution is the normalized histogram multiplied by the probability density
    probability_distribution = np.zeros_like(
        normalized_hist)   # Dimensions: 1/distance**2
    for i in range(probability_distribution.shape[0]):
        probability_distribution[i, :] = normalized_hist[i, :]*p_den_disp

    return probability_distribution


def make_prob_dist_plot(args, plotpath, probability_distribution, initial_disp_edge, final_disp_edge, initial_disp, b, m, lin_fit, **_):
    # Yildiz final disp vs initial disp line
    yildiz_line = [(0.6*x)+9.1 for x in np.asarray(initial_disp)]

    plt.figure('Probability Distribution to Match Yildiz', figsize=(8, 6))
    plt.pcolor(initial_disp_edge, final_disp_edge,
               probability_distribution)
    plt.plot(initial_disp, lin_fit,
             label='Model: y = ({:.2}) + ({:.1})x'.format(b, m), linestyle=":", color='C1')
    plt.plot(initial_disp, yildiz_line,
             label='Experiment: y = (9.1) + (0.6)x', linestyle=":", color='Red')
    plt.xlabel('Initial displacement (nm)')
    plt.ylabel('Final displacement (nm)')
    plt.xlim(-48, 48)
    plt.ylim(-48, 48)
    plt.colorbar().set_label(r'Probability (1/nm$^2$)')
    plt.legend()
    plt.subplots_adjust(top=0.90, bottom=0.10)
    # plt.title('kb = {0:.2e}, kstk = {1:.2e}, cb = {2}, cm = {3}, ct = {4}, eqb = {5}, eqmpre = {6}, eqmpost = {7}, C = {8}'.format(args.k_b,
    # args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C), fontsize=7)
    plt.savefig(plotpath+u+'final_disp_probability_distribution_' +
                params_string+'.pdf')


def make_filtered_prob_dist_plot(args, plotpath, probability_distribution, initial_disp_edge, final_disp_edge, initial_disp, final_disp_bin_width, **_):
    filtered_probability_distribution = probability_distribution * \
        1.0    # Dimensions: 1/distance**2

    # Filter steps where final and initial displacements are within 4 nm of each other
    init_L = initial_disp[int(len(initial_disp)/2):]
    min_init_L = int(min(init_L))
    index_subtract = 1
    while min_init_L < 2:
        index_subtract += 1
        init_L = init_L[1:]*1.0
        min_init_L = int(min(init_L))
    for i in range(len(probability_distribution)):
        for j in range(len(probability_distribution[i])):
            if i == j:
                filtered_probability_distribution[i, j] = 0.0
                if i >= index_subtract and i <= len(probability_distribution)-index_subtract-1:
                    filtered_probability_distribution[i -
                                                      index_subtract:i+index_subtract, j] = 0.0

    # Sum of filtered probability distribution
    filt_probability_distribution_sum = integrate_2d(
        filtered_probability_distribution, final_disp_bin_width, final_disp_bin_width)

    # Normalize the filtered historgram
    filtered_probability_distribution /= filt_probability_distribution_sum

    # Intercept & Slope for Best-fit Filtered
    b_filt, m_filt, lin_fit_filt = least_squares(
        filtered_probability_distribution, initial_disp, initial_disp, final_disp_bin_width, final_disp_bin_width)

    plt.figure('Filtered Probability Distribution to Match Yildiz')
    plt.pcolor(initial_disp_edge, final_disp_edge,
               filtered_probability_distribution)
    plt.plot(initial_disp, lin_fit_filt,
             label='Model: y = ({:.3}) + ({:.3})x'.format(b_filt, m_filt), linestyle=":", color='r')
    plt.xlabel('Initial displacement (nm)')
    plt.ylabel('Final displacement (nm)')
    plt.colorbar()
    plt.legend()
    # plt.title('kb = {0:.2e}, kstk = {1:.2e}, cb = {2}, cm = {3}, ct = {4}, eqb = {5}, eqmpre = {6}, eqmpost = {7}, C = {8}'.format(args.k_b,
    # args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C))
    plt.savefig(
        plotpath+u+'filtered_final_disp_probability_distribution_'+params_string+'.pdf')


def make_step_length_plots(args, plotpath, probability_distribution, initial_disp_edge, final_disp_edge, initial_disp, final_disp_bin_width, **_):
    # STEP LENGTH (s) PLOTS CALCULATION
    num_disp = len(probability_distribution)
    initial_disp_index = np.arange(0, num_disp)
    # step length probability density axis for 1d hist
    s_den = np.zeros((2*num_disp-1))
    # step length axis (or step length bin center)
    s_arr = np.zeros_like(s_den)
    ds = np.zeros(2*num_disp-1)
    ds[:num_disp] = final_disp_bin_width
    ds[num_disp:] = final_disp_bin_width[1:]
    s_arr[0] = -np.sum(final_disp_bin_width)+ds[0]  # set step length values
    s_arr[-1] = np.sum(final_disp_bin_width)-ds[-1]
    for i in range(1, len(ds)-1):
        s_arr[i] = s_arr[i-1] + ds[i]

    # 2D probability distribution of step length
    s_bin_edges = np.zeros(len(s_arr)+1)
    for i in range(1, len(s_bin_edges)-1):
        s_bin_edges[i] = (s_arr[i-1] + s_arr[i])*0.5
    s_bin_edges[0] = 2*s_arr[0] - s_bin_edges[1]
    s_bin_edges[-1] = 2*s_arr[-1] - s_bin_edges[-2]

    # step length histogram from Yildiz 2012 figure 1B top head-labelled panel
    yildiz_hist_fig_1B = np.array(
        [[-40, 0], [-39, 0], [-38, 0], [-37, 1], [-36, 0], [-35, 1], [-34, 1], [-33, 2], [-31, 2],
         [-30, 3], [-29, 1], [-28, 1], [-27, 4], [-26, 4], [-25,
                                                            2], [-24, 3], [-23, 4], [-22, 0], [-21, 4], [-20, 3],
         [-19, 3], [-18, 5], [-17, 3], [-16, 3], [-15,
                                                  7], [-14, 5], [-13, 7], [-12, 12], [-11, 16],
         [-10, 14], [-9, 20], [-8, 14], [-7, 10], [-6,
                                                   9], [-5, 11], [-4, 8], [-3, 2], [-2, 0], [-1, 0],
         [0, 0], [1, 0], [2, 0], [3, 0], [4, 6],
         [5, 7], [6, 12], [7, 20], [8, 19], [9, 22], [10, 30], [11, 34], [
             12, 26], [13, 21], [14, 23], [15, 22], [16, 30], [17, 29],
         [18, 23], [19, 22], [20, 26], [21, 12], [22, 21], [23, 16], [24, 7], [25, 8], [
             26, 7], [27, 8], [28, 5], [29, 9], [30, 7], [31, 8], [32, 6], [33, 2],
         [34, 2], [35, 9], [36, 4], [37, 9], [38, 5], [39, 1], [40, 2], [41, 1], [42, 3], [
             43, 2], [44, 4], [45, 4], [46, 1], [47, 1], [48, 0], [49, 0], [50, 0], [51, 0],
         [52, 0], [53, 0], [54, 0], [55, 0]])
    yildiz_step_lengths_fig_1B = []
    for i in range(len(yildiz_hist_fig_1B)):
        yildiz_step_lengths_fig_1B.extend(
            [yildiz_hist_fig_1B[i, 0]]*yildiz_hist_fig_1B[i, 1])
    yildiz_step_lengths_fig_1B = np.array(yildiz_step_lengths_fig_1B)
    step_length_bin_center_fig_1B = 1.0*yildiz_hist_fig_1B[:, 0]
    step_length_count_fig_1B = 1.0*yildiz_hist_fig_1B[:, 1]

    norm_length_min = 7.75
    if args.k_stk == 9e99:
        norm_length_min = 1

    norm_1B = np.sum(
        step_length_count_fig_1B[step_length_bin_center_fig_1B > norm_length_min])
    yildiz_normalized_prob_fig_1B = step_length_count_fig_1B/norm_1B

    step_length_bin_edges_fig_1B = np.zeros(
        2*len(step_length_bin_center_fig_1B))
    yildiz_normalized_prob_edges_fig_1B = np.zeros(
        2*len(step_length_bin_center_fig_1B))
    for i in range(len(step_length_bin_center_fig_1B)):
        step_length_bin_edges_fig_1B[2 *
                                     i] = step_length_bin_center_fig_1B[i]-0.5
        step_length_bin_edges_fig_1B[2*i +
                                     1] = step_length_bin_center_fig_1B[i]+0.5
        yildiz_normalized_prob_edges_fig_1B[2 *
                                            i] = yildiz_normalized_prob_fig_1B[i]
        yildiz_normalized_prob_edges_fig_1B[2 *
                                            i+1] = yildiz_normalized_prob_fig_1B[i]

    # load data from Yildiz 2012 Fig 3.A scatter plot
    yildiz_IF_data = np.loadtxt(
        "../data/yildiz_2012_if_scatter_coordinates.txt", delimiter=", ")

    center_px_x = 157  # from data file
    center_px_y = 128
    right_center_px_x = 311
    right_center_nm_x = 56
    bottom_center_px_y = 254
    bottom_center_nm_y = 56

    yildiz_data_nm = yildiz_IF_data + \
        np.array((-center_px_x, -center_px_y))  # shift origin to (0, 0)
    yildiz_data_nm = yildiz_data_nm*np.array((1.0, -1.0))  # invert the y
    yildiz_data_nm = yildiz_data_nm*np.array((right_center_nm_x / (right_center_px_x-center_px_x),
                                              bottom_center_nm_y / (bottom_center_px_y-center_px_y)))  # rescale to nm

    yildiz_step_lengths_fig_3A = yildiz_data_nm[:, 1]
    step_length_count_fig_3A = np.zeros_like(step_length_count_fig_1B)

    for i in range(len(yildiz_step_lengths_fig_3A)):
        bin = int(
            np.round(yildiz_step_lengths_fig_3A[i]) - step_length_bin_center_fig_1B[0])
        if bin == len(step_length_count_fig_3A):
            bin = bin-1
        step_length_count_fig_3A[bin] += 1
    norm_3A = np.sum(
        step_length_count_fig_3A[step_length_bin_center_fig_1B > norm_length_min])
    yildiz_normalized_prob_fig_3A = step_length_count_fig_3A/norm_3A

    step_length_bin_edges_fig_3A = np.zeros_like(step_length_bin_edges_fig_1B)
    yildiz_normalized_prob_edges_fig_3A = np.zeros_like(
        yildiz_normalized_prob_edges_fig_1B)
    for i in range(len(step_length_bin_center_fig_1B)):
        step_length_bin_edges_fig_3A[2 *
                                     i] = step_length_bin_center_fig_1B[i]-0.5
        step_length_bin_edges_fig_3A[2*i +
                                     1] = step_length_bin_center_fig_1B[i]+0.5
        yildiz_normalized_prob_edges_fig_3A[2 *
                                            i] = yildiz_normalized_prob_fig_3A[i]
        yildiz_normalized_prob_edges_fig_3A[2 *
                                            i+1] = yildiz_normalized_prob_fig_3A[i]


    # FOR 2D Plot
    s_initial_disp_edge, s_length_edge = np.meshgrid(
        initial_disp_edge[0], s_bin_edges)
    # a 1D array giving final displacement bin width
    s_bin_width = s_bin_edges[1:] - s_bin_edges[:-1]

    s_probability_distribution = np.zeros((2*num_disp-1, num_disp))
    initial_disp_center = 0.5*(initial_disp_edge[0][:-1]  + initial_disp_edge[0][1:])
    step_length_center = 0.5*(s_bin_edges[1:]+ s_bin_edges[:-1])
    for i_initial in initial_disp_index:
        for i_step_length in range(2*num_disp-1):
            this_step_length = step_length_center[i_step_length]
            final_displacement = this_step_length + initial_disp_center[i_initial]
            i_final = np.argmin(abs(initial_disp_center - final_displacement))
            s_probability_distribution[i_step_length, i_initial] = probability_distribution[i_final, i_initial]

    # Intercept & Slope for Best-fit of step length
    s_b, s_m, s_lin_fit = least_squares(
        s_probability_distribution, initial_disp, s_arr, final_disp_bin_width, s_bin_width)

    print('step length prob distribution 2d integral: ', integrate_2d(
        s_probability_distribution, final_disp_bin_width, s_bin_width))

    # Yildiz step length vs initial disp line
    yildiz_line = [(-0.4*x)+9.1 for x in np.asarray(initial_disp)]

    spd_max = np.max(s_probability_distribution)

    def percentile(prob):
        return s_probability_distribution[s_probability_distribution>prob].sum()/s_probability_distribution.sum()
    def find_percentile(fraction):
        pmax = 1
        pmin = 0
        while pmax - pmin > 1e-6:
            pmid = 0.5*(pmin + pmax)
            if percentile(pmid) > fraction:
                pmin = pmid
            else:
                pmax = pmid
        return pmid
    quantiles = {
        # find_percentile(0.01): '1%',
        find_percentile(0.05): '5%',
        find_percentile(0.25): '25%',
        find_percentile(0.5): '50%',
        find_percentile(0.9): '90%',
    }

    s_initial_disp_center, s_length_center = np.meshgrid(
    initial_disp_center,
    step_length_center)

    # 2D hist step_length
    plt.figure('Step length probability distribution', figsize=(8, 6))
    contours = plt.contour(s_initial_disp_center, s_length_center,
                           s_probability_distribution, levels=sorted(list(quantiles.keys())), colors='white')
    plt.clabel(contours, quantiles.keys(), fmt=quantiles, inline=True, fontsize=7)
    plt.pcolor(s_initial_disp_edge, s_length_edge,
               s_probability_distribution)
    # plt.pcolor(initial_disp_edge, step_length_edge, s_probability_distribution)
    plt.plot(initial_disp, s_lin_fit,
             label='Model: y = ({:.3}) + ({:.3})x'.format(s_b, s_m), linestyle=":", color='C1')
    plt.plot(initial_disp, yildiz_line,
             label='Experiment: y = (9.1) + (-0.4)x', linestyle=":", color='Red')
    plt.xlabel('Initial displacement (nm)')
    plt.ylabel('step length (nm)')
    plt.gca().set_aspect('equal')
    plt.ylim(-56, 56)
    plt.xlim(-56, 56)
    plt.colorbar().set_label(r'Probability (1/nm$^2$)')
    plt.legend(loc='lower left')
    # plt.title('kb = {0:.2e}, kstk = {1:.2e}, cb = {2}, cm = {3}, ct = {4}, eqb = {5}, eqmpre = {6}, eqmpost = {7}, C = {8}'.format(args.k_b,
    #         args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C))
    plt.savefig(plotpath+u+'step_length_probability_distribution_' +
                params_string+'.pdf')


    s_den = np.zeros_like(step_length_center)
    for i in range(len(step_length_center)):
        s_den[i] = integrate_1d(s_probability_distribution[i], final_disp_bin_width)
    s_den_norm = np.sum(s_den[step_length_center > norm_length_min])*(step_length_center[1]-step_length_center[0])

    # 1D hist step length
    plt.figure('Probability Density of Step Length', figsize=(9, 5))
    plt.fill_between(step_length_center, 0*s_den, s_den/s_den_norm,
                     label='Model', color='C1')
    plt.plot(step_length_bin_edges_fig_3A,
             yildiz_normalized_prob_edges_fig_3A, color='C0', alpha=0.8)
    plt.fill_between(step_length_bin_edges_fig_3A, 0*yildiz_normalized_prob_edges_fig_3A,
                     yildiz_normalized_prob_edges_fig_3A, label='Experiment', color='C2', alpha=0.5)
    # plt.plot(step_length_bin_edges_fig_1B, yildiz_normalized_prob_edges_fig_1B, color='C0', alpha = 0.8)
    # plt.fill_between(step_length_bin_edges_fig_1B, 0*yildiz_normalized_prob_edges_fig_1B,
    #             yildiz_normalized_prob_edges_fig_1B, label='Experiment Fig 1B', color='C2', alpha=0.5)
    plt.xlabel('Step Length (nm)')
    plt.ylabel('Probability Density (1/nm)')
    plt.xlim(-50, 65)
    if args.k_stk == 9e99:
        plt.ylim(0, 0.08)
    else:
        plt.ylim(0)
    plt.legend()
    plt.subplots_adjust(top=0.95)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    # plt.title('kb = {0:.2e}, kstk = {1:.2e}, cb = {2}, cm = {3}, ct = {4}, eqb = {5}, eqmpre = {6}, eqmpost = {7}, C = {8}'.format(args.k_b,
    # args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C), fontsize=7)
    plt.savefig(plotpath+u+'step_length_1d_probability_density_' +
                params_string+'.pdf')

def make_ob_time_plot(args, plotpath, time_hists, avg_affinity_time, **_):
    max_time = time_hists['max_time']   # 10 mu s
    increment = time_hists['increment']  # 10ns
    time_bin_center = np.arange(increment/2, max_time, increment, dtype=float)
    x_upper_lim = 0.5e-6
    x_increment = 0.1e-6
    decimal = 1
    if args.k_stk == 9e99:
        x_upper_lim = 0.05e-6
        x_increment = 0.01e-6
        decimal = 2
    ind = 1
    for j in time_bin_center:
        if j < x_upper_lim:
            ind += 1
        else:
            break

    iL = [num for num in time_hists.keys() if isinstance(num, (int, float))]

    # OB Time Plot
    for i in iL:
        if len(time_hists[i]) == 0:
            print("Do not have time data for initial distance ", i)
        else:
            if i == 16.0:

                print('Remaining probability outside of x range for ',
                      i, ', ', np.sum(time_hists[i][ind:]))
                print('Remaining probability outside of x range for -',
                      i, ', ', np.sum(time_hists[-i][ind:]))
                print('Remaining total prob outside of x range, ', (np.sum(
                    time_hists[i][ind:])+np.sum(time_hists[-i][ind:]))/2)

                plt.figure('Onebound Time plot for iL = {}'.format(
                    i), figsize=(7, 5))
                # plt.fill_between(time_bin_center,0*time_hists[8], time_hists[8], label='Model')
                # print(time_bin_center)
                plt.bar(time_bin_center, time_hists[i], width=increment,
                        align='center', label='Leading', alpha=0.8, color='xkcd:blue')
                plt.bar(time_bin_center, time_hists[-i], width=increment,
                        align='center', label='Trailing', alpha=0.6, color='xkcd:green')
                plt.xlabel(r'One-bound time ($\mu s$)')
                plt.ylabel('Probability')
                plt.axvline(avg_affinity_time[i], color='Red', ls='--',
                            linewidth=2, label='Average affinity time')
                # plt.axvline(avg_affinity_time[-i], color='Red', ls='--')
                xticks = np.arange(0.0, x_upper_lim+x_increment, x_increment)
                labels = np.arange(
                    0.0, (x_upper_lim+x_increment)*1e6, x_increment*1e6)
                plt.xlim(0, x_upper_lim)
                plt.xticks(xticks, np.around(labels, decimals=decimal))
                plt.legend()
                plt.subplots_adjust(top=0.93)
                plt.gca().spines['right'].set_visible(False)
                plt.gca().spines['top'].set_visible(False)
                # plt.title('kb = {0:.2e}, kstk = {1:.2e}, cb = {2}, cm = {3}, ct = {4}, eqb = {5}, eqmpre = {6}, eqmpost = {7}, C = {8}'.format(args.k_b,
                # args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C))
                plt.savefig(
                    plotpath+u+'ob_time_probability_density_'+params_string+'.pdf')


def make_bothbound_plots(args, plotpath, bb_L, bb_P_trailing, bb_avg_t, bb_lt, bb_tt, **_):
    # Bothbound PLOTS
    # Prob Lagging vs Initial L plot
    plt.figure('Prob lagging vs init L', figsize=(9, 4))

    yildiz_displacements = [10, 20, 30, 40, 50]
    yildiz_lagging_fractions = [0.525, 0.545, 0.61, 0.59, 0.67]
    yildiz_lagging_uncertainty = [0.06, 0.04, 0.035, 0.045, 0.075]
    plt.errorbar(yildiz_displacements, yildiz_lagging_fractions, yerr=yildiz_lagging_uncertainty,
                 label="Experiment", fmt='o-', c='C0', linestyle='', capsize=3)

    plt.plot(bb_L, bb_P_trailing, label='Model', color='C1')
    plt.xlim(0, 55)
    plt.xlabel('Binding domain separation (nm)')
    plt.ylabel('P(lagging step)')
    plt.legend()
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.subplots_adjust(top=0.92)
    # plt.title('C = {}'.format(args.C))
    plt.savefig(
        plotpath+'prob_lagging_vs_init_L_{}.pdf'.format(args.C), bbox_inches='tight')

    # Bothbound time plot
    plt.figure('BB time plot')
    plt.plot(bb_L, bb_avg_t, color='C0')
    plt.xlabel('initial L (nm)')
    plt.ylabel('Average time (s)')
    plt.xlim(0)
    plt.ylim(0, 0.3)
    plt.legend()
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    # plt.title('C = {}'.format(args.C))
    plt.savefig(plotpath+'bb_time_{}.pdf'.format(args.C))

    # # plt.figure('BB time plot')
    # plt.plot(bb_L, bb_lt, color='red')
    # plt.xlabel('initial L (nm)')
    # plt.ylabel('Average time (s)')
    # plt.xlim(0)
    # # plt.ylim(0, 0.3)
    # plt.legend()
    # plt.gca().spines['right'].set_visible(False)
    # plt.gca().spines['top'].set_visible(False)
    # # plt.title('C = {}'.format(args.C))
    # plt.savefig(plotpath+'bb_time_{}.pdf'.format(args.C))
    #
    # # plt.figure('BB time plot')
    # plt.plot(bb_L, bb_tt, color='green')
    # plt.xlabel('initial L (nm)')
    # plt.ylabel('Average time (s)')
    # plt.xlim(0)
    # # plt.ylim(0, 0.3)
    # plt.legend()
    # plt.gca().spines['right'].set_visible(False)
    # plt.gca().spines['top'].set_visible(False)
    # # plt.title('C = {}'.format(args.C))
    # plt.savefig(plotpath+'bb_time_{}.pdf'.format(args.C))


def make_trajectory_plot(args, plotpath, bb_P_trailing, bb_rate_total, hist, initial_disp, **_):
    L = 15.0
    num_steps = 400
    bx_position = np.array([[0.0, L]])
    bb_step_time = np.array([0.0])
    j = 0
    seed_for_kub_2 = 2
    np.random.seed(seed_for_kub_2)

    while j < num_steps:
        L = np.abs(bx_position[j, 1] - bx_position[j, 0])
        bb_time = np.random.exponential(1/bb_rate_total[int(L)-1])

        trailing_leading_rand = np.random.random()
        if trailing_leading_rand < bb_P_trailing[int(L)-1]:
            init_disp = -L
            am_trailing = True
        else:
            init_disp = L
            am_trailing = False

        ind = find_nearest(initial_disp, init_disp)

        disp_rand = np.random.random()
        disp_prob = 0
        for i in range(len(hist[:, ind])):
            disp_prob += hist[:, ind][i]
            if disp_rand < disp_prob:
                final_disp = initial_disp[i]
                break

        if am_trailing:
            if bx_position[j, 0] < bx_position[j, 1]:
                new_bx_position = np.array(
                    [final_disp-init_disp, 0]) + bx_position[j]
            else:
                new_bx_position = np.array(
                    [0, final_disp-init_disp]) + bx_position[j]
        else:
            if bx_position[j, 0] > bx_position[j, 1]:
                new_bx_position = np.array(
                    [final_disp-init_disp, 0]) + bx_position[j]
            else:
                new_bx_position = np.array(
                    [0, final_disp-init_disp]) + bx_position[j]

        new_L = np.abs(new_bx_position[1] - new_bx_position[0])

        j += 1

        bb_step_time = np.append(bb_step_time, bb_time)
        bx_position = np.concatenate(
            (bx_position, np.array([new_bx_position])))

    # print(bx_position)
    cumulative_time = np.append(
        np.cumsum(bb_step_time), np.cumsum(bb_step_time)[-1])
    vline_x = cumulative_time[1:-1]

    yildiz_data = np.loadtxt(
        "../data/yildiz_fig2A_tracking_data.txt", dtype=np.float64)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    plt.hlines(bx_position[:, 1], cumulative_time[:-1],
               cumulative_time[1:], color='blue')
    plot1 = plt.hlines(
        bx_position[:, 0], cumulative_time[:-1], cumulative_time[1:], color='red')
    plot2 = plt.hlines(bx_position[0, 0], cumulative_time[0],
                       cumulative_time[1], color='red', linestyle='--')
    plt.vlines(vline_x, bx_position[:-1, 1], bx_position[1:, 1], color='blue')
    plt.vlines(vline_x, bx_position[:-1, 0], bx_position[1:, 0], color='red')
    # plt.plot(yildiz_data[0], yildiz_data[1], '--',
    #          markersize=0, color="blue", alpha=0.8)
    # plt.plot(yildiz_data[0], yildiz_data[2], '--',
    #          markersize=0, color="red", alpha=0.8)
    plt.plot(yildiz_data[0], yildiz_data[1]+80, '--',
             markersize=0, color="blue", alpha=0.8)
    plt.plot(yildiz_data[0], yildiz_data[2]+80, '--',
             markersize=0, color="red", alpha=0.8)

    x1, x2, y1, y2 = [6, 9, 190, 320]
    box_v = np.array([x1, x2])
    box_h = np.array([y1, y2])
    plt.hlines(box_h, box_v[0], box_v[1], linewidth=0.5)
    plt.vlines(box_v, box_h[0], box_h[1], linewidth=0.5)

    # plt.yticks(np.append(initial_disp, 0))
    plt.ylim(np.amin(bx_position)-50, 720)
    # plt.ylim(np.amin(bx_position), 500)
    plt.xlim(0, right=24)
    plt.xlabel('Time (s)')
    plt.ylabel('On-axis position (nm)')
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.subplots_adjust(top=0.85, bottom=0.05)
    plt.legend([plot1, plot2], ["Model", "Experiment"], loc='upper left')
    plt.grid()

    ax2 = fig.add_axes([0.6, 0.08, 0.3, 0.35])
    ax2.hlines(bx_position[:, 1], cumulative_time[:-1],
               cumulative_time[1:], color='blue')
    ax2.hlines(bx_position[:, 0], cumulative_time[:-1],
               cumulative_time[1:], color='red')
    ax2.vlines(vline_x, bx_position[:-1, 1], bx_position[1:, 1], color='blue')
    ax2.vlines(vline_x, bx_position[:-1, 0], bx_position[1:, 0], color='red')
    ax2.axhline(300, color='grey', linewidth=1, alpha=0.5)
    ax2.plot(yildiz_data[0], yildiz_data[1]+80, '--',
             markersize=0, color="blue", alpha=0.8)
    ax2.plot(yildiz_data[0], yildiz_data[2]+80, '--',
             markersize=0, color="red", alpha=0.8)
    ax2.set_xlim(x1, x2)
    ax2.set_ylim(y1, y2)
    ax2.get_xaxis().set_ticks([])
    ax2.get_yaxis().set_ticks([])

    plt.savefig(plotpath+u+'trajectory_plot_' +
                params_string+'.pdf', bbox_inches='tight')

    print('Avg velocity Yildiz: ', yildiz_data[1][-1]/yildiz_data[0][-1])


def bug_checking_plots(args, plotpath, initial_disp_edge, final_disp_edge, normalized_hist, **_):
    plt.figure('From Data')
    plt.pcolor(initial_disp_edge, final_disp_edge, normalized_hist)
    plt.xlabel('Initial displacement (nm)')
    plt.ylabel('Final displacement (nm)')
    # plt.title('kb = {0:.2e}, kstk = {1:.2e}, cb = {2}, cm = {3}, ct = {4}, eqb = {5}, eqmpre = {6}, eqmpost = {7}, C = {8}'.format(args.k_b,
    # args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C))

    plt.colorbar()
    # plt.savefig(plotpath+'2dhist_initL_vs_finalL_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))

    # # Plot L to L probability density
    # plt.figure('prob density')
    # plt.legend(loc='best')
    # plt.xlabel('L')
    # plt.ylabel('probability per L')
    #
    # # Obtain L to L probability density
    # for i in range(len(P)):
    #     plt.plot(final_L, p_den_L, label=f'i is {i}')
    # # plt.savefig(plotpath+'L_to_L_prob_density_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))
    #
    # plt.figure('p_den_disp')
    # plt.plot(initial_disp, p_den_disp)
    # plt.xlabel('displacement')
    # plt.ylabel('probability density')
    # # plt.savefig(plotpath+'Probability_density_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))


def main():
    plotpath = '../plots/mc_plots/'
    if not path.exists(plotpath):
        mkdir(plotpath)
    args = get_cli_arguments()
    global u
    u = ''
    if args.underMT == False:
        u = 'u_'
    global params_string
    data_params_string = "{0:.2e}_{1:.2e}_{2}_{3}_{4}_{5}_{6}_{7}_{8}".format(
        args.k_b, args.k_stk, args.cb, args.cm, args.ct, args.eqb, args.eqmpre, args.eqmpost, args.C)
    params_string = str(args.k_ub) + '_' + data_params_string

    plt.rcParams.update({'font.size': 14})

    # args.C = -0.35
    plotting_data_file = "../data/mc_plotting_data/mc_plotting_data_" + \
        u + data_params_string + ".npz"
    bothbound_data_file = "../data/mc_bb_data/bb_exp-unbinding-constant_{}.npz".format(
        args.C)
    assert(path.exists(bothbound_data_file)), "Bothbound data missing. Need to run monte_carlo_simulation_bb.py with params exp-ub-const = {}".format(
        params.for_simulation['exp-unbinding-constant'])
    initial_disp, hist, normalized_hist, time_hists, avg_affinity_time = get_onebound_data(
        args, plotting_data_file)
    bb_L, bb_P_leading, bb_P_trailing, bb_rate_total, bb_avg_t, bb_lt, bb_tt = get_bothbound_data(
        args, bothbound_data_file)
    initial_disp_edge, final_disp_edge, final_L_bin_width, final_disp_bin_width, initial_L = make_bins_and_edges(
        initial_disp)
    probability_distribution = make_probability_distribution(**locals())
    # Sum of probability distribution
    probability_distribution_sum = integrate_2d(
        probability_distribution, final_disp_bin_width, final_disp_bin_width)
    print('final displacement prob distribution 2d integral: ',
          probability_distribution_sum)

    # Linear regression of probability distribution plot
    b, m, lin_fit = least_squares(probability_distribution, initial_disp,
                                  initial_disp, final_disp_bin_width, final_disp_bin_width)

    make_prob_dist_plot(**locals())
    make_step_length_plots(**locals())
    make_ob_time_plot(**locals())
    make_trajectory_plot(**locals())
    if args.extra == False:
        make_filtered_prob_dist_plot(**locals())
        make_bothbound_plots(**locals())
        bug_checking_plots(**locals())
    if args.plot == True:
        plt.show()


if __name__ == "__main__":
    params = importlib.import_module("params")
    main()
