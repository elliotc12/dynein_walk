import os
from os import path, mkdir
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

    for j in range(len(arr)):
        for i in range(len(arr[j])):
            x_mean += arr[j,i]*x[i]*dx[i]*dy[j]
            y_mean += arr[j,i]*y[j]*dx[i]*dy[j]
            x2_mean += arr[j,i]*x[i]**2*dx[i]*dy[j]
            xy_mean += arr[j,i]*y[j]*x[i]*dx[i]*dy[j]

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

def get_cli_arguments():
    parser = argparse.ArgumentParser(description = 'Script to generate various plots from Monte Carlo data.')
    parser.add_argument("-k", "--k_b", type=float, help="Binding const", default=params.for_simulation['k_b'])
    parser.add_argument("-s", "--k_stk", type=float, help="Sticky const", default=params.for_simulation['k_stk'])
    parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
    return parser.parse_args()

def get_onebound_data(args, plotting_data_file):
    # Onebound Data
    mc_data = np.load(plotting_data_file, allow_pickle=True)
    return mc_data['initial_disp'], mc_data['hist'], mc_data['normalized_hist']

def get_bothbound_data(args, bothbound_data_file):
    # Bothbound Data
    mc_bb_data = np.load(bothbound_data_file, allow_pickle=True)
    bb_L = mc_bb_data['L']
    bb_rate_leading = mc_bb_data['rate_leading']*params.for_simulation['k_ub']
    bb_rate_trailing = mc_bb_data['rate_trailing']*params.for_simulation['k_ub']
    bb_P_leading = bb_rate_leading/(bb_rate_leading+bb_rate_trailing)
    bb_P_trailing = bb_rate_trailing/(bb_rate_leading+bb_rate_trailing)
    bb_avg_t = 1/(bb_rate_leading+bb_rate_trailing)
    return bb_L, bb_P_leading, bb_P_trailing, bb_avg_t

def make_bins_and_edges(initial_disp):
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

def make_probability_distribution(hist, normalized_hist, bb_P_leading, bb_P_trailing, initial_L, final_L_bin_width, **_):
    # Transition Matrix
    T = np.matrix(hist)     # Dimensionless

    # Probability Column Vector for L
    P = np.matrix(np.zeros((int(len(T)/2),1)))      # Dimensionless column vector

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
        P[:,:]= 0
        P[i] = 1

        P_L_to_L = (L_to_L(T, P_ub_leading, P_ub_trailing)**num_steps)*P    # Dimensionless 2D Array that only has 1 column vector
        P_L_to_L = np.array(P_L_to_L)[:,0]                      # convert to a dimensionless 1D array from a column vector
        norm_const = 1/((P_L_to_L*final_L_bin_width).sum())     # Dimensions: 1/distance, sum of (P_L_to_L flat * bin width of both axis)
        p_den_L = P_L_to_L*norm_const                           # dimensions 1/distance, a probability density

    p_den_disp = L_to_initial_displacement(P_ub_leading, P_ub_trailing).dot(p_den_L)   # Dimensions: 1/distance

    # Probability Distribution is the normalized histogram multiplied by the probability density
    probability_distribution = np.zeros_like(normalized_hist)   # Dimensions: 1/distance**2
    for i in range(probability_distribution.shape[0]):
        probability_distribution[i,:] = normalized_hist[i,:]*p_den_disp

    return probability_distribution

def make_prob_dist_plot(args, probability_distribution, initial_disp_edge, final_disp_edge, initial_disp, b, m, lin_fit, **_):
    plt.figure('Probability Distribution to Match Yildiz')
    plt.pcolor(initial_disp_edge, final_disp_edge, probability_distribution)
    plt.plot(initial_disp, lin_fit, label='Model: y = ({:.3}) + ({:.3})x'.format(b,m), linestyle=":", color='r')
    plt.xlabel('initial displacement (nm)')
    plt.ylabel('final displacement (nm)')
    plt.colorbar()
    plt.legend()
    plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(args.k_b, args.k_stk))
    # plt.savefig(plotpath+'final_disp_probability_distribution_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))


def make_filtered_prob_dist_plot(args, probability_distribution, initial_disp_edge, final_disp_edge, initial_disp, final_disp_bin_width, **_):
    filtered_probability_distribution = probability_distribution*1.0    # Dimensions: 1/distance**2

    # Filter steps where final and initial displacements are within 4 nm of each other
    # FIXME use actual distances here!!!
    for i in range(len(probability_distribution)):
        for j in range(len(probability_distribution[i])):
            if i == j:
                filtered_probability_distribution[i, j] = 0.0
                if i >=4 and i <= 95:
                    filtered_probability_distribution[i-4:i+4,j] = 0.0

    # Sum of filtered probability distribution
    filt_probability_distribution_sum = integrate_2d(filtered_probability_distribution, final_disp_bin_width, final_disp_bin_width)

    # Normalize the filtered historgram
    filtered_probability_distribution /= filt_probability_distribution_sum

    # Intercept & Slope for Best-fit Filtered
    b_filt, m_filt, lin_fit_filt = least_squares(filtered_probability_distribution, initial_disp, initial_disp, final_disp_bin_width, final_disp_bin_width)

    plt.figure('Filtered Probability Distribution to Match Yildiz')
    plt.pcolor(initial_disp_edge, final_disp_edge, filtered_probability_distribution)
    plt.plot(initial_disp, lin_fit_filt, label='Model: y = ({:.3}) + ({:.3})x'.format(b_filt,m_filt), linestyle=":", color='r')
    plt.xlabel('initial displacement (nm)')
    plt.ylabel('final displacement (nm)')
    plt.colorbar()
    plt.legend()
    plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(args.k_b, args.k_stk))
    # plt.savefig(plotpath+'filtered_final_disp_probability_distribution_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))

def make_step_length_plots(args, probability_distribution, initial_disp_edge, final_disp_edge, initial_disp, final_disp_bin_width, **_):
    step_length_edge = final_disp_edge - initial_disp_edge
    plt.figure('Weird step length prob distribution')
    plt.pcolor(initial_disp_edge, step_length_edge, probability_distribution)
    plt.xlabel('initial displacement (nm)')
    plt.ylabel('final displacement (nm)')
    plt.colorbar()
    plt.legend()
    plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(args.k_b, args.k_stk))

    # STEP LENGTH (s) PLOTS CALCULATION
    num_disp = len(probability_distribution)
    initial_disp_index = np.arange(0,num_disp)
    s_den = np.zeros((2*num_disp-1))        # step length probability density axis for 1d hist
    s_arr = np.zeros_like(s_den)            # step length axis (or step length bin center)
    ds = np.zeros(2*num_disp-1)
    ds[:num_disp] = final_disp_bin_width
    ds[num_disp:] = final_disp_bin_width[1:]
    s_arr[0] = -np.sum(final_disp_bin_width)+ds[0] # set step length values
    s_arr[-1] = np.sum(final_disp_bin_width)-ds[-1]
    for i in range(1,len(ds)-1):
        s_arr[i] = s_arr[i-1] + ds[i]

    for i in initial_disp_index:
        s_range1 = np.arange(i,num_disp)
        s_range2 = np.arange(0,num_disp - i)
        s_current1 = np.zeros(np.shape(probability_distribution))
        s_current2 = np.zeros(np.shape(probability_distribution))
        s_current1[s_range1, s_range2] = probability_distribution[s_range1, s_range2]
        s_current2[s_range2, s_range1] = probability_distribution[s_range2, s_range1]
        s_den[i+len(probability_distribution)-1] = integrate_2d(s_current1, final_disp_bin_width, final_disp_bin_width) # negative step length
        s_den[len(probability_distribution)-i-1] = integrate_2d(s_current2, final_disp_bin_width, final_disp_bin_width) # positive step length
    s_den = s_den/ds # dimensions probability/length
    print('1d step length integral: ', integrate_1d(s_den, ds))

    # 2D probability distribution of step length
    s_bin_edges = np.zeros(len(s_arr)+1)
    for i in range(1,len(s_bin_edges)-1):
        s_bin_edges[i] = (s_arr[i-1] + s_arr[i])*0.5
    s_bin_edges[0] = 2*s_arr[0] - s_bin_edges[1]
    s_bin_edges[-1] = 2*s_arr[-1] - s_bin_edges[-2]

    # step length histogram from Yildiz 2012 figure 1B top head-labelled panel
    yildiz_step_lengths = np.concatenate(([-37]*1, [-35]*1, [-34]*1, [-33]*2, [-31]*2, [-30]*3, [-29]*1, [-28]*1, [-27]*4, [-26]*4, [-25]*2, [-24]*3, [-23]*4, [-21]*4, [-20]*3,
                                          [-19]*3, [-18]*5, [-17]*3, [-16]*3, [-15]*7, [-14]*5, [-13]*7, [-12]*12, [-11]*16, [-10]*14, [-9]*20, [-8]*14, [-7]*10, [-6]*9, [-5]*11,
                                          [-4]*8, [-3]*2, [4]*6, [5]*7, [6]*12, [7]*20, [8]*19, [9]*22, [10]*30, [11]*34, [12]*26, [13]*21, [14]*23, [15]*22, [16]*30, [17]*29,
                                          [18]*23, [19]*22, [20]*26, [21]*12, [22]*21, [23]*16, [24]*7, [25]*8, [26]*7, [27]*8, [28]*5, [29]*9, [30]*7, [31]*8, [32]*6, [33]*2,
                                          [34]*2, [35]*9, [36]*4, [37]*9, [38]*5, [39]*1, [40]*2, [41]*1, [42]*3, [43]*2, [44]*4, [45]*4, [46]*1, [47]*1))

    # 1D hist step length
    plt.figure('Probability Density of Step Length')
    plt.fill_between(s_arr,0*s_den, s_den, label='Model', color='C1')
    plt.hist(yildiz_step_lengths, alpha=0.5, label='Experiment', density=True, stacked=True, color='C0')
    plt.xlabel('Step Length (nm)')
    plt.ylabel('Probability Density (1/nm)')
    plt.legend()
    plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(args.k_b, args.k_stk))
    # plt.savefig(plotpath+'Probability_density_step_length_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))

    s_initial_disp_edge, s_final_disp_edge = np.meshgrid(initial_disp_edge[0], s_bin_edges)
    s_bin_width = s_bin_edges[1:] - s_bin_edges[:-1]    # a 1D array giving final displacement bin width

    s_probability_distribution = np.zeros((2*num_disp-1,num_disp))
    for j in initial_disp_index:
        s_probability_distribution[initial_disp_index+num_disp-1-j,j] = probability_distribution[initial_disp_index,j]*final_disp_bin_width
        s_probability_distribution[:,j] /= s_bin_width

    # Intercept & Slope for Best-fit of step length
    s_b, s_m, s_lin_fit = least_squares(s_probability_distribution, initial_disp, s_arr, final_disp_bin_width, s_bin_width)

    print('step length prob distribution 2d integral: ', integrate_2d(s_probability_distribution, final_disp_bin_width, s_bin_width))

    # Yildiz step length vs initial disp line
    yildiz_line = [(-0.4*x)+9.1 for x in np.asarray(initial_disp)]

    # 2D hist step_length
    plt.figure('Step length probability distribution')
    plt.pcolor(s_initial_disp_edge, s_final_disp_edge, s_probability_distribution)
    plt.plot(initial_disp, s_lin_fit, label='Model: y = ({:.3}) + ({:.3})x'.format(s_b,s_m), linestyle=":", color='r')
    plt.plot(initial_disp, yildiz_line, label='Experiment: y = (9.1) + (-0.4)x', linestyle=":", color='b')
    plt.xlabel('initial displacement (nm)')
    plt.ylabel('step length (nm)')
    plt.ylim(-50,50)
    plt.colorbar()
    plt.legend()
    plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(args.k_b, args.k_stk))
    # plt.savefig(plotpath+'step_length_probability_distribution_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))

def make_bothbound_plots(args, bb_L, bb_P_trailing, bb_avg_t, **_):
    # Bothbound PLOTS
    # Prob Lagging vs Initial L plot
    plt.figure('Prob lagging vs init L')

    yildiz_displacements = [10, 20, 30, 40, 50]
    yildiz_lagging_fractions = [0.525, 0.545, 0.61, 0.59, 0.67]
    yildiz_lagging_uncertainty = [0.06, 0.04, 0.035, 0.045, 0.075]
    plt.errorbar(yildiz_displacements, yildiz_lagging_fractions, yerr=yildiz_lagging_uncertainty, label="Experiment", fmt='o-', c='C0', linestyle='', capsize=3)

    plt.scatter(bb_L, bb_P_trailing, label='Model',color='C1')
    plt.xlabel('Binding domain separation (nm)')
    plt.ylabel('P(lagging step)')
    plt.legend()
    plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(args.k_b, args.k_stk))
    # plt.savefig(plotpath+'prob_lagging_vs_init_L_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))

    # Bothbound time plot
    plt.figure('BB time plot')
    plt.plot(bb_L, bb_avg_t,color='C0')
    plt.xlabel('initial L (nm)')
    plt.ylabel('Average time (s)')
    plt.legend()
    plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(args.k_b, args.k_stk))
    # plt.savefig(plotpath+'bb_time_{0:.2e}_{1:.2e}.pdf'.format(float(args.k_b), float(args.k_stk)))


def bug_checking_plots(args, initial_disp_edge, final_disp_edge, normalized_hist, **_):
    plt.figure('From Data')
    plt.pcolor(initial_disp_edge, final_disp_edge, normalized_hist)
    plt.xlabel('initial displacement (nm)')
    plt.ylabel('final displacement (nm)')
    plt.title('kb = {0:.2e}, kstk = {1:.2e}'.format(args.k_b, args.k_stk))
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
    plotting_data_file = "../data/mc_plotting_data/mc_plotting_data_{0:.2e}_{1:.2e}.npz".format(args.k_b, args.k_stk)
    bothbound_data_file = "../data/mc_bb_data/bb_exp-unbinding-constant_{}.npz".format(args.C)
    assert(path.exists(bothbound_data_file)), "Bothbound data missing. Need to run monte_carlo_simulation_bb.py with params exp-ub-const = {}".format(params.for_simulation['exp-unbinding-constant'])
    initial_disp, hist, normalized_hist = get_onebound_data(args, plotting_data_file)
    bb_L, bb_P_leading, bb_P_trailing, bb_avg_t = get_bothbound_data(args, bothbound_data_file)
    initial_disp_edge, final_disp_edge, final_L_bin_width, final_disp_bin_width, initial_L = make_bins_and_edges(initial_disp)
    probability_distribution = make_probability_distribution(**locals())
    # Sum of probability distribution
    probability_distribution_sum = integrate_2d(probability_distribution, final_disp_bin_width, final_disp_bin_width)
    print('final displacement prob distribution 2d integral: ', probability_distribution_sum)

    # Linear regression of probability distribution plot
    b, m, lin_fit = least_squares(probability_distribution, initial_disp, initial_disp, final_disp_bin_width, final_disp_bin_width)

    make_prob_dist_plot(**locals())
    # make_filtered_prob_dist_plot(args, probability_distribution, initial_disp_edge, final_disp_edge, initial_disp, final_disp_bin_width)
    make_step_length_plots(**locals())
    make_bothbound_plots(**locals())
    bug_checking_plots(**locals())



if __name__ == "__main__":
    params = importlib.import_module("params")
    main()
    # plt.show()
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
