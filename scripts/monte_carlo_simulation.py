import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy.constants
import sys
sys.path.append("../data")
import importlib
import argparse
import subprocess
import bb_energy_distribution

"""
Monte Carlo simulation for dynein taking a step
"""

def run_onebound(bba, bma, uma, uba, state, k):
        """
        Runs onebound.cpp with bb configuration and params.py
        """
        global seed
        seed += 1 # use a different seed every time.  ugh, global variables!
        print('running with inputs', bba, bma, uma, uba, state, k)
        process = subprocess.Popen(['../onebound',
                                    str(k_b),
                                    str(params.for_simulation['cb']),
                                    str(params.for_simulation['cm']),
                                    str(params.for_simulation['ct']),
                                    str(params.for_simulation['ls']),
                                    str(params.for_simulation['lt']),
                                    str(params.for_simulation['rt']),
                                    str(params.for_simulation['rm']),
                                    str(params.for_simulation['rb']),
                                    str(seed),
                                    str(dt),
                                    str(params.for_simulation['eqb']),
                                    str(params.for_simulation['eqmpre']),
                                    str(params.for_simulation['eqmpost']),
                                    str(params.for_simulation['eqt']),
                                    str(params.for_simulation['force']),
                                    str(params.for_simulation['exp-unbinding-constant']),
                                    str(bba), str(bma), str(uma), str(uba), str(state), str(k),
        ], stdout=subprocess.PIPE)
        (output, err) = process.communicate()
        exit_code = process.wait()
        output_data = eval(output.decode('utf-8'))
        # print('onebound for a step (which kind?) gives:\n%s\nEND OUTPUT' % output.decode('utf-8'))
        # print('data', output_data)
        assert(exit_code == 0);
        return output_data


def collect_bothbound_data(k, self, P, state, nma, fma, prob):
        """
        Collect bothbound statistics
        """
        P_arr.append(P)
        # Storing bb angles
        angles['nma'].append(nma)
        angles['fma'].append(fma)

        # Sum calculation for averages
        r_t['x_avg'] += self.r_t[0]*P
        r_t['y_avg'] += self.r_t[1]*P
        r_nm['x_avg'] += self.r_nm[0]*P
        r_nm['y_avg'] += self.r_nm[1]*P
        r_fm['x_avg'] += self.r_fm[0]*P
        r_fm['y_avg'] += self.r_fm[1]*P
        E['avg'] += self.E_total*P
        prob_unbinding['avg'] +=  prob*P

        # Storing data for histograms
        r_t['x'].append(self.r_t[0])
        r_t['y'].append(self.r_t[1])
        r_nm['x'].append(self.r_nm[0])
        r_nm['y'].append(self.r_nm[1])
        r_fm['x'].append(self.r_fm[0])
        r_fm['y'].append(self.r_fm[1])
        E['bb'].append(self.E_total)
        # FIXME These are wrong I'm pretty sure
        if state == 0:
            # NEARBOUND State - Leading step
            prob_unbinding['leading'].append(prob)
            if args.bb == True:
                data_file_bb_leading.write("{0:f}\t{1:f}\n".format(prob_unbinding['leading'][k[0]-len(prob_unbinding['trailing'])-1],
                prob_unbinding['unbinding'][k[0]]))
        else:
            # FARBOUND State - Trailing step
            prob_unbinding['trailing'].append(prob)
            if args.bb == True:
                data_file_bb_trailing.write("{0:f}\t{1:f}\n".format(prob_unbinding['trailing'][k[0]-len(prob_unbinding['leading'])-1],
                prob_unbinding['unbinding'][k[0]]))

def collect_onebound_data(k, state, bba, bma, uma, uba, L, step_data):
        """
        Call run_onebound function and collect onebound statistics
        """
        print('\n\nbothbound angles ',bba, bma, uma, uba, state)
        step = run_onebound(bba, bma, uma, uba, state, k[0])

        final_data['L'].append(step['L'])          # Final L array
        step_data['t'].append(step['t'])           # Specific type of step time array
        final_data['t'].append(step['t'])          # Onebound time array

        if state == 0:
            # NEARBOUND State - Leading step ddata
            print('leading stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            step_data['L'].append(step['L'])
            data['final_L'].append(step['L'])
            step_data['step_length'].append(step['L']-L)
            final_data['step_length'].append(step['L']-L)         # Contains both steps data

            # Storing final motor angles
            if step['L'] < 0:   #FIXME uma != nma !!! need to calculate new bb nma
                final_data['nma'].append(step['uma'])
                final_data['fma'].append(step['bma'])
            else:
                final_data['nma'].append(step['bma'])
                final_data['fma'].append(step['uma'])

            if args.ob == True:
                data_file_ob_leading.write("{0:f}\t{1:e}\n".format(final_data['L'][k[0]],final_data['t'][k[0]]))
        else:
            # FARBOUND State - Trailing step data
            print('trailing stepped with final displacement %g after time %g \n' % (step['L'], step['t']))
            step_data['L'].append(step['L'])
            data['final_L'].append(step['L'])
            step_data['step_length'].append(step['L']+L)
            final_data['step_length'].append(step['L']+L)         # Contains both steps data

            # Storing final motor angles
            if step['L'] > 0:
                final_data['nma'].append(step['bma'])
                final_data['fma'].append(step['uma'])
            else:
                final_data['nma'].append(step['uma'])
                final_data['fma'].append(step['bma'])

            if args.ob == True:
                data_file_ob_trailing.write("{0:f}\t{1:e}\n".format(final_data['L'][k[0]],final_data['t'][k[0]]))
        k[0]+=1


def plot_bb_before_step(self, dynein_color_nm, dynein_color_fm):
        """
        Plot just the figure of dynein for the both bound configuration given an
        array of motor angles and an initial displacement before the step.
        """
        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        x_coords_nm = [self.r_nb[0],
                        self.r_nm[0],
                        self.r_t[0]]

        y_coords_nm = [self.r_nb[1],
                        self.r_nm[1],
                        self.r_t[1]]

        x_coords_fm = [self.r_t[0],
                        self.r_fm[0],
                        self.r_fb[0]]

        y_coords_fm = [self.r_t[1],
                        self.r_fm[1],
                        self.r_fb[1]]


        ax.plot(x_coords_nm, y_coords_nm, color= dynein_color_nm, linewidth=3)
        ax.plot(x_coords_fm, y_coords_fm, color= dynein_color_fm, linewidth=3)
        ax.plot([-6, 30], [0, 0], color = 'black', linestyle='-', linewidth=3)
        ax.axis('off')
        ax.axis('equal')
        ax.legend()


def plot_bb_after_step(nbx, nby, nmx, nmy, tx, ty, fmx, fmy, fbx, fby, dynein_color_nm, dynein_color_fm):
        """
        Plot just the figure of dynein for the both bound configuration given an
        array of motor angles and an initial displacement after the step.
        """
        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        x_coords_nm = [nbx,
                        nmx,
                        tx]

        y_coords_nm = [nby,
                        nmy,
                        ty]

        x_coords_fm = [tx,
                        fmx,
                        fbx]

        y_coords_fm = [ty,
                        fmy,
                        fby]


        ax.plot(x_coords_nm, y_coords_nm, color= dynein_color_nm, linewidth=3)
        ax.plot(x_coords_fm, y_coords_fm, color= dynein_color_fm, linewidth=3)
        ax.plot([-6, 30], [0, 0], color = 'black', linestyle='-', linewidth=3)
        ax.axis('off')
        ax.axis('equal')
        ax.legend()


params = importlib.import_module("params")

parser = argparse.ArgumentParser()
# parser.add_argument("-L", "--L", type=float, help="displacement in nm", default=8)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=100)
parser.add_argument("-k", "--kb", type=float, help="Manually set the binding rate", default=params.for_simulation['k_b'])
parser.add_argument("-t", "--dt", type=float, help="Manually set the dt", default=params.for_simulation['dt'])
parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
parser.add_argument("--bb", type=bool, help="Collect Bothbound data", default=True)
parser.add_argument("--ob", type=bool, help="Colelct Onebound data", default=True)
args = parser.parse_args()


k_b = args.kb        # Binding Rate Constant
dt = args.dt         # Time Step
C = args.C         # exponential binding constant from paper_params.py April 12
data = {'init_L': [], 'final_L': []}

for L in range(0, 24, 2):
    # L = args.L           # Initial Length
    N = args.N           # Count
    Z = 0                # Partition Function
    k = [0]              # Dynein Count & RNG Seed

    max_unbinding = 1
    b = 11.82733524          # thermodynamic beta from default_parameters.h
    eqb_angle = params.for_simulation['eqb']
    if bb_energy_distribution.eq_in_degrees:
            eqb_angle = eqb_angle*np.pi/180

    rate_unbinding_leading = []                 # Leading (Far) Unbinding Rates
    rate_unbinding_trailing = []                # Trailing (Near) Unbinding Rates

    max_rate_trailing = 0
    max_rate_leading = 0

    # Bothbound Data
    P_arr = []
    angles = {'nma': [],'fma': []}          # Pair of Angles
    r_t = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}     # Tail data
    r_nm = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}    # Near motor data
    r_fm = {'x': [], 'y': [], 'x_avg': 0, 'y_avg': 0}    # Far motor data
    E = {'bb': [], 'avg': 0}               # Energy data
    rate_unbinding = {      # Unbinding rate
            'trailing': [],
            'leading': []
    }
    prob_unbinding = {      # Unbinding probability
            'trailing': [],
            'leading': [],
            'unbinding': [],
            'avg': 0
    }

    # Onebound Data
    trailing_data = {   # Trailing step data
            'L': [],
            't': [],
            'step_length': [],
    }
    leading_data = {    # Leading step data
            'L': [],
            't': [],
            'step_length': [],
    }
    final_data = {
            'nma': [],
            'fma': [],
            'L': [],
            't': [],
            'step_length': [],
            'L_avg': 0,
            't_avg': 0,
            'step_length_avg': 0,
    }

    seed = 0
    np.random.seed(0)

    # Creating Data File for Specific L
    if args.bb == True:
        data_file_bb_trailing = open("../data/mc_data_kb_{0:e}/mc_bb_trailing_data_{1}_{2}_{3}_{4}.txt".format(k_b, int(L), k_b, dt, N), "w")
        data_file_bb_trailing.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}, N-{3}, C-{4}********\n\n\n".format(L,
                        k_b, dt, N, C))
        data_file_bb_trailing.write("unbinding prob\t cumulative unbinding prob\n")
        data_file_bb_leading = open("../data/mc_data_kb_{0:e}/mc_bb_trailing_data_{1}_{2}_{3}_{4}.txt".format(k_b, int(L), k_b, dt, N), "w")
        data_file_bb_leading.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}, N-{3}, C-{4}********\n\n\n".format(L,
                        k_b, dt, N, C))
        data_file_bb_leading.write("unbinding prob\t cumulative unbinding prob\n")
    if args.ob == True:
        data_file_ob_trailing = open("../data/mc_data_kb_{0:e}/mc_ob_trailing_data_{1}_{2}_{3}_{4}.txt".format(k_b, int(L), k_b, dt, N), "w")
        data_file_ob_trailing.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}, N-{3}, C-{4}********\n\n\n".format(L,
                        k_b, dt, N, C))
        data_file_ob_trailing.write("final L\t t\n")
        data_file_ob_leading = open("../data/mc_data_kb_{0:e}/mc_ob_leading_data_{1}_{2}_{3}_{4}.txt".format(k_b, int(L), k_b, dt, N), "w")
        data_file_ob_leading.write("#********mc_data: L-{0}, k_b-{1}, dt-{2}, N-{3}, C-{4}********\n\n\n".format(L,
                        k_b, dt, N, C))
        data_file_ob_leading.write("final L\t t\n")

    while Z < N:
            # Making random motor angles
            nma = np.random.uniform(0, 2*np.pi)
            fma = np.random.uniform(0, 2*np.pi)

            dynein = bb_energy_distribution.DyneinBothBound(nma, fma, params, L)

            # Checking if energy is nan
            if np.isnan(dynein.E_total) == True:
                    continue
            else:
                    # Calculating partition function
                    P = np.exp(-b*dynein.E_total)
                    Z += P

                    rate_trailing = np.exp(C*(dynein.nba - eqb_angle))
                    rate_leading = np.exp(C*(dynein.fba - eqb_angle))
                    rate_unbinding['trailing'].append(rate_trailing)
                    rate_unbinding['leading'].append(rate_leading)
                    max_rate_leading = max(rate_leading, max_rate_leading)
                    max_rate_trailing = max(rate_trailing, max_rate_trailing)

                    prob_trailing = P*rate_trailing
                    prob_leading = P*rate_leading

                    # print("prob_leading: ", prob_leading)
                    # print("prob_trailing: ", prob_trailing)

                    new_nma = nma-(np.pi-dynein.nba)
                    new_fma = fma-(np.pi-dynein.fba)

                    if np.random.random() < prob_trailing: # FIXME need to normalize this a tad so it is never > 1.
                            # FARBOUND State
                            state = 1
                            data['init_L'].append(-L)

                            if (prob_trailing+prob_leading) > max_unbinding:
                                max_unbinding = prob_trailing+prob_leading
                            prob_unbinding['unbinding'].append(prob_trailing+prob_leading)
                            collect_bothbound_data(k, dynein, P, state, nma, fma, prob_trailing)


                            collect_onebound_data(k, state, dynein.fba, new_fma, new_nma, dynein.nba,
                                                    L, trailing_data)

                            # plot_bb_before_step(dynein, 'red', 'blue')
                            # plt.savefig('../plots/mc_plots/trailing_{}a_before_step.png'.format(k), transparent=False)

                            # plot_bb_after_step(step['ubx'], step['uby'], step['umx'], step['umy'],
                            #                 step['tx'], step['ty'], step['bmx'], step['bmy'],
                            #                 step['bbx'], step['bby'], 'red', 'blue')
                            # plt.savefig('../plots/mc_plots/trailing_{}b_after_step.png'.format(k), transparent=False)
                            # plt.show()


                    if np.random.random() < prob_leading:
                            # NEARBOUND State
                            state = 0
                            data['init_L'].append(L)

                            if (prob_trailing+prob_leading) > max_unbinding:
                                max_unbinding = prob_trailing+prob_leading
                            prob_unbinding['unbinding'].append(prob_trailing+prob_leading)
                            collect_bothbound_data(k, dynein, P, state, nma, fma, prob_leading)


                            collect_onebound_data(k, state, dynein.nba, new_nma, new_fma, dynein.fba,
                                                    L, leading_data)

                            # plot_bb_before_step(dynein, 'red', 'blue')
                            # plt.savefig('../plots/mc_plots/leading_{}a_before_step.png'.format(k), transparent=False)

                            # plot_bb_after_step(step['bbx'], step['bby'], step['bmx'], step['bmy'],
                            #                 step['tx'], step['ty'], step['umx'], step['umy'],
                            #                 step['ubx'], step['uby'], 'red', 'blue')
                            # plt.savefig('../plots/mc_plots/leading_{}b_after_step.png'.format(k), transparent=False)
                            # plt.show()



    print("FINAL DISPLACEMENTS: {0} \n".format(final_data['L']))
    for i in range(len(final_data['L'])):
        final_data['L_avg'] += final_data['L'][i]*P_arr[i]
        final_data['t_avg'] += final_data['t'][i]*P_arr[i]
        final_data['step_length_avg'] += final_data['step_length'][i]*P_arr[i]

    # What to collect and output or visualize?

    ### Bothbound data
    # Mean angles while bothbound? (no stepping required)
    # Mean motor/tail domain locations


    ### Stepping data (separately for leading/trailing)
    # Final displacement (mean/histogram/list)
    # Onebound time (mean/histogram/list)
    # Rate of stepping

    # print("rate_unbinding_leading: ", rate_unbinding_leading)
    # print("rate_unbinding_trailing: ", rate_unbinding_trailing)
    # print('max_rate_trailing', max_rate_trailing)
    # print('max_rate_leading', max_rate_leading)

    ### What to export, and in what format?
    # Histograms of final displacements?

    # Averages
    tx_avg = r_t['x_avg']/Z          # Tail x
    ty_avg = r_t['y_avg']/Z          # Tail y
    nmx_avg = r_nm['x_avg']/Z        # Near motor x
    nmy_avg = r_nm['y_avg']/Z        # Near Motor y
    fmx_avg = r_fm['x_avg']/Z        # Far motor x
    fmy_avg = r_fm['y_avg']/Z        # Far Motor y
    E_avg = E['avg']/Z          # Average energy
    prob_unbinding_avg = prob_unbinding['avg']/Z
    final_L_avg = final_data['L_avg']/Z
    step_length_avg = final_data['step_length_avg']/Z
    obt_avg = final_data['t_avg']/Z

    print("BOTHBOUND AVERAGES")
    print("Avg Tail x:", tx_avg)
    print("Avg Tail y:", ty_avg)
    print("Avg nmx:", nmx_avg)
    print("Avg nmy:", nmy_avg)
    print("Avg fmx:", fmx_avg)
    print("Avg fmy:", fmy_avg)
    print("Avg E:", E_avg)
    print("Avg prob_unbinding:", prob_unbinding_avg)
    print("Avg Final Displacement:", final_L_avg)
    print("Avg Step Length:", step_length_avg)
    print("Avg ob time:", obt_avg)

    data_file_bb.write("\n\n Average Unbinding Prob: {0:f}".format(prob_unbinding_avg))
    # data_file_ob.write("\n\n Average Final Displacement: {0:f}\nAverage Step Length: {1:f}\nAverage time: {2:f}".format(final_L_avg,
    #                     step_length_avg, obt_avg))

    if args.bb == True:
        data_file_bb_trailing.close()
        data_file_bb_leading.close()
    if args.ob == True:
        data_file_ob_trailing.close()
        data_file_ob_leading.close()

    # plot_hist(L, k_b, dt, N)

plt.hist2d(data['init_L'], data['final_L'], bins=(range(-30,30), 60), cmap=plt.cm.jet)
plt.show()

def make_hist(ax, stacked_hist, data, data0, bin, Label, Label0, tof, Color, Color0, Title, xlabel):
    ax.hist(data, bins=bin, alpha=0.5, label=Label, normed=tof, stacked=True, color=Color)
    if stacked_hist == True:
        ax.hist(data0, bins=bin, alpha=0.5, label=Label0, normed=tof, stacked=True, color=Color0)
    ax.legend(loc="upper right")
    ax.set_title(Title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")

def plot_hist(L, k_b, dt, N):
    fig0 = plt.figure(0, figsize=(12,8))
    gs0 = gridspec.GridSpec(2, 21)
    ax0 = fig0.add_subplot(gs0[0, 0:10])
    ax1 = fig0.add_subplot(gs0[1, 0:10])
    ax2 = fig0.add_subplot(gs0[0, 11:21])
    ax3 = fig0.add_subplot(gs0[1, 11:21])

    separate_step_hist = make_hist(ax0, True, trailing_data['L'], leading_data['L'], 30,
                        "Trailing Step", "Leading Step", True, "C0", "C1",
                        "Initial Displacement {0}nm".format(int(L)), "Final Displacement (nm)")
    step_hist = make_hist(ax1, False, final_data['step_length'], None, 50,
                        None, None, True, "C3", None,
                        "", "Step Length (nm)")
    separate_time_hist = make_hist(ax2, True, trailing_data['t'], leading_data['t'], 30,
                        "Trailing time", "Leading time", False, "C0", "C1",
                        "k_b: {0:e}".format(k_b), "time (s)")
    time_hist = make_hist(ax3, False, final_data['t'], None, 50,
                        None, None, False, "C3", None,
                        "", "time (s)")
    plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_onebound_length_time.pdf'.format(int(L), k_b, dt, N), transparent=False)


    fig6 = plt.figure(6, figsize=(6,8))
    ax12 = fig6.add_subplot(gs2[0,:])
    ax13 = fig6.add_subplot (gs2[1,:])
    prob_separate_hist = make_hist(ax12, True, prob_unbinding['trailing'], prob_unbinding['leading'], 30,
                        "Trailing", "Leading", True, "C0", "C1",
                        "Unbinding Probabilities", "Probability")
    prob_hist = make_hist(ax13, False, prob_unbinding['unbinding'], None, 30,
                        "Probabilities", None, True, "C0", None,
                        "Cumulative Unbinding Probabilities", "Probability")
    plt.savefig('../plots/mc_plots/mc_{0}_{1:e}_{2}_{3}_bothbound_unbinding_prob.pdf'.format(int(L), k_b, dt, N), transparent=False)



    # plt.show()
