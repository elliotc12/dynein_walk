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
Monte Carlo simulation for just Bothbound dynein
"""

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
            Z['leading'] += P
            prob_unbinding['leading'].append(prob)
            prob_unbinding['leading_avg'] += prob*P
            relative_prob['leading'].append(prob/prob_unbinding['cumulative'][k[0]-len(prob_unbinding['trailing'])-1])
        else:
            # FARBOUND State - Trailing step
            Z['trailing'] += P
            prob_unbinding['trailing'].append(prob)
            prob_unbinding['trailing_avg'] += prob*P
            relative_prob['trailing'].append(prob/prob_unbinding['cumulative'][k[0]-len(prob_unbinding['leading'])-1])
        k[0]+=1


params = importlib.import_module("params")

parser = argparse.ArgumentParser()
# parser.add_argument("-L", "--L", type=float, help="displacement in nm", default=8)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=100)
parser.add_argument("-u", "--kub", type=float, help="Manually set the unbinding const", default=params.for_simulation['k_ub'])
parser.add_argument("-C", "--C", type=float, help="Exponential unbinding constant", default=params.for_simulation['exp-unbinding-constant'])
args = parser.parse_args()

prob_avg = {'trailing': [], 'leading': [], 'L': [], '1': [], '2': []}
relative_prob = {'trailing': [], 'leading': []}

for L in range(1, 52, 2):
    # L = args.L           # Initial Length
    prob_avg['L'].append(L)
    N = args.N           # Count
    steps = 0
    C =  args.C          # exponential binding constant from paper_params.py April 12
    k_ub = args.kub
    Z = {'Z': 0, 'main': 0, 'trailing': 0, 'leading': 0}                # Partition Function
    k = [0]              # Dynein Count & RNG Seed

    max_unbinding = 1
    b = 1#1.82733524          # thermodynamic beta from default_parameters.h
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
            'leading': [],
            'cumulative': []
    }
    prob_unbinding = {      # Unbinding probability
            'trailing': [],
            'leading': [],
            'cumulative': [],
            'trailing_avg': 0,
            'leading_avg': 0,
            'cumulative_avg': 0,
            'avg': 0
    }

    seed = 0
    np.random.seed(0)

    nma_Z = np.linspace(0, 2*np.pi, np.sqrt(N))
    fma_Z = np.linspace(0, 2*np.pi, np.sqrt(N))
    NMA, FMA = np.meshgrid(nma_Z, fma_Z)
    dynein_Z = bb_energy_distribution.DyneinBothBound(NMA, FMA, params, L)
    big_P = dynein_Z.P[~np.isnan(dynein_Z.P)]

    Z['Z'] = sum(big_P)

    print('Z', Z['Z'])

    while steps < N:
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

                    # FIXME!
                    rate_trailing = np.exp(C*(dynein.nba - eqb_angle))
                    rate_leading = np.exp(C*(dynein.fba - eqb_angle))
                    max_rate_leading = max(rate_leading, max_rate_leading)
                    max_rate_trailing = max(rate_trailing, max_rate_trailing)
                    cumulative_rate = rate_trailing+rate_leading

                    prob_trailing = rate_trailing/10 #(P/Z['Z'])
                    prob_leading =  rate_leading/10 #(P/Z['Z'])

                    # print('prob trailing: ', prob_trailing)
                    # print('prob leading: ', prob_leading)

                    new_nma = nma-(np.pi-dynein.nba)
                    new_fma = fma-(np.pi-dynein.fba)

                    if np.random.random() < prob_trailing: # FIXME need to normalize this a tad so it is never > 1.
                            # FARBOUND State
                            state = 1
                            rate_unbinding['trailing'].append(rate_trailing)
                            rate_unbinding['cumulative'].append(cumulative_rate)
                            prob_unbinding['cumulative'].append(prob_trailing+prob_leading)
                            collect_bothbound_data(k, dynein, P, state, nma, fma, prob_trailing)
                            Z['main'] += P
                            prob_unbinding['cumulative_avg'] += (prob_trailing+prob_leading)*P
                            steps += 1

                    if np.random.random() < prob_leading:
                            # NEARBOUND State
                            state = 0
                            rate_unbinding['leading'].append(rate_leading)
                            rate_unbinding['cumulative'].append(cumulative_rate)
                            prob_unbinding['cumulative'].append(prob_trailing+prob_leading)
                            collect_bothbound_data(k, dynein, P, state, nma, fma, prob_leading)
                            Z['main'] += P
                            prob_unbinding['cumulative_avg'] += (prob_trailing+prob_leading)*P
                            steps += 1



    prob_unbinding_avg = prob_unbinding['avg']/Z['main']
    prob_unbinding_leading_avg = prob_unbinding['leading_avg']/Z['leading']
    prob_unbinding_trailing_avg = prob_unbinding['trailing_avg']/Z['trailing']
    prob_unbinding_cumulative_avg = prob_unbinding['cumulative_avg']/Z['main']


    print('BOTHBOUND AVERAGES')
    print('Prob unbinding trailing: ', prob_unbinding['trailing'])
    print('Prob unbinding leading: ', prob_unbinding['leading'])
    # print('P: ', P_arr)
    print('AVG trailing: ', np.mean(prob_unbinding['trailing']))
    print('AVG leading: ', np.mean(prob_unbinding['leading']))
    # print('AVG P: ', np.mean(P_arr))

    # print("Avg prob_unbinding:", prob_unbinding_avg)
    # print("Z trailing:", Z['trailing'])
    # print("Z leading:", Z['leading'])
    # print("Z Main:", Z['main'])
    # print("Avg rel trailing prob_unbinding:", prob_unbinding_trailing_avg)
    # print("Avg rel leading prob_unbinding:", prob_unbinding_leading_avg)

    prob_avg['trailing'].append(np.mean(relative_prob['trailing']))
    prob_avg['1'].append(np.mean(prob_unbinding['trailing']))
    prob_avg['leading'].append(np.mean(relative_prob['leading']))
    prob_avg['2'].append(np.mean(prob_unbinding['leading']))


def make_hist(ax, stacked_hist, data, data0, bin, Label, Label0, tof, Color, Color0, Title, xlabel):
    ax.hist(data, bins=bin, alpha=0.5, label=Label, normed=tof, stacked=True, color=Color)
    if stacked_hist == True:
        ax.hist(data0, bins=bin, alpha=0.5, label=Label0, normed=tof, stacked=True, color=Color0)
    ax.legend(loc="upper right")
    ax.set_title(Title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")

fig6 = plt.figure(6, figsize=(6,8))
gs2 = gridspec.GridSpec(2,1)
ax12 = fig6.add_subplot(gs2[0,:])
ax13 = fig6.add_subplot (gs2[1,:])

ax12.scatter(prob_avg['L'], prob_avg['trailing'])
ax12.set_title("Trailing Unbinding Prob vs. Initial L")
ax12.set_ylabel("Avg Unbinding Prob")

ax13.scatter(prob_avg['L'], prob_avg['leading'])
ax13.set_title("Leading Unbinding Prob vs. Initial L")
ax13.set_xlabel("Initial L")
ax13.set_ylabel("Avg Unbinding Prob")

fig7 = plt.figure(7, figsize=(6,8))
ax14 = fig7.add_subplot(gs2[0,:])
ax15 = fig7.add_subplot (gs2[1,:])

ax14.scatter(prob_avg['L'], prob_avg['1'])
ax14.set_title("Trailing Unbinding Prob vs. Initial L")
ax14.set_ylabel("Avg Unbinding Prob")

ax15.scatter(prob_avg['L'], prob_avg['2'])
ax15.set_title("Leading Unbinding Prob vs. Initial L")
ax15.set_xlabel("Initial L")
ax15.set_ylabel("Avg Unbinding Prob")

# prob_separate_hist = make_hist(ax12, True, prob_unbinding['rel_trailing'], prob_unbinding['rel_leading'], 30,
#                     "Trailing", "Leading", True, "C0", "C1",
#                     "Unbinding Probabilities", "Probability")
# prob_hist = make_hist(ax13, False, prob_unbinding['unbinding'], None, 30,
#                     "Probabilities", None, True, "C0", None,
#                     "Cumulative Unbinding Probabilities", "Probability")
# plt.savefig('../plots/mc_plots/mc_bb_unbinding_prob_{0}_{1}.pdf'.format(N, C), transparent=False)

plt.show()
