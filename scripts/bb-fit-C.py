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
import scipy.optimize as optimization

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

def func(x, a, b):
    return a*x + b

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
# parser.add_argument("-L", "--L", type=float, help="displacement in nm", default=8)
parser.add_argument("-N", "--N", type=float, help="how many steps to do", default=100)
parser.add_argument("-u", "--kub", type=float, help="Manually set the unbinding const", default=params.for_simulation['k_ub'])
parser.add_argument("-C", "--C", type=float, help="Initial Guess", default=params.for_simulation['exp-unbinding-constant'])
args = parser.parse_args()

params.for_simulation['exp-unbinding-constant'] = args.C
params.for_simulation['k_ub'] = args.kub

Ls = np.arange(5, 51, 2)
delta = 0.000001

yildiz_displacements = [10.0, 20.0, 30.0, 40.0, 50.0]
yildiz_lagging_fractions = [0.525, 0.545, 0.61, 0.59, 0.67]
yildiz_lagging_uncertainty = [0.06, 0.04, 0.035, 0.045, 0.075]
yildiz_fit, pcov = optimization.curve_fit(func, yildiz_displacements, yildiz_lagging_fractions, sigma=yildiz_lagging_uncertainty)

while True:
    rate_leading = []
    rate_trailing = []
    for L in Ls:
        steps = 0
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

        seed = 0
        np.random.seed(0)

        nma_Z = np.linspace(0, 2*np.pi, np.sqrt(args.N))
        fma_Z = np.linspace(0, 2*np.pi, np.sqrt(args.N))
        NMA, FMA = np.meshgrid(nma_Z, fma_Z)
        dynein_states = bb_energy_distribution.DyneinBothBound(NMA, FMA, params, L)
        valid_states = ~np.isnan(dynein_states.E_total)
        Z = dynein_states.P[valid_states].sum()
        dynein_states.P = dynein_states.P/Z # just normalize the probabilities

        rate_trailing.append((dynein_states.rate_trailing[valid_states]*dynein_states.P[valid_states]).sum())
        rate_leading.append((dynein_states.rate_leading[valid_states]*dynein_states.P[valid_states]).sum())

    rate_trailing = np.array(rate_trailing)
    rate_leading = np.array(rate_leading)
    rel_trailing = rate_trailing/(rate_leading + rate_trailing)
    rel_leading = rate_leading/(rate_leading + rate_trailing)
    curve_fit, ppcov = optimization.curve_fit(func, Ls, rel_trailing)

    # print('{}\t{}'.format(yildiz_fit[0], curve_fit[0]))

    if abs(curve_fit[0]-yildiz_fit[0]) < delta:
        break
    elif (curve_fit[0]-yildiz_fit[0]) < 0:
        params.for_simulation['exp-unbinding-constant'] -= delta
    elif (curve_fit[0]-yildiz_fit[0]) > 0:
        params.for_simulation['exp-unbinding-constant'] += delta

print('Best C: ', params.for_simulation['exp-unbinding-constant'])

def make_hist(ax, stacked_hist, data, data0, bin, Label, Label0, tof, Color, Color0, Title, xlabel):
    ax.hist(data, bins=bin, alpha=0.5, label=Label, normed=tof, stacked=True, color=Color)
    if stacked_hist == True:
        ax.hist(data0, bins=bin, alpha=0.5, label=Label0, normed=tof, stacked=True, color=Color0)
    ax.legend(loc="upper right")
    ax.set_title(Title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")

fig6 = plt.figure(6, figsize=(8,6))
gs1 = gridspec.GridSpec(1,1)
# gs2 = gridspec.GridSpec(2,1)
ax12 = fig6.add_subplot(gs1[:,:])
# ax13 = fig6.add_subplot (gs2[1,:])


ax12.scatter(Ls, rel_trailing)
ax12.errorbar(yildiz_displacements, yildiz_lagging_fractions, yerr=yildiz_lagging_uncertainty, label="Experiment", fmt='o-', c='C0', markersize=4, linestyle='', capsize=1, elinewidth=0.3, markeredgewidth=0.3)
ax12.plot(np.array(yildiz_displacements), func(np.array(yildiz_displacements), *np.array(yildiz_fit)), 'g--', label='Experiment')
ax12.plot(np.array(Ls), func(np.array(Ls), *np.array(curve_fit)), 'r--', label='Best fit')
ax12.set_title("Trailing Unbinding Prob vs. Initial L (C=={})".format(params.for_simulation['exp-unbinding-constant']))
ax12.set_ylabel("Relative Unbinding Prob")
ax12.legend()

# ax13.scatter(Ls, rate_leading/(rate_leading + rate_trailing))
# ax13.set_title("Leading Unbinding Prob vs. Initial L")
# ax13.set_xlabel("Initial L")
# ax13.set_ylabel("Relative Unbinding Prob")

# fig7 = plt.figure(7, figsize=(6,8))
# ax14 = fig7.add_subplot(gs2[0,:])
# ax15 = fig7.add_subplot (gs2[1,:])

# ax14.scatter(Ls, rate_trailing)
# ax14.set_title("Trailing Unbinding Prob vs. Initial L")
# ax14.set_ylabel("Avg Unbinding Prob")

# ax15.scatter(Ls, rate_leading)
# ax15.set_title("Leading Unbinding Prob vs. Initial L")
# ax15.set_xlabel("Initial L")
# ax15.set_ylabel("Avg Unbinding Prob")

plt.show()
