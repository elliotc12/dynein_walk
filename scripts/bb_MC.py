import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants
import sys
sys.path.append("../data")
import importlib
import argparse
import bb_energy_distribution

params = importlib.import_module("params")

parser = argparse.ArgumentParser()
parser.add_argument("-L", type=float, help="displacement in nm", required=True)
args = parser.parse_args()

angles = [[] for i in range(2)]		# Pair of angles

rate_unbinding_leading = []			# Leading (Far) Unbinding Rates
rate_unbinding_trailing = []		# Trailing (Near) Unbinding Rates


C =  0.1 	# exponential binding constant FIXME from parameters

Z = 0		# Partition Function
N = 100	 	# Count
L = args.L 	# Length

# These are sums of partition function
r_tx = 0		# Tail x position
r_ty = 0		# Tail y position
r_nmx = 0		# Near motor x position
r_fmx = 0		# Far motor x position
r_fbx = 0		# Far bound x position
E_avg = 0		# Energy Average

b = 1			# thermodynamic beta

while Z < N:
	# Making random motor angles
	nma = np.random.uniform(0, 2*np.pi)
	fma = np.random.uniform(0, 2*np.pi)
	angles[0].append(nma)
	angles[1].append(fma)

	dynein = bb_energy_distribution.DyneinBothBound(nma, fma, params, L)

	# Checking if energy is nan
	if np.isnan(dynein.E_total) == True:
		continue
	else:
		# Calculating partition function
		P = np.exp(-b*dynein.E_total)
		Z += P

		# Calculation for averages
		r_tx += dynein.r_t[0]*P
		r_ty += dynein.r_t[1]*P
		r_nmx += dynein.r_nm[0]*P
		r_fmx += dynein.r_fm[0]*P
		r_fbx += dynein.r_fb[0]*P
		E_avg += dynein.E_total*P

		rate_trailing = np.exp(C*(dynein.nba - params.for_simulation['eqb']))
		rate_leading = np.exp(C*(dynein.fba - params.for_simulation['eqb']))
		rate_unbinding_trailing.append(rate_trailing)
		rate_unbinding_leading.append(rate_leading)

		prob_trailing = P*rate_trailing
		prob_leading = P*rate_leading

		print("P:", P)
		print("prob_leading: ", prob_leading)
		print("prob_trailing: ", prob_trailing)
		if np.random.random() < prob_trailing:
		        print("I ought to simulate this")
		if np.random.random() < prob_leading:
		        print("I ought to simulate this")

print("rate_unbinding_leading: ", rate_unbinding_leading)
print("rate_unbinding_trailing: ", rate_unbinding_trailing)

tx = r_tx/Z		# Tail x array
ty = r_ty/Z		# Tail y array
nmx = r_nmx/Z	# Near motor array
fmx = r_fmx/Z	# Far motor array
fbx = r_fbx/Z	# Far bound array
E_avg_arr = E_avg/Z	# Average energy array

print("Avg Tail x:", tx)
print("Avg Tail y:", ty)
print("Avg nmx:", nmx)
print("Avg fmx:", fmx)
print("Avg fbx:", fbx)
print("Avg E:", E_avg_arr)

# print("Unbinding Rates: ", rate_ub)
# plt.title("Unbinding Rates")
# plt.hist(rate_ub, bins = 100, ec = 'black')
# plt.xlabel("Unbinding Rates")
# plt.show()


