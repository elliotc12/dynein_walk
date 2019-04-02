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


L = args.L

# Initialize arrays for histograms
angles = [[] for i in range(2)]		# Pair of angles
tx = []		# Tail x array
ty = []		# Tail y array
nmx = []	# Near motor array
fmx = []	# Far motor array
fbx = []	# Far bound array
E_avg_arr = []	# Average energy array


Z = 0		# Partition Function
N = 100	 	# Count

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

		# Array of averages
		tx.append(r_tx/Z)
		ty.append(r_ty/Z)
		nmx.append(r_nmx/Z)
		fmx.append(r_fmx/Z)
		fbx.append(r_fbx/Z)
		E_avg_arr.append(E_avg/Z)


print("Angles: ", angles)
plt.title("Average energies")
plt.hist(E_avg_arr, bins = 100, ec = 'black')
plt.xlabel("Energy")
plt.show()



