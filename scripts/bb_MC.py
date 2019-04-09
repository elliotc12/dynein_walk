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


class DyneinOneBound:
	"""
	Class for dynein one bound angles right after both bound state.
	"""


	def __init__(self, nma, fma, params, L=8, x=0):
		dyneinbb = bb_energy_distribution.DyneinBothBound(nma, fma, params, L)

		self.bba = dyneinbb.nba
		self.bma = np.pi/2 + np.arcsin((dyneinbb.r_nm[0]-dyneinbb.r_t[0])/(params.for_simulation['lt']))
		self.uba = dyneinbb.fba
		self.uma = np.pi/2 + np.arcsin((dyneinbb.r_fm[0]-dyneinbb.r_t[0])/(params.for_simulation['lt']))
			

		self.E_bba = bb_energy_distribution.spring_energy(self.bba, params.for_simulation['eqb'], params.for_simulation['cb'])
		self.E_bma = bb_energy_distribution.spring_energy(self.bma, params.for_simulation['eqmpre'], params.for_simulation['cm'])
		self.E_uba = bb_energy_distribution.spring_energy(self.uba, params.for_simulation['eqb'], params.for_simulation['cb'])
		self.E_uma = bb_energy_distribution.spring_energy(self.uma, params.for_simulation['eqmpre'], params.for_simulation['cm'])
		self.E_total = self.E_bba+self.E_bma+self.E_uba+self.E_uma



# Initialize arrays for histograms
angles = [[] for i in range(2)]		# Pair of angles
tx = []		# Tail x array
ty = []		# Tail y array
nmx = []	# Near motor array
fmx = []	# Far motor array
fbx = []	# Far bound array
E_avg_arr = []	# Average energy array
k_ub = []		# Unbinding Rates


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

	dyneinbb = bb_energy_distribution.DyneinBothBound(nma, fma, params, L)

	# Checking if energy is nan
	if np.isnan(dyneinbb.E_total) == True:
		continue
	else:
		# Calculating partition function
		P = np.exp(-b*dyneinbb.E_total)
		Z += P

		# Calculation for averages
		r_tx += dyneinbb.r_t[0]*P
		r_ty += dyneinbb.r_t[1]*P
		r_nmx += dyneinbb.r_nm[0]*P
		r_fmx += dyneinbb.r_fm[0]*P
		r_fbx += dyneinbb.r_fb[0]*P
		E_avg += dyneinbb.E_total*P

		# Array of averages
		tx.append(r_tx/Z)
		ty.append(r_ty/Z)
		nmx.append(r_nmx/Z)
		fmx.append(r_fmx/Z)
		fbx.append(r_fbx/Z)
		E_avg_arr.append(E_avg/Z)

		# One bound dynein immediately after both bound
		dyneinob = DyneinOneBound(nma, fma, params, L)

		dG = dyneinob.E_total - dyneinbb.E_total 
		print("dG: ", dG)
		k_ub.append(np.exp(-b*dG))




print("Unbinding Rates: ",k_ub )
plt.title("Unbinding Rates")
plt.hist(k_ub, bins = 100, ec = 'black')
plt.xlabel("Unbinding Rates")
plt.show()


