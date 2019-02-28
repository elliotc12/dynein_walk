import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../data")
import importlib

params = importlib.import_module("params")

nma = np.linspace(0, 2*np.pi, 500)
fma = np.linspace(0, 2*np.pi, 500)


def spring_energy(theta, theta_eq, c): #domain angle, equilibrium angle, spring constant
    """
    Calculate the Hooke's law energy for the spring at a particular domain
    """
    return 0.5*c*(theta-theta_eq)**2

def get_bb_total_energy(nma, fma, L=16): # length in nm, near motor angle, far motor angle
    """ Calculate total energy for both bound dynein given motor angles and binding domain displacement."""

    # calculate all of the angles 
    Lt = params.for_simulation['lt']
    Ls = params.for_simulation['ls']
    Ln = np.sqrt(Lt**2+Ls**2-2*Lt*Ls*np.cos(nma)) # squared lengths
    Lf = np.sqrt(Lt**2+Ls**2-2*Ls*Lt*np.cos(fma))

    na = np.arccos((Ln**2+Ls**2-Lt**2)/(2*Ln*Ls))
    fa = np.arccos((Lf**2+Ls**2-Lt**2)/(2*Lf*Ls))

    nba = np.arccos((L**2+Lf**2-Ln**2)/(2*L*Lf))-na
    fba = np.pi-np.arccos((L**2+Ln**2-Lf**2)/(2*L*Ln))-fa

    r_nm = [Ls*np.cos(nba), Ls*np.sin(nba)]
    r_fm = [-Ls*np.cos(fba), Ls*np.sin(fba)]

    Lm = np.sqrt((r_nm[0]-r_fm[0])**2+(r_nm[1]-r_fm[1])**2)

    ta = np.arccos(1-(Lm**2)/(2*Lt**2))

    # calculate all of the energies
    E_t = spring_energy(ta, params.for_simulation['eqt'], params.for_simulation['ct'])
    E_nm = spring_energy(nma, params.for_simulation['eqmpost'], params.for_simulation['cm'])
    E_fm = spring_energy(fma, params.for_simulation['eqmpost'], params.for_simulation['cm'])
    E_nb = spring_energy(nba, params.for_simulation['eqb'], params.for_simulation['cb'])
    E_fb = spring_energy(fba, params.for_simulation['eqb'], params.for_simulation['cb'])

    E_total = E_t+E_nm+E_fm+E_nb+E_fb

    return E_total

def remove_nans(E_total, nma, fma):
    nan_locations = np.isnan(E_total)
    print(row_nans)


def find_energy_extema(E_total):
    """
    Return indices of the maximum and minimum energy value
    """
    # need to be careful about nan values
    max_coords = np.where(E_total==np.nanmax(E_total))
    min_coords = np.where(E_total==np.nanmin(E_total))
    # return a list with the correct indices
    maximum = [max_coords[0][0], max_coords[1][0]]
    minimum = [min_coords[0][0], min_coords[1][0]]
    print(maximum)
    return maximum, minimum


def plot_bb_energy_distribution(nma, fma, L=16, extrema=True):
    NMA, FMA = np.meshgrid(nma, fma)
    E_total = get_bb_total_energy(NMA, FMA, L)

    fig = plt.figure()
    ax = plt.gca()
    ax.pcolor(NMA, FMA, E_total)
    ax.set_xlabel('Near motor angle')
    ax.set_ylabel('Far motor angle')
    ax.set_title('Total Energy distribution for L={0}'.format(L))

    if extrema == True:
        extrema = np.asarray(find_energy_extema(E_total))
        x = nma[extrema[:,0]]
        y = fma[extrema[:, 1]]
        ax.scatter(x,y, color='r')

    return fig, ax
if __name__ == "__main__":
   fig1, ax1 = plot_bb_energy_distribution(nma, fma, L=8, extrema=True)
   fig2, ax2 = plot_bb_energy_distribution(nma, fma, L=16)
   fig3, ax3 = plot_bb_energy_distribution(nma, fma, L=24)
   plt.show()





























