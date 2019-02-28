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

def get_energy_landscape(nma, fma, L=16): # length in nm, near motor angle, far motor angle
    """ Calculate total energy for both bound dynein given motor angles and binding domain displacement."""
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

    E_t = spring_energy(ta, params.for_simulation['eqt'], params.for_simulation['ct'])
    print(E_t)

def remove_nans():
    pass



if __name__ == "__main__":
    get_energy_landscape(nma, fma)
