import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../data")
import importlib

params = importlib.import_module("params")

nma = np.linspace(0, 2*np.pi, 500)
fma = np.linspace(0, 2*np.pi, 500)


def spring_energy(theta, theta_eq, c, eq_in_degrees=True): #domain angle, equilibrium angle, spring constant
    """
    Calculate the Hooke's law energy for the spring at a particular domain
    """
    if eq_in_degrees == True:
        theta_eq = theta_eq*np.pi/180
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
    r_fm = [L-Ls*np.cos(fba), Ls*np.sin(fba)]

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
    print('min', np.nanmin(E_total))
    # return a list with the correct indices
    maximum = [max_coords[0][0], max_coords[1][0]]
    minimum = [min_coords[0][0], min_coords[1][0]]
    print('energy at min', E_total[min_coords[0][0], min_coords[1][0]])
    return maximum, minimum

def get_bb_coordinates(nma, fma, L, x=0):
    """
    Calculate positions of each domain for the bothbound state given motor angles, 
    """
    # calculate all of the angles 
    Lt = params.for_simulation['lt']
    Ls = params.for_simulation['ls']
    Ln = np.sqrt(Lt**2+Ls**2-2*Lt*Ls*np.cos(nma)) # squared lengths
    Lf = np.sqrt(Lt**2+Ls**2-2*Ls*Lt*np.cos(fma))

    na = np.arccos((Ln**2+Ls**2-Lt**2)/(2*Ln*Ls))
    fa = np.arccos((Lf**2+Ls**2-Lt**2)/(2*Lf*Ls))

    nba = np.arccos((L**2+Lf**2-Ln**2)/(2*L*Lf))-na
    fba = np.pi-np.arccos((L**2+Ln**2-Lf**2)/(2*L*Ln))-fa

    # coordinates
    r_nb = [x, 0]
    r_fb = [x+L, 0]

    r_nm = [x+Ls*np.cos(nba), Ls*np.sin(nba)]
    r_fm = [x+L-Ls*np.cos(fba), Ls*np.sin(fba)]

    r_t = [x+Ln*np.cos(na+nba), Ln*np.sin(na+nba)]
    r_t = [x+Lf*np.cos(fa+fba), Lf*np.sin(fa+fba)]

    return r_nb, r_fb, r_nm, r_fm, r_t

def is_bb_crazy(nma, fma, L):
    r_nb, r_fb, r_nm, r_fm, r_t = get_bb_coordinates(nma, fma, L)
    return r_nm[1] < 0 or r_fm[1] < 0 or r_t[1] < 0

def plot_bb_energy_distribution(nma, fma, L=16, extrema=True):
    NMA, FMA = np.meshgrid(nma, fma)
    E_total = get_bb_total_energy(NMA, FMA, L)
    for i in range(NMA.shape[0]):
        for j in range(NMA.shape[1]):
            if is_bb_crazy(NMA[i,j], FMA[i,j], L):
                E_total[i,j] = np.nan

    fig = plt.figure()
    if extrema == False:
        ax1 = fig.add_subplot(1, 1, 1)
    else:
        ax1 = fig.add_subplot(1,2,1)
    ax1.pcolor(NMA, FMA, E_total)
    ax1.set_xlabel('Near motor angle')
    ax1.set_ylabel('Far motor angle')
    ax1.set_title('Total Energy distribution for L={0}'.format(L))
    ax1.set_xlim(0-0.1, 2*np.pi+0.1)
    ax1.set_ylim(0-0.1, 2*np.pi+0.1)

    if extrema == True:
        import dynein.draw.balls as balls
        # plot extrema as red dots on graph
        extrema = np.asarray(find_energy_extema(E_total))
        x = nma[extrema[:,1]]
        y = fma[extrema[:, 0]]
        print(x,y)
        print('craziness', is_bb_crazy(x[0],y[0],L))
        print('craziness', is_bb_crazy(x[1],y[1],L))
        ax1.scatter(x,y, color='r')

        #add another graph to show configuration at extrema
        r_nb1, r_fb1, r_nm1, r_fm1, r_t1 = get_bb_coordinates(x[0], y[0], L)
        r_nb2, r_fb2, r_nm2, r_fm2, r_t2 = get_bb_coordinates(x[1], y[1], L)
        ax2 = fig.add_subplot(1,2,2)

        x_coords1 = [r_nb1[0], r_nm1[0], r_t1[0], r_fm1[0], r_fb1[0]]
        y_coords1 = [r_nb1[1], r_nm1[1], r_t1[1], r_fm1[1], r_fb1[1]]

        x_coords2 = [r_nb2[0], r_nm2[0], r_t2[0], r_fm2[0], r_fb2[0]]
        # shift over second dynein to make sure they don't overlap on graph
        Delta = 3.5*L
        x_coords2 = [elem + Delta for elem in x_coords2]
        y_coords2 = [r_nb2[1], r_nm2[1], r_t2[1], r_fm2[1], r_fb2[1]]

        balls.draw(x_coords1, y_coords1, alpha=1)
        balls.draw(x_coords2, y_coords2, alpha=1)
        plt.axhline(0)

        ax2.axis('off')
        ax2.axis('equal')
        return fig, [ax1, ax2]
    else:
        return fig, ax1
if __name__ == "__main__":
   fig1, ax1 = plot_bb_energy_distribution(nma, fma, L=8, extrema=True)
   fig2, ax2 = plot_bb_energy_distribution(nma, fma, L=16)
   fig3, ax3 = plot_bb_energy_distribution(nma, fma, L=24)
   plt.show()





























