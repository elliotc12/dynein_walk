import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../data")
import importlib

params = importlib.import_module("params")


eq_in_degrees=True

def spring_energy(theta, theta_eq, c): #domain angle, equilibrium angle, spring constant
    """
    Calculate the Hooke's law energy for the spring at a particular domain
    """
    if eq_in_degrees == True:
        theta_eq = theta_eq*np.pi/180
    return 0.5*c*(theta-theta_eq)**2

def unbinding_rate(expconst, theta, theta_eq):
    """
    Calculate unbinding rates for Jins poster
    """
    if eq_in_degrees == True:
        theta_eq = theta_eq*np.pi/180
    return np.exp(expconst*(theta-theta_eq))
    




class DyneinBothBound:
    """
    Class for all dynein parameters in the both bound state.
    """


    def __init__(self, nma, fma, params, L=8, x=0):
        self.nma = nma
        self.fma = fma
        self.L = L

        self.Lt = params.for_simulation['lt']
        self.Ls = params.for_simulation['ls']

        self.Ln = np.sqrt(self.Ls**2+self.Lt**2-2*self.Ls*self.Lt*np.cos(2*np.pi-self.nma))
        self.Lf = np.sqrt(self.Ls**2+self.Lt**2-2*self.Ls*self.Lt*np.cos(2*np.pi-self.fma))

        # calculate angles from Ln, Lf to microtubule
        self.na = np.arccos((self.Ln**2+self.L**2-self.Lf**2)/(2*self.Ln*self.L))
        self.fa = np.pi-np.arccos((self.Lf**2+self.L**2-self.Ln**2)/(2*self.Lf*self.L))

        # calculate small triangle angles
        I think there may be a bug here, since we may assume here that the triangle is on one side of the thingy
        self.nsa = np.arccos((self.Ln**2+self.Ls**2-self.Lt**2)/(2*self.Ln*self.Ls))
        self.fsa = np.arccos((self.Lf**2+self.Ls**2-self.Lt**2)/(2*self.Lf*self.Ls))

        # calculate binding domain angles
        self.nba = self.na-self.nsa
        self.fba = self.fa-self.fsa

        # calculate positions
        self.r_nb = np.array([x*np.ones_like(self.nma), np.zeros_like(self.nma)])

        self.r_fb = np.array([self.r_nb[0]+self.L, self.r_nb[1]])


        self.r_nm = self.r_nb + np.array([self.Ls*np.cos(self.nba), self.Ls*np.sin(self.nba)])
        self.r_fm = self.r_fb + np.array([self.Ls*np.cos(self.fba), self.Ls*np.sin(self.fba)])

        self.r_t = self.r_nb + np.array([self.Ln*np.cos(self.na), self.Ln*np.sin(self.na)])

        # calculate distance between motors
        # NOTE: np.sqrt() is already element-wise 
        self.Lm = np.sqrt((self.r_nm[0]-self.r_fm[0])**2 + (self.r_nm[1]-self.r_fm[1])**2)
        # calculate tail angle
        self.ta = np.arccos(1- self.Lm**2/self.Lt**2)


        # calculate all of the energies
        self.E_t = spring_energy(self.ta, params.for_simulation['eqt'], params.for_simulation['ct'])
        self.E_nm = spring_energy(self.nma, params.for_simulation['eqmpost'], params.for_simulation['cm'])
        self.E_fm = spring_energy(self.fma, params.for_simulation['eqmpost'], params.for_simulation['cm'])
        self.E_nb = spring_energy(self.nba, params.for_simulation['eqb'], params.for_simulation['cb'])
        self.E_fb = spring_energy(self.fba, params.for_simulation['eqb'], params.for_simulation['cb'])

        self.E_total = self.E_t+self.E_nm+self.E_fm+self.E_nb+self.E_fb

        # Hokey trick to make result a nan if y coordinates are negative
        self.E_total += 0*np.sqrt(self.r_nm[1]) + 0*np.sqrt(self.r_fm[1]) + 0*np.sqrt(self.r_t[1])

        # added unbinding prob. for Jin's poster
        self.P = np.exp(-self.E_total)
        self.rate_trailing = unbinding_rate(params.for_simulation['exp-unbinding-constant'], self.nba, params.for_simulation['eqb'])
        self.rate_leading = unbinding_rate(params.for_simulation['exp-unbinding-constant'], self.fba, params.for_simulation['eqb'])
        self.prob_trailing = self.P*self.rate_trailing
        self.prob_leading = self.P*self.rate_leading
        # self.prob = []
        # if self.prob_trailing > self.prob_leading:
        #     self.prob.append(self.prob_trailing)
        # else:
        #     self.prob.append(self.prob_leading)

    def find_energy_extrema(self):
        """
        Return indices of the maximum and minimum energy value
        """

        # make sure that we actually have an array of energies to check
        assert(self.nma.shape != (1,))

        # need to be careful about nan values. Use np.nanmax and np.nanmin to ignore nans
        max_coords = np.where(self.E_total==np.nanmax(self.E_total))
        min_coords = np.where(self.E_total==np.nanmin(self.E_total))

        # return a list with the correct indices
        maximum = [max_coords[0][0], max_coords[1][0]]
        minimum = [min_coords[0][0], min_coords[1][0]]

        nm_max = maximum[1]
        fm_max = maximum[0]
        nm_min = minimum[1]
        fm_min = minimum[0]

        return nm_max, fm_max, nm_min, fm_min


    def plot_bb_energy_distribution(self):
        """Plot the energy distribution for the both bound configuration given an
        array of motor angles and an initial displacement.
        """

        # make sure that we actually have an array of energies to check
        assert(self.nma.shape!= (1,))
        fig = plt.figure()

        # make contourf graph
        ax1 = fig.add_subplot(1, 2, 1)
        energyPlot = ax1.contourf(self.nma, self.fma, self.E_total, levels=100)
        contour = ax1.contour(self.nma, self.fma, self.E_total, np.arange(1, 30, 1), colors='w', linewidth=0)
        ax1.set_xlabel(r'$\theta_{nm}$')
        ax1.set_ylabel(r'$\theta_{fm}$')
        ax1.set_title('Total Energy distribution for L={0}'.format(self.L))
        ax1.set_xlim(0-0.1, 2*np.pi+0.1)
        ax1.set_ylim(0-0.1, 2*np.pi+0.1)
        cb = plt.colorbar(energyPlot)
        cb.set_label(r"Energy [$k_BT$]")
        cb.set_ticks(np.arange(0, 31, 5))
        cb.add_lines(contour)
        # find the extrema
        j_max, i_max, j_min, i_min = self.find_energy_extrema()
        ax1.scatter(self.nma[i_max,j_max], self.fma[i_max,j_max], color='red')
        ax1.scatter(self.nma[i_min, j_min], self.fma[i_min, j_min], color='orange')

        ax2 = fig.add_subplot(1, 2, 2)
        x_coords_min = [self.r_nb[0,i_min, j_min],
                        self.r_nm[0,i_min, j_min],
                        self.r_t[0,i_min, j_min],
                        self.r_fm[0,i_min, j_min],
                        self.r_fb[0,i_min, j_min]]
        y_coords_min = [self.r_nb[1,i_min, j_min],
                        self.r_nm[1,i_min, j_min],
                        self.r_t[1,i_min, j_min],
                        self.r_fm[1,i_min, j_min],
                        self.r_fb[1,i_min, j_min]]
        x_coords_max = [self.r_nb[0,i_max, j_max],
                        self.r_nm[0,i_max, j_max],
                        self.r_t[0,i_max, j_max],
                        self.r_fm[0,i_max, j_max],
                        self.r_fb[0,i_max, j_max]]
        y_coords_max = [self.r_nb[1,i_max, j_max],
                        self.r_nm[1,i_max, j_max],
                        self.r_t[1,i_max, j_max],
                        self.r_fm[1,i_max, j_max],
                        self.r_fb[1,i_max, j_max]]

        shift = 40
        x_coords_max = [elem + shift for elem in x_coords_max]

        ax2.plot(x_coords_min, y_coords_min, color='orange', label="min")
        ax2.plot(x_coords_max, y_coords_max, color='red', label="max")
        ax2.axis('off')
        ax2.axis('equal')
        ax2.legend()



if __name__ == "__main__":
    num_points = 500
    nma = np.linspace(0, 2*np.pi, num_points)
    fma = np.linspace(0, 2*np.pi, num_points)
    NMA, FMA = np.meshgrid(nma, fma)
    dynein_8 = DyneinBothBound(NMA, FMA, params, L=8)
    dynein_16 = DyneinBothBound(NMA, FMA, params, L=16)
    dynein_24 = DyneinBothBound(NMA, FMA, params, L=24)
    dynein_32 = DyneinBothBound(NMA, FMA, params, L=32)

    dynein_8.plot_bb_energy_distribution()
    dynein_16.plot_bb_energy_distribution()
    dynein_24.plot_bb_energy_distribution()
    dynein_32.plot_bb_energy_distribution()

    plt.show()

































