import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../data")
import importlib
import bb_energy_distribution

N = 100000
circle = 2*np.pi

params = importlib.import_module("params")
ls = params.for_simulation['ls']
lt = params.for_simulation['lt']

zero_bound_dynein = np.random.uniform(0, 2*np.pi, (N,4))

# thetas are relative angles
theta_0_arr = []
theta_1_arr = []
theta_2_arr = []
theta_3_arr = []

for i in range(N):
    phi_0 = zero_bound_dynein[i,0]
    phi_1 = zero_bound_dynein[i,1]
    phi_2 = zero_bound_dynein[i,2]
    phi_3 = zero_bound_dynein[i,3]
    theta_0 = 1*phi_0
    if phi_1 > theta_0:
        theta_1 = phi_1-theta_0
        if phi_2 > phi_1:
            theta_2 = phi_2-phi_1
            if phi_3 > phi_2:
                theta_3 = phi_3-phi_2
            else:
                theta_3 = (circle+phi_3)-phi_2
        else:
            theta_2 = (circle+phi_2)-phi_1
            if phi_3 > phi_2:
                theta_3 = (circle+phi_3)-(theta_2+theta_1+theta_0)
            else:
                theta_3 = (2*circle+phi_3)-(theta_2+theta_1+theta_0)
    else:
        theta_1 = (circle+phi_1)-theta_0
        if phi_2 > phi_1:
            theta_2 = (circle+phi_2)-(theta_1+theta_0)
            if phi_3 > phi_2:
                theta_3 = phi_3-phi_2
            else:
                theta_3 = (circle+phi_3)-phi_2
        else:
            theta_2 = (2*circle+phi_2)-(theta_1+theta_0)
            if phi_3 > phi_2:
                theta_3 = (circle+phi_3)-(phi_1+theta_2)
            else:
                theta_3 = (3*circle+phi_3)-(theta_2+theta_1+theta_0)
    theta_0_arr.append(theta_0)
    theta_1_arr.append(theta_1)
    theta_2_arr.append(theta_2)
    theta_3_arr.append(theta_3)
rel_angles = np.transpose(np.array([theta_0_arr, theta_1_arr, theta_2_arr, theta_3_arr]))

r_nb = np.zeros([2, N])
r_nm = r_nb + np.array([ls*np.cos(rel_angles[:,0]), ls*np.sin(rel_angles[:,0])])
# r_t = r_nm +
# r_fm
# r_fb
print(r_nm)

for j in range(4):
    plt.figure()
    plt.hist(zero_bound_dynein[:,j], label='angle {}'.format(j), density=False, stacked=True, histtype='bar')
    plt.title('Random angle {}'.format(j))
    plt.xlabel('angles')
    plt.legend()

for j in range(4):
    plt.figure()
    plt.hist(rel_angles[:,j], label='angle {}'.format(j), density=False, stacked=True, histtype='bar')
    plt.title('Relative angle {}'.format(j))
    plt.xlabel('angles')
    plt.legend()
plt.show()
