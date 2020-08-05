import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../data")
import importlib
# import bb_energy_distribution

N = 100000
circle = 2*np.pi
want_L = 8
err = 1

params = importlib.import_module("params")
Ls = params.for_simulation['ls']
Lt = params.for_simulation['lt']

zero_bound_dynein = np.random.uniform(0, 2*np.pi, (N,4))

# thetas are relative angles
theta_0_arr = []
theta_1_arr = []
theta_2_arr = []
theta_3_arr = []
bb_0 = []
bb_1 = []
bb_2 = []
bb_3 = []
bb_4 = []
bb_nba = []
bb_nma = []
bb_ta = []
bb_fma = []

for i in range(N):
    phi_0 = zero_bound_dynein[i,0]
    phi_1 = zero_bound_dynein[i,1]
    phi_2 = zero_bound_dynein[i,2]
    phi_3 = zero_bound_dynein[i,3]
    x1 = Ls*np.cos(phi_0)
    y1 = Ls*np.sin(phi_0)
    x2 = x1 + Lt*np.cos(phi_1)
    y2 = y1 + Lt*np.sin(phi_1)
    x3 = x2 + Lt*np.cos(phi_2)
    y3 = y2 + Lt*np.sin(phi_2)
    x4 = x3 + Ls*np.cos(phi_3)
    y4 = y3 + Ls*np.sin(phi_3)
    L = np.sqrt(x4**2 + y4**2)

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
    # print('L: 0', L)
    if L < want_L+err and L > want_L-err:
        bb_0.append(theta_0)
        bb_1.append(theta_1)
        bb_2.append(theta_2)
        bb_3.append(theta_3)
        bb_4.append(theta_3)
        bb_nba.append(theta_0)
        if theta_1 > np.pi:
            bb_nma.append(theta_1-np.pi)
        if theta_1 < np.pi:
            bb_nma.append(theta_1+np.pi)
        if theta_2 > np.pi:
            bb_ta.append(theta_2-np.pi)
        if theta_2 < np.pi:
            bb_ta.append(theta_2+np.pi)
        if theta_3 > np.pi:
            bb_fma.append(3*np.pi-theta_3)
        if theta_3 < np.pi:
            bb_fma.append(np.pi-theta_3)

rel_angles = np.transpose(np.array([theta_0_arr, theta_1_arr, theta_2_arr, theta_3_arr]))


# for j in range(4):
#     plt.figure()
#     plt.hist(zero_bound_dynein[:,j], label='angle {}'.format(j), density=False, stacked=True, histtype='bar')
#     plt.title('Random angle {}'.format(j))
#     plt.xlabel('angles')
#     plt.legend()
#
# for j in range(4):
#     plt.figure()
#     plt.hist(rel_angles[:,j], label='angle {}'.format(j), density=False, stacked=True, histtype='bar')
#     plt.title('Relative angle {}'.format(j))
#     plt.xlabel('angles')
#     plt.legend()

plt.figure()
plt.hist(bb_nma, label='nm angle', density=False, stacked=True, histtype='bar')
plt.title('Both Bound Near Motor angle')
plt.xlabel('angles')
plt.legend()

plt.figure()
plt.hist(bb_fma, label='fm angle', density=False, stacked=True, histtype='bar')
plt.title('Both Bound Far Motor angle')
plt.xlabel('angles')
plt.legend()

plt.figure()
plt.scatter(bb_1, bb_3)
plt.title('Both Bound angle distribution')
plt.xlabel(r'$\theta_{1}$')
plt.ylabel(r'$\theta_{3}$')
plt.legend()

plt.figure()
plt.scatter(bb_nma, bb_fma)
plt.title('Both Bound motor angle distribution')
plt.xlabel(r'$\theta_{nm}$')
plt.ylabel(r'$\theta_{fm}$')
plt.legend()
plt.show()
