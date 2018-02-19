#!/usr/binpython3

import os
import dynein.run as run
import glob
import numpy as np


k_b = 300000000.0
k_ub = 50
runtime = 0.5
exp_binding_constant = 0

if os.path.exists('parameter_search.py'):
    os.chdir('../')
os.system("make generate_stepping_data")


basename = run.sim(**{"k_b": k_b,
                      "k_ub": k_ub,
                      "cb": 0.1,
                      "cm": 0.4,
                      "ct": 0.2,
                      "ls": 10.49,
                      "lt": 23.8,
                      "eqb": 120,
                      "eqmpre": 200,
                      "eqmpost": 224,
                      "eqt": 0,
                      "dt": 1e-10,
                      "label": "paramSearch",
                      "seed": 1,
                      "runtime": runtime,
                      "framerate": 1e-10,
                      "crash-movie": False,
                      "nomovie": True,
                      "exp_binding_constant": exp_binding_constant})


dataFile = glob.blob("data/stepping_data_paramSearch*.txt")
if len(dataFile) is not 1:
    print("Something went wrong. Make sure data was generated and that old paramSearch files have been deleted or renamed")

step_times = []
onebound_times = []
bothbound_times = []
step_lengths = []

data = np.loadtxt(dataFile)
print("Data loaded...")
bind_times = np.array(data[:,1])
unbind_times = np.array(data[:,0])
near_positions = np.around(np.array(data[:,2]), decimals=7)
far_positions = np.around(np.array(data[:,3]), decimals=7)
near_step_lens = near_positions[1:] - near_positions[:-1]
far_step_lens = far_positions[1:] - far_positions[:-1]

onebound_times = bind_times[1:]-unbind_times[1:]
bothbound_times = unbind_times[1:]-bind_times[:-1]
step_lengths = near_step_lens + far_step_lens
step_times = onebound_times + bothbound_times

print("Calculating...")
max_ob_t = np.amax(onebound_times)
min_ob_t = np.amin(onebound_times)
max_bb_t = np.amax(bothbound_times)
min_bb_t = np.amin(bothbound_times)
max_nb_step = np.amax(near_step_lens)
min_nb_step = np.amin(near_step_lens)
max_fb_step = np.amax(far_step_lens)
min_fb_step = np.amin(far_step_lens)
total_steps = len(bind_times)
nb_disp = near_positions[-1]-near_positions[0]
fb_disp = far_positions[-1]-far_positions[0]

print("Writing to testedParameters.txt")

with open("data/testedParameters.txt", "a") as file:
    file.write(k_b, k_ub, runtime, exp_binding_constant, max_ob_t, min_ob_t, max_bb_t, min_bb_t, max_nb_step, min_nb_step, max_fb_step, min_fb_step, total_steps, nb_disp, fb_disp)

print("All done.")
