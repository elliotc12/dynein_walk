#!/usr/bin/python
import simrunner
import numpy as np

for cm in np.arange(0.1, 5, 0.2):
    simrunner.run_sim(**{"k_b":  500, "k_ub": 2e11, "cb": 2.4, "cm": cm, "ct": 2.4, "dt": 1e-10, "seed": 1, "label": "test"})

simrunner.run_sim(k_b=5, k_ub=2, dt=1e-8, cb=3, label="test")

for run in simrunner.read_csv('example_run_csv.csv'):
    simrunner.run_sim(run)
