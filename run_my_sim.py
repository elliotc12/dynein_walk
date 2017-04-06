import simrunner
import numpy as np

for cm in np.arange(0.1, 5, 0.2):
    run_sim({"k_b":  500, "k_ub": 2e11, "cb": cb, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "seed": 1})

run_sim(k_b=5, silliness="Elliott")

for run in read_csv('foo.csv'):
    run_sim(run)
