#!/user/bin/python

import simrunner
import numpy as np

l = 'autocorr'

simrunner.run_sim(**{"k_b": 10e20, "k_ub": 10e-10, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-12, "label": l, "constant-write": True, "runtime": 1e-4})

simrunner.run_sim(**{"k_b": 10e-10, "k_ub": 10e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-12, "label": l, "constant-write": True, "runtime": 1e-4})
