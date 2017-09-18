#!/user/bin/python
import simrunner
import numpy as np

l = "autocorrelationData"


simrunner.run_sim(**{"k_b": 10e-10, "k_ub": 10e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-10, "label": l,
                     "constant-write": True, "runtime": 5e-4})


simrunner.run_sim(**{"k_b": 10e-10, "k_ub": 10e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-11, "label": l,
                     "constant-write": True, "runtime": 5e-4})


simrunner.run_sim(**{"k_b": 10e-10, "k_ub": 10e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-12, "label": l,
                     "constant-write": True, "runtime": 5e-4})