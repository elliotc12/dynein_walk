#!/user/bin/python
import dynein.run as run
import numpy as np


l = "autocorrelationData"


basename = run.sim(**{"k_b": 10e-10, "k_ub": 10e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-10, "label": l, "constant-write": True, "runtime": 5e-4})
print basename


