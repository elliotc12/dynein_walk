#!/usr/bin/python
import simrunner
import numpy as np

l = "4-6-tiny-rate-square"

simrunner.run_sim(**{"k_b":  1300, "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  1400, "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  1200, "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  1100, "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b":  1300, "k_ub": 3e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  1400, "k_ub": 3e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  1200, "k_ub": 3e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  1100, "k_ub": 3e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
