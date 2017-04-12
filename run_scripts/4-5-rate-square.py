#!/usr/bin/python3

import simrunner
import numpy as np

l = "4-5-rate-square"

simrunner.run_sim(**{"k_b":  500, "k_ub": 2e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 5e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 2e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 5e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 2e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 5e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 2e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 5e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 2e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 5e-10, "label": l})

simrunner.run_sim(**{"k_b":  500, "k_ub": 5e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 5e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 5e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 5e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 5e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b":  500, "k_ub": 8e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 8e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 8e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 8e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 8e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b":  500, "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b":  500, "k_ub": 4e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 4e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 4e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 4e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 4e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b":  500, "k_ub": 7e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 7e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 7e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 7e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 7e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
