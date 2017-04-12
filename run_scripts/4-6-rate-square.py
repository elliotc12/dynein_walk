#!/usr/bin/python
import simrunner
import numpy as np

l = "4-6-rate-square"

simrunner.run_sim(**{"k_b":  100, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  300, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  500, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1500, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1700, "k_ub": 7e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b":  100, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  300, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  500, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1500, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1700, "k_ub": 2e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b":  100, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  300, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  500, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  700, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b":  900, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1100, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1300, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1500, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1700, "k_ub": 7e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
