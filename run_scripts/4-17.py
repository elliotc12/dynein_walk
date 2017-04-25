#!/usr/bin/python
import simrunner
import numpy as np

l = "4-23"

simrunner.run_sim(**{"k_b": 1e-6, "k_ub":  1e5, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6, "k_ub":  1e6, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6, "k_ub":  1e7, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6,  "k_ub": 1e8, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6,  "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6,  "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6,  "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6,  "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6,  "k_ub": 1e13, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b": 1e-12, "k_ub":  1e5, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-12, "k_ub":  1e6, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-12, "k_ub":  1e7, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-12,  "k_ub": 1e8, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-12,  "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-12,  "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-12,  "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-12,  "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-12,  "k_ub": 1e13, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b": 1e-16,  "k_ub": 1e16, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-15,  "k_ub": 1e15, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-14,  "k_ub": 1e14, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
