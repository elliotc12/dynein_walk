#!/usr/bin/python
import simrunner
import numpy as np

l = "4-23"

simrunner.run_sim(**{"k_b": 1e-5, "k_ub":  1e2, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-5, "k_ub":  1e3, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-5, "k_ub":  1e4, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-5, "k_ub":  1e5, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-5, "k_ub":  1e6, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b": 1e-6, "k_ub":  1e2, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6, "k_ub":  1e3, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-6, "k_ub":  1e4, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b": 1e-7, "k_ub":  1e2, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-7, "k_ub":  1e3, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-7, "k_ub":  1e4, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-7, "k_ub":  1e5, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-7, "k_ub":  1e6, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
