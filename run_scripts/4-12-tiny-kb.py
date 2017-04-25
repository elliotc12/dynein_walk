#!/usr/bin/python
import simrunner
import numpy as np

l = "4-12-tiny-kb"

simrunner.run_sim(**{"k_b": 1e-3, "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-2, "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-1, "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e0,  "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e1,  "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 5e1,  "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e2,  "k_ub": 1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b": 1e-3, "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-2, "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-1, "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e0,  "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e1,  "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 5e1,  "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e2,  "k_ub": 1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b": 1e-3, "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-2, "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-1, "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e0,  "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e1,  "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 5e1,  "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e2,  "k_ub": 1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})

simrunner.run_sim(**{"k_b": 1e-3, "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-2, "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e-1, "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e0,  "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e1,  "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 5e1,  "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
simrunner.run_sim(**{"k_b": 1e2,  "k_ub": 1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l})
