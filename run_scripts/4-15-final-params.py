#!/usr/bin/python
import random
import simrunner
import numpy as np

l = "4-21"

simrunner.run_sim(**{"k_b": 1e-4, "k_ub":  1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e-4, "k_ub":  1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e-4, "k_ub":  1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e-4, "k_ub":  1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})

simrunner.run_sim(**{"k_b": 1e-2, "k_ub":  1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e-2, "k_ub":  1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e-2, "k_ub":  1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e-2, "k_ub":  1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})

simrunner.run_sim(**{"k_b": 1e0, "k_ub":  1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e0, "k_ub":  1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e0, "k_ub":  1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e0, "k_ub":  1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})

simrunner.run_sim(**{"k_b": 1e2, "k_ub":  1e9, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e2, "k_ub":  1e10, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e2, "k_ub":  1e11, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
simrunner.run_sim(**{"k_b": 1e2, "k_ub":  1e12, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l, "seed": random.randint(0,100)})
