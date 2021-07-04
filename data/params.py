#!/usr/bin/python3

for_simulation= {"k_b": 5.5e10, # 3.4e9, previously 3e6
                 "k_ub": 2.5,
                 "k_stk": 1e8,  # NEW: Sticky Rate
                 "cb": 0.0, # OLD: 0.08,
                 "cm": 1.0, # OLD: 1.05,
                 "ct": 1.0, # OLD: 0.36,
                 "ls": 20.75,
                 "lt": 23.0,
                 "rt": 6.0,
                 "rm": 5.5,
                 "rb": 1.75,
                 "seed": 1,
                 "dt": 1e-13, # 1e-10,
                 "eqb" :120.0,
                 "eqmpre": 197.0,
                 "eqmpost": 242.0,
                 "eqt": 0,
                 "force": 0,
                 "exp-unbinding-constant": -0.35, # -0.5,
                 "runtime": 0.1,
                 "crash-movie": False,
                 "nomovie": True,
                 "angle-logging-mode": False,
                 "long-angle-logging-mode": False,
                 "label": "test",
                 "boltzmann-constant": 2.726097017e-4, # kB in ATP energies per K (from default_parameters.h)
                 "T": 310.15  # temperature in K
}
# other good combos of ct/cm/cb: 0.03/1.01/0.2 and 0.23/1.36/0.09
