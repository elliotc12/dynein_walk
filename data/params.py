#!/usr/bin/python3

for_simulation= {"k_b": 5.5e9, # 3.4e9, # previously 3e6
                 "k_ub": 30,
                 "k_stk": 1e6,  # NEW: Sticky Rate
                 "cb": 0.08,
                 "cm": 1.05,
                 "ct": 0.36,
                 "ls": 20.75,
                 "lt": 23.0,
                 "rt": 6.0,
                 "rm": 5.5,
                 "rb": 1.75,
                 "seed": 1,
                 "dt": 1e-13, # 1e-10,
                 "eqb" :120,
                 "eqmpre": 197,
                 "eqmpost": 242,
                 "eqt": 0,
                 "force": 0,
                 "exp-unbinding-constant": -0.3563,
                 "runtime": 0.1,
                 "crash-movie": False,
                 "nomovie": True,
                 "angle-logging-mode": False,
                 "long-angle-logging-mode": False,
                 "label": "test"
}
# other good combos of ct/cm/cb: 0.03/1.01/0.2 and 0.23/1.36/0.09
