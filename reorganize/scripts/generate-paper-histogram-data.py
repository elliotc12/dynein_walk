#!/usr/bin/python3
import argparse
import os, sys
import numpy as np
import dynein.run as run

runtime = 2

os.system("make generate_stepping_data")
os.system("rm data/stepping_data_paperhisto_*")

# parser = argparse.ArgumentParser()
# parser.add_argument("-s", "--seed", type=int, help="Manually set the seed value, default=1")
# args = parser.parse_args()
# if args.seed:
#     seed = args.seed
# else:
#     seed = 1

label = "paperhisto"
seeds = [1, 2, 3, 4, 5]

params = {"k_b": 1e9,
          "k_ub": 100,
          "cb": 0.1,
          "cm": 0.5,
          "ct": 0.2,
          "ls": 10.49, # from urnavicius 2015 (paper.bib)
          "lt": 23.8,  # from urnavicius 2015
          "eqb": 120,  # from redwine 2012 supplemental
          "eqmpre": 200, # from burgess 2002, 360-160
          "eqmpost": 224, # from burgess 2002, 360-136
          "eqt": 0,
          "dt": 1e-10, "label": label, "seed": seed, "runtime": runtime,
          "framerate": 1e-5, "crash-movie": False,
          "no-slurm": True}

for s in seeds:
    params["seed"] = s
    run.sim(**params)

os.system("rm data/stepping_movie_config_paperhisto*.txt" % basename)
os.system("rm data/stepping_config_paperhisto*.txt" % basename)
