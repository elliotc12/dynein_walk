#!/usr/bin/python3
import argparse
import os, sys
import numpy as np
import dynein.run as run

runtime = 1e-4
l = "generic"

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", type=int, help="Manually set the seed value, default=1")
args = parser.parse_args()
if args.seed:
    seed = args.seed
else:
    seed = 1

basename  = run.sim(**{"k_b": 1e9,
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
                       "dt": 1e-10, "label": l, "seed": seed, "runtime": runtime,
                       "framerate": 1e-5, "crash-movie": False,
                       "no-slurm": True})

os.unlink("data/stepping_movie_config_%s.txt" % basename)
os.unlink("data/stepping_config_%s.txt" % basename)
