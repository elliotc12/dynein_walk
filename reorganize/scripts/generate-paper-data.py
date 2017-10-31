#!/usr/bin/python3
import os, sys
import numpy as np
import dynein.run as run

runtime = 1e-2
l = "paper"
if 'long' in sys.argv:
    runtime = 100e-3
    l = 'long_paper'

basename  = run.sim(**{"k_b": 1e14,
                       "k_ub": 1e4,
                       "cb": 5.0,
                       "cm": 2.0,
                       "ct": 0.2,
                       "ls": 22.1,
                       "lt": 10.0,
                       "eqb": 116,
                       "eqmpre": 224,
                       "eqmpost": 160,
                       "eqt": 0,
                       "dt": 1e-10, "label": l, "seed": 1, "runtime": runtime,
                       "framerate": 1e-8, "crash-movie": False,
                       "no-slurm": True})

os.rename("data/stepping_movie_data_%s.txt" % (basename), "data/%s_movie_data.txt" % l)
os.rename("data/stepping_data_%s.txt" % (basename), "data/%s_stepping_data.txt" % l)
os.rename("data/stepping_parameters_%s.tex" % (basename), "data/%s_stepping_parameters.tex" % l)

os.unlink("data/stepping_movie_config_%s.txt" % basename)
os.unlink("data/stepping_config_%s.txt" % basename)
