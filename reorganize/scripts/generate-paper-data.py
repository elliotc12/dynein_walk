#!/usr/bin/python3
import os, sys
import numpy as np
import dynein.run as run

runtime = 1e-3
l = "paper"
if 'long' in sys.argv:
    runtime = 100e-3
    l = 'long_paper'

basename  = run.sim(**{"k_b": 1e15,
                       "k_ub": 1e7,
                       "cb": 0.5,
                       "cm": 0.7,
                       "ct": 0.2,
                       "ls": 10.49, # from urnavicius 2015 (paper.bib)
                       "lt": 23.8,  # from urnavicius 2015
                       "eqb": 120,  # from redwine 2012 supplemental
                       "eqmpre": 200, # from burgess 2002, 360-160
                       "eqmpost": 224, # from burgess 2002, 360-136
                       "eqt": 0,
                       "dt": 1e-10, "label": l, "seed": 1, "runtime": runtime,
                       "framerate": 1e-8, "crash-movie": False,
                       "no-slurm": True})

os.rename("data/stepping_movie_data_%s.txt" % (basename), "data/%s_movie_data.txt" % l)
os.rename("data/stepping_data_%s.txt" % (basename), "data/%s_stepping_data.txt" % l)
os.rename("data/stepping_parameters_%s.tex" % (basename), "data/%s_stepping_parameters.tex" % l)

os.unlink("data/stepping_movie_config_%s.txt" % basename)
os.unlink("data/stepping_config_%s.txt" % basename)
