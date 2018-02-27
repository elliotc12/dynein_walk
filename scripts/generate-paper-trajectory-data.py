#!/usr/bin/python3
import os, sys
import numpy as np
import dynein.run as run

runtime = 1e-4
framerate = 1e-7

if (runtime / framerate > 1e5):
    print("Error: runtime/framerate > 1e5; this would result in a file larger than 1e5 lines. This is too big for git; please shorten.")
    exit(1)

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
                       "dt": 1e-10, "label": "papertrajecto", "seed": 1, "runtime": runtime,
                       "framerate": framerate, "crash-movie": False, "nomovie": False})

os.rename("data/stepping_movie_data_%s.txt" % (basename), "data/paper_trajectory_movie_data.txt")
os.rename("data/stepping_data_%s.txt" % (basename), "data/paper_trajectory_stepping_data.txt")
os.rename("data/stepping_parameters_%s.tex" % (basename), "data/paper_trajectory_stepping_parameters.tex")

os.unlink("data/stepping_movie_config_%s.txt" % basename)
os.unlink("data/stepping_config_%s.txt" % basename)
