#!/usr/bin/python3
import os
import numpy as np
import dynein.run as run

l = "paper"

basename  = run.sim(**{"k_b": 1e14,
                       "k_ub": 1e5,
                       "cb": 2.0,
                       "cm": 2.0,
                       "ct": 0.2,
                       "ls": 22.1,
                       "lt": 10.0,
                       "eqb": 116,
                       "eqmpre": 224,
                       "eqmpost": 160,
                       "eqt": 0,
                       "dt": 1e-10, "label": l, "seed": 1, "runtime": 1e-3,
                       "framerate": 1e-8, "constant-write": True,
                       "no-slurm": True})

os.rename("data/stepping_movie_data_%s.txt" % (basename), "data/paper_movie_data.txt")
os.rename("data/stepping_data_%s.txt" % (basename), "data/paper_stepping_data.txt")
os.rename("data/stepping_parameters_%s.tex" % (basename), "data/paper_stepping_parameters.tex")

os.unlink("data/stepping_movie_config_%s.txt" % basename)
os.unlink("data/stepping_config_%s.txt" % basename)
