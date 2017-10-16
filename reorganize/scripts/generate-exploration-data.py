#!/usr/bin/python3
import os
import numpy as np
import dynein.run as run

l = "exploration"

basename  = run.run(**{"k_b": 1e15,
                       "k_ub": 1e15,
                       "cb": 0.4,
                       "cm": 1.3,
                       "ct": 0.5,
                       "dt": 1e-10, "label": l, "seed": 1, "runtime": 5e-3,
                       "framerate": 1e-8, "constant-write": True,
                       "no-slurm": True})

os.rename("data/stepping_movie_data_%s.txt" % (basename), "data/exploration_movie_data.txt")
os.rename("data/stepping_data_%s.txt" % (basename), "data/exploration_stepping_data.txt")
os.rename("data/stepping_parameters_%s.tex" % (basename), "data/exploration_stepping_parameters.tex")

os.unlink("data/stepping_movie_config_%s.txt" % basename)
os.unlink("data/stepping_config_%s.txt" % basename)
