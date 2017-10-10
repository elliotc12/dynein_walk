#!/usr/bin/python3
import simrunner, os
import numpy as np

l = "exploration"

simrunner.avoid_roundqueue = True

basename = simrunner.run_sim(**{"k_b": 35*1e9, "k_ub": 5e9, "cb": 2, "cm": 2, "ct": 1, "dt": 1e-10, "label": l,
                                "seed": 1,
                                "runtime": 1e-2,
                                "framerate": 1e-8,
                                "constant-write": True,
                                "no-slurm": True})

os.rename("data/stepping_movie_data_%s.txt" % (basename), "data/exploration_movie_data.txt")
os.rename("data/stepping_data_%s.txt" % (basename), "data/exploration_stepping_data.txt")
os.rename("data/stepping_parameters_%s.tex" % (basename), "data/exploration_stepping_parameters.tex")

os.unlink("data/stepping_movie_config_%s.txt" % basename)
os.unlink("data/stepping_config_%s.txt" % basename)
