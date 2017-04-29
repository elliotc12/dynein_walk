#!/usr/bin/python
import simrunner, os
import numpy as np

l = "thesis"

simrunner.avoid_slurm = True

basename = simrunner.run_sim(**{"k_b": 1e5, "k_ub":  1e6, "cb": 2.4, "cm": 2.4, "ct": 2.4, "dt": 1e-10, "label": l,
                                "runtime": 1e-3,
                                "constant-write": True})

print basename

os.rename("../data/stepping_movie_data_%s.txt" % (basename),
          "../data/thesis_movie.txt")
os.rename("../data/stepping_data_%s.txt" % (basename),
          "../data/thesis_data.txt")
os.unlink("../data/stepping_movie_config_%s.txt" % basename)
os.unlink("../data/stepping_config_%s.txt" % basename)
