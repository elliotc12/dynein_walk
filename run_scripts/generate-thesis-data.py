#!/usr/bin/python3
import simrunner, os
import numpy as np

l = "thesis"

simrunner.avoid_slurm = True

basename = simrunner.run_sim(**{"k_b": 1e18, "k_ub": 1e20, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-10, "label": l,
                                "seed": 8,
                                "runtime": 1e-3,
                                "constant-write": True,
                                "no-slurm": True})

#this one has a bothbound at the end
# basename = simrunner.run_sim(**{"k_b": 1e18, "k_ub": 1e7, "cb": 2.0, "cm": 2.0, "ct": 1.0, "dt": 1e-10, "label": l,
#                                 "seed": 10, # 8 was good
#                                 "runtime": 1e-4,
#                                 "constant-write": True,
#                                 "no-slurm": True})

# basename = simrunner.run_sim(**{"k_b": 1e18, "k_ub": 1e11, "cb": 1.5, "cm": 1.5, "ct": 1.0, "dt": 1e-10, "label": l,
#                                 "seed": 10, # 8 was good
#                                 "runtime": 1e-4,
#                                 "constant-write": True,
#                                 "no-slurm": True})

print(basename)

os.rename("../data/stepping_movie_data_%s.txt" % (basename), "../data/thesis_movie_data.txt")
os.rename("../data/stepping_data_%s.txt" % (basename), "../data/thesis_stepping_data.txt")
os.rename("../data/stepping_parameters_%s.tex" % (basename), "../data/thesis_stepping_parameters.tex")

os.unlink("../data/stepping_movie_config_%s.txt" % basename)
os.unlink("../data/stepping_config_%s.txt" % basename)
