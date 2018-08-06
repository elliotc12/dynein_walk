#!/usr/bin/python3
import os, sys
import numpy as np
import dynein.run as run

sys.path.insert(0, os.getcwd() + "/data/")
import paper_params as params

runtime = 1e-5
framerate = 1e-10

os.system("python3 scripts/generate-stepping-data.py" \
          + " --ls " + str(params.ls)\
          + " --lt " + str(params.lt)\
          + " --kub " + str(params.trajectory_k_ub)\
          + " --kb  " + str(params.trajectory_k_b)\
          + " --cb " + str(params.cb)\
          + " --cm " + str(params.cm)\
          + " --ct " + str(params.ct)\
          + " --eqb " + str(params.eqb)\
          + " --eqt " + str(params.eqt)\
          + " --eqmpre " + str(params.eqmpre)\
          + " --eqmpost " + str(params.eqmpost)\
          + " --runtime " + str(runtime)\
          + " --framerate " + str(framerate)\
          + " --movie"\
          + " --label papertrajecto"\
          + " --renametrajectory")
