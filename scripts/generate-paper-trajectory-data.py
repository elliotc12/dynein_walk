#!/usr/bin/python3
import os, sys
import numpy as np
import dynein.run as run

runtime = 1e-5
framerate = 1e-10

os.system("python3 scripts/generate-stepping-data.py" \
          + " --kub 1e9"\
          + " --kb  1e25"\
          + " --cb  0.1"\
          + " --cm  0.5"\
          + " --ct  0.2"\
          + " --runtime " + str(runtime)\
          + " --framerate " + str(framerate)\
          + " --movie"\
          + " --label papertrajecto"\
          + " --renametrajectory")
