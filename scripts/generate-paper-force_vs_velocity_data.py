#!/usr/bin/python3

import os
import dynein.run as run
import glob
import numpy as np
import argparse
import subprocess
import sys

sys.path.insert(0, os.getcwd() + "/data/")
import paper_params as params

runtime = 2

for f in [4, 0, -1, -2, -3, -4, -5, -6, -7, -8, -10]:
    for s in [1, 2]:
        basename = "paper_force__F-%s,s-%s" % (str(f), str(s))
        os.system("rq run --job-name paperforce-F-{},s-{}".format(f, s) \
                  + " python3 scripts/generate-stepping-data.py" \
                  + " --kub " + str(params.k_ub) \
                  + " --kb " + str(params.k_b) \
                  + " --cb " + str(params.cb) \
                  + " --cm " + str(params.cm) \
                  + " --ct " + str(params.ct) \
                  + " --ls " + str(params.ls) \
                  + " --lt " + str(params.lt) \
                  + " --eqb " + str(params.eqb) \
                  + " --eqt " + str(params.eqt) \
                  + " --eqmpre " + str(params.eqmpre) \
                  + " --eqmpost " + str(params.eqmpost) \
                  + " --rt " + str(params.radius_t) \
                  + " --rm " + str(params.radius_m) \
                  + " --rb " + str(params.radius_b) \
                  + " --seed " + str(s) \
                  + " --force {}".format(f) \
                  + " --unbindingconst " + str(params.cexp)\
                  + " --label paperforce_F-{},s-{}".format(f, s)\
                  + " --renameforce"\
                  + " --runtime {}".format(runtime))
