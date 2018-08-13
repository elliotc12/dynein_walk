#!/usr/bin/python3
import argparse
import os, sys
import numpy as np
import dynein.run as run

sys.path.insert(0, os.getcwd() + "/data/")
import paper_params as params

runtime = 5
seeds = [1, 2, 3, 4, 5, 6, 7, 8]

for sim in sims:
    for s in seeds:
        os.system("rq run --job-name paperexponentialhisto-" + str(s) \
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
                  + " --unbindingconst " + str(params.cexp)\
                  + " --label paperexponentialhisto-" + str(s)\
                  + " --renameexponential"\
                  + " --runtime " + str(runtime))
