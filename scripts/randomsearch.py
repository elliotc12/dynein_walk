#!/usr/bin/python3

import os
import dynein.run as run
import glob
import numpy as np
import argparse

import sys
sys.path.insert(0, os.getcwd() + "/data/")
import paper_params as params

if os.path.exists('parameter_search.py'):
    os.chdir('../')
os.system("make generate_stepping_data")

runtime = 1e-1
kb = 1e9
kub = 1e9

def get_rand_spring(minn, maxx):
    p = 0
    x = 0
    while (np.random.random() > p):
        x = np.random.random()*(maxx-minn) + minn
        p = 2*(maxx - x ) / (maxx - minn)**2
    return x

for _ in range(5):
    cm = get_rand_spring(0.01, 2)
    cb = get_rand_spring(0.01, 2)
    ct = get_rand_spring(0.01, 2)

    os.system("rq run --job-name randomsearch"\
              + " python3 scripts/generate-stepping-data.py" \
              + " --kub " + str(kub) \
              + " --kb " + str(kb) \
              + " --cb " + "{:.2f}".format(cb) \
              + " --cm " + "{:.2f}".format(cm) \
              + " --ct " + "{:.2f}".format(ct) \
              + " --ls " + str(params.ls) \
              + " --lt " + str(params.lt) \
              + " --eqb " + str(params.eqb) \
              + " --eqt " + str(params.eqt) \
              + " --eqmpre " + str(params.eqmpre) \
              + " --eqmpost " + str(params.eqmpost) \
              + " --rt " + str(params.radius_t) \
              + " --rm " + str(params.radius_m) \
              + " --rb " + str(params.radius_b) \
              + " --seed " + str(1) \
              + " --unbindingconst " + str(params.cexp)\
              + " --label randomsearch"\
              + " --runtime " + str(runtime))
