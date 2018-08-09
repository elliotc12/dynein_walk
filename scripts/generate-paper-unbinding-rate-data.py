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

write_rate = 1e5
runtime = 1e-4
seed = 1

for L in [1, 5, 10, 15, 20, 25, 30, 35, 40]:
    basename = "paper_unbinding_probability__L-%s,s-%s" % (str(L), seed)

    cmd = ["./simulate_unbinding_rates",
           "--label", "%s" % "paper_unbinding_probability",
           "--k_b", "%g" % float(params.k_b),
           "--k_ub", "%g" % float(params.k_ub),
           "--c", "%g" % float(params.cexp),
           "--cb", "%g" % float(params.cb),
           "--cm", "%g" % float(params.cm),
           "--ct", "%g" % float(params.ct),
           "--ls", str(params.ls),
           "--lt", str(params.lt),
           "--eqb", str(params.eqb),
           "--eqmpre", str(params.eqmpre),
           "--eqmpost", str(params.eqmpost),
           "--eqt", str(params.eqt),
           "--write_rate", str(write_rate),
           "--runtime", str(runtime),
           "--seed", str(seed),
           "--dt", "1e-10",
           "--L", "%g" % float(L)]

    if not os.path.exists('runlogs'):
        os.makedirs('runlogs')
    out = open('runlogs/' + basename + '.out', 'w')

    print("Running: ", " ".join(cmd), out)
    out.flush()
    process_object = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE)
    err = process_object.communicate()[1]
    if (err != b''):
        print("\n##################################",
              "\nSimulation exited in error: \n\n",
              err.decode("utf-8"),
              "\n##################################\n\n")
