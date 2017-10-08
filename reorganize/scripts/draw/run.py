#/usr/bin/python
import subprocess, os
import numpy as np

def sim(**run):

    if 'label' in run:
      basename = "%s__k_b-%g,k_ub-%g,cb-%g,cm-%g,ct-%g,dt-%g" % (str(run["label"]), run["k_b"], run["k_ub"],
                                                                 run["cb"], run["cm"], run["ct"], run["dt"])
    else:
      basename = "k_b-%g,k_ub-%g,cb-%g,cm-%g,ct-%g,dt-%g" % (str(run["label"]), run["k_b"], run["k_ub"],
                                                             run["cb"], run["cm"], run["ct"], run["dt"])

    cmd = ["./generate_stepping_data"]

    for key in ["ls", "lt", "k_b", "k_ub", "cb", "cm", "ct", "T", "dt", "label", "seed", "runtime", "movie"]:
        if key in run:
            cmd.extend(["--"+key, str(run[key])])
    for key in ["nomovie", "onebound-debugging", "constant-write"]:
        if key in run:
            cmd.extend(["--"+key])
            
    print "Running: " + " ".join(cmd)
    subprocess.call(cmd)
    return basename
