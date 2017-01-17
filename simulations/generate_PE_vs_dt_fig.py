#!/usr/bin/python2.7

import numpy as np
import subprocess

assert(subprocess.call("make histogram-stuff", shell=True) == 0)

have_slurm = True

try:
    subprocess.check_call("squeue > /dev/null", shell=True)
except (OSError, subprocess.CalledProcessError):
    print "Not using slurm..."
    have_slurm = False

custom_runs = []

custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 5000, "k_ub": 50000, "T": 310.15,
                    "cb": 0.8,
                    "cm": 0.8,
                    "ct": 0.8,
                    "dt": 1e-13,
                    "movie": True,
                    "constant-write": True})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 5000, "k_ub": 50000, "T": 310.15,
                    "cb": 0.8,
                    "cm": 0.8,
                    "ct": 0.8,
                    "dt": 1e-12,
                    "movie": True,
                    "constant-write": True})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 5000, "k_ub": 50000, "T": 310.15,
                    "cb": 0.8,
                    "cm": 0.8,
                    "ct": 0.8,
                    "dt": 1e-11,
                    "movie": True,
                    "constant-write": True})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 5000, "k_ub": 50000, "T": 310.15,
                    "cb": 0.8,
                    "cm": 0.8,
                    "ct": 0.8,
                    "dt": 1e-10,
                    "movie": True,
                    "constant-write": True})

label = "make-PE-dt-figure-"

runtime = 1e-8

for run in custom_runs:
    cmd = ["srun"] if have_slurm else []
    cmd.extend([
        "./generate_stepping_data",
        "--Ls",  str(run["ls"]),
        "--Lt",  str(run["lt"]),
        "--k_b", str(run["k_b"]),
        "--k_ub",str(run["k_ub"]),
        "--cb",  str(run["cb"]),
        "--cm",  str(run["cm"]),
        "--ct",  str(run["ct"]),
        "--T",   str(run["T"]),
        "--label", label + str(run['dt']),
        "--runtime", str(runtime),
    ])
    if (run["movie"]):
        cmd.extend(["--movie"])
    if (run["constant-write"]):
        cmd.extend(["--constant-write"])

    basename = '%s__ls-%.3g,lt-%.3g,k_b-%s,k_ub-%s,cb-%s,cm-%s,ct-%s,T-%s' % (label + str(run['dt']), run['ls'], run['lt'], run["k_b"], run["k_ub"], run["cb"], run["cm"], run["ct"], run['T'])

    out = open('runlogs/' + basename + '.out', 'w') # eventually replace with /dev/null
    subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
    print "Running: ", ' '.join(cmd)

PEs = []
dts = []

for run in custom_runs:
    basename = '%s__ls-%.4g,lt-%.4g,k_b-%s,k_ub-%s,cb-%s,cm-%s,ct-%s,T-%s' % (label + str(run['dt']), run['ls'], run['lt'], run["k_b"], run["k_ub"], run["cb"], run["cm"], run["ct"], run['T'])
    data_file = 'data/stepping_movie_data_' + basename + '.txt'

    file_busy = True
    while file_busy:
        try:
            fd = open(data_file, 'a')
            fd.close()
            file_busy = False
        except IOError, message:
            print "file still busy"
            
    data = np.genfromtxt(data_file, delimiter="\t", invalid_raise=False)
    PE = np.average(data[:,4])
    PEs.extend([PE])
    dts.extend([run['dt']])

print PEs
print dts
