#!/usr/bin/python2.7
from __future__ import print_function
import numpy as np
import random
import subprocess, os
import sys

os.system('mkdir -p runlogs data')

assert(subprocess.call("make histogram-stuff", shell=True) == 0)

have_slurm = True

try:
    subprocess.check_call("squeue > /dev/null", shell=True)
except (OSError, subprocess.CalledProcessError):
    print("Not using slurm...")
    have_slurm = False

atp_in_kJ_per_mol = 30.5

binding_energy_high_affinity_kJ_mol = 71;
binding_energy_high_affinity_atp = binding_energy_high_affinity_kJ_mol / atp_in_kJ_per_mol;

random.seed() # seed from a random source
custom_runs = []

custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 1e-11,
                    "seed": random.randint(0,100)})
custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 3e-11,
                    "seed": random.randint(0,100)})
custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 6e-11,
                    "seed": random.randint(0,100)})
custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 8e-11,
                    "seed": random.randint(0,100)})
custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 1e-10,
                    "seed": random.randint(0,100)})
custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 2e-10,
                    "seed": random.randint(0,100)})
custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 3e-10,
                    "seed": random.randint(0,100)})
custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 5e-10,
                    "seed": random.randint(0,100)})
custom_runs.append({"k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 2e-11,
                    "seed": random.randint(0,100)})

if len(custom_runs) == 0:  # if no custom_runs specified above, load them from a file
        cmdArgs = sys.argv[1:]
        file = cmdArgs[0]
        try:
                values = np.loadtxt(file, delimiter=',', comments='#', dtype='string')
                if values.ndim != 1:
                        for i in range(0, values.shape[0]):
                                custom_runs.append({"ls":float(values[i,0]), "lt":float(values[i,1]), "k_b":float(values[i,2]), "k_ub":float(values[i,3]),
                                                    "T":float(values[i,4]), "cb":float(values[i,5]), "cm":float(values[i,6]), "ct":float(values[i,7])})
                else:
                        custom_runs.append({"ls":float(values[0]), "lt":float(values[1]), "k_b":float(values[2]), "k_ub":float(values[3]),
                                            "T":float(values[4]), "cb":float(values[5]), "cm":float(values[6]), "ct":float(values[7])})
        except:
                print("usage: %s filename", sys.argv[0])

runtime = 0.3
label = "test-nan-corrections"

for run in custom_runs:
        custom_label = label+str(run["dt"]) #edit this to uniquely name each simulation

        # cmd = ["srun"] if have_slurm else []
        cmd = []
        cmd.extend(["echo"])
        cmd.extend(["nice", "-19"])
        cmd.extend(["./generate_stepping_data"])
        cmd.extend(["--runtime", str(runtime)])
        cmd.extend(["--label", str(custom_label)])

        for key in ["ls", "lt", "k_b", "k_ub", "cb", "cm", "ct", "T", "dt", "custom-label"]:
            if key in run:
                cmd.extend(["--"+key, str(run[key])])

        if "movie" in run and run["movie"]:
            cmd.extend(["--movie"])
        if "onebound-debugging" in run and run["onebound-debugging"]:
            cmd.extend(["--onebound-debugging"])
        if "constant-write" in run and run["constant-write"]:
            cmd.extend(["--constant-write"])

        basename = '%s__k_b-%s,k_ub-%s,c-%s,dt-%s' % (str(custom_label), str(run["k_b"]), str(run["k_ub"]), str(run["cb"]), str(run["dt"]))
        out = open('runlogs/' + basename + '.out', 'w')
        subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
        print("Running: %s" % " ".join(cmd))
