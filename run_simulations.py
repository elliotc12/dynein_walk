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

# --- adding code to pull parameters from cmd line specified txt file - John - 3/27 --- #
# --- the code will now pull the values from runVars.txt file specified at terminal
custom_runs = []
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
	print "incorrect format --- could not gather simulation variables"
print "\n"

custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 0.1, "k_ub": 2e11, "T": 310.15,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "movie": False})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 1, "k_ub": 2e11, "T": 310.15,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "movie": False})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 2, "k_ub": 2e11, "T": 310.15,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "movie": False})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 4, "k_ub": 2e11, "T": 310.15,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "movie": False})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 6, "k_ub": 2e11, "T": 310.15,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "movie": False})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 8, "k_ub": 2e11, "T": 310.15,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "movie": False})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11, "T": 310.15,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "movie": False})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 15, "k_ub": 2e11, "T": 310.15,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "movie": False})


ls_min = 22.1 # nm
ls_max = 22.1 # nm
ls_num = 1

lt_min = 11.15 # nm
lt_max = 11.15 # nm
lt_num = 1

k_b_min = 5000 # s^-1
k_b_max = 10000 # s^-1
k_b_num = 1
>>>>>>> 498419f611caf31ae6f13616b613442642752bc0

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

<<<<<<< HEAD
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
=======
cm_min = 0.01*binding_energy_high_affinity_atp # s^-1
cm_max = 0.5*binding_energy_high_affinity_atp # s^-1
cm_num = 10

ct_min = .1*binding_energy_high_affinity_atp # s^-1
ct_max = .1*binding_energy_high_affinity_atp # s^-1
ct_num = 1

T_min = 310.15 # K
T_max = 310.15 # K
T_num = 1

label = "debug-unbinding-3-13"

ls_range = np.linspace(ls_min, ls_max, num=ls_num)
lt_range = np.linspace(lt_min, lt_max, num=lt_num)
k_b_range = np.linspace(k_b_min, k_b_max, num=k_b_num)
T_range = np.linspace(T_min, T_max, num=T_num)
cb_range = np.linspace(cb_min, cb_max, num=cb_num)
cm_range = np.linspace(cm_min, cm_max, num=cm_num)
ct_range = np.linspace(ct_min, ct_max, num=ct_num)

runtime = 0.33

if len(custom_runs) != 0:
    for run in custom_runs:
        cmd = ["srun"] if have_slurm else []
        cmd.extend([
            "nice", "-19",
            "./generate_stepping_data",
            "--Ls",  str(run["ls"]),
            "--Lt",  str(run["lt"]),
            "--k_b", str(run["k_b"]),
            "--k_ub",str(run["k_ub"]),
            "--cb",  str(run["cb"]),
            "--cm",  str(run["cm"]),
            "--ct",  str(run["ct"]),
            "--T",   str(run["T"]),
            "--label", label,
            "--runtime", str(runtime),
        ])
        if (run["movie"]):
            cmd.extend(["--movie"])
        if ("onebound-debugging" in run and run["onebound-debugging"]):
            cmd.extend(["--onebound-debugging"])

        basename = '%s__ls-%.3g,lt-%.3g,k_b-%s,k_ub-%s,cb-%s,cm-%s,ct-%s,T-%s' % (label, run['ls'], run['lt'], run["k_b"], run["k_ub"], run["cb"], run["cm"], run["ct"], run['T'])
        out = open('runlogs/' + basename + '.out', 'w')
        subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
        print "Running: ", ' '.join(cmd)
else:
    for permutation in [{"ls": ls,"lt": lt,"k_b": k_b, "T": T, "cb": cb, "cm": cm, "ct": ct} for ls in ls_range for lt in lt_range
                    for k_b in k_b_range for T in T_range for cb in cb_range for cm in cm_range for ct in ct_range]:
        cmd = ["srun"] if have_slurm else []
        cmd.extend([
            "./generate_stepping_data",
            "--Ls", str(permutation["ls"]),
            "--Lt", str(permutation["lt"]),
            "--k_b", str(permutation["k_b"]),
            "--cb", str(permutation["cb"]),
            "--cm", str(permutation["cm"]),
            "--ct", str(permutation["ct"]),
            "--T", str(permutation["T"]),
            "--label", label,
            "--runtime", str(runtime),
            "--movie"
        ])
        print "Running: ", ' '.join(cmd)

        basename = '%s__ls-%.3g,lt-%.3g,k_b-%s,cb-%s,cm-%s,ct-%s,T-%s' % (label, permutation['ls'], permutation['lt'], permutation["k_b"],
                                                                      permutation["cb"], permutation["cm"], permutation["ct"], permutation['T'])
        out = open('runlogs/' + basename + '.out', 'w')
        subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
>>>>>>> 498419f611caf31ae6f13616b613442642752bc0
