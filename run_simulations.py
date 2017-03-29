#!/usr/bin/python2.7
import numpy as np
import subprocess, os

os.system('mkdir -p runlogs data')

assert(subprocess.call("make histogram-stuff", shell=True) == 0)

have_slurm = True

try:
    subprocess.check_call("squeue > /dev/null", shell=True)
except (OSError, subprocess.CalledProcessError):
    print "Not using slurm..."
    have_slurm = False

atp_in_kJ_per_mol = 30.5

binding_energy_high_affinity_kJ_mol = 71;
binding_energy_high_affinity_atp = binding_energy_high_affinity_kJ_mol / atp_in_kJ_per_mol;

custom_runs = []

custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 1e-11})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 3e-11})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 6e-11})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 8e-11})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 1e-10})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 2e-10})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 3e-10})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 5e-10})
custom_runs.append({"ls": 22.1, "lt": 11.15, "k_b": 10, "k_ub": 2e11,
                    "cb": 2.4,
                    "cm": 2.4,
                    "ct": 2.4,
                    "dt": 2e-11})

ls_min = 22.1 # nm
ls_max = 22.1 # nm
ls_num = 1

lt_min = 11.15 # nm
lt_max = 11.15 # nm
lt_num = 1

k_b_min = 5000 # s^-1
k_b_max = 10000 # s^-1
k_b_num = 1

cb_min = .1*binding_energy_high_affinity_atp # s^-1
cb_max = .1*binding_energy_high_affinity_atp # s^-1
cb_num = 1

cm_min = 0.01*binding_energy_high_affinity_atp # s^-1
cm_max = 0.5*binding_energy_high_affinity_atp # s^-1
cm_num = 10

ct_min = .1*binding_energy_high_affinity_atp # s^-1
ct_max = .1*binding_energy_high_affinity_atp # s^-1
ct_num = 1

T_min = 310.15 # K
T_max = 310.15 # K
T_num = 1

ls_range = np.linspace(ls_min, ls_max, num=ls_num)
lt_range = np.linspace(lt_min, lt_max, num=lt_num)
k_b_range = np.linspace(k_b_min, k_b_max, num=k_b_num)
T_range = np.linspace(T_min, T_max, num=T_num)
cb_range = np.linspace(cb_min, cb_max, num=cb_num)
cm_range = np.linspace(cm_min, cm_max, num=cm_num)
ct_range = np.linspace(ct_min, ct_max, num=ct_num)

runtime = 0.3
label = "test-nan-corrections"

if len(custom_runs) != 0:
    for run in custom_runs:
        custom_label = label+str(run["dt"]) #edit this to uniquely name each simulation

        cmd = ["srun"] if have_slurm else []
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
