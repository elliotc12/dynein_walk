import subprocess, os

def have_slurm():
    try:
        subprocess.check_call("squeue > /dev/null", shell=True)
    except (OSError, subprocess.CalledProcessError):
        print("Not using slurm...")
        return False
    return True

def read_csv(fname):
    values = np.loadtxt(fname, delimiter=',', comments='#', dtype='string')
    custom_runs = []
    if values.ndim != 1:
        for i in range(0, values.shape[0]):
            custom_runs.append({"ls":float(values[i,0]), "lt":float(values[i,1]), "k_b":float(values[i,2]), "k_ub":float(values[i,3]),
                                "T":float(values[i,4]), "cb":float(values[i,5]), "cm":float(values[i,6]), "ct":float(values[i,7])})
        else:
            custom_runs.append({"ls":float(values[0]), "lt":float(values[1]), "k_b":float(values[2]), "k_ub":float(values[3]),
                                "T":float(values[4]), "cb":float(values[5]), "cm":float(values[6]), "ct":float(values[7])})
    return custom_runs


def run_sim(**run = {}, k_b=None):
    cmd = ["srun"] if have_slurm() else []
    cmd.extend(["nice", "-19"])
    cmd.extend(["./generate_stepping_data"])
    cmd.extend(["--runtime", str(runtime)])
    cmd.extend(["--label", str(custom_label)])

    for key in ["ls", "lt", "k_b", "k_ub", "cb", "cm", "ct", "T", "dt", "custom-label"]:
        if key in run:
            cmd.extend(["--"+key, str(run[key])])
    if k_b is not None:
        cmd.extend(['--k_b', k_b])

    if "movie" in run and run["movie"]:
        cmd.extend(["--movie"])
    if "onebound-debugging" in run and run["onebound-debugging"]:
        cmd.extend(["--onebound-debugging"])
    if "constant-write" in run and run["constant-write"]:
        cmd.extend(["--constant-write"])

    basename = '%s__k_b-%s,k_ub-%s,c-%s,dt-%s' % (str(custom_label), str(run["k_b"]), str(run["k_ub"]), str(run["cb"]), str(run["dt"]))
    out = open('runlogs/' + basename + '.out', 'w')
    subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
    print("Running: %s", " ".join(cmd))
