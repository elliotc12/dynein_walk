#/usr/bin/python
import subprocess, os
import numpy as np

def read_csv(fname):
    values = np.loadtxt(fname, delimiter=',', comments='#', dtype='string')
    custom_runs = []
    if values.ndim != 1:
        for i in range(0, values.shape[0]):
            custom_runs.append({"ls":float(values[i,0]), "lt":float(values[i,1]),
                                "k_b":float(values[i,2]), "k_ub":float(values[i,3]),
                                "T":float(values[i,4]), "cb":float(values[i,5]),
                                "cm":float(values[i,6]), "ct":float(values[i,7])})
    else:
        custom_runs.append({"ls":float(values[0]), "lt":float(values[1]),
                            "k_b":float(values[2]), "k_ub":float(values[3]),
                            "T":float(values[4]), "cb":float(values[5]),
                            "cm":float(values[6]), "ct":float(values[7])})
    return custom_runs


# ask professor round / elliott about this latex_format stuff - 6/26/17
def latex_format(x):
    if isinstance(x, float) or isinstance(x, int):
        x = '{:g}'.format(x)
        if 'e+0' in x:
            m,e = x.split('e+0')
            if m == '1':
                return r'10^{'+e+'}'
            return m + r'\times 10^{' + e+ '}'
        if 'e+' in x:
            m,e = x.split('e+')
            if m == '1':
                return r'10^{'+e+'}'
            return m + r'\times 10^{' + e+ '}'
        if 'e-0' in x:
            m,e = x.split('e-0')
            if m == '1':
                return r'10^{-'+e+'}'
            return m + r'\times 10^{-' + e+ '}'
        if 'e' in x:
            m,e = x.split('e')
            if m == '1':
                return r'10^{'+e+'}'
            return m + r'\times 10^{' + e+ '}'
    if isinstance(x, str):
        x = x.replace('-', '_')
    return x

def run(**run):
    if 'label' in run:
      basename = "%s__k_b-%g,k_ub-%g,cb-%g,cm-%g,ct-%g,dt-%g" % (str(run["label"]), run["k_b"], run["k_ub"],
                                                                 run["cb"], run["cm"], run["ct"], run["dt"])
    else:
      basename = "k_b-%g,k_ub-%g,cb-%g,cm-%g,ct-%g,dt-%g" % (str(run["label"]), run["k_b"], run["k_ub"],
                                                             run["cb"], run["cm"], run["ct"], run["dt"])

    cmd = ["./generate_stepping_data"]

    for key in ["ls", "lt", "k_b", "k_ub", "cb", "cm", "ct", "T", "dt", "label", "seed", "runtime", "movie", "framerate"]:
        if key in run:
            cmd.extend(["--"+key, str(run[key])])
    for key in ["nomovie", "onebound-debugging", "constant-write"]:
        if key in run:
            cmd.extend(["--"+key])

    with open("data/stepping_parameters_%s.tex" % basename, "w") as f:
        for k,v in sorted(run.items()):
            if k == "label":
                f.write(r'\newcommand\runlabel{%s}' % (latex_format(v)) + '\n')
            else:
                f.write(r'\newcommand\%s{%s}' % (latex_format(k).replace("_",""), latex_format(v)) + '\n')
    out = open('runlogs/' + basename + '.out', 'w')
    process_object = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
    print("Running: ", " ".join(cmd))
    if "no-slurm" in run:
        process_object.wait()
    return basename
