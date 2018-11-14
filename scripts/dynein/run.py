#/usr/bin/python
import subprocess, os
import numpy as np

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

def sim(**run):
    if 'label' in run:
      basename = "%s__k_b-%g,k_ub-%g,c-%g,cb-%g,cm-%g,ct-%g,ls-%g,lt-%g,seed-%g,dt-%g" % \
                 (str(run["label"]), run["k_b"], run["k_ub"], run["exp-unbinding-constant"], run["cb"], run["cm"],
                  run["ct"], run["ls"], run["lt"], run["seed"], run["dt"])
    else:
      basename = "k_b-%g,k_ub-%g,c-%g,cb-%g,cm-%g,ct-%g,dt-%g" % (str(run["label"]), run["k_b"], run["k_ub"],
                                                             run["exp-unbinding-constant"], run["cb"], run["cm"], run["ct"], run["dt"])

    cmd = ["./generate_stepping_data"]

    for key in run:
        if key in ["nomovie", "onebound-debugging", "crash-movie", "full-gibbs-transitions", "angle-logging-mode"]:
            if run[key] == True:
                cmd.extend(["--"+key])
        else:
            cmd.extend(["--"+key, str(run[key])])

    if not os.path.exists('data'):
        os.makedirs('data')

    with open("data/stepping_parameters_%s.tex" % basename, "w") as f:
        for k,v in sorted(run.items()):
            if k == "label":
                f.write(r'\newcommand\runlabel{%s}' % (latex_format(v)) + '\n')
            else:
                f.write(r'\newcommand\%s{%s}' % (latex_format(k).replace("_",""), latex_format(v)) + '\n')

    #os.makedirs('runlogs', exist_ok=True) # ensure runlogs directory exists
    # if not os.path.exists('runlogs'):
    #     os.makedirs('runlogs')
    # out = open('runlogs/' + basename + '.out', 'w')
    # print("Running: ", " ".join(cmd), out)
    # out.flush()
    # process_object = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE)
    process_object = subprocess.Popen(cmd, stderr=subprocess.PIPE)
    err = process_object.communicate()[1]

    if (err != b''):
        print("\n##################################",
              "\nSimulation exited in error: \n\n",
              err.decode("utf-8"),
              "\n##################################\n\n")
    
    return basename
