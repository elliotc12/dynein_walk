#/usr/bin/python
import subprocess, os
import numpy as np

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
        if key in ["nomovie", "onebound-debugging", "crash-movie", "full-gibbs-transitions", "angle-logging-mode", "long-angle-logging-mode"]:
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

    # add to simulation to runlogs
    os.makedirs('runlogs', exist_ok=True) # ensure runlogs directory exists
    out = open('runlogs/' + basename + '.out', 'w')
    print("Running: ", " ".join(cmd))
    out.flush()

    # run the simulation
    rc = subprocess.run(cmd, stdout=out, stderr=out).returncode
    out.flush()

    if (rc != 0):
        print("\n##################################",
              "\nSimulation exited in error! See {}".format('runlogs/' + basename + '.out'),
              "\n##################################\n\n")

    return basename


def sim2(path_to_best_params=None, verbose=False, **run):
    # detect operating system using strategy found at
    # https://stackoverflow.com/questions/1854/python-what-os-am-i-running-on
    import platform
    if platform.system() == 'Darwin':
        os.system("make generate_stepping_data")
    elif platform.system() == "Linux":
        if verbose: print("You are running linux. Fac should handle everything...")
    else:
        if verbose: print("Your os is not recognized... Beware!")

    # we know now that code should be built. First, try using best_params from data folder
    # we can do this recursively
    if path_to_best_params!=None:
        # make sure file exists
        assert(os.path.isfile(path_to_best_params))

        # find end of path
        path = os.path.split(path_to_best_params)
        import sys
        sys.path.append(path[0])
        import importlib
        params = importlib.import_module(path[1][:-3]) #make sure to chop off .py for import
        sim2(**params.for_simulation)

    else:


        # to run the command, find the distance to root of directory where simulation executable lives
        current_dir = os.getcwd()
        folders = current_dir.split('/')
        root_index = folders.index('dynein_walk')
        distance_to_root = int(len(folders)-(root_index+1))

        cmd = [distance_to_root*"../"+"generate_stepping_data"]

        for key in run:
            if key in ["nomovie", "onebound-debugging", "crash-movie", "full-gibbs-transitions", "angle-logging-mode", "long-angle-logging-mode"]:
                if run[key] == True:
                    cmd.extend(["--"+key])
            else:
                cmd.extend(["--"+key, str(run[key])])


        # run the simulation
        print(cmd)
        rc = subprocess.run(cmd).returncode





if __name__ == "__main__":
    sim2("../../data/params.py", verbose=True)






































