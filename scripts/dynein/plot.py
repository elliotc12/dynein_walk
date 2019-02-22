"""Functions to make plots for dynein simulation data
"""



import numpy as np
import matplotlib.pyplot as plt
import sys, os
import importlib
import glob

def find_data_file(paramsfile):
    # make sure file exists
    assert(os.path.isfile(paramsfile))

    # find end of path
    path = os.path.split(paramsfile)
    sys.path.append(path[0])
    params = importlib.import_module(path[1][:-3]) #make sure to chop off .py for import

    # create name for file that matches simulation output.
    if 'label' in params.for_simulation:
        search_string = "%s__k_b-%g,k_ub-%g,c-%g,cb-%g,cm-%g,ct-%g,ls-%g,lt-%g,seed-%g,dt-%g" % \
                (str(params.for_simulation["label"]), params.for_simulation["k_b"], params.for_simulation["k_ub"],
                 params.for_simulation["exp-unbinding-constant"], params.for_simulation["cb"], params.for_simulation["cm"],
                 params.for_simulation["ct"], params.for_simulation["ls"], params.for_simulation["lt"],
                 params.for_simulation["seed"], params.for_simulation["dt"])
    else:
         search_string = "k_ub-%g,c-%g,cb-%g,cm-%g,ct-%g,ls-%g,lt-%g,seed-%g,dt-%g" % \
                (params.for_simulation["k_b"], params.for_simulation["k_ub"], params.for_simulation["exp-unbinding-constant"],
                 params.for_simulation["cb"], params.for_simulation["cm"], params.for_simulation["ct"],
                 params.for_simulation["ls"], params.for_simulation["lt"], params.for_simulation["seed"],
                 params.for_simulation["dt"])

    data_string = os.path.dirname(paramsfile)+'/*'+search_string+'.txt'
    print(data_string)
    data_file = glob.glob(data_string)
    print(data_file)
    return data_file

if __name__ == "__main__":
    find_data_file("../../data/params.py")

