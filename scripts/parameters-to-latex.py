#! /usr/bin/env python3

import numpy as np
import sys, os, importlib

sys.path.insert(0, os.getcwd() + "/data/")
params = importlib.import_module("params")

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
    return x

parameters = {
"k_b" :               params.for_simulation['k_b']*0.01, #kb in nm/s
"k_ub" :              params.for_simulation['k_ub'],
"k_stk" :             params.for_simulation['k_stk'],
"trajectory_k_b" :    9.5e9, # params.trajectory_k_b,
"trajectory_k_ub" :   1e20, # params.trajectory_k_ub,
"cexp" :              params.for_simulation['exp-unbinding-constant'],
"cb" :                params.for_simulation['cb'],
"cm" :                params.for_simulation['cm'],
"ct" :                params.for_simulation['ct'],
"eqb" :               params.for_simulation['eqb'],
"eqmpre" :            params.for_simulation['eqmpre'],
"eqmpost" :           params.for_simulation['eqmpost'],
"eqt" :               params.for_simulation['eqt'],
"ls" :                params.for_simulation['ls'],
"lt" :                params.for_simulation['lt'],
"radius_t" :          params.for_simulation['rt'],
"radius_m" :          params.for_simulation['rm'],
"radius_b" :          params.for_simulation['rb'],
"MT_binding_distance":            0.01, #nm MICROTUBULE_BINDING_DISTANCE from default_parameters.h
}

with open("data/paper_params.tex", "w") as f:
    for key in parameters:
        val = parameters[key]
        f.write(r'\newcommand\%s{%s}' % (latex_format(key).replace("_",""), latex_format(val)) + '\n')
