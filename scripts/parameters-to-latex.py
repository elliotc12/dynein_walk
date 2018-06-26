#! /usr/bin/env python3

import numpy as np
import sys, os

sys.path.insert(0, os.getcwd() + "/data/")
import paper_params as params

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

parameters = {
"k_b" :               params.k_b,
"k_ub" :              params.k_ub,
"trajectory_k_b" :    params.trajectory_k_b,
"trajectory_k_ub" :   params.trajectory_k_ub,
"c" :                 params.c,
"cb" :                params.cb,
"cm" :                params.cm,
"ct" :                params.ct,
"eqb" :               params.eqb,
"eqmpre" :            params.eqmpre,
"eqmpost" :           params.eqmpost,
"eqt" :               params.eqt,
"ls" :                params.ls,
"lt" :                params.lt,
"radius_t" :          params.radius_t,
"radius_m" :          params.radius_m,
"radius_b" :          params.radius_b
}

with open("data/paper_params.tex", "w") as f:
    for key in parameters:
        val = parameters[key]
        f.write(r'\newcommand\%s{%s}' % (latex_format(key).replace("_",""), latex_format(val)) + '\n')
