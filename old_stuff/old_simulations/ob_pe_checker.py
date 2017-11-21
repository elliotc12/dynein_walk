#!/usr/bin/python2.7

import math
import numpy
import sys

kb = 1.3806e-5
eq_bba = 0.6 * math.pi
eq_bma = 0.6 * math.pi
eq_ta =  0.3 * math.pi
eq_uma = 0.6 * math.pi

fname_suffix = sys.argv[1]

config_txt = open("data/ob_config_" + fname_suffix + ".txt", 'r').read()

ct_st = config_txt.index("ct: ") + 4
ct_end = ct_st + config_txt[ct_st:].index('\n')
ct = float(config_txt[ct_st:ct_end])

cm_st = config_txt.index("cm: ") + 4
cm_end = cm_st + config_txt[cm_st:].index('\n')
cm = float(config_txt[cm_st:cm_end])

cb_st = config_txt.index("cb: ") + 4
cb_end = cb_st + config_txt[cb_st:].index('\n')
cb = float(config_txt[cb_st:cb_end])

T_st = config_txt.index("T: ") + 3
T_end = T_st + config_txt[T_st:].index('\n')
T = float(config_txt[T_st:T_end])

for (ang, c, eq) in [("bba", cb, eq_bba), ("bma", cm, eq_bma), ("ta", ct, eq_ta), ("uma", cm, eq_uma)]:
    pe_data = numpy.loadtxt("data/ob_" + ang + "_pe_" + fname_suffix + ".txt", skiprows=1)
    ang_data = numpy.loadtxt("data/ob_" + ang + "_angle_" + fname_suffix + ".txt", skiprows=1)
    for i in range(len(pe_data[:,1])):
        a = ang_data[i,1] - eq
        pe = pe_data[i,1]
        print 0.5*c*a*a - pe*(0.5*kb*T)
