#!/usr/bin/python3

import os
import dynein.run as run
import glob
import numpy as np
import argparse

def latex_format(x):
    if isinstance(x, float) or isinstance(x, int):
        x = '{:e}'.format(x)
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
    # if isinstance(x, str):
    #     x = x.replace('-', '_')
    return x

parser = argparse.ArgumentParser(description="script to generate unbinding probabilities of bb dynein")

parser.add_argument('-k_b', '--binding', dest='k_b', action='store', type=float,
                    default=1e8, help="pre-exponential binding constant", metavar='')
parser.add_argument('-k_ub', '--unbinding', dest='k_ub', action='store', type=float,
                    default=100, help="pre-exponential unbinding constant", metavar='')
parser.add_argument('-t', '--runtime', dest='runtime', action='store', type=float,
                    default=1.0, help='total runtime for simulation in seconds', metavar='')
parser.add_argument('-exp', '--exp-unbinding-constant', dest='exp_unbinding_constant',
                    action='store', type=float, default=0.0, help="exponential unbinding constant", metavar='')

parser.add_argument('-s', '--seed', dest ='seed', action='store', type=float, default=1.0, help ="random number seed", metavar='')
parser.add_argument('-cb', '--cb', dest ='cb', action='store', type=float, default=0.1, help ="cb", metavar='')
parser.add_argument('-cm', '--cm', dest ='cm', action='store', type=float, default=0.4, help ="cm", metavar='')
parser.add_argument('-ct', '--ct', dest ='ct', action='store', type=float, default=0.2, help ="ct", metavar='')

parser.add_argument('-f', '--logfile', dest='logfile', action='store', type=str, default='data/parameterSearch/testedParameters.csv', help=".txt file for logging stepping statistics", metavar='')

parser.add_argument('-l', '--label', dest='label', action='store', type=str, default='default', help="label for run", metavar='')

parser.add_argument('-w', '--writerate', dest='write_rate', action='store', type=str, default=1e6, help="writes per second", metavar='')

args = parser.parse_args()

if os.path.exists('run-unbinding-rate-simulations.py'):
    os.chdir('../')
os.system("make simulate_unbinding_rates")

if not os.path.exists('data/unbinding-probability/'):
        os.makedirs('data/unbinding-probability/')

for L in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]:
    basename = "%s__L-%s,s-%s" % (args.label, str(L), args.seed)

    cmd = ["rq run ./simulate_unbinding_rates --label %s --k_b %g --k_ub %g --c %g "\
           "--cb %g --cm %g --ct %g --ls 10.49 --lt 23.8 --eqb 120 --eqmpre 200 --eqmpost 224 --eqt 0 --write_rate %g "\
           "--runtime %g --seed %g --dt 1e-10 --L %g" %
           (str(args.label), float(args.k_b), float(args.k_ub), float(args.exp_unbinding_constant), float(args.cb),
            float(args.cm), float(args.ct), float(args.write_rate), float(args.runtime), args.seed, float(L))]

    if not os.path.exists('runlogs'):
        os.makedirs('runlogs')
    out = open('runlogs/' + basename + '.out', 'w')

    print("Running: ", " ".join(cmd), out)
    out.flush()
    # process_object = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE)
    # err = process_object.communicate()[1]
    # if (err != b''):
    #     print("\n##################################",
    #           "\nSimulation exited in error: \n\n",
    #           err.decode("utf-8"),
    #           "\n##################################\n\n")

with open("data/unbinding-probability/parameters_%s.tex" % args.label, "w") as f:
    f.write(r'\newcommand\%s{%s}' % (latex_format("label").replace("_",""), latex_format(args.label)) + '\n')
    f.write(r'\newcommand\%s{%s}' % (latex_format("kb").replace("_",""), latex_format(args.k_b)) + '\n')
    f.write(r'\newcommand\%s{%s}' % (latex_format("kub").replace("_",""), latex_format(args.k_ub)) + '\n')
    f.write(r'\newcommand\%s{%s}\n' %(latex_format("c").replace("_",""), latex_format(args.exp_unbinding_constant)))
    f.write(r'\newcommand\%s{%s}' % (latex_format("cb").replace("_",""), latex_format(args.cb))+ '\n')
    f.write(r'\newcommand\%s{%s}' % (latex_format("cm").replace("_",""), latex_format(args.cm))+ '\n')
    f.write(r'\newcommand\%s{%s}' % (latex_format("ct").replace("_",""), latex_format(args.ct))+ '\n')
    f.write(r'\newcommand\%s{%s}' % (latex_format("ls").replace("_",""), latex_format(10.49))+ '\n')
    f.write(r'\newcommand\%s{%s}' % (latex_format("lt").replace("_",""), latex_format(23.8))+ '\n')
    f.write(r'\newcommand\%s{%s}' % (latex_format("w_rate").replace("_",""), latex_format(args.write_rate))+ '\n')
