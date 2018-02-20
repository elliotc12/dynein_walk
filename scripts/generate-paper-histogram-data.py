#!/usr/bin/python3
import argparse
import os, sys
import numpy as np
import dynein.run as run

runtime = 1e-1

os.system("make generate_stepping_data")

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", type=int, help="Manually set the seed value, default=1", default=1)
parser.add_argument("-b", "--cb", type=float, help="Manually set the cb value, default=0.1", default=0.1)
parser.add_argument("-m", "--cm", type=float, help="Manually set the cm value, default=0.4", default=0.4)
parser.add_argument("-t", "--ct", type=float, help="Manually set the ct value, default=0.2", default=0.2)
parser.add_argument("-l", "--label", type=str, help="Manually set the label", default="paperhisto")
parser.add_argument("-k", "--kb", type=float, help="Manually set the binding rate", default=3e9)
parser.add_argument("-u", "--kub", type=float, help="Manually set the unbinding rate", default=1000)
parser.add_argument("-p", "--notpaper", help="Not making paper plot", action="store_true")
args = parser.parse_args()

basename = run.sim(**{"k_b": args.kb, # 6e9?
                      "k_ub": args.kub, # 1000?
                      "cb": args.cb, # 0.1?
                      "cm": args.cm, # 0.4?
                      "ct": args.ct, # 0.2?
                      "ls": 10.49, # from urnavicius 2015 (paper.bib)
                      "lt": 23.8,  # from urnavicius 2015
                      "eqb": 120,  # from redwine 2012 supplemental
                      "eqmpre": 200, # from burgess 2002, 360-160
                      "eqmpost": 224, # from burgess 2002, 360-136
                      "eqt": 0,
                      "dt": 1e-10, "label": args.label, "seed": args.seed, "runtime": runtime,
                      "framerate": 1e-10, "crash-movie": False, "nomovie": True})

if not args.notpaper:
    os.system("mv data/stepping_data_%s.txt data/paper_histogram_stepping_data-%s.txt" % (basename, args.seed))
    os.system("mv data/stepping_parameters_%s.tex data/paper_histogram_stepping_parameters.tex" % basename)

os.system("rm data/stepping_movie_config_%s.txt" % basename)
os.system("rm data/stepping_movie_data_%s.txt" % basename)
os.system("rm data/stepping_config_%s.txt" % basename)
