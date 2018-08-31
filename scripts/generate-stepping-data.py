#!/usr/bin/python3
import argparse
import os, sys
import numpy as np
import dynein.run as run

os.system("make generate_stepping_data")

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", type=int, help="Manually set the seed value, default=1", default=1)
parser.add_argument("-b", "--cb", type=float, help="Manually set the cb value, default=0.1", default=0.1)
parser.add_argument("-m", "--cm", type=float, help="Manually set the cm value, default=0.4", default=0.4)
parser.add_argument("-t", "--ct", type=float, help="Manually set the ct value, default=0.2", default=0.2)
parser.add_argument("-k", "--kb", type=float, help="Manually set the binding rate", default=3e9)
parser.add_argument("-u", "--kub", type=float, help="Manually set the unbinding rate", default=1000)
parser.add_argument("-y", "--ls", type=float, help="", default=10)
parser.add_argument("-w", "--lt", type=float, help="", default=10)

parser.add_argument("-F", "--force", type=float, help="", default=0)

parser.add_argument("-ee", "--rt", type=float, help="", default=8)
parser.add_argument("-ff", "--rm", type=float, help="", default=11)
parser.add_argument("-gg", "--rb", type=float, help="", default=3.5)

parser.add_argument("-aa", "--eqb", type=float, help="", default=10)
parser.add_argument("-bb", "--eqt", type=float, help="", default=10)
parser.add_argument("-cc", "--eqmpre", type=float, help="", default=10)
parser.add_argument("-dd", "--eqmpost", type=float, help="", default=10)

parser.add_argument("-x", "--unbindingconst", type=float, help="Manually set the unbinding const", default=0.0)
parser.add_argument("-r", "--runtime", type=float, help="Manually set the runtime value, default=0.1", default=0.1)
parser.add_argument("-f", "--framerate", type=float, help="Manually set the frame rate, default=1e-10", default=1)
parser.add_argument("-l", "--label", type=str, help="Manually set the label", default="default")
parser.add_argument("-v", "--movie", help="Movie flag", action="store_true")
parser.add_argument("-n", "--renametrajectory", help="Rename outputs to paper trajectory files", action="store_true")
parser.add_argument("-o", "--renamepaper", help="Rename outputs to paper exponential histogram files", action="store_true")
parser.add_argument("-R", "--renameforce", help="Rename outputs to force files", action="store_true")
args = parser.parse_args()

if args.movie:
    if (args.runtime / args.framerate > 1e5):
        print("Error: runtime/framerate > 1e5; this would result in a movie file larger than 1e5 lines. Please shorten.")
        exit(1)

basename = run.sim(**{"k_b": args.kb,
                      "k_ub": args.kub,
                      "cb": args.cb, # 0.1?
                      "cm": args.cm, # 0.4?
                      "ct": args.ct, # 0.2?
                      "ls": args.ls,
                      "lt": args.lt,
                      "eqb": args.eqb,  # from redwine 2012 supplemental
                      "eqmpre": args.eqmpre, # from burgess 2002, 360-160
                      "eqmpost": args.eqmpost, # from burgess 2002, 360-136
                      "eqt": args.eqt,
                      "rt": args.rt,
                      "rm": args.rm,
                      "rb": args.rb,
                      "force": args.force,
                      "exp-unbinding-constant": args.unbindingconst,
                      "dt": 1e-10, "label": args.label, "seed": args.seed, "runtime": args.runtime,
                      "framerate": args.framerate, "crash-movie": False, "nomovie": not args.movie})

if args.renamepaper:
    os.system("mv data/stepping_data_%s.txt data/paper_main_stepping_data-%s.txt" % (basename, args.seed))
    os.system("mv data/stepping_parameters_%s.tex data/paper_main_stepping_parameters.tex" % basename)

if args.renameforce:
    os.system("mv data/stepping_data_{}.txt data/paper_force_stepping_data-F-{},s-{}.txt".format(basename, args.force, args.seed))
    os.system("rm -f data/stepping_parameters_{}.tex".format(basename))
    os.system("rm -f data/stepping_movie_data_{}.tex".format(basename))
    os.system("rm -f data/stepping_parameters_{}.tex".format(basename))

if args.renametrajectory:
    os.rename("data/stepping_movie_data_%s.txt" % (basename), "data/paper_trajectory_movie_data.txt")
    os.rename("data/stepping_data_%s.txt" % (basename), "data/paper_trajectory_stepping_data.txt")
    os.rename("data/stepping_parameters_%s.tex" % (basename), "data/paper_trajectory_stepping_parameters.tex")
    os.unlink("data/stepping_movie_config_%s.txt" % basename)

os.system("rm data/stepping_config_%s.txt" % basename)
