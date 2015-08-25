#! /usr/bin/env python2

import subprocess
import sys

plot = sys.argv[1]

if (plot == "pentagon"):
    subprocess.call(["./walk", ".6", ".2", ".8", ".4"])
    subprocess.call(["./plot.py"])

if (plot == "foot-wiggle"):
    subprocess.call(["./walk", ".6", ".2", ".8", ".7"])
    subprocess.call(["./plot.py"])

if (plot == "mega-wiggle"):
    subprocess.call(["./walk", "1.0", ".7", ".1", ".2"])
    subprocess.call(["./plot.py"])
