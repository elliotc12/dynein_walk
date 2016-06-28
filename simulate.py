#! /usr/bin/env python2

import random
import subprocess
import sys

speed = sys.argv[1]
length = sys.argv[2]

if len(sys.argv) > 3:
    plot = sys.argv[3]
else:
    plot = "natural"

if speed == "veryfast":
    rate = 150
elif speed == "fast":
    rate = 20
elif speed == "normal":
    rate = 10
elif speed == "slow":
    rate = 1
else:
    rate = int(sys.argv[1])

if length == "verylong":
    runtime = 100000
elif length == "long":
    runtime = 10000
elif length == "normal":
    runtime = 5000
elif length == "short":
    runtime = 50
else:
    runtime = int(sys.argv[2])

if len(sys.argv) >= 5:
    if sys.argv[4] == "loop":
        flag = "loop"
    elif sys.argv[4] == "step":
        flag = "step"
    elif sys.argv[4] == "save":
        if len(sys.argv) < 6:
            print "Usage: make save NAME=savename"
            exit(-1)
        flag = "save" + "/" + sys.argv[5]
else:
    flag = ""

def verbose_run(command):
    cmd = ' '.join(command)
    subprocess.check_call(command)

if (plot == "natural"):
    verbose_run(["./plot", str(runtime), "0", "0", "0", "0"])
    subprocess.call(["./plot.py", "speed=" + str(rate), flag])

if (plot == "pretty"):
    verbose_run(["./plot", str(runtime), "0", "0", "0.5", "0"])
    subprocess.call(["./plot.py", "speed=" + str(rate), flag])

if (plot == "random"):
    random.seed()
    subprocess.check_call(["./plot", str(runtime), str(2*random.random()), \
        str(2*random.random()), str(2*random.random()), str(2*random.random())])
    subprocess.call(["./plot.py", "speed=" + str(rate), flag])
