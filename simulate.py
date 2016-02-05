#! /usr/bin/env python2

import random
import subprocess
import sys

speed = sys.argv[1]
length = sys.argv[2]
plot = sys.argv[3]

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
    runtime = 1000
else:
    runtime = int(sys.argv[2])


if len(sys.argv) >= 5 and sys.argv[4] == "loop":
        loop = True
else:
    loop = False

if (plot == "natural"):
    subprocess.call(["./walk", str(runtime), "0", "0", "0", "0"])
    subprocess.call(["./plot.py", "speed=" + str(rate), "loop" if loop else ""])

if (plot == "random"):
    random.seed()
    subprocess.call(["./walk", str(runtime), str(2*random.random()), \
        str(2*random.random()), str(2*random.random()), str(2*random.random())])
    subprocess.call(["./plot.py", "speed=" + str(rate), "loop" if loop else ""])
