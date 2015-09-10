#! /usr/bin/env python2

import random
import subprocess
import sys

speed = sys.argv[1]
length = sys.argv[2]
plot = sys.argv[3]

dt = 0.1

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
    runtime = 5000
elif length == "long":
    runtime = 1000
elif length == "normal":
    runtime = 500
elif length == "short":
    runtime = 100
else:
    runtime = int(sys.argv[2])


if len(sys.argv) >= 5 and sys.argv[4] == "loop":
        loop = True
else:
    loop = False

if (plot == "pentagon"):
    subprocess.call(["./walk", str(dt), str(runtime), ".6", ".2", ".8", ".4"])
    subprocess.call(["./plot.py", "speed=" + str(rate), "loop" if loop else ""])

if (plot == "foot-wiggle"):
    subprocess.call(["./walk", str(dt), str(runtime), ".6", ".2", ".8", ".7"])
    subprocess.call(["./plot.py", "speed=" + str(rate), "loop" if loop else ""])

if (plot == "mega-wiggle"):
    subprocess.call(["./walk", str(dt), str(runtime), "1.0", ".7", ".1", ".2"])
    subprocess.call(["./plot.py", "speed=" + str(rate), "loop" if loop else ""])

if (plot == "random"):
    random.seed()
    subprocess.call(["./walk", str(dt), str(runtime), str(2*random.random()), \
        str(2*random.random()), str(2*random.random()), str(2*random.random())])
    subprocess.call(["./plot.py", "speed=" + str(rate), "loop" if loop else ""])
