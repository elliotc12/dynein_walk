#!/usr/bin/python3

import numpy as np
import argparse


parser = argparse.ArgumentParser(description="Generate drunk vs sober stepping stats from stepping data")

parser.add_argument('-d', '--datafile', dest ='datafile', action='store', type=str,
                    help = "stepping data file")

args = parser.parse_args()

step_times = []
onebound_times = []
bothbound_times = []
step_lengths = []

data = np.loadtxt(args.datafile)
print("Data loaded...")
bind_times = np.array(data[:, 1])
unbind_times = np.array(data[:, 0])
near_positions = np.around(np.array(data[:, 2]), decimals=7)
far_positions = np.around(np.array(data[:, 3]), decimals=7)
near_step_lens = near_positions[1:] - near_positions[:-1]
far_step_lens = far_positions[1:] - far_positions[:-1]


assert len(near_positions) == len(far_positions)
Steps = []

for i in range(1, len(near_step_lens)):
    thisStep = np.array([near_positions[i], far_positions[i]])
    prevStep = np.array([near_positions[i-1], far_positions[i-1]])

    difference = prevStep-thisStep

    # notation: 0 for stride, 1 for stutter

    if difference[0] != 0.0 and difference[1] == 0.0:
        # near foot stepped#
        if difference[0] < 0:
            # negative displacement
            if prevStep[0] < prevStep[1]:  # check which foot stepped
                Steps.append([i, 1])
            else:
                Steps.append([i, 0])

        else:
            # positive displacement
            if prevStep[0] < prevStep[1]:
                Steps.append([i, 0])
            else:
                Steps.append([i, 1])  

    elif difference[0] == 0.0 and difference[1] != 0.0:
        # far foot stepped #
        if difference[1] < 0:
            # negative displacement
            if prevStep[1] < prevStep[0]:
                Steps.append([i, 1])
            else:
                Steps.append([i, 0])
        else:
            # positive displacement
            if prevStep[1] < prevStep[0]:
                Steps.append([i, 0])
            else:
                Steps.append([i, 1])
    
    elif difference[0] == 0.0 and difference[1] == 0.0:
        # nothing happened #
        print("NO STEP! current: {0} previous: {1}\n".format(thisStep, prevStep))

    else:
        print("Can not discern which stepped. Diff: {0}\n".format(difference))


#  now we can total everything up

stutter = 0
stride = 0
total = 0
for step in Steps:
    if step[1] == 1:
        stutter = stutter + 1
    elif step[1] == 0:
        stride = stride + 1
    else:
        pass
    total = total + 1

print("\nstutter: {0}, stride: {1}, total: {2}\n".format(stutter, stride, total))
    
    



