#!/usr/bin/env python
import os

os.system("rm data/stepping_data_paperhisto_*")

seeds = [1, 2, 3, 4, 5]

for s in seeds:
    os.system("rq run generate-paper-histogram-data.py --seed " + str(s))
