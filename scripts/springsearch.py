#!/usr/bin/env python
import os

os.system("mkdir -p plots/springsearch")

seeds = [1, 2, 3, 4]

sims = []
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.01"})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.1"})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.2"})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.3"})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.4"})

simnum = 0

for sim in sims:
    for s in seeds:
        os.system("rq run --job-name springsearch-" + str(simnum) \
                  + " python3 scripts/generate-paper-histogram-data.py" \
                  + " --cb " + sim["cb"] \
                  + " --cm " + sim["cm"] \
                  + " --ct " + sim["ct"] \
                  + " --seed " + str(s) \
                  + " --notpaper" \
                  + " --label springsearch-" + str(simnum))
    simnum += 1
