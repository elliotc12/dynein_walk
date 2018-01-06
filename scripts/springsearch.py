#!/usr/bin/env python
import os

os.system("mkdir -p plots/springsearch")

seeds = [1, 2, 3, 4]

sims = []
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.05", "num" : 5})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.8" , "num" : 6})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.12" , "num" : 7})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.15" , "num" : 8})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.6" , "num" : 9})

for sim in sims:
    for s in seeds:
        os.system("rq run --job-name springsearch-" + str(sim["num"]) \
                  + " python3 scripts/generate-paper-histogram-data.py" \
                  + " --cb " + sim["cb"] \
                  + " --cm " + sim["cm"] \
                  + " --ct " + sim["ct"] \
                  + " --seed " + str(s) \
                  + " --notpaper" \
                  + " --label springsearch-" + str(sim["num"]))
