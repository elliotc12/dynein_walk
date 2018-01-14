#!/usr/bin/env python
import os

os.system("mkdir -p plots/springsearch")

seeds = [1, 2, 3, 4]

sims = []
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.6", "kub" : "1000",  "num" : 19})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.6", "kub" : "3000",  "num" : 20})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.6", "kub" : "5000",  "num" : 21})
sims.append({"cb" : "0.1", "cm" : "0.4", "ct" : "0.6", "kub" : "10000", "num" : 22})

for sim in sims:
    for s in seeds:
        os.system("rq run --job-name springsearch-" + str(sim["num"]) \
                  + " python3 scripts/generate-paper-histogram-data.py" \
                  + " --kub " + sim["kub"] \
                  + " --cb " + sim["cb"] \
                  + " --cm " + sim["cm"] \
                  + " --ct " + sim["ct"] \
                  + " --seed " + str(s) \
                  + " --notpaper" \
                  + " --label springsearch-" + str(sim["num"]))
