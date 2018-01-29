#!/usr/bin/env python
import os

os.system("mkdir -p plots/springsearch")

seeds = [1, 2, 3, 4]

sims = []
sims.append({"cb" : "0.1", "cm" : "1.5", "ct" : "0.6", "kb" : "1e14", "kub" : "420",  "num" : 32})
sims.append({"cb" : "0.1", "cm" : "1.5", "ct" : "0.6", "kb" : "5e13", "kub" : "420",  "num" : 32})

for sim in sims:
    for s in seeds:
        os.system("rq run --job-name springsearch-" + str(sim["num"]) \
                  + " python3 scripts/generate-paper-histogram-data.py" \
                  + " --kub " + sim["kub"] \
                  + " --kb " + sim["kb"] \
                  + " --cb " + sim["cb"] \
                  + " --cm " + sim["cm"] \
                  + " --ct " + sim["ct"] \
                  + " --seed " + str(s) \
                  + " --notpaper" \
                  + " --label springsearch-" + str(sim["num"]))
