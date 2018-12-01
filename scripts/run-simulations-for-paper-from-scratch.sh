#!/bin/sh

set -ev

rq run -J what-does-this-do python3 scripts/generate-stepping-data.py --ls 20.75 --lt 23 --kub 30 --kb 3000000.0 --cb 0.08 --cm 1.05 --ct 0.36 --eqb 120 --eqt 0 --eqmpre 197 --eqmpost 242 --unbindingconst -0.35 --runtime 2 --framerate 1e-10 --seed 1 --label paper_stroke_angles_1 --renameangles --longanglemode

# X elliotc12 bingley    2h:10m  0:32    1 python3 scripts/generate-stepping-data.py --ls 20.75 --lt 23 --kub 30 --kb 3000000.0 --cb 0.08 --cm 1.05 --ct 0.36 --eqb 120 --eqt 0 --eqmpre 197 --eqmpost 242 --unbindingconst -0.35 --runtime 2 --framerate 1e-10 --seed 5 --label paper_stroke_angles_5 --renameangles --longanglemode
# X elliotc12 bingley    2h:10m  0:32    1 python3 scripts/generate-stepping-data.py --ls 20.75 --lt 23 --kub 30 --kb 3000000.0 --cb 0.08 --cm 1.05 --ct 0.36 --eqb 120 --eqt 0 --eqmpre 197 --eqmpost 242 --unbindingconst -0.35 --runtime 2 --framerate 1e-10 --seed 4 --label paper_stroke_angles_4 --renameangles --longanglemode
# X elliotc12 bingley    2h:10m  0:32    1 python3 scripts/generate-stepping-data.py --ls 20.75 --lt 23 --kub 30 --kb 3000000.0 --cb 0.08 --cm 1.05 --ct 0.36 --eqb 120 --eqt 0 --eqmpre 197 --eqmpost 242 --unbindingconst -0.35 --runtime 2 --framerate 1e-10 --seed 3 --label paper_stroke_angles_3 --renameangles --longanglemode
# X elliotc12 bingley    2h:10m  0:02    1 python3 scripts/generate-stepping-data.py --ls 20.75 --lt 23 --kub 30 --kb 3000000.0 --cb 0.08 --cm 1.05 --ct 0.36 --eqb 120 --eqt 0 --eqmpre 197 --eqmpost 242 --unbindingconst -0.35 --runtime 2 --framerate 1e-10 --seed 1 --label paper_stroke_angles_1 --renameangles --longanglemode
# F elliotc12 darcy      26:20   0:08    1 python3 scripts/generate-stepping-data.py --ls 20.75 --lt 23 --kub 30 --kb 3000000.0 --cb 0.08 --cm 1.05 --ct 0.36 --eqb 120 --eqt 0 --eqmpre 197 --eqmpost 242 --unbindingconst -0.35 --runtime 2 --framerate 1e-10 --seed 2 --label paper_stroke_angles_2 --renameangles --longanglemode
