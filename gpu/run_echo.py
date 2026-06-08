#!/usr/bin/env python3
# #I3 handshake proof: read the input cloud, emit 3 fixed detections near its
# centroid. No model, no GPU -- proves R -> python -> CSV -> scorer end to end.
import sys, numpy as np, laspy
inp, out = sys.argv[1], sys.argv[2]
las = laspy.read(inp)
cx, cy = float(np.median(las.x)), float(np.median(las.y))
det = np.array([[cx, cy, 20.0], [cx + 5, cy, 18.0], [cx, cy + 5, 15.0]])
np.savetxt(out, det, header="x y z", comments="")
