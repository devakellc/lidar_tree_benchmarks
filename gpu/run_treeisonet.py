#!/usr/bin/env python3
# #M7 headless driver (apex-only). treeLoc -> postPeakExtraction -> tops CSV in
# UTM with ABSOLUTE Z (R converts to height-above-ground via ground_dtm.tif).
# Usage: run_treeisonet.py <input.laz> <out.csv> <loc.pth> <loc.json> [voxel] [conf]
#   voxel <= 0 -> checkpoint native [0.1,0.1,0.2]; >0 -> isotropic override.
import os, sys, numpy as np, laspy
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "TreeAIBox"))
from modules.treeisonet.treeLoc import treeLoc, postPeakExtraction

inp, out, loc_pth, loc_cfg = sys.argv[1:5]
voxel = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
conf  = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5
las = laspy.read(inp)
try:
    epsg = las.header.parse_crs().to_epsg()
except Exception:
    epsg = None
if epsg == 3857:
    sys.exit("ERROR: input is EPSG:3857 (web-mercator); reproject to metric UTM")
pcd = np.transpose([np.asarray(las.x), np.asarray(las.y), np.asarray(las.z)]).astype(float)
pmin = pcd.min(0); pcd[:, :3] -= pmin
cr = np.array([voxel, voxel, voxel]) if voxel > 0 else np.zeros(3)
preds = np.atleast_2d(treeLoc(loc_cfg, pcd, loc_pth, use_cuda=True,
                              if_stem=False, custom_resolution=cr))
sel = preds[preds[:, 3] > conf] if preds.size else preds   # col 3 = confidence
n = sel.shape[0]
if n == 0:                                                 # no peaks -> 0 tops
    np.savetxt(out, np.empty((0, 3)), header="x y z", comments=""); sys.exit(0)
tops = sel[:, :3] if n == 1 else \
       np.atleast_2d(postPeakExtraction(sel, K=min(5, n)))[:, :3]
# TreeLoc gives a tree LOCATION (z ~ stem base); snap z to the local canopy max
# within `snap_r` so the apex height matches the scorer's height gate.
snap_r = 2.0
kd = cKDTree(pcd[:, :2])
for i in range(tops.shape[0]):
    idx = kd.query_ball_point(tops[i, :2], snap_r)
    if idx:
        tops[i, 2] = pcd[idx, 2].max()
tops = tops + pmin[:3]                                      # back to UTM + height frame
np.savetxt(out, tops, header="x y z", comments="")
print(f"wrote {len(tops)} tops -> {out}")
