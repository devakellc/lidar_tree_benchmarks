#!/usr/bin/env python3
# #M7 crown variant (#20): treeLoc -> postPeakExtraction seeds -> treeOff offset
# net -> per-point instance ids. Writes the labeled CANOPY points (height >= hmin
# in the normalized frame, so ground returns don't inflate the crown footprint)
# as a CSV (x y z crown_id) in UTM. R reduces these to apexes (reduce_instances)
# and crown diameters (crown_diameter_table) for the crown-diameter benchmark.
#
# Usage: run_treeisonet_crowns.py <input.laz> <out.csv> <loc.pth> <loc.json>
#            <off.pth> <off.json> [voxel] [conf] [hmin]
import os, sys, numpy as np, laspy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "TreeAIBox"))
from modules.treeisonet.treeLoc import treeLoc, postPeakExtraction
from modules.treeisonet.treeOff import treeOff

inp, out, loc_pth, loc_cfg, off_pth, off_cfg = sys.argv[1:7]
voxel = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
conf  = float(sys.argv[8]) if len(sys.argv) > 8 else 0.22
hmin  = float(sys.argv[9]) if len(sys.argv) > 9 else 2.0
HEADER = "x y z crown_id"
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
sel = preds[preds[:, 3] > conf] if preds.size else preds
if sel.shape[0] < 2:                                       # no seeds -> no crowns
    open(out, "w").write(HEADER + "\n"); sys.exit(0)
tops = np.atleast_2d(postPeakExtraction(sel, K=min(5, sel.shape[0])))[:, :3]
ids = np.asarray(treeOff(off_cfg, pcd, tops, off_pth, use_cuda=True))  # per-point 1..N (0 = unassigned)
xyz = pcd[:, :3] + pmin[:3]                                # back to UTM + height frame
keep = (ids > 0) & (pcd[:, 2] >= hmin)                     # assigned canopy points only
if not keep.any():
    open(out, "w").write(HEADER + "\n"); sys.exit(0)
arr = np.column_stack([xyz[keep], ids[keep].astype(int)])
np.savetxt(out, arr, header=HEADER, comments="",
           fmt=["%.3f", "%.3f", "%.3f", "%d"])
print(f"wrote {int(keep.sum())} canopy pts in {len(np.unique(ids[keep]))} crowns -> {out}")
