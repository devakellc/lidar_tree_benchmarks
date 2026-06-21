#!/usr/bin/env python
# Non-interactive Treeiso runner: input LAZ -> per-point tree-instance labels in
# a "treeiso" int32 extra dim. Runs init/intermediate/final segs on ONE
# decimation (upstream main() mismatches dec1/dec2 indices). Ground (class 2) is
# removed first (Treeiso expects vegetation only).
import sys, os, numpy as np, laspy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import treeiso as T

inp, outp = sys.argv[1], sys.argv[2]
las = laspy.read(inp)
x = np.asarray(las.x, float); y = np.asarray(las.y, float); z = np.asarray(las.z, float)
cls = np.asarray(las.classification)
veg = cls != 2
if veg.sum() < 50: veg = np.ones(len(x), bool)
xv, yv, zv = x[veg], y[veg], z[veg]
work = np.column_stack([xv, yv, zv]); work = work - work.mean(axis=0)
uidx, inv = T.decimate_pcd(work[:, :3], T.PR_DECIMATE_RES1)
pdec = work[uidx]
init = np.asarray(T.init_segs(pdec[:, :3]))
inter = np.asarray(T.intermediate_segs(np.column_stack([pdec[:, :3], init])))
final = np.asarray(T.final_segs(np.column_stack([pdec[:, :3], init, inter])))
labels = (final[inv] + 1).astype(np.int32)  # 1-based: 0 is the reader's "unassigned"
out = laspy.LasData(las.header); out.points = las.points[veg].copy()
out.add_extra_dim(laspy.ExtraBytesParams(name="treeiso", type="int32", description="treeiso instance"))
out.treeiso = labels
out.write(outp)
print("treeiso: %d veg points, %d instances -> %s" % (len(labels), len(np.unique(labels)), outp))
