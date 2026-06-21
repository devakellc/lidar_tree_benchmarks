#!/usr/bin/env python
# Seeded SAM2Point arm (#P3): a LAZ clip + CHM-VWF apex prompts -> per-point crown
# instance labels in a `sam2point` int32 extra dim. Adapts SAM2Point's main.py to
# (a) read a LiDAR LAZ (ground removed), (b) normalize to a unit cube, (c) feed
# EACH world-coordinate apex as a 3-D point prompt, and (d) assemble the per-prompt
# voxel masks into per-point instances (first-claim wins). LiDAR has no RGB, so the
# voxel feature is the normalized height replicated to 3 channels.
#   python run_sam2point_arm.py --input clip.laz --prompts apex.csv --output out.laz
import sys, os, argparse, numpy as np, laspy, torch
sys.path.insert(0, "/workspace/SAM2Point")
os.chdir("/workspace/SAM2Point")                      # seg_point writes ./video frames
from sam2point.voxelizer import Voxelizer
import segment as S

ap = argparse.ArgumentParser()
ap.add_argument("--input", required=True); ap.add_argument("--prompts", required=True)
ap.add_argument("--output", required=True)
ap.add_argument("--voxel_size", type=float, default=0.02)   # in the [0,1] normalized cube
ap.add_argument("--mode", default="nearest"); ap.add_argument("--theta", type=float, default=0.5)
ap.add_argument("--dataset", default="NEON"); ap.add_argument("--prompt_type", default="point")
ap.add_argument("--sample_idx", type=int, default=0)
args = ap.parse_args()
torch.autocast(device_type="cuda", dtype=torch.bfloat16).__enter__()
if torch.cuda.get_device_properties(0).major >= 8:
    torch.backends.cuda.matmul.allow_tf32 = True; torch.backends.cudnn.allow_tf32 = True

las = laspy.read(args.input)
cls = np.asarray(las.classification); veg = cls != 2
if veg.sum() < 50: veg = np.ones(len(las.x), bool)
P = np.column_stack([np.asarray(las.x)[veg], np.asarray(las.y)[veg],
                     np.asarray(las.z)[veg]]).astype(np.float64)
mn = P.min(0); scale = (P.max(0) - mn).max()          # uniform scale -> preserve aspect
if scale <= 0: scale = 1.0
Pn = (P - mn) / scale
pr = np.loadtxt(args.prompts, delimiter=",", ndmin=2)[:, :3]
prn = (pr - mn) / scale
col = np.repeat(Pn[:, 2:3], 3, axis=1).astype(np.float64)   # normalized height as pseudo-RGB

vox = Voxelizer(voxel_size=args.voxel_size, clip_bound=None)
locs, feats, _, inds = vox.voxelize(Pn, col, Pn[:, :1].astype(int))
point_locs = locs[inds].astype(int)
labels = np.zeros(Pn.shape[0], np.int32)
for i, p in enumerate(prn):
    a = argparse.Namespace(**vars(args)); a.prompt_idx = 0; a.sample_idx = i
    try:
        mask = S.seg_point(locs, feats, np.array([p], dtype=np.float64), a)
        m = np.asarray(mask)
        pm = m[point_locs[:, 0], point_locs[:, 1], point_locs[:, 2]].astype(bool)
    except Exception as e:
        sys.stderr.write("prompt %d failed: %s\n" % (i, e)); continue
    labels[(labels == 0) & pm] = i + 1                # first-claim, 1-based
print("sam2point: %d veg points, %d/%d prompts produced masks, %d instances" %
      (len(labels), len(np.unique(labels[labels > 0])), len(prn), len(np.unique(labels[labels > 0]))))
out = laspy.LasData(las.header); out.points = las.points[veg].copy()
out.add_extra_dim(laspy.ExtraBytesParams(name="sam2point", type="int32", description="sam2point instance"))
out.sam2point = labels
out.write(args.output)
