#!/usr/bin/env python
"""Stage one point cloud as FORAINetV2 instance-data for FF3D inference, WITHOUT
batch_load (which needs segmentator/open3d). Writes the *_vert/_sem/_ins/_bbox
npy that converter_forainetv2.ForAINetV2Data.get_infos consumes. Labels are dummy
(test/zero-shot): sem=0, ins=0, empty (0,7) bboxes -> gt_num=0. Centering matches
load_forainetv2_data.export (no-blue path)."""
import os
import numpy as np
import laspy

SCAN = os.environ.get("SCAN", "neonsjer054")
SRC = os.environ.get("SRC", "/payload/test.laz")
CLIP_R = float(os.environ.get("CLIP_R", "0"))   # 0 = no clip
ROOT = "/payload/ForestFormer3D/data/ForAINetV2"
inst = os.path.join(ROOT, "forainetv2_instance_data")
meta = os.path.join(ROOT, "meta_data")
os.makedirs(inst, exist_ok=True)
os.makedirs(meta, exist_ok=True)

las = laspy.read(SRC)
xyz = np.vstack([las.x, las.y, las.z]).T.astype(np.float64)
print(f"raw: {xyz.shape[0]} pts  extent x={xyz[:,0].ptp():.1f} y={xyz[:,1].ptp():.1f} z={xyz[:,2].ptp():.1f}")

if CLIP_R > 0:
    cx, cy = xyz[:, 0].mean(), xyz[:, 1].mean()
    keep = np.hypot(xyz[:, 0] - cx, xyz[:, 1] - cy) <= CLIP_R
    xyz = xyz[keep]
    print(f"clipped to {CLIP_R} m radius -> {xyz.shape[0]} pts")

# center: subtract mean x,y and min z (load_forainetv2_data.export, no-blue)
xyz[:, 0] -= xyz[:, 0].mean()
xyz[:, 1] -= xyz[:, 1].mean()
xyz[:, 2] -= xyz[:, 2].min()
pts = xyz.astype(np.float32)
n = pts.shape[0]

np.save(f"{inst}/{SCAN}_vert.npy", pts)
np.save(f"{inst}/{SCAN}_sem_label.npy", np.zeros(n, np.int64))
np.save(f"{inst}/{SCAN}_ins_label.npy", np.zeros(n, np.int64))
np.save(f"{inst}/{SCAN}_unaligned_bbox.npy", np.zeros((0, 7), np.float32))
np.save(f"{inst}/{SCAN}_aligned_bbox.npy", np.zeros((0, 7), np.float32))
np.save(f"{inst}/{SCAN}_axis_align_matrix.npy", np.eye(4))

with open(f"{meta}/test_list.txt", "w") as f:
    f.write(SCAN + "\n")
open(f"{meta}/train_list.txt", "w").close()
open(f"{meta}/val_list.txt", "w").close()
print(f"staged {SCAN}: {n} pts -> {inst}")
