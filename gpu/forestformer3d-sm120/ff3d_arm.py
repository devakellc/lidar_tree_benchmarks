#!/usr/bin/env python
"""FF3D benchmark arm: stage every cylinder clip in <in_dir> as a ForAINetV2
scene, run FF3D once over all of them, and write ONE merged UTM LAZ <out_laz>
(point_source_id = per-cylinder instance id, user_data = cylinder/block index,
instance score in an extra dim). Run from the FF3D repo root.

Usage: python ff3d_arm.py <in_dir> <out_laz> <ckpt>
"""
import os, sys, glob, shutil
import numpy as np
import laspy

_orig = __import__("torch").load
def _load(*a, **k):
    k.setdefault("weights_only", False); return _orig(*a, **k)
__import__("torch").load = _load

IN_DIR, OUT_LAZ, CKPT = sys.argv[1], sys.argv[2], sys.argv[3]
# This driver is mounted from outside the repo, so sys.path[0] is its own dir,
# not the repo root. Add the repo root (cwd, set by ff3d_entry.sh `cd $repo`) so
# the repo-local `oneformer3d` package + tools/ helpers import.
sys.path.insert(0, os.getcwd()); sys.path.insert(0, os.path.join(os.getcwd(), "tools"))
ROOT = "data/ForAINetV2"                       # relative to repo CWD (ephemeral)
INST = os.path.join(ROOT, "forainetv2_instance_data")
META = os.path.join(ROOT, "meta_data")
for d in (INST, META):
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d, exist_ok=True)

# 1. stage each cylinder as a scene; remember its centering offset for UTM restore
clips = sorted(glob.glob(os.path.join(IN_DIR, "cyl_*.laz")))
offsets, scans = {}, []
for clip in clips:
    scan = os.path.splitext(os.path.basename(clip))[0]        # cyl_000
    las = laspy.read(clip)
    xyz = np.vstack([las.x, las.y, las.z]).T.astype(np.float64)
    if xyz.shape[0] < 50:                                     # too sparse -> skip
        continue
    off = np.array([xyz[:, 0].mean(), xyz[:, 1].mean(), xyz[:, 2].min()])
    offsets[scan] = off
    pts = (xyz - off).astype(np.float32)
    n = pts.shape[0]
    np.save(f"{INST}/{scan}_vert.npy", pts)
    np.save(f"{INST}/{scan}_sem_label.npy", np.zeros(n, np.int64))
    np.save(f"{INST}/{scan}_ins_label.npy", np.zeros(n, np.int64))
    np.save(f"{INST}/{scan}_unaligned_bbox.npy", np.zeros((0, 7), np.float32))
    np.save(f"{INST}/{scan}_aligned_bbox.npy", np.zeros((0, 7), np.float32))
    np.save(f"{INST}/{scan}_axis_align_matrix.npy", np.eye(4))
    scans.append(scan)
with open(f"{META}/test_list.txt", "w") as f:
    f.write("\n".join(scans) + "\n")
open(f"{META}/train_list.txt", "w").close(); open(f"{META}/val_list.txt", "w").close()
if not scans:
    print("no usable cylinders"); laspy.LasData(laspy.LasHeader(point_format=3,
        version="1.2")).write(OUT_LAZ); print("ARM_DONE"); sys.exit(0)

# 2. build points/*.bin + the test pkl IN-PROCESS, updating ONLY the test pkl.
# The stock tools/create_data_forainetv2.py also runs update_pkl_infos on the
# train/val pkls, which UnboundLocalErrors on our empty train/val splits (#27
# note). Replicate just create_info_file + the test-pkl update.
from mmdet3d.utils import register_all_modules
register_all_modules()
from converter_forainetv2 import create_info_file
from update_infos_to_v2 import update_pkl_infos
create_info_file(ROOT, "forainetv2", ROOT, workers=4)
update_pkl_infos("forainetv2", out_dir=ROOT,
                 pkl_path=os.path.join(ROOT, "forainetv2_oneformer3d_infos_test.pkl"))

# 3. one runner; iterate all scenes; collect UTM-restored labelled points
from mmengine.config import Config, ConfigDict
from mmengine.runner import Runner
import oneformer3d  # noqa: F401
cfg = Config.fromfile("configs/oneformer3d_qs_radius16_qp300_2many.py")
cfg.work_dir = "./work_dirs/arm"; cfg.load_from = CKPT
if cfg.model.get("test_cfg") is None:
    cfg.model.test_cfg = ConfigDict()
cfg.model.test_cfg["output_dir"] = cfg.work_dir
runner = Runner.from_cfg(cfg); runner.load_or_resume(); model = runner.model.eval()
import torch
def npy(x): return x.detach().cpu().numpy() if torch.is_tensor(x) else np.asarray(x)

XS, YS, ZS, INST_ID, BLK, SCORE = [], [], [], [], [], []
for block_i, data in enumerate(runner.test_dataloader):
    scan = scans[block_i]                                     # dataloader is ordered
    with torch.no_grad():
        res = model.test_step(data, epoch=0)
    seg = res[0].pred_pts_seg
    per_point = npy(seg.pts_instance_mask[1]).astype(np.int64)  # (N,)
    pts = npy(data["inputs"]["points"][0] if isinstance(
        data["inputs"]["points"], list) else data["inputs"]["points"])[:, :3]
    n = min(len(pts), len(per_point))
    off = offsets[scan]
    XS.append(pts[:n, 0] + off[0]); YS.append(pts[:n, 1] + off[1])
    ZS.append(pts[:n, 2] + off[2]); INST_ID.append(per_point[:n])
    BLK.append(np.full(n, block_i, np.int64))
    # per-instance score broadcast to points (diagnostic only)
    scores = npy(seg.instance_scores).astype(np.float32)
    sc = np.zeros(n, np.float32)
    for k, s in enumerate(scores):
        sc[per_point[:n] == (k + 1)] = s                      # ids are 1-based
    SCORE.append(sc)
    print(f"{scan}: {n} pts, {len(np.unique(per_point[per_point>0]))} instances")

X = np.concatenate(XS); Y = np.concatenate(YS); Z = np.concatenate(ZS)
inst = np.concatenate(INST_ID); blk = np.concatenate(BLK); score = np.concatenate(SCORE)

# 4. write the merged UTM LAZ
h = laspy.LasHeader(point_format=3, version="1.2")
h.offsets = np.array([X.min(), Y.min(), Z.min()]); h.scales = [0.001, 0.001, 0.001]
las = laspy.LasData(h)
las.x, las.y, las.z = X, Y, Z
las.point_source_id = np.clip(inst, 0, 65535).astype(np.uint16)
las.user_data = np.clip(blk, 0, 255).astype(np.uint8)
las.add_extra_dim(laspy.ExtraBytesParams(name="ff3d_score", type=np.float32))
las.ff3d_score = score
las.write(OUT_LAZ)
print(f"wrote {OUT_LAZ}: {len(X)} pts from {len(scans)} cylinders"); print("ARM_DONE")
