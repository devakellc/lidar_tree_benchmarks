#!/usr/bin/env python
"""Run FF3D on the staged plot and extract predicted tree instances.
pred_pts_seg.pts_instance_mask = [ (K,N) bool masks , (N,) per-point instance ids ].
Reports the model's own per-point assignment + a score-threshold sweep, and saves
a viewable instance-labelled cloud (.laz, instance id in point_source_id)."""
import os
import numpy as np
import torch

_orig = torch.load
def _load(*a, **k):
    k.setdefault("weights_only", False)
    return _orig(*a, **k)
torch.load = _load

from mmengine.config import Config, ConfigDict
from mmengine.runner import Runner
import oneformer3d  # noqa: F401

cfg = Config.fromfile("configs/oneformer3d_qs_radius16_qp300_2many.py")
cfg.work_dir = "./work_dirs/infer"
cfg.load_from = "work_dirs/clean_forestformer/epoch_3000_fix.pth"
if cfg.model.get("test_cfg") is None:
    cfg.model.test_cfg = ConfigDict()
cfg.model.test_cfg["output_dir"] = cfg.work_dir

runner = Runner.from_cfg(cfg)
runner.load_or_resume()
model = runner.model.eval()

def np_(x):
    return x.detach().cpu().numpy() if torch.is_tensor(x) else np.asarray(x)

for data in runner.test_dataloader:
    with torch.no_grad():
        res = model.test_step(data, epoch=0)
    seg = res[0].pred_pts_seg

    masks = np_(seg.pts_instance_mask[0]).astype(bool)      # (K, N)
    per_point = np_(seg.pts_instance_mask[1]).astype(np.int64)  # (N,)
    scores = np_(seg.instance_scores).astype(float)         # (K,)
    labels = np_(seg.instance_labels).astype(np.int64)      # (K,)
    K, N = masks.shape

    print("\n================ FF3D PREDICTION (NEON SJER plot, zero-shot) ================")
    print(f"points: {N:,}   raw instance candidates: {K}")
    print(f"instance score range: {scores.min():.3f} .. {scores.max():.3f}  (cfg score_th=0.4)")
    # model's own per-point assignment
    final_ids = np.unique(per_point[per_point > 0])
    print(f"model per-point assignment: {len(final_ids)} distinct tree ids, "
          f"{(per_point > 0).mean()*100:.1f}% of points assigned to a tree")
    # threshold sweep over candidate scores
    for thr in (0.05, 0.10, 0.15, 0.20):
        sel = scores >= thr
        sizes = masks[sel].sum(1)
        sizes = sizes[sizes >= 50]   # ignore specks
        print(f"  score>={thr:.2f}: {sel.sum():3d} candidates, "
              f"{len(sizes):3d} with >=50 pts (median {int(np.median(sizes)) if len(sizes) else 0} pts)")

    # save the model's per-point instances as a viewable LAZ + npy
    pts = np_(data["inputs"]["points"][0] if isinstance(data["inputs"]["points"], list)
              else data["inputs"]["points"])[:, :3]
    n = min(len(pts), len(per_point))
    out = np.column_stack([pts[:n], per_point[:n]]).astype(np.float32)
    np.save(os.path.join(cfg.work_dir, "pred_instances.npy"), out)
    try:
        import laspy
        h = laspy.LasHeader(point_format=3, version="1.2")
        las = laspy.LasData(h)
        las.x, las.y, las.z = pts[:n, 0], pts[:n, 1], pts[:n, 2]
        las.point_source_id = np.clip(per_point[:n], 0, 65535).astype(np.uint16)
        las.write(os.path.join(cfg.work_dir, "pred_instances.laz"))
        print(f"saved: {cfg.work_dir}/pred_instances.laz  (instance id in point_source_id)")
    except Exception as e:
        print("laz save skipped:", e)
    print(f"saved: {cfg.work_dir}/pred_instances.npy  (x,y,z,instance_id)  {out.shape}")
print("INFER_DONE")
