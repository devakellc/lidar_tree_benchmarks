#!/usr/bin/env python3
"""Smoke: ME + PyG + torch_points3d import on sm_120."""
import os
import sys

import torch
import MinkowskiEngine as ME
import torch_scatter
import torch_sparse
from torch_points_kernels import region_grow  # noqa: F401

print("torch", torch.__version__, "cuda", torch.version.cuda)
print("ME", ME.__version__)
print("arch", torch.cuda.get_arch_list())
print("device", torch.cuda.get_device_name(0))

sat_home = os.environ.get("SAT_HOME", "/opt/segment-any-tree")
sys.path.insert(0, sat_home)
import torch_points3d  # noqa: E402

print("torch_points3d", getattr(torch_points3d, "__version__", "ok"))
print("SAT_HOME", sat_home)
print("weights", os.path.exists(os.path.join(sat_home, "model_file/PointGroup-PAPER.pt")))

feats = torch.randn(10, 3, device="cuda")
coords = torch.randint(0, 10, (10, 4), dtype=torch.int32, device="cuda")
st = ME.SparseTensor(feats, coords)
conv = ME.MinkowskiConvolution(3, 8, kernel_size=3, dimension=3).cuda()
out = conv(st)
print("MinkowskiConvolution ok", out.F.shape)
