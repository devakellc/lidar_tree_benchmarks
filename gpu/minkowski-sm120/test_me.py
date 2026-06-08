#!/usr/bin/env python3
import torch
import MinkowskiEngine as ME

print("torch", torch.__version__, "cuda", torch.version.cuda)
print("ME", ME.__version__)
print("arch", torch.cuda.get_arch_list())
print("device", torch.cuda.get_device_name(0))

feats = torch.randn(10, 3, device="cuda")
coords = torch.randint(0, 10, (10, 4), dtype=torch.int32, device="cuda")
st = ME.SparseTensor(feats, coords)
print("SparseTensor ok", st.F.shape, st.F.device)

# tiny conv to exercise CUDA kernels
conv = ME.MinkowskiConvolution(3, 8, kernel_size=3, dimension=3).cuda()
out = conv(st)
print("MinkowskiConvolution ok", out.F.shape, out.F.device)
