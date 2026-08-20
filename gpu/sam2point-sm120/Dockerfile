# SAM2Point (#P3) seed->refine arm on Blackwell (sm_120). Builds FROM the proven
# MinkowskiEngine sm_120 base (torch 2.7.0 / cu128, TORCH_CUDA_ARCH_LIST=12.0)
# the SAT/FF3D arms already use, so no PyTorch-on-Blackwell rebuild is needed.
# SAM2Point (Apache-2.0) voxelizes a point cloud, slices it into a "video", and
# runs SAM2's video predictor from a 3D point prompt. The SAM2 connected-
# components CUDA ext is left UNCOMPILED (optional; falls back to CPU postproc).
ARG ME_BASE=me-sm120-test:latest
FROM ${ME_BASE}
ENV TORCH_CUDA_ARCH_LIST=12.0 OMP_NUM_THREADS=8 DEBIAN_FRONTEND=noninteractive
# open3d / opencv need GL + glib at import time even headless.
RUN apt-get update && apt-get install -y --no-install-recommends \
      libgl1-mesa-glx libglib2.0-0 libgomp1 wget git && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir \
      plotly open3d opencv-python-headless hydra-core "numpy<2.0" Pillow tqdm \
      iopath matplotlib scipy laspy lazrs
WORKDIR /workspace
RUN git clone --depth 1 https://github.com/ZiyuGuo99/SAM2Point.git
WORKDIR /workspace/SAM2Point
# SAM2 hiera-large checkpoint (~900 MB) from the official FAIR mirror.
RUN cd checkpoints && wget -q https://dl.fbaipublicfiles.com/segment_anything_2/072824/sam2_hiera_large.pt
ENV PYTHONPATH=/workspace/SAM2Point
