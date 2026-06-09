# ForestFormer3D inference image for RTX 5090 (sm_120), torch 2.7.0 / cu128.
# Consolidates the working sequence; ff3d_freeze.txt is the authoritative version
# lock. FROM the proven MinkowskiEngine-sm120 base (see ../minkowski-sm120), which
# is pytorch/pytorch:2.7.0-cuda12.8-cudnn9-devel + ME 0.5.4 built for arch 12.0.
FROM me-sm120-test:latest

ENV FORCE_CUDA=1 MMCV_WITH_OPS=1 TORCH_CUDA_ARCH_LIST=12.0 OMP_NUM_THREADS=8
ENV PIP="pip install --no-cache-dir --no-deps"

# open3d / opencv runtime libs (slim base lacks them)
RUN apt-get update && apt-get install -y --no-install-recommends \
      libgl1 libgomp1 libx11-6 libxext6 libxrender1 libsm6 libglib2.0-0 \
      libusb-1.0-0 git ninja-build cmake build-essential \
    && rm -rf /var/lib/apt/lists/*

# numpy<2 track + sparse-conv backends (spconv backbone). cumm-cu128 MUST come
# from rathaROG (PyPI is ABI-stale); ccimport/pccm must stay current for cumm.
RUN $PIP "numpy==1.26.4" "numba==0.59.1" "llvmlite==0.42.0" && \
    pip install --no-cache-dir --no-deps --index-strategy unsafe-best-match \
      --extra-index-url https://ratharog.github.io/cumm-spconv/ \
      "cumm-cu128==0.9.1" "spconv-cu128==2.4.1" && \
    $PIP -U "ccimport==0.4.4" "pccm==0.4.16"

# OpenMMLab stack (FF3D pins), pure-python, --no-deps
RUN $PIP "mmengine==0.7.3" "mmdet==3.0.0" "mmsegmentation==1.0.0" && \
    $PIP "git+https://github.com/open-mmlab/mmdetection3d.git@22aaa47fdb53ce1870ff92cb7e3f96ae38d17f61"

# mmcv 2.0.0 _ext from source — torch 2.7 headers need C++17 (mmcv hardcodes c++14)
RUN git clone --depth 1 --branch v2.0.0 https://github.com/open-mmlab/mmcv.git /opt/mmcv && \
    sed -i 's/c++14/c++17/g' /opt/mmcv/setup.py && \
    cd /opt/mmcv && MAX_JOBS="$(nproc)" pip install --no-cache-dir --no-deps -v .

# compiled geometry ops for sm_120
RUN $PIP --no-build-isolation "torch-scatter==2.1.2" "torch-cluster==1.6.3" \
      "torch-points-kernels==0.7.0"

# support + viz + dataset-SDK deps that mmdet3d/oneformer3d import eagerly
RUN $PIP "opencv-python-headless==4.7.0.72" "open3d==0.19.0" plotly tenacity \
      "laspy[lazrs]" && \
    pip install --no-cache-dir dash && $PIP "numpy==1.26.4" && \
    $PIP "lyft-dataset-sdk==0.0.8" "nuscenes-devkit==1.1.10" Pillow descartes fire \
      addict yapf termcolor "Shapely==1.8.5" "scipy==1.10.1" terminaltables \
      pyquaternion "pandas==2.0.1" matplotlib scikit-learn trimesh plyfile \
      tensorboard pycocotools cachetools

# Runtime: mount the FF3D repo at /workspace with replace_mmdetection_files/ and
# apply the four code patches (ff3d_repo.patch) + the mmengine/mmdet3d file swaps:
#   cp replace_mmdetection_files/loops.py       <site>/mmengine/runner/
#   cp replace_mmdetection_files/base_model.py  <site>/mmengine/model/base_model/
#   cp replace_mmdetection_files/transforms_3d.py <site>/mmdet3d/datasets/transforms/
#   git -C /workspace apply ff3d_repo.patch     # test.py + config (drop tensorboard)
WORKDIR /workspace
CMD ["bash"]
