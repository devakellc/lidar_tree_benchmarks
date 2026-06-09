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

# numpy<2 track + sparse-conv backends (spconv backbone). cumm-cu128/spconv-cu128
# MUST come from rathaROG (PyPI is ABI-stale), so pin it as the primary --index-url
# (--no-deps means only those two wheels are fetched); ccimport/pccm must stay
# current for cumm. (--index-strategy is a uv-only flag; plain pip rejects it.)
RUN $PIP "numpy==1.26.4" "numba==0.59.1" "llvmlite==0.42.0" && \
    pip install --no-cache-dir --no-deps \
      --index-url https://ratharog.github.io/cumm-spconv/ \
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

# support + viz + dataset-SDK deps that mmdet3d/oneformer3d import eagerly.
# Pin to ff3d_freeze.txt: unpinned, these resolve too-new and break the mm-stack
# at runtime (e.g. yapf>=0.40 drops FormatCode(verify=) that mmengine 0.7.3's
# cfg.pretty_text calls). dash stays unpinned (not in the hot path).
RUN $PIP "opencv-python-headless==4.7.0.72" "open3d==0.19.0" \
      "plotly==6.8.0" "tenacity==9.1.4" "laspy[lazrs]" && \
    pip install --no-cache-dir dash && $PIP "numpy==1.26.4" && \
    $PIP "lyft-dataset-sdk==0.0.8" "nuscenes-devkit==1.1.10" "Pillow==11.0.0" \
      "descartes==1.1.0" "fire==0.7.1" "addict==2.4.0" "yapf==0.33.0" \
      "termcolor==2.3.0" "Shapely==1.8.5" "scipy==1.10.1" terminaltables \
      pyquaternion "pandas==2.0.1" "matplotlib==3.5.2" "scikit-learn==1.2.2" \
      "trimesh==3.21.6" "plyfile==1.0.2" "tensorboard==2.15.1" \
      "pycocotools==2.0.6" "cachetools==5.3.0"

# Transitive deps the --no-deps installs above skip but that fail at IMPORT:
# pybind11 + nvidia-arch (spconv/cumm), rich (mmengine + whole mm-stack), and the
# leaf deps of matplotlib/sklearn/laspy. Versions per ff3d_freeze.txt. Installed
# WITH deps so rich's own deps (markdown-it-py) come along; none touch the pinned
# GPU wheels.
RUN pip install --no-cache-dir \
      "pybind11==2.10.4" "nvidia-arch==7.1.0" "rich==13.3.5" "joblib==1.2.0" \
      "threadpoolctl==3.1.0" "cycler==0.11.0" "kiwisolver==1.4.4" \
      "pyparsing==3.0.9" "lark==1.1.5" "lazrs==0.8.1" "packaging==23.1" \
      "portalocker==3.2.0" "zipp==3.17.0" "python-dateutil==2.8.2"

# Bake the mmengine/mmdet3d site-packages swaps (FF3D's replace_mmdetection_files)
# into the image so they are present for every --rm run. SITE is conda's
# site-packages. ff3d_repo.patch (which edits the mounted REPO, not site-packages)
# is applied at container start by ff3d_entry.sh — see the FF3D arm (#M8).
COPY replace_mmdetection_files/ /opt/replace_mmdetection_files/
RUN SITE="$(python -c 'import site; print(site.getsitepackages()[0])')" && \
    cp /opt/replace_mmdetection_files/loops.py        "$SITE/mmengine/runner/loops.py" && \
    cp /opt/replace_mmdetection_files/base_model.py   "$SITE/mmengine/model/base_model/base_model.py" && \
    cp /opt/replace_mmdetection_files/transforms_3d.py "$SITE/mmdet3d/datasets/transforms/transforms_3d.py"

# Runtime: the #M8 arm mounts the FF3D repo + weights and runs ff3d_entry.sh,
# which applies ff3d_repo.patch (tools/test.py spconv-permute + tensorboard drop)
# to the mounted repo idempotently. The site-packages swaps above are already baked.
WORKDIR /workspace
CMD ["bash"]
