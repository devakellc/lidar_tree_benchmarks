# SegmentAnyTree (#M6) on Blackwell — extends the MinkowskiEngine sm_120 base.
# Build the base first:  bash gpu/build_me_sm120.sh
ARG ME_BASE=me-sm120-test
FROM ${ME_BASE}

ARG SAT_REF=main
ARG PYG_INDEX=https://data.pyg.org/whl/torch-2.7.0+cu128.html

ENV PYTHONUNBUFFERED=1 \
    SAT_HOME=/opt/segment-any-tree \
    PYTHONPATH=/opt/segment-any-tree

RUN apt-get update && apt-get install -y --no-install-recommends \
    unzip wget curl \
    && rm -rf /var/lib/apt/lists/*

# PyG stack for torch 2.7 / cu128 (cp311 in the pytorch devel image). The
# compiled companions come from the cu128 index; torch-geometric is pinned to
# 2.3.1 — the last release where Data.keys is a property (modern PyG made it a
# method, breaking torch_points3d's `set(data.keys)`). Verified to import on
# torch 2.7 and to carry SAT inference end-to-end (#M6 bring-up, 2026-06-09).
RUN pip install --no-cache-dir \
    pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
    -f ${PYG_INDEX} \
 && pip install --no-cache-dir "torch-geometric==2.3.1"

# torch-points-kernels pins numpy<1.20; patch and compile for sm_120.
RUN git clone --depth 1 https://github.com/nicolas-chaulet/torch-points-kernels.git /tmp/tpk \
 && sed -i 's/"numpy<1.20"/"numpy>=1.21"/' /tmp/tpk/setup.py \
 && cd /tmp/tpk \
 && FORCE_CUDA=1 pip install --no-cache-dir .

# hdbscan: imported by the PointGroup model. pyarrow: required by modern dask's
# dataframe, which the merge step (merge_pt_ss_is.py) uses. Both were missing and
# blocked inference/merge respectively (#M6 bring-up).
RUN pip install --no-cache-dir \
    hydra-core==1.3.2 omegaconf==2.3.0 h5py plyfile laspy \
    numba scikit-image gdown tqdm addict \
    matplotlib pandas tensorboard wandb torchnet pytorch-metric-learning \
    jaklas pykdtree dask hdbscan pyarrow

WORKDIR /opt
# Checkout (not --branch) so SAT_REF may be a branch, tag, OR a pinned commit SHA
# (gpu/SegmentAnyTree.commit). A full clone is needed to resolve an arbitrary SHA.
RUN git clone https://github.com/SmartForest-no/SegmentAnyTree.git segment-any-tree \
 && git -C segment-any-tree checkout ${SAT_REF}

WORKDIR ${SAT_HOME}

# Weights: prefer bind-mount from gpu/store/segmentanytree at run time; fetch
# a copy at build for import smoke if the mirror file is absent.
RUN mkdir -p model_file \
 && curl -fL --retry 3 \
    https://media.githubusercontent.com/media/SmartForest-no/SegmentAnyTree/main/model_file/PointGroup-PAPER.pt \
    -o model_file/PointGroup-PAPER.pt

COPY test_imports.py /opt/test_imports.py
CMD ["python", "/opt/test_imports.py"]
