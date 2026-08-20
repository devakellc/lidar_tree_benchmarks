# Vendored Treeiso (classical, unsupervised 3-D tree isolation) — #P5

`treeiso.py` and `cut_pursuit_L0.py` are **vendored verbatim** (MIT, see
[`LICENSE`](LICENSE)) from
[truebelief/artemis_treeiso](https://github.com/truebelief/artemis_treeiso)
(the Feb-2025 pure-Python build). Treeiso is the graph-cut (cut-pursuit)
unsupervised instance segmenter used as the non-learned 3-D baseline in
FGI-EMIT; its errors are decorrelated from CHM/local-maximum detectors, so it
adds diversity to the #P1 consensus pool.

Please cite:

> Xi, Z.; Hopkinson, C. 3D Graph-Based Individual-Tree Isolation (Treeiso) from
> Terrestrial Laser Scanning Point Clouds. *Remote Sens.* **2022**, 14, 6116.
> <https://doi.org/10.3390/rs14236116>
>
> Landrieu, L.; Obozinski, G. Cut Pursuit: Fast Algorithms to Learn Piecewise
> Constant Functions on General Weighted Graphs. *SIAM J. Imaging Sci.* 2017.

## `run_treeiso.py` (added here, not upstream)

A non-interactive runner: `python run_treeiso.py <input.laz> <output.laz>`. It
removes ground (Classification 2), centers, and runs the three Treeiso stages
(`init_segs` → `intermediate_segs` → `final_segs`) on **one** consistent
decimation, writing per-point instance ids to a `treeiso` int32 extra dim on the
vegetation points (absolute coordinates preserved). The upstream `main()` is
not used — it mismatches its two decimation index maps (`dec_inverse_idx` vs
`dec_inverse_idx2`), raising an IndexError; the single-decimation path here
avoids that and is appropriate for plot-sized ALS clips (the 0.05 m decimation
is a near-no-op at ALS density).

## Environment

CPU only (no GPU). Create a dedicated env and install the pinned deps:

```sh
conda create -n treeiso python=3.11 -y
conda run -n treeiso pip install "numpy<2.0" numpy-indexed laspy scipy \
  scikit-image PyMaxflow lazrs
```

`scripts/detect_treeiso_sweep.R` calls this runner per (plot, rung) via
`PYTHON=~/miniconda3/envs/treeiso/bin/python` (override with `PYTHON=`).
