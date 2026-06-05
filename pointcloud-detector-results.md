# Point-Cloud Detectors vs CHM-VWF at Native Density — Results

*Addresses [issue #6](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/6):
"Add a point-cloud detector arm (Li 2012) at native density vs CHM-VWF". The
density-ladder sweep ([density-ladder-sweep-results.md](density-ladder-sweep-results.md))
is CHM-VWF only. A 2.5-D canopy height model sees only the top surface, so the
one place high point density **should** pay off is a point-cloud detector able
to resolve sub-dominant apexes a CHM cannot. This run tests that directly. Code:
[scripts/detect_pc_sweep.R](scripts/detect_pc_sweep.R). Last run: 2026-06-05.*

---

## TL;DR

- At **native density** (measured 4.2–17.8 first-returns/m², median 11.5) across
  **46 plots** at SJER + SOAP + TEAK (699 live field stems), three point-cloud
  detectors were run on the **same** normalized clip as the CHM-VWF baseline and
  scored against the **same** stems by crown class.
- **Only Li 2012 moves understory recall meaningfully**: pooled
  intermediate+suppressed recall **0.257 vs 0.162** for CHM-VWF (**+0.095**),
  and overall recall 0.495 vs 0.419 (**+0.076**) — but it **pays in precision**
  (0.266 vs 0.336) and runs **~37× slower** than point-cloud `lmf`.
- **lidR `lmf`-on-points and lasR point `local_maximum` are byte-identical**
  (same allometry, same points): understory **0.190** (+0.029) but overall
  recall slightly **worse** (0.386, −0.033) — more apexes, mostly false.
- **The understory floor is real and largely occlusion.** The best detector
  still recovers only **27 of 105** intermediate+suppressed stems. No
  point-cloud method gets understory recall above ~0.26 at any site. Native
  point density does not buy back the stems hidden beneath the dominant canopy.

---

## Method

Per plot, one native clip is produced by `prepare_clip(..., rung=NA)` (clip +
TIN-normalize using the existing ground class; no decimation). The same clip
feeds all four detectors:

| Detector | What it is |
| --- | --- |
| `chm_vwf` | Baseline: `detect_lasr(file, res=0.5, a=0.10, frdens)` — pit-filled CHM + variable-window `local_maximum_raster`. |
| `lidr_lmf_pc` | `lidR::locate_trees(las, lmf(ws, hmin=2, shape="circular"))` on the **point cloud**. |
| `lidr_li2012` | `lidR::segment_trees(las, li2012(dt1=1.5, dt2=2, R=2, Zu=15, hmin=2, speed_up=10))`; apex = max-Z point per `treeID`. |
| `lasr_lmax_pc` | lasR point-based `local_maximum(ws, min_height=2)` (**not** `local_maximum_raster`). |

All windows use the **same** variable-window allometry as the CHM path:
`ws <- ws_factory(0.10)` (slope a=0.10, clamped [3,5]); `hmin`/`min_height` = 2;
`shape="circular"` for `lmf`. Scoring is identical to the sweep: `score_plot`
with `tol_xy=4`, plot core via `plot_half(plotType)` (tower ±20 m, distributed
±10 m), global nearest-distance greedy 1:1 matching with the height-consistency
gate. Plots kept if ≥6 live stems (mirrors `run_sweep.R`). Pooling is by summed
counts (`sum(TP)/sum(n_ref)`), never a mean of per-plot rates; per-class TP is
recovered as `round(rec_class × n_class)` (reuses `analyze_sweep.R`'s `pool()`).

## Pooled results (all sites, 46 plots, 699 stems)

Recall and per-crown-class recall by detector:

| Detector | n_det | recall | precision | F1 | dominant | codominant | intermediate | suppressed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `chm_vwf` | 798 | 0.419 | 0.336 | 0.373 | 0.639 | 0.387 | 0.156 | 0.222 |
| `lidr_lmf_pc` | 855 | 0.386 | 0.285 | 0.328 | 0.617 | 0.324 | 0.188 | 0.222 |
| `lidr_li2012` | 1189 | 0.495 | 0.266 | 0.346 | 0.722 | 0.447 | 0.250 | 0.333 |
| `lasr_lmax_pc` | 855 | 0.386 | 0.285 | 0.328 | 0.617 | 0.324 | 0.188 | 0.222 |

Class reference counts (pooled): dominant 180, codominant 367,
intermediate 96, suppressed 9.

## The understory question

Combined understory = intermediate + suppressed (105 stems pooled). This is the
metric issue #6 cares about — the stems a surface model structurally cannot see.

| Detector | understory recall | understory TP / 105 | Δ recall vs CHM-VWF | Δ understory vs CHM-VWF |
| --- | ---: | ---: | ---: | ---: |
| `chm_vwf` | 0.162 | 17 | — | — |
| `lidr_lmf_pc` | 0.190 | 20 | −0.033 | +0.029 |
| `lidr_li2012` | 0.257 | 27 | +0.076 | +0.095 |
| `lasr_lmax_pc` | 0.190 | 20 | −0.033 | +0.029 |

Per-site understory recall (the structure gradient SJER → SOAP → TEAK):

| Site | n_und | chm_vwf | lidr_lmf_pc | lidr_li2012 | lasr_lmax_pc |
| --- | ---: | ---: | ---: | ---: | ---: |
| SJER | 2 | 0.50 | 1.00 | 1.00 | 1.00 |
| SOAP | 45 | 0.200 | 0.267 | 0.333 | 0.267 |
| TEAK | 58 | 0.121 | 0.103 | 0.172 | 0.103 |

SJER has only 2 understory stems (open oak savanna — almost no sub-canopy), so
its 1.00 is two trees, not a trend. The load-bearing sites are SOAP (mixed
conifer) and TEAK (red fir): Li 2012 is the only detector that improves
understory recall at **both**, and even there it tops out at 0.333 / 0.172.

## Runtime (per native plot clip, 46 plots, 6 workers)

| Detector | median s | mean s | max s | total s |
| --- | ---: | ---: | ---: | ---: |
| `chm_vwf` | 0.15 | 0.17 | 0.37 | 8.0 |
| `lidr_lmf_pc` | 0.04 | 0.05 | 0.14 | 2.2 |
| `lidr_li2012` | 1.50 | 1.75 | 4.53 | 80.7 |
| `lasr_lmax_pc` | 1.23 | 1.44 | 5.44 | 66.2 |

Li 2012 costs **~37× the lidR point-`lmf`** mean time. The lasR point
`local_maximum` is also ~28× the lidR `lmf` here (its per-call `exec`/EPT setup
dominates on these small clips), yet returns the **identical** result set — so
it buys nothing over `lmf` on a per-plot basis.

## Why `lidr_lmf_pc` ≡ `lasr_lmax_pc`

Both apply the same `a·h+3` window (clamped [3,5]), `min_height=2`, to the same
normalized points; a local-maximum search over identical points with an
identical window returns identical apexes. Verified: identical apex count and TP
on **every** plot (max recall difference 0.000). They are two implementations of
the same operator, not two methods — included per the issue's request to check
**all** named detectors, but they should be read as one result.

## Verdict: occlusion floor

**The understory floor is genuinely occlusion-limited, not detector-limited.**
At native density the best detector (Li 2012, true 3D segmentation) lifts pooled
understory recall from 0.162 to 0.257 — it recovers **10 extra** of 105
sub-canopy stems by splitting clusters in the point cloud — but the remaining
~75% are simply not separable: their returns are sparse and entangled with the
dominant crowns above them. Point-`lmf` (= lasR `local_maximum`) adds 3 stems
and *loses* overall recall, because more apexes on the same surface are mostly
commission. The 2.5-D CHM is not the binding constraint for these stems; canopy
occlusion is. Li 2012's gain is real but small and bought with precision (530 →
873 pooled core false positives) and a ~37× runtime hit — a wall-to-wall
deployment trade-off, not a free lunch.

**Bottom line for the sweep:** adding a point-cloud arm does not overturn the
sweep's finding that detection accuracy is structure- and occlusion-bound rather
than density-bound. It refines it: high density + a 3D segmenter recovers a
*thin* slice of understory the CHM misses, but the floor is set by what the
LiDAR can physically see under the canopy.
