# NEON Deep-Model Benchmark — Design & Issue Decomposition

**Status.** Approved design (brainstorming output). Ready for implementation
planning. Last updated 2026-06-07.

**Scope.** Extend the existing density-first CHM-VWF benchmark to evaluate the
3D tree-segmentation models from
[`docs/deep-research-report.md`](../../deep-research-report.md) on the NEON
dataset — zero-shot, scored as detections through the existing greedy-match
harness, across the full density ladder, SOAP first.

---

## 1. Locked decisions

These framing decisions were settled during brainstorming and bound everything
below:

| Decision | Choice | Consequence |
|---|---|---|
| Model set | Whatever is *actually runnable* (feasibility-gated) | A triage gate decides the set; see §2. |
| Compute | Local CUDA GPU available | Deep models are on the table; single GPU ⇒ serial runs. |
| Scoring basis | **Reduce predictions to detections**, reuse the existing harness | Each model collapses to a detection table `(x, y, apex height)` scored by `score_plot`/`greedy_match` unchanged. Crown diameter is a **separate optional side metric** (the scorer ignores it). Fully comparable to the current CHM-VWF / Li 2012 arms. No new ground truth. |
| Model application | **Zero-shot pretrained transfer** | No fine-tuning (NEON has no per-point instance labels). |
| Density axis | **Full density ladder** (native ~20 / 8 / 4 / 2 / 1 pts/m²) | Directly tests the report's density-robustness claims against the repo's decimation thesis. |
| Sites | **SOAP first**, others follow | SOAP 2021 is the one site that natively clears 8 pts/m². |
| Sequencing | **Hybrid B→A** | Walking skeleton on one model to de-risk, extract the shared bridge, then parallel fan-out. |

### Density framing (corrected)

The rung labels are **all-return** point density (the unit `homogenize()`
targets). Measured SOAP 2021 native is **~20 pts/m² all-return / ~11.8
first-returns/m²** — pulse density is roughly half the all-return number. Deep
models voxelize **all returns**, so the operative figure for them is the
all-return ~20.

This matters: at native ~20 the deep models are **within or near their tested
envelope** (above ForestFormer3D's and SegmentAnyTree's floor of 10 pts/m²,
inside AMS3D's 4–20 comfort zone). Feasibility risk is concentrated at the
decimated **4 / 2 / 1** rungs. A second risk survives even at matched 20 pts/m²
and is **structural, not density-driven**: NEON discrete-return ALS is
crown-surface-dominated with few sub-canopy returns, whereas these models
trained on dense ULS decimated down. So "20 pts/m² of NEON" is *easier on
density, harder on vertical structure* — and that structural gap is the crux of
the report's understory-advantage claim.

---

## 2. Feasibility verdict (the runnable set)

Each profile below was produced by a research fan-out and **adversarially
verified against primary sources**. Full per-model profiles, env pins, and
source URLs live in §8.

### Run (the fan-out)

| Model | Family | Code / weights | Zero-shot at native ~20 | Rung scope | Role |
|---|---|---|---|---|---|
| **AMS3D** (`crownsegmentr`, CRAN R+C++) | classical mean-shift | public; no weights (classical) | **likely** (comfort 4–20) | full ladder | **Walking skeleton** — pure R, no GPU/Docker, normalized input matches `prepare_clip`. |
| **SegmentAnyTree** | PointGroup 3D CNN | public (BSD); weights inside 8.6 GB Docker image | risky (floor 10; in-envelope at 20) | full ladder | Its headline claim *is* density-robustness — the low rungs test it directly; cheap to run. |
| **TreeisoNet-ALS** (`NRCan/TreeAIBox`) | SegFormer voxel | public; `.pth` on a personal server (mirror!) | risky (domain: reclamation regrowth) | full ladder | Risk is domain not density; all rungs informative. Run coarsened voxels + both checkpoints. |
| **ForestFormer3D** | transformer panoptic | public (CC-BY-NC/GPL tangle); Zenodo checkpoint | risky (floor 10; in-envelope at 20) | **native + 8 only** | Report's #1 leader; native clears its floor. Heaviest env ⇒ don't burn collapse rungs. |
| **lidRplugins** (`lmfauto`/`multichm`/`ptree`) | classical CHM + point | public; classical | yes (true zero-shot) | full ladder | **Fair LiDAR-native, density-responsive competitor** to CHM-VWF. |

All scored against the existing **CHM-VWF** and **Li 2012** baselines.

### Defer (revisit in Phase 5)

- **ForAINet** — needs stem/branch returns NEON lacks; zero-shot verdict
  **refuted**. Very-high env effort (torch 1.9 / CUDA 11.1 / MinkowskiEngine).
- **HFC** — public repo is an **empty stub**; code requires a signed email
  agreement (~1 week lead). Start that email early if we want it.

### Drop (architecturally impossible on NEON ALS)

- **TreeLearn** — trunk-clustering at ~1 pt/(0.1 m)³; authors explicitly
  exclude ALS.
- **Dersch graph-cut + stem detection** — **no public code**, and needs
  trunk returns NEON ALS does not capture under closed canopy.

### DeepForest (optional, Phase 5)

RGB RetinaNet box detector — **density-invariant** (never reads the point
cloud), **near-in-domain on NEON** (the prebuilt model was trained on 22 NEON
sites, so *not* a clean zero-shot transfer), and needs **NEON AOP RGB + the CHM
for height** (boxes carry no `z`). Valuable as a flat reference line, but a
different kind of arm; kept optional with the confound flagged loudly.

---

## 3. Architecture

The benchmark is **one new data path bolted onto the existing scorer, which
stays untouched.** Everything funnels into `score_plot` / `greedy_match`
(both already in `sweep_lib.R`) exactly as the Li 2012 arm already does.
Pooling is **not** yet a shared function — `pool()` is copy-pasted in four
scripts (`analyze_sweep.R`, `detect_pc_sweep.R`, `temporal_sensitivity.R`,
`calval_split.R`), so #B2 *canonicalizes* it rather than reusing it. Three
layers.

### Data flow

```text
NEON SOAP LiDAR  +  field stems (existing ground truth)
   │
   │  ① frozen-clip provider  (decimate raw ONCE w/ stable seed keyed by site/plot/rung; reused by ALL arms)
   │      ├─ clip_rawground.laz    → SegmentAnyTree, TreeisoNet, FF3D       (ground retained)
   │      ├─ clip_normalized.laz   → CHM-VWF, Li2012, AMS3D, lidRplugins    (derived from the SAME decimated set)
   │      ├─ ground_dtm.tif (TIN)  → raw→height transform for raw-ground arms' apex z
   │      └─ manifest.json         → seed, point counts, pdens/frdens, source tile ids
   │
   ├─ R-native arms ───────────────────────────────► read clip in R, emit instances
   └─ GPU arms ─► ② LAZ→PLY ─► [Docker per model] ─► instance PLY/LAS
                                                            │
   per-instance point groups ──── ③ reducer ──► data.frame(x, y, z) ──┘   (lowercase, UTM, normalized z)
                                                            │
                                          ④ conformance harness (assert contract before any real run)
                                                            │
                       score_plot → greedy_match  (UNCHANGED)  ──► long-form CSV (model × rung × plot × class)
                                                            │
                       ⑤ pool() + per-(plot,rung) equal-set guard ──► per-class / per-height-band tables + curves
```

### Layer 1 — the bridge (core deliverable, proven by the AMS3D skeleton)

- **① Frozen-clip provider** — the single most important new component.
  Decimation today is *non-deterministic* (`decimate_points(homogenize())` with
  no `set.seed`), and `prepare_clip()` decimates inline, normalizes, and writes
  **only** the normalized LAZ ([`sweep_lib.R:90`](../../../scripts/sweep_lib.R)).
  **Acceptance criterion:** decimate the raw clip **exactly once**, with a stable
  seed keyed by `site/plot/rung`, then derive every variant from that one
  decimated raw point set — never a second stochastic decimation. Feed the
  **identical** point subset to every arm or cross-arm deltas are partly RNG
  noise. It persists:
  - `clip_rawground.laz` — ground retained, for the deep trunk/offset arms that
    need ground as a semantic class.
  - `clip_normalized.laz` — derived from the *same* decimated set, for the
    CHM/AMS3D/lidRplugins arms. Feeding the wrong variant is the **#1
    silent-wrong-result bug** (a trunk model fed a flattened cloud looks like a
    "density" failure when it was really starved of ground).
  - `ground_dtm.tif` (the TIN/DTM) — the raw→height-above-ground transform, so
    apexes computed on `clip_rawground.laz` can be normalized **back to height**
    before scoring (see ③); guards against absolute elevations leaking into the
    height gate as `z`.
  - `manifest.json` — seed, per-variant point counts, `pdens`/`frdens` (native
    `pdens` drives the no-upsampling guard), and source-tile identity, so a run
    is reproducible and auditable.
- **③ Reducer** — the universal collapse `instance point-group → (x, y, z)`.
  Output is **exactly `x,y,z`** (the only columns the scorer reads); one
  **standardized apex rule** across all arms. For raw-ground arms, the apex `z`
  is converted to **height above ground** via the persisted `ground_dtm.tif`
  before it enters the scorer — never an absolute elevation. Crown diameter is
  **not** part of this contract: it is a separate optional side metric (§ Layer
  3 / #B2), so the reducer stays a clean `data.frame(x,y,z)`.
- **④ Conformance harness** — a synthetic one-tree clip; every arm must return
  a base `data.frame` with exactly lowercase `x, y, z`, UTM coords, normalized
  height `z`, non-NULL on empty, *before* any real run. (`sf::st_coordinates`
  returns capital `X/Y/Z` that must be renamed — easy to get wrong, so assert.)

### Layer 2 — GPU/infra (unblocks the heavy arms; never touches the scorer)

- **Containerized per-model runner** — the four deep stacks have *mutually
  incompatible* CUDA / torch / MinkowskiEngine / spconv pins and cannot share an
  env. One Docker image per arm + a shared bind-mount I/O dir + a thin R
  `system2()` wrapper (drop clip in → read CSV out). Single GPU ⇒ arms run
  serially.
- **Weights/image mirror + checksums** — TreeisoNet's weights live on one
  researcher's personal site, ForAINet's on a single Dropbox link,
  SegmentAnyTree's only inside an 8.6 GB image. Mirror + SHA256 + record source
  URLs on day one or the benchmark is not reproducible after link-rot.

### Layer 3 — scoring (reuse + four additions)

- Reuse `greedy_match` / `score_plot` **unchanged** (both already in
  `sweep_lib.R`): ±20 m tower / ±10 m distributed cores, the `[0.5·az, az+8]`
  height-consistency gate, per-crown-class and per-height-band recall.
- **Canonicalize `pool()`** — it is *not* in `sweep_lib.R`; it is duplicated in
  `analyze_sweep.R`, `detect_pc_sweep.R`, `temporal_sensitivity.R`, and
  `calval_split.R` (sum-TP / sum-n_ref). #B2 extracts **one** canonical pooler
  (with the equal-set guard below) into the shared lib so the benchmark doesn't
  create a fifth subtly-different copy; the analysis scripts can migrate to it
  opportunistically.
- **Per-(plot,rung) equal-set guard** — arms fail at *different* rungs (one
  OOMs at native, another collapses at 1), so the drop-set must be computed per
  `(plot, rung)`, not per plot, or a rung's cross-arm table silently compares
  different plot populations.
- **Per-crown-class + per-height-band recall is the primary table**, not pooled
  F1 — at these densities pooled F1 is dominated by overstory and hides the
  understory signal that is the entire scientific point.
- **Zero-shot protocol ledger** — record every non-default knob (AMS3D
  allometry ratios, TreeisoNet voxel resolution, any HFC-style intensity
  threshold). These are *de facto* tuning and must be declared, not hidden.

**Crown diameter — separate optional side metric (not in the detection
contract).** `score_plot` reads only `x,y,z`. Crown-diameter RMSE is scored
exactly as `crown_metrics_sweep.R` already does it — joining field
`maxCrownDiameter`/`ninetyCrownDiameter` from the NEON `vst` rds and matching on
shared tree-tops — as an **add-on table** keyed off the same detections, with
the per-crown point-count floor applied there (a caliper over 3 points at
1 pt/m² is garbage). It never enters the primary reducer or the greedy match.

### Reducer output contract (exact, from the repo)

`score_plot(stems, det, tol_xy, core_cx, core_cy, core_half)` expects `det` to
be a **base `data.frame`** with exactly:

- `x` — treetop easting, plot UTM (EPSG:32610/32611), same CRS as `stems$E`.
- `y` — treetop northing, same CRS as `stems$N`.
- `z` — apex **height above ground** in metres (normalized Z, **not**
  elevation), same datum as `stems$height`; used for both the height gate and
  height RMSE.
- Empty case returns `data.frame(x=numeric(), y=numeric(), z=numeric())`, never
  `NULL`. No `treeID`/`crown_class` column is read — crown class comes from the
  **stem side** only.

Coordinates **must be metric UTM**, not EPSG:3857 (Web Mercator inflates
distance ~1.32× at this latitude and breaks `tol_xy`). The integration template
to copy is `detect_pc_sweep.R`.

---

## 4. Issue decomposition

Sequencing follows the approved **B→A** shape. Each issue is one script (+
shared lib) and, where it produces findings, one result-doc section — matching
how issues #3–#8 are already tracked. Research-surfaced gaps are folded in as
**acceptance criteria**, not separate issues.

### Phase 0 — Gate

- **#A0 · Triage + standing protocol.** Record the run/defer/drop verdict
  table, the zero-shot protocol (which per-arm knobs may be set from prior
  literature vs. forbidden from SOAP-tuning), and the weights-mirroring policy.
  Deliverable: `docs/model-benchmark-plan.md`. *Deps: none.*

### Phase 1 — Walking skeleton (B), then extract the layer (A)

- **#B1 · AMS3D arm end-to-end on SOAP** (native → full ladder).
  `detect_ams3d_sweep.R` patterned on `detect_pc_sweep.R`. Proves the whole
  bridge in-language, zero GPU. *Absorbs:* reducer v1 (`x,y,z` only), conformance
  harness v1, AMS3D allometry recorded in the ledger, and crown-diameter scoring
  as a **separate** side table (not in the reducer). *Deps: #A0.*
  **← de-risking spike.**
- **#B2 · Extract the shared bridge** into `model_bench_lib.R`:
  - **Frozen-clip provider** — decimate raw **once** per `(site, plot, rung)`
    with a stable keyed seed; persist `clip_rawground.laz`, the
    `clip_normalized.laz` derived from that *same* decimated set, the
    `ground_dtm.tif` raw→height transform, and a `manifest.json` (seed, counts,
    `pdens`/`frdens`, source tiles). Native `pdens` drives the no-upsampling
    guard.
  - **Universal reducer** — `instance → data.frame(x,y,z)`; one apex rule;
    raw-ground apexes normalized to height-above-ground via `ground_dtm.tif`.
  - **Conformance harness** — synthetic one-tree contract assertion.
  - **Canonical `pool()` + per-(plot,rung) equal-set guard** — one extracted
    pooler (today duplicated in four scripts), plus per-class/height-band
    primary tables and the zero-shot ledger.
  - **Crown-diameter side metric** — the optional add-on table (with the
    point-count floor), kept out of the detection path.

  AMS3D re-runs unchanged on top as the regression check. *Deps: #B1.*
  **← A-layer hinge; all arms depend on this.**

### Phase 2 — GPU/infra scaffolding (parallelizable)

- **#I3 · Containerized runner contract.** One Docker image per GPU arm +
  shared bind-mount + R `system2()` wrapper. A no-op "echo a fixed CSV" arm
  proves the R↔container↔scorer handshake before any model. *Deps: #B2.*
- **#I4 · I/O bridge.** LAZ→PLY preserving intensity/return · instance
  PLY/LAS → `data.frame(x,y,z)` parser · CRS/units round-trip acceptance test
  (known stem → apex within tol of UTM; `z` = normalized height). *Deps: #B2.*
- **#I5 · Weights/image mirror + checksums.** SegmentAnyTree image ·
  TreeisoNet `.pth` (treeloc+treeoff, both ALS-reclamation **and**
  UAV-mixedwood) · FF3D Zenodo checkpoint → project store, SHA256 + source URL.
  *Deps: none — start at t0, parallel to Phase 1.*

### Phase 3 — GPU model arms (fan-out; each *Deps: #B2, #I3, #I4, #I5*)

- **#M6 · SegmentAnyTree arm**, full ladder. Docker; raw-with-ground; smoke on
  SOAP native → ladder.
- **#M7 · TreeisoNet-ALS arm**, full ladder. **Coarsened voxels (~0.8 m /
  2.0 m)**, not native 10 cm; run **both** checkpoints (reclamation +
  UAV-mixedwood) as variants.
- **#M8 · ForestFormer3D arm**, **native + 8 only**. Heaviest env (their Docker
  as-is); cross-block instance dedup verified before the reducer.

### Phase 4 — Competitor + headline analysis

- **#C9 · lidRplugins competitor arm**, full ladder. `lmfauto`/`multichm`
  (CHM-native) + `ptree` (point-native); same `ws_factory`/`hmin` discipline as
  CHM-VWF for a fair fight. R-native ⇒ *Deps: #B2 only* — can land alongside
  Phase 2.
- **#R10 · Cross-model density-ladder result doc.**
  `results/model-benchmark-results.md`: per-class + per-height-band recall as
  the primary tables; density-robustness curves; head-to-head vs CHM-VWF + Li
  2012; zero-shot-ledger appendix; decimation-as-simulation and
  discrete-return-vs-ULS caveats. Written incrementally as arms land.
  *Deps: the arms.*

### Phase 5 — Extensions (deferred / optional)

- **#E11 · Add SJER + TEAK** (note their native is already ~4–6 pts/m²).
- **#E12 · DeepForest RGB reference line** (AOP RGB + CHM-for-height; flag the
  NEON-trained confound).
- **#E13 · Revisit deferred** — ForAINet; **start HFC's signed-agreement email
  early** (~1 week lead).
- **#E14 · Optional external point-IoU side-arm** on decimated FOR-instance /
  FGI-EMIT, clearly labeled non-NEON.

### Critical path

```text
#A0 → #B1 → #B2 → (#I3 ∥ #I4 ∥ #I5) → (#M6 ∥ #M7 ∥ #M8) → #R10
                   #I5 starts at t0 (parallel)
                   #C9 joins off #B2
```

First real cross-arm result lands after **#B1 + #B2 + #C9** (AMS3D +
lidRplugins + existing baselines; #C9 depends on the extracted bridge #B2) —
before any Docker work. **10 core issues (#A0–#R10) + 4 deferred
(#E11–#E14).**

---

## 5. Gap → issue mapping

Research-surfaced risks and where each is handled, so none is lost:

| Gap / risk | Handled in |
|---|---|
| Normalized-vs-raw-ground input mismatch (silent-wrong-result #1) | #B2 frozen-clip two variants from one decimated set |
| Non-deterministic decimation breaks cross-arm comparability | #B2 single keyed-seed decimation + manifest |
| Raw-ground apex elevations leaking into the height gate as `z` | #B2 persisted `ground_dtm.tif` raw→height transform |
| Instance→apex reduction ambiguity | #B1/#B2 standardized apex rule |
| Crown-diameter from sparse points unstable | #B2 crown-diameter side metric (point-count floor, out of detection path) |
| `pool()` duplicated in 4 scripts (no canonical pooler) | #B2 extract one canonical `pool()` |
| Per-model env incompatibility | #I3 containerized runner |
| Single-GPU runtime / OOM | #I3 (serial) + per-arm chunk levers in #M6–#M8 |
| Weights link-rot / single points of failure | #I5 mirror + checksum |
| R↔Python format conversion, field preservation | #I4 I/O bridge |
| CRS / units handoff (UTM, normalized z) | #I4 round-trip test + #B2 `ground_dtm.tif` |
| Cylinder/block tiling + cross-block instance dedup | #M8 (and #M6/#M7 as needed) |
| Per-(plot,rung) equal-set guard | #B2 canonical pooler |
| Understory bias hidden by pooled F1 | #B2 per-class/height-band primary tables |
| Zero-shot leakage via per-arm tuning | #A0 protocol + #B2 ledger |
| Missing fair LiDAR-native competitor | #C9 lidRplugins arm |
| Scorer-contract conformance | #B2 conformance harness |

---

## 6. Out of scope / non-goals

- **No fine-tuning** of any model (zero-shot only).
- **No new per-point ground truth** on NEON; native point-IoU/coverage/mIoU
  scoring is *not* done on NEON (only via the optional external side-arm #E14).
- **No wall-to-wall production pipeline** — plot-scale clips only, matching the
  existing sweep.
- **No engine runtime benchmark** — wall-clock is recorded for budgeting, not
  as a headline metric.

## 7. Open questions / caveats carried forward

- **Compute budget.** Single GPU + full ladder × N SOAP plots × deep inference
  (FF3D wants A100-class; TreeLearn-style peaks aside, expect multi-GB VRAM and
  long serial runs). #I3 should record a per-arm VRAM/wall-clock estimate that
  can gate the SOAP-first decision.
- **Expected outcome is honest about failure.** The runnable deep arms are
  predicted to lose to CHM-VWF at 4/2/1 pts/m²; the benchmark's value is
  *quantifying the degradation curve and the understory gap*, not crowning a
  winner. #R10 must frame it that way.
- **Decimation ≠ native low-density.** Carry the existing repo caveat; the
  native QL2 cross-check (#4) is the precedent.

---

## 8. References

## Repo internals

- [`docs/deep-research-report.md`](../../deep-research-report.md) — the model
  source.
- [`docs/dataset-research-and-sweep-plan.md`](../../dataset-research-and-sweep-plan.md)
  — NEON density-ladder rationale.
- `scripts/sweep_lib.R` (`greedy_match`, `plot_half`, `prepare_clip`,
  `score_plot`), `scripts/detect_pc_sweep.R` (the arm template),
  `scripts/crown_metrics_sweep.R` (crown-diameter scoring),
  `scripts/neon_ground_truth.R` (stem schema, crown class).

## Models (official code / weights, verified)

- ForestFormer3D — <https://github.com/bxiang233/ForestFormer3d>;
  inference <https://github.com/bxiang233/FF3D_inference>; weights
  <https://zenodo.org/records/16742708>; ICCV 2025; arXiv 2506.16991.
- SegmentAnyTree — <https://github.com/SmartForest-no/SegmentAnyTree>; Docker
  `maciekwielgosz/segment-any-tree`; arXiv 2401.15739; RSE 2024.
- TreeisoNet / TreeAIBox — <https://github.com/NRCan/TreeAIBox>; weights
  <https://github.com/NRCan/TreeAIBox/releases/tag/v1.0>; Sci. Remote Sensing
  2025 (S266739322500002X).
- ForAINet — <https://github.com/prs-eth/ForAINet>; arXiv 2312.15084.
- HFC (ITS-HFC) — <https://github.com/Geo3DSmart/ITS-HFC> (agreement-gated);
  T&F 10.1080/17538947.2024.2356124.
- TreeLearn — <https://github.com/ecker-lab/TreeLearn>; arXiv 2309.08471.
- AMS3D — `crownsegmentr`
  <https://cran.r-project.org/web/packages/crownsegmentr/>;
  <https://github.com/nikoknapp/MeanShiftR>; Ferraz et al. 2016,
  doi:10.1016/j.rse.2016.05.028.
- DeepForest — <https://deepforest.readthedocs.io>;
  <https://github.com/weecology/DeepForest>; Weinstein et al. 2020,
  doi:10.1111/2041-210X.13472.
- lidRplugins — <https://github.com/Jean-Romain/lidRplugins>
  (`lmfauto`, `multichm`, `ptree`).
- Dersch et al. 2021 (no public code) — ISPRS J. P&RS,
  doi:10.1016/j.isprsjprs.2020.10.011.

## External labelled datasets (optional point-IoU side-arm)

- FOR-instance — <https://zenodo.org/records/8287792> (CC-BY-4.0).
- FOR-instanceV2 — <https://zenodo.org/records/16742708> (GPL-3.0).
- FGI-EMIT — <https://doi.org/10.5281/zenodo.19351234> (CC-BY-NC-SA-4.0,
  released 2026-04-11).
