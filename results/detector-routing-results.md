# Per-cell detector routing (#P2)

The repo's only routing knob today is measured density (the `frdens`/`pdens`
guards in `sweep_lib`), and the SJER→SOAP→TEAK structure gradient is described
only descriptively. But the density ladder already shows the crossovers a router
would key on. This study quantifies the achievable meta-pipeline by **selecting
a detector per cell** from cheap deploy-time features instead of shipping one
fixed arm — and asks how much of that headroom a simple learned policy actually
captures.

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/route_detectors.R SITES=SOAP,SJER,TEAK
# -> work/neon/<SITE>/router_policy.csv (per-cell features + oracle label + held-out prediction)
```

## What this is

The per-cell × per-arm scored ladder is assembled from the per-arm `*_results.csv`
(CHM-VWF + multichm at the canonical density-derived `chm_res`, `vwf_a = 0.10`;
SegmentAnyTree, native Li 2012, native+8 ForestFormer3D as scored), restricted to
the **equal set** where the three laddered arms are all present (229 cells, three
sites, five density rungs). Deploy-time structure features are computed over each
cell's `clip_normalized.laz` with `lidR` — `rumple_index`, canopy cover, height
CV, gap fraction, mean canopy height — plus measured `frdens`/`pdens`; **no field
data**. Each cell is labelled with its argmax-F1 arm and an interpretable `rpart`
CART is fit and evaluated **leave-one-plot-out** (a plot's rungs are not
independent). Four policies are scored by the canonical `pool()` (sum counts,
never average rates): every single arm, the fixed best-single arm, the **oracle**
(per-cell argmax), and the **learned** router (held-out predictions). The routing
helpers `oracle_pick`/`select_policy_rows` (in
[`route_lib.R`](../scripts/route_lib.R)) are unit-tested
(`tests/testthat/test-route-lib.R`).

## Generated tables

### Routing policies (pooled, equal-set ladder, 229 cells)

| policy | n_ref | F1 | recall | understory recall |
|---|--:|--:|--:|--:|
| single: chm_vwf | 3489 | 0.353 | 0.314 | 0.133 |
| single: forestformer3d | 699 | 0.279 | 0.369 | 0.190 |
| single: multichm | 3489 | 0.417 | 0.550 | 0.377 |
| single: segmentanytree | 3489 | 0.373 | 0.404 | 0.238 |
| **FIXED-BEST (multichm)** | 3489 | 0.417 | 0.550 | 0.377 |
| LEARNED router (LOPO) | 3489 | 0.420 | 0.529 | 0.345 |
| **ORACLE (per-cell best)** | 3489 | **0.479** | 0.563 | **0.404** |

Learned vs fixed-best ΔF1 = **+0.002**; oracle headroom ΔF1 = **+0.062**.

### Per-arm pooled F1 by density rung (the crossover the oracle exploits)

| rung | chm_vwf | multichm | segmentanytree | oracle picks |
|---|--:|--:|--:|---|
| native | 0.365 | 0.412 | **0.443** | segmentanytree > multichm |
| 8 | 0.346 | 0.414 | **0.425** | segmentanytree > multichm |
| 4 | 0.367 | **0.430** | 0.397 | multichm > chm_vwf |
| 2 | 0.338 | **0.421** | 0.308 | multichm > chm_vwf |
| 1 | 0.343 | **0.409** | 0.102 | multichm > chm_vwf |

### Learned routing table (CART, fit on all cells)

The tree keys on density and structure exactly as hypothesized — high `pdens`
routes to SegmentAnyTree (understory regime), low `pdens` to multichm/chm_vwf:

```text
root: multichm
├─ pdens < 4.5
│  ├─ mean_ht < 14.2
│  │  ├─ mean_ht < 9.1 → chm_vwf  (tall-gap / short stands)
│  │  └─ mean_ht ≥ 9.1 → multichm (rumple ≥ 6 → chm_vwf)
│  └─ mean_ht ≥ 14.2 → multichm   (tall, sparse → flattest arm wins)
└─ pdens ≥ 4.5
   ├─ height_cv < 0.46 → multichm
   └─ height_cv ≥ 0.46 → segmentanytree (rumple ≥ 9.2 → multichm)
```

## Readings

- **The density crossover is real and the oracle exploits it.** SegmentAnyTree is
  the best arm at native/8 pts/m² (F1 0.443/0.425) but collapses by 1 pt/m²
  (0.102), while multichm is the flattest arm and the strongest low-density
  baseline. A per-cell oracle that switches on this gains **+0.062 F1** and, more
  importantly for the meta-pipeline, **+0.027 understory recall** (0.404 vs
  multichm's 0.377) over the best fixed arm.
- **But a simple learned router barely captures it (+0.002 F1) and loses
  understory.** The CART recovers the right *structure* — it splits on `pdens`
  near the 4–8 pts/m² crossover and routes dense cells to SegmentAnyTree — yet
  its pooled F1 is statistically indistinguishable from just shipping
  multichm, and its understory recall (0.345) is **below** multichm's (0.377). The
  per-cell F1 gaps near the crossover are small and noisy, multichm is a strong
  flat default, and 229 cells across four classes is thin for leave-one-plot-out
  generalization.
- **Strategic conclusion for the meta-pipeline: prefer fusion over selection.**
  Routing's realizable gain here is marginal, whereas the fusion union (#P1) lifts
  understory recall +0.15 over the best single arm on the same cells. Detector
  *selection* leaves most of its (already modest) oracle headroom on the table;
  detector *fusion* is the more promising direction. The router is still useful
  as a density-gated default (ship multichm, escalate to SegmentAnyTree where
  `pdens ≥ ~5` and structure is rough), but it is not the meta-pipeline's win.
- **multichm is the single best default**, consistent with the model benchmark:
  its flat density response makes it the safest fixed choice when only one arm can
  ship.

## Caveats

- **Labels are pooled-F1 argmax**, so the oracle and router optimize F1; on the
  meta-pipeline's priority metric (understory recall) the oracle still wins but
  the learned router does not. A router trained to maximize understory recall
  directly might fare differently.
- **Equal-set restriction**: only cells where CHM-VWF, multichm, and SegmentAnyTree
  all scored enter the study; ForestFormer3D (native+8) and Li 2012 (native) are
  candidate picks only where present, so their low pick share partly reflects
  limited rung coverage, not only quality.
- **Cheap features only** (`rumple`, cover, height CV, gap, mean height, density):
  a richer feature set or a confidence-aware router (consuming the #P4 calibrated
  scores) could narrow the oracle gap; this study measures what deploy-time
  geometry alone buys.
- **F1 is coverage-limited** (the #V4 isolated-FP finding), so absolute F1 is a
  lower bound; the policy *ordering* is the result.
