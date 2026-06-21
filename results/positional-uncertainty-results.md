# Monte-Carlo stem-position uncertainty (#V3)

Every leaderboard delta the router (#P2) would switch arms on is a single point
estimate against stem coordinates that carry real uncertainty — `neon_ground_
truth.R` stores `pos_unc` (NEON coordinate uncertainty + 0.3 m TruPulse
rangefinder) per stem, yet `score_plot` ignores it. This arm puts **confidence
bands** on recall / precision / F1 (and per-crown-class recall) by re-scoring every
ladder arm — and the #P1 fusion union/layered arm — under K reproducible draws of
the field positions, answering two questions: which arm-vs-arm gaps survive
ground-truth jitter, and whether fusion is *more stable* than any single arm.

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
# CORES=1 is required: prep_cell -> detect_lasr uses lasR exec, which deadlocks under fork
Rscript scripts/mc_positional_uncertainty.R SITES=SOAP,SJER,TEAK K=200
# -> work/neon/<SITE>/positional_uncertainty.csv (median + p05/p95 per arm x tol)
```

## What this is

Detections — including the fused union/layered sets — are materialized **once**
per native cell (the expensive step). Each of K=200 draws then perturbs only the
stems: `perturb_positions` (in [`sweep_lib.R`](../scripts/sweep_lib.R),
unit-tested) draws each stem's (E,N) from an independent 2-D Gaussian with
per-stem σ = its `pos_unc`, re-runs `greedy_match`/`score_plot` against the fixed
detections, and pools by the canonical summed-count rule. Draw *k* uses
`seed = seed_for(site, plot, "native") + k`, so the bands reproduce exactly. The
report gives the median and 5th/95th percentiles over the 200 draws.

## Generated tables

Native density, three sites, K=200 draws, matching tol = 4 m (with a 3/4/5 m
co-sweep on SOAP).

### F1 bands — median [p05, p95]

| arm | SOAP | SJER | TEAK |
|---|---|---|---|
| chm_vwf | 0.374 [0.362, 0.388] | 0.319 [0.301, 0.338] | 0.361 [0.352, 0.372] |
| multichm | 0.438 [0.425, 0.450] | 0.354 [0.347, 0.365] | 0.413 [0.403, 0.424] |
| li2012 | 0.356 [0.345, 0.369] | 0.270 [0.260, 0.280] | 0.368 [0.360, 0.378] |
| segmentanytree | 0.459 [0.446, 0.472] | 0.299 [0.287, 0.312] | 0.476 [0.466, 0.484] |
| forestformer3d | 0.267 [0.259, 0.278] | 0.241 [0.230, 0.252] | 0.303 [0.293, 0.315] |
| fusion_union | 0.323 [0.316, 0.331] | 0.225 [0.218, 0.237] | 0.413 [0.405, 0.420] |
| fusion_layered | 0.331 [0.323, 0.338] | 0.225 [0.218, 0.234] | 0.425 [0.416, 0.433] |

### F1 band width (p95 − p05) — smaller = more stable

| arm | SOAP | SJER | TEAK | mean |
|---|--:|--:|--:|--:|
| chm_vwf | 0.026 | 0.037 | 0.020 | 0.028 |
| multichm | 0.025 | 0.018 | 0.020 | 0.021 |
| li2012 | 0.023 | 0.021 | 0.018 | 0.020 |
| segmentanytree | 0.026 | 0.025 | 0.018 | 0.023 |
| forestformer3d | 0.020 | 0.022 | 0.023 | 0.022 |
| **fusion_union** | **0.014** | 0.019 | 0.016 | **0.016** |
| **fusion_layered** | 0.016 | 0.016 | 0.017 | **0.017** |

### Gap survival: SegmentAnyTree vs multichm F1 (do the bands separate?)

| site | SegmentAnyTree | multichm | gap survives jitter? |
|---|---|---|---|
| SOAP | 0.459 [0.446, 0.472] | 0.438 [0.425, 0.450] | **NO — bands overlap** |
| SJER | 0.299 [0.287, 0.312] | 0.354 [0.347, 0.365] | YES (multichm > SAT) |
| TEAK | 0.476 [0.466, 0.484] | 0.413 [0.403, 0.424] | YES (SAT > multichm) |

### SOAP matching-tol co-sweep (tol and pos_unc trade off)

| tol (m) | fusion_union | multichm | segmentanytree |
|--:|---|---|---|
| 3 | 0.285 [0.277, 0.293] | 0.369 [0.357, 0.382] | 0.397 [0.384, 0.414] |
| 4 | 0.323 [0.316, 0.331] | 0.438 [0.425, 0.450] | 0.459 [0.446, 0.472] |
| 5 | 0.349 [0.342, 0.356] | 0.477 [0.467, 0.490] | 0.508 [0.497, 0.518] |

## Readings

- **Fusion is the most stable arm under field jitter.** The #P1 union/layered
  consensus has the **tightest F1 bands of any arm** (mean width 0.016–0.017 vs
  0.020–0.028 for the single arms), at every site. Consensus over multiple
  detectors averages out the per-stem matching noise that moves any single arm's
  score, so the fused metric is the most reproducible — exactly the property that
  makes a fusion recommendation defensible against ground-truth uncertainty.
- **The SegmentAnyTree-vs-multichm gap is site-specific and only sometimes
  real.** On SOAP the two arms' F1 bands **overlap** (0.459 [0.446, 0.472] vs
  0.438 [0.425, 0.450]) — the apparent "SAT beats multichm at native" edge does
  *not* survive stem jitter, so a router switching SOAP-native to SAT would be
  acting on noise. On SJER multichm cleanly wins and on TEAK SAT cleanly wins
  (bands fully separated). The honest routing rule: switch only where the bands
  separate (SJER, TEAK), default otherwise (SOAP).
- **The matching tolerance and `pos_unc` trade off, but the ordering holds.**
  Widening tol 3→5 m lifts every arm's F1 (more jittered stems still find their
  apex) without materially widening the bands or reordering the arms — so the
  stability and gap-survival conclusions are not artifacts of the 4 m default.
- **Fusion's value is recall and stability, not F1.** Consistent with #P1, the
  fused F1 sits below the best single arm (coverage-limited precision), but its
  bands are the tightest and its understory recall the highest (#P1). The bands
  make the trade explicit: fusion buys reproducible recall, not a higher F1 point
  estimate.

## Caveats

- **Independent per-stem Gaussian jitter** at σ = `pos_unc`. Real NEON error has
  a shared per-plot georeferencing component (correlated across stems) the i.i.d.
  draw omits; correlated error would widen the bands somewhat, so these are a
  lower bound on uncertainty. The relative arm ordering (the result) is robust to
  this.
- **Only stem positions are perturbed** — detections are fixed. This isolates
  ground-truth uncertainty from detector stochasticity (the question the issue
  asks); detector-side variance is out of scope.
- **F1/precision are coverage-limited** (the #V4 isolated-FP finding), so absolute
  bands are lower bounds; the band *widths* and *overlaps* are the result.
- **Native density only**; the driver accepts `RUNG=`/`K=`/`TOLS=` to extend.
