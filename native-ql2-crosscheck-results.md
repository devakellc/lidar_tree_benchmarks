# Native USGS 3DEP Cross-Check of the Density-Ladder Sweep — Results

*Cross-check for [issue #4](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/4).
The density-ladder sweep ([density-ladder-sweep-results.md](density-ladder-sweep-results.md))
**thins** a dense NEON cloud to simulate sparse acquisitions
("decimation-as-simulation", a stated caveat). This note pulls the **native**
USGS 3DEP cloud over the same NEON plots through the public entwine EPT catalog,
runs the **same** CHM-VWF detection pipeline natively, and asks: does decimation
faithfully predict native-sparse behaviour? Tooling:
[`scripts/native_ql2_crosscheck.R`](scripts/native_ql2_crosscheck.R) +
[`scripts/ept_discovery.R`](scripts/ept_discovery.R). Last run: 2026-06-05.*

---

## TL;DR (honest acquisition status first)

- **No public *native ~2-pulse QL2* cloud covers these NEON plots.** Every USGS
  3DEP project in the entwine public catalog
  (`raw.githubusercontent.com/hobuinc/usgs-lidar/.../resources.geojson`) that
  *covers* SOAP / SJER / TEAK is **far denser than QL2 locally**: measured
  native first-return (pulse) density is **~44 / ~5.5 / ~31 pulses/m²** (median
  over plots), not ~2. The intended "pull a native ~2-pulse cloud" cannot be
  done from public data here — reported truthfully rather than faked.
- **The QL designation in a project name is a floor, not the local density.**
  SJER is covered by `CA_FEMAR9Fresno_2_2019` — the `_2_` is the USGS **QL2**
  tag — yet its measured first-return density over the SJER plots is **2.5–8.5**
  (median 5.5), i.e. mostly above the 2 ppsm QL2 nominal.
- **Because no native sparse cloud is available, the equivalence test is a
  cross-source decimation check, made unit- and resolution-consistent.** Both
  arms thin to **all-return density `pdens = 2`** (the sweep's `homogenize`
  unit — *not* a pulse target), and both detect at **`chm_res = 0.5`**:
  - `native_dec2` — the *native non-NEON* survey decimated to `pdens = 2`.
  - `neon_dec2` — the cached *NEON* `rung==2, chm_res==0.5, vwf_a==0.10` arm.
  Same pipeline, same all-return target, same CHM grid, **different source
  sensor**. With both at `res = 0.5` the residual delta is genuinely
  sensor/return-structure, not CHM resolution.
- **`pdens = 2` all-return ⇒ `frdens ≈ 1.4–1.6` first-return.** Thinning *all*
  returns to 2 pts/m² leaves the *first-return* (pulse) density below 2 — so
  neither arm is "2 pulses/m² QL2". Measured first-return density per arm:

  | Site | native_dec2 frdens (med) | neon_dec2 frdens (med) |
  |------|:------------------------:|:----------------------:|
  | SOAP | 1.46 | 1.46 |
  | SJER | 1.69 | 1.60 |
  | TEAK | 1.35 | 1.43 |

  Both arms land at the same sub-2 pulse density, so the two-arm comparison is
  unit-consistent. We no longer call either arm "2 pulses/m² QL2".
- **Verdict — with CHM resolution matched, decimation is faithful overall and
  the residual error is a modest overstory bias.** Pooled overall recall now
  agrees within **0.05 at every site** (SOAP Δ = −0.05, SJER Δ = −0.01, TEAK
  Δ = −0.003). The largest remaining gaps are in the **dominant** class (SOAP
  −0.10, TEAK −0.10), where the native survey decimated to `pdens = 2` still
  under-detects the *NEON* cloud decimated to the same target. This is much
  smaller than the pre-fix figure: most of the previously reported −0.06…−0.09
  overall and −0.21…−0.23 dominant gaps were a **CHM-resolution confound**
  (`native_dec2` had been detecting at `res = 1.0` while `neon_dec2` pooled at
  `res = 0.5`), not a sensor difference.
- **Native-full context:** at full native density these surveys detect at
  overall recall 0.49 (SOAP), 0.41 (SJER), 0.34 (TEAK) — consistent with the
  sweep's native rung and the SJER→SOAP→TEAK structure gradient.

---

## 1. EPT discovery (real, reproducible)

[`scripts/ept_discovery.R`](scripts/ept_discovery.R) downloads the public 3DEP
entwine boundary index (`usgs.entwine.io/data/boundaries.json` is **404** as of
2026-06; the live mirror is
`raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson`,
~8.7 MB, WGS84), reprojects each NEON plot centroid (UTM 11N / EPSG:32611) into
WGS84, and point-in-polygons against every project footprint. Covering projects
found (1 km buffer around all plot centroids per site):

| Site | Covering public EPT project(s) | Resolved (preferred) EPT | EPT SRS |
|------|--------------------------------|--------------------------|---------|
| SOAP | `CA_SierraNevada_14_B22` | `.../CA_SierraNevada_14_B22/ept.json` | EPSG:3857 |
| SJER | `CA_FEMAR9Fresno_2_2019` (QL2), `CA_SierraNevada_11_B22` | `.../CA_FEMAR9Fresno_2_2019/ept.json` | EPSG:3857 |
| TEAK | `CA_SierraNevada_14_B22` | `.../CA_SierraNevada_14_B22/ept.json` | EPSG:3857 |

SOAP and TEAK are covered by **only one** public project — the 2022
`CA_SierraNevada_14_B22` high-density survey. There is **no QL2 project** over
them in the catalog. SJER has two: discovery prefers the QL2-tagged
`CA_FEMAR9Fresno_2_2019` (name score bonus for the `_2_` QL2 tag), the
genuine candidate.

All EPTs are **EPSG:3857** (Web Mercator). Distances there are inflated ~1.32×
at this latitude, so the pipeline reprojects to **UTM 11N (EPSG:32611)** before
any metric step — the AOI box corners are transformed UTM→3857 to form the
`readers.ept` `bounds`, and `filters.reprojection out_srs=EPSG:32611` (mirroring
[`scripts/extract.json`](scripts/extract.json)) brings the points back to UTM
before detection.

### Provenance-gated LAZ cache

Each per-plot LAZ is pulled once and cached under `neon/<SITE>/ql2/<plot>.laz`
with a sidecar manifest `<plot>.laz.json` recording **`ept_url` + `pad` +
`outcrs`**. The cache is reused only when all three match; a different EPT, AOI
pad, or reprojection target re-pulls. A pure detection-parameter change (e.g.
this run's CHM-resolution match) does **not** trigger a slow EPT re-pull. This
run reused all cached LAZ — **no EPT re-pull occurred** (manifests for the
pre-existing LAZ were backfilled after verifying the EPT/pad/outcrs were
unchanged).

## 2. Acquisition status — native density actually obtained

For each eligible plot (≥6 mapped live-tree stems in the plot core; core ±20 m
tower / ±10 m distributed, mirroring `run_sweep.R`) a ±45 m AOI was pulled via
`pdal pipeline` (readers.ept → drop class 7/18 + withheld → reproject → laz) and
its first-return density measured directly from the laz.

| Site | EPT project | Plots pulled | Native first-return ppsm (min / median / max) | ~2 ppsm? |
|------|-------------|:-----------:|:---------------------------------------------:|:--------:|
| SOAP | CA_SierraNevada_14_B22 | 16 | 32.9 / **43.9** / 66.1 | no (≈22×) |
| SJER | CA_FEMAR9Fresno_2_2019 | 8 | 2.5 / **5.5** / 8.5 | partial (3 plots ≈2.5) |
| TEAK | CA_SierraNevada_14_B22 | 19 | 21.4 / **31.4** / 43.8 | no (≈16×) |

The **only** plots that land near 2 ppsm are three SJER distributed plots —
`SJER_046` (2.5), `SJER_054` (2.7), `SJER_008` (2.8). Everything else is a
high-density survey. All three sites are therefore classified `dense` (median
native first-return density above the QL2 band). This is the honest reason the
"native QL2" comparison is reframed below as a **decimation cross-source check**
rather than a true native sparse run.

## 3. Native-vs-decimated comparison (pooled, by crown class)

Pooling matches the sweep exactly (`analyze_sweep.R`):

- **Recall** = `sum(TP) / sum(n_ref)`; per-class TP recovered as
  `round(rec_class × n_class)`.
- **Precision** = `sum(tp_core) / sum(n_det)`, where `tp_core =
  round(precision × n_det)` is recovered **per row**. The per-row precision
  numerator is core-only, but `TP` includes boundary matches outside the core;
  pooling `sum(TP)/sum(n_det)` would mix a boundary-inclusive numerator with a
  core-only denominator and **inflate precision**. This was fixed — corrected
  precision is reported below.

Three engines per site, **all detecting at `chm_res = 0.5`** for the two
decimated arms:

- `native_full` — detection at the survey's **native** density (context); its
  CHM resolution is derived from the native first-return density (`res = 0.25`
  for the dense surveys).
- `native_dec2` — the native non-NEON cloud decimated to **all-return
  `pdens = 2`**, CHM **pinned to `res = 0.5`** to match `neon_dec2`.
- `neon_dec2` — the **cached** NEON `rung==2, chm_res==0.5, vwf_a==0.10` arm
  from `sweep_results.csv`, same plots.

The equivalence test compares **`native_dec2` vs `neon_dec2`** (matched
all-return target, matched CHM grid, different source sensor).
`decimate_points(homogenize)` is random, so `*_dec2` recall jitters ±1–2 TP
run-to-run; `native_full` is deterministic.

### SOAP — mixed conifer (native median 44 ppsm; 16 plots, 222 stems)

| engine | overall | dominant | codominant | intermediate | suppressed | precision |
|--------|:------:|:--------:|:----------:|:------------:|:----------:|:---------:|
| native_full | **0.486** | 0.507 | 0.552 | 0.263 | 0.429 | 0.281 |
| native_dec2 (pdens→2, res=0.5) | **0.288** | 0.348 | 0.314 | 0.158 | 0.000 | 0.410 |
| neon_dec2 (pdens→2, res=0.5) | **0.333** | 0.449 | 0.324 | 0.158 | 0.000 | 0.466 |
| **Δ (native_dec2 − neon_dec2)** | **−0.045** | **−0.101** | −0.010 | 0.000 | 0.000 | −0.056 |

### SJER — open oak savanna (native median 5.5 ppsm; 8 plots, 71 stems)

| engine | overall | dominant | codominant | intermediate | suppressed | precision |
|--------|:------:|:--------:|:----------:|:------------:|:----------:|:---------:|
| native_full | **0.408** | 0.696 | 0.500 | 0.500 | — | 0.292 |
| native_dec2 (pdens→2, res=0.5) | **0.394** | 0.696 | 0.417 | 0.500 | — | 0.300 |
| neon_dec2 (pdens→2, res=0.5) | **0.408** | 0.652 | 0.500 | 0.500 | — | 0.337 |
| **Δ (native_dec2 − neon_dec2)** | **−0.014** | +0.043 | −0.083 | 0.000 | — | −0.037 |

### TEAK — red fir / high structure (native median 31 ppsm; 19 plots, 392 stems)

| engine | overall | dominant | codominant | intermediate | suppressed | precision |
|--------|:------:|:--------:|:----------:|:------------:|:----------:|:---------:|
| native_full | **0.342** | 0.679 | 0.271 | 0.107 | 0.500 | 0.334 |
| native_dec2 (pdens→2, res=0.5) | **0.222** | 0.452 | 0.171 | 0.054 | 0.500 | 0.400 |
| neon_dec2 (pdens→2, res=0.5) | **0.224** | 0.548 | 0.142 | 0.054 | 0.500 | 0.414 |
| **Δ (native_dec2 − neon_dec2)** | **−0.003** | **−0.095** | +0.029 | 0.000 | 0.000 | −0.014 |

### By height band (recall; short <8 m, mid 8–15 m, tall ≥15 m)

| Site | engine | tall | mid | short |
|------|--------|:----:|:---:|:-----:|
| SOAP | native_dec2 | 0.397 | 0.230 | 0.270 |
| SOAP | neon_dec2 | 0.534 | 0.241 | 0.257 |
| SJER | native_dec2 | 1.000 | 0.800 | 0.444 |
| SJER | neon_dec2 | 1.000 | 0.800 | 0.444 |
| TEAK | native_dec2 | 0.424 | 0.184 | 0.071 |
| TEAK | neon_dec2 | 0.456 | 0.136 | 0.084 |

## 4. Equivalence verdict

**Does decimation faithfully predict native-sparse behaviour?**

*Caveat up front:* this could not be tested against a true native ~2-pulse cloud,
because none is publicly available over these sites (§1–2). It was tested
against the strongest available proxy — a *different* native sensor decimated to
the same all-return `pdens = 2`, at the same CHM resolution.

- **Overall: faithful within 0.05 recall at every site.** With the CHM grid
  matched, pooled overall recall agrees to **−0.045 (SOAP), −0.014 (SJER),
  −0.003 (TEAK)**. Decimation reproduces both the qualitative density response
  (overstory holds up, understory/short trees collapse) and the absolute pooled
  recall to within a few points.
- **The residual error is a modest, structured overstory bias.** The remaining
  gaps concentrate in the **dominant** class: SOAP −0.10, TEAK −0.10 (SJER
  +0.04, where the native cloud is already near the target). The likely cause is
  **sensor/return-structure differences** — all-return homogenisation equalises
  point count but not footprint, scan pattern, or multi-return distribution, so
  the resulting CHM apexes differ slightly for the tallest trees. The height-band
  view agrees: the tall band carries the SOAP/TEAK gap (−0.14 / −0.03), while
  mid/short are within ±0.05.
- **Most of the previously reported gap was a resolution confound, not the
  sensor.** Before this fix `native_dec2` selected `res = 1.0` from its decimated
  first-return density (`frdens < 4`) while `neon_dec2` was pooled at `res = 0.5`
  — the two arms ran at *different* CHM resolutions. Matching them collapses the
  overall delta from −0.06…−0.09 to ≤0.05 and the dominant-class delta from
  −0.21…−0.23 to ~−0.10. The honest reading is that the coarser CHM, not the
  source sensor, drove most of the old under-detection.
- **Precision is also corrected.** Recovering core-only `tp_core` per row lowers
  the pooled precision relative to the previous (inflated) `sum(TP)/sum(n_det)`:
  e.g. TEAK `native_full` 0.40→0.33, SOAP `native_full` 0.31→0.28. The
  native_dec2 vs neon_dec2 precision gap is small at every site (≤0.06).

**Bottom line for the sweep's caveat:** with resolution and density units matched,
the density-ladder's decimation-as-simulation is a **faithful predictor of pooled
sparse-cloud recall (within ~0.05) and of the qualitative density response**, with
a **residual ~0.10 overstory (dominant/tall) optimism** when a much denser cloud
is thinned rather than acquired sparse. Treat the sweep's absolute sparse-rung
recall as **slightly optimistic for the overstory**, but otherwise a reasonable
stand-in. A true public native ~2-pulse validation over these NEON sites is
**not currently possible** from the entwine catalog. The earlier, larger gaps in
this note were an artefact of an unmatched CHM resolution between the two arms,
now corrected.

## 5. Reproduce

```sh
export CLAUDE_JOB_DIR=/path/to/work        # NEON ground truth + sweep_results live here

# 1. discover covering EPT projects (writes neon/<SITE>/ql2/ept_candidates.csv)
Rscript scripts/ept_discovery.R SITES=SOAP,SJER,TEAK

# 2. pull native AOIs via PDAL, detect, score, compare to decimated-2 rung
Rscript scripts/native_ql2_crosscheck.R SITES=SOAP,SJER,TEAK
#   per-site:   neon/<SITE>/ql2/ql2_detect_results.csv  (+ <plot>.laz + .laz.json)
#   combined:   neon/native_ql2_vs_decimated.csv
# Override a site's EPT:  EPT_SJER=https://.../ept.json
# Tune the QL2 band classifier: QL2_LO=1.5 QL2_HI=4
```

All `*.laz` / `*.csv` outputs are gitignored — regenerate them. The cache is
provenance-gated (`<plot>.laz.json` records ept_url/pad/outcrs); a detection-only
re-run reuses the cached LAZ. Numbers above are from the 2026-06-05 run; `*_dec2`
rows carry ±1–2 TP of decimation randomness.
