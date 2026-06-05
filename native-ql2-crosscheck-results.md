# Native USGS 3DEP Cross-Check of the Density-Ladder Sweep — Results

*Cross-check for [issue #4](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/4).
The density-ladder sweep ([density-ladder-sweep-results.md](density-ladder-sweep-results.md))
**thins** a dense NEON cloud to simulate sparse acquisitions
("decimation-as-simulation", a stated caveat). This note pulls the **native**
USGS 3DEP cloud over the same NEON plots through the public entwine EPT catalog,
runs the **same** CHM-VWF detection pipeline natively, and asks: does decimation
faithfully predict native-sparse (QL2, ~2 pulses/m²) behaviour? Tooling:
[`scripts/native_ql2_crosscheck.R`](scripts/native_ql2_crosscheck.R) +
[`scripts/ept_discovery.R`](scripts/ept_discovery.R). Last run: 2026-06-05.*

---

## TL;DR (honest acquisition status first)

- **No public *native QL2* cloud exists over these NEON plots.** Every USGS 3DEP
  project in the entwine public catalog
  (`raw.githubusercontent.com/hobuinc/usgs-lidar/.../resources.geojson`) that
  *covers* SOAP / SJER / TEAK is **far denser than QL2 locally**: measured
  native first-return (pulse) density is **~44 / ~5.5 / ~31 pulses/m²** (median
  over plots), not ~2. The intended "pull native ~2 ppsm" cannot be done from
  public data here — reported truthfully rather than faked.
- **The QL designation in a project name is a floor, not the local density.**
  SJER is covered by `CA_FEMAR9Fresno_2_2019` — the `_2_` is the USGS **QL2**
  tag — yet its measured first-return density over the SJER plots is **2.5–8.5**
  (median 5.5), i.e. mostly above the 2 ppsm QL2 nominal.
- **Because true native-QL2 is unavailable, the equivalence test becomes a
  cross-source decimation check:** decimate the *native non-NEON* survey to
  2 ppsm (`native_dec2`) and compare it, per crown class, to the cached
  NEON-decimated-2 rung (`neon_dec2`, `rung==2, chm_res=0.5, vwf_a=0.10`). Same
  pipeline, same target density, **different source sensor**.
- **Verdict — decimation is *directionally* faithful but *not* numerically
  equivalent, and the error is structured.** Pooled overall recall agrees at
  SJER (Δ ≈ −0.04) but the *native cloud decimated to 2* under-detects relative
  to *NEON decimated to 2* at the two dense sites: SOAP Δ = −0.09, TEAK
  Δ = −0.06 overall, with the **largest gaps in the overstory** (dominant trees:
  SOAP −0.23, TEAK −0.25; tall band similar). Decimation predicts the *shape* of
  the density response but can be off by 6–25 recall points depending on source
  sensor and crown class.
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

All EPTs are **EPSG:3857** (Web Mercator). Distances there are inflated ~1.32× at
this latitude, so the pipeline reprojects to **UTM 11N (EPSG:32611)** before any
metric step — the AOI box corners are transformed UTM→3857 to form the
`readers.ept` `bounds`, and `filters.reprojection out_srs=EPSG:32611` (mirroring
[`scripts/extract.json`](scripts/extract.json)) brings the points back to UTM
before detection.

## 2. Acquisition status — native density actually obtained

For each eligible plot (≥6 mapped live-tree stems in the plot core; core ±20 m
tower / ±10 m distributed, mirroring `run_sweep.R`) a ±45 m AOI was pulled via
`pdal pipeline` (readers.ept → drop class 7/18 + withheld → reproject → laz) and
its first-return density measured directly from the laz.

| Site | EPT project | Plots pulled | Native first-return ppsm (min / median / max) | QL2 (~2)? |
|------|-------------|:-----------:|:---------------------------------------------:|:---------:|
| SOAP | CA_SierraNevada_14_B22 | 16 | 33.0 / **43.9** / 67.5 | no (≈22× QL2) |
| SJER | CA_FEMAR9Fresno_2_2019 | 8 | 2.5 / **5.5** / 8.5 | partial (3 plots ≈2.5) |
| TEAK | CA_SierraNevada_14_B22 | 19 | 21.4 / **31.4** / 43.9 | no (≈16× QL2) |

The **only** plots that land near the QL2 target are three SJER distributed
plots — `SJER_046` (2.5), `SJER_054` (2.7), `SJER_008` (2.8). Everything else is
a high-density survey. This is the honest reason the "native QL2" comparison is
reframed below as a **decimation cross-source check** rather than a true native
QL2 run.

## 3. Native-vs-decimated comparison (pooled, by crown class)

Pooling matches the sweep exactly: **sum counts, never average rates** — site
recall = `sum(TP) / sum(n_ref)`; per-class TP recovered as
`round(rec_class × n_class)`. Three engines per site:

- `native_full` — detection at the survey's **native** density (context).
- `native_dec2` — the native non-NEON cloud **decimated to 2 ppsm** (matched to
  the rung; the cross-source stand-in for native QL2).
- `neon_dec2` — the **cached** NEON-decimated-2 rung from `sweep_results.csv`
  (`rung==2, chm_res=0.5, vwf_a=0.10`), same plots.

The equivalence test compares **`native_dec2` vs `neon_dec2`** (matched target
density, different source sensor). `decimate_points(homogenize)` is random, so
`*_dec2` recall jitters ±1–2 TP run-to-run; `native_full` is deterministic.

### SOAP — mixed conifer (native median 44 ppsm; 16 plots, 222 stems)

| engine | overall | dominant | codominant | intermediate | suppressed |
|--------|:------:|:--------:|:----------:|:------------:|:----------:|
| native_full | **0.491** | 0.507 | 0.562 | 0.263 | 0.429 |
| native_dec2 (→2) | **0.248** | 0.217 | 0.286 | 0.184 | 0.143 |
| neon_dec2 (→2) | **0.333** | 0.449 | 0.324 | 0.158 | 0.000 |
| **Δ (native_dec2 − neon_dec2)** | **−0.085** | **−0.232** | −0.038 | +0.026 | +0.143 |

### SJER — open oak savanna (native median 5.5 ppsm; 8 plots, 71 stems)

| engine | overall | dominant | codominant | intermediate | suppressed |
|--------|:------:|:--------:|:----------:|:------------:|:----------:|
| native_full | **0.408** | 0.696 | 0.500 | 0.500 | — |
| native_dec2 (→2) | **0.366** | 0.652 | 0.417 | 0.500 | — |
| neon_dec2 (→2) | **0.408** | 0.652 | 0.500 | 0.500 | — |
| **Δ (native_dec2 − neon_dec2)** | **−0.042** | 0.000 | −0.083 | 0.000 | — |

### TEAK — red fir / high structure (native median 31 ppsm; 19 plots, 392 stems)

| engine | overall | dominant | codominant | intermediate | suppressed |
|--------|:------:|:--------:|:----------:|:------------:|:----------:|
| native_full | **0.342** | 0.679 | 0.271 | 0.107 | 0.500 |
| native_dec2 (→2) | **0.163** | 0.298 | 0.133 | 0.054 | 0.500 |
| neon_dec2 (→2) | **0.224** | 0.548 | 0.142 | 0.054 | 0.500 |
| **Δ (native_dec2 − neon_dec2)** | **−0.066** | **−0.214** | −0.009 | 0.000 | 0.000 |

### By height band (recall; short <8 m, mid 8–15 m, tall ≥15 m)

| Site | engine | tall | mid | short |
|------|--------|:----:|:---:|:-----:|
| SOAP | native_dec2 | 0.362 | 0.161 | 0.230 |
| SOAP | neon_dec2 | 0.534 | 0.241 | 0.257 |
| SJER | native_dec2 | 1.000 | 0.500 | 0.556 |
| SJER | neon_dec2 | 1.000 | 0.800 | 0.444 |
| TEAK | native_dec2 | 0.304 | 0.107 | 0.058 |
| TEAK | neon_dec2 | 0.456 | 0.136 | 0.084 |

## 4. Equivalence verdict

**Does decimation faithfully predict native-sparse behaviour?**

*Caveat up front:* this could not be tested against a true native ~2 ppsm cloud,
because none is publicly available over these sites (§1–2). It was tested
against the strongest available proxy — a *different* native sensor decimated to
the same 2 ppsm.

- **Directionally yes, numerically no.** Decimation reproduces the *qualitative*
  density response (overstory holds up, understory/short trees collapse, the
  SJER→SOAP→TEAK ordering is preserved) but the absolute recall can differ by
  **6–9 points overall** and **20–25 points in the dominant/tall classes**
  between two clouds decimated to the *same* density.
- **The error is one-directional and structured:** at both dense sites the
  *native survey* decimated to 2 under-detects relative to the *NEON* cloud
  decimated to 2, concentrated in the **overstory/tall** trees (SOAP dominant
  −0.23, TEAK dominant −0.21; tall band −0.17 / −0.15). The likely cause is
  **sensor/return-structure differences** — pulse-density homogenisation
  equalises first-return count but not footprint, scan pattern, or
  multi-return distribution, so the resulting CHM apexes differ.
- **At the genuinely sparse site (SJER, native ~5.5, three plots ≈2.5)
  decimation is most faithful** (overall Δ −0.04; dominant Δ 0.00). Where the
  native cloud is already near the target, the residual sensor gap is small.

**Bottom line for the sweep's caveat:** the density-ladder's decimation-as-
simulation is a **reasonable predictor of the trend and of overstory behaviour
to within ~0.05 recall when the native cloud is itself near the target**, but it
**systematically over-states sparse-cloud recall for dominant/tall trees by up
to ~0.2** when a much denser cloud is thinned. Treat the sweep's absolute
sparse-rung recall as an **optimistic** estimate for the overstory, not a
faithful native-QL2 value. A true public native-QL2 validation over these NEON
sites is **not currently possible** from the entwine catalog.

## 5. Reproduce

```sh
export CLAUDE_JOB_DIR=/path/to/work        # NEON ground truth + sweep_results live here

# 1. discover covering EPT projects (writes neon/<SITE>/ql2/ept_candidates.csv)
Rscript scripts/ept_discovery.R SITES=SOAP,SJER,TEAK

# 2. pull native AOIs via PDAL, detect, score, compare to decimated-2 rung
Rscript scripts/native_ql2_crosscheck.R SITES=SOAP,SJER,TEAK
#   per-site:   neon/<SITE>/ql2/ql2_detect_results.csv  (+ <plot>.laz)
#   combined:   neon/native_ql2_vs_decimated.csv
# Override a site's EPT:  EPT_SJER=https://.../ept.json
# Tune the QL2 band classifier: QL2_LO=1.5 QL2_HI=4
```

All `*.laz` / `*.csv` outputs are gitignored — regenerate them. Numbers above are
from the 2026-06-05 run; `*_dec2` rows carry ±1–2 TP of decimation randomness.
