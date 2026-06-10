# NEON benchmark sites — field data & LiDAR acquisitions

Reference for the three **Domain D17 (Pacific Southwest)** sites used in the
density-ladder sweep and follow-on analyses. All sites share **UTM zone 11N /
EPSG:32611** for NEON woody-vegetation and AOP LiDAR products.

**Field ground truth:** NEON Woody Plant Vegetation Structure `DP1.10098.001`
(mapped stems >10 cm DBH in 20×20 m distributed and 40×40 m tower plots).
**Airborne LiDAR (primary):** NEON discrete-return point cloud `DP1.30003.001`
(AOP tiles via `neon_download_lidar.R`, restricted to 1 km tiles overlapping
field plots). **Cross-check LiDAR:** USGS 3DEP public EPT projects
([`ept_discovery.R`](../scripts/ept_discovery.R), issue #4).

Acquisition months below are from the NEON Data Portal API
(`GET /api/v0/sites/{SITE}`, product `DP1.30003.001`) as of 2026-06-10.

---

## Site overview

| Code | Site name | Forest type (benchmark role) | Lat, lon | Benchmark plots / live stems† |
|------|-----------|------------------------------|----------|-------------------------------|
| **SJER** | San Joaquin Experimental Range | Open oak / foothill-pine woodland (open-canopy gradient) | 37.11°N, 119.73°W | 8 / 71 |
| **SOAP** | Soaproot Saddle | Mixed conifer/deciduous (**anchor** site) | 37.03°N, 119.26°W | 18 / 232 |
| **TEAK** | Lower Teakettle | Red fir / subalpine conifer (dense-canopy gradient) | 37.01°N, 119.01°W | 20 / 396 |

†Live, mapped tree stems in the sweep ground truth (`ground_truth_stems.csv`),
paired with the **2021** NEON AOP acquisition (±4 yr nearest field measurement).

---

## NEON AOP LiDAR (`DP1.30003.001`)

| Site | All portal acquisition months | Benchmark acquisition | Sensor (benchmark year) | Native all-return pts/m² | Native first-return pts/m² |
|------|------------------------------|----------------------|-------------------------|:------------------------:|:--------------------------:|
| SJER | 2013-06, 2017-03, 2018-03, 2019-03, **2021-03**, 2023-04, 2024-04 | **2021-03** | Optech Galaxy Prime (~20 pts/m² class) | 16.3 | 9.0 |
| SOAP | 2013-06, 2017-07, 2018-06, 2019-06, **2021-07**, 2023-06, 2023-07, 2024-06, 2026-04 | **2021-07** | Optech Galaxy Prime (~20 pts/m² class) | 18.2 | 11.9 |
| TEAK | 2013-06, 2017-06, 2018-06, 2019-06, **2021-07**, 2023-07, 2024-06 | **2021-07** | Optech Galaxy Prime (~20 pts/m² class) | 19.2 | 11.9 |

**Sensor timeline (NEON airborne).** Optech Gemini era (2013–2020) yields
~4–6 pts/m² at these sites; **2021+ Galaxy Prime** is the first acquisition
that clears the repository's >8 pts/m² design threshold. Pre-2021 site-years
remain on the portal but are not used in the benchmark pipeline.

**Download:** `Rscript scripts/neon_download_lidar.R SITE=<CODE> YEAR=2021`
(after `neon_ground_truth.R`).

---

## USGS 3DEP LiDAR (native QL2 cross-check, issue #4)

Public entwine EPT projects covering each site (preferred project per
[`native-ql2-crosscheck-results.md`](../results/native-ql2-crosscheck-results.md)):

| Site | Covering EPT project(s) | Acquisition (from project name) | Median native first-return pts/m² over plots | Notes |
|------|-------------------------|---------------------------------|:--------------------------------------------:|-------|
| SOAP | `CA_SierraNevada_14_B22` | **2022** (Sierra Nevada block B22) | ~44 | Only one public project; no QL2-tagged alternative |
| SJER | `CA_FEMAR9Fresno_2_2019`, `CA_SierraNevada_11_B22` | **2019** (preferred QL2), 2022 (alt.) | ~5.5 | Discovery prefers the QL2-tagged 2019 project |
| TEAK | `CA_SierraNevada_14_B22` | **2022** (Sierra Nevada block B22) | ~31 | Same high-density block as SOAP |

EPT tiles are stored in **EPSG:3857** (Web Mercator); metric detection reprojects
to UTM 11N before CHM construction ([`extract.json`](../scripts/extract.json)
/ `native_ql2_crosscheck.R`).

---

## Field measurement timing vs. 2021 LiDAR

Ground truth pairs each stem with the `apparentindividual` record **nearest
2021 within ±4 yr** (`neon_ground_truth.R`: `meas_year`, `dist21`). Exact-2021
field coverage (issue #5):

| Site | Mapped live-tree stems | Measured in 2021 | Exact-2021 share |
|------|:----------------------:|:----------------:|:----------------:|
| TEAK | 483 | 233 | 48% |
| SOAP | 268 | 52 | 19% |
| SJER | 113 | 0 | 0% |

SJER field stems were mostly measured in **2022 and 2024**, not 2021; treat its
sweep metrics as carrying the full ±4 yr temporal slack.
