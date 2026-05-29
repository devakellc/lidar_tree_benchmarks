# lasR EPT acquisition sweep results

Benchmarks of `reader_rectangles(AOI) + summarise()` on remote USGS 3DEP EPT,
run on an EC2 instance (16 usable OpenMP threads, 40 GB RAM) against
`CA_CarrHirzDeltaFires_2_2019`.

**lasR:** `0.21.0` pre-devel (commit `366b956a` and earlier for “old default” rows)  
**Pipeline:** EPT spatial clip + point count (`summarise()`), no CHM/detection  
**Endpoint:** `https://s3-us-west-2.amazonaws.com/usgs-lidar-public/CA_CarrHirzDeltaFires_2_2019/ept.json`

---

## AOIs tested

| Label | Bbox (EPSG:3857) | Approx. side | Points | Notes |
|-------|------------------|--------------|--------|-------|
| Small (original) | `-13578452,4988399,-13577792,4989059` | 660 × 660 m | 5.7M | ~25 ha; matches `extract.json` |
| 5× small | `-13578860,4987991,-13577384,4989467` | 1476 × 1476 m | 31.3M | ~5× area vs small |
| 4× 5× (sweep AOI) | `-13579688,4987253,-13576556,4990205` | 3132 × 2952 m | 112.5M | Used for full parameter sweep |

All coordinates are Web Mercator (EPSG:3857). Distances and density-derived
parameters are inflated ~1.32× vs UTM at this latitude.

---

## Full parameter sweep (112.5M points, 2 reps each)

Script: `scripts/sweep_lasr_ept_params.R`  
CSV: `/tmp/lidar_bench/sweep_lasr_ept_params.csv`

Sorted by median wall time (fastest first):

| Config | Workers | `LASR_EPT_PARTITIONS` | `LASR_EPT_PREFETCH` | `VSI_CACHE_SIZE` | Chunks | Median (s) | Min–Max (s) | Mpts/s | vs seq |
|--------|---------|----------------------|---------------------|------------------|--------|------------|-------------|--------|--------|
| **par16_parts32** | 16 | 32 | default (4) | default (64 MB) | 36 | **47.3** | 47.2–47.4 | 2.38 | **2.93×** |
| par16_parts36 | 16 | 36 | default | default | 36 | 47.7 | 47.3–48.0 | 2.36 | 2.90× |
| par16_parts32_cache1024 | 16 | 32 | default | 1024 MB | 36 | 48.1 | 47.9–48.3 | 2.34 | 2.88× |
| par16_parts32_prefetch32 | 16 | 32 | 32 | default | 36 | 48.4 | 48.4–48.4 | 2.33 | 2.86× |
| par16_parts32_prefetch16 | 16 | 32 | 16 | default | 36 | 48.7 | 48.1–49.4 | 2.31 | 2.84× |
| par16_parts32_cache256 | 16 | 32 | default | 256 MB | 36 | 47.9 | 47.2–48.5 | 2.35 | 2.89× |
| par16_parts32_prefetch8 | 16 | 32 | 8 | default | 36 | 47.8 | 47.7–47.9 | 2.35 | 2.89× |
| par8_cache256 | 8 | default | default | 256 MB | 36 | 51.9 | 51.7–52.1 | 2.17 | 2.66× |
| par8_cache1024 | 8 | default | default | 1024 MB | 36 | 52.0 | 51.7–52.3 | 2.16 | 2.66× |
| par8_prefetch16 | 8 | default | 16 | default | 36 | 52.2 | 52.1–52.3 | 2.15 | 2.65× |
| par8_prefetch8 | 8 | default | 8 | default | 36 | 52.4 | 51.9–52.8 | 2.15 | 2.64× |
| **par8_default** | 8 | default | default | default | 36 | **52.6** | 52.3–52.8 | 2.14 | **2.63×** |
| par8_prefetch32 | 8 | default | 32 | default | 36 | 52.6 | 52.4–52.8 | 2.14 | 2.63× |
| par16_default | 16 | default (`4×workers=64`) | default | default | **100** | 67.3 | 67.1–67.4 | 1.67 | 2.06× |
| par16_parts64 | 16 | 64 | default | default | **100** | 67.3 | 67.1–67.5 | 1.67 | 2.06× |
| sequential | 1 | — | — | — | 1 | 138.4 | 135.7–141.0 | 0.81 | — |

All configs returned **112,519,502 points**.

### Sweep takeaways

1. **Partition count dominates.** Old default `4 × workers` with 16 threads → target 64 → **100 chunks** and ~67 s. Setting `LASR_EPT_PARTITIONS=32` → **36 chunks** and ~47 s (~30% faster than par16 default, ~10% faster than par8).
2. **`LASR_EPT_PREFETCH` (8/16/32)** and **`VSI_CACHE_SIZE` (256 MB / 1 GB)** moved results by ≤1 s on this AOI — not material vs partition tuning.
3. **Best config:** `concurrent_files(16)` + `LASR_EPT_PARTITIONS=32` (or lasR pre-devel default after `366b956a`).

---

## Ad-hoc runs (single rep, varying AOI size)

### Small AOI (~5.7M points)

| Mode | Wall (s) | Chunks | vs seq |
|------|----------|--------|--------|
| sequential | 8.1 | 1 | — |
| concurrent_files(8) | 11.6 | 64 | 0.70× (slower) |
| concurrent_files(16) | 9.4 | 64 | 0.86× (slower) |

Parallel hurt on this tiny AOI: partition overhead > download/compute gain.

### 5× AOI (~31.3M points)

| Mode | Wall (s) | Chunks | vs seq |
|------|----------|--------|--------|
| sequential | 43–47 | 1 | — |
| concurrent_files(8) | 24–25 | 36 | ~1.7× |
| concurrent_files(16), old default | 36–68 | 100 | 0.6–1.2× |

With old default, 16 workers often created 100 chunks and could segfault under load.

### 4× 5× AOI (~112.5M points, single rep)

| Mode | Wall (s) | Chunks | vs seq |
|------|----------|--------|--------|
| sequential | 143.5 | 1 | — |
| concurrent_files(8) | 52.8 | 36 | 2.72× |
| concurrent_files(16), old default | 68.5 | 100 | 2.09× |

---

## After code default change (`default_ept_auto_partitions = 32`)

Commit `366b956a` on lasR `pre-devel`: auto-partition target is **32** when
`LASR_EPT_PARTITIONS` is unset (was `4 × concurrent_files` workers).

Verification on 112.5M-point AOI, **no env override**:

| Mode | Wall (s) | Chunks | Notes |
|------|----------|--------|-------|
| concurrent_files(16) | 52.3 | 36 | Matches old manual `PARTITIONS=32`; ~5 s slower than sweep median (network variance) |

Unit test on bundled fixture: auto chunk count matches `cpp_ept_partition_inspect(..., 32)`.

---

## How lasR chooses chunk count (relevant to these results)

When `concurrent_files(N)` runs on EPT with no explicit `chunk`/`LASR_EPT_PARTITIONS`:

| lasR version | Auto target | par8 chunks (typical) | par16 chunks (typical) |
|--------------|-------------|----------------------|------------------------|
| Before `366b956a` | `4 × N` | 32 → ~36 occupied | 64 → ~100 occupied |
| After `366b956a` | **32** (fixed) | ~36 | ~36 |

Env override: `LASR_EPT_PARTITIONS=<1..4096>` always wins.

Other knobs (defaults in parentheses):

- `LASR_EPT_PREFETCH` — tile prefetch depth per `EPTio` reader (**4**)
- `VSI_CACHE_SIZE` — GDAL shared cache (**64 MiB**)
- `Concurrent points` — **1** for `summarise()` (no inner parallelism)

---

## Recommended settings

```bash
# After lasR pre-devel >= 366b956a (default target = 32):
export CLAUDE_JOB_DIR=/path/to/work
Rscript scripts/detect_lasr_ept_aoi.R   # uses concurrent_files via LASR_EPT_PARALLEL

# Older lasR or explicit override:
export LASR_EPT_PARTITIONS=32
# concurrent_files(16) on large AOIs; concurrent_files(8) is a safe fallback
```

For AOIs ≲30M points, **sequential** or **concurrent_files(8)** is often enough;
parallel overhead can exceed benefit below ~25–30M points on this network path.

---

## Reproduce

```bash
export CLAUDE_JOB_DIR=/tmp/lidar_bench && mkdir -p "$CLAUDE_JOB_DIR"

# Full sweep (~30 min on 112M-point AOI)
export LASR_SWEEP_REPS=2
Rscript scripts/sweep_lasr_ept_params.R

# Quick seq vs par8 vs par16
export LASR_BENCH_REPS=3
export LASR_BENCH_FILE_CORES=16
Rscript scripts/bench_lasr_ept_acquisition.R
```
