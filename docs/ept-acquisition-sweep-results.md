# lasR EPT acquisition parameter sweep

Benchmarks run on **USGS 3DEP** `CA_CarrHirzDeltaFires_2_2019` remote EPT  
(`https://s3-us-west-2.amazonaws.com/usgs-lidar-public/CA_CarrHirzDeltaFires_2_2019/ept.json`)

Pipeline timed: `reader_rectangles(AOI) + summarise()` with `-drop_class 7 18 -drop_withheld`.

Host: **16 OpenMP threads** (`available_threads=16`), lasR **0.21.0** (pre-devel build).

---

## AOI sizes

| Label | Bbox (EPSG:3857) | Approx. area | Points |
|-------|------------------|--------------|--------|
| Original | `-13578452,4988399,-13577792,4989059` | ~25 ha | 5.7M |
| 5× test | `-13578860,4987991,-13577384,4989467` | ~5× original | 31.3M |
| 4× larger | `-13579688,4987253,-13576556,4990205` | ~4× 5× test | 112.5M |

---

## Small AOI (~5.7M pts)

| Mode | Median wall | Chunks | vs seq |
|------|-------------|--------|--------|
| sequential | 8.2 s | 1 | — |
| concurrent_files(8) | 11.6 s | 64 | 0.70× (slower) |
| concurrent_files(16) | 11.7 s | 64 | 0.70× (slower) |

Parallel hurt on this small AOI: partition overhead dominates; network-bound S3 reads.

---

## Medium AOI (~31M pts, 5× original)

| Mode | Median wall | Chunks | vs seq |
|------|-------------|--------|--------|
| sequential | 43 s | 1 | — |
| concurrent_files(8) | 25.1 s | 36 | 1.7× |
| concurrent_files(16) default | 36.2 s | 100 | 1.2× |

With **`LASR_EPT_PARTITIONS=32`** and 16 workers: **~47 s**, **36 chunks** (matches par8 chunk count).

---

## Large AOI (~112M pts, 4× medium) — full parameter sweep

2 reps per config unless noted. All returned **112,519,502 points**.

| Config | Workers | Partitions | Prefetch | Cache (MB) | Median (s) | Chunks | Mpts/s | vs seq |
|--------|---------|------------|----------|------------|------------|--------|--------|--------|
| sequential | 1 | default | default | default | **138.4** | 1 | 0.81 | — |
| **par8_default** | 8 | default (32) | default | default | **52.6** | 36 | 2.14 | **2.63×** |
| par16_default | 16 | default (64) | default | default | 67.3 | 100 | 1.67 | 2.06× |
| **par16_parts32** | 16 | **32** | default | default | **47.3** | 36 | 2.38 | **2.93×** |
| par16_parts36 | 16 | 36 | default | default | 47.7 | 36 | 2.36 | 2.90× |
| par16_parts64 | 16 | 64 | default | default | 67.3 | 100 | 1.67 | 2.06× |
| par8_prefetch8 | 8 | default | 8 | default | 52.4 | 36 | 2.15 | 2.64× |
| par8_prefetch16 | 8 | default | 16 | default | 52.2 | 36 | 2.15 | 2.65× |
| par8_prefetch32 | 8 | default | 32 | default | 52.6 | 36 | 2.14 | 2.63× |
| par16_parts32_prefetch8 | 16 | 32 | 8 | default | 47.8 | 36 | 2.35 | 2.90× |
| par16_parts32_prefetch16 | 16 | 32 | 16 | default | 48.7 | 36 | 2.31 | 2.84× |
| par16_parts32_prefetch32 | 16 | 32 | 32 | default | 48.4 | 36 | 2.33 | 2.86× |
| par8_cache256 | 8 | default | default | 256 | 51.9 | 36 | 2.17 | 2.66× |
| par8_cache1024 | 8 | default | default | 1024 | 52.0 | 36 | 2.16 | 2.66× |
| par16_parts32_cache256 | 16 | 32 | default | 256 | 47.9 | 36 | 2.35 | 2.90× |
| par16_parts32_cache1024 | 16 | 32 | default | 1024 | 48.1 | 36 | 2.34 | 2.88× |

Raw CSV: run `scripts/sweep_lasr_ept_params.R` → `$CLAUDE_JOB_DIR/sweep_lasr_ept_params.csv`

---

## Key findings

1. **Partition count dominates.** Default `4 × workers` at 16 cores → target 64 → **100 chunks**; forcing **32** → **36 chunks** and ~30% faster than default par16.

2. **Best config on large AOI:** `concurrent_files(16)` + partition target **32** → **47.3 s** (2.93× vs sequential). Beats par8 default (52.6 s) by ~10%.

3. **Prefetch (`LASR_EPT_PREFETCH`) and `VSI_CACHE_SIZE`:** negligible effect (~0–1 s) on this workload.

4. **Not OOM.** Peak RSS ~1 GB at par16; failures at high chunk counts were heap corruption / segfaults (intermittent), not memory exhaustion.

5. **Code change (lasR pre-devel):** default auto-partition target changed from `4 × concurrent_files` to fixed **32** when `LASR_EPT_PARTITIONS` is unset (`366b956a` on pre-devel, `78282bb6` on `parallel-ept-acquisition`).

---

## Reproduce

```bash
export CLAUDE_JOB_DIR=/tmp/lidar_bench
export LASR_EPT_AOI_BBOX="-13579688,4987253,-13576556,4990205"
export LASR_SWEEP_REPS=2
Rscript scripts/sweep_lasr_ept_params.R
```

Single comparison:

```bash
export LASR_EPT_AOI_BBOX="-13579688,4987253,-13576556,4990205"
export LASR_BENCH_FILE_CORES=16
Rscript scripts/bench_lasr_ept_acquisition.R
```

Override partition target:

```bash
export LASR_EPT_PARTITIONS=32
```
