# GPU Runner: Docker Backend + LAZ↔PLY Bridge — Design

**Status.** Approved design (brainstorming output). Ready for implementation
planning. Last updated 2026-06-08. Implements
[issue #19](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/19).

**Scope.** Build the two pieces of GPU-arm infrastructure that the
TreeisoNet-first plan
([`2026-06-08-gpu-arm-infra-m7-first.md`](../plans/2026-06-08-gpu-arm-infra-m7-first.md))
explicitly **deferred** to the Docker-bound arms: a **Docker backend** in
`scripts/model_runner.R` and a **LAZ↔PLY bridge** in `scripts/io_bridge.R`.
Both are prerequisites for #M6 (SegmentAnyTree) and #M8 (ForestFormer3D), which
ingest/emit PLY and run inside a CUDA container. This issue ships the
**infrastructure plus unit tests only** — not the model arms themselves.

---

## 1. Locked decisions

These were settled during brainstorming and bound everything below:

| Decision | Choice | Consequence |
|---|---|---|
| PLY I/O location | **Hand-rolled binary PLY in base R**, inside `io_bridge.R` | No new deps; the bridge stays unit-testable in pure R (no venv/Docker), matching the repo's GPU-free-bridge property. Model-specific PLY dialects stay a per-arm concern. |
| Coordinates | **float64 absolute UTM**, no offset | `property double x/y/z` round-trips UTM losslessly, so no offset bookkeeping. Models that internally re-center (FF3D, SAT both do) are unaffected. |
| Docker output | **`reader=` injection** on one runner | Default reader = the same `x y z` CSV parse as the venv arm; #M6/#M8 pass `read_instances_ply`/`read_instances_laz` for the labeled-cloud path. One runner, not two. |
| Bind mounts | **Identity mounts** (`-v dir:dir`) | In-container paths == host paths ⇒ zero path translation; the entrypoint receives real paths valid inside the container. |
| Crash contract | **Identical to the venv backend** | Non-zero exit OR missing output → `NULL`; wrong-schema → `NULL`; valid header, no rows → 0-row frame. Never a fake 0-row from a stale/partial file. |

---

## 2. Component: `io_bridge.R` LAZ↔PLY bridge

A minimal **binary little-endian** PLY reader/writer (base R `readBin`/
`writeBin`/`writeChar`), plus an ascii read fallback. The header lists
properties by name, so the canonical Hugues-Thomas reader that FF3D uses (and
the `plyfile`/`plyutils` family in general) can parse the output.

### Functions

- **`laz_to_ply(laz_path, ply_path, fields = character())`** — `readLAS` →
  write `element vertex N` with `property double x/y/z` plus any requested
  `fields`. A small name→PLY-dtype map keeps LAS-native widths
  (`Intensity`→`ushort`, `ReturnNumber`/`NumberOfReturns`/`Classification`→
  `uchar`, `gpstime`→`double`); unknown numeric fields fall back to `float`.
  Returns the point count invisibly.
- **`read_ply(ply_path)`** — generic reader → base `data.frame`, one column per
  PLY property (named by the header). Supports the standard PLY dtype set and an
  ascii fallback.
- **`ply_to_laz(ply_path, laz_path)`** — the inverse: `read_ply` → `LAS` →
  `writeLAS`. Recognized fields map back to LAS columns.
- **`read_instances_ply(ply_path, id_field = "treeID", x = "x", y = "y",
  z = "z")`** — `read_ply` → `instances_to_det` → apex `det(x, y, z)`. The
  "parse the instance-labeled PLY output back" deliverable, mirroring the
  existing `read_instances_laz` (which already covers a labeled-LAS output).

`instances_to_det`, `reduce_instances`, and `assert_detection_contract` are
reused unchanged: `id == 0` OR `NA` is unassigned and dropped, apex = max-Z per
id.

## 3. Component: `model_runner.R` Docker backend

**`run_docker_arm(image, input, out_csv, extra = character(),
cmd = character(), mounts = NULL, gpus = "all", docker = "docker",
timeout = 1800, label = NULL, reader = NULL)`**

- Unlink any stale `out_csv` first (never read a stale file).
- **Identity-mount** `dirname(input)` + `dirname(out_csv)`, plus every entry in
  `mounts` (host weight/config dirs), as `-v dir:dir`.
- Assemble `docker run --rm [--gpus all] -v … <image> [cmd…] <input>
  <out_csv> <extra…>` and `system2` it with the same `tryCatch` +
  exit-status discipline as `run_python_arm`.
- Default reader parses the `x y z` CSV exactly like the venv arm (so a
  CSV-emitting container behaves identically). A supplied `reader(out_csv)` is
  used instead for the labeled-cloud path.
- `gpus = NULL` omits `--gpus` (CPU/test). `cmd` overrides the in-container
  command (default: rely on the image's own entrypoint/CMD).

### DRY refactor

Extract the CSV-parse-with-contract block (valid `x y z` → `det`; wrong-schema
→ `NULL`; header-only → 0-row) into one helper `.read_detection_csv()`, used by
**both** `run_python_arm` and `run_docker_arm`. Behavior is unchanged and
guarded by the existing `test-model-runner.R` cases.

---

## 4. Tests (pure-R, existing testthat style)

### `tests/testthat/test-io-bridge.R` (additions)

- `laz_to_ply` → `read_ply` round-trips X/Y/Z within tol and carries
  `Intensity` + `ReturnNumber`.
- `laz_to_ply` → `ply_to_laz` → `readLAS` round-trips XYZ (the inverse).
- `read_instances_ply` reduces a labeled PLY (`treeID`) to per-id apexes,
  dropping id 0.
- A float64 absolute-UTM losslessness guard (UTM 11N easting/northing survive
  the PLY round-trip to sub-mm), the PLY analogue of the existing CRS/units LAS
  round-trip test.

### `tests/testthat/test-model-runner.R` (additions)

A tiny **fake `docker`** script (a shell stub that strips `run`, `--rm`,
`--gpus <v>`, repeated `-v <v>`, and the image token, then `exec`s the rest)
drives `run_docker_arm`, reusing the same arm scripts as the venv tests:

- valid CSV → `x y z` det;
- non-zero container exit → `NULL` even with a stale CSV present;
- missing output → `NULL`;
- valid header, no rows → 0-row frame;
- `reader=` injection routes a non-CSV output through a custom reader.

The fake `docker` exercises the real argument assembly (mounts, `--gpus`, image,
command) without needing a daemon or the GPU.

---

## 5. Scope boundary (non-goals)

- The #M6/#M8 **headless drivers** and container entrypoints, real weights, and
  GPU inference runs — those arms' own issues.
- Any **model-specific PLY dialect** (e.g. FF3D's dummy `semantic_seg`/`treeID`
  test columns, SAT's `PredInstance` extra-dim naming) — handled where each arm
  is wired, by passing `fields=`/`id_field=`.
- Folding a Docker arm into `analyze_model_benchmark.R` or any result doc.

## 6. Doc hygiene (in-scope cleanup)

- Refresh the now-stale "added in their issue / added with #M6" comments at the
  top of `io_bridge.R` and `model_runner.R` to state the backend/bridge now
  exist.
- Tick the `gpu/segmentanytree-sm120/README.md` next-steps lines for "Wire
  Docker backend in `scripts/model_runner.R`" and "LAZ→PLY conversion".
