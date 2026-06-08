# GPU Runner: Docker Backend + LAZ↔PLY Bridge — Design

**Status.** Approved design (brainstorming output, revised after a Codex
design review). Ready for implementation. Last updated 2026-06-08. Implements
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

These were settled during brainstorming and a follow-up Codex design review:

| Decision | Choice | Consequence |
|---|---|---|
| PLY I/O location | **Hand-rolled binary PLY in base R**, inside `io_bridge.R` | No new deps; the bridge stays unit-testable in pure R (no venv/Docker), matching the repo's GPU-free-bridge property. |
| Coordinate contract | **`coord_type` + `offset`**, offset returned and restored | float64 absolute is lossless but **float32 + absolute UTM northing (~4.1e6) loses ~0.5 m** (float32 ulp). So `float` forces a local offset; the offset is returned on write and re-added on read. See §2.1. |
| PLY schema | **Configurable property names + dtype + constant columns** | SAT writes lowercase `x/y/z/intensity` + `semantic_seg`/`treeID`; the writer must express those, not only LAS-native `Intensity`. The *mechanism* is here; each model's exact recipe stays in its arm. |
| Docker output | **Contract-owned `reader=`** on one runner | Default reader = the same `x y z` CSV parse as the venv arm; #M6/#M8 pass `read_instances_ply`/`read_instances_laz`. The runner owns failure/missing/exception/contract; the reader returns `NULL` on schema failure (skip cell), 0-row only when legitimately empty. |
| Bind mounts | **Identity mounts** (`-v abs:abs`), absolute paths only | In-container paths == host paths ⇒ zero translation. `normalizePath` + `shQuote` as in `run_python_arm` (paths with spaces/parens). |
| Crash contract | **Identical to the venv backend** | Non-zero exit OR missing output → `NULL`; wrong-schema → `NULL`; valid header, no rows → 0-row frame. Never a fake 0-row from a stale/partial file. |

---

## 2. Component: `io_bridge.R` LAZ↔PLY bridge

A minimal **binary little-endian** PLY reader/writer (base R `readBin`/
`writeBin`/`writeChar`), **vertex-scalar only** — no ascii branch, no `face`
element, no list properties (none of the models emit them). The header lists
properties by name and ends with `end_header\n`; the body is written
**interleaved per vertex in header-property order** (not column blocks).

### 2.1 Coordinate contract (load-bearing)

`property double` preserves absolute UTM losslessly, but `property float` does
not — a UTM northing near 4.1e6 has a float32 ulp of ~0.5 m. Because SAT runs
UTM→local preprocessing and FF3D recenters then casts to float32, the bridge
must make the offset explicit:

- `laz_to_ply(..., coord_type = "double")` → `offset = c(0, 0, 0)`, absolute
  coords (default; lossless).
- `laz_to_ply(..., coord_type = "float")` → `offset = floor(apply(XYZ, 2, min))`
  by default, so stored coords are small and float32-safe.
- `laz_to_ply` **returns the `offset`** (length-3 numeric). The arm holds it and
  passes it to the reader, which **adds it back** so detections score in the
  input UTM frame. The arm's in-container driver is required to emit output in
  the same frame it received (it MUST NOT silently re-offset its output); the
  documented contract is the alternative to per-model offset archaeology.

### 2.2 Functions

- **`laz_to_ply(laz_path, ply_path, fields = character(), props = NULL,
  coord_type = c("double", "float"), offset = NULL)`** — `readLAS` → subtract
  `offset` → write `element vertex N` with `x/y/z` of `coord_type` plus any
  `fields`. `props` is an optional named spec that controls **output property
  name, dtype, and constant columns** (e.g. `intensity` from `Intensity`,
  `semantic_seg = 1L`, `treeID = 0L`) so a model-shaped PLY is expressible. A
  small default name→PLY-dtype map keeps LAS-native widths
  (`Intensity`→`ushort`, `ReturnNumber`/`NumberOfReturns`/`Classification`→
  `uchar`, `gpstime`→`double`); unknown numerics fall back to `float`. Returns
  `list(n, offset)` invisibly.
- **`read_ply(ply_path, offset = c(0, 0, 0))`** — binary vertex-scalar reader →
  base `data.frame`, one column per PLY property (named by the header); adds
  `offset` to `x/y/z`. Supports the standard scalar dtype set.
- **`ply_to_laz(ply_path, laz_path, offset = c(0, 0, 0))`** — the inverse
  converter (issue asks for "and inverse"): `read_ply` → `LAS` → `writeLAS`,
  recognized fields mapped back to LAS columns. Minimal; vertex-scalar only.
- **`read_instances_ply(ply_path, id_field = "treeID", offset = c(0, 0, 0),
  x = "x", y = "y", z = "z")`** — `read_ply` → strict reduce → apex
  `det(x, y, z)`. **Strict**: if `id_field` is absent from the PLY it returns
  `NULL` (schema failure → the runner skips the cell), never a fake 0-row.
  Mirrors `read_instances_laz`, which gets the same strict id-field check (an
  empty cloud is still a legitimate 0-row; a missing id field is `NULL`).

`instances_to_det`, `reduce_instances`, and `assert_detection_contract` are
reused unchanged for the reduction itself (`id == 0` OR `NA` is unassigned and
dropped, apex = max-Z per id). The strict-vs-empty distinction is added in the
model-output readers, leaving the general `instances_to_det` utility untouched.

## 3. Component: `model_runner.R` Docker backend

**`run_docker_arm(image, input, out_csv, extra = character(),
cmd = character(), mounts = NULL, gpus = "all", docker = "docker",
timeout = 1800, label = NULL, reader = NULL)`**

- Unlink any stale `out_csv` first (never read a stale file).
- **Identity-mount** `normalizePath(dirname(input))` +
  `normalizePath(dirname(out_csv))`, plus every entry in `mounts` (host
  weight/config dirs), as `-v abs:abs`. Relative paths are normalized to
  absolute so `-v` is always valid.
- Assemble `docker run --rm [--gpus all] -v … <image> [cmd…] <input>
  <out_csv> <extra…>`, `shQuote` every token, and `system2` it with the same
  `tryCatch` + exit-status discipline as `run_python_arm`. `gpus = NULL` omits
  `--gpus`; `cmd` overrides the in-container command (default: the image's own
  entrypoint/CMD).
- **The runner owns the contract.** On non-zero exit / missing output → `NULL`.
  Otherwise the reader is applied inside `tryCatch` (a throwing reader → `NULL`);
  a `NULL` reader result → `NULL` (schema failure, skip); a non-NULL result is
  validated with `assert_detection_contract`. Default reader = the `x y z` CSV
  parse (so a CSV-emitting container behaves identically to the venv arm).

### DRY refactor

Extract the CSV-parse-with-contract block (valid `x y z` → `det`; wrong-schema
→ `NULL`; header-only → 0-row) into one helper `.read_detection_csv()`, used by
**both** `run_python_arm` and `run_docker_arm` as the default reader. Behavior
is unchanged and guarded by the existing `test-model-runner.R` cases.

---

## 4. Tests (pure-R, existing testthat style)

### `tests/testthat/test-io-bridge.R` (additions)

- `laz_to_ply` (double) → `read_ply` round-trips X/Y/Z within tol and carries
  `Intensity` + `ReturnNumber` at their LAS-native widths.
- A float64 absolute-UTM losslessness guard (UTM 11N easting/northing survive to
  sub-mm) — the PLY analogue of the existing CRS/units LAS round-trip.
- `coord_type = "float"` with `offset` round-trips: stored coords are small,
  and `read_ply(offset=)` restores absolute UTM within the float32 tolerance
  (and **without** the offset the error is ~0.5 m — proving the contract earns
  its keep).
- `props` emits a **SAT-shaped PLY** (lowercase `x/y/z/intensity`, constant
  `semantic_seg`/`treeID`) — asserts header property names/dtypes.
- `read_instances_ply` reduces a labeled PLY (`treeID`) to per-id apexes,
  dropping id 0; **returns `NULL`** when `id_field` is absent (strict).
- `laz_to_ply` → `ply_to_laz` → `readLAS` round-trips XYZ (the inverse).
- **Skippable** cross-reader test: `skip_if` Python `plyfile` is unavailable,
  else read the emitted PLY with `plyfile` and assert X/Y/Z + a field match —
  breaking the self-round-trip circularity.

### `tests/testthat/test-model-runner.R` (additions)

A tiny **fake `docker`** script that **records its argv** to a file and then
emulates the daemon (strips `run`, `--rm`, `--gpus <v>`, repeated `-v <v>`, and
the image token, then `exec`s the rest), reusing the same arm scripts as the
venv tests:

- valid CSV → `x y z` det, **and the recorded argv equals**
  `run --rm --gpus all -v host:host … image cmd input out extra`;
- `gpus = NULL` omits `--gpus`; repeated `mounts` each appear as `-v`;
- an `input` path **with spaces** survives (shQuote + identity mount);
- non-zero container exit → `NULL` even with a stale CSV present;
- missing output → `NULL`; valid header, no rows → 0-row frame;
- a `reader=` that returns `NULL` (schema failure) → `NULL`; a throwing
  `reader=` → `NULL`; a valid custom `reader=` → its det.

The fake `docker` exercises the real argument assembly without a daemon or GPU;

## M6/#M8 still owe one real Docker smoke when they wire their drivers

---

### 5. Scope boundary (non-goals)

- The #M6/#M8 **headless drivers** and container entrypoints, real weights, and
  GPU inference runs — those arms' own issues.
- Each model's **exact PLY recipe** (which constant columns, which `id_field`,
  float32-vs-float64 choice) — expressed *through* `props`/`coord_type`/
  `id_field` where each arm is wired. The bridge ships the mechanism + one
  SAT-shaped example test, not a per-model preset.
- Folding a Docker arm into `analyze_model_benchmark.R` or any result doc.

### 6. Doc hygiene (in-scope cleanup)

- Refresh the now-stale "added in their issue / added with #M6" comments at the
  top of `io_bridge.R` and `model_runner.R` to state the backend/bridge now
  exist.
- Tick the `gpu/segmentanytree-sm120/README.md` next-steps lines for "Wire
  Docker backend in `scripts/model_runner.R`" and "LAZ→PLY conversion".
