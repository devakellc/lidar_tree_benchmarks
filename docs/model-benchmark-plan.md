# NEON Model Benchmark — Triage & Standing Protocol (#A0)

Decision record for the deep-model benchmark. Source of truth for the runnable
set and the rules every arm must follow. Derived from
[the design spec](superpowers/specs/2026-06-07-neon-deep-model-benchmark-design.md).

## Runnable set

| Model | Disposition | Rung scope | Input variant |
|---|---|---|---|
| AMS3D (`crownsegmentr`) | run — walking skeleton | full ladder | normalized |
| lidRplugins (`lmfauto`/`multichm`/`ptree`) | run — competitor | full ladder | normalized |
| SegmentAnyTree | run | full ladder | raw-with-ground |
| TreeisoNet-ALS | run | full ladder | raw-with-ground |
| ForestFormer3D | run | native + 8 only | raw-with-ground |
| ForAINet | defer | — | — |
| HFC | defer (agreement) | — | — |
| TreeLearn | drop | — | — |
| Dersch graph-cut | drop | — | — |
| DeepForest (RGB) | optional reference | n/a (density-invariant) | RGB + CHM |

Baselines: existing CHM-VWF and Li 2012 arms.

## Zero-shot protocol

- No fine-tuning. Published weights / classical defaults only.
- Per-arm knobs MAY be set from prior literature/allometry; they MUST NOT be
  tuned on SOAP scoring outcomes.
- Every non-default knob is recorded in the run's ledger (the
  `*_ledger.csv` written next to each arm's results), with the value and its
  literature source.
- AMS3D knobs to record: `crown_diameter_to_tree_height`,
  `crown_length_to_tree_height`, `segment_crowns_only_above`.

## Weights-mirror policy

- Mirror every model weight/image to a project store before first use; record
  SHA256 + source URL. Upstream links (TreeisoNet personal server, ForAINet
  Dropbox) are single points of failure.
- Tracked in the GPU-infra plan (#I5), not here.

## Density framing

Rungs are all-return pts/m². SOAP native ≈ 20 all-return / ≈ 12 first-return.
Deep models clear their floor (10) at native; feasibility risk is the 4/2/1
rungs and the discrete-return-vs-ULS structural gap, not native sparsity.
