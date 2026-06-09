#!/usr/bin/env python3
# #M6 SegmentAnyTree headless driver. Wraps the upstream folder-in/folder-out
# run_inference.sh into the single-file (input -> output) contract that
# scripts/model_runner.R::run_docker_arm expects.
#
# Input : ONE SAT-shaped PLY (lowercase x y z, ABSOLUTE UTM, raw-with-ground).
#         SAT self-localizes XY (pipeline_utm2local) and does NOT reproject or
#         normalize height -- feed UTM, never EPSG:3857.
# Output: the merged per-point LAS the pipeline writes under final_results/,
#         carrying the PredInstance extra dim (0 = non-tree, 1..N = tree id) at
#         absolute UTM Z. The R arm reads PredInstance -> reduce_instances
#         (max-Z per id) -> det_to_agl(ground_dtm.tif).
#
# Usage: run_segmentanytree.py <input.ply> <output.las>
#
# Runs INSIDE the sat-sm120-test container; uses only stdlib so it can be
# bind-mounted from the host (no rebuild) and executed by the container python3.
import os
import sys
import glob
import shutil
import subprocess
import tempfile

SAT_HOME = os.environ.get("SAT_HOME", "/opt/segment-any-tree")
# run_inference.sh hardcodes SCRIPT_DIR and conf/eval.yaml hardcodes
# checkpoint_dir, both to the upstream image layout /home/nibio/mutable-
# outside-world. Our sm_120 image clones to SAT_HOME instead, so reconcile both
# with a single symlink rather than patching upstream files.
LEGACY = "/home/nibio/mutable-outside-world"

# Upstream run_inference.sh runs eval.py with an ABSOLUTE --config-name, invalid
# under Hydra 1.3 (it strips the leading '/' and fails). Hydra 1.3 resolves the
# primary config by NAME across the search path, and the in-tree conf/ provider
# outranks an extra --config-dir, so it would load the UNMODIFIED conf/eval.yaml
# (default `data.fold`) instead of modify_eval.py's edited copy in DEST_DIR.
# Overwrite conf/eval.yaml with the modified copy, then load by bare name. Safe:
# the container is --rm, so mutating the in-tree config does not persist.
_EVAL_BAD = 'python3 eval.py --config-name "$DEST_DIR/eval.yaml"'
_EVAL_FIX = ('cp "$DEST_DIR/eval.yaml" "$SCRIPT_DIR/conf/eval.yaml" && '
             'python3 eval.py --config-name eval')


def _patch_eval_config(src_sh, dst_sh):
    with open(src_sh) as fh:
        s = fh.read()
    if _EVAL_BAD not in s:
        sys.exit("ERROR: run_inference.sh no longer contains the expected eval.py "
                 "call (%r); the Hydra --config-name patch needs review." % _EVAL_BAD)
    with open(dst_sh, "w") as fh:
        fh.write(s.replace(_EVAL_BAD, _EVAL_FIX))


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: run_segmentanytree.py <input.ply> <output.las>")
    inp = os.path.abspath(sys.argv[1])
    out = os.path.abspath(sys.argv[2])
    if not inp.lower().endswith(".ply"):
        sys.exit("ERROR: SegmentAnyTree ingests PLY; got %s" % inp)
    if not os.path.exists(inp):
        sys.exit("ERROR: input not found: %s" % inp)

    if not os.path.exists(LEGACY):
        os.makedirs(os.path.dirname(LEGACY), exist_ok=True)
        os.symlink(SAT_HOME, LEGACY)

    work = tempfile.mkdtemp(prefix="sat_")
    src = os.path.join(work, "input")
    dst = os.path.join(work, "results")
    os.makedirs(src)
    os.makedirs(dst)
    # Stable, '-'-free basename: fix_naming_of_input_files.py rewrites '-'->'_',
    # and the final_results basename is echoed from the input stem.
    shutil.copy(inp, os.path.join(src, "clip.ply"))

    env = dict(os.environ)
    # Put the compat dir FIRST so its usercustomize.py auto-applies the
    # stack-modernization shims (omegaconf missing-attr, torch.load weights_only,
    # multiprocessing spawn, cross-device indexing) to the `python3 eval.py` that
    # run_inference.sh spawns -- upstream files stay unmodified. With the image's
    # pinned deps (torch-geometric 2.3.1, hdbscan, pyarrow) these carry inference
    # end-to-end (#M6 bring-up, verified 2026-06-09).
    compat = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sat_compat")
    env["PYTHONPATH"] = "%s:%s:%s" % (compat, SAT_HOME, env.get("PYTHONPATH", ""))
    # run_inference.sh <source_dir> <dest_dir> <clean_output_dir>. Patch the one
    # broken line first: it calls eval.py with an ABSOLUTE --config-name, which
    # Hydra 1.3 strips to a relative path and fails to find. The documented way
    # to load a config from an out-of-tree dir is --config-dir + a bare name;
    # --config-dir ADDS to the search path, so the `visualization` default still
    # resolves from conf/. Patch a copy (the source image stays untouched).
    patched = os.path.join(work, "run_inference.sh")
    _patch_eval_config(os.path.join(SAT_HOME, "run_inference.sh"), patched)
    subprocess.run(
        ["bash", patched, src, dst, "true"],
        cwd=SAT_HOME, env=env, check=True,
    )

    # The merge writes final_results/<stem>.las; fall back to any .las under dst
    # (e.g. .laz vs .las extension drift upstream), preferring final_results.
    res = sorted(glob.glob(os.path.join(dst, "final_results", "*.la[sz]")))
    if not res:
        res = sorted(glob.glob(os.path.join(dst, "**", "*.la[sz]"), recursive=True))
    if not res:
        tree = []
        for root, _, files in os.walk(dst):
            for f in files:
                tree.append(os.path.join(root, f))
        sys.exit("ERROR: no .las/.laz under %s\ntree:\n%s" % (dst, "\n".join(tree)))
    if len(res) > 1:
        sys.stderr.write("WARN: %d result LAS; using %s\n" % (len(res), res[0]))
    shutil.copy(res[0], out)
    print("wrote %s from %s" % (out, res[0]))


if __name__ == "__main__":
    main()
