#!/usr/bin/env python3
# Batch variant of run_segmentanytree.py. Runs the upstream SAT folder-in /
# folder-out inference once for a directory of PLY clips and copies every merged
# final LAS/LAZ to the requested output directory.
import glob
import os
import shutil
import subprocess
import sys
import tempfile

from run_segmentanytree import SAT_HOME, LEGACY, _patch_eval_config


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: run_segmentanytree_batch.py <input_dir> <output_dir>")
    inp = os.path.abspath(sys.argv[1])
    out = os.path.abspath(sys.argv[2])
    if not os.path.isdir(inp):
        sys.exit("ERROR: input directory not found: %s" % inp)
    plys = sorted(glob.glob(os.path.join(inp, "*.ply")))
    if not plys:
        sys.exit("ERROR: no .ply files under %s" % inp)
    os.makedirs(out, exist_ok=True)

    if not os.path.exists(LEGACY):
        os.makedirs(os.path.dirname(LEGACY), exist_ok=True)
        os.symlink(SAT_HOME, LEGACY)

    work = tempfile.mkdtemp(prefix="sat_batch_")
    src = os.path.join(work, "input")
    dst = os.path.join(work, "results")
    os.makedirs(src)
    os.makedirs(dst)
    for ply in plys:
        name = os.path.basename(ply).replace("-", "_")
        shutil.copy(ply, os.path.join(src, name))

    env = dict(os.environ)
    compat = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sat_compat")
    env["PYTHONPATH"] = "%s:%s:%s" % (compat, SAT_HOME, env.get("PYTHONPATH", ""))
    patched = os.path.join(work, "run_inference.sh")
    _patch_eval_config(os.path.join(SAT_HOME, "run_inference.sh"), patched)
    subprocess.run(["bash", patched, src, dst, "true"], cwd=SAT_HOME, env=env, check=True)

    res = sorted(glob.glob(os.path.join(dst, "final_results", "*.la[sz]")))
    if not res:
        res = sorted(glob.glob(os.path.join(dst, "**", "*.la[sz]"), recursive=True))
    if not res:
        tree = []
        for root, _, files in os.walk(dst):
            for f in files:
                tree.append(os.path.join(root, f))
        sys.exit("ERROR: no .las/.laz under %s\ntree:\n%s" % (dst, "\n".join(tree)))
    for path in res:
        name = os.path.basename(path)
        stem, ext = os.path.splitext(name)
        if stem.endswith("_out"):
            name = stem[:-4] + ext
        shutil.copy(path, os.path.join(out, name))
    print("wrote %d LAS/LAZ files to %s" % (len(res), out))


if __name__ == "__main__":
    main()
