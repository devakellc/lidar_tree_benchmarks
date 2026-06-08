#!/usr/bin/env bash
# #I5 mirror: fetch model weights/images into a local store and verify SHA256
# against a VERSIONED manifest (gpu/CHECKSUMS.sha256). First run records the
# full hash; later runs hard-fail on any mismatch (manual verify against the
# manifest -- the store layout is not `sha256sum -c`-shaped). GitHub-release
# assets go through authenticated `gh` (the unauthenticated github.com download
# path 504s here); other sources use curl.
#
# Usage:
#   bash gpu/mirror_weights.sh        # TreeisoNet only (the #M7 detection arm)
#   bash gpu/mirror_weights.sh ALL    # + SegmentAnyTree (#M6) + ForestFormer3D (#M8)
set -euo pipefail
WHICH="${1:-treeisonet}"
HERE="$(cd "$(dirname "$0")" && pwd)"
STORE="$HERE/store"; MAN="$HERE/CHECKSUMS.sha256"
mkdir -p "$STORE/treeaibox" "$STORE/segmentanytree" "$STORE/forestformer3d"
touch "$MAN"
verify() {  # dest sha_prefix(sanity; "" to skip)
  local dest="$1" pre="$2" base got exp
  base="$(basename "$dest")"; got="$(sha256sum "$dest" | cut -d' ' -f1)"
  if [ -n "$pre" ]; then case "$got" in "$pre"*) : ;;
    *) echo "FAIL prefix $base: $got !~ $pre"; exit 1;; esac; fi
  if grep -q "  $base\$" "$MAN"; then
    exp="$(grep "  $base\$" "$MAN" | cut -d' ' -f1)"
    [ "$got" = "$exp" ] || { echo "FAIL manifest $base: $got != $exp"; exit 1; }
  else echo "$got  $base" >> "$MAN"; fi
  echo "OK $base ($got)"
}
gh_get() {  # asset dest sha_prefix  (TreeAIBox release asset via authenticated gh)
  local asset="$1" dest="$2" dir; dir="$(dirname "$dest")"
  if [ ! -f "$dest" ]; then echo "GET $asset"
    gh release download v1.0 --repo NRCan/TreeAIBox --pattern "$asset" --dir "$dir" --clobber
    [ "$dir/$asset" = "$dest" ] || mv "$dir/$asset" "$dest"
  fi
  verify "$dest" "$3"
}
url_get() {  # url dest sha_prefix  (curl for non-release sources)
  local url="$1" dest="$2"
  [ -f "$dest" ] || { echo "GET $(basename "$dest")"; curl -fL --retry 3 "$url" -o "$dest"; }
  verify "$dest" "$3"
}
# --- TreeisoNet ALS (the #M7 detection arm; always fetched) ---
gh_get "treeisonet_als_reclamation_treeloc_esegformer3D_128_10cm_GPU4GB.pth" \
       "$STORE/treeaibox/als_treeloc.pth" c6e15a21
gh_get "treeisonet_als_reclamation_treeoff_esegformer3D_128_10cm_GPU4GB.pth" \
       "$STORE/treeaibox/als_treeoff.pth" 6cf5d890
if [ "$WHICH" = "ALL" ]; then
  url_get "https://media.githubusercontent.com/media/SmartForest-no/SegmentAnyTree/main/model_file/PointGroup-PAPER.pt" \
          "$STORE/segmentanytree/PointGroup-PAPER.pt" 0b4d74b4
  url_get "https://zenodo.org/api/records/16742708/files/clean_forestformer.zip/content" \
          "$STORE/forestformer3d/clean_forestformer.zip" ""
fi
echo "manifest -> $MAN"; cat "$MAN"
# NOTE: TreeisoNet config JSONs are copied from the cloned repo in
# gpu/setup_treeisonet_env.sh (upstream names vary in punctuation -> glob).
# The SegmentAnyTree Docker image is pulled in #M6.
