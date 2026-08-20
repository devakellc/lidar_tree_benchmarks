#!/usr/bin/env python
# #X1 DeepForest runner: an RGB GeoTIFF tile -> tree-crown box centroids in the
# tile's CRS (UTM), written as x,y,score,w_m,h_m CSV. The density-INVARIANT
# optical arm: RGB detection does not degrade with LiDAR point density.
#   python run_deepforest.py <tile.tif> <out.csv> [patch_size] [patch_overlap]
import sys, numpy as np, pandas as pd, rasterio
from rasterio.transform import xy
from deepforest import main

tile, out = sys.argv[1], sys.argv[2]
patch = int(sys.argv[3]) if len(sys.argv) > 3 else 400
overlap = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05

m = main.deepforest()
m.load_model("weecology/deepforest-tree")
boxes = m.predict_tile(path=tile, patch_size=patch, patch_overlap=overlap)
if boxes is None or len(boxes) == 0:
    pd.DataFrame(columns=["x", "y", "score", "w_m", "h_m"]).to_csv(out, index=False)
    print("deepforest: 0 boxes -> %s" % out); sys.exit(0)

with rasterio.open(tile) as src:
    T = src.transform; res = abs(T.a)                   # m / pixel

# DeepForest predict_tile returns boxes in PIXEL space (geometry or xmin..ymax),
# NOT georeferenced -- so always map pixel box centroids through the raster
# transform to the tile CRS. xmin/xmax are columns (x-pixel), ymin/ymax rows.
g = getattr(boxes, "geometry", None)
if g is not None and g.notna().any():
    bn = boxes.geometry.bounds
    pminx, pminy, pmaxx, pmaxy = bn.minx.values, bn.miny.values, bn.maxx.values, bn.maxy.values
else:
    pminx, pminy = boxes.xmin.values, boxes.ymin.values
    pmaxx, pmaxy = boxes.xmax.values, boxes.ymax.values
col = (pminx + pmaxx) / 2.0; row = (pminy + pmaxy) / 2.0
xs, ys = xy(T, row, col)                                 # pixel -> CRS x,y
xs, ys = np.asarray(xs, float), np.asarray(ys, float)
w_m = (pmaxx - pminx) * res; h_m = (pmaxy - pminy) * res

score = boxes.score.values if "score" in boxes else np.ones(len(xs))
pd.DataFrame({"x": xs, "y": ys, "score": score, "w_m": w_m, "h_m": h_m}).to_csv(out, index=False)
print("deepforest: %d crown boxes -> %s" % (len(xs), out))
