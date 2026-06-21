#!/usr/bin/env python
# #X3 Detectree2 runner: a (small, plot-sized) georeferenced RGB crop -> tree-crown
# POLYGON centroids + equivalent crown diameter, in the crop's CRS, as
# x,y,score,d_eq,area CSV. Detectree2 is a Detectron2 Mask R-CNN crown segmenter;
# its masks (vs DeepForest's boxes) let the optical modality contribute crown
# WIDTH, not just detection. CPU inference (Detectron2 built CPU-only against
# torch 2.12). NOTE: the pretrained weights are tropical-forest-trained, so
# CA-conifer transfer is a tested unknown -- run on plot crops, not whole tiles
# (each 40 m sub-tile is a full Mask R-CNN pass, too slow for a 1 km2 mosaic).
#   python run_detectree2.py <crop.tif> <out.csv> <model.pth>
import sys, os, shutil, numpy as np, pandas as pd
from detectree2.preprocessing.tiling import tile_data
from detectree2.models.train import setup_cfg
from detectree2.models.predict import predict_on_data
from detectree2.models.outputs import project_to_geojson, stitch_crowns, clean_crowns
from detectron2.engine import DefaultPredictor

crop, out_csv, model = sys.argv[1], sys.argv[2], sys.argv[3]
work = sys.argv[4] if len(sys.argv) > 4 else "/tmp/dt2_work"
if os.path.isdir(work): shutil.rmtree(work)
os.makedirs(work, exist_ok=True)

EMPTY = lambda: pd.DataFrame(columns=["x", "y", "score", "d_eq", "area"]).to_csv(out_csv, index=False)
try:
    tile_data(crop, work, buffer=10, tile_width=40, tile_height=40)   # CRS metres
    cfg = setup_cfg(update_model=model)
    cfg.MODEL.DEVICE = "cpu"
    preds = os.path.join(work, "predictions") + "/"
    geo = os.path.join(work, "predictions_geo") + "/"
    predict_on_data(work + "/", out_folder=preds, predictor=DefaultPredictor(cfg))
    project_to_geojson(work + "/", preds, geo)
    crowns = stitch_crowns(geo, 1)
    if crowns is None or len(crowns) == 0:
        EMPTY(); print("detectree2: 0 crowns -> %s" % out_csv); sys.exit(0)
    crowns = clean_crowns(crowns, iou_threshold=0.7, confidence=0.2, verbose=False)
    crowns = crowns[crowns.geometry.notna() & crowns.geometry.is_valid & ~crowns.geometry.is_empty]
except Exception as e:
    EMPTY(); sys.stderr.write("detectree2 failed: %s\n" % e); print("detectree2: error"); sys.exit(0)

if len(crowns) == 0:
    EMPTY(); print("detectree2: 0 crowns -> %s" % out_csv); sys.exit(0)
c = crowns.geometry.centroid
area = crowns.geometry.area.values                       # m^2 (CRS in metres)
d_eq = 2.0 * np.sqrt(area / np.pi)
conf = crowns["Confidence_score"].values if "Confidence_score" in crowns else np.ones(len(area))
pd.DataFrame({"x": c.x.values, "y": c.y.values, "score": conf, "d_eq": d_eq,
             "area": area}).to_csv(out_csv, index=False)
print("detectree2: %d crown polygons -> %s" % (len(area), out_csv))
