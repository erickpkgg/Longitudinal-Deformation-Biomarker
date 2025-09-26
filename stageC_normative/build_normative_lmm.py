
import argparse, pandas as pd
from pathlib import Path
ap=argparse.ArgumentParser()
ap.add_argument("--jobs",required=True); ap.add_argument("--roi_csv",required=True); ap.add_argument("--out_dir",required=True)
args=ap.parse_args()
roi=pd.read_csv(args.roi_csv); out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
cn=roi[roi["DX"]=="CN"].copy()
stats=cn.groupby("ROI")["value"].agg(mu_CN="mean", sd_CN="std").reset_index()
stats.to_csv(out/"normative_models.csv",index=False); print("WROTE",out/"normative_models.csv")
