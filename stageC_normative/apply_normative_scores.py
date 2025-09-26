
import argparse, pandas as pd, numpy as np
from pathlib import Path
ap=argparse.ArgumentParser()
ap.add_argument("--jobs",required=True); ap.add_argument("--roi_csv",required=True); ap.add_argument("--models",required=True); ap.add_argument("--out_csv",required=True)
args=ap.parse_args()
roi=pd.read_csv(args.roi_csv); models=pd.read_csv(args.models)
df=roi.merge(models,on="ROI",how="left"); df["z"]=(df["value"]-df["mu_CN"])/df["sd_CN"]
comp=df.groupby(["Subject","JobID"]).agg(IDE_global=("z","mean")).reset_index()
meta=roi[["Subject","JobID","DX","DataSource"]].drop_duplicates()
out=comp.merge(meta,on=["Subject","JobID"],how="left")
out.to_csv(args.out_csv,index=False); print("WROTE",args.out_csv)
