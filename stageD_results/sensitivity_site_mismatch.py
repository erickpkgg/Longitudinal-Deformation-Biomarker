
import argparse, pandas as pd, numpy as np
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument("--composites",required=True); ap.add_argument("--jobs",required=True); ap.add_argument("--out_dir",required=True)
args=ap.parse_args()
comp=pd.read_csv(args.composites); jobs=pd.read_csv(args.jobs)
df=comp.merge(jobs[["Subject","JobID","DX","SITE"]],on=["Subject","JobID"],how="left")
cn=df[df["DX"]=="CN"].copy(); cn["SITE_next"]=cn.groupby("Subject")["SITE"].shift(-1); cn["SITE_mismatch"]=(cn["SITE"]!=cn["SITE_next"]).astype(int)
cn.to_csv(Path(args.out_dir)/"cn_site_mismatch.csv",index=False); print("WROTE",Path(args.out_dir)/"cn_site_mismatch.csv")
