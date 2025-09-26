
import argparse, pandas as pd
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument("--composites",required=True); ap.add_argument("--jobs",required=True); ap.add_argument("--out_dir",required=True)
args=ap.parse_args()
comp=pd.read_csv(args.composites); jobs=pd.read_csv(args.jobs)
comp.merge(jobs[["Subject","JobID","DX","SITE"]],on=["Subject","JobID"],how="left").to_csv(Path(args.out_dir)/"composites_merged.csv",index=False)
print("WROTE", Path(args.out_dir)/"composites_merged.csv")
