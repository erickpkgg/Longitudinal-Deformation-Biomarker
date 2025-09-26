
import argparse, pandas as pd
from pathlib import Path
ap=argparse.ArgumentParser()
ap.add_argument("--jobs",required=True); ap.add_argument("--reg_root",required=True); ap.add_argument("--out",required=True)
args=ap.parse_args()
jobs=pd.read_csv(args.jobs); reg=Path(args.reg_root)
def ok(jid): return (reg/jid/"ANTs_Warped.nii.gz").exists()
df=jobs.assign(has_warped=jobs["JobID"].apply(ok))
Path(args.out).mkdir(parents=True,exist_ok=True)
df.to_csv(Path(args.out)/"qc_pairs_summary.csv",index=False)
df[df["has_warped"]].to_csv(Path(args.out)/"qc_pairs_summary_passlist.csv",index=False)
df[~df["has_warped"]].to_csv(Path(args.out)/"qc_pairs_summary_rejects.csv",index=False)
jobs.merge(df[["JobID","has_warped"]],on="JobID").query("has_warped").drop(columns=["has_warped"]).to_csv(Path(args.out)/"jobs_adjacent_registration.csv",index=False)
print("WROTE passlist", df["has_warped"].sum())
