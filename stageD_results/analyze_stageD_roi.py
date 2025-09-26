
import argparse, pandas as pd, numpy as np
from pathlib import Path
from scipy.stats import ttest_ind
def welch(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return ttest_ind(a,b,equal_var=False).pvalue, (np.nanmean(a)-np.nanmean(b))
ap=argparse.ArgumentParser(); ap.add_argument("--roi_csv",required=True); ap.add_argument("--out_dir",required=True)
args=ap.parse_args(); df=pd.read_csv(args.roi_csv); out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
for tgt,ref in [("AD","CN"),("MCI","CN"),("PD","CN"),("PPD","CN")]:
    sub=df[df["DX"].isin([tgt,ref])]; rows=[]
    for roi,g in sub.groupby("ROI"):
        a=g[g["DX"]==tgt]["value"]; b=g[g["DX"]==ref]["value"]; p,d=welch(a,b); rows.append({"contrast":f"{tgt}_vs_{ref}","ROI":roi,"p":p,"delta":d})
    pd.DataFrame(rows).sort_values("p").to_csv(out/f"roi_effects_{tgt}_vs_{ref}.csv",index=False)
print("WROTE",out)
