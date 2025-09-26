
import argparse, pandas as pd, numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
ap=argparse.ArgumentParser(); ap.add_argument("--composites",required=True); ap.add_argument("--jobs",required=True); ap.add_argument("--out_dir",required=True)
args=ap.parse_args(); comp=pd.read_csv(args.composites); jobs=pd.read_csv(args.jobs)
df=comp.merge(jobs[["Subject","JobID","DX","SITE"]],on=["Subject","JobID"],how="left"); out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
order=["CN","MCI","AD","PD","PPD"]; data=[df[df["DX"]==dx]["IDE_global"].dropna().values for dx in order if (df["DX"]==dx).any()]; labels=[dx for dx in order if (df["DX"]==dx).any()]
plt.figure(figsize=(6,4)); plt.boxplot(data,labels=labels,showfliers=False); plt.ylabel("IDE_global"); plt.title("IDE_global por diagnóstico"); plt.tight_layout(); plt.savefig(out/"fig_box_IDE_by_DX.png",dpi=300); plt.close()
cn=df[df["DX"]=="CN"].copy(); cn["SITE_next"]=cn.groupby("Subject")["SITE"].shift(-1); cn["SITE_mismatch"]=np.where(cn["SITE"]==cn["SITE_next"],"same-site","mismatch")
plt.figure(figsize=(6,4)); plt.boxplot([cn[cn['SITE_mismatch']=='same-site']['IDE_global'],cn[cn['SITE_mismatch']=='mismatch']['IDE_global']],labels=["same-site","mismatch"],showfliers=False)
plt.ylabel("IDE_global"); plt.title("CN: IDE_global por SITE_mismatch"); plt.tight_layout(); plt.savefig(out/"fig_box_CN_IDE_by_site_mismatch.png",dpi=300); plt.close()
print("WROTE",out.as_posix())
