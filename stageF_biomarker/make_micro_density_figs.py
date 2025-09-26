
import argparse, numpy as np, pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
def find_col(cols, cands):
    low={c.lower():c for c in cols}
    for x in cands:
        if x.lower() in low: return low[x.lower()]
    for c in cols:
        for x in cands:
            if x.lower() in c.lower(): return c
    return None
def site_centering(df, site="SITE", score="IDE_global", dx="DX", min_cn=30):
    df=df.copy(); 
    if site not in df.columns: df["IDE_siteZ"]=df[score]; return df
    cn=df[df[dx]=="CN"][score].values
    mu_all=float(np.nanmean(cn)) if len(cn)>0 else 0.0; sd_all=float(np.nanstd(cn,ddof=1)) if len(cn)>1 else 1.0
    if sd_all<=1e-6: sd_all=1.0
    stats={}
    for s, g in df.groupby(site):
        cn=g[g[dx]=="CN"][score].values
        if len(cn)>=min_cn and np.isfinite(cn).sum()>=10:
            mu=float(np.nanmean(cn)); sd=float(np.nanstd(cn,ddof=1)); sd=sd if sd>1e-6 else 1.0
        else: mu,sd=mu_all,sd_all
        stats[s]=(mu,sd)
    df["IDE_siteZ"]=df.apply(lambda r: (r[score]-stats.get(r.get(site),(mu_all,sd_all))[0])/stats.get(r.get(site),(mu_all,sd_all))[1],axis=1); return df
def aggregate_subject(df, score="IDE_siteZ", dx="DX"):
    return df.groupby("Subject").agg(score=(score,"mean"), DX=(dx, lambda x: x.mode().iat[0] if len(x.mode())>0 else x.iloc[0])).reset_index()
def plot_density(pos,neg,thr,out_png,title):
    plt.figure(figsize=(5,4))
    plt.hist(neg,bins=50,density=True,alpha=0.5,label="CN",histtype='stepfilled')
    plt.hist(pos,bins=50,density=True,alpha=0.5,label="Caso",histtype='stepfilled')
    if thr is not None and np.isfinite(thr): plt.axvline(thr,ls="--",lw=1.5,label=f"Umbral={thr:.3f}")
    plt.xlabel("IDE_global (z por sitio)"); plt.ylabel("Densidad"); plt.title(title); plt.legend(); plt.tight_layout(); plt.savefig(out_png,dpi=300); plt.close()
ap=argparse.ArgumentParser()
ap.add_argument("--composites",required=True); ap.add_argument("--jobs",required=True); ap.add_argument("--youdens",required=True); ap.add_argument("--out_dir",required=True)
args=ap.parse_args()
comp=pd.read_csv(args.composites); jobs=pd.read_csv(args.jobs, parse_dates=["A_AcqDate","B_AcqDate"], dayfirst=True)
if "DeltaYears" not in jobs.columns: jobs["DeltaYears"]=(jobs["B_AcqDate"]-jobs["A_AcqDate"]).dt.days/365.25
meta=jobs[["Subject","JobID","DeltaYears","DX","SITE"]].copy()
meta["DX"]=(meta["DX"].astype(str).str.upper().str.strip().str.replace(r"\.$","",regex=True).str.replace(r"\s+","",regex=True).replace({"PRODROMAL":"PPD"}))
df=comp.merge(meta,on=["Subject","JobID"],how="left"); df=df[df["DeltaYears"]>=0.5].copy(); df=site_centering(df,"SITE","IDE_global","DX",30); subj=aggregate_subject(df,"IDE_siteZ","DX")
y=pd.read_csv(args.youdens); y=y[y["cohort"]=="ALL"].set_index("contrast")
for con,label in [("AD vs CN","AD"),("MCI vs CN","MCI"),("PD vs CN","PD"),("PPD vs CN","PPD")]:
    if con not in y.index: continue
    thr=float(y.loc[con,"threshold"]) if "threshold" in y.columns else np.nan
    sub=subj[subj["DX"].isin([label,"CN"])]; pos=sub.loc[sub["DX"]==label,"score"].values; neg=sub.loc[sub["DX"]=="CN","score"].values
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    plot_density(pos,neg,thr, Path(args.out_dir)/f"density_{label}_vs_CN_subjects.png", f"Densidad — {label} vs CN (sujetos)")
print("WROTE", args.out_dir)
