
import argparse, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
def find_col(cols, cands):
    low={c.lower():c for c in cols}
    for x in cands:
        if x.lower() in low: return low[x.lower()]
    for c in cols:
        for x in cands:
            if x.lower() in c.lower(): return c
    return None
def auc_mw(y,s):
    y=np.asarray(y).astype(int); s=np.asarray(s).astype(float)
    pos=s[y==1]; neg=s[y==0]
    if len(pos)==0 or len(neg)==0: return np.nan
    order=np.argsort(s); ranks=np.empty_like(order); ranks[order]=np.arange(1,len(s)+1)
    rpos=ranks[y==1].sum(); U=rpos-len(pos)*(len(pos)+1)/2.0
    return float(U/(len(pos)*len(neg)))
def cindex_ordinal_fast(labels,scores,levels):
    rank={lv:i for i,lv in enumerate(levels)}
    y=np.array([rank.get(l,np.nan) for l in labels],float); s=np.asarray(scores,float)
    ok=np.isfinite(y)&np.isfinite(s); y=y[ok].astype(int); s=s[ok]
    n0=(y==0).sum(); n1=(y==1).sum(); n2=(y==2).sum(); total=n0*n1+n0*n2+n1*n2
    if total==0: return np.nan
    order=np.argsort(s,kind='mergesort'); ys=y[order]; ss=s[order]
    conc=0; c0=c1=c2=0; i=0; n=len(ys)
    while i<n:
        j=i
        while j<n and ss[j]==ss[i]: j+=1
        g=ys[i:j]; conc += int((g==1).sum())*c0; conc += int((g==2).sum())*(c0+c1)
        c0 += int((g==0).sum()); c1 += int((g==1).sum()); c2 += int((g==2).sum()); i=j
    return conc/total
def bootstrap_metric(func, y, s, n=300, seed=17):
    rng=np.random.default_rng(seed); vals=[]
    for _ in range(n):
        idx=rng.integers(0,len(y),len(y)); vals.append(func(y[idx],s[idx]))
    a=np.array(vals,float); lo,hi=np.nanpercentile(a,[2.5,97.5]); return float(np.nanmean(a)),float(lo),float(hi)
ap=argparse.ArgumentParser()
ap.add_argument("--composites",required=True); ap.add_argument("--jobs",required=True); ap.add_argument("--out_dir",required=True)
ap.add_argument("--axes",required=True,help='Ej.: "ADNI:CN<MCI<AD,PPMI:CN<PPD<PD"')
ap.add_argument("--boot",type=int,default=300); ap.add_argument("--cindex_boot",type=int,default=150); ap.add_argument("--min_per_class",type=int,default=30)
args=ap.parse_args()
comp=pd.read_csv(args.composites); jobs=pd.read_csv(args.jobs)
dx_col=find_col(jobs.columns,["DX","Diagnosis","Condition"]); coh_col=find_col(jobs.columns,["DataSource","Cohort","Dataset"])
jobs["_DX"]=jobs[dx_col].astype(str).str.upper().str.strip().str.replace(r"\.$","",regex=True).str.replace(r"\s+","",regex=True).replace({"PRODROMAL":"PPD"})
jobs["_Cohort"]=jobs[coh_col].astype(str)
df=comp.merge(jobs[["Subject","JobID","_DX","_Cohort"]],on=["Subject","JobID"],how="left").rename(columns={"_DX":"DX","_Cohort":"Cohort"})
axes={}
for chunk in [s.strip() for s in args.axes.split(",") if s.strip()]:
    name, order = chunk.split(":",1)
    levels=[x.strip().upper() for x in order.split("<") if x.strip()]
    levels=["PPD" if x in ("PPD","PRODROMAL") else x for x in levels]
    axes[name.strip()] = levels
trend=[]; cidx=[]; ovr=[]
for coh,levels in axes.items():
    sub=df[(df["Cohort"]==coh) & (df["DX"].isin(levels))].copy()
    if any(sub["DX"].value_counts().reindex(levels,fill_value=0)<args.min_per_class):
        trend.append({"cohort":coh,"levels":"<".join(levels),"rho_s":np.nan,"p_s":np.nan,"n":len(sub)})
        cidx.append({"cohort":coh,"levels":"<".join(levels),"c_index":np.nan,"ci_lo":np.nan,"ci_hi":np.nan,"n":len(sub)})
        continue
    rank={lv:i for i,lv in enumerate(levels)}; sub["rank"]=sub["DX"].map(rank).astype(int)
    from scipy.stats import spearmanr
    rho,p=spearmanr(sub["rank"].values, sub["IDE_global"].values, nan_policy='omit')
    trend.append({"cohort":coh,"levels":"<".join(levels),"rho_s":float(rho),"p_s":float(p),"n":len(sub)})
    ci = cindex_ordinal_fast(sub["DX"].values, sub["IDE_global"].values, levels)
    y_ord=sub["rank"].values; s=sub["IDE_global"].values
    b,lo,hi=bootstrap_metric(lambda yy,ss: cindex_ordinal_fast(sub["DX"].values, sub["IDE_global"].values, levels), y_ord, s, n=args.cindex_boot)
    cidx.append({"cohort":coh,"levels":"<".join(levels),"c_index":float(ci),"ci_lo":float(lo),"ci_hi":float(hi),"n":len(sub)})
    # OVR
    for lv in levels:
        y=(sub["DX"].values==lv).astype(int); s=sub["IDE_global"].values.astype(float)
        a=auc_mw(y,s); b,lo,hi=bootstrap_metric(auc_mw,y,s,n=args.boot)
        ovr.append({"cohort":coh,"levels":"<".join(levels),"class":lv,"AUC":float(a),"ci_lo":float(lo),"ci_hi":float(hi)})
    macro=np.mean([r["AUC"] for r in ovr[-len(levels):]]); ovr.append({"cohort":coh,"levels":"<".join(levels),"class":"MACRO","AUC":float(macro),"ci_lo":np.nan,"ci_hi":np.nan})
out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
pd.DataFrame(trend).to_csv(out/"ordinal_trend.csv",index=False)
pd.DataFrame(cidx).to_csv(out/"ordinal_cindex.csv",index=False)
pd.DataFrame(ovr).to_csv(out/"ovr_auc.csv",index=False)
print("WROTE", out.as_posix())
