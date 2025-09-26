
import argparse, numpy as np, pandas as pd
from pathlib import Path
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
    r=ranks[y==1].sum(); U=r - len(pos)*(len(pos)+1)/2.0
    return float(U/(len(pos)*len(neg)))
def boot_ci(y,s,n=1000,seed=13):
    rng=np.random.default_rng(seed); vals=[]
    for _ in range(n): idx=rng.integers(0,len(y),len(y)); vals.append(auc_mw(y[idx],s[idx]))
    a=np.array(vals); lo,hi=np.nanpercentile(a,[2.5,97.5]); return float(np.nanmean(a)),float(lo),float(hi)
def youden(y,s):
    thr=np.unique(s[np.isfinite(s)]); best=-1; BT=np.nan; M=None
    for t in thr:
        pred=(s>=t).astype(int); tp=((pred==1)&(y==1)).sum(); tn=((pred==0)&(y==0)).sum(); fp=((pred==1)&(y==0)).sum(); fn=((pred==0)&(y==1)).sum()
        sens=tp/(tp+fn) if (tp+fn)>0 else np.nan; spec=tn/(tn+fp) if (tn+fp)>0 else np.nan
        if np.isfinite(sens) and np.isfinite(spec) and sens+spec-1>best:
            best=sens+spec-1; BT=t; ppv=tp/(tp+fp) if (tp+fp)>0 else np.nan; npv=tn/(tn+fn) if (tn+fn)>0 else np.nan
            M={"sens":sens,"spec":spec,"ppv":ppv,"npv":npv,"tp":int(tp),"tn":int(tn),"fp":int(fp),"fn":int(fn)}
    return float(BT), M
def site_centering(df, site="SITE", score="IDE_global", dx="DX", min_cn=30):
    df=df.copy()
    if site not in df.columns: df["IDE_siteZ"]=df[score]; return df
    cn_all=df[df[dx]=="CN"][score].values
    mu_all=float(np.nanmean(cn_all)) if len(cn_all)>0 else 0.0
    sd_all=float(np.nanstd(cn_all,ddof=1)) if len(cn_all)>1 else 1.0
    if sd_all<=1e-6: sd_all=1.0
    stats={}
    for s, g in df.groupby(site):
        cn=g[g[dx]=="CN"][score].values
        if len(cn)>=min_cn and np.isfinite(cn).sum()>=10:
            mu=float(np.nanmean(cn)); sd=float(np.nanstd(cn,ddof=1)); sd=sd if sd>1e-6 else 1.0
        else: mu,sd=mu_all,sd_all
        stats[s]=(mu,sd)
    df["IDE_siteZ"]=df.apply(lambda r: (r[score]-stats.get(r.get(site),(mu_all,sd_all))[0])/stats.get(r.get(site),(mu_all,sd_all))[1],axis=1)
    return df
def aggregate_subject(df, score="IDE_siteZ", dx="DX"):
    return df.groupby("Subject").agg(score=(score,"mean"), DX=(dx, lambda x: x.mode().iat[0] if len(x.mode())>0 else x.iloc[0])).reset_index()
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--composites",required=True); ap.add_argument("--jobs",required=True); ap.add_argument("--out_dir",required=True)
    ap.add_argument("--contrasts",default="AD:CN,MCI:CN,PD:CN,PPD:CN"); ap.add_argument("--min_dt",type=float,default=0.5); ap.add_argument("--boot",type=int,default=1000)
    args=ap.parse_args()
    comp=pd.read_csv(args.composites); jobs=pd.read_csv(args.jobs, parse_dates=["A_AcqDate","B_AcqDate"], dayfirst=True)
    if "DeltaYears" not in jobs.columns: jobs["DeltaYears"]=(jobs["B_AcqDate"]-jobs["A_AcqDate"]).dt.days/365.25
    meta=jobs[["Subject","JobID","DeltaYears","DX","DataSource","SITE"]].copy()
    meta["DX"]=(meta["DX"].astype(str).str.upper().str.strip().str.replace(r"\.$","",regex=True).str.replace(r"\s+","",regex=True).replace({"PRODROMAL":"PPD"}))
    df=comp.merge(meta,on=["Subject","JobID"],how="left"); df=df[df["DeltaYears"]>=args.min_dt].copy()
    df=site_centering(df, site="SITE", score="IDE_global", dx="DX", min_cn=30)
    out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
    # pair-level
    auc_rows=[]; thr_rows=[]
    for con in [c for c in args.contrasts.split(",") if c.strip()]:
        tgt,ref=[x.strip().upper() for x in con.split(":")]
        sub=df[df["DX"].isin([tgt,ref])]; y=(sub["DX"]==tgt).astype(int).values; s=sub["IDE_siteZ"].values.astype(float)
        auc=auc_mw(y,s); b,lo,hi=boot_ci(y,s,n=args.boot); thr,met=youden(y,s)
        auc_rows.append({"contrast":f"{tgt} vs {ref}","cohort":"ALL","AUC":auc,"CI95_lo":lo,"CI95_hi":hi,"n_tgt":int((y==1).sum()),"n_ref":int((y==0).sum())})
        thr_rows.append({"contrast":f"{tgt} vs {ref}","cohort":"ALL","threshold":thr,**met})
    pd.DataFrame(auc_rows).to_csv(out/"auc_pairs.csv",index=False); pd.DataFrame(thr_rows).to_csv(out/"youdens_pairs.csv",index=False)
    # subject-level
    subj=aggregate_subject(df, score="IDE_siteZ", dx="DX"); auc_rows=[]; thr_rows=[]
    for con in [c for c in args.contrasts.split(",") if c.strip()]:
        tgt,ref=[x.strip().upper() for x in con.split(":")]
        sub=subj[subj["DX"].isin([tgt,ref])]; y=(sub["DX"]==tgt).astype(int).values; s=sub["score"].values.astype(float)
        auc=auc_mw(y,s); b,lo,hi=boot_ci(y,s,n=args.boot); thr,met=youden(y,s)
        auc_rows.append({"contrast":f"{tgt} vs {ref}","cohort":"ALL","AUC":auc,"CI95_lo":lo,"CI95_hi":hi,"n_tgt":int((y==1).sum()),"n_ref":int((y==0).sum())})
        thr_rows.append({"contrast":f"{tgt} vs {ref}","cohort":"ALL","threshold":thr,**met})
    pd.DataFrame(auc_rows).to_csv(out/"auc_subjects.csv",index=False); pd.DataFrame(thr_rows).to_csv(out/"youdens_subjects.csv",index=False)
    print("DONE", out.as_posix())
if __name__=="__main__": main()
