
import argparse, pandas as pd, numpy as np
from pathlib import Path
def auc_mw(y,s):
    y=np.asarray(y).astype(int); s=np.asarray(s).astype(float)
    pos=s[y==1]; neg=s[y==0]
    if len(pos)==0 or len(neg)==0: return np.nan
    order=np.argsort(s); ranks=np.empty_like(order); ranks[order]=np.arange(1,len(s)+1)
    r_pos=ranks[y==1].sum(); U=r_pos - len(pos)*(len(pos)+1)/2.0
    return float(U/(len(pos)*len(neg)))
def boot_ci(y,s,n=1000,seed=13):
    rng=np.random.default_rng(seed); vals=[]
    for _ in range(n):
        idx=rng.integers(0,len(y),len(y)); vals.append(auc_mw(y[idx],s[idx]))
    a=np.array(vals); lo,hi=np.nanpercentile(a,[2.5,97.5]); return float(np.nanmean(a)),float(lo),float(hi)
def youden(y,s):
    thr=np.unique(s[np.isfinite(s)]); best=-1; bestT=np.nan; stats=None
    for t in thr:
        pred=(s>=t).astype(int); tp=((pred==1)&(y==1)).sum(); tn=((pred==0)&(y==0)).sum(); fp=((pred==1)&(y==0)).sum(); fn=((pred==0)&(y==1)).sum()
        sens=tp/(tp+fn) if (tp+fn)>0 else np.nan; spec=tn/(tn+fp) if (tn+fp)>0 else np.nan
        if np.isfinite(sens) and np.isfinite(spec) and sens+spec-1>best:
            best=sens+spec-1; bestT=t; ppv=tp/(tp+fp) if (tp+fp)>0 else np.nan; npv=tn/(tn+fn) if (tn+fn)>0 else np.nan
            stats={"sens":sens,"spec":spec,"ppv":ppv,"npv":npv,"tp":int(tp),"tn":int(tn),"fp":int(fp),"fn":int(fn)}
    return float(bestT), stats
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--composites",required=True); ap.add_argument("--jobs",required=True)
    ap.add_argument("--out_dir",required=True); ap.add_argument("--contrasts",default="AD:CN,MCI:CN,PD:CN,PPD:CN")
    ap.add_argument("--boot",type=int,default=1000)
    args=ap.parse_args()
    comp=pd.read_csv(args.composites); jobs=pd.read_csv(args.jobs)
    df=comp.merge(jobs[["Subject","JobID","DX","DataSource"]],on=["Subject","JobID"],how="left")
    out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
    auc_rows=[]; thr_rows=[]
    for c in [x.strip() for x in args.contrasts.split(",") if x.strip()]:
        tgt,ref=[t.strip().upper() for t in c.split(":")]
        sub=df[df["DX"].isin([tgt,ref])].copy(); y=(sub["DX"]==tgt).astype(int).values; s=sub["IDE_global"].values.astype(float)
        auc=auc_mw(y,s); b,lo,hi=boot_ci(y,s,n=args.boot)
        thr,stats=youden(y,s); auc_rows.append({"contrast":f"{tgt} vs {ref}","cohort":"ALL","AUC":auc,"CI95_lo":lo,"CI95_hi":hi,"n_tgt":int((y==1).sum()),"n_ref":int((y==0).sum())})
        thr_rows.append({"contrast":f"{tgt} vs {ref}","cohort":"ALL","threshold":thr,**stats})
    pd.DataFrame(auc_rows).to_csv(out/"auc_overall.csv",index=False); pd.DataFrame(thr_rows).to_csv(out/"thresholds_youden.csv",index=False)
    print("WROTE", out.as_posix())
if __name__=="__main__": main()
