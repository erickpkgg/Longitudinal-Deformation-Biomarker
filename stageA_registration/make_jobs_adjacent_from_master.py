
import argparse, pandas as pd
from pathlib import Path
def find(cols, cands):
    low={c.lower():c for c in cols}
    for x in cands:
        if x.lower() in low: return low[x.lower()]
    for c in cols:
        for x in cands:
            if x.lower() in c.lower(): return c
    return None
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--master",required=True)
    ap.add_argument("--out_csv",required=True)
    args=ap.parse_args()
    p=Path(args.master)
    df=pd.read_excel(p, engine="openpyxl") if p.suffix.lower() in (".xlsx",".xls") else pd.read_csv(p)
    subj=find(df.columns,["Subject","RID","ID"]); acq=find(df.columns,["Acq Date","AcqDate"])
    age=find(df.columns,["Age"]); sex=find(df.columns,["Gender","Sex"]); site=find(df.columns,["SITE","Site"])
    cohort=find(df.columns,["DataSource","Cohort","Dataset"]); dx=find(df.columns,["Condition","DX","Diagnosis"])
    ppath=find(df.columns,["Preprocessed Path","Path","ImagePath"])
    assert all([subj,acq,age,sex,site,cohort,dx,ppath]), "Faltan columnas en master"
    df=df.rename(columns={subj:"Subject",acq:"AcqDate",age:"Age",sex:"Gender",site:"SITE",cohort:"DataSource",dx:"DX",ppath:"PreprocPath"})
    df["DX"]=(df["DX"].astype(str).str.upper().str.strip().str.replace(r"\.$","",regex=True).str.replace(r"\s+","",regex=True).replace({"PRODROMAL":"PPD","PPD ":"PPD"}))
    df["AcqDate"]=pd.to_datetime(df["AcqDate"],dayfirst=True,errors="coerce")
    df=df.sort_values(["Subject","AcqDate"])
    rows=[]
    for sid,g in df.groupby("Subject"):
        g=g.dropna(subset=["AcqDate"]).sort_values("AcqDate")
        for i in range(len(g)-1):
            A=g.iloc[i]; B=g.iloc[i+1]
            rows.append({"JobID":f"{sid}__{i+1:02d}","Subject":sid,
                         "A_AcqDate":A["AcqDate"],"B_AcqDate":B["AcqDate"],
                         "DeltaYears":(B["AcqDate"]-A["AcqDate"]).days/365.25,
                         "AgeA":A["Age"],"AgeB":B["Age"],"Gender":A["Gender"],
                         "SITE":A["SITE"],"DataSource":A["DataSource"],"DX":A["DX"],
                         "A_PreprocPath":A["PreprocPath"],"B_PreprocPath":B["PreprocPath"]})
    jobs=pd.DataFrame(rows).dropna(subset=["A_PreprocPath","B_PreprocPath"])
    jobs.to_csv(args.out_csv,index=False); print("WROTE",args.out_csv,len(jobs))
if __name__=="__main__": main()
