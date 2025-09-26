
import argparse, pandas as pd, numpy as np, nibabel as nib
from pathlib import Path
ap=argparse.ArgumentParser()
ap.add_argument("--jobs",required=True); ap.add_argument("--out_root",required=True)
ap.add_argument("--atlas_nii",required=True); ap.add_argument("--labels_csv",required=True); ap.add_argument("--brain_mask",required=True)
ap.add_argument("--out_csv",required=True)
args=ap.parse_args()
jobs=pd.read_csv(args.jobs)
atlas=nib.load(args.atlas_nii); atlas_data=atlas.get_fdata().astype(int)
labels=pd.read_csv(args.labels_csv)
mask=nib.load(args.brain_mask).get_fdata()>0
rows=[]
for _,r in jobs.iterrows():
    jid=r["JobID"]
    # En tu flujo real, aquí leerías el mapa LOG-JAC. Este es un placeholder con la imagen deformada.
    p=Path(args.out_root)/"reg"/jid/"ANTs_Warped.nii.gz"
    if not p.exists(): continue
    img=nib.load(str(p)).get_fdata()
    for _,lr in labels.iterrows():
        lab=int(lr[0]); name=str(lr[1])
        m=(atlas_data==lab)&mask
        val=float(np.nanmean(img[m])) if m.any() else np.nan
        rows.append({"JobID":jid,"Subject":r["Subject"],"DX":r["DX"],"DataSource":r.get("DataSource",""),"ROI":name,"value":val})
pd.DataFrame(rows).to_csv(args.out_csv,index=False); print("WROTE",args.out_csv)
