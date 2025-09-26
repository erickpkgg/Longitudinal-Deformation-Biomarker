
import argparse, pandas as pd
from pathlib import Path
import ants
def register(fixed,moving,out_prefix):
    f=ants.image_read(str(fixed)); m=ants.image_read(str(moving))
    reg=ants.registration(f,m,type_of_transform="SyN")
    ants.image_write(reg["warpedmovout"], out_prefix+"Warped.nii.gz")
    if "fwdtransforms" in reg:
        for i,t in enumerate(reg["fwdtransforms"]):
            ants.write_transform(t, out_prefix+f"{i+1}Warp.h5")
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--jobs",required=True); ap.add_argument("--out_root",required=True)
    ap.add_argument("--skip_if_done",action="store_true")
    args=ap.parse_args()
    jobs=pd.read_csv(args.jobs)
    out_root=Path(args.out_root); out_root.mkdir(parents=True,exist_ok=True)
    rows=[]
    for _,r in jobs.iterrows():
        jid=r["JobID"]; d=out_root/jid; d.mkdir(parents=True,exist_ok=True)
        out_prefix=str((d/"ANTs_").as_posix())
        warped=d/"ANTs_Warped.nii.gz"
        if args.skip_if_done and warped.exists():
            rows.append((jid,"SKIP")); continue
        try:
            register(r["A_PreprocPath"], r["B_PreprocPath"], out_prefix); rows.append((jid,"OK"))
        except Exception as e:
            rows.append((jid,f"ERR:{e}"))
    pd.DataFrame(rows,columns=["JobID","status"]).to_csv(out_root/"registration_results.csv",index=False)
    print("DONE",out_root)
if __name__=="__main__": main()
