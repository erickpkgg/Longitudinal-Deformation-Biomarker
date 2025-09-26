
import argparse, os, time
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument("--temp_dir",default=None); ap.add_argument("--older_than_min",type=int,default=30); ap.add_argument("--sleep",type=int,default=120)
args=ap.parse_args()
temp=Path(args.temp_dir or os.environ.get("TEMP",".")).resolve()
print("Monitoring",temp)
while True:
    now=time.time(); n=0
    for p in temp.rglob("*.nii*"):
        try:
            if now-p.stat().st_mtime>args.older_than_min*60: p.unlink(); n+=1
        except: pass
    print("Removed",n,"old NIfTI; sleeping",args.sleep,"s"); time.sleep(args.sleep)
