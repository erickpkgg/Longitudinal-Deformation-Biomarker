
import argparse, pandas as pd
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument("--roi_dir",required=True); ap.add_argument("--topk",type=int,default=5); ap.add_argument("--out_dir",required=True)
args=ap.parse_args(); out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
for con in ["AD_vs_CN","MCI_vs_CN","PD_vs_CN","PPD_vs_CN"]:
    p=Path(args.roi_dir)/f"roi_effects_{con}.csv"
    if not p.exists(): continue
    df=pd.read_csv(p).sort_values("p").head(args.topk)
    lines=["\\begin{tabular}{l r r}","\\toprule","ROI & Δmean & p \\\\","\\midrule"]
    for _,r in df.iterrows(): lines.append(f"{r['ROI']} & {r['delta']:.3f} & {r['p']:.2e} \\\\")
    lines+=["\\bottomrule","\\end{tabular}"]; (out/f"table_top{args.topk}_{con}.tex").write_text("\\n".join(lines),encoding="utf-8")
print("WROTE",out.as_posix())
