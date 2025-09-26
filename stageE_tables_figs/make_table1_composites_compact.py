
import argparse, pandas as pd
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument("--composites",required=True); ap.add_argument("--out_tex",required=True)
args=ap.parse_args(); df=pd.read_csv(args.composites)
lines=["\\begin{tabular}{l r r}","\\toprule","DX & n & IDE$_{global}$ (mean) \\\\","\\midrule"]
for dx,g in df.groupby("DX"): lines.append(f"{dx} & {len(g)} & {g['IDE_global'].mean():.3f} \\\\")
lines+=["\\bottomrule","\\end{tabular}"]; Path(args.out_tex).write_text("\\n".join(lines),encoding="utf-8"); print("WROTE",args.out_tex)
