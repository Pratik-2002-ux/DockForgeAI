#!/usr/bin/env python3
import argparse,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=argparse.ArgumentParser();p.add_argument("--seeds",nargs="+",default=["2026","2027","2028","2029","2030"]);p.add_argument("--bootstrap-iterations",default="1000");p.add_argument("--overwrite",action="store_true");a=p.parse_args();py=sys.executable
 commands=[[py,str(ROOT/"scripts/model_v2/build_ki_regression_dataset.py"),"--overwrite"],[py,str(ROOT/"scripts/model_v2/build_ki_regression_features.py")],[py,str(ROOT/"scripts/model_v2/train_ki_regression.py"),"--seeds",*a.seeds,"--bootstrap-iterations",str(a.bootstrap_iterations),"--overwrite"],[py,str(ROOT/"scripts/model_v2/generate_ki_regression_report.py"),"--overwrite"]]
 for c in commands:subprocess.run(c,cwd=ROOT,check=True)
 return 0
if __name__=="__main__":raise SystemExit(main())
