#!/usr/bin/env python3
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.publication_claims import claims_for_phases
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--completed-phases",nargs="+",type=int,default=[1]);p.add_argument("--output",default="reports/publication_claims.json");a=p.parse_args();path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(claims_for_phases(set(a.completed_phases)),indent=2)+"\n")
