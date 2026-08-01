#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.xai_reports.final_report import generate_xai_report
if __name__=="__main__":print(generate_xai_report(ROOT))
