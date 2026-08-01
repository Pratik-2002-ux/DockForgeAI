#!/usr/bin/env python
"""Regenerate reports from existing validation artifacts."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.validation.validation_report import build_scorecard,generate_figures,generate_reports
if __name__=="__main__": build_scorecard(); generate_figures(); generate_reports()
