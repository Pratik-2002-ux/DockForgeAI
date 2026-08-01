#!/usr/bin/env python3
"""Generate DockForgeAI final PDF and HTML reports."""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path: sys.path.insert(0, str(PROJECT_ROOT))

from src.reports.final_report import generate_final_reports  # noqa: E402


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    group = result.add_mutually_exclusive_group()
    group.add_argument("--pdf-only", action="store_true"); group.add_argument("--html-only", action="store_true")
    result.add_argument("--overwrite", action="store_true"); result.add_argument("--target-id", default="1HSG")
    result.add_argument("--ligand-id", default="MK1"); result.add_argument("--output-prefix")
    result.add_argument("--open", action="store_true", dest="open_report"); result.add_argument("--strict", action="store_true")
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    log_dir = PROJECT_ROOT / "logs" / "reporting"; log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=log_dir / f"report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log", level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    try:
        manifest = generate_final_reports(PROJECT_ROOT, args.target_id, args.ligand_id, args.output_prefix, not args.html_only, not args.pdf_only, args.overwrite, args.strict)
        reports = [entry for entry in manifest["generated_files"] if entry["path"].endswith((".pdf", ".html"))]
        print(f"Generated {len(reports)} report(s) in {manifest['runtime_seconds']:.3f} s")
        for entry in reports: print(f"- {entry['path']} ({entry['size_bytes']} bytes)")
        if manifest["warnings"]: print(f"Warnings: {len(manifest['warnings'])} (see manifest)")
        if args.open_report and reports:
            subprocess.run(["open", str(PROJECT_ROOT / reports[0]["path"])], check=False, timeout=10)
        return 0
    except Exception:
        logging.exception("Final report generation failed")
        print("Report generation failed; see logs/reporting/ for the traceback.", file=sys.stderr)
        return 1


if __name__ == "__main__": raise SystemExit(main())
