#!/usr/bin/env python3
"""Run the Phase 2 redocking workflow end to end."""

from __future__ import annotations

import sys
import traceback
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.docking.input_validator import DockingInputError, validate_docking_inputs
from src.docking.pose_exporter import export_poses
from src.docking.redocking_validator import validate_redocking
from src.docking.report_generator import generate_reports
from src.docking.result_parser import parse_vina_result
from src.docking.rmsd import calculate_rmsd
from src.docking.vina_runner import VinaExecutionError, run_vina


class RedockingExecutionError(RuntimeError):
    """Raised when the overall redocking workflow fails."""


def main(argv: list[str] | None = None) -> int:
    """Execute the complete redocking workflow."""

    argv = list(argv or sys.argv[1:])
    overwrite = "--overwrite" in argv
    skip_docking = "--skip-docking" in argv
    threshold = 2.0
    for index, arg in enumerate(argv):
        if arg == "--threshold" and index + 1 < len(argv):
            threshold = float(argv[index + 1])

    try:
        validation = validate_docking_inputs()
        settings = validation["settings"]
        docking_settings = settings.get("docking", {})
        threshold = float(docking_settings.get("redocking_rmsd_threshold_angstrom", threshold))
        run_dir = PROJECT_ROOT / "data" / "processed" / "docking_runs" / "1HSG_MK1"
        pose_dir = PROJECT_ROOT / "data" / "processed" / "docking_poses" / "1HSG_MK1"
        table_dir = PROJECT_ROOT / "data" / "processed" / "docking_tables"
        result_dir = PROJECT_ROOT / "data" / "processed" / "docking_results"
        report_dir = PROJECT_ROOT / "data" / "processed" / "docking_reports"
        log_dir = PROJECT_ROOT / "logs" / "docking"
        run_dir.mkdir(parents=True, exist_ok=True)
        pose_dir.mkdir(parents=True, exist_ok=True)
        table_dir.mkdir(parents=True, exist_ok=True)
        result_dir.mkdir(parents=True, exist_ok=True)
        report_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        output_pdbqt = run_dir / "1HSG_MK1_redocked.pdbqt"
        stdout_path = run_dir / "vina_stdout.txt"
        stderr_path = run_dir / "vina_stderr.txt"
        execution_metadata_path = run_dir / "vina_execution.json"
        score_table_path = table_dir / "1HSG_MK1_vina_scores.csv"
        parser_metadata_path = pose_dir / "parser_metadata.json"
        pose_sdf = pose_dir / "1HSG_MK1_redocked_poses.sdf"
        export_report_path = pose_dir / "pose_export_report.json"
        analysis_table_path = table_dir / "1HSG_MK1_redocking_analysis.csv"
        rmsd_metadata_path = pose_dir / "rmsd_metadata.json"
        summary_path = result_dir / "1HSG_MK1_redocking_summary.json"
        report_md_path = report_dir / "1HSG_MK1_redocking_report.md"
        report_json_path = report_dir / "1HSG_MK1_redocking_report.json"
        preparation_report_path = PROJECT_ROOT / "data" / "processed" / "preparation_reports" / "1HSG_receptor_preparation.json"
        box_report_path = PROJECT_ROOT / "data" / "processed" / "docking_boxes" / "1HSG_MK1_box.json"
        config_path = PROJECT_ROOT / "config" / "docking" / "1HSG_MK1_redocking.txt"

        if not skip_docking:
            print("Running Vina docking...")
            execution_metadata = run_vina(
                config_path=config_path,
                output_path=output_pdbqt,
                log_dir=run_dir,
                overwrite=overwrite,
            )
            execution_metadata_path.write_text(__import__("json").dumps(execution_metadata, indent=2), encoding="utf-8")
        else:
            if not output_pdbqt.exists():
                raise RedockingExecutionError("Docking output missing and --skip-docking was specified")
            execution_metadata = {
                "command": f"vina --config {config_path}",
                "executable": "vina",
                "vina_version": "skipped",
                "runtime_seconds": 0.0,
                "output_path": str(output_pdbqt),
            }

        print("Parsing docking output...")
        parse_vina_result(
            vina_output_path=output_pdbqt,
            output_dir=pose_dir,
            score_table_path=score_table_path,
            metadata_path=parser_metadata_path,
        )

        print("Exporting docked poses...")
        export_poses(
            pose_pdbqts=[pose_dir / f"pose_{index:03d}.pdbqt" for index in range(1, 11)],
            output_sdf=pose_sdf,
            reference_sdf=validation["reference_sdf"],
            report_path=export_report_path,
        )

        print("Calculating RMSD...")
        calculate_rmsd(
            reference_sdf=validation["reference_sdf"],
            pose_sdf=pose_sdf,
            score_table_path=score_table_path,
            output_table_path=analysis_table_path,
            metadata_path=rmsd_metadata_path,
            threshold=threshold,
        )

        print("Evaluating redocking success...")
        validate_redocking(analysis_table_path, summary_path, threshold)

        print("Generating reports...")
        generate_reports(
            analysis_csv=analysis_table_path,
            summary_path=summary_path,
            execution_metadata=execution_metadata,
            report_md_path=report_md_path,
            report_json_path=report_json_path,
            preparation_report_path=preparation_report_path,
            box_report_path=box_report_path,
            config_path=config_path,
            threshold=threshold,
        )

        print("Redocking completed successfully.")
        return 0
    except DockingInputError as error:
        print(f"Input validation failed: {error}")
        return 1
    except VinaExecutionError as error:
        print(f"Vina execution failed: {error}")
        return 1
    except Exception as error:
        log_path = PROJECT_ROOT / "logs" / "docking" / "redocking_traceback.log"
        log_path.write_text(traceback.format_exc(), encoding="utf-8")
        print(f"Redocking workflow failed: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
