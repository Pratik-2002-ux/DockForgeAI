"""Run sequential AutoDock Vina screening for prepared demo ligands."""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
import time
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]

RECEPTOR = (
    PROJECT_ROOT
    / "data/processed/pdbqt/receptors/"
    / "1HSG_chains-A-B_receptor.pdbqt"
)

LIGAND_DIRECTORY = (
    PROJECT_ROOT
    / "data/processed/screening_library/pdbqt"
)

RUN_DIRECTORY = (
    PROJECT_ROOT
    / "data/processed/screening_runs"
)

OUTPUT_TABLE = (
    PROJECT_ROOT
    / "data/processed/screening_tables/"
    / "screening_scores.csv"
)

REDOCKING_CONFIG = (
    PROJECT_ROOT
    / "config/docking/1HSG_MK1_redocking.txt"
)


class ScreeningError(RuntimeError):
    """Raised when demo virtual screening cannot be completed."""


def parse_config(path: Path) -> dict[str, str]:
    """Read simple key-value entries from a Vina configuration file."""

    if not path.exists():
        raise FileNotFoundError(f"Vina configuration not found: {path}")

    values: dict[str, str] = {}

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()

        if not line or line.startswith("#") or "=" not in line:
            continue

        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()

    required = {
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
    }

    missing = required - values.keys()

    if missing:
        raise ScreeningError(
            f"Missing docking-box settings: {sorted(missing)}"
        )

    return values


def parse_vina_output(output_path: Path) -> list[float]:
    """Extract affinity scores from a Vina multi-model PDBQT."""

    if not output_path.exists():
        return []

    pattern = re.compile(
        r"REMARK VINA RESULT:\s+"
        r"(-?\d+(?:\.\d+)?)"
    )

    scores: list[float] = []

    for line in output_path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        match = pattern.search(line)

        if match:
            scores.append(float(match.group(1)))

    return scores


def write_ligand_config(
    ligand_path: Path,
    output_path: Path,
    config_path: Path,
    box: dict[str, str],
) -> None:
    """Write one Vina configuration for a prepared ligand."""

    content = f"""receptor = {RECEPTOR.relative_to(PROJECT_ROOT)}
ligand = {ligand_path.relative_to(PROJECT_ROOT)}
out = {output_path.relative_to(PROJECT_ROOT)}

center_x = {box["center_x"]}
center_y = {box["center_y"]}
center_z = {box["center_z"]}

size_x = {box["size_x"]}
size_y = {box["size_y"]}
size_z = {box["size_z"]}

exhaustiveness = 4
num_modes = 3
energy_range = 3
cpu = 4
seed = 2026
"""

    config_path.write_text(content, encoding="utf-8")


def run_screening() -> list[dict[str, object]]:
    """Run Vina sequentially for all prepared demo ligands."""

    vina = shutil.which("vina")

    if vina is None:
        raise ScreeningError(
            "AutoDock Vina was not found in the active environment."
        )

    if not RECEPTOR.exists():
        raise FileNotFoundError(
            f"Receptor PDBQT not found: {RECEPTOR}"
        )

    ligand_paths = sorted(
        (
            path
            for path in LIGAND_DIRECTORY.glob("*.pdbqt")
            if not path.stem.endswith("_out")
        ),
        key=lambda path: path.stem,
    )

    if not ligand_paths:
        raise ScreeningError(
            f"No ligand PDBQT files found in {LIGAND_DIRECTORY}"
        )

    box = parse_config(REDOCKING_CONFIG)

    RUN_DIRECTORY.mkdir(parents=True, exist_ok=True)
    OUTPUT_TABLE.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []

    print(f"Vina executable: {vina}")
    print(f"Prepared ligands: {len(ligand_paths)}")
    print()

    for ligand_path in ligand_paths:
        compound_id = ligand_path.stem
        compound_run_directory = RUN_DIRECTORY / compound_id
        compound_run_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        output_path = (
            compound_run_directory
            / f"{compound_id}_screened.pdbqt"
        )
        config_path = compound_run_directory / "config.txt"
        stdout_path = compound_run_directory / "vina_stdout.txt"
        stderr_path = compound_run_directory / "vina_stderr.txt"

        write_ligand_config(
            ligand_path=ligand_path,
            output_path=output_path,
            config_path=config_path,
            box=box,
        )

        command = [
            vina,
            "--config",
            str(config_path.relative_to(PROJECT_ROOT)),
        ]

        print(f"Docking compound {compound_id}...")

        started = time.perf_counter()

        completed = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            check=False,
            timeout=1800,
        )

        runtime_seconds = time.perf_counter() - started

        stdout_path.write_text(
            completed.stdout or "",
            encoding="utf-8",
        )

        stderr_path.write_text(
            completed.stderr or "",
            encoding="utf-8",
        )

        scores = parse_vina_output(output_path)

        if completed.returncode == 0 and scores:
            sorted_scores = sorted(scores)
            top_three = sorted_scores[:3]

            row = {
                "compound_id": compound_id,
                "best_affinity": min(scores),
                "mean_top3_affinity": round(
                    sum(top_three) / len(top_three),
                    3,
                ),
                "median_affinity": sorted_scores[
                    len(sorted_scores) // 2
                ],
                "pose_count": len(scores),
                "status": "success",
                "runtime_seconds": round(runtime_seconds, 3),
                "error_message": "",
                "output_path": str(
                    output_path.relative_to(PROJECT_ROOT)
                ),
            }

            print(
                f"  SUCCESS: best score "
                f"{row['best_affinity']} kcal/mol"
            )
        else:
            error_message = (
                completed.stderr.strip()
                or completed.stdout.strip()
                or "Vina produced no parseable poses."
            )

            row = {
                "compound_id": compound_id,
                "best_affinity": "",
                "mean_top3_affinity": "",
                "median_affinity": "",
                "pose_count": 0,
                "status": "failed",
                "runtime_seconds": round(runtime_seconds, 3),
                "error_message": error_message,
                "output_path": str(
                    output_path.relative_to(PROJECT_ROOT)
                ),
            }

            print(
                f"  FAILED: inspect {stderr_path}"
            )

        rows.append(row)

    return rows


def write_results(rows: list[dict[str, object]]) -> None:
    """Write screening results to CSV."""

    fieldnames = [
        "compound_id",
        "best_affinity",
        "mean_top3_affinity",
        "median_affinity",
        "pose_count",
        "status",
        "runtime_seconds",
        "error_message",
        "output_path",
    ]

    with OUTPUT_TABLE.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    """Run the direct demo screening workflow."""

    try:
        rows = run_screening()
        write_results(rows)
    except (
        FileNotFoundError,
        ScreeningError,
        subprocess.TimeoutExpired,
    ) as error:
        print(f"\nScreening failed: {error}")
        return 1

    successful = sum(
        row["status"] == "success"
        for row in rows
    )

    failed = len(rows) - successful

    print()
    print("Screening summary")
    print("-" * 50)
    print(f"Total ligands: {len(rows)}")
    print(f"Successful:    {successful}")
    print(f"Failed:        {failed}")
    print(f"Score table:   {OUTPUT_TABLE}")

    return 0 if successful > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
