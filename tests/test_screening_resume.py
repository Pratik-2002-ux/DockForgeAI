from pathlib import Path

from src.screening.batch_runner import run_batch


def test_run_batch_respects_existing_outputs(tmp_path: Path) -> None:
    config_template = tmp_path / "config.txt"
    config_template.write_text("receptor = x\nligand = {ligand}\n", encoding="utf-8")
    output_dir = tmp_path / "runs"
    output_dir.mkdir(parents=True, exist_ok=True)
    existing_output = output_dir / "A" / "output.pdbqt"
    existing_output.parent.mkdir(parents=True, exist_ok=True)
    existing_output.write_text("dummy", encoding="utf-8")

    result = run_batch([{"compound_id": "A"}], config_template=config_template, output_dir=output_dir, timeout_seconds=10, resume=True)

    assert result["results"][0]["status"] == "skipped"
