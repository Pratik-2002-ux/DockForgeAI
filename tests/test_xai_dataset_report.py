import json
from pathlib import Path
from src.datasets.dataset_report import generate_dataset_reports, write_table
from src.datasets.decoy_matcher import FEATURES

def row(label): return {"binary_label":label,"source_database":"x","scaffold_id":"s","dataset_split":"train","compound_id":"c","pchembl_value":6,"standard_activity_value_nm":1,**{key:1 for key in FEATURES}}
def test_class_balance_reporting(tmp_path: Path) -> None:
    records=[row("active"),row("inactive_or_decoy")]; generate_dataset_reports(records,tmp_path,[records[0]],[records[1]])
    assert json.loads((tmp_path/"dataset_summary.json").read_text())["dataset_size"]==2
def test_empty_exclusion_table_written(tmp_path: Path) -> None:
    output=write_table(tmp_path/"empty.csv",[]); assert output.is_file()
