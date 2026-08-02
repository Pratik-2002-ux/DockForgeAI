import json
from pathlib import Path

from scripts.model_v2.audit_v2_readiness import _available


ROOT = Path(__file__).resolve().parents[2]


def test_locked_external_is_not_inside_v2_training_data():
    source = (ROOT / "scripts/model_v2/audit_v2_readiness.py").read_text()
    assert '"training_use_forbidden": True' in source


def test_required_backend_detection():
    status = _available("scikit-learn")
    assert status["available"] is True
    assert status["version"]


def test_generated_readiness_gate_blocks_training():
    path = ROOT / "data/model_v2/readiness/readiness_report.json"
    if not path.exists():
        return
    report = json.loads(path.read_text())
    assert report["training_allowed"] is False
    assert report["locked_external"]["training_use_forbidden"] is True
