from pathlib import Path
def test_v2_directories_exist():
 for p in ["src/models_v2","app/components/model_v2","scripts/model_v2","data/model_v2","reports/model_v2","docs/model_v2"]: assert Path(p).exists()
def test_v2_app_is_standalone(): assert "app.streamlit_app" not in Path("app/components/model_v2/streamlit_app.py").read_text()
def test_v1_model_not_referenced(): assert "xai_models" not in "".join(p.read_text() for p in Path("src/models_v2").glob("*.py"))
