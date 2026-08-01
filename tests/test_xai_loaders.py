from app.components.xai_loaders import load_xai_csv,load_xai_json
def test_missing_xai_json_safe(): assert load_xai_json("missing.json")=={}
def test_missing_xai_csv_safe(): assert load_xai_csv("missing.csv").empty
