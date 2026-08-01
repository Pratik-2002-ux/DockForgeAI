import pandas as pd,pytest
from src.features.descriptor_features import descriptor_frame,descriptor_quality
def data(): return pd.DataFrame([{"canonical_isomeric_smiles":"CCO","molecular_weight":46,"exact_molecular_weight":46,"logp":0,"tpsa":20,"hbd":1,"hba":1,"rotatable_bonds":0,"ring_count":0,"aromatic_ring_count":0,"formal_charge":0,"heavy_atom_count":3,"fraction_csp3":1,"qed":.4,"pains_alert_count":0,"brenk_alert_count":0}])
def test_descriptor_features_finite(): assert descriptor_frame(data()).notna().all().all()
def test_descriptor_quality_reports_constants(): assert descriptor_quality(descriptor_frame(data()))["constant_columns"]
