"""Activity harmonization and conflict-aware multi-source curation."""
from __future__ import annotations
import numpy as np,pandas as pd
from src.datasets.structure_standardizer import standardize_structure,StructureStandardizationError
SUPPORTED={"KI","KD","IC50","EC50"}
UNIT_TO_NM={"NM":1.,"UM":1000.,"ΜM":1000.,"MM":1_000_000.,"PM":.001}
def to_nm(value,unit):
 factor=UNIT_TO_NM.get(str(unit).strip().upper())
 if factor is None: raise ValueError(f"Unsupported activity unit: {unit}")
 result=float(value)*factor
 if not np.isfinite(result) or result<=0: raise ValueError("Invalid activity value")
 return result
def p_activity(nm): return float(9-np.log10(float(nm)))
def assay_quality_score(confidence,direct_target=True,standardized=True): return float(np.clip(.6*float(confidence)/10+.25*direct_target+.15*standardized,0,1))
def curate_records(frame:pd.DataFrame,active_p=6.,inactive_p=5.,conflict_delta=1.):
 rows=[]; excluded=[]
 for i,r in frame.iterrows():
  try:
   kind=str(r.measurement_type).upper()
   if kind not in SUPPORTED: raise ValueError("unsupported_measurement_type")
   structure=standardize_structure(r.smiles); nm=to_nm(r.value,r.unit); pa=p_activity(nm)
   label="active" if pa>=active_p else "inactive" if pa<inactive_p else None
   if label is None: raise ValueError("ambiguous_activity_range")
   quality=assay_quality_score(r.get("assay_confidence",5),bool(r.get("direct_target",True)),True); source_weight=float(r.get("source_weight",1.)); rows.append({**r.to_dict(),**structure,"standard_value_nm":nm,"p_activity":pa,"binary_label":label,"assay_quality_score":quality,"record_weight":quality*source_weight})
  except (ValueError,StructureStandardizationError) as exc: excluded.append({**r.to_dict(),"exclusion_reason":str(exc)})
 curated=pd.DataFrame(rows)
 if not curated.empty:
  conflicts=curated.groupby("inchikey_connectivity").p_activity.agg(lambda x:x.max()-x.min()).ge(conflict_delta); conflict_ids=set(conflicts[conflicts].index); excluded.extend(curated[curated.inchikey_connectivity.isin(conflict_ids)].assign(exclusion_reason="conflicting_measurements").to_dict("records")); curated=curated[~curated.inchikey_connectivity.isin(conflict_ids)].sort_values(["assay_quality_score","record_weight"],ascending=False).drop_duplicates("canonical_isomeric_smiles")
 return curated.reset_index(drop=True),pd.DataFrame(excluded)
