"""Curated physicochemical and topological descriptor features."""
from __future__ import annotations
import math
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, GraphDescriptors, Lipinski, rdMolDescriptors

BASE = ["molecular_weight","exact_molecular_weight","logp","tpsa","hbd","hba","rotatable_bonds","ring_count","aromatic_ring_count","formal_charge","heavy_atom_count","fraction_csp3","qed","pains_alert_count","brenk_alert_count"]
EXTRA = ["MolMR","LabuteASA","BalabanJ","BertzCT","HallKierAlpha","Kappa1","Kappa2","Kappa3","Chi0","Chi1","Chi2n","Chi3n","NumAliphaticRings","NumSaturatedRings","NumHeterocycles","NumAmideBonds","NumBridgeheadAtoms","NumSpiroAtoms"]

def descriptor_frame(dataset: pd.DataFrame) -> pd.DataFrame:
    """Return finite, ordered descriptors without metadata or labels."""
    rows=[]
    for smiles in dataset["canonical_isomeric_smiles"]:
        mol=Chem.MolFromSmiles(smiles)
        if mol is None: raise ValueError(f"Invalid standardized SMILES: {smiles}")
        rows.append({"MolMR":Descriptors.MolMR(mol),"LabuteASA":rdMolDescriptors.CalcLabuteASA(mol),"BalabanJ":GraphDescriptors.BalabanJ(mol),"BertzCT":GraphDescriptors.BertzCT(mol),"HallKierAlpha":GraphDescriptors.HallKierAlpha(mol),"Kappa1":GraphDescriptors.Kappa1(mol),"Kappa2":GraphDescriptors.Kappa2(mol),"Kappa3":GraphDescriptors.Kappa3(mol),"Chi0":GraphDescriptors.Chi0(mol),"Chi1":GraphDescriptors.Chi1(mol),"Chi2n":GraphDescriptors.Chi2n(mol),"Chi3n":GraphDescriptors.Chi3n(mol),"NumAliphaticRings":rdMolDescriptors.CalcNumAliphaticRings(mol),"NumSaturatedRings":rdMolDescriptors.CalcNumSaturatedRings(mol),"NumHeterocycles":rdMolDescriptors.CalcNumHeterocycles(mol),"NumAmideBonds":rdMolDescriptors.CalcNumAmideBonds(mol),"NumBridgeheadAtoms":rdMolDescriptors.CalcNumBridgeheadAtoms(mol),"NumSpiroAtoms":rdMolDescriptors.CalcNumSpiroAtoms(mol)})
    frame=pd.concat([dataset[BASE].astype(float).reset_index(drop=True),pd.DataFrame(rows)],axis=1)
    bad=[column for column in frame if not frame[column].map(math.isfinite).all()]
    if bad: raise ValueError(f"Non-finite descriptors: {bad}")
    return frame

def descriptor_quality(frame: pd.DataFrame) -> dict:
    constants=[c for c in frame if frame[c].nunique(dropna=False)<=1]
    correlation=frame.corr().abs(); pairs=[{"feature_a":a,"feature_b":b,"absolute_correlation":float(correlation.loc[a,b])} for i,a in enumerate(frame) for b in frame[i+1:] if correlation.loc[a,b]>=0.95]
    return {"row_count":len(frame),"feature_count":len(frame.columns),"constant_columns":constants,"missing_rate":{c:float(frame[c].isna().mean()) for c in frame},"highly_correlated_pairs":pairs}
