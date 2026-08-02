"""Representation robustness transformations."""
from __future__ import annotations
import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
def perturb_descriptors(values,scale=.01,seed=2026): return np.asarray(values)+np.random.default_rng(seed).normal(0,scale,np.asarray(values).shape)
def randomized_smiles(smiles,count=5,seed=2026):
 mol=Chem.MolFromSmiles(smiles); return [Chem.MolToSmiles(mol,canonical=False,doRandom=True) for _ in range(count)]
def canonical_tautomer(smiles): return Chem.MolToSmiles(rdMolStandardize.TautomerEnumerator().Canonicalize(Chem.MolFromSmiles(smiles)),isomericSmiles=True)
def remove_stereochemistry(smiles):
 mol=Chem.MolFromSmiles(smiles); Chem.RemoveStereochemistry(mol); return Chem.MolToSmiles(mol)
