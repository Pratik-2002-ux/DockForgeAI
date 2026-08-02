"""Molecular vector and graph representations."""
from __future__ import annotations
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import AllChem,MACCSkeys
from src.features.descriptor_features import descriptor_frame
def fingerprints(smiles,kind="morgan",n_bits=2048):
 rows=[]; gen=AllChem.GetMorganGenerator(radius=3 if kind=="ecfp6" else 2,fpSize=n_bits,includeChirality=True)
 for value in smiles:
  mol=Chem.MolFromSmiles(value); fp=MACCSkeys.GenMACCSKeys(mol) if kind=="maccs" else gen.GetFingerprint(mol); arr=np.zeros(fp.GetNumBits(),np.uint8); DataStructs.ConvertToNumpyArray(fp,arr); rows.append(arr)
 return np.asarray(rows)
def descriptors(frame): return descriptor_frame(frame)
def graph_from_smiles(smiles):
 mol=Chem.MolFromSmiles(smiles); atoms=np.asarray([[a.GetAtomicNum(),a.GetTotalDegree(),a.GetFormalCharge(),int(a.GetIsAromatic()),int(a.GetChiralTag())] for a in mol.GetAtoms()],dtype=float); edges=[]; bonds=[]
 for b in mol.GetBonds():
  i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx(); feature=[float(b.GetBondTypeAsDouble()),int(b.GetIsAromatic()),int(b.IsInRing()),int(b.GetStereo())]; edges.extend([[i,j],[j,i]]); bonds.extend([feature,feature])
 return {"atom_features":atoms,"edge_index":np.asarray(edges,dtype=int).T if edges else np.empty((2,0),int),"bond_features":np.asarray(bonds,dtype=float)}
