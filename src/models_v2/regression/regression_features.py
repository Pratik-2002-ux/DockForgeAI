from __future__ import annotations
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import AllChem,Descriptors,MACCSkeys
DESCRIPTORS=[("MolWt",Descriptors.MolWt),("MolLogP",Descriptors.MolLogP),("TPSA",Descriptors.TPSA),("HBD",Descriptors.NumHDonors),("HBA",Descriptors.NumHAcceptors),("RotB",Descriptors.NumRotatableBonds),("Rings",Descriptors.RingCount),("FracCSP3",Descriptors.FractionCSP3),("QED",Descriptors.qed)]
def matrices(smiles):
 mols=[Chem.MolFromSmiles(s) for s in smiles]; desc=np.asarray([[fn(m) for _,fn in DESCRIPTORS] for m in mols],float)
 def fp(radius,bits=2048):
  arr=[]
  for m in mols:
   v=AllChem.GetMorganGenerator(radius=radius,fpSize=bits).GetFingerprint(m); a=np.zeros(bits,dtype=np.uint8); DataStructs.ConvertToNumpyArray(v,a); arr.append(a)
  return np.asarray(arr)
 m2=fp(2);m3=fp(3);mac=[]
 for m in mols:
  a=np.zeros(167,dtype=np.uint8);DataStructs.ConvertToNumpyArray(MACCSkeys.GenMACCSKeys(m),a);mac.append(a)
 return {"descriptors":desc,"morgan_r2":m2,"morgan_r3":m3,"maccs":np.asarray(mac),"fused_r2":np.hstack([desc,m2]),"fused_r3":np.hstack([desc,m3])}
