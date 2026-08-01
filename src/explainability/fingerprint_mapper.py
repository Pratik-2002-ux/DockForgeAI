from __future__ import annotations
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
def map_morgan_bits(dataset,bit_ids,radius=2,n_bits=2048):
 results=[]; generator=AllChem.GetMorganGenerator(radius=radius,fpSize=n_bits)
 for bit in bit_ids:
  examples=[]
  for _,row in dataset.iterrows():
   mol=Chem.MolFromSmiles(row.canonical_isomeric_smiles); output=AllChem.AdditionalOutput();output.AllocateBitInfoMap();fp=generator.GetFingerprint(mol,additionalOutput=output);info=output.GetBitInfoMap()
   if fp[bit]:
    atom,rad=info[bit][0]; env=Chem.FindAtomEnvironmentOfRadiusN(mol,rad,atom); fragment=Chem.MolFragmentToSmiles(mol,atomsToUse=list({a for bond in env for a in (mol.GetBondWithIdx(bond).GetBeginAtomIdx(),mol.GetBondWithIdx(bond).GetEndAtomIdx())}) or [atom],canonical=True); examples.append({"compound_id":row.compound_id,"atom_index":atom,"radius":rad,"fragment_smiles":fragment})
   if len(examples)>=5: break
  results.append({"bit":bit,"occurrence_examples":len(examples),"representative_compound_ids":";".join(e["compound_id"] for e in examples),"fragment_smiles":";".join(sorted({e["fragment_smiles"] for e in examples})),"collision_warning":"Hashed bits may represent multiple environments"})
 return results
