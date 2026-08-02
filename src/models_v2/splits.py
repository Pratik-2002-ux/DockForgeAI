"""Leakage-safe scaffold, temporal, assay, and source splits."""
from __future__ import annotations
import numpy as np,pandas as pd
def scaffold_split(frame,seed=2026,fractions=(.7,.15,.15)):
 rng=np.random.default_rng(seed); groups=frame.groupby("murcko_scaffold_smiles").size().sort_values(ascending=False); buckets=[[],[],[]]; counts=np.zeros(3); targets=np.asarray(fractions)*len(frame)
 for scaffold,n in groups.items():
  ratios=np.divide(counts,targets,out=np.zeros(3),where=targets>0); choices=np.flatnonzero(ratios==ratios.min()); bucket=int(rng.choice(choices)); buckets[bucket].append(scaffold); counts[bucket]+=n
 mapping={s:name for values,name in zip(buckets,("train","validation","test")) for s in values}; return frame.murcko_scaffold_smiles.map(mapping)
def time_split(frame,fractions=(.7,.15,.15)):
 order=pd.to_datetime(frame.assay_date,errors="coerce").sort_values().index; out=pd.Series(index=frame.index,dtype=str); a=int(len(frame)*fractions[0]); b=int(len(frame)*(fractions[0]+fractions[1])); out.loc[order[:a]]="train"; out.loc[order[a:b]]="validation"; out.loc[order[b:]]="test"; return out
def leave_group_out(frame,column): return {str(group):(frame[column].ne(group),frame[column].eq(group)) for group in sorted(frame[column].dropna().unique())}
def assert_no_scaffold_leakage(frame,split): return int(pd.DataFrame({"s":frame.murcko_scaffold_smiles,"split":split}).groupby("s").split.nunique().gt(1).sum())==0
