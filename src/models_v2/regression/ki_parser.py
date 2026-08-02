from __future__ import annotations
import math
from src.models_v2.data.measurement_parser import parse_measurement
def parse_ki(value,units="nM"):
 rel,num=parse_measurement(value); factors={"nm":1.,"um":1000.,"µm":1000.,"μm":1000.,"pm":.001,"mm":1_000_000.}; key=str(units).lower()
 if num is None or key not in factors:return rel,None
 return rel,num*factors[key]
def pki_from_nm(value):
 value=float(value)
 if not math.isfinite(value) or value<=0: raise ValueError("Ki must be positive and finite")
 return 9-math.log10(value)
