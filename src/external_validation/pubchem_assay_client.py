"""Rate-limited clients for official PubChem and NCBI services."""
from __future__ import annotations
import hashlib, json, time
from io import StringIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import pandas as pd

EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUG="https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUG_VIEW="https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

class PubChemRequestError(RuntimeError): pass

class PubChemAssayClient:
    def __init__(self,cache_dir:Path,offline=False,timeout=45,retries=4,rate_limit=0.34,user_agent="DockForgeAI/1.0 PubChem external validation"):
        self.cache_dir=Path(cache_dir); self.cache_dir.mkdir(parents=True,exist_ok=True); self.offline=offline; self.timeout=timeout; self.retries=retries; self.rate_limit=rate_limit; self.user_agent=user_agent; self.last_request=0.
    def _path(self,url): return self.cache_dir/(hashlib.sha256(url.encode()).hexdigest()+".cache")
    def get(self,url,refresh=False):
        path=self._path(url)
        if path.exists() and not refresh: return path.read_bytes()
        if self.offline: raise PubChemRequestError(f"Offline cache miss: {url}")
        for attempt in range(self.retries):
            try:
                delay=self.rate_limit-(time.monotonic()-self.last_request)
                if delay>0: time.sleep(delay)
                with urlopen(Request(url,headers={"User-Agent":self.user_agent,"Accept":"application/json,text/csv"}),timeout=self.timeout) as response: payload=response.read()
                self.last_request=time.monotonic(); path.write_bytes(payload); return payload
            except Exception as exc:
                if attempt+1==self.retries: raise PubChemRequestError(f"Request failed after {self.retries} attempts: {url}") from exc
                time.sleep(2**attempt)
    def get_json(self,url,refresh=False): return json.loads(self.get(url,refresh).decode())
    def search_assays(self,term,retmax=200,refresh=False):
        url=f"{EUTILS}/esearch.fcgi?"+urlencode({"db":"pcassay","term":term,"retmode":"json","retmax":retmax}); result=self.get_json(url,refresh)["esearchresult"]; return [int(x) for x in result.get("idlist",[])],int(result.get("count",0)),url
    def assay_summaries(self,aids,refresh=False,batch_size=100):
        rows=[]
        for start in range(0,len(aids),batch_size):
            values=aids[start:start+batch_size]; url=f"{EUTILS}/esummary.fcgi?"+urlencode({"db":"pcassay","id":",".join(map(str,values)),"retmode":"json"}); result=self.get_json(url,refresh).get("result",{}); rows.extend(result[str(a)] for a in values if str(a) in result)
        return rows
    def assay_description(self,aid,refresh=False): return self.get_json(f"{PUG}/assay/aid/{aid}/description/JSON",refresh)
    def assay_results(self,aid,refresh=False): return pd.read_csv(StringIO(self.get(f"{PUG}/assay/aid/{aid}/CSV",refresh).decode()))
    def compound_properties(self,cids,refresh=False,batch_size=100):
        properties="CanonicalSMILES,IsomericSMILES,InChI,InChIKey,MolecularFormula"
        rows=[]
        for start in range(0,len(cids),batch_size):
            chunk=','.join(map(str,cids[start:start+batch_size])); url=f"{PUG}/compound/cid/{chunk}/property/{properties}/JSON"; rows.extend(self.get_json(url,refresh).get("PropertyTable",{}).get("Properties",[]))
        return pd.DataFrame(rows)
