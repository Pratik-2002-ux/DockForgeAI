import json
from pathlib import Path
import pytest
from src.datasets.chembl_client import ChEMBLClient, ChEMBLClientError

class Response:
    def __init__(self, payload): self.payload = payload
    def raise_for_status(self): return None
    def json(self): return self.payload
class Session:
    def __init__(self, pages): self.pages=iter(pages); self.headers={}; self.calls=0
    def get(self,*args,**kwargs): self.calls+=1; return Response(next(self.pages))

def test_pagination_and_cache(tmp_path: Path) -> None:
    session=Session([{"activities":[{"id":1}],"page_meta":{"next":"x"}},{"activities":[{"id":2}],"page_meta":{"next":None}}])
    assert len(ChEMBLClient(tmp_path,page_size=1,session=session).fetch_activities()) == 2 and session.calls == 2
def test_offline_without_cache_fails(tmp_path: Path) -> None:
    with pytest.raises(ChEMBLClientError): ChEMBLClient(tmp_path).fetch_activities(offline=True)
