from typing import Optional
from pydantic import BaseModel

class RunNeighborsRequest(BaseModel):
    project_id: str
    dataset_id: str
    source_snapshot_id: Optional[str] = None
    n_neighbors: int = 15
    n_pcs: int = 30
