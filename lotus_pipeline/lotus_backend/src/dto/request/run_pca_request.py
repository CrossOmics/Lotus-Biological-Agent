from typing import Optional

from pydantic import BaseModel

class RunPCARequest(BaseModel):
    project_id: str
    dataset_id: str
    source_snapshot_id: Optional[str] = None
    n_comps: int = 50
    svd_solver: str = 'arpack'