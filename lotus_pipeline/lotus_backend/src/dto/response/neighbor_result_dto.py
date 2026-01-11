from typing import Optional
from pydantic import BaseModel


class NeighborsResultDTO(BaseModel):
    snapshot_path: str
    snapshot_id: str
    n_connectivities: int
