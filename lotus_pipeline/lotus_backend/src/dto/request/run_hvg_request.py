from pydantic import BaseModel, Field
from typing import Optional, Literal


class RunHVGRequest(BaseModel):
    project_id: str
    dataset_id: str
    source_snapshot_id: Optional[str] = None
    n_top_genes: int = Field(2000, ge=500, le=10000, description="Target number of HVGs")
    flavor: Literal['seurat', 'cell_ranger', 'seurat_v3'] = 'seurat'
    target_sum: float = 1e4
