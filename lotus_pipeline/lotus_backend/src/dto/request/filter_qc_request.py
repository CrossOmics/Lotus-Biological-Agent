from pydantic import BaseModel, Field
from typing import Optional


class FilterQCRequest(BaseModel):
    """
    Request DTO for applying QC filters.
    """
    project_id: str = Field(..., description="Project Business ID")
    dataset_id: str = Field(..., description="Source Dataset Business ID (Root for this snapshot)")

    # Filter thresholds
    min_genes: int = Field(200, description="Minimum number of genes per cell")
    min_cells: int = Field(3, description="Minimum number of cells per gene")

    # Optional filters
    max_counts: Optional[int] = Field(None, description="Maximum total counts per cell")
    pct_mt_max: Optional[float] = Field(None, description="Maximum mitochondrial percentage (0-100)")
    pct_hb_max: Optional[float] = Field(None, description="Maximum hemoglobin percentage (0-100)")