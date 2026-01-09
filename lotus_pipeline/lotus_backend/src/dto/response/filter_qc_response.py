from pydantic import BaseModel, Field

class FilterQCResponse(BaseModel):
    """
    Response DTO after applying filtering.
    """
    snapshot_id: str = Field(..., description="Unique ID of the generated snapshot")
    snapshot_path: str = Field(..., description="Path to the saved H5AD snapshot")
    n_obs_remaining: int = Field(..., description="Number of cells remaining after filtering")
    n_vars_remaining: int = Field(..., description="Number of genes remaining after filtering")