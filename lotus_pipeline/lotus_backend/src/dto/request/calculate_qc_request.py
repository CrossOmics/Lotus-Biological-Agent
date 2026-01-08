from pydantic import BaseModel, Field
from typing import Optional, Dict


class CalculateQCRequest(BaseModel):
    """
    Request DTO for QC calculation parameters.
    """
    project_id: str = Field(..., description="Unique project business ID")
    dataset_id: str = Field(..., description="Unique dataset business ID")

    # Optional parameters with defaults
    organism: Optional[str] = Field("Human", description="Organism type (Human/Mouse) for gene prefix detection")
    custom_prefixes: Optional[Dict[str, str]] = Field(None, description="Custom gene prefixes (e.g., {'mt': 'MT-'})")