from dataclasses import dataclass
from typing import Optional


@dataclass
class CreateProjectRequest():
    """Request body for creating a new project."""
    project_name: str
    local_file_path: str  # Absolute path to the raw data on user's disk
    organism: Optional[str] = None
    tissue_type: Optional[str] = None
    description: Optional[str] = None
