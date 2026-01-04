from dataclasses import dataclass


@dataclass
class ProjectResponse():
    """Standard response after project creation."""
    project_id: str
    project_name: str
    dataset_id: str
    dataset_name: str
    workspace_path: str
