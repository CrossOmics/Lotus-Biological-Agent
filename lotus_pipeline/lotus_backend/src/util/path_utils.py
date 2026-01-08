from pathlib import Path
from typing import Optional
from infrastructure.filesystem.constants.filesystem_constants import USER_PROJECT_ROOT


def get_project_relative(project_id: str, subspace: str, filename: Optional[str] = None) -> str:
    """
    Build a POSIX-style relative_key for saving files inside a project workspace.

    Example:
        get_project_relative("projects", "proj123", "qc", "plot.pdf")
        -> "projects/proj123/qc/plot.pdf"
    """
    parts = [USER_PROJECT_ROOT, project_id, subspace]
    if filename:
        parts.append(filename)

    return Path(*parts).as_posix()
