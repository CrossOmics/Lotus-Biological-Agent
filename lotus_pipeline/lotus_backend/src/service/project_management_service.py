import uuid
import time
from pathlib import Path
from typing import Optional, List

from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO
from infrastructure.database.model.project_meta_model import ProjectMeta


class ProjectManagementService:
    """
    Service layer for Project Metadata management.
    Handles business logic such as ID generation, path resolution, and bridging
    data access to the controller layer.
    """

    def __init__(self):
        pass

    def create_project(
            self,
            project_name: str,
            organism: Optional[str] = None,
            tissue_type: Optional[str] = None,
            description: Optional[str] = None
    ) -> Optional[ProjectMeta]:
        """
        Creates a new project. Generates a unique Project ID and assigns a workspace path.

        Args:
            project_name: User-defined name.
            organism: Biological organism (e.g., Human).
            tissue_type: Tissue origin (e.g., Liver).
            description: Optional description.

        Returns:
            ProjectMeta: The created project object, or None if failed.
        """
        # 1. Generate Unique Business ID
        # Format: p_<timestamp>_<short_uuid>
        timestamp = int(time.time())
        short_uuid = uuid.uuid4().hex[:5]
        project_id = f"p_{timestamp}_{short_uuid}"

        # 2. Determine Project Root Path
        # Logic: USER_PROJECT_ROOT / project_id
        # We store the relative path string (compatible with pathlib logic later)
        from infrastructure.filesystem.constants.filesystem_constants import USER_PROJECT_ROOT
        project_path_obj = Path(USER_PROJECT_ROOT) / project_id
        project_path_str = project_path_obj.as_posix()

        # 3. Persist to Database via DAO
        project = ProjectMetaDAO.create_project(
            project_id=project_id,
            project_name=project_name,
            project_path=project_path_str,
            organism=organism,
            tissue_type=tissue_type,
            description=description
        )

        return project

    def get_project_by_id(self, project_id: str) -> Optional[ProjectMeta]:
        """
        Retrieves a project by its business ID (e.g., 'p_2025...').
        """
        return ProjectMetaDAO.get_project_by_project_id(project_id)

    def get_project_by_pk(self, pk_id: int) -> Optional[ProjectMeta]:
        """
        Retrieves a project by its database primary key.
        """
        return ProjectMetaDAO.get_project_by_id(pk_id)

    def list_all_projects(self) -> List[ProjectMeta]:
        """
        Returns a list of all projects, ordered by last update.
        """
        return ProjectMetaDAO.get_all_projects()

    def update_project_details(
            self,
            pk_id: int,
            name: Optional[str] = None,
            description: Optional[str] = None,
            organism: Optional[str] = None,
            tissue_type: Optional[str] = None
    ) -> bool:
        """
        Updates project metadata.
        """
        return ProjectMetaDAO.update_project(
            pk_id=pk_id,
            project_name=name,
            description=description,
            organism=organism,
            tissue_type=tissue_type
        )

    def delete_project(self, pk_id: int) -> bool:
        """
        Deletes a project record.

        Note: This only deletes the DB record. File cleanup should strictly
        be handled by an orchestration layer or a separate cleanup service
        to avoid data loss accidents.
        """
        return ProjectMetaDAO.delete_project(pk_id)
