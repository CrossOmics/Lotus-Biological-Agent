from datetime import datetime
from typing import List, Optional
from peewee import DoesNotExist, IntegrityError

from ..model.project_meta_model import ProjectMeta


class ProjectMetaDAO:
    """
    Data Access Object (DAO) for managing ProjectMeta database operations.
    Encapsulates all database interactions related to projects to ensure
    consistent data handling and abstraction from the business logic.
    """

    @staticmethod
    def create_project(
            project_id: str,
            project_name: str,
            organism: Optional[str] = None,
            tissue_type: Optional[str] = None,
            description: Optional[str] = None,
            ext_info: Optional[str] = None
    ) -> Optional[ProjectMeta]:
        """
        Creates a new ProjectMeta record in the database.

        Args:
            project_id (str): Unique identifier for the project (e.g., 'p_20250101_12345').
            project_name (str): Human-readable name for the project.
            organism (str, optional): Species information (e.g., 'Human').
            tissue_type (str, optional): Tissue origin (e.g., 'Liver').
            description (str, optional): Project description.
            ext_info (str, optional): Additional metadata in JSON string format.

        Returns:
            ProjectMeta: The created project instance if successful.
            None: If creation fails (e.g., due to duplicate project_id).
        """
        try:
            # Current time for both create_time and update_time
            now = datetime.utcnow()

            project = ProjectMeta.create(
                project_id=project_id,
                project_name=project_name,
                organism=organism,
                tissue_type=tissue_type,
                description=description,
                ext_info=ext_info,
                create_time=now,
                update_time=now
            )
            return project
        except IntegrityError as e:
            print(f"[Error] Failed to create project '{project_id}': {e}")
            return None
        except Exception as e:
            print(f"[Error] Unexpected error creating project: {e}")
            return None

    @staticmethod
    def get_project_by_id(pk_id: int) -> Optional[ProjectMeta]:
        """
        Retrieves a single project by its primary key ID.

        Args:
            pk_id (int): The auto-increment primary key of the project.

        Returns:
            ProjectMeta: The found project instance.
            None: If no record is found.
        """
        try:
            return ProjectMeta.get_by_id(pk_id)
        except DoesNotExist:
            print(f"[Warning] Project with Primary Key {pk_id} not found.")
            return None

    @staticmethod
    def get_project_by_project_id(project_id: str) -> Optional[ProjectMeta]:
        """
        Retrieves a single project by its unique business string identifier.

        Args:
            project_id (str): The unique string identifier (e.g., 'p_2025...').

        Returns:
            ProjectMeta: The found project instance.
            None: If no record is found.
        """
        try:
            return ProjectMeta.get(ProjectMeta.project_id == project_id)
        except DoesNotExist:
            print(f"[Warning] Project with ID '{project_id}' not found.")
            return None

    @staticmethod
    def get_all_projects() -> List[ProjectMeta]:
        """
        Retrieves all projects from the database, ordered by update time (most recent first).

        Returns:
            List[ProjectMeta]: A list of all project instances.
        """
        try:
            return list(ProjectMeta.select().order_by(ProjectMeta.update_time.desc()))
        except Exception as e:
            print(f"[Error] Failed to retrieve all projects: {e}")
            return []

    @staticmethod
    def filter_projects(
            organism: Optional[str] = None,
            tissue_type: Optional[str] = None
    ) -> List[ProjectMeta]:
        """
        Filters projects based on organism and/or tissue type.

        Args:
            organism (str, optional): Filter by species (exact match).
            tissue_type (str, optional): Filter by tissue type (exact match).

        Returns:
            List[ProjectMeta]: A list of matching project instances ordered by recent updates.
        """
        try:
            query = ProjectMeta.select()

            if organism:
                query = query.where(ProjectMeta.organism == organism)

            if tissue_type:
                query = query.where(ProjectMeta.tissue_type == tissue_type)

            return list(query.order_by(ProjectMeta.update_time.desc()))
        except Exception as e:
            print(f"[Error] Failed to filter projects: {e}")
            return []

    @staticmethod
    def update_project(
            pk_id: int,
            project_name: Optional[str] = None,
            description: Optional[str] = None,
            organism: Optional[str] = None,
            tissue_type: Optional[str] = None,
            ext_info: Optional[str] = None
    ) -> bool:
        """
        Updates an existing project record. Automatically updates the 'update_time' field.
        Only updates fields that are provided (not None).

        Args:
            pk_id (int): Primary key of the project to update.
            project_name (str, optional): New name.
            description (str, optional): New description.
            organism (str, optional): New organism.
            tissue_type (str, optional): New tissue type.
            ext_info (str, optional): New ext info JSON string.

        Returns:
            bool: True if the update was successful, False otherwise.
        """
        try:
            project = ProjectMeta.get_by_id(pk_id)

            # Helper flag to check if we actually need to save and update timestamp
            changed = False

            if project_name is not None:
                project.project_name = project_name
                changed = True

            if description is not None:
                project.description = description
                changed = True

            if organism is not None:
                project.organism = organism
                changed = True

            if tissue_type is not None:
                project.tissue_type = tissue_type
                changed = True

            if ext_info is not None:
                project.ext_info = ext_info
                changed = True

            if changed:
                project.update_time = datetime.utcnow()
                project.save()
                return True

            return True  # No changes needed, but operation was valid

        except DoesNotExist:
            print(f"[Error] Cannot update: Project ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to update project {pk_id}: {e}")
            return False

    @staticmethod
    def delete_project(pk_id: int) -> bool:
        """
        Deletes a project record from the database.

        Args:
            pk_id (int): Primary key of the project to delete.

        Returns:
            bool: True if deletion was successful, False otherwise.
        """
        try:
            project = ProjectMeta.get_by_id(pk_id)
            project.delete_instance()
            print(f"[Info] Project ID {pk_id} deleted successfully.")
            return True
        except DoesNotExist:
            print(f"[Warning] Cannot delete: Project ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to delete project {pk_id}: {e}")
            return False
