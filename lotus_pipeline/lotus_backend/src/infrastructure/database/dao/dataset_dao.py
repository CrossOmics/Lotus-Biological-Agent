from datetime import datetime
from typing import List, Optional, Dict, Any
from peewee import DoesNotExist, IntegrityError

from ..model.dataset_model import Dataset
from ..model.project_meta_model import ProjectMeta


class DatasetDAO:
    """
    Data Access Object (DAO) for managing Dataset database operations.
    Handles interaction with the 'dataset' table, linking raw files to projects.
    """

    @staticmethod
    def create_dataset(
            project: ProjectMeta,
            dataset_id: str,
            dataset_name: str,
            ext_info: Optional[Dict[str, Any]] = None
    ) -> Optional[Dataset]:
        """
        Creates a new Dataset record.

        Args:
            project (ProjectMeta): The parent project instance this dataset belongs to.
            dataset_id (str): Unique business identifier (e.g., 'dataset_20250101_name_1234').
            dataset_name (str): Display name for the dataset.
            ext_info (dict, optional): Dictionary containing file metadata (e.g., file path, size).

        Returns:
            Dataset: The created dataset instance if successful.
            None: If creation fails (e.g., integrity error).
        """
        try:
            dataset = Dataset.create(
                project_primary_id=project,
                dataset_id=dataset_id,
                dataset_name=dataset_name,
                uploaded_time=datetime.utcnow(),
                ext_info=ext_info or {}
            )
            return dataset
        except IntegrityError as e:
            print(f"[Error] Failed to create dataset '{dataset_id}': {e}")
            return None
        except Exception as e:
            print(f"[Error] Unexpected error creating dataset: {e}")
            return None

    @staticmethod
    def get_dataset_by_id(pk_id: int) -> Optional[Dataset]:
        """
        Retrieves a single dataset by its primary key ID.

        Args:
            pk_id (int): The auto-increment primary key.

        Returns:
            Dataset: The found dataset instance or None.
        """
        try:
            return Dataset.get_by_id(pk_id)
        except DoesNotExist:
            print(f"[Warning] Dataset with PK {pk_id} not found.")
            return None

    @staticmethod
    def get_dataset_by_business_id(dataset_id: str) -> Optional[Dataset]:
        """
        Retrieves a dataset by its unique business string ID.

        Args:
            dataset_id (str): The unique string identifier.

        Returns:
            Dataset: The found dataset instance or None.
        """
        try:
            return Dataset.get(Dataset.dataset_id == dataset_id)
        except DoesNotExist:
            print(f"[Warning] Dataset with business ID '{dataset_id}' not found.")
            return None

    @staticmethod
    def get_datasets_by_project(project_pk: int) -> List[Dataset]:
        """
        Retrieves all datasets linked to a specific project.

        Args:
            project_pk (int): The primary key of the project.

        Returns:
            List[Dataset]: List of dataset instances, ordered by upload time.
        """
        try:
            query = (Dataset
                     .select()
                     .where(Dataset.project_primary_id == project_pk)
                     .order_by(Dataset.uploaded_time.desc()))
            return list(query)
        except Exception as e:
            print(f"[Error] Failed to retrieve datasets for project {project_pk}: {e}")
            return []

    @staticmethod
    def update_dataset(
            pk_id: int,
            dataset_name: Optional[str] = None,
            ext_info_update: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        Updates an existing dataset record.

        Args:
            pk_id (int): Primary key of the dataset to update.
            dataset_name (str, optional): New display name.
            ext_info_update (dict, optional): Dictionary of metadata to merge/replace.

        Returns:
            bool: True if successful, False otherwise.
        """
        try:
            dataset = Dataset.get_by_id(pk_id)

            if dataset_name is not None:
                dataset.dataset_name = dataset_name

            # Since ext_info is a JSONField, we can update it directly
            if ext_info_update is not None:
                # Assuming complete replacement for simplicity,
                # logic could be added to merge dictionaries if needed.
                dataset.ext_info = ext_info_update

            dataset.save()
            return True
        except DoesNotExist:
            print(f"[Error] Cannot update: Dataset ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to update dataset {pk_id}: {e}")
            return False

    @staticmethod
    def delete_dataset(pk_id: int) -> bool:
        """
        Deletes a dataset record from the database.

        Args:
            pk_id (int): Primary key of the dataset to delete.

        Returns:
            bool: True if deletion was successful, False otherwise.
        """
        try:
            dataset = Dataset.get_by_id(pk_id)
            dataset.delete_instance()
            print(f"[Info] Dataset ID {pk_id} deleted successfully.")
            return True
        except DoesNotExist:
            print(f"[Warning] Cannot delete: Dataset ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to delete dataset {pk_id}: {e}")
            return False
