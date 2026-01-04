from datetime import datetime
from typing import List, Optional, Dict, Any, Union
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
            project: Union[ProjectMeta, str],
            dataset_id: str,
            dataset_name: str,
            dataset_path: str,
            ext_info: Optional[Dict[str, Any]] = None
    ) -> Optional[Dataset]:
        """
        Creates a new Dataset record.

        Args:
            project (ProjectMeta | str): The parent project instance or project_id string.
            dataset_id (str): Unique business identifier.
            dataset_name (str): Display name for the dataset.
            dataset_path (str): Relative file path in storage.
            ext_info (dict, optional): Dictionary containing file metadata.

        Returns:
            Dataset: The created dataset instance if successful.
        """
        try:
            dataset = Dataset.create(
                project_id=project,  # Matches the ForeignKeyField name in Model
                dataset_id=dataset_id,
                dataset_name=dataset_name,
                dataset_path=dataset_path,  # [Fix] Insert the path
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
                     .where(Dataset.project_id == project_pk)
                     .order_by(Dataset.uploaded_time.desc()))
            return list(query)
        except Exception as e:
            print(f"[Error] Failed to retrieve datasets for project {project_pk}: {e}")
            return []

    @staticmethod
    def update_dataset(
            pk_id: int,
            dataset_name: Optional[str] = None,
            dataset_path: Optional[str] = None,  # [New Field] Allow update
            ext_info_update: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        Updates an existing dataset record.
        """
        try:
            dataset = Dataset.get_by_id(pk_id)

            if dataset_name is not None:
                dataset.dataset_name = dataset_name

            if dataset_path is not None:
                dataset.dataset_path = dataset_path

            if ext_info_update is not None:
                # Merge or replace logic for JSON
                current_info = dataset.ext_info or {}
                current_info.update(ext_info_update)
                dataset.ext_info = current_info

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
