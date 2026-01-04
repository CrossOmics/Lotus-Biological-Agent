import uuid
import os
from datetime import datetime
from pathlib import Path
from typing import Optional, Union

from loguru import logger

# Infrastructure & Constants
from infrastructure.filesystem.constants.filesystem_constants import RAW_DATA
from infrastructure.filesystem.storage import AssetStorage

# Database
from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.database.model.project_meta_model import ProjectMeta
from infrastructure.database.model.dataset_model import Dataset

from lotus.io import standardize_load


class DatasetService:
    """
    Service layer for managing Dataset lifecycle (Ingestion, Persistence, and DB Records).
    """

    def __init__(self):
        self.storage = AssetStorage()

    def import_dataset_from_local(
            self,
            local_file_path: str,
            project: Union[ProjectMeta, int],
            dataset_name: Optional[str] = None
    ) -> Optional[Dataset]:
        """
        Orchestrates the user upload process:
        1. Validates & Converts raw data to AnnData (Memory).
        2. Defines a secure relative path in the workspace.
        3. Saves the standardized .h5ad file.
        4. Inserts the record into the Dataset table via DAO.

        Args:
            local_file_path (str): The absolute path chosen by the user (e.g., "C:/Downloads/raw.csv").
            project (ProjectMeta | int): The project object or its Primary Key ID to link this dataset to.
            dataset_name (str, optional): User-defined name. Defaults to the original filename.

        Returns:
            Dataset: The created database model instance, or None if failed.
        """

        # Validate Local File
        local_path_obj = Path(local_file_path)
        if not local_path_obj.exists():
            logger.error(f"[Service] Local file not found: {local_file_path}")
            raise FileNotFoundError(f"File not found: {local_file_path}")

        # Ingest & Convert to Standard AnnData (in memory)
        try:
            adata = standardize_load(str(local_path_obj))
        except Exception as e:
            logger.error(f"[Service] Failed to load or standardize file: {e}")
            raise ValueError(f"Invalid data format: {e}")

        # 2Prepare Metadata & Naming Conventions
        file_uuid = uuid.uuid4().hex[:5]
        original_name = local_path_obj.stem
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Business Logic: Construct unique IDs and Filenames
        # Naming convention: {timestamp}_{original_name}_{uuid}.h5ad
        filename = f"{timestamp}_{original_name}_{file_uuid}.h5ad"

        # Dataset Business ID
        dataset_business_id = f"ds_{timestamp}_{file_uuid}"

        # Determine Display Name
        final_dataset_name = dataset_name if dataset_name else original_name

        # Persist File to Workspace (Storage Layer)
        if isinstance(project, ProjectMeta):
            project_str_id = project.project_id
        else:
            project_str_id = project

        try:
            saved_relative_path = self.storage.save_anndata_project(
                adata=adata,
                project_id=project_str_id,
                file_type=RAW_DATA,
                file_name=filename
            )
        except Exception as e:
            logger.error(f"[Service] Storage write failed: {e}")
            raise e

        logger.info(f"[Service] Inserting dataset record into DB: {final_dataset_name}")
        dataset_record = DatasetDAO.create_dataset(
            project=project_str_id,
            dataset_id=dataset_business_id,
            dataset_name=final_dataset_name,
            dataset_path=saved_relative_path
        )

        if dataset_record:
            logger.success(f"[Service] Dataset imported successfully. ID: {dataset_record.id}")
            return dataset_record
        else:
            logger.error("[Service] Database insertion failed.")
            return None
