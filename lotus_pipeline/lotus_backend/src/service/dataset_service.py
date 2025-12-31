import uuid
from datetime import datetime
from pathlib import Path

from infrastructure.filesystem.constants.filesystem_constants import RAW_DATA
from infrastructure.filesystem.storage import AssetStorage
from lotus.io import standardize_load
from loguru import logger


class DatasetService:
    def __init__(self):
        self.storage = AssetStorage()

    def import_dataset_from_local(self, local_file_path: str, project_id: str, dataset_name: str = None) -> dict:
        """
        Orchestrates the user upload process:
        1. Validates & Converts raw data to AnnData (Memory).
        2. Defines a secure relative path in the workspace.
        3. Saves the standardized .h5ad file.
        4. Returns metadata for DB insertion.

        Args:
            local_file_path: The absolute path chosen by the user (e.g., "C:/Downloads/raw.csv")
            dataset_name: User-defined name (e.g., "Patient 001")
            project_id: user project id
        """

        # converts CSV/MTX to Standard AnnData object in memory
        logger.info(f"[Service] Ingesting file: {local_file_path}")
        adata = standardize_load(local_file_path)
        # Determine Storage Path (Relative Key)
        file_uuid = uuid.uuid4().hex[:5]
        original_name = Path(local_file_path).name
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # anndata filename key convention: {timestamp}_{original file name}_{last 5 chars of UUID}.h5ad
        filename = f"{timestamp}_{original_name}_{file_uuid}.h5ad"

        # Persist to Workspace
        saved_path = self.storage.save_anndata_project(adata, project_id, RAW_DATA, filename)

        # Return Metadata (Ready for Database Insert)
        # TODO: replace plain hashmap with dataclass object
        return {
            "name": dataset_name or original_name,
            "relative_path": saved_path,
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
        }
