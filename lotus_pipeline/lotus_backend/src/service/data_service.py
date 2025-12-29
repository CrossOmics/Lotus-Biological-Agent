import shutil
import uuid
from datetime import datetime
from pathlib import Path
from infrastructure.filesystem.storage import AssetStorage
from lotus.io import standardize_load


class DatasetService:
    def __init__(self):
        self.storage = AssetStorage()

    def import_dataset_from_local(self, local_file_path: str, dataset_name: str = None) -> dict:
        """
        Orchestrates the user upload process:
        1. Validates & Converts raw data to AnnData (Memory).
        2. Defines a secure relative path in the workspace.
        3. Saves the standardized .h5ad file.
        4. Returns metadata for DB insertion.

        Args:
            local_file_path: The absolute path chosen by the user (e.g., "C:/Downloads/raw.csv")
            dataset_name: User-defined name (e.g., "Patient 001")
        """

        # converts CSV/MTX -> Standard AnnData object in memory
        print(f"[Service] Ingesting file: {local_file_path}")
        adata = standardize_load(local_file_path)

        # Determine Storage Path (Relative Key)
        file_uuid = uuid.uuid4().hex[:8]
        original_name = Path(local_file_path).stem
        timestamp = datetime.now().strftime("%Y%m%d")

        # relative key convention: raw_data/ timestamp_original file name_last 8 chars of UUID.h5ad
        # example: "raw_data/20231229_MyData_a1b2.h5ad"
        relative_key = f"raw_data/{timestamp}_{original_name}_{file_uuid}.h5ad"

        # Persist to Workspace
        saved_path = self.storage.save_anndata(adata, relative_key)

        # Return Metadata (Ready for Database Insert)
        # TODO: replace plain hashmap with dataclass object
        return {
            "name": dataset_name or original_name,
            "relative_path": saved_path,
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "file_size_bytes": Path(self.storage.resolve(saved_path)).stat().st_size
        }
