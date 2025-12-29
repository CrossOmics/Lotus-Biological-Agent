import shutil
from pathlib import Path
from anndata import AnnData
from lotus.io import read_h5ad
# Import the singleton instance from the context module
from ..workspace_context import workspace_path_manager


class AssetStorage:
    """
    Handles physical file I/O operations for the biological assets.
    The Service Layer uses this class to save/load data without knowing absolute paths.
    """

    def save_anndata(self, adata: AnnData, relative_key: str) -> str:
        """
        Saves an AnnData object to the file system using a relative key.

        Args:
            adata (AnnData): The object to save.
            relative_key (str): The intended logical path (e.g., 'analysis/run_01.h5ad').

        Returns:
            str: The standardized POSIX relative path string to store in SQLite.
        """
        # 1. Resolve the absolute path using the Context Manager
        target_abs_path = workspace_path_manager.resolve(relative_key)

        # 2. Ensure parent directory exists (e.g., if key is 'analysis/batch1/run.h5ad')
        target_abs_path.parent.mkdir(parents=True, exist_ok=True)

        # 3. Write the file (using gzip to save disk space)
        print(f"[IO] Saving AnnData to: {target_abs_path}")
        adata.write_h5ad(target_abs_path, compression='gzip')

        # 4. Return POSIX path (forward slashes) for database consistency across OS
        return Path(relative_key).as_posix()

    def load_anndata(self, relative_key: str) -> AnnData:
        """
        Loads an AnnData object from the workspace.

        Args:
            relative_key (str): The relative path stored in SQLite.

        Returns:
            AnnData: The loaded data object.
        """
        source_abs_path = workspace_path_manager.resolve(relative_key)

        if not source_abs_path.exists():
            raise FileNotFoundError(f"Asset file missing at: {source_abs_path}")

        print(f"[IO] Loading AnnData from: {source_abs_path}")
        return read_h5ad(source_abs_path)

    def import_external_file(self, source_path_str: str) -> str:
        """
        Copies an external file (e.g., from Downloads) into the Workspace 'raw_data' folder.
        """
        src = Path(source_path_str)
        if not src.exists():
            raise FileNotFoundError("Source file not found.")

        # Construct a relative key for the internal workspace
        # e.g., raw_data/pbmc68k.h5ad
        relative_key = f"raw_data/{src.name}"

        # Resolve destination
        dest_abs_path = workspace_path_manager.resolve(relative_key)

        # Perform the copy
        shutil.copy2(src, dest_abs_path)

        return relative_key
