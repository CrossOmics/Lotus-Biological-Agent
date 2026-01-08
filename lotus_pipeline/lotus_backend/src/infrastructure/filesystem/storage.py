import json
import pickle
from pathlib import Path
from typing import Union, Any

from anndata import AnnData
from lotus.io import read_h5ad
from .constants.filesystem_constants import USER_PROJECT_ROOT
# Import the singleton instance from the context module
from ..workspace_context import workspace_path_manager


class AssetStorage:
    """
    Handles physical file I/O operations for biological assets.
    The Service Layer uses this class to save/load data without managing absolute paths directly.
    """

    def save_anndata_project(self, adata: AnnData, project_id: str, file_type: str, file_name: str) -> str:
        """
        Constructs the relative path and saves the AnnData object.

        Args:
            adata (AnnData): The adata object to save.
            project_id (str): The unique business ID of the current project (e.g., 'p_2025...').
            file_type (str): Sub-folder name (e.g., 'raw_data', 'analysis').
            file_name (str): The filename without extension.

        Returns:
            str: The standardized POSIX relative path key for database storage.
        """
        # Construct the logical relative path using USER_PROJECT_ROOT
        relative_path_obj = Path(USER_PROJECT_ROOT) / project_id / file_type / file_name

        # Convert to POSIX style (forward slashes) for DB consistency
        relative_key = relative_path_obj.as_posix()

        return self.save_anndata(adata, relative_key)

    def save_anndata(self, adata: AnnData, relative_key: str) -> str:
        """
        Saves an AnnData object to the file system using a relative key.

        Args:
            adata (AnnData): The object to save.
            relative_key (str): The intended logical path.

        Returns:
            str: The standardized POSIX relative path string.
        """
        # 1. Resolve the absolute path using the Context Manager
        target_abs_path = workspace_path_manager.resolve(relative_key)

        # 2. Ensure parent directory exists
        target_abs_path.parent.mkdir(parents=True, exist_ok=True)

        # 3. Write the file (using gzip to save disk space)
        print(f"[IO] Saving AnnData to: {target_abs_path}")
        try:
            adata.write_h5ad(target_abs_path, compression='gzip')
        except Exception as e:
            # Clean up if write fails to avoid corrupt files
            if target_abs_path.exists():
                target_abs_path.unlink()
            raise IOError(f"Failed to write AnnData file: {e}")

        # 4. Return POSIX path
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

    def save_file(self, content: Union[str, bytes, dict, Any], relative_key: str, **kwargs) -> str:
        """
        Generic method to save arbitrary files to the workspace.
        Supports: JSON, Strings, Bytes, Matplotlib Figures, PIL Images, and Pickle objects.

        Args:
            content: The data object to save.
            relative_key: The logical path (e.g., 'qc/violin_plot.pdf').
            **kwargs: Additional arguments passed to the underlying save method (e.g., indent=2 for JSON).

        Returns:
            str: The POSIX relative path.
        """
        target_abs_path = workspace_path_manager.resolve(relative_key)
        target_abs_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"[IO] Saving generic file to: {target_abs_path}")

        try:
            # --- 1. Matplotlib Figure Support (Priority for Bio-Science) ---
            # Duck typing: Check if object has 'savefig' method (e.g., plt.figure)
            if hasattr(content, 'savefig'):
                # Default to tight bounding box for cleaner plots if not specified
                save_kwargs = {'bbox_inches': 'tight'}
                save_kwargs.update(kwargs) # Allow override
                content.savefig(target_abs_path, **save_kwargs)

            # --- 2. PIL Image Support ---
            # Duck typing: Check if object has 'save' method (but is not bytes/str)
            elif hasattr(content, 'save') and not isinstance(content, (bytes, str)):
                content.save(target_abs_path, **kwargs)

            # --- 3. JSON Dictionary ---
            elif isinstance(content, (dict, list)) and target_abs_path.suffix == '.json':
                with open(target_abs_path, 'w', encoding='utf-8') as f:
                    # Allow passing 'indent' in kwargs, default to 2
                    indent = kwargs.get('indent', 2)
                    json.dump(content, f, indent=indent)

            # --- 4. Binary Data (Bytes) ---
            # Handles raw image bytes, PDF streams, etc.
            elif isinstance(content, bytes):
                with open(target_abs_path, 'wb') as f:
                    f.write(content)

            # --- 5. Text Data (String) ---
            elif isinstance(content, str):
                with open(target_abs_path, 'w', encoding='utf-8') as f:
                    f.write(content)

            # --- 6. Fallback: Python Pickle ---
            else:
                with open(target_abs_path, 'wb') as f:
                    pickle.dump(content, f)

        except Exception as e:
            # Clean up partial files on failure
            if target_abs_path.exists():
                target_abs_path.unlink()
            raise IOError(f"Failed to save generic file ({relative_key}): {e}")

        return Path(relative_key).as_posix()
