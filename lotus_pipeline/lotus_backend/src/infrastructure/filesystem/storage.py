import json
import pickle
from pathlib import Path
from typing import Union, Any, Optional, List

from anndata import AnnData
from lotus.io import read_h5ad
from .constants.filesystem_constants import USER_PROJECT_ROOT
# Import the singleton instance from the context module
from ..workspace_context import workspace_path_manager
from loguru import logger

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
            # 1. Matplotlib Figure Support (Priority for Bio-Science)
            # Duck typing: Check if object has 'savefig' method (e.g., plt.figure)
            if hasattr(content, 'savefig'):
                # Default to tight bounding box for cleaner plots if not specified
                save_kwargs = {'bbox_inches': 'tight'}
                save_kwargs.update(kwargs) # Allow override
                content.savefig(target_abs_path, **save_kwargs)

            # 2. PIL Image Support
            # Duck typing: Check if object has 'save' method (but is not bytes/str)
            elif hasattr(content, 'save') and not isinstance(content, (bytes, str)):
                content.save(target_abs_path, **kwargs)

            # 3. JSON Dictionary
            elif isinstance(content, (dict, list)) and target_abs_path.suffix == '.json':
                with open(target_abs_path, 'w', encoding='utf-8') as f:
                    # Allow passing 'indent' in kwargs, default to 2
                    indent = kwargs.get('indent', 2)
                    json.dump(content, f, indent=indent)

            # 4. Binary Data (Bytes)
            # Handles raw image bytes, PDF streams, etc.
            elif isinstance(content, bytes):
                with open(target_abs_path, 'wb') as f:
                    f.write(content)

            # 5. Text Data (String)
            elif isinstance(content, str):
                with open(target_abs_path, 'w', encoding='utf-8') as f:
                    f.write(content)

            # 6. Fallback: Python Pickle
            else:
                with open(target_abs_path, 'wb') as f:
                    pickle.dump(content, f)

        except Exception as e:
            # Clean up partial files on failure
            if target_abs_path.exists():
                target_abs_path.unlink()
            raise IOError(f"Failed to save generic file ({relative_key}): {e}")

        return Path(relative_key).as_posix()

    def save_incremental_anndata(
            self,
            adata_source: AnnData,
            obs_cols: List[str],
            var_cols: List[str],
            relative_key: str
    ) -> Optional[str]:
        """
        Creates and saves a lightweight AnnData object containing only specific
        obs/var columns (Incremental Cache).

        Note: The X matrix is NOT saved (X=None). Indices are preserved automatically.

        Args:
            adata_source: The source AnnData object containing the data.
            obs_cols: List of column names from adata.obs to cache.
            var_cols: List of column names from adata.var to cache.
            relative_key: The destination path key.

        Returns:
            str: The saved path if successful, None otherwise.
        """
        try:
            # Create a lightweight object.
            # AnnData automatically preserves the index (barcodes/genes),
            # which is crucial for alignment during merging.
            adata_cache = AnnData(
                X=None,
                obs=adata_source.obs[obs_cols].copy(),
                var=adata_source.var[var_cols].copy()
            )

            logger.debug(
                f"[IO] Saving incremental cache to {relative_key} (Obs: {len(obs_cols)}, Var: {len(var_cols)})")
            return self.save_anndata(adata_cache, relative_key)

        except Exception as e:
            logger.warning(f"[IO] Failed to save incremental cache: {e}")
            return None

    def load_and_merge_anndata(self, adata_target: AnnData, relative_key: str) -> bool:
        """
        Loads an incremental cache file and merges its obs/var columns into the target AnnData.

        Strategies:
        1. Checks Index Alignment (Barcodes/Genes must match).
        2. Calculates column differences (only merges new or updated columns).
        3. Assigns columns individually to avoid AnnData structural errors.

        Args:
            adata_target: The main AnnData object in memory (modified in-place).
            relative_key: Path to the cache file.

        Returns:
            bool: True if merge was successful, False if file missing or index mismatch.
        """
        try:
            logger.debug(f"[IO] Attempting to merge cache from: {relative_key}")

            # 1. Load the cache
            # load_anndata handles path resolution and existence checks internally
            try:
                adata_cache = self.load_anndata(relative_key)
            except FileNotFoundError:
                logger.debug("[IO] Cache miss (file not found).")
                return False

            # 2. Safety Check: Index Alignment
            if not adata_target.obs_names.equals(adata_cache.obs_names):
                logger.warning(f"[IO] Cache index mismatch (Obs) for {relative_key}. Aborting merge.")
                return False

            # Note: We usually trust var_names (genes) to be consistent if obs aligns,
            # but strictly speaking, we could check var_names too.

            # 3. Calculate Delta (Difference)
            # We only merge columns that exist in cache but not in target (or to update them).
            new_obs_cols = adata_cache.obs.columns.difference(adata_target.obs.columns)
            new_var_cols = adata_cache.var.columns.difference(adata_target.var.columns)

            # 4. Safe Merge (Column-by-Column Assignment)
            # This is the most robust way to update AnnData.obs/var without breaking internal links.
            if len(new_obs_cols) > 0:
                for col in new_obs_cols:
                    adata_target.obs[col] = adata_cache.obs[col]

            if len(new_var_cols) > 0:
                for col in new_var_cols:
                    adata_target.var[col] = adata_cache.var[col]

            logger.info(
                f"[IO] Incremental Merge Success! Added {len(new_obs_cols)} obs cols, {len(new_var_cols)} var cols.")
            return True

        except Exception as e:
            logger.warning(f"[IO] Error during incremental merge: {e}")
            return False