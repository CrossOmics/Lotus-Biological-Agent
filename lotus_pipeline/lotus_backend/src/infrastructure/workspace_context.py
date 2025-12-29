from pathlib import Path
from typing import Optional


class WorkspaceContext:
    """
    Singleton class to manage the application's runtime workspace.
    It acts as a bridge between the relative paths stored in the database
    and the absolute file system paths on the user's machine (Mac/Win).
    """
    _instance = None
    _root_path: Optional[Path] = None

    def __new__(cls):
        """Standard Singleton implementation."""
        if cls._instance is None:
            cls._instance = super(WorkspaceContext, cls).__new__(cls)
        return cls._instance

    def initialize(self, workspace_root: str) -> None:
        """
        Initializes the workspace root directory.
        This should be called when the desktop app starts or when the user selects a project folder.

        Args:
            workspace_root (str): The absolute path string chosen by the user.
                                  e.g., 'D:/Research/Project_Alpha' or '/Users/Tony/Project_Alpha'
        """
        path = Path(workspace_root).resolve()

        # Ensure the root directory exists
        if not path.exists():
            raise FileNotFoundError(f"Workspace root does not exist: {path}")

        self._root_path = path

        # Automatically ensure the standard directory structure exists
        self._ensure_directories()

        print(f"[System] Workspace initialized at: {self._root_path}")

    def _ensure_directories(self):
        """Creates the standard sub-folders if they don't exist."""
        dirs = ["raw_data", "analysis", "cache", "objects"]
        for d in dirs:
            (self.root / d).mkdir(exist_ok=True)

    @property
    def root(self) -> Path:
        """Accessor for the root path. Raises error if not initialized."""
        if self._root_path is None:
            raise RuntimeError("WorkspaceContext is not initialized! Call initialize() first.")
        return self._root_path

    def resolve(self, relative_path_key: str) -> Path:
        """
        CORE METHOD: Converts a relative DB key to a system absolute Path.

        Args:
            relative_path_key (str): The POSIX-style relative path from the DB.
                                     e.g., 'raw_data/sample_01.h5ad'

        Returns:
            Path: A pathlib object representing the absolute path on the current OS.
        """
        # Convert string to Path object
        clean_rel = Path(relative_path_key)

        # Combine root + relative path
        abs_path = (self.root / clean_rel).resolve()

        # Security Check: Prevent Path Traversal attacks (e.g., "../../windows/system32")
        # The resolved path MUST be inside the workspace root.
        if not str(abs_path).startswith(str(self.root)):
            raise ValueError(f"Security Alert: Attempted to access path outside workspace: {relative_path_key}")

        return abs_path


# Global singleton instance to be imported by other modules
workspace_path_manager = WorkspaceContext()
