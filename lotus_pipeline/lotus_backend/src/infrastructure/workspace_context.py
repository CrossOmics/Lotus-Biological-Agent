from pathlib import Path
from typing import Optional


class WorkspaceContext:
    """
    Singleton to manage the application's workspace root.

    It maps POSIX-style relative paths stored in the database
    to absolute OS-specific paths (Windows / macOS).
    """

    _instance = None
    _root_path: Optional[Path] = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    # Public API
    def initialize(self, workspace_root: Optional[str] = None) -> None:
        """
        Initialize the workspace root.

        Priority:
        1. User-specified workspace_root
        2. Default development/test data_root

        Args:
            workspace_root: Absolute path chosen by the user, or None.
        """
        if workspace_root is None:
            path = self._get_default_root()
            path.mkdir(parents=True, exist_ok=True)
            source = "DEFAULT"
        else:
            path = Path(workspace_root).expanduser().resolve()
            source = "USER"

        if not path.exists():
            raise FileNotFoundError(f"Workspace root does not exist: {path}")

        self._root_path = path
        self._ensure_directories()

        print(f"[System] Workspace initialized ({source}) at: {self._root_path}")

    @property
    def root(self) -> Path:
        if self._root_path is None:
            raise RuntimeError("WorkspaceContext not initialized. Call initialize() first.")
        return self._root_path

    def resolve(self, relative_path_key: str) -> Path:
        """
        Convert a POSIX-style relative DB path to an absolute filesystem path.
        """
        rel_path = Path(relative_path_key)
        abs_path = (self.root / rel_path).resolve()

        # Security: ensure path is within workspace
        if not abs_path.is_relative_to(self.root):
            raise ValueError(
                f"Security violation: path escapes workspace: {relative_path_key}"
            )

        return abs_path

    def _ensure_directories(self) -> None:
        for name in ("raw_data", "analysis", "cache", "objects"):
            (self.root / name).mkdir(exist_ok=True)

    def _get_default_root(self) -> Path:
        """
        Resolve the default data_root used for development or testing.
        """
        current_file = Path(__file__).resolve()
        project_root = current_file.parents[3]
        default_root = project_root / "data_root"

        return default_root.resolve()


# Global singleton instance
workspace_path_manager = WorkspaceContext()
