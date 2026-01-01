from pathlib import Path
from typing import Optional
import threading

from infrastructure.filesystem.constants.filesystem_constants import (
    USER_PROJECT_ROOT,
    CACHE_PATH,
    ARCHIVE_PATH
)


class WorkspaceContext:
    """
    Singleton to manage the application's workspace root.

    Automatically initializes with default workspace on first access.
    Thread-safe implementation.
    """

    _instance = None
    _lock = threading.Lock()  # 线程安全
    _root_path: Optional[Path] = None
    _is_initialized: bool = False

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        pass

    # Public API

    def initialize(self, workspace_root: Optional[str] = None, force: bool = False) -> None:
        """
        Explicitly initialize or re-initialize the workspace root.

        Args:
            workspace_root: Absolute path chosen by the user, or None for default.
            force: If True, allow re-initialization even if already initialized.

        Raises:
            RuntimeError: If already initialized and force=False.
        """
        if self._is_initialized and not force:
            raise RuntimeError(
                f"WorkspaceContext already initialized at: {self._root_path}. "
                "Use force=True to re-initialize."
            )

        with self._lock:
            self._initialize_internal(workspace_root)

    @property
    def root(self) -> Path:
        """
        Get workspace root. Auto-initializes with default path if not yet initialized.
        """
        if not self._is_initialized:
            with self._lock:
                if not self._is_initialized:  # Double-check locking
                    print("[System] Auto-initializing workspace with default path...")
                    self._initialize_internal(workspace_root=None)

        return self._root_path

    @property
    def user_project_root(self) -> Path:
        """User project space root folder."""
        return self.root / USER_PROJECT_ROOT

    @property
    def cache_root(self) -> Path:
        """Cache folder path."""
        return self.root / CACHE_PATH

    @property
    def archive_root(self) -> Path:
        """Archive folder path."""
        return self.root / ARCHIVE_PATH

    def resolve(self, relative_path_key: str) -> Path:
        """
        Convert a POSIX-style relative DB path to an absolute filesystem path.

        Args:
            relative_path_key: Relative path stored in database (POSIX format).

        Returns:
            Absolute OS-specific path.

        Raises:
            ValueError: If path escapes workspace (security check).
        """
        rel_path = Path(relative_path_key)
        abs_path = (self.root / rel_path).resolve()

        # Security: ensure path is within workspace
        if not abs_path.is_relative_to(self.root):
            raise ValueError(
                f"Security violation: path escapes workspace: {relative_path_key}"
            )

        return abs_path

    def is_initialized(self) -> bool:
        """Check if workspace has been initialized."""
        return self._is_initialized

    def reset(self) -> None:
        """Reset workspace context (mainly for testing)."""
        with self._lock:
            self._root_path = None
            self._is_initialized = False
            print("[System] WorkspaceContext reset.")

    # Internal Methods
    def _initialize_internal(self, workspace_root: Optional[str]) -> None:
        """
        Internal initialization logic.

        This method is NOT thread-safe by itself - caller must hold _lock.
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
        self._is_initialized = True

        print(f"[System] Workspace initialized ({source}) at: {self._root_path}")

    def _ensure_directories(self) -> None:
        """Create essential subdirectories if they don't exist."""
        for name in (USER_PROJECT_ROOT, CACHE_PATH, ARCHIVE_PATH):
            (self._root_path / name).mkdir(parents=True, exist_ok=True)

    def _get_default_root(self) -> Path:
        """
        Resolve the default data_root used for development or testing.
        """
        current_file = Path(__file__).resolve()
        project_root = current_file.parents[3]
        default_root = project_root / "data_root"

        return default_root.resolve()


# Global Singleton Instance
workspace_path_manager = WorkspaceContext()
