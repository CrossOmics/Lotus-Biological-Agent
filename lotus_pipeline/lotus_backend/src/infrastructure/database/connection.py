# infrastructure/database/connection.py
"""
Database connection management module for FastAPI.
Implements thread-local SQLite connections with proper encapsulation.
"""

import threading
from pathlib import Path
from typing import Optional

from peewee import SqliteDatabase, Proxy
from infrastructure.workspace_context import workspace_path_manager

# Global Database Proxy
database_proxy = Proxy()


class DatabaseConnection:
    """
    Thread-safe database connection manager for SQLite.
    """

    _thread_local = threading.local()

    def __init__(self, default_db_path: Optional[Path] = None):
        """
        Initialize database connection manager.

        Args:
            default_db_path: Default path to SQLite database file.
        """
        self.default_db_path = default_db_path

    def get_connection(self, db_path: Optional[Path] = None) -> SqliteDatabase:
        """
        Get or create a thread-local database connection.
        """
        path_to_use = db_path or self.default_db_path or (workspace_path_manager.root / "lotus.db")

        if not hasattr(self._thread_local, 'connection'):
            self._thread_local.connection = SqliteDatabase(
                path_to_use,
                pragmas={
                    'journal_mode': 'wal',
                    'cache_size': -1024 * 128,
                    'foreign_keys': 1,
                    'synchronous': 1
                }
            )
        elif db_path and db_path != Path(self._thread_local.connection.database):
            self.close_connection()
            return self.get_connection(db_path)

        return self._thread_local.connection

    def initialize_proxy(self, db_path: Optional[Path] = None) -> None:
        """
        Create global proxy, bind it to the real database. Invoke this function before creating tables.
        """
        db = self.get_connection(db_path)
        database_proxy.initialize(db)
        print(f"[DB] Proxy initialized with database: {db.database}")

    def close_connection(self) -> None:
        """Close and remove the current thread's database connection."""
        if hasattr(self._thread_local, 'connection'):
            if not self._thread_local.connection.is_closed():
                self._thread_local.connection.close()
            delattr(self._thread_local, 'connection')

    def connection_exists(self) -> bool:
        """Check if the current thread has an active connection."""
        return (hasattr(self._thread_local, 'connection') and
                not self._thread_local.connection.is_closed())

    def get_current_path(self) -> Optional[Path]:
        """Get the database path for the current thread's connection."""
        if hasattr(self._thread_local, 'connection'):
            return Path(self._thread_local.connection.database)
        return None


# Global Database Manager
default_db_manager: Optional[DatabaseConnection] = None


def get_default_db_manager() -> DatabaseConnection:
    """
    Lazily create default DatabaseConnection after workspace is initialized.
    """
    global default_db_manager

    if default_db_manager is None:
        default_db_manager = DatabaseConnection(
            workspace_path_manager.root / "lotus.db"
        )

    return default_db_manager


def get_db() -> SqliteDatabase:
    """
    FastAPI / service layer dependency.
    
    Returns:
        SqliteDatabase: Thread-local database connection
    """
    return get_default_db_manager().get_connection()
