"""
Database connection management module.
Handles SQLite database initialization and connection lifecycle.
"""

from pathlib import Path
from peewee import SqliteDatabase
import threading

# Thread-local storage for connection management
_db_state = threading.local()


class DatabaseManager:
    """
    Singleton database manager for SQLite connections.
    Ensures only one database instance exists across the application.
    """

    _instance = None
    _lock = threading.Lock()

    def __new__(cls, db_path: Path = None):
        """
        Implement singleton pattern with thread-safe initialization.

        Args:
            db_path: Path to SQLite database file.
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, db_path: Path = None):
        """
        Initialize database connection.

        Args:
            db_path: Path to SQLite database file. Defaults to './lotus.db'
        """
        if not hasattr(self, '_initialized'):
            self.db_path = db_path or Path('./lotus.db')

            # Initialize SQLite database with customized settings
            self.db = SqliteDatabase(
                self.db_path,
                pragmas={
                    'journal_mode': 'wal',  # Write-Ahead Logging for better concurrency
                    'cache_size': -1024 * 128,  # 128MB cache
                    'foreign_keys': 1,  # Enable foreign key constraints
                    'ignore_check_constraints': 0,
                    'synchronous': 0  # Faster writes (trade-off: less durability)
                }
            )

            self._initialized = True
            print(f"[DB] Database initialized at: {self.db_path.absolute()}")

    def connect(self):
        """Open database connection."""
        if self.db.is_closed():
            self.db.connect()
            print("[DB] Database connected")

    def close(self):
        """Close database connection."""
        if not self.db.is_closed():
            self.db.close()
            print("[DB] Database closed")

    def get_connection(self) -> SqliteDatabase:
        """
        Get the current database connection.

        Returns:
            SqliteDatabase: Active database connection.
        """
        return self.db

    def __enter__(self):
        """Context manager entry: open connection."""
        self.connect()
        return self.db

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: close connection."""
        self.close()


# Global database instance (initialized on import)
db_manager = DatabaseManager()
db = db_manager.get_connection()
