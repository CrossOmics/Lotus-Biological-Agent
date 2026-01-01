from pathlib import Path

from infrastructure.database.connection import DatabaseConnection, database_proxy
from infrastructure.database.model.dataset_model import Dataset
from infrastructure.database.model.project_meta_model import ProjectMeta
from infrastructure.database.model.analysis_snapshots_model import AnalysisSnapshot
from infrastructure.workspace_context import workspace_path_manager


def initialize_database(db_path: Path):
    """
    Initialize the database schema for the current workspace.
    """

    # 1. Create database manager
    db_manager = DatabaseConnection(db_path)

    # 2. Bind the proxy to the actual database instance
    db_manager.initialize_proxy()

    # 3. Acquire database connection
    db = db_manager.get_connection()

    # 4. List of models to create
    models = [
        ProjectMeta,
        Dataset,
        AnalysisSnapshot
    ]

    try:
        print(f"[DB] Initializing database at: {db_path}")

        # 5. Ensure the connection is open
        if db.is_closed():
            db.connect(reuse_if_open=True)

        # 6. Create tables if they do not exist
        db.create_tables(models, safe=True)

        print("[DB] Tables created successfully:")
        for model in models:
            print(f"  - {model._meta.table_name}")

    except Exception as e:
        print(f"[DB] Error creating tables: {e}")
        raise

    finally:
        # 7. Close the connection
        if not db.is_closed():
            db.close()
            print("[DB] Database connection closed.")


if __name__ == "__main__":
    # Resolve workspace root directory
    workspace_root = workspace_path_manager.root
    db_file = workspace_root / "lotus.db"

    # Initialize database schema
    initialize_database(db_file)
