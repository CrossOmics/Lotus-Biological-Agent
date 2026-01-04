from infrastructure.database.connection import get_default_db_manager
from infrastructure.database.model.project_meta_model import ProjectMeta
from infrastructure.database.model.dataset_model import Dataset
from infrastructure.database.model.analysis_snapshots_model import AnalysisSnapshot


def db_setup(ensure_schema: bool = True):
    """
    Handles the full database initialization workflow:
    connecting, binding the proxy, and creating tables.
    """
    # 1. Acquire the database manager
    manager = get_default_db_manager()

    # 2. Bind the proxy to the actual database instance
    manager.initialize_proxy()

    # 3. Get the database connection and ensure it is open
    db = manager.get_connection()
    if db.is_closed():
        db.connect()

    # 4. Centralized table creation logic
    if ensure_schema:
        # Explicitly list all models that require table creation
        models_to_create = [
            ProjectMeta,
            Dataset,
            AnalysisSnapshot
        ]
        db.create_tables(models_to_create, safe=True)
        print(f"[DB Setup] Verified tables: {[m.__name__ for m in models_to_create]}")

    return db
