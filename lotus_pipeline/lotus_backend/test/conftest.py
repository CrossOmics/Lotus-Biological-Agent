import pytest

from infrastructure.filesystem.logging import setup_logging
from infrastructure.workspace_context import workspace_path_manager
from infrastructure.database.connection import get_default_db_manager


@pytest.fixture(scope="session", autouse=True)
def initialize_test_environment():
    """
    This fixture runs EXACTLY ONCE per test session (before all tests).
    It sets up the necessary environment (Context + Logging) for the entire suite.

    pytest applies it automatically.
    """

    print("\n[Test Config] Setting up Global Test Environment...")

    # Initialize the Workspace Context
    workspace_path_manager.initialize()

    # Initialize Logging
    setup_logging(debug_mode=True)

    print(f"[Test Config] Context and Logging initialized at: {workspace_path_manager.root}")

    # Initialize the database connection
    # Retrieve the global DB manager (default points to lotus.db)
    db_manager = get_default_db_manager()
    # Initialize the proxy and bind it to the actual lotus.db
    db_manager.initialize_proxy()
    # Get the database connection object
    db = db_manager.get_connection()
    # Connect and ensure the schema exists
    if db.is_closed():
        db.connect()

    # Run the unit tests
    yield

    db.close()
    print("\n[Dev Info] Test session finished: connection closed, data preserved.")

    # (Optional): Clean up after all tests are done
    # print("\n[Test Config] Tearing down Test Environment...")
    # shutil.rmtree(TEST_ROOT_DIR)
