import pytest
import shutil
from pathlib import Path

# Import your infrastructure components
from infrastructure.filesystem.logging import setup_logging
from infrastructure.workspace_context import workspace_path_manager


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

    yield

    # (Optional): Clean up after all tests are done
    # print("\n[Test Config] Tearing down Test Environment...")
    # shutil.rmtree(TEST_ROOT_DIR)
