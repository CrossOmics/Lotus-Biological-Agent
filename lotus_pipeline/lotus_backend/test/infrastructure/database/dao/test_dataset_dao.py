import pytest
import os
import time
from pathlib import Path
from datetime import datetime

# Infrastructure imports
from infrastructure.database.connection import get_default_db_manager
from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO
from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.database.model.project_meta_model import ProjectMeta
from infrastructure.database.model.dataset_model import Dataset

# Define a specific path for this test suite
TEST_DB_PATH = Path("./test_dataset_cicd.db")


# --- Fixture for CI/CD Style Isolation ---
@pytest.fixture(scope="function", autouse=True)
def setup_database():
    """
    CI/CD Style Setup:
    1. Creates a fresh, isolated database file for EACH test function.
    2. Creates tables.
    3. Runs the test.
    4. Immediately drops tables and deletes the file to clean up.
    """
    # 1. Initialize Proxy to point to our specific test file
    db_manager = get_default_db_manager()
    db_manager.initialize_proxy(TEST_DB_PATH)

    # 2. Get connection and ensure it's open
    test_db = db_manager.get_connection(TEST_DB_PATH)
    if test_db.is_closed():
        test_db.connect()

    # 3. Create tables (Fresh state)
    # Note: We must create BOTH tables to satisfy Foreign Key constraints
    test_db.create_tables([ProjectMeta, Dataset], safe=True)

    yield  # Run the test case

    # 4. Teardown (Cleanup)
    # Close connection first
    test_db.close()

    # Remove the physical file to ensure no data persistence between tests
    if TEST_DB_PATH.exists():
        try:
            os.remove(TEST_DB_PATH)
        except PermissionError:
            pass


# --- Test Cases ---

def test_create_dataset_success():
    """
    Test creating a dataset linked to a valid project.
    Verifies that all fields (including Foreign Keys, JSON, and Path) are persisted correctly.
    """
    # 1. Setup: Create a parent project first (Foreign Key requirement)
    timestamp = int(time.time())
    project_id = f"p_ds_test_{timestamp}"

    project = ProjectMetaDAO.create_project(
        project_id=project_id,
        project_name="Parent Project for Dataset",
        project_path=f"/tmp/{project_id}"
    )

    assert project is not None, "Prerequisite project creation failed"

    # 2. Action: Create the dataset with the required path
    dataset_business_id = f"ds_{timestamp}_raw"
    metadata = {"file_size": 102400, "file_type": "h5ad", "source": "lab_upload"}
    # [Update]: Added dataset_path
    relative_file_path = "raw_data/pbmc_test.h5ad"

    dataset = DatasetDAO.create_dataset(
        project=project,
        dataset_id=dataset_business_id,
        dataset_name="PBMC Raw Data",
        dataset_path=relative_file_path,  # <--- Added required arg
        ext_info=metadata
    )

    # 3. Validation: Check if the returned object matches input
    assert dataset is not None
    assert dataset.dataset_id == dataset_business_id
    assert dataset.dataset_name == "PBMC Raw Data"
    assert dataset.dataset_path == relative_file_path  # <--- Verify path persistence

    # Verify Foreign Key relationship
    assert dataset.project_id == project

    # Verify JSON storage
    assert dataset.ext_info == metadata
    assert dataset.ext_info['file_type'] == "h5ad"
    # Verify timestamp generation
    assert isinstance(dataset.uploaded_time, datetime)


def test_get_dataset_by_business_id():
    """
    Test retrieving a dataset using its unique business string ID.
    """
    # 1. Setup
    timestamp = int(time.time())
    project_id = f"p_get_{timestamp}"

    project = ProjectMetaDAO.create_project(
        project_id=project_id,
        project_name="Get Test Project",
        project_path=f"/tmp/{project_id}"
    )

    target_id = f"ds_target_{timestamp}"
    # [Update]: Added dataset_path
    DatasetDAO.create_dataset(
        project,
        target_id,
        "Target Dataset",
        dataset_path="raw_data/target.h5ad"
    )

    # 2. Action
    retrieved_ds = DatasetDAO.get_dataset_by_business_id(target_id)

    # 3. Validation
    assert retrieved_ds is not None
    assert retrieved_ds.dataset_id == target_id
    assert retrieved_ds.dataset_name == "Target Dataset"
    assert retrieved_ds.dataset_path == "raw_data/target.h5ad"  # <--- Verify path

    # Verify we can access parent project details through the relationship
    assert retrieved_ds.project_id.project_id == project.project_id


def test_update_dataset_info():
    """
    Test updating a dataset's name, path, and metadata (JSON).
    """
    # 1. Setup
    timestamp = int(time.time())
    project_id = f"p_upd_{timestamp}"

    project = ProjectMetaDAO.create_project(
        project_id=project_id,
        project_name="Update Test Project",
        project_path=f"/tmp/{project_id}"
    )

    # [Update]: Added initial dataset_path
    dataset = DatasetDAO.create_dataset(
        project,
        f"ds_upd_{timestamp}",
        "Original Name",
        dataset_path="raw_data/original.h5ad",
        ext_info={"status": "raw"}
    )

    # 2. Action: Update name, path, and replace metadata
    new_metadata = {"status": "processed", "obs_count": 5000}

    # [Update]: Testing path update capability as well
    success = DatasetDAO.update_dataset(
        pk_id=dataset.id,
        dataset_name="Renamed Dataset",
        dataset_path="analysis/moved_dataset.h5ad",  # <--- Updating path
        ext_info_update=new_metadata
    )

    # 3. Validation
    assert success is True

    # Reload from DB to verify persistence
    updated_ds = DatasetDAO.get_dataset_by_id(dataset.id)
    assert updated_ds.dataset_name == "Renamed Dataset"
    assert updated_ds.dataset_path == "analysis/moved_dataset.h5ad"  # <--- Verify update
    assert updated_ds.ext_info == new_metadata
    assert updated_ds.ext_info['status'] == "processed"


def test_delete_dataset_lifecycle():
    """
    Test the full lifecycle: Create -> Exists -> Delete -> Gone.
    """
    # 1. Setup
    timestamp = int(time.time())
    project_id = f"p_del_{timestamp}"

    project = ProjectMetaDAO.create_project(
        project_id=project_id,
        project_name="Delete Test Project",
        project_path=f"/tmp/{project_id}"
    )

    # [Update]: Added dataset_path
    dataset = DatasetDAO.create_dataset(
        project,
        f"ds_del_{timestamp}",
        "To Delete",
        dataset_path="raw_data/delete_me.h5ad"
    )

    # Confirm existence
    assert DatasetDAO.get_dataset_by_id(dataset.id) is not None

    # 2. Action: Delete
    success = DatasetDAO.delete_dataset(dataset.id)

    # 3. Validation
    assert success is True
    # Verify retrieval fails
    assert DatasetDAO.get_dataset_by_id(dataset.id) is None
    # Verify business ID retrieval fails
    assert DatasetDAO.get_dataset_by_business_id(dataset.dataset_id) is None


def test_create_duplicate_dataset_id_fails():
    """
    Test that creating a dataset with a duplicate 'dataset_id' is rejected by the DB constraints.
    """
    # 1. Setup
    timestamp = int(time.time())
    project_id = f"p_dup_{timestamp}"

    project = ProjectMetaDAO.create_project(
        project_id=project_id,
        project_name="Dup Test Project",
        project_path=f"/tmp/{project_id}"
    )

    dup_id = f"ds_unique_{timestamp}"

    # Create first instance
    # [Update]: Added dataset_path
    DatasetDAO.create_dataset(
        project,
        dup_id,
        "First Instance",
        dataset_path="raw_data/1.h5ad"
    )

    # 2. Action: Try creating second instance with same ID
    # [Update]: Added dataset_path
    result2 = DatasetDAO.create_dataset(
        project,
        dup_id,
        "Second Instance",
        dataset_path="raw_data/2.h5ad"
    )

    # 3. Validation: Should return None due to IntegrityError
    assert result2 is None
