import pytest
import time
from datetime import datetime

from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO
from infrastructure.database.dao.dataset_dao import DatasetDAO


def test_create_dataset_success():
    """
    Test creating a dataset linked to a valid project.
    Verifies that all fields (including Foreign Keys and JSON) are persisted correctly.
    """
    # 1. Setup: Create a parent project first (Foreign Key requirement)
    # Using a timestamp to ensure uniqueness in the persistent DB
    timestamp = int(time.time())
    project_id = f"p_ds_test_{timestamp}"
    project = ProjectMetaDAO.create_project(project_id, "Parent Project for Dataset")

    assert project is not None, "Prerequisite project creation failed"

    # 2. Action: Create the dataset
    dataset_business_id = f"ds_{timestamp}_raw"
    metadata = {"file_size": 102400, "file_type": "h5ad", "source": "lab_upload"}

    dataset = DatasetDAO.create_dataset(
        project=project,
        dataset_id=dataset_business_id,
        dataset_name="PBMC Raw Data",
        ext_info=metadata
    )

    # 3. Validation: Check if the returned object matches input
    assert dataset is not None
    assert dataset.dataset_id == dataset_business_id
    assert dataset.dataset_name == "PBMC Raw Data"
    # Verify Foreign Key relationship
    assert dataset.project_primary_id == project
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
    project = ProjectMetaDAO.create_project(f"p_get_{timestamp}", "Get Test Project")
    target_id = f"ds_target_{timestamp}"

    DatasetDAO.create_dataset(project, target_id, "Target Dataset")

    # 2. Action
    retrieved_ds = DatasetDAO.get_dataset_by_business_id(target_id)

    # 3. Validation
    assert retrieved_ds is not None
    assert retrieved_ds.dataset_id == target_id
    assert retrieved_ds.dataset_name == "Target Dataset"
    # Verify we can access parent project details through the relationship
    assert retrieved_ds.project_primary_id.project_id == project.project_id


def test_update_dataset_info():
    """
    Test updating a dataset's name and metadata (JSON).
    """
    # 1. Setup
    timestamp = int(time.time())
    project = ProjectMetaDAO.create_project(f"p_upd_{timestamp}", "Update Test Project")
    dataset = DatasetDAO.create_dataset(
        project,
        f"ds_upd_{timestamp}",
        "Original Name",
        ext_info={"status": "raw"}
    )

    # 2. Action: Update name and replace metadata
    new_metadata = {"status": "processed", "obs_count": 5000}
    success = DatasetDAO.update_dataset(
        pk_id=dataset.id,
        dataset_name="Renamed Dataset",
        ext_info_update=new_metadata
    )

    # 3. Validation
    assert success is True

    # Reload from DB to verify persistence
    updated_ds = DatasetDAO.get_dataset_by_id(dataset.id)
    assert updated_ds.dataset_name == "Renamed Dataset"
    assert updated_ds.ext_info == new_metadata
    assert updated_ds.ext_info['status'] == "processed"


def test_delete_dataset_lifecycle():
    """
    Test the full lifecycle: Create -> Exists -> Delete -> Gone.
    """
    # 1. Setup
    timestamp = int(time.time())
    project = ProjectMetaDAO.create_project(f"p_del_{timestamp}", "Delete Test Project")
    dataset = DatasetDAO.create_dataset(project, f"ds_del_{timestamp}", "To Delete")

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
    project = ProjectMetaDAO.create_project(f"p_dup_{timestamp}", "Dup Test Project")
    dup_id = f"ds_unique"

    # Create first instance
    DatasetDAO.create_dataset(project, dup_id, "First Instance")

    # 2. Action: Try creating second instance with same ID
    result2 = DatasetDAO.create_dataset(project, dup_id, "Second Instance")

    # 3. Validation: Should return None due to IntegrityError
    assert result2 is None