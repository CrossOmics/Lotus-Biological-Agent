import time
from datetime import datetime
import pytest
import os
from pathlib import Path

# Infrastructure Imports
from infrastructure.database.connection import get_default_db_manager
from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO
from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO

# Models
from infrastructure.database.model.project_meta_model import ProjectMeta
from infrastructure.database.model.dataset_model import Dataset
from infrastructure.database.model.analysis_snapshots_model import AnalysisSnapshot

TEST_DB_PATH = Path("./test_snapshots.db")


@pytest.fixture(scope="function", autouse=True)
def setup_database():
    db_manager = get_default_db_manager()
    db_manager.initialize_proxy(TEST_DB_PATH)

    test_db = db_manager.get_connection(TEST_DB_PATH)
    if test_db.is_closed():
        test_db.connect()

    test_db.create_tables([ProjectMeta, Dataset, AnalysisSnapshot], safe=True)

    yield

    test_db.close()
    if TEST_DB_PATH.exists():
        try:
            os.remove(TEST_DB_PATH)
        except PermissionError:
            pass


def test_create_snapshot_success():
    """
    Test creating a snapshot linked to a valid dataset.
    Verifies that JSON parameters, foreign keys, and the new snapshot_path are stored correctly.
    """
    # 1. Setup Hierarchy
    timestamp = int(time.time())

    project = ProjectMetaDAO.create_project(
        project_id=f"p_snap_{timestamp}",
        project_name="Snapshot Parent Project",
        project_path=f"projects/p_snap_{timestamp}"
    )

    dataset = DatasetDAO.create_dataset(
        project=project,
        dataset_id=f"ds_snap_{timestamp}",
        dataset_name="Snapshot Parent Dataset",
        dataset_path=f"raw_data/ds_snap_{timestamp}.h5ad"
    )

    # 2. Action: Create Snapshot
    snap_id = f"snap_base_{timestamp}"
    params = {"resolution": 0.5, "n_neighbors": 15}
    dummy_path = f"snapshots/{snap_id}.h5ad"

    snapshot = AnalysisSnapshotsDAO.create_snapshot(
        dataset_id=dataset,
        snapshot_id=snap_id,
        branch_name="Initial Clustering",
        snapshot_path=dummy_path,
        params_json=params,
        user_notes="First run looks good"
    )

    # 3. Validation
    assert snapshot is not None
    assert snapshot.snapshot_id == snap_id
    assert snapshot.snapshot_path == dummy_path  # <--- Validate path
    assert snapshot.branch_name == "Initial Clustering"
    assert snapshot.dataset_id == dataset
    assert snapshot.params_json == params


def test_get_snapshot_by_business_id():
    """
    Test retrieving a snapshot by its unique business string ID.
    """
    # 1. Setup
    timestamp = int(time.time())

    project = ProjectMetaDAO.create_project(
        f"p_get_snap_{timestamp}",
        "Get Snapshot Project",
        f"projects/p_get_{timestamp}"
    )
    dataset = DatasetDAO.create_dataset(
        project,
        f"ds_get_snap_{timestamp}",
        "Dataset",
        f"raw_data/ds_get_{timestamp}.h5ad"
    )

    target_id = f"snap_target_{timestamp}"

    AnalysisSnapshotsDAO.create_snapshot(
        dataset,
        target_id,
        "Target Snapshot",
        snapshot_path="snapshots/target.h5ad"
    )

    # 2. Action
    retrieved_snap = AnalysisSnapshotsDAO.get_snapshot_by_business_id(target_id)

    # 3. Validation
    assert retrieved_snap is not None
    assert retrieved_snap.snapshot_id == target_id
    assert retrieved_snap.snapshot_path == "snapshots/target.h5ad"  # <--- Validate path
    assert retrieved_snap.dataset_id.dataset_id == dataset.dataset_id


def test_update_snapshot_json_merge_logic():
    """
    Test updating a snapshot, specifically verifying the JSON merge logic.
    """
    # 1. Setup
    timestamp = int(time.time())

    project = ProjectMetaDAO.create_project(
        f"p_upd_snap_{timestamp}",
        "Update Project",
        f"projects/p_upd_{timestamp}"
    )
    dataset = DatasetDAO.create_dataset(
        project,
        f"ds_upd_snap_{timestamp}",
        "Dataset",
        f"raw_data/ds_upd_{timestamp}.h5ad"
    )

    initial_params = {"method": "leiden", "resolution": 0.5}

    snapshot = AnalysisSnapshotsDAO.create_snapshot(
        dataset,
        f"snap_upd_{timestamp}",
        "Original Name",
        snapshot_path="snapshots/original.h5ad",
        params_json=initial_params
    )

    # 2. Action
    new_params = {"resolution": 1.0, "min_dist": 0.3}

    success = AnalysisSnapshotsDAO.update_snapshot(
        pk_id=snapshot.id,
        branch_name="Renamed Branch",
        params_update=new_params
    )

    # 3. Validation
    assert success is True
    updated_snap = AnalysisSnapshotsDAO.get_snapshot_by_id(snapshot.id)
    assert updated_snap.branch_name == "Renamed Branch"

    current_json = updated_snap.params_json
    assert current_json['method'] == "leiden"
    assert current_json['resolution'] == 1.0
    assert current_json['min_dist'] == 0.3


def test_delete_snapshot():
    """
    Test deleting a single snapshot.
    """
    # 1. Setup
    timestamp = int(time.time())

    project = ProjectMetaDAO.create_project(
        f"p_del_snap_{timestamp}",
        "Delete Project",
        f"projects/p_del_{timestamp}"
    )
    dataset = DatasetDAO.create_dataset(
        project,
        f"ds_del_snap_{timestamp}",
        "Dataset",
        f"raw_data/ds_del_{timestamp}.h5ad"
    )

    snapshot = AnalysisSnapshotsDAO.create_snapshot(
        dataset,
        f"snap_del_{timestamp}",
        "To Delete",
        snapshot_path="snapshots/del.h5ad"
    )

    assert AnalysisSnapshotsDAO.get_snapshot_by_id(snapshot.id) is not None

    # 2. Action
    success = AnalysisSnapshotsDAO.delete_snapshot(snapshot.id)

    # 3. Validation
    assert success is True
    assert AnalysisSnapshotsDAO.get_snapshot_by_id(snapshot.id) is None


def test_delete_snapshots_by_dataset_batch():
    """
    Test batch deletion of all snapshots belonging to a specific dataset.
    """
    # 1. Setup
    timestamp = int(time.time())

    project = ProjectMetaDAO.create_project(
        f"p_batch_{timestamp}",
        "Batch Project",
        f"projects/p_batch_{timestamp}"
    )
    dataset = DatasetDAO.create_dataset(
        project,
        f"ds_batch_{timestamp}",
        "Dataset to Empty",
        f"raw_data/ds_batch_{timestamp}.h5ad"
    )

    id_list = []
    # Create 3 snapshots
    for i in range(3):
        result = AnalysisSnapshotsDAO.create_snapshot(
            dataset,
            f"snap_batch_{i}_{timestamp}",
            f"Snap {i}",
            snapshot_path=f"snapshots/batch_{i}.h5ad"
        )
        id_list.append(result.snapshot_id)

    # Verify creation
    initial_list = AnalysisSnapshotsDAO.get_snapshots_by_dataset(dataset.dataset_id)
    assert len(initial_list) == 3

    # 2. Action
    deleted_count = AnalysisSnapshotsDAO.delete_snapshots_by_dataset(dataset.dataset_id)

    # 3. Validation
    assert deleted_count == 3
    remaining = AnalysisSnapshotsDAO.get_snapshots_by_dataset(dataset.dataset_id)
    assert len(remaining) == 0


def test_create_duplicate_snapshot_id_fails():
    """
    Test that creating a snapshot with a duplicate 'snapshot_id' is rejected.
    """
    # 1. Setup
    timestamp = int(time.time())

    project = ProjectMetaDAO.create_project(
        f"p_dup_snap_{timestamp}",
        "Dup Project",
        f"projects/p_dup_{timestamp}"
    )
    dataset = DatasetDAO.create_dataset(
        project,
        f"ds_dup_snap_{timestamp}",
        "Dataset",
        f"raw_data/ds_dup_{timestamp}.h5ad"
    )

    dup_id = f"snap_unique_{timestamp}"

    # Create first
    AnalysisSnapshotsDAO.create_snapshot(
        dataset,
        dup_id,
        "First",
        snapshot_path="snapshots/1.h5ad"
    )

    # 2. Action: Try creating second with same ID
    result = AnalysisSnapshotsDAO.create_snapshot(
        dataset,
        dup_id,
        "Second",
        snapshot_path="snapshots/2.h5ad"
    )

    # 3. Validation
    assert result is None