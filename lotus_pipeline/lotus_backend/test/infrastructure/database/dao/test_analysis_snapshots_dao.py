import pytest
import time
from datetime import datetime

from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO
from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO


def test_create_snapshot_success():
    """
    Test creating a snapshot linked to a valid dataset.
    Verifies that JSON parameters and foreign keys are stored correctly.
    """
    # 1. Setup Hierarchy: Project -> Dataset
    timestamp = int(time.time())
    project = ProjectMetaDAO.create_project(f"p_snap_{timestamp}", "Snapshot Parent Project")
    dataset = DatasetDAO.create_dataset(project, f"ds_snap_{timestamp}", "Snapshot Parent Dataset")

    # 2. Action: Create Snapshot
    snap_id = f"snap_base_{timestamp}"
    params = {"resolution": 0.5, "n_neighbors": 15}

    snapshot = AnalysisSnapshotsDAO.create_snapshot(
        dataset_primary_id=dataset,
        snapshot_id=snap_id,
        branch_name="Initial Clustering",
        params_json=params,
        user_notes="First run looks good"
    )

    # 3. Validation
    assert snapshot is not None
    assert snapshot.snapshot_id == snap_id
    assert snapshot.branch_name == "Initial Clustering"
    # Verify Foreign Key
    assert snapshot.dataset_primary_id == dataset
    # Verify JSON storage
    assert snapshot.params_json == params
    assert snapshot.params_json['resolution'] == 0.5
    # Verify timestamps
    assert isinstance(snapshot.create_time, datetime)


def test_get_snapshot_by_business_id():
    """
    Test retrieving a snapshot by its unique business string ID.
    """
    # 1. Setup
    timestamp = int(time.time())
    project = ProjectMetaDAO.create_project(f"p_get_snap_{timestamp}", "Get Snapshot Project")
    dataset = DatasetDAO.create_dataset(project, f"ds_get_snap_{timestamp}", "Dataset")

    target_id = f"snap_target_{timestamp}"
    AnalysisSnapshotsDAO.create_snapshot(dataset, target_id, "Target Snapshot")

    # 2. Action
    retrieved_snap = AnalysisSnapshotsDAO.get_snapshot_by_business_id(target_id)

    # 3. Validation
    assert retrieved_snap is not None
    assert retrieved_snap.snapshot_id == target_id
    assert retrieved_snap.branch_name == "Target Snapshot"
    # Verify traversal to parent dataset
    assert retrieved_snap.dataset_primary_id.dataset_id == dataset.dataset_id


def test_update_snapshot_json_merge_logic():
    """
    Test updating a snapshot, specifically verifying the JSON merge logic.
    Existing keys should remain, new keys added, and overlapping keys updated.
    """
    # 1. Setup
    timestamp = int(time.time())
    project = ProjectMetaDAO.create_project(f"p_upd_snap_{timestamp}", "Update Project")
    dataset = DatasetDAO.create_dataset(project, f"ds_upd_snap_{timestamp}", "Dataset")

    # Initial Params: {'method': 'leiden', 'resolution': 0.5}
    initial_params = {"method": "leiden", "resolution": 0.5}
    snapshot = AnalysisSnapshotsDAO.create_snapshot(
        dataset,
        f"snap_upd_{timestamp}",
        "Original Name",
        params_json=initial_params
    )

    # 2. Action: Update with new params
    # Change 'resolution' (update), Add 'min_dist' (new), Leave 'method' alone (preserve)
    new_params = {"resolution": 1.0, "min_dist": 0.3}

    success = AnalysisSnapshotsDAO.update_snapshot(
        pk_id=snapshot.id,
        branch_name="Renamed Branch",
        params_update=new_params
    )

    # 3. Validation
    assert success is True

    # Reload from DB
    updated_snap = AnalysisSnapshotsDAO.get_snapshot_by_id(snapshot.id)

    assert updated_snap.branch_name == "Renamed Branch"

    # Verify JSON Merge
    current_json = updated_snap.params_json
    assert current_json['method'] == "leiden"  # Preserved
    assert current_json['resolution'] == 1.0  # Updated
    assert current_json['min_dist'] == 0.3  # Added


def test_delete_snapshot():
    """
    Test deleting a single snapshot.
    """
    # 1. Setup
    timestamp = int(time.time())
    project = ProjectMetaDAO.create_project(f"p_del_snap_{timestamp}", "Delete Project")
    dataset = DatasetDAO.create_dataset(project, f"ds_del_snap_{timestamp}", "Dataset")

    snapshot = AnalysisSnapshotsDAO.create_snapshot(dataset, f"snap_del_{timestamp}", "To Delete")

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
    project = ProjectMetaDAO.create_project(f"p_batch_{timestamp}", "Batch Project")
    dataset = DatasetDAO.create_dataset(project, f"ds_batch_{timestamp}", "Dataset to Empty")

    id_list = []
    # Create 3 snapshots
    for i in range(3):
        result = AnalysisSnapshotsDAO.create_snapshot(dataset, f"snap_batch_{i}_{timestamp}", f"Snap {i}")
        id_list.append(result.snapshot_id)

    # Verify creation
    initial_list = AnalysisSnapshotsDAO.get_snapshots_by_dataset(dataset.id)
    returned_ids = {snap.snapshot_id for snap in initial_list}

    assert set(id_list).issubset(returned_ids)

    # 2. Action: Delete all for this dataset
    deleted_count = AnalysisSnapshotsDAO.delete_snapshots_by_dataset(dataset.id)

    # 3. Validation
    assert deleted_count == 3
    remaining = AnalysisSnapshotsDAO.get_snapshots_by_dataset(dataset.id)
    remaining_ids = {snap.snapshot_id for snap in remaining}
    assert set(id_list).isdisjoint(remaining_ids)


def test_create_duplicate_snapshot_id_fails():
    """
    Test that creating a snapshot with a duplicate 'snapshot_id' is rejected.
    """
    # 1. Setup
    timestamp = int(time.time())
    project = ProjectMetaDAO.create_project(f"p_dup_snap_{timestamp}", "Dup Project")
    dataset = DatasetDAO.create_dataset(project, f"ds_dup_snap_{timestamp}", "Dataset")

    dup_id = f"snap_unique_{timestamp}"

    # Create first
    AnalysisSnapshotsDAO.create_snapshot(dataset, dup_id, "First")

    # 2. Action: Try creating second with same ID
    result = AnalysisSnapshotsDAO.create_snapshot(dataset, dup_id, "Second")

    # 3. Validation
    assert result is None
