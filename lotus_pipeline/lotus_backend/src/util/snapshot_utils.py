import uuid
from datetime import datetime
from typing import Optional


def resolve_source_snapshot(dataset_id: str, snapshot_id: Optional[str], target_branch: str, snapshot_dao):
    """
    Helper to find the specific source snapshot or auto-discover the latest one.
    """
    if snapshot_id:
        snapshot = snapshot_dao.get_snapshot_by_business_id(snapshot_id)
        if not snapshot:
            raise ValueError(f"Source snapshot {snapshot_id} not found.")
        return snapshot

    # Auto-discovery
    latest_snap = snapshot_dao.get_latest_snapshot(dataset_id, branch_name=target_branch)
    if not latest_snap:
        raise ValueError(f"No preceding '{target_branch}' snapshot found for dataset {dataset_id}. Cannot proceed.")

    return latest_snap


def generate_snapshot_id() -> str:
    """Generates a unique snapshot ID."""
    timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    unique_suffix = uuid.uuid4().hex[:4]
    return f"snap_{timestamp_str}_{unique_suffix}"
