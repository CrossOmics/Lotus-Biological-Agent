from datetime import datetime
from typing import List, Optional, Dict, Any, Union
from peewee import DoesNotExist, IntegrityError

from ..model.analysis_snapshots_model import AnalysisSnapshot
from ..model.dataset_model import Dataset


class AnalysisSnapshotsDAO:
    def __init__(self):
        pass
    """
    Data Access Object (DAO) for managing AnalysisSnapshot database operations.
    """
    @staticmethod
    def get_latest_snapshot(dataset_id: str, branch_name: str) -> Optional[AnalysisSnapshot]:
        """
        Retrieves the most recently created snapshot for a specific dataset and branch.

        Args:
            dataset_id: The business ID of the dataset.
            branch_name: The specific branch/step name to filter by (e.g., 'QC Filtered').

        Returns:
            The latest AnalysisSnapshot object, or None if no match found.
        """
        try:
            return (AnalysisSnapshot
                    .select()
                    .where(
                (AnalysisSnapshot.dataset_id == dataset_id) &
                (AnalysisSnapshot.branch_name == branch_name)
            )
                    .order_by(AnalysisSnapshot.create_time.desc())  # Newest first
                    .first())
        except Exception as e:
            print(f"[DB Error] Error fetching latest snapshot: {e}")
            return None

    @staticmethod
    def create_snapshot(
            dataset_id: Union[Dataset, str],
            snapshot_id: str,
            branch_name: str,
            snapshot_path: str,
            params_json: Optional[Dict[str, Any]] = None,
            thumbnail_json: Optional[Dict[str, Any]] = None,
            user_notes: Optional[str] = None
    ) -> Optional[AnalysisSnapshot]:
        """
        Creates a new AnalysisSnapshot record in the database.
        """
        try:
            snapshot = AnalysisSnapshot.create(
                dataset_id=dataset_id,
                snapshot_id=snapshot_id,
                branch_name=branch_name,
                snapshot_path=snapshot_path,
                params_json=params_json or {},
                thumbnail_json=thumbnail_json or {},
                user_notes=user_notes,
                create_time=datetime.utcnow()
            )
            return snapshot
        except IntegrityError as e:
            print(f"[Error] Failed to create snapshot '{snapshot_id}': {e}")
            return None
        except Exception as e:
            print(f"[Error] Unexpected error creating snapshot: {e}")
            return None

    # ... (Rest of the methods: get_snapshot_by_id, get_snapshot_by_business_id, etc. remain unchanged)

    @staticmethod
    def get_snapshot_by_id(pk_id: int) -> Optional[AnalysisSnapshot]:
        try:
            return AnalysisSnapshot.get_by_id(pk_id)
        except DoesNotExist:
            return None

    @staticmethod
    def get_snapshot_by_business_id(snapshot_id: str) -> Optional[AnalysisSnapshot]:
        try:
            return AnalysisSnapshot.get(AnalysisSnapshot.snapshot_id == snapshot_id)
        except DoesNotExist:
            return None

    @staticmethod
    def get_snapshots_by_dataset(dataset_id: str) -> List[AnalysisSnapshot]:
        try:
            query = (AnalysisSnapshot
                     .select()
                     .where(AnalysisSnapshot.dataset_id == dataset_id)
                     .order_by(AnalysisSnapshot.create_time.desc()))
            return list(query)
        except Exception as e:
            print(f"[Error] Failed to retrieve snapshots for dataset {dataset_id}: {e}")
            return []

    @staticmethod
    def update_snapshot(
            pk_id: int,
            branch_name: Optional[str] = None,
            user_notes: Optional[str] = None,
            end_time: Optional[datetime] = None,
            params_update: Optional[Dict[str, Any]] = None
    ) -> bool:
        try:
            snapshot = AnalysisSnapshot.get_by_id(pk_id)

            if branch_name is not None:
                snapshot.branch_name = branch_name

            if user_notes is not None:
                snapshot.user_notes = user_notes

            if end_time is not None:
                snapshot.end_time = end_time

            if params_update is not None:
                current_params = snapshot.params_json or {}
                current_params.update(params_update)
                snapshot.params_json = current_params

            snapshot.save()
            return True
        except DoesNotExist:
            print(f"[Error] Cannot update: Snapshot ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to update snapshot {pk_id}: {e}")
            return False

    @staticmethod
    def delete_snapshot(pk_id: int) -> bool:
        try:
            snapshot = AnalysisSnapshot.get_by_id(pk_id)
            snapshot.delete_instance()
            return True
        except DoesNotExist:
            return False
        except Exception as e:
            print(f"[Error] Failed to delete snapshot {pk_id}: {e}")
            return False

    @staticmethod
    def delete_snapshots_by_dataset(dataset_id: str) -> int:
        try:
            query = AnalysisSnapshot.delete().where(AnalysisSnapshot.dataset_id == dataset_id)
            deleted_count = query.execute()
            return deleted_count
        except Exception as e:
            print(f"[Error] Failed to batch delete snapshots: {e}")
            return 0


