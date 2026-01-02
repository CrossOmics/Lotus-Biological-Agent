from datetime import datetime
from typing import List, Optional, Dict, Any, Union
from peewee import DoesNotExist, IntegrityError

from ..model.analysis_snapshots_model import AnalysisSnapshot
from ..model.dataset_model import Dataset


class AnalysisSnapshotsDAO:
    """
    Data Access Object (DAO) for managing AnalysisSnapshot database operations.
    Encapsulates all database interactions related to analysis snapshots to ensure
    consistent data handling and abstraction from the business logic.
    """

    @staticmethod
    def create_snapshot(
            dataset_primary_id: Union[Dataset, int],
            snapshot_id: str,
            branch_name: str,
            params_json: Optional[Dict[str, Any]] = None,
            thumbnail_json: Optional[Dict[str, Any]] = None,
            user_notes: Optional[str] = None
    ) -> Optional[AnalysisSnapshot]:
        """
        Creates a new AnalysisSnapshot record in the database.

        Args:
            dataset_primary_id (int): The source dataset id this snapshot belongs to.
            snapshot_id (str): Unique business identifier (e.g., 'snap_20250101_abcd').
            branch_name (str): Human-readable name for the analysis step (e.g., 'QC Filtered').
            params_json (dict, optional): Dictionary of Scanpy parameters used.
            thumbnail_json (dict, optional): Dictionary mapping visualization keys to file paths.
            user_notes (str, optional): Initial notes from the user.

        Returns:
            AnalysisSnapshot: The created snapshot instance if successful.
            None: If the creation fails (e.g., due to duplicate snapshot_id).
        """
        try:
            snapshot = AnalysisSnapshot.create(
                dataset_primary_id=dataset_primary_id,
                snapshot_id=snapshot_id,
                branch_name=branch_name,
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

    @staticmethod
    def get_snapshot_by_id(pk_id: int) -> Optional[AnalysisSnapshot]:
        """
        Retrieves a single snapshot by its primary key ID.

        Args:
            pk_id (int): The auto-increment primary key of the snapshot.

        Returns:
            AnalysisSnapshot: The found snapshot instance.
            None: If no record is found.
        """
        try:
            return AnalysisSnapshot.get_by_id(pk_id)
        except DoesNotExist:
            print(f"[Warning] Snapshot with Primary Key {pk_id} not found.")
            return None

    @staticmethod
    def get_snapshot_by_business_id(snapshot_id: str) -> Optional[AnalysisSnapshot]:
        """
        Retrieves a single snapshot by its unique business string ID.

        Args:
            snapshot_id (str): The unique string identifier (e.g., 'snap_2025...').

        Returns:
            AnalysisSnapshot: The found snapshot instance.
            None: If no record is found.
        """
        try:
            return AnalysisSnapshot.get(AnalysisSnapshot.snapshot_id == snapshot_id)
        except DoesNotExist:
            print(f"[Warning] Snapshot with ID '{snapshot_id}' not found.")
            return None

    @staticmethod
    def get_snapshots_by_dataset(dataset_id: int) -> List[AnalysisSnapshot]:
        """
        Retrieves all snapshots associated with a specific dataset.

        Args:
            dataset_id (int): The primary key of the parent dataset.

        Returns:
            List[AnalysisSnapshot]: A list of snapshot instances ordered by creation time.
        """
        try:
            query = (AnalysisSnapshot
                     .select()
                     .where(AnalysisSnapshot.dataset_primary_id == dataset_id)
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
        """
        Updates an existing snapshot record. Only updates fields that are provided (not None).

        Args:
            pk_id (int): Primary key of the snapshot to update.
            branch_name (str, optional): New name for the branch.
            user_notes (str, optional): New notes.
            end_time (datetime, optional): Timestamp when the analysis finished.
            params_update (dict, optional): New parameters to merge or replace existing ones.

        Returns:
            bool: True if the update was successful, False otherwise.
        """
        try:
            snapshot = AnalysisSnapshot.get_by_id(pk_id)

            if branch_name is not None:
                snapshot.branch_name = branch_name

            if user_notes is not None:
                snapshot.user_notes = user_notes

            if end_time is not None:
                snapshot.end_time = end_time

            # Handle JSON field updates: merge new params into existing ones
            if params_update is not None:
                # Load existing params (default to empty dict if None)
                current_params = snapshot.params_json or {}

                # Update the dictionary (this modifies current_params in-place)
                #    - New keys are added
                #    - Existing keys are updated with new values
                current_params.update(params_update)

                # Assign back to the model field
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
        """
        Deletes a snapshot record from the database.

        Args:
            pk_id (int): Primary key of the snapshot to delete.

        Returns:
            bool: True if deletion was successful, False otherwise.
        """
        try:
            snapshot = AnalysisSnapshot.get_by_id(pk_id)
            snapshot.delete_instance()
            print(f"[Info] Snapshot ID {pk_id} deleted successfully.")
            return True
        except DoesNotExist:
            print(f"[Warning] Cannot delete: Snapshot ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to delete snapshot {pk_id}: {e}")
            return False

    @staticmethod
    def delete_snapshots_by_dataset(dataset_id: int) -> int:
        """
        Deletes all snapshots associated with a specific dataset (Batch delete).

        Args:
            dataset_id (int): Primary key of the dataset.

        Returns:
            int: The number of records deleted.
        """
        try:
            query = AnalysisSnapshot.delete().where(AnalysisSnapshot.dataset_primary_id == dataset_id)
            deleted_count = query.execute()
            print(f"[Info] Deleted {deleted_count} snapshots for dataset {dataset_id}.")
            return deleted_count
        except Exception as e:
            print(f"[Error] Failed to batch delete snapshots: {e}")
            return 0
