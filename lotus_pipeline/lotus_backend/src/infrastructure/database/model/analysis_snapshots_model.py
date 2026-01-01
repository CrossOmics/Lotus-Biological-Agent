from peewee import (
    AutoField,
    CharField,
    DateTimeField,
    TextField,
    ForeignKeyField
)
from datetime import datetime

from playhouse.sqlite_ext import JSONField

from .base_model import BaseModel
from .dataset_model import Dataset


class AnalysisSnapshot(BaseModel):
    """
    Analysis Snapshot table for tracking analysis history and branching.

    Represents a specific state node in the analysis pipeline (e.g., QC result,
    Clustering result). Supports a tree structure to enable 'undo/redo' and
    'branching' analysis workflows.

    Attributes:
        id: Auto-increment primary key
        snapshot_id: Unique business identifier (e.g., snap_20250101_xxxx)
        dataset_primary_id: Foreign key to the source Dataset
        branch_name: Human-readable name for this analysis branch/step
        params_json: JSON string storing Scanpy parameters used
        thumbnail_json: JSON string mapping visualization keys to file paths
        create_time: Timestamp of run start
        end_time: Timestamp of run completion
        user_notes: Manual annotations for this experiment result
        ext_info: Reserved field
    """

    # Primary Key
    id = AutoField(
        primary_key=True,
        help_text="Auto-increment primary key"
    )

    # Business Logic ID
    snapshot_id = CharField(
        max_length=128,
        unique=True,
        index=True,
        null=False,
        help_text="Unique snapshot identifier (e.g., snap_<timestamp>_<4 gits random>)"
    )

    # Foreign Keys
    # Link to the raw dataset (Root of the analysis)
    dataset_primary_id = ForeignKeyField(
        Dataset,
        backref='snapshots',
        on_delete='CASCADE',
        help_text="Reference to the raw dataset"
    )

    # Metadata
    branch_name = CharField(
        max_length=255,
        null=False,
        help_text="Name of this analysis branch or step (e.g., 'Leiden Res 0.5')"
    )

    # Scanpy Technical parameters (JSON Field)
    params_json = JSONField(
        null=True,
        help_text="Scanpy parameters (QC thresholds, resolution, etc.) in JSON format"
    )

    thumbnail_json = JSONField(
        null=True,
        help_text="Map of visualization names to relative file paths in JSON format"
    )

    # Temporal Fields
    create_time = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Experiment start timestamp"
    )

    end_time = DateTimeField(
        null=True,
        help_text="Experiment completion timestamp"
    )

    user_notes = TextField(
        null=True,
        help_text="User's evaluation or notes for this result"
    )

    ext_info = JSONField(
        null=True,
        help_text="Reserved field for additional metadata"
    )

    class Meta:
        table_name = 'analysis_snapshots'
        indexes = (
            (('dataset_primary_id', 'create_time'), False),
        )

    def __repr__(self) -> str:
        return (
            f"<Snapshot id={self.id} "
            f"code='{self.snapshot_id}' "
            f"name='{self.branch_name}'>"
        )
