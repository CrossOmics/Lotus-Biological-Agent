from peewee import (
    AutoField,
    CharField,
    DateTimeField,
    TextField,
    ForeignKeyField
)
from datetime import datetime
from .base_model import BaseModel
from .project_meta_model import ProjectMeta
from playhouse.sqlite_ext import JSONField


class Dataset(BaseModel):
    """
    Dataset table for managing raw single-cell data files.

    Links uploaded data files to projects and tracks their physical
    storage location and basic biological statistics.

    Attributes:
        id: Auto-increment primary key
        project: Foreign key to ProjectMeta
        name: User-defined display name for the dataset
        relative_path: Path to .h5ad file relative to workspace root (raw_data/...)
        file_size_bytes: Size of the stored file in bytes
        n_obs: Number of observations (cells)
        n_vars: Number of variables (genes)
        created_at: ISO8601 timestamp of upload/ingestion
        ext_info: Reserved JSON field for extensibility
    """

    # Primary Key
    id = AutoField(
        primary_key=True,
        help_text="Auto-increment primary key"
    )

    # dataset business id
    dataset_id = CharField(
        max_length=128,
        null=False,
        help_text="dataset business id, naming rule: dataset_<timestamp>_<original_name>_<4_digits_random_suffix>"
    )
    # Foreign Key links to project id
    project_primary_id = ForeignKeyField(
        ProjectMeta,
        backref='datasets',
        on_delete='CASCADE',
        help_text="Reference to the parent project"
    )

    # Display Name
    dataset_name = CharField(
        max_length=255,
        null=False,
        help_text="User-defined dataset display name"
    )

    # Temporal Fields
    uploaded_time = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Upload timestamp"
    )

    # Extensibility
    ext_info = JSONField(
        null=True,
        help_text="Reserved field for additional metadata (JSON format)"
    )

    class Meta:
        table_name = 'dataset'
        indexes = (
            # Index for quick lookup by project
            (('project_primary_id',), False),
        )

    def __repr__(self) -> str:
        return (
            f"<Dataset id={self.id} "
            f"name='{self.name}' "
            f"project_id={self.project_id}>"
        )
