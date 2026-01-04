from peewee import (
    AutoField,
    CharField,
    TextField,
    DateTimeField
)
from datetime import datetime


from .base_model import BaseModel


class ProjectMeta(BaseModel):
    """
    Project table for single-cell analysis projects.

    Stores metadata about research projects including organism,
    tissue type, and temporal information.

    Attributes:
        id: Auto-increment primary key
        project_id: Unique project identifier (UUID recommended)
        project_name: Human-readable project name displayed in UI
        create_time: ISO8601 timestamp of project creation
        update_time: ISO8601 timestamp of last modification
        description: Optional project description
        organism: Species information (Human/Mouse)
        tissue_type: Tissue origin (PBMC/Brain/Liver/etc.)
        ext_info: Reserved JSON field for extensibility
    """

    # Primary key
    id = AutoField(
        primary_key=True,
        help_text="Auto-increment primary key"
    )

    # Unique project identifier
    project_id = CharField(
        max_length=64,
        unique=True,
        index=True,
        null=False,
        help_text="Unique project identifier (p_timestamp_5 digits random number)"
    )

    # Display name
    project_name = CharField(
        max_length=255,
        null=False,
        help_text="Project name displayed in user interface"
    )

    # Project Root folder
    project_path = TextField(
        null=False,
        help_text="Project root folder"
    )

    # Temporal fields (ISO8601 format: 2025-01-30T18:20:00Z)
    create_time = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Project creation timestamp in ISO8601 format"
    )

    update_time = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Last update timestamp in ISO8601 format"
    )

    # Description
    description = TextField(
        null=True,
        help_text="Optional project description"
    )

    # Biological metadata
    organism = CharField(
        max_length=50,
        null=True,
        index=True,
        help_text="Species (e.g., Human, Mouse)"
    )

    tissue_type = CharField(
        max_length=100,
        null=True,
        index=True,
        help_text="Tissue type (e.g., PBMC, Brain, Liver)"
    )

    # Extensibility field
    ext_info = TextField(
        null=True,
        help_text="Reserved field for additional metadata (JSON format)"
    )

    class Meta:
        table_name = 'project_meta'
        indexes = (
            # Composite index for common query patterns
            (('organism', 'tissue_type'), False),
        )

    def __repr__(self) -> str:
        """String representation for debugging."""
        return (
            f"<Project id={self.id} "
            f"project_id='{self.project_id}' "
            f"name='{self.project_name}'>"
        )

    def __str__(self) -> str:
        """Human-readable string representation."""
        return f"{self.project_name} ({self.project_id})"
