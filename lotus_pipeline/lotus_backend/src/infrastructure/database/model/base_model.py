"""
ORM models for single-cell analysis project management.
Defines the Project table schema and business logic.
"""

from peewee import (
    Model,

)

from ..connection import db


class BaseModel(Model):
    """
    Base model for all database tables.
    Provides common configuration and utility methods.
    """

    class Meta:
        database = db
        legacy_table_names = False

