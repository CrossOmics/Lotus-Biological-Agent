from peewee import Model
from infrastructure.database.connection import database_proxy


class BaseModel(Model):
    """
    Base model for all database tables.

    All models inherit from this class to share the same database connection.
    """

    class Meta:
        database = database_proxy
