from django.db.backends.postgresql_psycopg2.base import DatabaseWrapper as Psycopg2DatabaseWrapper
from django_rdkit.db.backends.postgresql_psycopg2.creation import DatabaseCreation
#from django_rdkit.db.backends.postgresql_psycopg2.introspection import DatabaseIntrospection
#from django_rdkit.db.backends.postgresql_psycopg2.operations import DatabaseOperations
from django_rdkit.db.backends.postgresql_psycopg2.schema import DatabaseSchemaEditor

class DatabaseWrapper(Psycopg2DatabaseWrapper):
    def __init__(self, *args, **kwargs):
        super(DatabaseWrapper, self).__init__(*args, **kwargs)
        self.creation = DatabaseCreation(self)
        # self.ops = DatabaseOperations(self)
        # self.introspection = DatabaseIntrospection(self)

    def schema_editor(self, *args, **kwargs):
        "Returns a new instance of this backend's SchemaEditor"
        return DatabaseSchemaEditor(self, *args, **kwargs)

