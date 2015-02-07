from __future__ import unicode_literals

from django.db.backends.postgresql_psycopg2.base import DatabaseWrapper as Psycopg2DatabaseWrapper

#from .creation import DatabaseCreation
#from .introspection import DatabaseIntrospection
#from .operations import DatabaseOperations
from .schema import DatabaseSchemaEditor

class DatabaseWrapper(Psycopg2DatabaseWrapper):

    SchemaEditorClass = DatabaseSchemaEditor

    def prepare_database(self):
        super(DatabaseWrapper, self).prepare_database()
        with self.cursor() as cursor:
            cursor.execute("CREATE EXTENSION IF NOT EXISTS rdkit")

