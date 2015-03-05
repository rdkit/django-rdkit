from django import VERSION as _django_version

from django.utils.functional import cached_property
from django.db.migrations.operations.base import Operation

if _django_version[:2] > (1, 7):
    from django.contrib.postgres.operations import CreateExtension
else:
    # this is just the same CreateExtension Operation backported
    # (actually, copied) from django 1.8 contrib.postgres, for those
    # still using django 1.7
    class CreateExtension(Operation):
        reversible = True

        def __init__(self, name):
            self.name = name

        def state_forwards(self, app_label, state):
            pass

        def database_forwards(self, app_label, schema_editor, 
                              from_state, to_state):
            schema_editor.execute("CREATE EXTENSION IF NOT EXISTS %s" % 
                                  self.name)

        def database_backwards(self, app_label, schema_editor, 
                               from_state, to_state):
            schema_editor.execute("DROP EXTENSION %s" % self.name)

        def describe(self):
            return "Creates extension %s" % self.name


class RDKitExtension(CreateExtension):

    def __init__(self):
        self.name = 'rdkit'

    def database_forwards(self, app_label, schema_editor, from_state, to_state):
        super(RDKitExtension, self).database_forwards(app_label, 
                                                      schema_editor, 
                                                      from_state, to_state)
        # Ensure the RDKit extension is fully loaded, because a subsequent
        # data migration would use the same connection
        with schema_editor.connection.cursor() as c:
            c.execute("SELECT mol_to_smiles('C'::mol)")


class GiSTIndex(Operation):

    reversible = True

    def __init__(self, model_name, name, index_name=None):
        self.model_name = model_name
        self.name = name
        self.index_name = index_name

    def state_forwards(self, app_label, state):
        pass

    def database_forwards(self, app_label, schema_editor, from_state, to_state):
        model = from_state.apps.get_model(app_label, self.model_name)
        field = model._meta.get_field(self.name)
        qn = schema_editor.quote_name

        sql_add_gist_index = "CREATE INDEX %(index)s ON %(table)s USING GIST (%(column)s)"

        index_name = (
            self.index_name if self.index_name is not None
            else '%s_%s_gist_idx' % (model._meta.db_table, field.column)
            )

        sql = sql_add_gist_index % {
            "index": qn(index_name),
            "table": qn(model._meta.db_table),
            "column": qn(field.column),
        }
        schema_editor.execute(sql)

    def database_backwards(self, app_label, schema_editor, from_state, to_state):
        model = from_state.apps.get_model(app_label, self.model_name)
        field = model._meta.get_field(self.name)
        qn = schema_editor.quote_name

        sql_remove_gist_index = "DROP INDEX %(index)s"

        index_name = (
            self.index_name if self.index_name is not None
            else '%s_%s_gist_idx' % (model._meta.db_table, field.column)
            )

        sql = sql_remove_gist_index % {
            "index": qn(index_name),
        }
        schema_editor.execute(sql)

    def describe(self):
        return "Creates GiST index on %s.%s" % (self.model_name, self.name)
   
