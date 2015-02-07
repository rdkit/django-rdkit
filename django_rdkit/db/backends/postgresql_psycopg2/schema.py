from __future__ import unicode_literals

from django.db.backends.postgresql_psycopg2.schema import DatabaseSchemaEditor as Psycopg2SchemaEditor

class DatabaseSchemaEditor(Psycopg2SchemaEditor):

    sql_add_chem_index = "CREATE INDEX %(index)s ON %(table)s USING GIST (%(column)s)"

    def __init__(self, *args, **kwargs):
        super(DatabaseSchemaEditor, self).__init__(*args, **kwargs)
        self.chem_sql = []

    def column_sql(self, model, field, include_default=False):
        from django_rdkit.db.models.fields import ChemField

        column_sql = super(DatabaseSchemaEditor, self).column_sql(model, field, include_default)

        if not isinstance(field, ChemField):
            return column_sql

        if field.chem_index:
            self.chem_sql.append(
                self.sql_add_chem_index % {
                    "index": self.quote_name('%s_%s_gist_idx' % (model._meta.db_table, field.column)),
                    "table": self.quote_name(model._meta.db_table),
                    "column": self.quote_name(field.column),
                }
            )
        return column_sql

    def create_model(self, model):
        super(DatabaseSchemaEditor, self).create_model(model)
        # create a gist index for ChemFields
        for sql in self.chem_sql:
            self.execute(sql)
        self.chem_sql = []

    def add_field(self, model, field):
        super(DatabaseSchemaEditor, self).add_field(model, field)
        #  create a gist index for ChemFields
        for sql in self.chem_sql:
            self.execute(sql)
        self.chem_sql = []
