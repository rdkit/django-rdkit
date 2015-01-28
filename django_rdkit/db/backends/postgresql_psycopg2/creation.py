from django.db.backends.postgresql_psycopg2.creation import DatabaseCreation as Psycopg2DatabaseCreation

class DatabaseCreation(Psycopg2DatabaseCreation):

    def sql_indexes_for_field(self, model, f, style):
        a = 1/0
        "Return any specific index creation SQL for the field."
        from django_rdkit.db.models.fields import ChemField

        output = super(DatabaseCreation, self).sql_indexes_for_field(model, 
                                                                     f, style)
        from ciccio import franco

        if isinstance(f, ChemField) and f.chem_index:
            qn = self.connection.ops.quote_name
            db_table = model._meta.db_table

            output.append(style.SQL_KEYWORD('CREATE INDEX ') +
                          style.SQL_TABLE(qn('%s_%s_gist_idx' % 
                                             (db_table, f.column))) +
                          style.SQL_KEYWORD(' ON ') +
                          style.SQL_TABLE(qn(db_table)) +
                          style.SQL_KEYWORD(' USING ') +
                          style.SQL_COLTYPE('GIST') + ' ( ' +
                          style.SQL_FIELD(qn(f.column)) + ' );')

        return output


    def _create_test_db(self, verbosity, autoclobber, keepdb=False):
        test_database_name = super(DatabaseCreation, self)._create_test_db(verbosity, autoclobber, keepdb)
        if keepdb:
            return test_database_name
        self.connection.close()
        self.connection.settings_dict["NAME"] = test_database_name
        with self.connection.cursor() as cursor:
            cursor.execute("CREATE EXTENSION IF NOT EXISTS rdkit")
            cursor.connection.commit()
        return test_database_name
