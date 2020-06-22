from django.contrib.postgres.operations import CreateExtension


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

   
