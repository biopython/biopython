"""Load biopyton objects into a BioSQL database for persistant storage.

This code makes it possible to store biopython objects in a relational
database and then retrieve them back. You shouldn't use any of the
classes in this module directly. Rather, call the load() method on
a database object.
"""
# standard modules
from time import gmtime, strftime

class DatabaseLoader:
    """Load a database with biopython objects.
    """
    def __init__(self, adaptor, dbid):
        """Initialize with connection information for the database.

        XXX Figure out what I need to load a database and document it.
        """
        self.adaptor = adaptor
        self.dbid = dbid
    
    def load_seqrecord(self, record):
        """Load a Biopython SeqRecord into the database.
        """
        bioentry_id = self._load_bioentry_table(record)
        # self._load_bioentry_date(record, bioentry_id)
        # self._load_bioentry_taxa(record, bioentry_id)

    def _next_bioentry_id(self):
        """Find a unique bioentry id for a record.
        """
        all_ids = self.adaptor.all_bioentry_ids()
        if len(all_ids) >= 1:
            max_id = max(all_ids)
        else:
            max_id = 0L

        new_id = max_id + 1
        assert new_id not in all_ids, "Failed to creat a unique id"
        return new_id
    
    def _load_bioentry_table(self, record):
        """Fill the bioentry table with sequence information.
        """
        # get the stuff we need to fill up the table
        bioentry_id = self._next_bioentry_id()
        accession, version = record.id.split('.')
        try:
            division = record.annotations["data_file_divison"]
        except KeyError:
            division = "No"
        # now do it 
        sql = r"INSERT INTO bioentry VALUES" \
              r" (%s, %s, %s, %s, %s, %s)"
        self.adaptor.execute_one(sql, (bioentry_id, self.dbid, record.name, 
                                       accession, version, division))
        return bioentry_id

    def _load_bioentry_date(self, record, bioentry_id):
        """Add the effective date of the entry into the database.
        """
        # dates are GenBank style, like:
        # 14-SEP-2000
        try:
            date = record.annotations["date"]
        except KeyError:
            # just use today's date
            date = strftime("%d-%b-%Y", gmtime())
        sql = r"INSERT INTO bioentry_date VALUES" \
              r" (%s, %s)" 
        self.adaptor.execute_one(sql, (bioentry_id, date))

    def _load_bioentry_taxa(self, record, bioentry_id):
        """Add taxa information to the database.
        """
        return None # XXX don't do anything right now
        try:
            # XXX this isn't right, we need taxa ids and other junk
            taxa = record.annotations["taxa"]
            sql = r"INSERT INTO bioentry_taxa VALUES" \
                  r" (%s, %s)" 
            self.adapter.execute_one(sql, (bioentry_id, taxa))
        except KeyError:
            pass
        
