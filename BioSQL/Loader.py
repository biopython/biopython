"""Load biopyton objects into a BioSQL database for persistant storage.

This code makes it possible to store biopython objects in a relational
database and then retrieve them back. You shouldn't use any of the
classes in this module directly. Rather, call the load() method on
a database object.
"""
# standard modules
from time import gmtime, strftime

# biopython
from Bio import Alphabet

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
        self._load_bioentry_date(record, bioentry_id)
        # self._load_bioentry_taxa(record, bioentry_id)
        self._load_biosequence(record, bioentry_id)
        self._load_bioentry_description(record, bioentry_id)
        for seq_feature in record.features:
            self._load_seqfeature(seq_feature, bioentry_id)
    
    def _load_bioentry_table(self, record):
        """Fill the bioentry table with sequence information.
        """
        # get the pertinet info and insert it
        accession, version = record.id.split('.')
        try:
            division = record.annotations["data_file_divison"]
        except KeyError:
            division = "No"
        sql = r"INSERT INTO bioentry (biodatabase_id, display_id, " \
              r"accession, entry_version, division) VALUES" \
              r" (%s, %s, %s, %s, %s)"
        self.adaptor.execute_one(sql, (self.dbid, record.name, 
                                       accession, version, division))
        # now retrieve the id for the bioentry
        sql = r"SELECT max(bioentry_id) FROM bioentry"
        results = self.adaptor.execute_one(sql, ())
        bioentry_id = results[0]

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

    def _load_biosequence(self, record, bioentry_id):
        """Load the biosequence table in the database.
        """
        accession, version = record.id.split(".")
        # determine the string representation of the alphabet
        if isinstance(record.seq.alphabet, Alphabet.DNAAlphabet):
            mol_type = "DNA"
        elif isinstance(record.seq.alphabet, Alphabet.RNAAlphabet):
            mol_type = "RNA"
        elif isinstance(record.seq.alphabet, Alphabet.ProteinAlphabet):
            mol_type = "PROTEIN"
        else:
            mol_type = "UNKNOWN"
        
        sql = r"INSERT INTO biosequence (bioentry_id, seq_version, " \
              r"biosequence_str, molecule) VALUES (%s, %s, %s, %s)"
        self.adaptor.execute_one(sql, (bioentry_id, version, record.seq.data,
                                       mol_type))

    def _load_bioentry_description(self, record, bioentry_id):
        """Load the description table.
        """
        sql = r"INSERT INTO bioentry_description VALUES (%s, %s)"
        self.adaptor.execute_one(sql, (bioentry_id, record.description))

    def _load_seqfeature(self, feature, bioentry_id):
        """Load a biopython SeqFeature into the database.
        """
        seqfeature_id = self._load_seqfeature_basic(feature.type, 
                                                    bioentry_id)
        self._load_seqfeature_locations(feature, seqfeature_id)
        self._load_seqfeature_qualifiers(feature.qualifiers, seqfeature_id)

    def _load_seqfeature_basic(self, feature_type, bioentry_id):
        """Load the first tables of a seqfeature and returns the id.

        This loads the "key" of the seqfeature (ie. CDS, gene) and
        the basic seqfeature table itself.
        """
        sql = r"INSERT INTO seqfeature_key (key_name) VALUES (%s)"
        self.adaptor.execute_one(sql, (feature_type))
        sql = r"SELECT max(seqfeature_key_id) FROM seqfeature_key"
        results = self.adaptor.execute_one(sql, ())                
        seqfeature_key_id = results[0]               
        
        # XXX This doesn't do source or rank yet, since I'm not
        # sure I understand them.
        sql = r"INSERT INTO seqfeature (bioentry_id, seqfeature_key_id) " \
              r"VALUES (%s, %s)"
        self.adaptor.execute_one(sql, (bioentry_id, seqfeature_key_id))
        sql = r"SELECT max(seqfeature_id) FROM seqfeature"
        results = self.adaptor.execute_one(sql, ())
        seqfeature_id = results[0]

        return seqfeature_id

    def _load_seqfeature_locations(self, feature, seqfeature_id):
        """Load all of the locations for a SeqFeature into tables.

        This adds the locations related to the SeqFeature into the
        seqfeature_location table. Fuzzies are not handled right now.
        For a simple location, ie (1..2), we have a single table row
        with seq_start = 1, seq_end = 2, location_rank = 1.

        For split locations, ie (1..2, 3..4, 5..6) we would have three
        row tables with:
            seq_start = 1, seq_end = 2, location_rank = 1
            seq_start = 3, seq_end = 4, location_rank = 2
            seq_start = 5, seq_end = 6, location_rank = 3
        """
        # two cases, a simple location or a split location
        if len(feature.sub_features) == 0: # simple location
            self._insert_seqfeature_location(feature, 1, seqfeature_id)
        else: # split location
            for feature_rank in range(len(feature.sub_features)):
                cur_feature = feature.sub_features[feature_rank]
                self._insert_seqfeature_location(cur_feature, feature_rank,
                                                 seqfeature_id)

    def _insert_seqfeature_location(self, feature, rank, seqfeature_id):
        """Add a location of a SeqFeature to the seqfeature_location table.
        """
        sql = r"INSERT INTO seqfeature_location (seqfeature_id, " \
              r"seq_start, seq_end, seq_strand, location_rank) " \
               r"VALUES (%s, %s, %s, %s, %s)"

        # hack for NOT NULL in strand -- we have None be the same as 0
        # for strand information
        if feature.strand is None:
            strand = 0
        else:
            strand = feature.strand
            
        self.adaptor.execute_one(sql, (seqfeature_id, 
            feature.location.nofuzzy_start, feature.location.nofuzzy_end, 
            strand, rank))

    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):
        """Insert the (key, value) pair qualifiers relating to a feature.

        Qualifiers should be a dictionary of the form:
            {key : [value1, value2]}
        """
        for qualifier_key in qualifiers.keys():
            # add the key to the appropriate table
            sql = r"INSERT INTO seqfeature_qualifier (qualifier_name) " \
                  r"VALUES (%s)"
            self.adaptor.execute_one(sql, (qualifier_key))
            sql = r"SELECT max(seqfeature_qualifier_id) FROM " \
                  r"seqfeature_qualifier"
            results = self.adaptor.execute_one(sql, ())
            seqfeature_qualifier_id = results[0]

            # now add all of the values to their table
            for qual_value_rank in range(len(qualifiers[qualifier_key])):
                qualifier_value = qualifiers[qualifier_key][qual_value_rank]
                sql = r"INSERT INTO seqfeature_qualifier_value VALUES" \
                      r" (%s, %s, %s, %s)"
                self.adaptor.execute_one(sql, (seqfeature_id,
                  seqfeature_qualifier_id, qual_value_rank, qualifier_value))
       
class DatabaseRemover:
    """Compliment the Loader functionality by fully removing a database.

    This probably isn't really useful for normal purposes, since you
    can just do a:
        DROP DATABASE db_name
    and then recreate the database. But, it's really useful for testing
    purposes.

    XXX I think this might be the worst optimized SQL in the history
    of the world. There is probably a much better way to do it.
    """
    def __init__(self, adaptor, dbid):
        """Initialize with a database id and adaptor connection.
        """
        self.adaptor = adaptor
        self.dbid = dbid

    def remove(self):
        """Remove everything related to the given database id.
        """
        bioentry_ids = self.adaptor.list_bioentry_ids(self.dbid)

        # now remove all the entries
        for bioentry_id in bioentry_ids:
            self._remove_bioentry_basic(bioentry_id)
            self._remove_bioentry_features(bioentry_id)
            self._remove_bioentry_metadata(bioentry_id)

        # finally remove the database
        sql = r"DELETE FROM biodatabase WHERE biodatabase_id = %s"
        self.adaptor.execute_one(sql, (self.dbid))

    def _remove_bioentry_basic(self, bioentry_id):
        """Remove basic stuff relating to a bioentry.
        """
        sql = r"DELETE FROM bioentry WHERE bioentry_id = %s"
        self.adaptor.execute_one(sql, (bioentry_id))
        sql = r"DELETE FROM biosequence WHERE bioentry_id = %s"
        self.adaptor.execute_one(sql, (bioentry_id))

    def _remove_bioentry_metadata(self, bioentry_id):
        """Remove all metadata relating to a bioentry.
        """
        sql = r"DELETE FROM bioentry_date WHERE bioentry_id = %s"
        self.adaptor.execute_one(sql, (bioentry_id))
        sql = r"DELETE FROM bioentry_description WHERE bioentry_id = %s"
        self.adaptor.execute_one(sql, (bioentry_id))

    def _remove_bioentry_features(self, bioentry_id):
        """Remove all feature relating to a bioentry.
        """
        # first find all the relevant seqfeature ids
        sql = r"SELECT seqfeature_id FROM seqfeature WHERE " \
              r" bioentry_id = %s"
        seqfeature_ids = self.adaptor.list_any_ids(sql, (bioentry_id))
        for seqfeature_id in seqfeature_ids:
            # location
            sql = r"DELETE FROM seqfeature_location WHERE seqfeature_id = %s"
            self.adaptor.cursor.execute(sql, (seqfeature_id))
            # qualifiers
            sql = r"SELECT seqfeature_qualifier_id FROM "\
                  r"seqfeature_qualifier_value WHERE seqfeature_id = %s"
            qualifier_ids = self.adaptor.list_any_ids(sql, (seqfeature_id))
            for qualifier_id in qualifier_ids:
                sql = r"DELETE FROM seqfeature_qualifier WHERE "\
                      r"seqfeature_qualifier_id = %s"
                self.adaptor.cursor.execute(sql, (qualifier_id))
            sql = r"DELETE FROM seqfeature_qualifier_value WHERE " \
                  r"seqfeature_id = %s"
            self.adaptor.cursor.execute(sql, (seqfeature_id))

        # delete the main seqfeature table and its key
        sql = r"SELECT seqfeature_key_id FROM seqfeature WHERE " \
              r"bioentry_id = %s"
        seqfeature_key_ids = self.adaptor.list_any_ids(sql, (bioentry_id))
        for seqfeature_key_id in seqfeature_key_ids:
            sql = r"DELETE FROM seqfeature_key WHERE seqfeature_key_id = %s"
            self.adaptor.execute_one(sql, (seqfeature_key_id))
        sql = r"DELETE FROM seqfeature WHERE bioentry_id = %s"
        self.adaptor.cursor.execute(sql, (bioentry_id))
        

