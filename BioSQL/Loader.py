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
        for seq_feature_num in range(len(record.features)):
            seq_feature = record.features[seq_feature_num]
            self._load_seqfeature(seq_feature, seq_feature_num, bioentry_id)

    def _get_ontology_id(self, term_name, term_description = ""):
        """Get the id that corresponds to any term in an ontology.

        This looks through the ontology table for a the given term. If it
        is not found, a new id corresponding to this ontology is created.
        In either case, the id corresponding to that term is returned, so
        that you can reference it in another table.
        """
        # try to get the ontology term
        sql = r"SELECT ontology_term_id FROM ontology_term " \
              r"WHERE term_name = %s"
        id_results = self.adaptor.execute_and_fetchall(sql, (term_name,))
        # something is wrong
        if len(id_results) > 1:
            raise ValueError("Multiple ontology ids for %s: %s" % 
                             term_name, id_results)
        # we already have the ontology term inserted
        elif len(id_results) == 1:
            return id_results[0][0]
        # we need to create it
        else:
            sql = r"INSERT INTO ontology_term (term_name, term_definition)" \
                  r"VALUES (%s, %s)"
            self.adaptor.execute_one(sql, (term_name, term_description))
            # recursively call this to give back the id
            return self._get_ontology_id(term_name, term_description)
   
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
        date_id = self._get_ontology_id("date", "Sequence date")
        sql = r"INSERT INTO bioentry_qualifier_value VALUES" \
              r" (%s, %s, %s)" 
        self.adaptor.execute_one(sql, (bioentry_id, date_id, date))

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
        descr_id = self._get_ontology_id("description", "Sequence description")
        sql = r"INSERT INTO bioentry_qualifier_value VALUES (%s, %s, %s)"
        self.adaptor.execute_one(sql, (bioentry_id, descr_id, 
                                       record.description))

    def _load_seqfeature(self, feature, feature_rank, bioentry_id):
        """Load a biopython SeqFeature into the database.
        """
        seqfeature_id = self._load_seqfeature_basic(feature.type, feature_rank,
                                                    bioentry_id)
        self._load_seqfeature_locations(feature, seqfeature_id)
        self._load_seqfeature_qualifiers(feature.qualifiers, seqfeature_id)

    def _load_seqfeature_basic(self, feature_type, feature_rank, bioentry_id):
        """Load the first tables of a seqfeature and returns the id.

        This loads the "key" of the seqfeature (ie. CDS, gene) and
        the basic seqfeature table itself.
        """
        seqfeature_key_id = self._get_ontology_id(feature_type)
        
        # XXX This doesn't do source yet, since I'm not sure I understand it.
        sql = r"INSERT INTO seqfeature (bioentry_id, seqfeature_key_id, " \
              r"seqfeature_rank) VALUES (%s, %s, %s)"
        self.adaptor.execute_one(sql, (bioentry_id, seqfeature_key_id,
                                       feature_rank))
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

        # convert biopython locations to the 1-based location system
        # used in bioSQL
        # XXX This could also handle fuzzies
        start = feature.location.nofuzzy_start + 1
        end = feature.location.nofuzzy_end 
            
        self.adaptor.execute_one(sql, (seqfeature_id, start, end, strand, rank))

    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):
        """Insert the (key, value) pair qualifiers relating to a feature.

        Qualifiers should be a dictionary of the form:
            {key : [value1, value2]}
        """
        for qualifier_key in qualifiers.keys():
            qualifier_key_id = self._get_ontology_id(qualifier_key)

            # now add all of the values to their table
            for qual_value_rank in range(len(qualifiers[qualifier_key])):
                qualifier_value = qualifiers[qualifier_key][qual_value_rank]
                sql = r"INSERT INTO seqfeature_qualifier_value VALUES" \
                      r" (%s, %s, %s, %s)"
                self.adaptor.execute_one(sql, (seqfeature_id,
                  qualifier_key_id, qual_value_rank, qualifier_value))
       
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
        self.adaptor.execute(sql, (self.dbid))

    def _remove_bioentry_basic(self, bioentry_id):
        """Remove basic stuff relating to a bioentry.
        """
        sql = r"DELETE FROM bioentry WHERE bioentry_id = %s"
        self.adaptor.execute(sql, (bioentry_id))
        sql = r"DELETE FROM biosequence WHERE bioentry_id = %s"
        self.adaptor.execute(sql, (bioentry_id))

    def _remove_bioentry_metadata(self, bioentry_id):
        """Remove all metadata relating to a bioentry.
        """
        sql = r"DELETE FROM bioentry_qualifier_value WHERE bioentry_id = %s"
        self.adaptor.execute(sql, (bioentry_id))

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
            sql = r"DELETE FROM seqfeature_qualifier_value WHERE " \
                  r"seqfeature_id = %s"
            self.adaptor.cursor.execute(sql, (seqfeature_id))

        sql = r"DELETE FROM seqfeature WHERE bioentry_id = %s"
        self.adaptor.cursor.execute(sql, (bioentry_id))
        

