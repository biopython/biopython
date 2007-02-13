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
from Bio.crc import crc64

try:
    enumerate
except NameError:
    import sys
    _indices = xrange(sys.maxint)
    def enumerate(s):
        return zip(_indices, s)

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
        self._load_biosequence(record, bioentry_id)
        self._load_comment(record, bioentry_id)
        references = record.annotations.get('references', ())
        for reference, rank in zip(references, range(len(references))):
            self._load_reference(reference, rank, bioentry_id)
        for seq_feature_num in range(len(record.features)):
            seq_feature = record.features[seq_feature_num]
            self._load_seqfeature(seq_feature, seq_feature_num, bioentry_id)

    def _get_ontology_id(self, name, definition=None):
        oids = self.adaptor.execute_and_fetch_col0(
            "SELECT ontology_id FROM ontology WHERE name = %s",
            (name,))
        if oids:
            return oids[0]
        self.adaptor.execute(
            "INSERT INTO ontology(name, definition) VALUES (%s, %s)",
            (name, definition))
        return self.adaptor.last_id("ontology")

    
    def _get_term_id(self,
                     name,
                     ontology_id=None,
                     definition=None,
                     identifier=None):
        """Get the id that corresponds to a term.

        This looks through the term table for a the given term. If it
        is not found, a new id corresponding to this term is created.
        In either case, the id corresponding to that term is returned, so
        that you can reference it in another table.

        The ontology_id should be used to disambiguate the term.
        """

        # try to get the term id
        sql = r"SELECT term_id FROM term " \
              r"WHERE name = %s"
        fields = [name]
        if ontology_id:
            sql += ' AND ontology_id = %s'
            fields.append(ontology_id)
        id_results = self.adaptor.execute_and_fetchall(sql, fields)
        # something is wrong
        if len(id_results) > 1:
            raise ValueError("Multiple term ids for %s: %r" % 
                             (name, id_results))
        elif len(id_results) == 1:
            return id_results[0][0]
        else:
            sql = r"INSERT INTO term (name, definition," \
                  r" identifier, ontology_id)" \
                  r" VALUES (%s, %s, %s, %s)"
            self.adaptor.execute(sql, (name, definition,
                                       identifier, ontology_id))
            return self.adaptor.last_id("term")

    def _add_dbxref(self, dbname, accession, version):
       """Insert a dbxref and return its id"""
       
       self.adaptor.execute(
           "INSERT INTO dbxref(dbname, accession, version)" \
           " VALUES (%s, %s, %s)", (dbname, accession, version))
       return self.adaptor.last_id("dbxref")
           
    def _get_taxon_id(self, record):
        """Get the id corresponding to a taxon.

        If the species isn't in the taxon table, it is created.
        """

        ncbi_taxon_id = record.annotations.get("ncbi_taxid")
        if not ncbi_taxon_id:
            # Try the hard way...
            for f in record.features:
                if f.type == 'source':
                    quals = getattr(f, 'qualifiers', {})
                    if "db_xref" in quals:
                        for db_xref in f.qualifiers["db_xref"]:
                            if db_xref.startswith("taxon:"):
                                ncbi_taxon_id = int(db_xref[6:])
                                break
                            if ncbi_taxon_id: break
        if ncbi_taxon_id:
            taxa = self.adaptor.execute_and_fetch_col0(
                "SELECT taxon_id FROM taxon WHERE ncbi_taxon_id = %s",
                (ncbi_taxon_id,))
            if taxa:
                return taxa[0]

        # Tough luck. Let's try the binomial
        if record.annotations["organism"]:
            taxa = self.adaptor.execute_and_fetch_col0(
                "SELECT taxon_id FROM taxon_name" \
                " WHERE name_class = 'scientific name' AND name = %s",
                (record.annotations["organism"],))
            if taxa:
                return taxa[0]


        # Last chance...
        if record.annotations["source"]:
            taxa = self.adaptor.execute_and_fetch_col0(
                "SELECT DISTINCT taxon_id FROM taxon_name" \
                " WHERE name = %s",
                (record.annotations["source"],))
            if len(taxa) > 1:
                raise ValueError("Taxa: %d species have name %r" % (
                    len(taxa),
                    record.annotations["source"]))
            if taxa:
                return taxa[0]

        # OK, let's try inserting the species.
        # Chances are we don't have enough information ...
        # Furthermore, it won't be in the hierarchy.

        lineage = []
        for c in record.annotations.get("taxonomy", []):
            lineage.append([None, None, c])
        if lineage:
            lineage[-1][1] = "genus"
        lineage.append([None, "species", record.annotations["organism"]])
        # XXX do we have them?
        if "subspecies" in record.annotations:
            lineage.append([None, "subspecies",
                            record.annotations["subspecies"]])
        if "variant" in record.annotations:
            lineage.append([None, "varietas",
                            record.annotations["variant"]])
        lineage[-1][0] = ncbi_taxon_id
        
        left_value = self.adaptor.execute_one(
            "SELECT MAX(left_value) FROM taxon")[0]
        if not left_value:
            left_value = 0
        left_value += 1
        
        # XXX -- Brad: Fixing this for now in an ugly way because
        # I am getting overlaps for right_values. I need to dig into this
        # more to actually understand how it works. I'm not sure it is
        # actually working right anyhow.
        right_start_value = self.adaptor.execute_one(
            "SELECT MAX(right_value) FROM taxon")[0]
        if not right_start_value:
            right_start_value = 0
        right_value = right_start_value + 2 * len(lineage) - 1

        parent_taxon_id = None
        for taxon in lineage:
            self.adaptor.execute(
                "INSERT INTO taxon(parent_taxon_id, ncbi_taxon_id, node_rank,"\
                " left_value, right_value)" \
                " VALUES (%s, %s, %s, %s, %s)", (parent_taxon_id,
                                                 taxon[0],
                                                 taxon[1],
                                                 left_value,
                                                 right_value))
            taxon_id = self.adaptor.last_id("taxon")
            self.adaptor.execute(
                "INSERT INTO taxon_name(taxon_id, name, name_class)" \
                "VALUES (%s, %s, 'scientific name')", (taxon_id, taxon[2]))
            left_value += 1
            right_value -= 1
            parent_taxon_id = taxon_id
        if "source" in record.annotations:
            self.adaptor.execute(
                "INSERT INTO taxon_name(taxon_id, name, name_class)" \
                "VALUES (%s, %s, 'common name')", (
                taxon_id, record.annotations["source"]))

        return taxon_id

    def _load_bioentry_table(self, record):
        """Fill the bioentry table with sequence information.
        """
        # get the pertinent info and insert it
        
        if record.id.find('.') >= 0: # try to get a version from the id
            accession, version = record.id.split('.')
            version = int(version)
        else: # otherwise just use a version of 0
            accession = record.id
            version = 0
            
#        taxon_id = self._get_taxon_id(record)
        #taxon_id = "0" # inserted this because the taxon population code is out of date
                       # with the tables
        identifier = record.annotations.get('gi')
        description = getattr(record, 'description', None)
        division = record.annotations.get("data_file_division", "UNK")
        
	# removed taxon_id field, as it was causing difficulties with the 
	# schema  - not inserting a value allows it to default to NULL, 
	# avoiding the foreign key constraint.
        sql = """
        INSERT INTO bioentry (
         biodatabase_id,
         name,
         accession,
         identifier,
         division,
         description,
         version)
        VALUES (
         %s,
         %s,
         %s,
         %s,
         %s,
         %s,
         %s)"""
        self.adaptor.execute(sql, (self.dbid,
                                   record.name, 
                                   accession,
                                   identifier,
                                   division,
                                   description,
                                   version))
        # now retrieve the id for the bioentry
        bioentry_id = self.adaptor.last_id('bioentry')

        return bioentry_id

    def _load_bioentry_date(self, record, bioentry_id):
        """Add the effective date of the entry into the database.
        """
        # dates are GenBank style, like:
        # 14-SEP-2000
        date = record.annotations.get("date",
                                      strftime("%d-%b-%Y", gmtime()).upper())
        annotation_tags_id = self._get_ontology_id("Annotation Tags")
        date_id = self._get_term_id("date_changed", annotation_tags_id)
        sql = r"INSERT INTO bioentry_qualifier_value" \
              r" (bioentry_id, term_id, value, rank)" \
              r" VALUES (%s, %s, %s, 1)" 
        self.adaptor.execute(sql, (bioentry_id, date_id, date))

    def _load_biosequence(self, record, bioentry_id):
        """Load the biosequence table in the database.
        """
        # determine the string representation of the alphabet
        if isinstance(record.seq.alphabet, Alphabet.DNAAlphabet):
            alphabet = "dna"
        elif isinstance(record.seq.alphabet, Alphabet.RNAAlphabet):
            alphabet = "rna"
        elif isinstance(record.seq.alphabet, Alphabet.ProteinAlphabet):
            alphabet = "protein"
        else:
            alphabet = "unknown"

        sql = r"INSERT INTO biosequence (bioentry_id, version, " \
              r"length, seq, alphabet) " \
              r"VALUES (%s, 0, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id,
                                   len(record.seq.data),
                                   record.seq.data,
                                   alphabet))

    def _load_comment(self, record, bioentry_id):
        # Assume annotations['comment'] is not a list
        comment = record.annotations.get('comment')
        if not comment:
            return
        comment = comment.replace('\n', ' ')
        
        sql = "INSERT INTO comment (bioentry_id, comment_text, rank)" \
              " VALUES (%s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, comment, 1))
        
    def _load_reference(self, reference, rank, bioentry_id):

        refs = None
        if reference.medline_id:
            refs = self.adaptor.execute_and_fetch_col0(
                "SELECT reference_id" \
                "  FROM reference JOIN dbxref USING (dbxref_id)" \
                " WHERE dbname = 'MEDLINE' AND accession = %s",
                (reference.medline_id,))
        if not refs and reference.pubmed_id:
            refs = self.adaptor.execute_and_fetch_col0(
                "SELECT reference_id" \
                "  FROM reference JOIN dbxref USING (dbxref_id)" \
                " WHERE dbname = 'PUBMED' AND accession = %s",
                (reference.pubmed_id,))
        if not refs:
            s = []
            for f in reference.authors, reference.title, reference.journal:
                s.append(f or "<undef>")
            crc = crc64("".join(s))
            refs = self.adaptor.execute_and_fetch_col0(
                "SELECT reference_id FROM reference" \
                  r" WHERE crc = %s", (crc,))
        if not refs:
            if reference.medline_id:
                dbxref_id = self._add_dbxref("MEDLINE",
                                             reference.medline_id, 0)
            elif reference.pubmed_id:
                dbxref_id = self._add_dbxref("PUBMED",
                                             reference.pubmed_id, 0)
            else:
                dbxref_id = None
            authors = reference.authors or None
            title =  reference.title or None
            journal = reference.journal or None
            self.adaptor.execute(
                "INSERT INTO reference (dbxref_id, location," \
                " title, authors, crc)" \
                " VALUES (%s, %s, %s, %s, %s)",
                (dbxref_id, journal, title,
                 authors, crc))
            reference_id = self.adaptor.last_id("reference")
        else:
            reference_id = refs[0]

        if reference.location:
            start = 1 + int(str(reference.location[0].start))
            end = int(str(reference.location[0].end))
        else:
            start = None
            end = None
        
        sql = "INSERT INTO bioentry_reference (bioentry_id, reference_id," \
              " start_pos, end_pos, rank)" \
              " VALUES (%s, %s, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, reference_id,
                                   start, end, rank + 1))
        
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
        ontology_id = self._get_ontology_id('SeqFeature Keys')
        seqfeature_key_id = self._get_term_id(feature_type,
                                              ontology_id = ontology_id)
        
        # XXX source is always EMBL/GenBank/SwissProt here; it should depend on
        # the record (how?)
        source_cat_id = self._get_ontology_id('SeqFeature Sources')
        source_term_id = self._get_term_id('EMBL/GenBank/SwissProt',
                                      ontology_id = source_cat_id)
        
        sql = r"INSERT INTO seqfeature (bioentry_id, type_term_id, " \
              r"source_term_id, rank) VALUES (%s, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, seqfeature_key_id,
                                   source_term_id, feature_rank + 1))
        seqfeature_id = self.adaptor.last_id('seqfeature')

        return seqfeature_id

    def _load_seqfeature_locations(self, feature, seqfeature_id):
        """Load all of the locations for a SeqFeature into tables.

        This adds the locations related to the SeqFeature into the
        seqfeature_location table. Fuzzies are not handled right now.
        For a simple location, ie (1..2), we have a single table row
        with seq_start = 1, seq_end = 2, location_rank = 1.

        For split locations, ie (1..2, 3..4, 5..6) we would have three
        row tables with:
            start = 1, end = 2, rank = 1
            start = 3, end = 4, rank = 2
            start = 5, end = 6, rank = 3
        """
        # two cases, a simple location or a split location
        if not feature.sub_features:    # simple location
            self._insert_seqfeature_location(feature, 1, seqfeature_id)
        else: # split location
            for rank, cur_feature in enumerate(feature.sub_features):
                self._insert_seqfeature_location(cur_feature,
                                                 rank + 1,
                                                 seqfeature_id)

    def _insert_seqfeature_location(self, feature, rank, seqfeature_id):
        """Add a location of a SeqFeature to the seqfeature_location table.
        """
        sql = r"INSERT INTO location (seqfeature_id, " \
              r"start_pos, end_pos, strand, rank) " \
              r"VALUES (%s, %s, %s, %s, %s)"

        # convert biopython locations to the 1-based location system
        # used in bioSQL
        # XXX This could also handle fuzzies
        start = feature.location.nofuzzy_start + 1
        end = feature.location.nofuzzy_end
        # Biopython uses None when we don't know strand information but
        # BioSQL requires something (non null) and sets this as zero
        # So we'll use the strand or 0 if Biopython spits out None
        strand = feature.strand or 0

        self.adaptor.execute(sql, (seqfeature_id, start, end, strand, rank))

    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):
        """Insert the (key, value) pair qualifiers relating to a feature.

        Qualifiers should be a dictionary of the form:
            {key : [value1, value2]}
        """
        tag_ontology_id = self._get_ontology_id('Annotation Tags')
        for qualifier_key in qualifiers.keys():
            qualifier_key_id = self._get_term_id(qualifier_key,
                                                 ontology_id = tag_ontology_id)

            # now add all of the values to their table
            for qual_value_rank in range(len(qualifiers[qualifier_key])):
                qualifier_value = qualifiers[qualifier_key][qual_value_rank]
                sql = r"INSERT INTO seqfeature_qualifier_value VALUES" \
                      r" (%s, %s, %s, %s)"
                self.adaptor.execute(sql, (seqfeature_id,
                  qualifier_key_id, qual_value_rank + 1, qualifier_value))
       
class DatabaseRemover:
    """Complement the Loader functionality by fully removing a database.

    This probably isn't really useful for normal purposes, since you
    can just do a:
        DROP DATABASE db_name
    and then recreate the database. But, it's really useful for testing
    purposes.

    YB: now use the cascaded deletions
    """
    def __init__(self, adaptor, dbid):
        """Initialize with a database id and adaptor connection.
        """
        self.adaptor = adaptor
        self.dbid = dbid

    def remove(self):
        """Remove everything related to the given database id.
        """
        sql = r"DELETE FROM bioentry WHERE biodatabase_id = %s"
        self.adaptor.execute(sql, (self.dbid,))
        sql = r"DELETE FROM biodatabase WHERE biodatabase_id = %s"
        self.adaptor.execute(sql, (self.dbid,))
