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

    def _get_ontology_id(self,
                         term_name,
                         term_description = None,
                         term_identifier = None,
                         category_id = 0):
        """Get the id that corresponds to any term in an ontology.

        This looks through the ontology table for a the given term. If it
        is not found, a new id corresponding to this ontology is created.
        In either case, the id corresponding to that term is returned, so
        that you can reference it in another table.

        The category_id can be needed to disambiguate the term:
        it will be used if != 0.
        """

        # try to get the ontology term
        sql = r"SELECT ontology_term_id FROM ontology_term " \
              r"WHERE term_name = %s"
        fields = [term_name]
        if category_id != 0:            # 'None' is legitimate
            sql += ' AND category_id '
            if category_id is None:
                sql += 'IS NULL'
            else:
                sql += '= %s'
                fields.append(category_id)
        id_results = self.adaptor.execute_and_fetchall(sql, fields)
        # something is wrong
        if len(id_results) > 1:
            raise ValueError("Multiple ontology ids for %s: %s" % 
                             (term_name, id_results))
        # we already have the ontology term inserted
        elif len(id_results) == 1:
            return id_results[0][0]
        # we need to create it
        else:
            # If no category_id specified, set it to null, as 0 isn't possible
            if category_id == 0: category_id = None
            
            sql = r"INSERT INTO ontology_term (term_name, term_definition," \
                  r" term_identifier, category_id)" \
                  r" VALUES (%s, %s, %s, %s)"
            self.adaptor.execute(sql, (term_name, term_description,
                                       term_identifier, category_id))
            return self.adaptor.last_id('ontology_term')
   
    def _get_taxon_id(self, record):
        """Get the id corresponding to a taxon.

        If the species isn't in the taxon table, it is created.
        
        The code to find the species in the record is brittle.
        """
        # Binomial and full lineage
        try:
            binomial = record.annotations["organism"]
        except KeyError:
            binomial = None

        # XXX no variant
        variant = '-'

        if binomial and variant:
            sql = "SELECT taxon_id FROM taxon WHERE binomial = %s" \
                  " AND variant = %s"
            taxa = self.adaptor.execute_and_fetchall(sql, (binomial, variant))
            if taxa:
                return taxa[0][0]

        # Didn't found the binomial/variant... Let's try with the taxon id
        ncbi_taxon_id = None
        for f in record.features:
            if (f.type == 'source' and getattr(f, 'qualifiers', None)
                and f.qualifiers.has_key('db_xref')):
                for db_xref in f.qualifiers['db_xref']:
                    if db_xref[:6] == 'taxon:':
                        ncbi_taxon_id = int(db_xref[6:])
                        break
            if ncbi_taxon_id: break

        if ncbi_taxon_id:
            sql = "SELECT taxon_id FROM taxon WHERE ncbi_taxon_id = %s"
            taxa = self.adaptor.execute_and_fetchall(sql, (ncbi_taxon_id,))
            if taxa:
                return taxa[0][0]

        # OK, so we're gonna try to insert the taxon
        
        # Common name
        try:
            common_name = record.annotations["source"]
        except KeyError:
            common_name = None

        # Full lineage
        try:
            full_lineage = record.annotations["taxonomy"]
            ante, last = binomial.split()
            if full_lineage[-1] == ante:
                full_lineage.append(last)
            full_lineage.reverse()
            full_lineage = ':'.join(full_lineage)
        except KeyError:
            full_lineage = None

        # Check for the NON NULLs
        if binomial == None or variant == None or full_lineage == None:
            return
        
        # Insert into the taxon table
        sql = r"INSERT INTO taxon (binomial, variant, common_name," \
              r" ncbi_taxon_id, full_lineage)" \
              r" VALUES (%s, %s, %s, %s, %s)"
        self.adaptor.execute(sql, (binomial, variant, common_name,
                                   ncbi_taxon_id, full_lineage))
        taxon_id = self.adaptor.last_id('taxon')

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
            
        taxon_id = self._get_taxon_id(record)
        identifier = record.annotations.get('gi')
        description = getattr(record, 'description', None)
        
        sql = r"INSERT INTO bioentry (biodatabase_id, taxon_id, display_id, " \
              r"accession, identifier, description, entry_version) VALUES" \
              r" (%s, %s, %s, %s, %s, %s, %s)"
        self.adaptor.execute(sql, (self.dbid, taxon_id, record.name, 
                                   accession, identifier, description,
                                   version))
        # now retrieve the id for the bioentry
        bioentry_id = self.adaptor.last_id('bioentry')

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
        sql = r"INSERT INTO bioentry_qualifier_value" \
              r" (bioentry_id, ontology_term_id, qualifier_value)" \
              r" VALUES (%s, %s, %s)" 
        self.adaptor.execute(sql, (bioentry_id, date_id, date))

    def _load_biosequence(self, record, bioentry_id):
        """Load the biosequence table in the database.
        """
        accession, version = record.id.split(".")
        version = int(version)
        # determine the string representation of the alphabet
        if isinstance(record.seq.alphabet, Alphabet.DNAAlphabet):
            alphabet = "dna"
        elif isinstance(record.seq.alphabet, Alphabet.RNAAlphabet):
            alphabet = "rna"
        elif isinstance(record.seq.alphabet, Alphabet.ProteinAlphabet):
            alphabet = "protein"
        else:
            alphabet = "unknown"
        
        try:
            division = record.annotations["data_file_division"]
        except KeyError:
            division = "UNK"

        sql = r"INSERT INTO biosequence (bioentry_id, seq_version, " \
              r"seq_length, biosequence_str, alphabet, division) " \
              r"VALUES (%s, %s, %s, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, version,
                                   len(record.seq.data),
                                   record.seq.data,
                                   alphabet, division))

    def _load_comment(self, record, bioentry_id):
        # Assume annotations['comment'] is not a list
        comment = record.annotations.get('comment')
        if not comment:
            return
        comment = comment.replace('\n', ' ')
        
        sql = "INSERT INTO comment (bioentry_id, comment_text, comment_rank)" \
              " VALUES (%s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, comment, 1))
        
    def _load_reference(self, reference, rank, bioentry_id):
        s = ''
        for f in reference.authors, reference.title, reference.journal:
            if f: s += f
            else: s += "<undef>"
        doc_id = crc64(s)
        
        # The UK is either the medline id or the CRC64 'docid'
        if reference.medline_id:
            sql = r"SELECT reference_id FROM reference" \
                  r" WHERE reference_medline = %s"
            refs = self.adaptor.execute_and_fetch_col0(sql,
                                                       (reference.medline_id,))
        else:
            sql = r"SELECT reference_id FROM reference" \
                  r" WHERE reference_docid = %s"
            refs = self.adaptor.execute_and_fetch_col0(sql, (doc_id,))

        if not len(refs):
            authors = reference.authors or None
            title =  reference.title or None
            journal = reference.journal or None
            medline_id = reference.medline_id or None
            sql = r"INSERT INTO reference (reference_location," \
                  r" reference_title, reference_authors, reference_medline," \
                  r" reference_docid)" \
                  r" VALUES (%s, %s, %s, %s, %s)"
            self.adaptor.execute(sql, (journal, title,
                                   authors, medline_id, doc_id))
            reference_id = self.adaptor.last_id('reference')
        else:
            reference_id = refs[0]

        if len(reference.location):
            start = 1 + int(str(reference.location[0].start))
            end = int(str(reference.location[0].end))
        else:
            start = None
            end = None
        
        sql = "INSERT INTO bioentry_reference (bioentry_id, reference_id," \
              " reference_start, reference_end, reference_rank)" \
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
        category_id = self._get_ontology_id('SeqFeature Keys')
        seqfeature_key_id = self._get_ontology_id(feature_type,
                                                  category_id = category_id)
        
        # XXX source is always EMBL/GenBank/SwissProt here; it should depend on
        # the record
        source_cat_id = self._get_ontology_id('SeqFeature Sources')
        source_id = self._get_ontology_id('EMBL/GenBank/SwissProt',
                                          category_id = source_cat_id)
        
        sql = r"INSERT INTO seqfeature (bioentry_id, ontology_term_id, " \
              r"seqfeature_source_id, seqfeature_rank) VALUES (%s, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, seqfeature_key_id,
                                   source_id, feature_rank + 1))
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
            
        self.adaptor.execute(sql, (seqfeature_id, start, end, strand, rank+1))

    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):
        """Insert the (key, value) pair qualifiers relating to a feature.

        Qualifiers should be a dictionary of the form:
            {key : [value1, value2]}
        """
        tag_category_id = self._get_ontology_id('Annotation Tags')
        for qualifier_key in qualifiers.keys():
            qualifier_key_id = self._get_ontology_id(qualifier_key,
                                                     category_id = tag_category_id)

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

    XXX I think this might be the worst optimized SQL in the history
    of the world. There is probably a much better way to do it.
    [The "right" way is of course to have FKs--YB]
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
        

