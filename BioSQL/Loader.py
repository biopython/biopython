# Copyright 2002 by Andrew Dalke.  All rights reserved.
# Revisions 2007-2009 copyright by Peter Cock.  All rights reserved.
# Revisions 2008 copyright by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that BioSQL (including the database schema and scripts) is
# available and licensed separately.  Please consult www.biosql.org

"""Load biopython objects into a BioSQL database for persistent storage.

This code makes it possible to store biopython objects in a relational
database and then retrieve them back. You shouldn't use any of the
classes in this module directly. Rather, call the load() method on
a database object.
"""
# standard modules
from time import gmtime, strftime

# biopython
from Bio import Alphabet
from Bio.SeqUtils.CheckSum import crc64
from Bio import Entrez
from Bio.Seq import UnknownSeq

from Bio._py3k import _is_int_or_long

class DatabaseLoader:
    """Object used to load SeqRecord objects into a BioSQL database."""
    def __init__(self, adaptor, dbid, fetch_NCBI_taxonomy=False):
        """Initialize with connection information for the database.

        Creating a DatabaseLoader object is normally handled via the
        BioSeqDatabase DBServer object, for example:

        from BioSQL import BioSeqDatabase
        server = BioSeqDatabase.open_database(driver="MySQLdb", user="gbrowse",
                         passwd = "biosql", host = "localhost", db="test_biosql")
        try:
            db = server["test"]
        except KeyError:
            db = server.new_database("test", description="For testing GBrowse")
        """
        self.adaptor = adaptor
        self.dbid = dbid
        self.fetch_NCBI_taxonomy = fetch_NCBI_taxonomy
    
    def load_seqrecord(self, record):
        """Load a Biopython SeqRecord into the database.
        """
        bioentry_id = self._load_bioentry_table(record)
        self._load_bioentry_date(record, bioentry_id)
        self._load_biosequence(record, bioentry_id)
        self._load_comment(record, bioentry_id)
        self._load_dbxrefs(record, bioentry_id)
        references = record.annotations.get('references', ())
        for reference, rank in zip(references, range(len(references))):
            self._load_reference(reference, rank, bioentry_id)
        self._load_annotations(record, bioentry_id)
        for seq_feature_num in range(len(record.features)):
            seq_feature = record.features[seq_feature_num]
            self._load_seqfeature(seq_feature, seq_feature_num, bioentry_id)

    def _get_ontology_id(self, name, definition=None):
        """Returns the identifier for the named ontology (PRIVATE).

        This looks through the onotology table for a the given entry name.
        If it is not found, a row is added for this ontology (using the
        definition if supplied).  In either case, the id corresponding to
        the provided name is returned, so that you can reference it in
        another table.
        """
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
        """Get the id that corresponds to a term (PRIVATE).

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
       """Insert a dbxref and return its id."""
       
       self.adaptor.execute(
           "INSERT INTO dbxref(dbname, accession, version)" \
           " VALUES (%s, %s, %s)", (dbname, accession, version))
       return self.adaptor.last_id("dbxref")
           
    def _get_taxon_id(self, record):
        """Get the taxon id for this record (PRIVATE).

        record - a SeqRecord object

        This searches the taxon/taxon_name tables using the
        NCBI taxon ID, scientific name and common name to find
        the matching taxon table entry's id.
        
        If the species isn't in the taxon table, and we have at
        least the NCBI taxon ID, scientific name or common name,
        at least a minimal stub entry is created in the table.

        Returns the taxon id (database key for the taxon table,
        not an NCBI taxon ID), or None if the taxonomy information
        is missing.

        See also the BioSQL script load_ncbi_taxonomy.pl which
        will populate and update the taxon/taxon_name tables
        with the latest information from the NCBI.
        """
        
        # To find the NCBI taxid, first check for a top level annotation
        ncbi_taxon_id = None
        if "ncbi_taxid" in record.annotations:
            #Could be a list of IDs.
            if isinstance(record.annotations["ncbi_taxid"],list):
                if len(record.annotations["ncbi_taxid"])==1:
                    ncbi_taxon_id = record.annotations["ncbi_taxid"][0]
            else:
                ncbi_taxon_id = record.annotations["ncbi_taxid"]
        if not ncbi_taxon_id:
            # Secondly, look for a source feature
            for f in record.features:
                if f.type == 'source':
                    quals = getattr(f, 'qualifiers', {})
                    if "db_xref" in quals:
                        for db_xref in f.qualifiers["db_xref"]:
                            if db_xref.startswith("taxon:"):
                                ncbi_taxon_id = int(db_xref[6:])
                                break
                if ncbi_taxon_id: break

        try:
            scientific_name = record.annotations["organism"][:255]
        except KeyError:
            scientific_name = None
        try:
            common_name = record.annotations["source"][:255]
        except KeyError:
            common_name = None
        # Note: The maximum length for taxon names in the schema is 255.
        # Cropping it now should help in getting a match when searching,
        # and avoids an error if we try and add these to the database.


        if ncbi_taxon_id:
            #Good, we have the NCBI taxon to go on - this is unambiguous :)
            #Note that the scientific name and common name will only be
            #used if we have to record a stub entry.
            return self._get_taxon_id_from_ncbi_taxon_id(ncbi_taxon_id,
                                                         scientific_name,
                                                         common_name)
        
        if not common_name and not scientific_name:
            # Nothing to go on... and there is no point adding
            # a new entry to the database.  We'll just leave this
            # sequence's taxon as a NULL in the database.
            return None

        # Next, we'll try to find a match based on the species name
        # (stored in GenBank files as the organism and/or the source).
        if scientific_name:
            taxa = self.adaptor.execute_and_fetch_col0(
                "SELECT taxon_id FROM taxon_name" \
                " WHERE name_class = 'scientific name' AND name = %s",
                (scientific_name,))
            if taxa:
                #Good, mapped the scientific name to a taxon table entry
                return taxa[0]

        # Last chance...
        if common_name:
            taxa = self.adaptor.execute_and_fetch_col0(
                "SELECT DISTINCT taxon_id FROM taxon_name" \
                " WHERE name = %s",
                (common_name,))
            #Its natural that several distinct taxa will have the same common
            #name - in which case we can't resolve the taxon uniquely.
            if len(taxa) > 1:
                raise ValueError("Taxa: %d species have name %r" % (
                    len(taxa),
                    common_name))
            if taxa:
                #Good, mapped the common name to a taxon table entry
                return taxa[0]

        # At this point, as far as we can tell, this species isn't
        # in the taxon table already.  So we'll have to add it.
        # We don't have an NCBI taxonomy ID, so if we do record just
        # a stub entry, there is no simple way to fix this later.
        #
        # TODO - Should we try searching the NCBI taxonomy using the
        # species name?
        #
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
                "VALUES (%s, %s, 'scientific name')", (taxon_id, taxon[2][:255]))
            #Note the name field is limited to 255, some SwissProt files
            #have a multi-species name which can be longer.  So truncate this.
            left_value += 1
            right_value -= 1
            parent_taxon_id = taxon_id
        if common_name:
            self.adaptor.execute(
                "INSERT INTO taxon_name(taxon_id, name, name_class)" \
                "VALUES (%s, %s, 'common name')", (
                taxon_id, common_name))

        return taxon_id

    def _fix_name_class(self, entrez_name):
        """Map Entrez name terms to those used in taxdump (PRIVATE).

        We need to make this conversion to match the taxon_name.name_class
        values used by the BioSQL load_ncbi_taxonomy.pl script.
        
        e.g.
        "ScientificName" -> "scientific name",
        "EquivalentName" -> "equivalent name",
        "Synonym" -> "synonym",
        """
        #Add any special cases here:
        #
        #known = {}
        #try:
        #    return known[entrez_name]
        #except KeyError:
        #    pass

        #Try automatically by adding spaces before each capital
        def add_space(letter):
            if letter.isupper():
                return " "+letter.lower()
            else:
                return letter
        answer = "".join([add_space(letter) for letter in entrez_name]).strip()
        assert answer == answer.lower()
        return answer

    def _get_taxon_id_from_ncbi_taxon_id(self, ncbi_taxon_id,
                                         scientific_name = None,
                                         common_name = None):
        """Get the taxon id for this record from the NCBI taxon ID (PRIVATE).

        ncbi_taxon_id - string containing an NCBI taxon id
        scientific_name - string, used if a stub entry is recorded
        common_name - string, used if a stub entry is recorded
        
        This searches the taxon table using ONLY the NCBI taxon ID
        to find the matching taxon table entry's ID (database key).
        
        If the species isn't in the taxon table, and the fetch_NCBI_taxonomy
        flag is true, Biopython will attempt to go online using Bio.Entrez
        to fetch the official NCBI lineage, recursing up the tree until an
        existing entry is found in the database or the full lineage has been
        fetched.

        Otherwise the NCBI taxon ID, scientific name and common name are
        recorded as a minimal stub entry in the taxon and taxon_name tables.
        Any partial information about the lineage from the SeqRecord is NOT
        recorded.  This should mean that (re)running the BioSQL script
        load_ncbi_taxonomy.pl can fill in the taxonomy lineage.

        Returns the taxon id (database key for the taxon table, not
        an NCBI taxon ID).
        """
        assert ncbi_taxon_id

        taxon_id = self.adaptor.execute_and_fetch_col0(
            "SELECT taxon_id FROM taxon WHERE ncbi_taxon_id = %s",
            (ncbi_taxon_id,))
        if taxon_id:
            #Good, we have mapped the NCBI taxid to a taxon table entry
            return taxon_id[0]

        # At this point, as far as we can tell, this species isn't
        # in the taxon table already.  So we'll have to add it.

        parent_taxon_id = None
        rank = "species"
        genetic_code = None
        mito_genetic_code = None
        species_names = []
        if scientific_name:
            species_names.append(("scientific name", scientific_name))
        if common_name:
            species_names.append(("common name", common_name))
        
        if self.fetch_NCBI_taxonomy:
            #Go online to get the parent taxon ID!
            handle = Entrez.efetch(db="taxonomy",id=ncbi_taxon_id,retmode="XML")
            taxonomic_record = Entrez.read(handle)
            if len(taxonomic_record) == 1:
                assert taxonomic_record[0]["TaxId"] == str(ncbi_taxon_id), \
                       "%s versus %s" % (taxonomic_record[0]["TaxId"],
                                         ncbi_taxon_id)
                parent_taxon_id = self._get_taxon_id_from_ncbi_lineage( \
                                            taxonomic_record[0]["LineageEx"])
                rank = taxonomic_record[0]["Rank"]
                genetic_code = taxonomic_record[0]["GeneticCode"]["GCId"]
                mito_genetic_code = taxonomic_record[0]["MitoGeneticCode"]["MGCId"]
                species_names = [("scientific name",
                                  taxonomic_record[0]["ScientificName"])]
                try:
                    for name_class, names in taxonomic_record[0]["OtherNames"].iteritems():
                        name_class = self._fix_name_class(name_class)
                        if not isinstance(names, list):
                            #The Entrez parser seems to return single entry
                            #lists as just a string which is annoying.
                            names = [names]
                        for name in names:
                            #Want to ignore complex things like ClassCDE entries
                            if isinstance(name, basestring):
                                species_names.append((name_class, name))
                except KeyError:
                    #OtherNames isn't always present,
                    #e.g. NCBI taxon 41205, Bromheadia finlaysoniana
                    pass
        else:
            pass
            # If we are not allowed to go online, we will record the bare minimum;
            # as long as the NCBI taxon id is present, then (re)running
            # load_ncbi_taxonomy.pl should fill in the taxonomomy lineage
            # (and update the species names).
            #
            # I am NOT going to try and record the lineage, even if it
            # is in the record annotation as a list of names, as we won't
            # know the NCBI taxon IDs for these parent nodes.

        self.adaptor.execute(
            "INSERT INTO taxon(parent_taxon_id, ncbi_taxon_id, node_rank,"\
            " genetic_code, mito_genetic_code, left_value, right_value)" \
            " VALUES (%s, %s, %s, %s, %s, %s, %s)", (parent_taxon_id,
                                                     ncbi_taxon_id,
                                                     rank,
                                                     genetic_code,
                                                     mito_genetic_code,
                                                     None,
                                                     None))
        taxon_id = self.adaptor.last_id("taxon")

        #Record the scientific name, common name, etc
        for name_class, name in species_names:
            self.adaptor.execute(
                "INSERT INTO taxon_name(taxon_id, name, name_class)" \
                " VALUES (%s, %s, %s)", (taxon_id, 
                                         name[:255], 
                                         name_class))
        return taxon_id

    def _get_taxon_id_from_ncbi_lineage(self, taxonomic_lineage):
        """This is recursive! (PRIVATE).

        taxonomic_lineage - list of taxonomy dictionaries from Bio.Entrez

        First dictionary in list is the taxonomy root, highest would be the species.
        Each dictionary includes:
        - TaxID (string, NCBI taxon id)
        - Rank (string, e.g. "species", "genus", ..., "phylum", ...)
        - ScientificName (string)
        (and that is all at the time of writing)

        This method will record all the lineage given, returning the the taxon id
        (database key, not NCBI taxon id) of the final entry (the species).
        """
        ncbi_taxon_id = taxonomic_lineage[-1]["TaxId"]

        #Is this in the database already?  Check the taxon table...
        taxon_id = self.adaptor.execute_and_fetch_col0(
            "SELECT taxon_id FROM taxon" \
            " WHERE ncbi_taxon_id=%s" % ncbi_taxon_id)
        if taxon_id:
            # we could verify that the Scientific Name etc in the database
            # is the same and update it or print a warning if not...
            if isinstance(taxon_id, list):
                assert len(taxon_id)==1
                return taxon_id[0]
            else:
                return taxon_id

        #We have to record this.
        if len(taxonomic_lineage) > 1:
            #Use recursion to find out the taxon id (database key) of the parent.
            parent_taxon_id = self._get_taxon_id_from_ncbi_lineage(taxonomic_lineage[:-1])
            assert _is_int_or_long(parent_taxon_id), repr(parent_taxon_id)
        else:
            parent_taxon_id = None

        # INSERT new taxon
        rank = taxonomic_lineage[-1].get("Rank", None)
        self.adaptor.execute(
                "INSERT INTO taxon(ncbi_taxon_id, parent_taxon_id, node_rank)"\
                " VALUES (%s, %s, %s)", (ncbi_taxon_id, parent_taxon_id, rank))
        taxon_id = self.adaptor.last_id("taxon")
        assert isinstance(taxon_id, int) or isinstance(taxon_id, long), repr(taxon_id)
        # ... and its name in taxon_name
        scientific_name = taxonomic_lineage[-1].get("ScientificName", None)
        if scientific_name:
            self.adaptor.execute(
                    "INSERT INTO taxon_name(taxon_id, name, name_class)" \
                    " VALUES (%s, %s, 'scientific name')", (taxon_id, 
                                                            scientific_name[:255]))
        return taxon_id


    def _load_bioentry_table(self, record):
        """Fill the bioentry table with sequence information (PRIVATE).

        record - SeqRecord object to add to the database.
        """
        # get the pertinent info and insert it
        
        if record.id.count(".") == 1: # try to get a version from the id
            #This assumes the string is something like "XXXXXXXX.123"
            accession, version = record.id.split('.')
            try:
                version = int(version)
            except ValueError:
                accession = record.id
                version = 0
        else: # otherwise just use a version of 0
            accession = record.id
            version = 0

        if "accessions" in record.annotations \
        and isinstance(record.annotations["accessions"],list) \
        and record.annotations["accessions"]:
            #Take the first accession (one if there is more than one)
            accession = record.annotations["accessions"][0]

        #Find the taxon id (this is not just the NCBI Taxon ID)
        #NOTE - If the species isn't defined in the taxon table,
        #a new minimal entry is created.
        taxon_id = self._get_taxon_id(record)

        if "gi" in record.annotations:
            identifier = record.annotations["gi"]
        else:
            identifier = record.id

        #Allow description and division to default to NULL as in BioPerl.
        description = getattr(record, 'description', None)
        division = record.annotations.get("data_file_division", None)
        
        sql = """
        INSERT INTO bioentry (
         biodatabase_id,
         taxon_id,
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
         %s,
         %s)"""
        #print self.dbid, taxon_id, record.name, accession, identifier, \
        #        division, description, version
        self.adaptor.execute(sql, (self.dbid,
                                   taxon_id,
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

        record - a SeqRecord object with an annotated date
        bioentry_id - corresponding database identifier
        """
        # dates are GenBank style, like:
        # 14-SEP-2000
        date = record.annotations.get("date",
                                      strftime("%d-%b-%Y", gmtime()).upper())
        if isinstance(date, list) : date = date[0]
        annotation_tags_id = self._get_ontology_id("Annotation Tags")
        date_id = self._get_term_id("date_changed", annotation_tags_id)
        sql = r"INSERT INTO bioentry_qualifier_value" \
              r" (bioentry_id, term_id, value, rank)" \
              r" VALUES (%s, %s, %s, 1)" 
        self.adaptor.execute(sql, (bioentry_id, date_id, date))

    def _load_biosequence(self, record, bioentry_id):
        """Record a SeqRecord's sequence and alphabet in the database (PRIVATE).

        record - a SeqRecord object with a seq property
        bioentry_id - corresponding database identifier
        """
        if record.seq is None:
            #The biosequence table entry is optional, so if we haven't
            #got a sequence, we don't need to write to the table.
            return
        
        # determine the string representation of the alphabet
        if isinstance(record.seq.alphabet, Alphabet.DNAAlphabet):
            alphabet = "dna"
        elif isinstance(record.seq.alphabet, Alphabet.RNAAlphabet):
            alphabet = "rna"
        elif isinstance(record.seq.alphabet, Alphabet.ProteinAlphabet):
            alphabet = "protein"
        else:
            alphabet = "unknown"

        if isinstance(record.seq, UnknownSeq):
            seq_str = None
        else:
            seq_str = str(record.seq)

        sql = r"INSERT INTO biosequence (bioentry_id, version, " \
              r"length, seq, alphabet) " \
              r"VALUES (%s, 0, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id,
                                   len(record.seq),
                                   seq_str,
                                   alphabet))

    def _load_comment(self, record, bioentry_id):
        """Record a SeqRecord's annotated comment in the database (PRIVATE).

        record - a SeqRecord object with an annotated comment
        bioentry_id - corresponding database identifier
        """
        comments = record.annotations.get('comment')
        if not comments:
            return
        if not isinstance(comments, list):
            #It should be a string then...
            comments = [comments]

        for index, comment in enumerate(comments):
            comment = comment.replace('\n', ' ')
            #TODO - Store each line as a separate entry?  This would preserve
            #the newlines, but we should check BioPerl etc to be consistent.
            sql = "INSERT INTO comment (bioentry_id, comment_text, rank)" \
                  " VALUES (%s, %s, %s)"
            self.adaptor.execute(sql, (bioentry_id, comment, index+1))
        
    def _load_annotations(self, record, bioentry_id):
        """Record a SeqRecord's misc annotations in the database (PRIVATE).

        The annotation strings are recorded in the bioentry_qualifier_value
        table, except for special cases like the reference, comment and
        taxonomy which are handled with their own tables.

        record - a SeqRecord object with an annotations dictionary
        bioentry_id - corresponding database identifier
        """
        mono_sql = "INSERT INTO bioentry_qualifier_value" \
                   "(bioentry_id, term_id, value)" \
                   " VALUES (%s, %s, %s)"
        many_sql = "INSERT INTO bioentry_qualifier_value" \
                   "(bioentry_id, term_id, value, rank)" \
                   " VALUES (%s, %s, %s, %s)"
        tag_ontology_id = self._get_ontology_id('Annotation Tags')
        for key, value in record.annotations.iteritems():
            if key in ["references", "comment", "ncbi_taxid", "date"]:
                #Handled separately
                continue
            term_id = self._get_term_id(key, ontology_id=tag_ontology_id)
            if isinstance(value, list) or isinstance(value, tuple):
                rank = 0
                for entry in value:
                    if isinstance(entry, str) or isinstance(entry, int):
                        #Easy case
                        rank += 1
                        self.adaptor.execute(many_sql, \
                                     (bioentry_id, term_id, str(entry), rank))
                    else:
                        pass
                        #print "Ignoring annotation '%s' sub-entry of type '%s'" \
                        #      % (key, str(type(entry)))
            elif isinstance(value, str) or isinstance(value, int):
                #Have a simple single entry, leave rank as the DB default
                self.adaptor.execute(mono_sql, \
                                     (bioentry_id, term_id, str(value)))
            else:
                pass
                #print "Ignoring annotation '%s' entry of type '%s'" \
                #      % (key, type(value))


    def _load_reference(self, reference, rank, bioentry_id):
        """Record a SeqRecord's annotated references in the database (PRIVATE).

        record - a SeqRecord object with annotated references
        bioentry_id - corresponding database identifier
        """

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
            #The location/journal field cannot be Null, so default
            #to an empty string rather than None:
            journal = reference.journal or ""
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
        """Load a biopython SeqFeature into the database (PRIVATE).
        """
        seqfeature_id = self._load_seqfeature_basic(feature.type, feature_rank,
                                                    bioentry_id)
        self._load_seqfeature_locations(feature, seqfeature_id)
        self._load_seqfeature_qualifiers(feature.qualifiers, seqfeature_id)

    def _load_seqfeature_basic(self, feature_type, feature_rank, bioentry_id):
        """Load the first tables of a seqfeature and returns the id (PRIVATE).

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
        """Load all of the locations for a SeqFeature into tables (PRIVATE).

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
        # TODO - Record an ontology for the locations (using location.term_id)
        # which for now as in BioPerl we leave defaulting to NULL.
        if feature.location_operator and feature.location_operator != "join":
            # e.g. order locations... we don't record "order" so it
            # will become a "join" on reloading. What does BioPerl do?
            import warnings
            warnings.warn("%s location operators are not fully supported" \
                          % feature.location_operator)
        
        # two cases, a simple location or a split location
        if not feature.sub_features:    # simple location
            self._insert_seqfeature_location(feature, 1, seqfeature_id)
        else: # split location
            for rank, cur_feature in enumerate(feature.sub_features):
                self._insert_seqfeature_location(cur_feature,
                                                 rank + 1,
                                                 seqfeature_id)

    def _insert_seqfeature_location(self, feature, rank, seqfeature_id):
        """Add a location of a SeqFeature to the seqfeature_location table (PRIVATE).

        TODO - Add location_operators to location_qualifier_value.
        """
        # convert biopython locations to the 1-based location system
        # used in bioSQL
        # XXX This could also handle fuzzies
        start = int(feature.location.start) + 1
        end = int(feature.location.end)

        # Biopython uses None when we don't know strand information but
        # BioSQL requires something (non null) and sets this as zero
        # So we'll use the strand or 0 if Biopython spits out None
        strand = feature.strand or 0

        # TODO - Record an ontology term for the location (location.term_id)
        # which for now like BioPerl we'll leave as NULL.
        # This might allow us to record "between" positions properly, but I
        # doesn't really see how it could work for before/after fuzzy positions
        loc_term_id = None

        if feature.ref:
            # sub_feature remote locations when they are in the same db as the current
            # record do not have a value for ref_db, which the SeqFeature object
            # stores as None. BioSQL schema requires a varchar and is not NULL 
            dbxref_id = self._get_dbxref_id(feature.ref_db or "", feature.ref)
        else:
            dbxref_id = None

        sql = r"INSERT INTO location (seqfeature_id, dbxref_id, term_id," \
              r"start_pos, end_pos, strand, rank) " \
              r"VALUES (%s, %s, %s, %s, %s, %s, %s)"
        self.adaptor.execute(sql, (seqfeature_id, dbxref_id, loc_term_id,
                                   start, end, strand, rank))

        """
        # See Bug 2677
        # TODO - Record the location_operator (e.g. "join" or "order")
        # using the location_qualifier_value table (which we and BioPerl
        # have historically left empty).
        # Note this will need an ontology term for the location qualifer
        # (location_qualifier_value.term_id) for which oddly the schema
        # does not allow NULL.
        if feature.location_operator:
            #e.g. "join" (common),
            #or "order" (see Tests/GenBank/protein_refseq2.gb)
            location_id = self.adaptor.last_id('location')
            loc_qual_term_id = None # Not allowed in BioSQL v1.0.1
            sql = r"INSERT INTO location_qualifier_value" \
                  r"(location_id, term_id, value)" \
                  r"VALUES (%s, %s, %s)"
            self.adaptor.execute(sql, (location_id, loc_qual_term_id,
                                       feature.location_operator))
        """

    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):
        """Insert the (key, value) pair qualifiers relating to a feature (PRIVATE).

        Qualifiers should be a dictionary of the form:
            {key : [value1, value2]}
        """
        tag_ontology_id = self._get_ontology_id('Annotation Tags')
        for qualifier_key in qualifiers:
            # Treat db_xref qualifiers differently to sequence annotation
            # qualifiers by populating the seqfeature_dbxref and dbxref
            # tables.  Other qualifiers go into the seqfeature_qualifier_value
            # and (if new) term tables.
            if qualifier_key != 'db_xref':
                qualifier_key_id = self._get_term_id(qualifier_key,
                                                  ontology_id=tag_ontology_id)
                # now add all of the values to their table
                entries = qualifiers[qualifier_key]
                if not isinstance(entries, list):
                    # Could be a plain string, or an int or a float.
                    # However, we exect a list of strings here.
                    entries = [entries]
                for qual_value_rank in range(len(entries)):
                    qualifier_value = entries[qual_value_rank]
                    sql = r"INSERT INTO seqfeature_qualifier_value "\
                          r" (seqfeature_id, term_id, rank, value) VALUES"\
                          r" (%s, %s, %s, %s)"
                    self.adaptor.execute(sql, (seqfeature_id,
                                               qualifier_key_id,
                                               qual_value_rank + 1,
                                               qualifier_value))
            else:
                # The dbxref_id qualifier/value sets go into the dbxref table
                # as dbname, accession, version tuples, with dbxref.dbxref_id
                # being automatically assigned, and into the seqfeature_dbxref
                # table as seqfeature_id, dbxref_id, and rank tuples
                self._load_seqfeature_dbxref(qualifiers[qualifier_key],
                                             seqfeature_id)


    def _load_seqfeature_dbxref(self, dbxrefs, seqfeature_id):
        """Add database crossreferences of a SeqFeature to the database (PRIVATE).

            o dbxrefs           List, dbxref data from the source file in the
                                format <database>:<accession>

            o seqfeature_id     Int, the identifier for the seqfeature in the
                                seqfeature table

            Insert dbxref qualifier data for a seqfeature into the
            seqfeature_dbxref and, if required, dbxref tables.
            The dbxref_id qualifier/value sets go into the dbxref table
            as dbname, accession, version tuples, with dbxref.dbxref_id
            being automatically assigned, and into the seqfeature_dbxref
            table as seqfeature_id, dbxref_id, and rank tuples
        """
        # NOTE - In older versions of Biopython, we would map the GenBank
        # db_xref "name", for example "GI" to "GeneIndex", and give a warning
        # for any unknown terms.  This was a long term maintainance problem,
        # and differed from BioPerl and BioJava's implementation.  See bug 2405
        for rank, value in enumerate(dbxrefs):
            # Split the DB:accession format string at colons.  We have to
            # account for multiple-line and multiple-accession entries
            try:
                dbxref_data = value.replace(' ','').replace('\n','').split(':')
                db = dbxref_data[0]
                accessions = dbxref_data[1:]
            except:
                raise ValueError("Parsing of db_xref failed: '%s'" % value)
            # Loop over all the grabbed accessions, and attempt to fill the
            # table
            for accession in accessions:
                # Get the dbxref_id value for the dbxref data
                dbxref_id = self._get_dbxref_id(db, accession)
                # Insert the seqfeature_dbxref data
                self._get_seqfeature_dbxref(seqfeature_id, dbxref_id, rank+1)
        
    def _get_dbxref_id(self, db, accession):
        """ _get_dbxref_id(self, db, accession) -> Int

            o db          String, the name of the external database containing
                          the accession number

            o accession   String, the accession of the dbxref data

            Finds and returns the dbxref_id for the passed data.  The method
            attempts to find an existing record first, and inserts the data
            if there is no record.
        """
        # Check for an existing record
        sql = r'SELECT dbxref_id FROM dbxref WHERE dbname = %s ' \
              r'AND accession = %s'
        dbxref_id = self.adaptor.execute_and_fetch_col0(sql, (db, accession))
        # If there was a record, return the dbxref_id, else create the
        # record and return the created dbxref_id
        if dbxref_id:
            return dbxref_id[0]
        return self._add_dbxref(db, accession, 0)

    def _get_seqfeature_dbxref(self, seqfeature_id, dbxref_id, rank):
        """ Check for a pre-existing seqfeature_dbxref entry with the passed
            seqfeature_id and dbxref_id.  If one does not exist, insert new
            data

        """
        # Check for an existing record
        sql = r"SELECT seqfeature_id, dbxref_id FROM seqfeature_dbxref " \
              r"WHERE seqfeature_id = %s AND dbxref_id = %s"
        result = self.adaptor.execute_and_fetch_col0(sql, (seqfeature_id,
                                                           dbxref_id))
        # If there was a record, return without executing anything, else create
        # the record and return
        if result:
            return result
        return self._add_seqfeature_dbxref(seqfeature_id, dbxref_id, rank)

    def _add_seqfeature_dbxref(self, seqfeature_id, dbxref_id, rank):
        """ Insert a seqfeature_dbxref row and return the seqfeature_id and
            dbxref_id
        """
        sql = r'INSERT INTO seqfeature_dbxref ' \
              '(seqfeature_id, dbxref_id, rank) VALUES' \
              r'(%s, %s, %s)'
        self.adaptor.execute(sql, (seqfeature_id, dbxref_id, rank))
        return (seqfeature_id, dbxref_id)

    def _load_dbxrefs(self, record, bioentry_id):
        """Load any sequence level cross references into the database (PRIVATE).

        See table bioentry_dbxref."""
        for rank, value in enumerate(record.dbxrefs):
            # Split the DB:accession string at first colon.
            # We have to cope with things like:
            # "MGD:MGI:892" (db="MGD", accession="MGI:892")
            # "GO:GO:123" (db="GO", accession="GO:123")
            #
            # Annoyingly I have seen the NCBI use both the style
            # "GO:GO:123" and "GO:123" in different vintages.
            assert value.count("\n")==0
            try:
                db, accession = value.split(':',1)
                db = db.strip()
                accession = accession.strip()
            except:
                raise ValueError("Parsing of dbxrefs list failed: '%s'" % value)
            # Get the dbxref_id value for the dbxref data
            dbxref_id = self._get_dbxref_id(db, accession)
            # Insert the bioentry_dbxref  data
            self._get_bioentry_dbxref(bioentry_id, dbxref_id, rank+1)

    def _get_bioentry_dbxref(self, bioentry_id, dbxref_id, rank):
        """ Check for a pre-existing bioentry_dbxref entry with the passed
            seqfeature_id and dbxref_id.  If one does not exist, insert new
            data

        """
        # Check for an existing record
        sql = r"SELECT bioentry_id, dbxref_id FROM bioentry_dbxref " \
              r"WHERE bioentry_id = %s AND dbxref_id = %s"
        result = self.adaptor.execute_and_fetch_col0(sql, (bioentry_id,
                                                           dbxref_id))
        # If there was a record, return without executing anything, else create
        # the record and return
        if result:
            return result
        return self._add_bioentry_dbxref(bioentry_id, dbxref_id, rank)

    def _add_bioentry_dbxref(self, bioentry_id, dbxref_id, rank):
        """ Insert a bioentry_dbxref row and return the seqfeature_id and
            dbxref_id
        """
        sql = r'INSERT INTO bioentry_dbxref ' \
              '(bioentry_id,dbxref_id,rank) VALUES ' \
              '(%s, %s, %s)'
        self.adaptor.execute(sql, (bioentry_id, dbxref_id, rank))
        return (bioentry_id, dbxref_id)
            
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

