"""Implementations of Biopython-like Seq objects on top of BioSQL.

This allows retrival of items stored in a BioSQL database using
a biopython-like Seq interface.
"""
class DBSeq:  # This implements the biopython Seq interface
    def __init__(self, primary_id, adaptor, alphabet, start, length):
        self.primary_id = primary_id
        self.adaptor = adaptor
        self.alphabet = alphabet
        self._length = length
        self.start = start

    def __getattr__(self, name):
        if name == "data":
            return self.tostring()
        raise AttributeError, name
    
    def __len__(self):
        return self._length
    
    def __getitem__(self, i):
        if i < 0:
            if -i >= self._length:
                raise IndexError(i)
            i = i + self._length
        elif i >= self._length:
            raise IndexError(i)

        return self.adaptor.get_subseq_as_string(self.primary_id,
                                                 self.start + i,
                                                 self.start + i + 1)
    def __getslice__(self, i, j):
        i = max(i, 0)
        i = min(i, self._length)
        j = max(j, 0)
        j = min(j, self._length)
        if i >= j:
            return self.__class__(self.primary_id, self.adaptor, self.alphabet,
                                  self.start, 0)
        return self.__class__(self.primary_id, self.adaptor, self.alphabet,
                              self.start + i, j - i)
    def tostring(self):
        return self.adaptor.get_subseq_as_string(self.primary_id,
                                                 self.start,
                                                 self.start + self._length)

class DBInternalSeq:
    """Intermediate database object for retrieving sequence and other info.

    You should never need to deal directly with this object except when
    debugging. This is basically the equivalent of bioperl's PrimarySeq
    """
    def __init__(self, primary_id, adaptor):
        self.primary_id = primary_id
        self.adaptor = adaptor

        self.name, self.id, _length, self.description, self.moltype = \
                         self.adaptor.execute_one(
            """select en.display_id, en.accession, length(bs.biosequence_str),
                      en.description, bs.alphabet
               from bioentry en, biosequence bs
               where bs.bioentry_id = en.bioentry_id and
                     bs.bioentry_id = %s""",
            (self.primary_id,))

        self._length = int(_length)
        
    def __getattr__(self, name):
        if name[:1] == '_':
            raise AttributeError, name
        if name == "seq":
            moltype = self.moltype.upper()
            from Bio.Alphabet import IUPAC
            if moltype == "DNA":
                alphabet = IUPAC.unambiguous_dna
            elif moltype == "RNA":
                alphabet = IUPAC.unambiguous_rna
            elif moltype == "protein":
                alphabet = IUPAC.protein
            else:
                raise AssertionError("Unknown moltype: %s" % self.moltype)
            seq = DBSeq(self.primary_id, self.adaptor, alphabet,
                        0, self._length)
            self.seq = seq
            return seq
        f = getattr(self, "_get_" + name, None)
        if f is None:
            raise AttributeError, name
        return f()

    def __len__(self):
        return self._length

def _get_ontology_terms(ontology_name, bioentry_id, adaptor):
    """Retrieve ontology values associated with the given id and name.
    """
    sql = r"SELECT ontology_term_id FROM ontology_term " \
          r"WHERE term_name = %s" 
    id_info = adaptor.execute_and_fetchall(sql, (ontology_name,))
    if id_info is None:
        return None
    
    ontology_id = id_info[0][0]

    sql = r"SELECT qualifier_value FROM bioentry_qualifier_value " \
          r"WHERE bioentry_id = %s AND ontology_term_id = %s"
    values = adaptor.execute_and_fetch_col0(sql, (bioentry_id,
                                                  ontology_id))
    return values

class DBLink:
    def __init__(self, database, primary_id, optional_id = None,
                 version = None, comment = None):
        self.database = database
        self.primary_id = primary_id
        self.optional_id = optional_id
        self.version = version
        self.comment = comment

class Reference:
    def __init__(self, title, location, authors, medline, start, end):
        self.title = title
        self.location = location
        self.authors = authors
        self.medline = medline
        self.start = start
        self.end = end
        
import re
_is_title = re.compile(r"[A-Z][^A-Z]*$")
_is_lower = re.compile(r"[^A-Z]+$")

class Species:
    def __init__(self, classification, common_name = None, ncbi_id = None,
                 sub_species = None, organelle = None):
        self.classification = classification
        self.common_name = common_name
        self.ncbi_id = ncbi_id
        self.sub_species = sub_species
        self.organelle = organelle

    def __str__(self):
        terms = self.classification[:]
        terms.reverse()
        return " ".join(terms)

    def _assert_classification(self, classification):
        assert len(classification) >= 2, "classification must contain at least the species and genus"
        self._assert_species(classification[0])
        for x in classification[1:]:
            self._assert_name(x)
            
    def _assert_species(self, name):
        assert name, "species name must be defined"
        assert _is_lower.match(name), \
               "%r not allowed in a classification (species name)" % name

    def _assert_name(self, name):
        assert name, "classification name must be defined"
        assert _is_title.match(name), \
               "%r not allowed in a classification (not species name)" % name

    def binomial(self, full = 1):
        """return the 'Genus species' and optional subspecies (if requested)

        If there is a subspecies and full is true (the default), then
        this returns the string containing "Genus species subspecies".
        Otherwise it returns the string containing "Genus species".
        """
        if full and self.sub_species:
            return " ".join( (self.classification[1],
                              self.classification[0],
                              self.sub_species) )
        else:
            return " ".join( (self.classification[1],
                              self.classification[0]) )
            
    def __getattr__(self, name):
        if name == "species":
            return self.classification[0]
        elif name == "genus":
            return self.classification[1]

        raise AttributeError, name

    def __setattr__(self, name, val):
        if name == "species":
            self._assert_species(val)
            self.classification[0] = val
        elif name == "genus":
            self._assert_name(val)
            self.classification[1] = val
        elif name == "classification":
            self._assert_classification(val)
            self.__dict__["classification"] = val
        else:
            self.__dict__[name] = val


class Annotation:
    def __init__(self, adaptor, primary_id):
        self.adaptor = adaptor
        self.primary_id = primary_id

    def __getattr__(self, name):
        if name[:1] == '_':
            raise AttributeError, name
        f = getattr(self, "_get_" + name, None)
        if f is None:
            raise AttributeError, name
        return f()
    
    # functions to make this more like a dictionary
    def __getitem__(self, key):
        if key in ["comments", "dblinks", "references"]:
            return getattr(self, key)
        else:
            raise KeyError("Unexpected item: %s" % key)

    def has_key(self, key):
        if key in ["comments", "dblinks", "references"]:
            return 1
        else:
            return 0

    def __delitem__(self, key):
        # don't really get rid of things
        pass

    def _get_comments(self):
        comments = self.adaptor.execute_and_fetch_col0(
            """select comment_text from comment where bioentry_id = %s""",
            (self.primary_id,))
        self.comments = comments
        return comments
        
    def _get_dblinks(self):
        dblink_info = self.adaptor.execute_and_fetchall(
            "SELECT dbname, accession FROM bioentry_dblink, dbxref" \
            " WHERE bioentry_dblink.dbxref_id = dbxref.dbxref_id" \
            " AND bioentry_id = %s""",
            (self.primary_id,))
        dblinks =  [DBLink(database, primary_id) for (database, primary_id)
                                                     in dblink_info]
        self.dblinks = dblinks
        return dblinks

    def _get_references(self):
        results = self.adaptor.execute_and_fetchall(
            """select reference_id, reference_start,reference_end
                       from bioentry_reference
                       where bioentry_id = %s
                       order by reference_rank""",
            (self.primary_id,))

        references = []
        for reference_id, start, end in results:
            authors, title, location, medline = \
                     self.adaptor.execute_one(
                """select reference_authors, reference_title,
                          reference_location, reference_medline
                       from reference
                       where reference_id = %s""",
                (reference_id,))
            references.append( Reference(title, location, authors, medline,
                                         start, end) )
        self.references = references
        return references

# bioperl-db             -->        biopython
#  seqfeature_key                 SeqFeature.type
#  seqfeature_source                XXXX
#  seqfeature_qualifier_value     SeqFeature.qualifier.values()
#  seqfeature_qualifier           SeqFeature.qualifier.keys()
#  seqfeature_location            SeqFeature in SeqFeature.sub_features

def load_seq_features(adaptor, primary_id):
    from Bio import SeqFeature
    
    # Get the seqfeature id list
    sql = r"SELECT seqfeature_id, seqfeature_rank, ontology_term_id " \
          r"FROM seqfeature WHERE bioentry_id = %s"
    results = adaptor.execute_and_fetchall(sql, (primary_id,))

    seq_feature_list = []

    for seqfeature_id, seqfeature_rank, seqfeature_key_id in results:
        # Get the name
        sql = r"SELECT term_name FROM ontology_term " \
              r"WHERE ontology_term_id = %s"
        seqfeature_key = adaptor.execute_one(sql, (seqfeature_key_id,))[0]
        
        # Get its qualifiers
        qualifiers = {}
        sql = r"SELECT ot.term_name, sqv.qualifier_value " \
              r"FROM ontology_term ot, seqfeature_qualifier_value sqv " \
              r"WHERE ot.ontology_term_id = sqv.ontology_term_id AND " \
              r"sqv.seqfeature_id = %s"
        results = adaptor.execute_and_fetchall(sql, (seqfeature_id,))

        for key, value in results:
            if qualifiers.has_key(key):
                qualifiers[key].append(value)
            else:
                qualifiers[key] = [value]

        # Get its locations
        locations = []
        results = adaptor.execute_and_fetchall(
            """select seqfeature_location_id, seq_start, seq_end, seq_strand
                     from seqfeature_location
                     where seqfeature_id = %s""",
            (seqfeature_id,))

        # convert to Python standard form
        for location_id, start, end, strand in results:
            locations.append( (location_id, start-1, end, strand) )

        # Get any remote reference information
        remote_results = adaptor.execute_and_fetchall("""
          SELECT seqfeature_location_id, accession, version
            FROM seqfeature_location sfl, dbxref drf
            WHERE drf.dbxref_id = sfl.dbxref_id AND
                  sfl.seqfeature_id = %s""",
                                                      (seqfeature_id,))
        # Do the merge locally
        lookup = {}
        for location_id, accession, version in remote_results:
            lookup[location_id] = (accession, version)

        feature = SeqFeature.SeqFeature()
        feature.type = seqfeature_key
        feature.qualifiers = qualifiers
        if len(locations) == 0:
            pass
        elif len(locations) == 1:
            location_id, start, end, strand = locations[0]
            accession, version = lookup.get(location_id, (None, None))

            feature.location = SeqFeature.FeatureLocation(start, end)
            feature.strand = strand
            feature.ref = accession
        else:
            min_start = locations[0][1]
            max_end = locations[0][2]
            for location in locations:
                location_id, start, end, strand = location
                accession, version = lookup.get(location_id, (None, None))
                min_start = min(min_start, start)
                max_end = max(max_end, end)

                subfeature = SeqFeature.SeqFeature()
                subfeature.type = seqfeature_key 
                subfeature.location_operator = "join"
                subfeature.location = SeqFeature.FeatureLocation(start, end)
                subfeature.strand = strand
                subfeature.ref = accession
                feature.sub_features.append(subfeature)
            feature.sub_features.sort(lambda x, y: cmp(x.location.start, y.location.start))
            feature.location = SeqFeature.FeatureLocation(min_start, max_end)
            feature.strand = feature.sub_features[0].strand

        seq_feature_list.append( (seqfeature_rank, feature.location.start.position, feature) )

    # Primary sort is on the feature's rank
    #  .. then on the start position
    #  .. then arbitrary on the feature's id (SeqFeature has no __cmp__)
    seq_feature_list.sort()

    # Get just the SeqFeature
    return [x[2] for x in seq_feature_list]

class DBSeqRecord:
    """BioSQL equivalent of the biopython SeqRecord object.
    """
    _forward_getattr = ['seq', 'description', 'name', 'id']

    def __init__(self, adaptor, primary_id):
        self.adaptor = adaptor
        self.primary_id = primary_id

        self.version, _length, self.division = \
                      self.adaptor.execute_one(
            """select en.entry_version, length(bs.biosequence_str), bs.division
                    from bioentry en, biosequence bs
                    where bs.bioentry_id = en.bioentry_id and
                          bs.bioentry_id = %s""",
            (self.primary_id,))
        self._length = int(_length)
        
    def __len__(self):
        return self._length

    def __getattr__(self, name):
        if name[:1] == "_":
            raise AttributeError, name
        if name in self._forward_getattr:
            return getattr(self.primary_seq, name)
        f = getattr(self, "_get_" + name, None)
        if f is None:
            raise AttributeError, name
        return f()

    def _get_primary_seq(self):
        primary_seq = DBInternalSeq(self.primary_id, self.adaptor)
        self.primary_seq = primary_seq
        return primary_seq

    # -- SeqRecord compatilibity attributes

    def _get_annotations(self):
        # XXX This should return a dictionary of annotations ala
        # the actual SeqRecord object
        annotation = Annotation(self.adaptor, self.primary_id)
        self.annotation = annotation
        return annotation

    def _get_features(self):
        features = load_seq_features(self.adaptor, self.primary_id)
        self.seq_features = features
        return self.seq_features
   
    # -- BioSQL only attributes
   
    def _get_dates(self):
        self.dates = _get_ontology_terms("date", self.primary_id, self.adaptor)
        return self.dates

    def _get_species(self):
        full_lineage, common_name, ncbi_id = self.adaptor.execute_one(
            """select tx.full_lineage, tx.common_name, tx.ncbi_taxon_id
                           from taxon tx, bioentry be
                           where tx.taxon_id = be.taxon_id and
                                 be.bioentry_id = %s""",
            (self.primary_id,))
        terms = full_lineage.split(":")
        species = Species(terms, common_name, ncbi_id)
        return species

    def _get_seq_version(self):
        version = self.accession + "." + self.version
        return version

    def _get_keywords(self):
        keywords = _get_ontology_term('Keywords', self.primary_id,
                                      self.adaptor)
        self.keywords = keywords
        return keywords
