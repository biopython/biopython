"""Implementations of Biopython-like Seq objects on top of BioSQL.

This allows retrival of items stored in a BioSQL database using
a biopython-like Seq interface.
"""

from Bio import SeqFeature

class DBSeq:  # This implements the biopython Seq interface
    def __init__(self, primary_id, adaptor, alphabet, start, length):
        self.primary_id = primary_id
        self.adaptor = adaptor
        self.alphabet = alphabet
        self._length = length
        self.start = start

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

    data = property(tostring)

def _retrieve_seq(adaptor, primary_id):
    seqs = adaptor.execute_and_fetchall(
        "SELECT alphabet, length(seq) FROM biosequence" \
        " WHERE bioentry_id = %s", (primary_id,))
    if seqs:
        moltype, length = seqs[0]
        from Bio.Alphabet import IUPAC
        if moltype == "dna":
            alphabet = IUPAC.unambiguous_dna
        elif moltype == "rna":
            alphabet = IUPAC.unambiguous_rna
        elif moltype == "protein":
            alphabet = IUPAC.protein
        else:
            raise AssertionError("Unknown moltype: %s" % moltype)
        seq = DBSeq(primary_id, adaptor, alphabet, 0, length)
        return seq
    else:
        return None

def _retrieve_dbxrefs(adaptor, primary_id):
    _dbxrefs = []
    dbxrefs = adaptor.execute_and_fetchall(
        "SELECT dbname, accession, version" \
        " FROM bioentry_dbxref join dbxref using (dbxref_id)" \
        " WHERE bioentry_id = %s" \
        " ORDER BY rank", (primary_id,))
    for dbname, accession, version in dbxrefs:
        if version and version != "0":
            v = "%s.%s" % (accession, version)
        else:
            v = accession
        _dbxrefs.append((dbname, v))
    return _dbxrefs

def _retrieve_features(adaptor, primary_id):
    sql = "SELECT seqfeature_id, type.name, rank" \
          " FROM seqfeature join term type on (type_term_id = type.term_id)" \
          " WHERE bioentry_id = %s" \
          " ORDER BY rank"
    results = adaptor.execute_and_fetchall(sql, (primary_id,))
    seq_feature_list = []
    for seqfeature_id, seqfeature_type, seqfeature_rank in results:
        # Get qualifiers
        qvs = adaptor.execute_and_fetchall(
            "SELECT name, value" \
            " FROM seqfeature_qualifier_value  join term using (term_id)" \
            " WHERE seqfeature_id = %s", (seqfeature_id,))
        qualifiers = {}
        for qv_name, qv_value in qvs:
            qualifiers.setdefault(qv_name, []).append(qv_value)
        # Get locations
        results = adaptor.execute_and_fetchall(
            "SELECT location_id, start_pos, end_pos, strand" \
            " FROM location" \
            " WHERE seqfeature_id = %s" \
            " ORDER BY rank", (seqfeature_id,))
        locations = []
        # convert to Python standard form
        for location_id, start, end, strand in results:
            if start:
                start -= 1
            locations.append( (location_id, start, end, strand) )
        # Get possible remote reference information
        remote_results = adaptor.execute_and_fetchall(
            "SELECT location_id, dbname, accession, version" \
            " FROM location join dbxref using (dbxref_id)" \
            " WHERE seqfeature_id = %s", (seqfeature_id,))
        lookup = {}
        for location_id, dbname, accession, version in remote_results:
            if version and version != "0":
                v = "%s.%s" % (accession, version)
            else:
                v = accession
            lookup[location_id] = (dbname, v)
        
        feature = SeqFeature.SeqFeature(type = seqfeature_type)
        feature.qualifiers = qualifiers
        if len(locations) == 0:
            pass
        elif len(locations) == 1:
            location_id, start, end, strand = locations[0]
            dbname, version = lookup.get(location_id, (None, None))

            feature.location = SeqFeature.FeatureLocation(start, end)
            feature.strand = strand
            feature.ref_db = dbname
            feature.ref = version
        else:
            min_start = locations[0][1]
            max_end = locations[0][2]
            sub_feature_list = []       # (start, sub feature) for sorting
            for location in locations:
                location_id, start, end, strand = location
                dbname, version = lookup.get(location_id, (None, None))
                min_start = min(min_start, start)
                max_end = max(max_end, end)

                subfeature = SeqFeature.SeqFeature()
                subfeature.type = seqfeature_type 
                subfeature.location_operator = "join"
                subfeature.location = SeqFeature.FeatureLocation(start, end)
                subfeature.strand = strand
                subfeature.ref_db = dbname
                subfeature.ref = version
                sub_feature_list.append((start, subfeature))
            sub_feature_list.sort()
            feature.sub_features = [sub_feature[1]
                                    for sub_feature in sub_feature_list]
            feature.location = SeqFeature.FeatureLocation(min_start, max_end)
            feature.strand = feature.sub_features[0].strand

        seq_feature_list.append(
            (seqfeature_rank, feature.location.start.position, feature) )

    # Primary sort is on the feature's rank
    #  .. then on the start position
    #  .. then arbitrary on the feature's id (SeqFeature has no __cmp__)
    seq_feature_list.sort()

    # Get just the SeqFeature
    return [x[2] for x in seq_feature_list]

def _retrieve_annotations(adaptor, primary_id, taxon_id):
    annotations = {}
    annotations.update(_retrieve_qualifier_value(adaptor, primary_id))
    annotations.update(_retrieve_reference(adaptor, primary_id))
    annotations.update(_retrieve_dbxref(adaptor, primary_id))
    annotations.update(_retrieve_taxon(adaptor, primary_id, taxon_id))
    return annotations

def _retrieve_qualifier_value(adaptor, primary_id):
    qvs = adaptor.execute_and_fetchall(
        "SELECT name, value" \
        " FROM bioentry_qualifier_value JOIN term USING (term_id)" \
        " WHERE bioentry_id = %s" \
        " ORDER BY term_id, rank", (primary_id,))
    qualifiers = {}
    for name, value in qvs:
        if name == "keyword": name = "keywords"
        elif name == "date_changed": name = "dates"
        elif name == "secondary_accession": name = "accessions"
        qualifiers.setdefault(name, []).append(value)
    return qualifiers

def _retrieve_reference(adaptor, primary_id):
    # XXX dbxref_qualifier_value
 
    refs = adaptor.execute_and_fetchall(
        "SELECT start_pos, end_pos, " \
        " location, title, authors," \
        " dbname, accession" \
        " FROM bioentry_reference" \
        " JOIN reference USING (reference_id)" \
        " LEFT JOIN dbxref USING (dbxref_id)" \
        " WHERE bioentry_id = %s" \
        " ORDER BY reference_id, rank", (primary_id,))
    references = []
    for start, end, location, title, authors, dbname, accession in refs:
        reference = SeqFeature.Reference()
        if start: start -= 1
        reference.location = [SeqFeature.FeatureLocation(start, end)]
        reference.authors = authors
        reference.title = title
        reference.journal = location
        if dbname == 'PUBMED':
            reference.pubmed_id = accession
        elif dbname == 'MEDLINE':
            reference.medline_id = accession
        references.append(reference)
    return {'references': references}

def _retrieve_dbxref(adaptor, primary_id):
    # XXX dbxref_qualifier_value

    refs = adaptor.execute_and_fetchall(
        "SELECT dbname, accession, version" \
        " FROM bioentry_dbxref JOIN dbxref USING (dbxref_id)" \
        " WHERE bioentry_id = %s" \
        " ORDER BY dbxref_id, rank", (primary_id,))
    dbxrefs = []
    for dbname, accession, version in refs:
        if version and version != "0":
            v = "%s.%s" % (accession, version)
        else:
            v = accession
        dbxrefs.append((dbname, v))
    return {'cross_references': dbxrefs}

def _retrieve_taxon(adaptor, primary_id, taxon_id):
    a = {}
    common_names = adaptor.execute_and_fetch_col0(
        "SELECT name FROM taxon_name WHERE taxon_id = %s" \
        " AND name_class = 'genbank common name'", (taxon_id,))
    if common_names:
        a['source'] = common_names[0]
    scientific_names = adaptor.execute_and_fetch_col0(
        "SELECT name FROM taxon_name WHERE taxon_id = %s" \
        " AND name_class = 'scientific name'", (taxon_id,))
    if scientific_names:
        a['organism'] = scientific_names[0]
    ncbi_taxids = adaptor.execute_and_fetch_col0(
        "SELECT ncbi_taxon_id FROM taxon WHERE taxon_id = %s", (taxon_id,))
    if ncbi_taxids and ncbi_taxids[0] and ncbi_taxids[0] != "0":
        a['ncbi_taxid'] = ncbi_taxids[0]

    # XXX 'no rank' or not?
    # If we keep them, perhaps add " AND parent.parent_taxon_id <> 1"
    taxonomy = adaptor.execute_and_fetch_col0(
        "SELECT taxon_name.name" \
        " FROM taxon t, taxon parent, taxon_name" \
        " WHERE parent.taxon_id=taxon_name.taxon_id" \
        " AND t.left_value BETWEEN parent.left_value AND parent.right_value" \
        " AND name_class='scientific name'" \
        " AND parent.node_rank <> 'no rank'" \
        " AND t.taxon_id = %s" \
        " ORDER BY parent.left_value", (taxon_id,))
    if taxonomy:
        a['taxonomy'] = taxonomy
    return a

class DBSeqRecord(object):
    """BioSQL equivalent of the biopython SeqRecord object.
    """

    def __init__(self, adaptor, primary_id):
        self._adaptor = adaptor
        self._primary_id = primary_id

        (self._biodatabase_id, self._taxon_id, self.name,
         accession, version, self._identifier,
         self._division, self.description) = self._adaptor.execute_one(
            "SELECT biodatabase_id, taxon_id, name, accession, version," \
            " identifier, division, description" \
            " FROM bioentry" \
            " WHERE bioentry_id = %s", (self._primary_id,))
        if version and version != "0":
            self.id = "%s.%s" % (accession, version)
        else:
            self.id = accession

    def __get_seq(self):
        if not hasattr(self, "_seq"):
            self._seq = _retrieve_seq(self._adaptor, self._primary_id)
        return self._seq
    def __set_seq(self, seq): self._seq = seq
    def __del_seq(self):      del self._seq
    seq = property(__get_seq, __set_seq, __del_seq, "Seq object")

    def __get_dbxrefs(self):
        if not hasattr(self,"_dbxrefs"):
            self._dbxrefs = _retrieve_dbxrefs(self._adaptor, self._primary_id)
        return self._dbxrefs
    def __set_dbxrefs(self, dbxrefs): self._dbxrefs = dbxrefs
    def __del_dbxrefs(self):      del self._dbxrefs
    dbxrefs = property(__get_dbxrefs, __set_dbxrefs, __del_dbxrefs, "Dbxrefs")

    def __get_features(self):
        if not hasattr(self, "_features"):
            self._features = _retrieve_features(self._adaptor,
                                                self._primary_id)
        return self._features
    def __set_features(self, features): self._features = features
    def __del_features(self):      del self._features
    features = property(__get_features, __set_features, __del_features,
                        "Features")

    def __get_annotations(self):
        if not hasattr(self, "_annotations"):
            self._annotations = _retrieve_annotations(self._adaptor,
                                                      self._primary_id,
                                                      self._taxon_id)
            if self._identifier:
                self._annotations["gi"] = self._identifier
            if self._division:
                self._annotations["data_file_division"] = self._division
        return self._annotations
    def __set_annotations(self, annotations): self._annotations = annotations
    def __del_annotations(self): del self._annotations
    annotations = property(__get_annotations, __set_annotations,
                           __del_annotations, "Annotations")
