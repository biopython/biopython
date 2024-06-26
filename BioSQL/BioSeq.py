# Copyright 2002 by Andrew Dalke.  All rights reserved.
# Revisions 2007-2016 copyright by Peter Cock.  All rights reserved.
# Revisions 2008-2009 copyright by Cymon J. Cox.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Note that BioSQL (including the database schema and scripts) is
# available and licensed separately.  Please consult www.biosql.org
"""Implementations of Biopython-like Seq objects on top of BioSQL.

This allows retrieval of items stored in a BioSQL database using
a biopython-like SeqRecord and Seq interface.

Note: Currently we do not support recording per-letter-annotations
(like quality scores) in BioSQL.
"""

from typing import List, Optional

from Bio.Seq import Seq, SequenceDataAbstractBaseClass
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from Bio import SeqFeature


class _BioSQLSequenceData(SequenceDataAbstractBaseClass):
    """Retrieves sequence data from a BioSQL database (PRIVATE)."""

    __slots__ = ("primary_id", "adaptor", "_length", "start")

    def __init__(self, primary_id, adaptor, start=0, length=0):
        """Create a new _BioSQLSequenceData object referring to a BioSQL entry.

        You wouldn't normally create a _BioSQLSequenceData object yourself,
        this is done for you when retrieving a DBSeqRecord object from the
        database, which creates a Seq object using a _BioSQLSequenceData
        instance as the data provider.
        """
        self.primary_id = primary_id
        self.adaptor = adaptor
        self._length = length
        self.start = start
        super().__init__()

    def __len__(self):
        """Return the length of the sequence."""
        return self._length

    def __getitem__(self, key):
        """Return a subsequence as a bytes or a _BioSQLSequenceData object."""
        if isinstance(key, slice):
            start, end, step = key.indices(self._length)
            size = len(range(start, end, step))
            if size == 0:
                return b""
        else:
            # Return a single letter as an integer (consistent with bytes)
            i = key
            if i < 0:
                i += self._length
                if i < 0:
                    raise IndexError(key)
            elif i >= self._length:
                raise IndexError(key)
            c = self.adaptor.get_subseq_as_string(
                self.primary_id, self.start + i, self.start + i + 1
            )
            return ord(c)

        if step == 1:
            if start == 0 and size == self._length:
                # Return the full sequence as bytes
                sequence = self.adaptor.get_subseq_as_string(
                    self.primary_id, self.start, self.start + self._length
                )
                return sequence.encode("ASCII")
            else:
                # Return a _BioSQLSequenceData with the start and end adjusted
                return _BioSQLSequenceData(
                    self.primary_id, self.adaptor, self.start + start, size
                )
        else:
            # Will have to extract the sequence because of the stride
            full = self.adaptor.get_subseq_as_string(
                self.primary_id, self.start + start, self.start + end
            )
            return full[::step].encode("ASCII")


def _retrieve_seq_len(adaptor, primary_id):
    # The database schema ensures there will be only one matching row
    seqs = adaptor.execute_and_fetchall(
        "SELECT length FROM biosequence WHERE bioentry_id = %s", (primary_id,)
    )
    if not seqs:
        return None
    if len(seqs) != 1:
        raise ValueError(f"Expected 1 response, got {len(seqs)}.")
    (given_length,) = seqs[0]
    return int(given_length)


def _retrieve_seq(adaptor, primary_id):
    # The database schema ensures there will be only one matching
    # row in the table.

    # If an undefined sequence was recorded, seq will be NULL,
    # but length will be populated.  This means length(seq)
    # will return None.
    seqs = adaptor.execute_and_fetchall(
        "SELECT alphabet, length, length(seq) FROM biosequence WHERE bioentry_id = %s",
        (primary_id,),
    )
    if not seqs:
        return
    if len(seqs) != 1:
        raise ValueError(f"Expected 1 response, got {len(seqs)}.")
    moltype, given_length, length = seqs[0]

    try:
        length = int(length)
        given_length = int(given_length)
        if length != given_length:
            raise ValueError(
                f"'length' differs from sequence length, {given_length}, {length}"
            )
        have_seq = True
    except TypeError:
        if length is not None:
            raise ValueError(f"Expected 'length' to be 'None', got {length}.")
        seqs = adaptor.execute_and_fetchall(
            "SELECT alphabet, length, seq FROM biosequence WHERE bioentry_id = %s",
            (primary_id,),
        )
        if len(seqs) != 1:
            raise ValueError(f"Expected 1 response, got {len(seqs)}.")
        moltype, given_length, seq = seqs[0]
        if seq:
            raise ValueError(f"Expected 'seq' to have a falsy value, got {seq}.")
        length = int(given_length)
        have_seq = False
        del seq
    del given_length

    if have_seq:
        data = _BioSQLSequenceData(primary_id, adaptor, start=0, length=length)
        return Seq(data)
    else:
        return Seq(None, length=length)


def _retrieve_dbxrefs(adaptor, primary_id):
    """Retrieve the database cross references for the sequence (PRIVATE)."""
    _dbxrefs = []
    dbxrefs = adaptor.execute_and_fetchall(
        "SELECT dbname, accession, version"
        " FROM bioentry_dbxref join dbxref using (dbxref_id)"
        " WHERE bioentry_id = %s"
        ' ORDER BY "rank"',
        (primary_id,),
    )
    for dbname, accession, version in dbxrefs:
        if version and version != "0":
            v = f"{accession}.{version}"
        else:
            v = accession
        _dbxrefs.append(f"{dbname}:{v}")
    return _dbxrefs


def _retrieve_features(adaptor, primary_id):
    sql = (
        'SELECT seqfeature_id, type.name, "rank"'
        " FROM seqfeature join term type on (type_term_id = type.term_id)"
        " WHERE bioentry_id = %s"
        ' ORDER BY "rank"'
    )
    results = adaptor.execute_and_fetchall(sql, (primary_id,))
    seq_feature_list = []
    for seqfeature_id, seqfeature_type, seqfeature_rank in results:
        # Get qualifiers [except for db_xref which is stored separately]
        qvs = adaptor.execute_and_fetchall(
            "SELECT name, value"
            " FROM seqfeature_qualifier_value  join term using (term_id)"
            " WHERE seqfeature_id = %s"
            ' ORDER BY "rank"',
            (seqfeature_id,),
        )
        qualifiers = {}
        for qv_name, qv_value in qvs:
            qualifiers.setdefault(qv_name, []).append(qv_value)
        # Get db_xrefs [special case of qualifiers]
        qvs = adaptor.execute_and_fetchall(
            "SELECT dbxref.dbname, dbxref.accession"
            " FROM dbxref join seqfeature_dbxref using (dbxref_id)"
            " WHERE seqfeature_dbxref.seqfeature_id = %s"
            ' ORDER BY "rank"',
            (seqfeature_id,),
        )
        for qv_name, qv_value in qvs:
            value = f"{qv_name}:{qv_value}"
            qualifiers.setdefault("db_xref", []).append(value)
        # Get locations
        results = adaptor.execute_and_fetchall(
            "SELECT location_id, start_pos, end_pos, strand"
            " FROM location"
            " WHERE seqfeature_id = %s"
            ' ORDER BY "rank"',
            (seqfeature_id,),
        )
        locations = []
        # convert to Python standard form
        # Convert strand = 0 to strand = None
        # re: comment in Loader.py:
        # Biopython uses None when we don't know strand information but
        # BioSQL requires something (non null) and sets this as zero
        # So we'll use the strand or 0 if Biopython spits out None
        for location_id, start, end, strand in results:
            if start:
                start -= 1
            if strand == 0:
                strand = None
            if strand not in (+1, -1, None):
                raise ValueError(
                    "Invalid strand %s found in database for "
                    "seqfeature_id %s" % (strand, seqfeature_id)
                )
            if start is not None and end is not None and end < start:
                import warnings
                from Bio import BiopythonWarning

                warnings.warn(
                    "Inverted location start/end (%i and %i) for "
                    "seqfeature_id %s" % (start, end, seqfeature_id),
                    BiopythonWarning,
                )

            # For SwissProt unknown positions (?)
            if start is None:
                start = SeqFeature.UnknownPosition()
            if end is None:
                end = SeqFeature.UnknownPosition()

            locations.append((location_id, start, end, strand))
        # Get possible remote reference information
        remote_results = adaptor.execute_and_fetchall(
            "SELECT location_id, dbname, accession, version"
            " FROM location join dbxref using (dbxref_id)"
            " WHERE seqfeature_id = %s",
            (seqfeature_id,),
        )
        lookup = {}
        for location_id, dbname, accession, version in remote_results:
            if version and version != "0":
                v = f"{accession}.{version}"
            else:
                v = accession
            # subfeature remote location db_ref are stored as a empty string
            # when not present
            if dbname == "":
                dbname = None
            lookup[location_id] = (dbname, v)

        feature = SeqFeature.SeqFeature(type=seqfeature_type)
        # Store the key as a private property
        feature._seqfeature_id = seqfeature_id
        feature.qualifiers = qualifiers
        if len(locations) == 0:
            pass
        elif len(locations) == 1:
            location_id, start, end, strand = locations[0]
            # See Bug 2677, we currently don't record the location_operator
            # For consistency with older versions Biopython, default to "".
            feature.location_operator = _retrieve_location_qualifier_value(
                adaptor, location_id
            )
            dbname, version = lookup.get(location_id, (None, None))
            feature.location = SeqFeature.SimpleLocation(start, end)
            feature.location.strand = strand
            feature.location.ref_db = dbname
            feature.location.ref = version
        else:
            locs = []
            for location in locations:
                location_id, start, end, strand = location
                dbname, version = lookup.get(location_id, (None, None))
                locs.append(
                    SeqFeature.SimpleLocation(
                        start, end, strand=strand, ref=version, ref_db=dbname
                    )
                )
            # Locations are typically in biological in order (see negative
            # strands below), but because of remote locations for
            # sub-features they are not necessarily in numerical order:
            strands = {_.strand for _ in locs}
            if len(strands) == 1 and -1 in strands:
                # Evil hack time for backwards compatibility
                # TODO - Check if BioPerl and (old) Biopython did the same,
                # we may have an existing incompatibility lurking here...
                locs = locs[::-1]
            feature.location = SeqFeature.CompoundLocation(locs, "join")
            # TODO - See Bug 2677 - we don't yet record location operator,
            # so for consistency with older versions of Biopython default
            # to assuming its a join.
        seq_feature_list.append(feature)
    return seq_feature_list


def _retrieve_location_qualifier_value(adaptor, location_id):
    value = adaptor.execute_and_fetch_col0(
        "SELECT value FROM location_qualifier_value WHERE location_id = %s",
        (location_id,),
    )
    try:
        return value[0]
    except IndexError:
        return ""


def _retrieve_annotations(adaptor, primary_id, taxon_id):
    annotations = {}
    annotations.update(_retrieve_alphabet(adaptor, primary_id))
    annotations.update(_retrieve_qualifier_value(adaptor, primary_id))
    annotations.update(_retrieve_reference(adaptor, primary_id))
    annotations.update(_retrieve_taxon(adaptor, primary_id, taxon_id))
    annotations.update(_retrieve_comment(adaptor, primary_id))
    return annotations


def _retrieve_alphabet(adaptor, primary_id):
    results = adaptor.execute_and_fetchall(
        "SELECT alphabet FROM biosequence WHERE bioentry_id = %s", (primary_id,)
    )
    if len(results) != 1:
        raise ValueError(f"Expected 1 response, got {len(results)}.")
    alphabets = results[0]
    if len(alphabets) != 1:
        raise ValueError(f"Expected 1 alphabet in response, got {len(alphabets)}.")
    alphabet = alphabets[0]
    if alphabet == "dna":
        molecule_type = "DNA"
    elif alphabet == "rna":
        molecule_type = "RNA"
    elif alphabet == "protein":
        molecule_type = "protein"
    else:
        molecule_type = None
    if molecule_type is not None:
        return {"molecule_type": molecule_type}
    else:
        return {}


def _retrieve_qualifier_value(adaptor, primary_id):
    qvs = adaptor.execute_and_fetchall(
        "SELECT name, value"
        " FROM bioentry_qualifier_value JOIN term USING (term_id)"
        " WHERE bioentry_id = %s"
        ' ORDER BY "rank"',
        (primary_id,),
    )
    qualifiers = {}
    for name, value in qvs:
        if name == "keyword":
            name = "keywords"
        # See handling of "date" in Loader.py
        elif name == "date_changed":
            name = "date"
        elif name == "secondary_accession":
            name = "accessions"
        qualifiers.setdefault(name, []).append(value)
    return qualifiers


def _retrieve_reference(adaptor, primary_id):
    # XXX dbxref_qualifier_value

    refs = adaptor.execute_and_fetchall(
        "SELECT start_pos, end_pos, "
        " location, title, authors,"
        " dbname, accession"
        " FROM bioentry_reference"
        " JOIN reference USING (reference_id)"
        " LEFT JOIN dbxref USING (dbxref_id)"
        " WHERE bioentry_id = %s"
        ' ORDER BY "rank"',
        (primary_id,),
    )
    references = []
    for start, end, location, title, authors, dbname, accession in refs:
        reference = SeqFeature.Reference()
        # If the start/end are missing, reference.location is an empty list
        if (start is not None) or (end is not None):
            if start is not None:
                start -= 1  # python counting
            reference.location = [SeqFeature.SimpleLocation(start, end)]
        # Don't replace the default "" with None.
        if authors:
            reference.authors = authors
        if title:
            reference.title = title
        reference.journal = location
        if dbname == "PUBMED":
            reference.pubmed_id = accession
        elif dbname == "MEDLINE":
            reference.medline_id = accession
        references.append(reference)
    if references:
        return {"references": references}
    else:
        return {}


def _retrieve_taxon(adaptor, primary_id, taxon_id):
    a = {}
    common_names = adaptor.execute_and_fetch_col0(
        "SELECT name FROM taxon_name WHERE taxon_id = %s"
        " AND name_class = 'genbank common name'",
        (taxon_id,),
    )
    if common_names:
        a["source"] = common_names[0]
    scientific_names = adaptor.execute_and_fetch_col0(
        "SELECT name FROM taxon_name WHERE taxon_id = %s"
        " AND name_class = 'scientific name'",
        (taxon_id,),
    )
    if scientific_names:
        a["organism"] = scientific_names[0]
    ncbi_taxids = adaptor.execute_and_fetch_col0(
        "SELECT ncbi_taxon_id FROM taxon WHERE taxon_id = %s", (taxon_id,)
    )
    if ncbi_taxids and ncbi_taxids[0] and ncbi_taxids[0] != "0":
        a["ncbi_taxid"] = ncbi_taxids[0]

    # Old code used the left/right values in the taxon table to get the
    # taxonomy lineage in one SQL command.  This was actually very slow,
    # and would fail if the (optional) left/right values were missing.
    #
    # The following code is based on a contribution from Eric Gibert, and
    # relies on the taxon table's parent_taxon_id field only (ignoring the
    # optional left/right values).  This means that it has to make a
    # separate SQL query for each entry in the lineage, but it does still
    # appear to be *much* faster.  See Bug 2494.
    taxonomy = []
    while taxon_id:
        name, rank, parent_taxon_id = adaptor.execute_one(
            "SELECT taxon_name.name, taxon.node_rank, taxon.parent_taxon_id"
            " FROM taxon, taxon_name"
            " WHERE taxon.taxon_id=taxon_name.taxon_id"
            " AND taxon_name.name_class='scientific name'"
            " AND taxon.taxon_id = %s",
            (taxon_id,),
        )
        if taxon_id == parent_taxon_id:
            # If the taxon table has been populated by the BioSQL script
            # load_ncbi_taxonomy.pl this is how top parent nodes are stored.
            # Personally, I would have used a NULL parent_taxon_id here.
            break

        taxonomy.insert(0, name)
        taxon_id = parent_taxon_id

    if taxonomy:
        a["taxonomy"] = taxonomy
    return a


def _retrieve_comment(adaptor, primary_id):
    qvs = adaptor.execute_and_fetchall(
        'SELECT comment_text FROM comment WHERE bioentry_id=%s ORDER BY "rank"',
        (primary_id,),
    )
    comments = [comm[0] for comm in qvs]
    # Don't want to add an empty list...
    if comments:
        return {"comment": comments}
    else:
        return {}


class DBSeqRecord(SeqRecord):
    """BioSQL equivalent of the Biopython SeqRecord object."""

    def __init__(self, adaptor, primary_id):
        """Create a DBSeqRecord object.

        Arguments:
         - adaptor - A BioSQL.BioSeqDatabase.Adaptor object
         - primary_id - An internal integer ID used by BioSQL

        You wouldn't normally create a DBSeqRecord object yourself,
        this is done for you when using a BioSeqDatabase object
        """
        self._adaptor = adaptor
        self._primary_id = primary_id

        (
            self._biodatabase_id,
            self._taxon_id,
            self.name,
            accession,
            version,
            self._identifier,
            self._division,
            self.description,
        ) = self._adaptor.execute_one(
            "SELECT biodatabase_id, taxon_id, name, accession, version,"
            " identifier, division, description"
            " FROM bioentry"
            " WHERE bioentry_id = %s",
            (self._primary_id,),
        )
        if version and version != "0":
            self.id = f"{accession}.{version}"
        else:
            self.id = accession
        # We don't yet record any per-letter-annotations in the
        # BioSQL database, but we should set this property up
        # for completeness (and the __str__ method).
        # We do NOT want to load the sequence from the DB here!
        length = _retrieve_seq_len(adaptor, primary_id)
        self._per_letter_annotations = _RestrictedDict(length=length)

    def __get_seq(self):
        if not hasattr(self, "_seq"):
            self._seq = _retrieve_seq(self._adaptor, self._primary_id)
        return self._seq

    def __set_seq(self, seq):
        # TODO - Check consistent with self._per_letter_annotations
        self._seq = seq

    def __del_seq(self):
        del self._seq

    seq = property(__get_seq, __set_seq, __del_seq, "Seq object")

    @property
    def dbxrefs(self) -> List[str]:
        """Database cross references."""
        if not hasattr(self, "_dbxrefs"):
            self._dbxrefs = _retrieve_dbxrefs(self._adaptor, self._primary_id)
        return self._dbxrefs

    @dbxrefs.setter
    def dbxrefs(self, value: List[str]) -> None:
        self._dbxrefs = value

    @dbxrefs.deleter
    def dbxrefs(self) -> None:
        del self._dbxrefs

    def __get_features(self):
        if not hasattr(self, "_features"):
            self._features = _retrieve_features(self._adaptor, self._primary_id)
        return self._features

    def __set_features(self, features):
        self._features = features

    def __del_features(self):
        del self._features

    features = property(__get_features, __set_features, __del_features, "Features")

    @property
    def annotations(self) -> SeqRecord._AnnotationsDict:
        """Annotations."""
        if not hasattr(self, "_annotations"):
            self._annotations = _retrieve_annotations(
                self._adaptor, self._primary_id, self._taxon_id
            )
            if self._identifier:
                self._annotations["gi"] = self._identifier
            if self._division:
                self._annotations["data_file_division"] = self._division
        return self._annotations

    @annotations.setter
    def annotations(self, value: Optional[SeqRecord._AnnotationsDict]) -> None:
        if value:
            self._annotations = value
        else:
            self._annotations = {}

    @annotations.deleter
    def annotations(self) -> None:
        del self._annotations
