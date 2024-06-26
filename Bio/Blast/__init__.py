# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# Revisions 2023 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Code to parse and store BLAST XML output, and to invoke the NCBI BLAST web server.

This module provides code to parse and store BLAST XML output, following its
definition in the associated BLAST XML DTD file:
https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd

This module also provides code to invoke the BLAST web server provided by NCBI.
https://blast.ncbi.nlm.nih.gov/

Variables:

    - email        Set the Blast email parameter (default is None).
    - tool         Set the Blast tool parameter (default is ``biopython``).

"""

import io
import textwrap
import time
import warnings
from collections import UserList
from urllib.parse import urlencode
from urllib.request import build_opener
from urllib.request import HTTPBasicAuthHandler
from urllib.request import HTTPPasswordMgrWithDefaultRealm
from urllib.request import install_opener
from urllib.request import Request
from urllib.request import urlopen
from xml.parsers import expat

import numpy as np

from Bio import BiopythonWarning
from Bio import StreamModeError
from Bio._utils import function_with_previous
from Bio.Align import Alignment
from Bio.Align import Alignments
from Bio.Blast import _writers

email = None
tool = "biopython"


NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

BLOCK = 2048  # default block size from expat


class NotXMLError(ValueError):
    """Failed to parse file as XML."""

    def __init__(self, message):
        """Initialize the class."""
        self.msg = message

    def __str__(self):
        """Return a string summary of the exception."""
        return (
            "Failed to parse the XML data (%s). Please make sure that the input data "
            "are in XML format." % self.msg
        )


class CorruptedXMLError(ValueError):
    """Corrupted XML."""

    def __init__(self, message):
        """Initialize the class."""
        self.msg = message

    def __str__(self):
        """Return a string summary of the exception."""
        return (
            "Failed to parse the XML data (%s). Please make sure that the input data "
            "are not corrupted." % self.msg
        )


class HSP(Alignment):
    """Stores an alignment of one query sequence against a target sequence.

    An HSP (High-scoring Segment Pair) stores the alignment of one query
    sequence segment against one target (hit) sequence segment. The
    ``Bio.Blast.HSP`` class inherits from the ``Bio.Align.Alignment`` class.

    In addition to the ``target`` and ``query`` attributes of a
    ``Bio.Align.Alignment``, a ``Bio.Blast.HSP`` object has the following
    attributes:

     - score:       score of HSP;
     - annotations: a dictionary that may contain the following keys:
                     - 'bit score': score (in bits) of HSP (float);
                     - 'evalue':    e-value of HSP (float);
                     - 'identity':  number of identities in HSP (integer);
                     - 'positive':  number of positives in HSP (integer);
                     - 'gaps':      number of gaps in HSP (integer);
                     - 'midline':   formatting middle line.

    A ``Bio.Blast.HSP`` object behaves the same as a `Bio.Align.Alignment``
    object and can be used as such. However, when printing a ``Bio.Blast.HSP``
    object, the BLAST e-value and bit score are included in the output (in
    addition to the alignment itself).

    See the documentation of ``Bio.Blast.Record`` for a more detailed
    explanation of how the information in BLAST records is stored in
    Biopython.
    """

    def __repr__(self):
        query = self.query
        target = self.target
        n, m = self.shape
        return f"<Bio.Blast.HSP target.id={target.id!r} query.id={query.id!r}; {n} rows x {m} columns>"

    def __str__(self):
        alignment_text = super().__str__()
        query = self.query
        target = self.target
        query_strand = (
            "Plus" if self.coordinates[1, 0] <= self.coordinates[1, -1] else "Minus"
        )
        target_strand = (
            "Plus" if self.coordinates[0, 0] <= self.coordinates[0, -1] else "Minus"
        )
        indent = " " * 8
        query_description = textwrap.fill(
            query.description,
            width=80,
            initial_indent=indent,
            subsequent_indent=indent,
        )
        target_description = textwrap.fill(
            target.description,
            width=80,
            initial_indent=indent,
            subsequent_indent=indent,
        )
        evalue = self.annotations["evalue"]
        bitscore = self.annotations["bit score"]
        score = self.score
        steps = np.diff(self.coordinates, 1)
        aln_span = sum(abs(steps).max(0))
        terms = []
        identity = self.annotations["identity"]
        identity_percentage = round(100.0 * identity / aln_span)
        identity_text = "Identities:%d/%d(%d%%)" % (
            identity,
            aln_span,
            identity_percentage,
        )
        terms.append(identity_text)
        positive = self.annotations.get("positive")
        if positive is not None:
            positive_percentage = round(100.0 * positive / aln_span)
            positive_text = "Positives:%d/%d(%d%%)" % (
                positive,
                aln_span,
                positive_percentage,
            )
            terms.append(positive_text)
        try:
            gaps = self.annotations["gaps"]
        except KeyError:
            pass
        else:
            gaps_percentage = round(100.0 * gaps / aln_span)
            gaps_text = "Gaps:%d.%d(%d%%)" % (gaps, aln_span, gaps_percentage)
            terms.append(gaps_text)
        counts_line = ",  ".join(terms)
        return """\
Query : %s Length: %d Strand: %s
%s
Target: %s Length: %d Strand: %s
%s

Score:%d bits(%d), Expect:%.1g,
%s

%s
""" % (
            query.id,
            len(query),
            query_strand,
            query_description,
            target.id,
            len(target),
            target_strand,
            target_description,
            bitscore,
            score,
            evalue,
            counts_line,
            alignment_text,
        )


class Hit(Alignments):
    """Stores a single BLAST hit of one single query against one target.

    The ``Bio.Blast.Hit`` class inherits from the ``Bio.Align.Alignments``
    class, which is a subclass of a Python list. The ``Bio.Blast.Hit`` class
    stores ``Bio.Blast.HSP`` objwcts, which inherit from
    ``Bio.Align.Alignment``. A ``Bio.Blast.Hit`` object is therefore
    effectively a list of ``Bio.Align.Alignment`` objects. Most hits consist of
    only 1 or a few Alignment objects.

    Each ``Bio.Blast.Hit`` object has a ``target`` attribute containing the
    following information:

     - target.id:          seqId of subject;
     - target.description: definition line of subject;
     - target.name:        accession of subject;
     - len(target.seq):    sequence length of subject.

    See the documentation of ``Bio.Blast.Record`` for a more detailed
    explanation of the information stored in the alignments contained in the
    ``Bio.Blast.Hit`` object.
    """

    def __getitem__(self, key):
        try:
            value = super().__getitem__(key)
        except IndexError:
            raise IndexError("index out of range") from None
        if isinstance(key, slice):
            hit = Hit(value)
            hit.target = self.target
            return hit
        else:
            return value

    def __repr__(self):
        target = self.target
        try:
            alignment = self[0]
        except IndexError:
            return f"<Bio.Blast.Hit target.id={target.id!r}; no hits>"
        query = alignment.query
        nhsps = len(self)
        if nhsps == 1:
            unit = "HSP"
        else:  # nhsps > 1
            unit = "HSPs"
        return f"<Bio.Blast.Hit target.id={target.id!r} query.id={query.id!r}; {nhsps} {unit}>"

    def __str__(self):
        """Return a human readable summary of the Hit object."""
        lines = []

        # set query id line
        query = self[0].query
        qid_line = "Query: %s" % query.id
        lines.append(qid_line)
        indent = " " * 7
        description_lines = textwrap.wrap(
            query.description,
            width=80,
            initial_indent=indent,
            subsequent_indent=indent,
        )
        lines.extend(description_lines)

        # set hit id line
        target = self.target
        hid_line = "  Hit: %s (length=%i)" % (target.id, len(target))
        lines.append(hid_line)
        description_lines = textwrap.wrap(
            target.description,
            width=80,
            initial_indent=indent,
            subsequent_indent=indent,
        )
        lines.extend(description_lines)

        # set hsp line and table
        lines.append(
            " HSPs: %s  %s  %s  %s  %s  %s"
            % ("-" * 4, "-" * 8, "-" * 9, "-" * 6, "-" * 15, "-" * 21)
        )
        pattern = "%11s  %8s  %9s  %6s  %15s  %21s"
        lines.append(
            pattern % ("#", "E-value", "Bit score", "Span", "Query range", "Hit range")
        )
        lines.append(pattern % ("-" * 4, "-" * 8, "-" * 9, "-" * 6, "-" * 15, "-" * 21))
        for idx, hsp in enumerate(self):
            # evalue
            evalue = format(hsp.annotations["evalue"], ".2g")
            # bit score
            bitscore = format(hsp.annotations["bit score"], ".2f")
            # alignment length
            steps = np.diff(hsp.coordinates, 1)
            aln_span = sum(abs(steps).max(0))
            # query region
            query_start = hsp.coordinates[1, 0]
            query_end = hsp.coordinates[1, -1]
            query_range = f"[{query_start}:{query_end}]"
            # max column length is 18
            query_range = (
                query_range[:13] + "~]" if len(query_range) > 15 else query_range
            )
            # hit region
            hit_start = hsp.coordinates[0, 0]
            hit_end = hsp.coordinates[0, -1]
            hit_range = f"[{hit_start}:{hit_end}]"
            hit_range = hit_range[:19] + "~]" if len(hit_range) > 21 else hit_range
            # append the hsp row
            lines.append(
                pattern % (idx, evalue, bitscore, aln_span, query_range, hit_range)
            )

        return "\n".join(lines)


class Record(list):
    """Stores the BLAST results for a single query.

    A ``Bio.Blast.Record`` object is a list of ``Bio.Blast.Hit`` objects, each
    corresponding to one hit for the query in the BLAST output.

    The ``Bio.Blast.Record`` object may have the following attributes:
     - query:   A ``SeqRecord`` object which may contain some or all of the
                following information:
                 - query.id:          SeqId of query;
                 - query.description: Definition line of query;
                 - len(query.seq):    Length of the query sequence.
     - stat:    A dictionary with summary statistics of the BLAST run. It may
                contain the following keys:
                 - 'db-num':    number of sequences in BLAST db (integer);
                 - 'db-len':    length of BLAST db (integer);
                 - 'hsp-len':   effective HSP length (integer);
                 - 'eff-space': effective search space (float);
                 - 'kappa':     Karlin-Altschul parameter K (float);
                 - 'lambda':    Karlin-Altschul parameter Lambda (float);
                 - 'entropy':   Karlin-Altschul parameter H (float).
     - message: Some (error?) information.

    Each ``Bio.Blast.Hit`` object has a ``target`` attribute containing the
    following information:

     - target.id:          seqId of subject;
     - target.description: definition line of subject;
     - target.name:        accession of subject;
     - len(target.seq):    sequence length of subject.

    The ``Bio.Blast.Hit`` class inherits from the ``Bio.Align.Alignments``
    class, which inherits from a Python list. In this list, the
    ``Bio.Blast.Hit`` object stores ``Bio.Blast.HSP`` objects, which inherit
    from the ``Bio.Align.Alignment`` class.  A ``Bio.Blast.Hit`` object is
    therefore effectively a list of alignment objects.

    Each HSP in a ``Bio.Blast.Hit`` object has the attributes ``target`` and
    ``query`` attributes, as usual for of a ``Bio.Align.Alignment`` object
    storing a pairwise alignment, pointing to a ``SeqRecord`` object
    representing the target and query, respectively.  For translated BLAST
    searches, the ``features`` attribute of the target or query may contain a
    ``SeqFeature`` of type CDS that stores the amino acid sequence region.  The
    ``qualifiers`` attribute of such a feature is a dictionary with  a single
    key 'coded_by'; the corresponding value specifies the nucleotide sequence
    region, in a GenBank-style string with 1-based coordinates, that encodes
    the amino acid sequence.

    Each ``Bio.Blast.HSP`` object has the following additional attributes:

     - score:       score of HSP;
     - annotations: a dictionary that may contain the following keys:
                     - 'bit score': score (in bits) of HSP (float);
                     - 'evalue':    e-value of HSP (float);
                     - 'identity':  number of identities in HSP (integer);
                     - 'positive':  number of positives in HSP (integer);
                     - 'gaps':      number of gaps in HSP (integer);
                     - 'midline':   formatting middle line.

    >>> from Bio import Blast
    >>> record = Blast.read("Blast/xml_2212L_blastx_001.xml")
    >>> record.query
    SeqRecord(seq=Seq(None, length=556), id='gi|1347369|gb|G25137.1|G25137', name='<unknown name>', description='human STS EST48004, sequence tagged site', dbxrefs=[])
    >>> record.stat
    {'db-num': 2934173, 'db-len': 1011751523, 'hsp-len': 0, 'eff-space': 0, 'kappa': 0.041, 'lambda': 0.267, 'entropy': 0.14}
    >>> len(record)
    78
    >>> hit = record[0]
    >>> type(hit)
    <class 'Bio.Blast.Hit'>
    >>> from Bio.Align import Alignments
    >>> isinstance(hit, Alignments)
    True
    >>> hit.target
    SeqRecord(seq=Seq(None, length=319), id='gi|12654095|gb|AAH00859.1|', name='AAH00859', description='Unknown (protein for IMAGE:3459481) [Homo sapiens]', dbxrefs=[])

    Most hits consist of only 1 or a few Alignment objects:

    >>> len(hit)
    1
    >>> alignment = hit[0]
    >>> type(alignment)
    <class 'Bio.Blast.HSP'>
    >>> alignment.score
    630.0
    >>> alignment.annotations
    {'bit score': 247.284, 'evalue': 1.69599e-64, 'identity': 122, 'positive': 123, 'gaps': 0, 'midline': 'DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERF GPYTDFTP TTE QKL EQAL TYPVNT ERW  IA AVPGR K+'}

    Target and query information are stored in the respective attributes of the
    alignment:

    >>> alignment.target
    SeqRecord(seq=Seq({155: 'DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDP...TKK'}, length=319), id='gi|12654095|gb|AAH00859.1|', name='AAH00859', description='Unknown (protein for IMAGE:3459481) [Homo sapiens]', dbxrefs=[])
    >>> alignment.query
    SeqRecord(seq=Seq('DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDP...XKE'), id='gi|1347369|gb|G25137.1|G25137', name='<unknown name>', description='human STS EST48004, sequence tagged site', dbxrefs=[])

    This was a BLASTX run, so the query sequence was translated:

    >>> len(alignment.target.features)
    0
    >>> len(alignment.query.features)
    1
    >>> feature = alignment.query.features[0]
    >>> feature
    SeqFeature(SimpleLocation(ExactPosition(0), ExactPosition(133)), type='CDS', qualifiers=...)
    >>> feature.qualifiers
    {'coded_by': 'gi|1347369|gb|G25137.1|G25137:1..399'}

    i.e., nucleotides 0:399 (in zero-based coordinates) encode the amino acids
    of the query in the alignment.

    For an alignment against the reverse strand, the location in the qualifier
    is shown as in this example:

    >>> record[72][0].query.features[0].qualifiers
    {'coded_by': 'complement(gi|1347369|gb|G25137.1|G25137:345..530)'}

    """

    def __init__(self):
        """Initialize the Record object."""
        self.query = None

    def __repr__(self):
        query = self.query
        try:
            query_id = query.id
        except AttributeError:
            query_id = "unknown"
        nhits = len(self)
        if nhits == 0:
            return f"<Bio.Blast.Record query.id={query_id!r}; no hits>"
        elif nhits == 1:
            return f"<Bio.Blast.Record query.id={query_id!r}; 1 hit>"
        else:
            return f"<Bio.Blast.Record query.id={query_id!r}; {nhits} hits>"

    def __str__(self):
        lines = []
        try:
            version = self.version
        except AttributeError:
            pass
        else:
            lines.append(f"Program: {version}")
        try:
            db = self.db
        except AttributeError:
            pass
        else:
            lines.append(f"     db: {db}")
        if self.query is not None:
            # self.query may be None with legacy Blast if there are no hits
            lines.append("  Query: %s (length=%d)" % (self.query.id, len(self.query)))
            indent = " " * 9
            description_lines = textwrap.wrap(
                self.query.description,
                width=80,
                initial_indent=indent,
                subsequent_indent=indent,
            )
            lines.extend(description_lines)
        if len(self) == 0:
            lines.append("   Hits: No hits found")
        else:
            lines.append("   Hits: %s  %s  %s" % ("-" * 4, "-" * 5, "-" * 58))
            pattern = "%13s  %5s  %s"
            lines.append(pattern % ("#", "# HSP", "ID + description"))
            lines.append(pattern % ("-" * 4, "-" * 5, "-" * 58))
            for idx, hit in enumerate(self):
                n = len(hit)  # Number of HSPs
                if idx < 30:
                    hid_line = "%s  %s" % (hit.target.id, hit.target.description)
                    if len(hid_line) > 58:
                        hid_line = hid_line[:55] + "..."
                    lines.append(pattern % (idx, len(hit), hid_line))
                elif idx > len(self) - 4:
                    hid_line = "%s  %s" % (hit.target.id, hit.target.description)
                    if len(hid_line) > 58:
                        hid_line = hid_line[:55] + "..."
                    lines.append(pattern % (idx, len(hit), hid_line))
                elif idx == 30:
                    lines.append("%14s" % "~~~")
        return "\n".join(lines)

    def __getitem__(self, key):
        try:
            value = super().__getitem__(key)
        except IndexError:
            raise IndexError("index out of range") from None
        except TypeError:
            if not isinstance(key, str):
                raise TypeError("key must be an integer, slice, or str") from None
            for hit in self:
                if hit.target.id == key:
                    return hit
            raise KeyError(key)
        else:
            if isinstance(key, slice):
                record = Record()
                record.extend(value)
                # Only store the query attribute, as the other attributes
                # pertain to the complete Blast record:
                try:
                    query = self.query
                except AttributeError:
                    pass
                else:
                    record.query = query
                # The following keys may be present if the record was created
                # by Blast.read:
                keys = (
                    "source",
                    "program",
                    "version",
                    "reference",
                    "db",
                    "param",
                    "mbstat",
                )
                for key in keys:
                    try:
                        value = getattr(self, key)
                    except AttributeError:
                        pass
                    else:
                        setattr(record, key, value)
                return record
            return value

    def keys(self):
        """Return a list of the target.id of each hit."""
        return [hit.target.id for hit in self]

    def __contains__(self, key):
        for hit in self:
            if hit.target.id == key:
                return True
        return False

    def index(self, key):
        """Return the index of the hit for which the target.id is equal to the key."""
        for i, hit in enumerate(self):
            if hit.target.id == key:
                return i
        raise ValueError(f"'{key}' not found")


class Records(UserList):
    """Stores the BLAST results of a single BLAST run.

    A ``Bio.Blast.Records`` object is an iterator. Iterating over it returns
    returns ``Bio.Blast.Record`` objects, each of which corresponds to one
    BLAST query.

    Common attributes of a ``Bio.Blast.Records`` object are

     - source:     The input data from which the ``Bio.Blast.Records`` object
                   was constructed.
     - program:    The specific BLAST program that was used (e.g., 'blastn').
     - version:    The version of the BLAST program (e.g., 'BLASTN 2.2.27+').
     - reference:  The literature reference to the BLAST publication.
     - db:         The BLAST database against which the query was run
                   (e.g., 'nr').
     - query:      A ``SeqRecord`` object which may contain some or all of the
                   following information:
                    - query.id:          SeqId of the query;
                    - query.description: Definition line of the query;
                    - query.seq:         The query sequence. The query sequence.
                                         The query sequence.
     - param:      A dictionary with the parameters used for the BLAST run.
                   You may find the following keys in this dictionary:
                    - 'matrix':       the scoring matrix used in the BLAST run
                                      (e.g., 'BLOSUM62') (string);
                    - 'expect':       threshold on the expected number of chance
                                      matches (float);
                    - 'include':      e-value threshold for inclusion in
                                      multipass model in psiblast (float);
                    - 'sc-match':     score for matching nucleotides (integer);
                    - 'sc-mismatch':  score for mismatched nucleotides
                                      (integer);
                    - 'gap-open':     gap opening cost (integer);
                    - 'gap-extend':   gap extension cost (integer);
                    - 'filter':       filtering options applied in the BLAST
                                      run (string);
                    - 'pattern':      PHI-BLAST pattern (string);
                    - 'entrez-query': Limit of request to Entrez query (string).
     - mbstat:     A dictionary with Mega BLAST search statistics.  As this
                   information is stored near the end of the XML file, this
                   attribute can only be accessed after the file has been read
                   completely (by iterating over the records until a
                   ``StopIteration`` is issued. This dictionary can contain the
                   same keys as the dictionary stored under the ``stat``
                   attribute of a ``Record`` object.

    >>> from Bio import Blast
    >>> path = "Blast/xml_2218_blastp_002.xml"

    In a script, you would use a ``with`` block, as in

    >>> with Blast.parse(path) as records:
    ...     print(records.source)
    ...
    Blast/xml_2218_blastp_002.xml

    to ensure that the file is closed at the end of the block.
    Here, we will simply do

    >>> records = Blast.parse("Blast/xml_2218_blastp_002.xml")

    so we can see the output of each command right away.

    >>> type(records)
    <class 'Bio.Blast.Records'>
    >>> records.source
    'Blast/xml_2218_blastp_002.xml'
    >>> records.program
    'blastp'
    >>> records.version
    'BLASTP 2.2.18+'
    >>> records.reference
    'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.'
    >>> records.db
    'gpipe/9606/Previous/protein'
    >>> records.param
    {'matrix': 'BLOSUM62', 'expect': 0.01, 'gap-open': 11, 'gap-extend': 1, 'filter': 'm L; R -d repeat/repeat_9606;'}

    Iterating over the records returns Bio.Blast.Record objects:

    >>> record = next(records)
    >>> type(record)
    <class 'Bio.Blast.Record'>
    >>> record.query.id
    'gi|585505|sp|Q08386|MOPB_RHOCA'
    >>> record = next(records)
    >>> type(record)
    <class 'Bio.Blast.Record'>
    >>> record.query.id
    'gi|129628|sp|P07175.1|PARA_AGRTU'
    >>> record = next(records)  # doctest:+ELLIPSIS
    Traceback (most recent call last):
    ...
    StopIteration

    You can also use the records as a list, for example by extracting a record
    by index, or by calling ``len`` or ``print`` on the records. The parser
    will then automatically iterate over the records and store them:

    >>> records = Blast.parse("Blast/wnts.xml")
    >>> record = records[3]  # this causes all records to be read in and stored
    >>> record.query.id
    'Query_4'
    >>> len(records)
    5

    After the records have been read in, you can still iterate over them:

    >>> for i, record in enumerate(records):
    ...     print(i, record.query.id)
    ...
    0 Query_1
    1 Query_2
    2 Query_3
    3 Query_4
    4 Query_5

    """  # noqa: RST201, RST203, RST301

    def __init__(self, source):
        """Initialize the Records object."""
        if isinstance(source, list):  # UserList API requirement
            self._records = source
            self._loaded = True
            self._index = 0
            return

        self.source = source
        try:
            stream = open(source, "rb")
        except TypeError:  # not a path, assume we received a stream
            if source.read(0) != b"":
                raise StreamModeError(
                    "BLAST output files must be opened in binary mode."
                ) from None
            stream = source
        self._stream = stream
        self._read_header()
        self._records = []  # for when we want to use Records as a list
        self._loaded = False

    def _read_header(self):
        from Bio.Blast._parser import XMLHandler

        stream = self._stream
        source = self.source
        try:  # context manager won't kick in until after parse returns
            parser = expat.ParserCreate(namespace_separator=" ")
            self._parser = parser

            handler = XMLHandler(parser)
            handler._records = self

            while True:
                data = stream.read(BLOCK)
                if data == b"":
                    try:
                        handler._parser
                    except AttributeError:
                        break
                    else:
                        raise ValueError(
                            f"premature end of XML file (after reading {parser.CurrentByteIndex} bytes)"
                        )
                try:
                    parser.Parse(data, False)
                except expat.ExpatError as e:
                    if parser.StartElementHandler:
                        # We saw the initial <!xml declaration, so we can be
                        # sure that we are parsing XML data. Most likely, the
                        # XML file is corrupted.
                        raise CorruptedXMLError(e) from None
                    else:
                        # We have not seen the initial <!xml declaration, so
                        # probably the input data is not in XML format.
                        raise NotXMLError(e) from None
                try:
                    self._cache
                except AttributeError:
                    pass
                else:
                    # We have finished reading the header
                    break
        except Exception:
            if stream is not source:
                stream.close()
            raise
        self._index = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        try:
            stream = self._stream
        except AttributeError:
            return
        if stream is not self.source:
            stream.close()
        del self._stream

    def __iter__(self):
        return self

    def __next__(self):
        if self._loaded is True:
            try:
                record = self._records[self._index]
            except IndexError:
                raise StopIteration from None
            self._index += 1
            return record
        try:
            cache = self._cache
        except AttributeError:
            raise StopIteration from None
        parser = self._parser
        stream = self._stream
        while True:
            try:
                record = self._cache.popleft()
            except IndexError:  # no record ready to be returned
                pass
            else:
                self._index += 1
                return record
            # Read in another block of data from the file.
            data = stream.read(BLOCK)
            if data == b"":
                del self._cache
                del self._parser
                if parser.StartElementHandler is not None:
                    raise ValueError(
                        f"premature end of XML file (after reading {parser.CurrentByteIndex} bytes)"
                    )
                raise StopIteration
            try:
                parser.Parse(data, False)
            except expat.ExpatError as e:
                raise CorruptedXMLError(e) from None

    def __getitem__(self, index):
        item = super().__getitem__(index)
        if index == slice(None, None, None):
            for key, value in self.__dict__.items():
                if key not in ("_stream", "_parser", "_cache"):
                    item.__dict__[key] = value
        return item

    @property
    def data(self):
        """Overrides the data attribute of UserList."""
        if self._loaded is False:
            # Read all records and store them
            index = self._index
            if index > 0:
                try:
                    self._stream.seek(0)
                except io.UnsupportedOperation:
                    raise ValueError(
                        "list-like access after iterating is supported only if the input data is seekable."
                    )
                self._read_header()
                self._index = 0
            for record in self:
                self._records.append(record)
            stream = self._stream
            if stream is not self.source:
                stream.close()
            del self._stream
            self._loaded = True
            self._index = index
        return self._records

    def __repr__(self):
        return f"<Bio.Blast.Records source={self.source!r} program={self.program!r} version={self.version!r} db={self.db!r}>"

    def __str__(self):
        text = """\
Program: %s
     db: %s""" % (
            self.version,
            self.db,
        )
        records = self[:]  # to ensure that the records are read in
        for record in self._records:
            text += "\n\n" + str(record)
        return text


def parse(source):
    """Parse an XML file containing BLAST output and return a Bio.Blast.Records object.

    This returns an iterator object; iterating over it returns Bio.Blast.Record
    objects one by one.

    The source can be a file stream or the path to an XML file containing the
    BLAST output. If a file stream, source  must be in binary mode. This allows
    the parser to detect the encoding from the XML file,and to use it to convert
    any text in the XML to the correct Unicode string. The qblast function in
    Bio.Blast returns a file stream in binary mode. For files, please use mode
    "rb" when opening the file, as in

    >>> from Bio import Blast
    >>> stream = open("Blast/wnts.xml", "rb")  # opened in binary mode
    >>> records = Blast.parse(stream)
    >>> for record in records:
    ...     print(record.query.id, record.query.description)
    ...
    Query_1 gi|195230749:301-1383 Homo sapiens wingless-type MMTV integration site family member 2 (WNT2), transcript variant 1, mRNA
    Query_2 gi|325053704:108-1166 Homo sapiens wingless-type MMTV integration site family, member 3A (WNT3A), mRNA
    Query_3 gi|156630997:105-1160 Homo sapiens wingless-type MMTV integration site family, member 4 (WNT4), mRNA
    Query_4 gi|371502086:108-1205 Homo sapiens wingless-type MMTV integration site family, member 5A (WNT5A), transcript variant 2, mRNA
    Query_5 gi|53729353:216-1313 Homo sapiens wingless-type MMTV integration site family, member 6 (WNT6), mRNA
    >>> stream.close()

    """
    return Records(source)


def read(source):
    """Parse an XML file containing BLAST output for a single query and return it.

    Internally, this function uses Bio.Blast.parse to obtain an iterator over
    BLAST records.  The function then reads one record from the iterator,
    ensures that there are no further records, and returns the record it found
    as a Bio.Blast.Record object. An exception is raised if no records are
    found, or more than one record is found.

    The source can be a file stream or the path to an XML file containing the
    BLAST output. If a file stream, source  must be in binary mode. This allows
    the parser to detect the encoding from the XML file,and to use it to convert
    any text in the XML to the correct Unicode string. The qblast function in
    Bio.Blast returns a file stream in binary mode. For files, please use mode
    "rb" when opening the file, as in

    >>> from Bio import Blast
    >>> stream = open("Blast/xml_21500_blastn_001.xml", "rb")  # opened in binary mode
    >>> record = Blast.read(stream)
    >>> record.query.id
    'Query_78041'
    >>> record.query.description
    'G26684.1 human STS STS_D11570, sequence tagged site'
    >>> len(record)
    11
    >>> stream.close()

    Use the Bio.Blast.parse function if you want to read a file containing
    BLAST output for more than one query.
    """
    with parse(source) as records:
        try:
            record = next(records)
        except StopIteration:
            raise ValueError("No BLAST output found.") from None
        try:
            next(records)
            raise ValueError("BLAST output for more than one query found.")
        except StopIteration:
            pass
    for key in ("source", "program", "version", "reference", "db", "param", "mbstat"):
        try:
            value = getattr(records, key)
        except AttributeError:
            pass
        else:
            setattr(record, key, value)
    return record


def write(records, destination, fmt="XML"):
    """Write BLAST records as an XML file, and return the number of records.

    Arguments:
     - records     - A ``Bio.Blast.Records`` object.
     - destination - File or file-like object to write to, or filename as
                     string.
                     The File object must have been opened for writing in
                     binary mode, and must be closed (or flushed) by the caller
                     after this function returns to ensure that all records are
                     written.
     - fmt         - string describing the file format to write
                     (case-insensitive).
                     Currently, only "XML" and "XML2" are accepted.

    Returns the number of records written (as an integer).
    """
    fmt = fmt.upper()
    if fmt == "XML":
        Writer = _writers.XMLWriter
    elif fmt == "XML2":
        Writer = _writers.XML2Writer
    else:
        raise ValueError(f"Unknown format {fmt}; expected 'XML' or 'XML2'")
    try:
        stream = open(destination, "wb")
    except TypeError:  # not a path, assume we received a stream
        try:
            destination.write(b"")
        except TypeError:
            # destination was opened in text mode
            raise StreamModeError(
                "File must be opened in binary mode for writing."
            ) from None
        stream = destination
    writer = Writer(stream)
    try:
        count = writer.write(records)
    finally:
        if stream is not destination:
            stream.close()
    return count


@function_with_previous
def qblast(
    program,
    database,
    sequence,
    url_base=NCBI_BLAST_URL,
    auto_format=None,
    composition_based_statistics=None,
    db_genetic_code=None,
    endpoints=None,
    entrez_query="(none)",
    expect=10.0,
    filter=None,
    gapcosts=None,
    genetic_code=None,
    hitlist_size=50,
    i_thresh=None,
    layout=None,
    lcase_mask=None,
    matrix_name=None,
    nucl_penalty=None,
    nucl_reward=None,
    other_advanced=None,
    perc_ident=None,
    phi_pattern=None,
    query_file=None,
    query_believe_defline=None,
    query_from=None,
    query_to=None,
    searchsp_eff=None,
    service=None,
    threshold=None,
    ungapped_alignment=None,
    word_size=None,
    short_query=None,
    alignments=500,
    alignment_view=None,
    descriptions=500,
    entrez_links_new_window=None,
    expect_low=None,
    expect_high=None,
    format_entrez_query=None,
    format_object=None,
    format_type="XML",
    ncbi_gi=None,
    results_file=None,
    show_overview=None,
    megablast=None,
    template_type=None,
    template_length=None,
    username="blast",
    password=None,
):
    """BLAST search using NCBI's QBLAST server or a cloud service provider.

    Supports all parameters of the old qblast API for Put and Get.

    Please note that NCBI uses the new Common URL API for BLAST searches
    on the internet (http://ncbi.github.io/blast-cloud/dev/api.html). Thus,
    some of the parameters used by this function are not (or are no longer)
    officially supported by NCBI. Although they are still functioning, this
    may change in the future.

    The Common URL API (http://ncbi.github.io/blast-cloud/dev/api.html) allows
    doing BLAST searches on cloud servers. To use this feature, please set
    ``url_base='http://host.my.cloud.service.provider.com/cgi-bin/blast.cgi'``
    and ``format_object='Alignment'``. For more details, please see
    https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=CloudBlast

    Some useful parameters:

     - program        blastn, blastp, blastx, tblastn, or tblastx (lower case)
     - database       Which database to search against (e.g. "nr").
     - sequence       The sequence to search.
     - ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
     - descriptions   Number of descriptions to show.  Def 500.
     - alignments     Number of alignments to show.  Def 500.
     - expect         An expect value cutoff.  Def 10.0.
     - matrix_name    Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
     - filter         "none" turns off filtering.  Default no filtering
     - format_type    "XML" (default), "HTML", "Text", "XML2", "JSON2",
                      or "Tabular".
     - entrez_query   Entrez query to limit Blast search
     - hitlist_size   Number of hits to return. Default 50
     - megablast      TRUE/FALSE whether to use MEga BLAST algorithm (blastn only)
     - short_query    TRUE/FALSE whether to adjust the search parameters for a
                      short query sequence. Note that this will override
                      manually set parameters like word size and e value. Turns
                      off when sequence length is > 30 residues. Default: None.
     - service        plain, psi, phi, rpsblast, megablast (lower case)

    This function does no checking of the validity of the parameters
    and passes the values to the server as is.  More help is available at:
    https://ncbi.github.io/blast-cloud/dev/api.html

    The http.client.HTTPResponse object returned by this function has the
    additional attributes rid and rtoe with the Request ID and Request Time Of
    Execution for this BLAST search.

    """
    programs = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]
    if program not in programs:
        raise ValueError(
            f"Program specified is {program}. Expected one of {', '.join(programs)}"
        )

    # SHORT_QUERY_ADJUST throws an error when using blastn (wrong parameter
    # assignment from NCBIs side).
    # Thus we set the (known) parameters directly:
    if short_query and program == "blastn":
        short_query = None
        # We only use the 'short-query' parameters for short sequences:
        if len(sequence) < 31:
            expect = 1000
            word_size = 7
            nucl_reward = 1
            filter = None
            lcase_mask = None
            warnings.warn(
                '"SHORT_QUERY_ADJUST" is incorrectly implemented (by NCBI) for blastn.'
                " We bypass the problem by manually adjusting the search parameters."
                " Thus, results may slightly differ from web page searches.",
                BiopythonWarning,
            )

    # Format the "Put" command, which sends search requests to qblast.
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node5.html on 9 July 2007
    # Additional parameters are taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node9.html on 8 Oct 2010
    # To perform a PSI-BLAST or PHI-BLAST search the service ("Put" and "Get" commands) must be specified
    # (e.g. psi_blast = NCBIWWW.qblast("blastp", "refseq_protein", input_sequence, service="psi"))
    parameters = {
        "AUTO_FORMAT": auto_format,
        "COMPOSITION_BASED_STATISTICS": composition_based_statistics,
        "DATABASE": database,
        "DB_GENETIC_CODE": db_genetic_code,
        "ENDPOINTS": endpoints,
        "ENTREZ_QUERY": entrez_query,
        "EXPECT": expect,
        "FILTER": filter,
        "GAPCOSTS": gapcosts,
        "GENETIC_CODE": genetic_code,
        "HITLIST_SIZE": hitlist_size,
        "I_THRESH": i_thresh,
        "LAYOUT": layout,
        "LCASE_MASK": lcase_mask,
        "MEGABLAST": megablast,
        "MATRIX_NAME": matrix_name,
        "NUCL_PENALTY": nucl_penalty,
        "NUCL_REWARD": nucl_reward,
        "OTHER_ADVANCED": other_advanced,
        "PERC_IDENT": perc_ident,
        "PHI_PATTERN": phi_pattern,
        "PROGRAM": program,
        # ('PSSM': pssm: - It is possible to use PSI-BLAST via this API?
        "QUERY": sequence,
        "QUERY_FILE": query_file,
        "QUERY_BELIEVE_DEFLINE": query_believe_defline,
        "QUERY_FROM": query_from,
        "QUERY_TO": query_to,
        # 'RESULTS_FILE': ...: - Can we use this parameter?
        "SEARCHSP_EFF": searchsp_eff,
        "SERVICE": service,
        "SHORT_QUERY_ADJUST": short_query,
        "TEMPLATE_TYPE": template_type,
        "TEMPLATE_LENGTH": template_length,
        "THRESHOLD": threshold,
        "UNGAPPED_ALIGNMENT": ungapped_alignment,
        "WORD_SIZE": word_size,
        "CMD": "Put",
    }

    if password is not None:
        # handle authentication for BLAST cloud
        password_mgr = HTTPPasswordMgrWithDefaultRealm()
        password_mgr.add_password(None, url_base, username, password)
        handler = HTTPBasicAuthHandler(password_mgr)
        opener = build_opener(handler)
        install_opener(opener)

    if url_base == NCBI_BLAST_URL:
        parameters.update({"email": email, "tool": tool})
    parameters = {key: value for key, value in parameters.items() if value is not None}
    message = urlencode(parameters).encode()
    request = Request(url_base, message, {"User-Agent": "BiopythonClient"})
    # Send off the initial query to qblast.
    # Note the NCBI do not currently impose a rate limit here, other
    # than the request not to make say 50 queries at once using multiple
    # threads.
    stream = urlopen(request)

    # Format the "Get" command, which gets the formatted results from qblast
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007
    rid, rtoe = _parse_qblast_ref_page(stream)
    parameters = {
        "ALIGNMENTS": alignments,
        "ALIGNMENT_VIEW": alignment_view,
        "DESCRIPTIONS": descriptions,
        "ENTREZ_LINKS_NEW_WINDOW": entrez_links_new_window,
        "EXPECT_LOW": expect_low,
        "EXPECT_HIGH": expect_high,
        "FORMAT_ENTREZ_QUERY": format_entrez_query,
        "FORMAT_OBJECT": format_object,
        "FORMAT_TYPE": format_type,
        "NCBI_GI": ncbi_gi,
        "RID": rid,
        "RESULTS_FILE": results_file,
        "SERVICE": service,
        "SHOW_OVERVIEW": show_overview,
        "CMD": "Get",
    }
    parameters = {key: value for key, value in parameters.items() if value is not None}
    message = urlencode(parameters).encode()

    # Poll NCBI until the results are ready.
    # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
    # 1. Do not contact the server more often than once every 10 seconds.
    # 2. Do not poll for any single RID more often than once a minute.
    # 3. Use the URL parameter email and tool, so that the NCBI
    #    can contact you if there is a problem.
    # 4. Run scripts weekends or between 9 pm and 5 am Eastern time
    #    on weekdays if more than 50 searches will be submitted.
    # --
    # Could start with a 10s delay, but expect most short queries
    # will take longer thus at least 70s with delay. Therefore,
    # start with 20s delay, thereafter once a minute.
    delay = 20  # seconds
    start_time = time.time()
    while True:
        current = time.time()
        wait = qblast.previous + delay - current
        if wait > 0:
            time.sleep(wait)
            qblast.previous = current + wait
        else:
            qblast.previous = current
        # delay by at least 60 seconds only if running the request against the public NCBI API
        if delay < 60 and url_base == NCBI_BLAST_URL:
            # Wasn't a quick return, must wait at least a minute
            delay = 60
        elapsed = time.time() - start_time
        if elapsed >= 600:
            warnings.warn(
                f"BLAST request {rid} is taking longer than 10 minutes, consider re-issuing it",
                BiopythonWarning,
            )
        request = Request(url_base, message, {"User-Agent": "BiopythonClient"})
        stream = urlopen(request)
        data = stream.peek()
        if format_type == "HTML" and b"<title>NCBI Blast:</title>" in data:
            continue
        elif data.startswith(b"<!DOCTYPE html"):
            continue
        else:
            break
    if format_type == "XML":
        assert data.startswith(b"<?xml ")
    elif format_type == "HTML":
        assert data.startswith(b"<!DOCTYPE html ")
    elif format_type in ("Text", "Tabular"):
        assert data.startswith(b"<p><!--\nQBlastInfoBegin")
    elif format_type in ("XML2", "JSON2"):
        assert data.startswith(b"PK\x03\x04")  # zipped file
    stream.rid = rid
    stream.rtoe = rtoe
    return stream


qblast.previous = 0


def _parse_qblast_ref_page(handle):
    """Extract a tuple of RID, RTOE from the 'please wait' page (PRIVATE).

    The NCBI FAQ pages use TOE for 'Time of Execution', so RTOE is probably
    'Request Time of Execution' and RID would be 'Request Identifier'.
    """
    s = handle.read().decode()
    i = s.find("RID =")
    if i == -1:
        rid = None
    else:
        j = s.find("\n", i)
        rid = s[i + len("RID =") : j].strip()

    i = s.find("RTOE =")
    if i == -1:
        rtoe = None
    else:
        j = s.find("\n", i)
        rtoe = s[i + len("RTOE =") : j].strip()

    if not rid and not rtoe:
        # Can we reliably extract the error message from the HTML page?
        # e.g.  "Message ID#24 Error: Failed to read the Blast query:
        #       Nucleotide FASTA provided for protein sequence"
        # or    "Message ID#32 Error: Query contains no data: Query
        #       contains no sequence data"
        #
        # This used to occur inside a <div class="error msInf"> entry:
        i = s.find('<div class="error msInf">')
        if i != -1:
            msg = s[i + len('<div class="error msInf">') :].strip()
            msg = msg.split("</div>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError(f"Error message from NCBI: {msg}")
        # In spring 2010 the markup was like this:
        i = s.find('<p class="error">')
        if i != -1:
            msg = s[i + len('<p class="error">') :].strip()
            msg = msg.split("</p>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError(f"Error message from NCBI: {msg}")
        # Generic search based on the way the error messages start:
        i = s.find("Message ID#")
        if i != -1:
            # Break the message at the first HTML tag
            msg = s[i:].split("<", 1)[0].split("\n", 1)[0].strip()
            raise ValueError(f"Error message from NCBI: {msg}")
        # We didn't recognise the error layout :(
        # print(s)
        raise ValueError(
            "No RID and no RTOE found in the 'please wait' page, "
            "there was probably an error in your request but we "
            "could not extract a helpful error message."
        )
    elif not rid:
        # Can this happen?
        raise ValueError(
            f"No RID found in the 'please wait' page. (although RTOE = {rtoe!r})"
        )
    elif not rtoe:
        # Can this happen?
        raise ValueError(
            f"No RTOE found in the 'please wait' page. (although RID = {rid!r})"
        )

    try:
        return rid, int(rtoe)
    except ValueError:
        raise ValueError(
            f"A non-integer RTOE found in the 'please wait' page, {rtoe!r}"
        ) from None


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
