# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# Copyright 2000 by Bertrand Frottier.  All rights reserved.
# Revisions 2005-2006 copyright Michiel de Hoon
# Revisions 2006-2009 copyright Peter Cock
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code to work with XML output from BLAST.

The BLAST XML DTD file is on the NCBI FTP site at:
ftp://ftp.ncbi.nlm.nih.gov/blast/documents/xml/NCBI_BlastOutput.dtd
"""

# fmt: off
# ruff: noqa
# flake8: noqa
# mypy: ignore-errors

import sys
from xml.parsers import expat
from collections import deque


from Bio import StreamModeError
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment, Alignments


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


class Record:
    pass


class Records:

    BLOCK = 2048  # default block size from expat

    _start_methods = {}
    _end_methods = {}

    def __init__(self, stream):
        parser = expat.ParserCreate()
        parser.XmlDeclHandler = self._xmlDeclHandler
        parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self._parser = parser
        self._stream = stream
        BLOCK = self.BLOCK
        while True:
            data = stream.read(BLOCK)
            try:
                parser.Parse(data, False)
            except expat.ExpatError as e:
                if parser.StartElementHandler:
                    # We saw the initial <!xml declaration, so we can be sure
                    # that we are parsing XML data. Most likely, the XML file
                    # is corrupted.
                    raise CorruptedXMLError(e) from None
                else:
                    # We have not seen the initial <!xml declaration, so
                    # probably the input data is not in XML format.
                    raise NotXMLError(e) from None
            try:
                self._cache
            except AttributeError:
                continue
            else:
                break

    def _start_blastoutput(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_program(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_version(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_reference(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_db(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_query_id(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_query_def(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_param(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        self.param = {}

    def _start_parameters_matrix(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters_expect(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters_sc_match(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters_sc_mismatch(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters_gap_open(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters_gap_extend(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters_filter(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_iterations(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        self._cache = deque()

    def _start_blastoutput_query_len(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_iteration(self, name, attrs):
        record = Record()
        self._record = record

    def _start_iteration_iter_num(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_iteration_query_id(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_iteration_query_def(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_iteration_query_len(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_iteration_hits(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        self._hits = []

    def _start_hit(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        self._hit = Alignments()

    def _start_hit_num(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_id(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_def(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_hsps(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_len(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_accession(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        self._hsp = {}

    def _start_hsp_num(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_bit_score(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_score(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_evalue(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_query_from(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_query_to(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_hit_from(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_hit_to(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_query_frame(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_hit_frame(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_identity(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_positive(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_gaps(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_align_len(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_qseq(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_hseq(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_midline(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_iteration_stat(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        self._record.stat = {}

    def _start_statistics_db_num(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_db_len(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_hsp_len(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_eff_space(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_kappa(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_lambda(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_entropy(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_blastoutput(self, name):
        assert self._characters.strip() == ""
        del self._characters
        del self._parser

    def _end_blastoutput_program(self, name):
        self.program = self._characters
        self._characters = ""

    def _end_blastoutput_version(self, name):
        self.version = self._characters
        self._characters = ""

    def _end_blastoutput_reference(self, name):
        self.reference = self._characters
        self._characters = ""

    def _end_blastoutput_db(self, name):
        self.db = self._characters
        self._characters = ""

    def _end_blastoutput_query_id(self, name):
        query_id = self._characters
        self.query = SeqRecord(None, query_id)
        self._characters = ""

    def _end_blastoutput_query_def(self, name):
        query_def = self._characters
        self.query.description = query_def
        self._characters = ""

    def _end_blastoutput_query_len(self, name):
        length = int(self._characters)
        self.query.seq = Seq(None, length=length)
        self._characters = ""

    def _end_blastoutput_param(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_parameters(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_parameters_matrix(self, name):
        self.param["matrix"] = self._characters
        self._characters = ""

    def _end_parameters_expect(self, name):
        self.param["expect"] = float(self._characters)
        self._characters = ""

    def _end_parameters_sc_match(self, name):
        self.param["sc-match"] = int(self._characters)
        self._characters = ""

    def _end_parameters_sc_mismatch(self, name):
        self.param["sc-mismatch"] = int(self._characters)
        self._characters = ""

    def _end_parameters_gap_open(self, name):
        self.param["gap-open"] = int(self._characters)
        self._characters = ""

    def _end_parameters_gap_extend(self, name):
        self.param["gap-extend"] = int(self._characters)
        self._characters = ""

    def _end_parameters_filter(self, name):
        self.param["filter"] = self._characters
        self._characters = ""

    def _end_blastoutput_iterations(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_iteration(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        self._cache.append(self._record)
        del self._record

    def _end_iteration_iter_num(self, name):
        self._record.num = int(self._characters)
        self._characters = ""

    def _end_iteration_query_id(self, name):
        query_id = self._characters
        self._record.query = SeqRecord(None, query_id)
        self._characters = ""

    def _end_iteration_query_def(self, name):
        query_def = self._characters
        self._record.query.description = query_def
        self._characters = ""

    def _end_iteration_query_len(self, name):
        length = int(self._characters)
        self._record.query.seq = Seq(None, length=length)
        self._characters = ""

    def _end_iteration_hits(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        hits = self._hits
        self._record.hits = hits
        del self._hits

    def _end_hit(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        hit = self._hit
        del self._hit
        self._hits.append(hit)

    def _end_hit_num(self, name):
        num = int(self._characters)
        if num != len(self._hits) + 1:
            raise ValueError(f"unexpected value in tag <Hit_num> (found {num}, expected {len(self._hits) + 1})")
        self._characters = ""

    def _end_hit_id(self, name):
        hit_id = self._characters
        self._hit.target = SeqRecord(None, hit_id)
        self._characters = ""

    def _end_hit_def(self, name):
        description = self._characters
        self._hit.target.description = description
        self._characters = ""

    def _end_hit_accession(self, name):
        accession = self._characters
        self._hit.target.name = accession
        self._characters = ""

    def _end_hit_len(self, name):
        length = int(self._characters)
        self._hit.target.seq = Seq(None, length=length)
        self._characters = ""

    def _end_hit_hsps(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_hsp_num(self, name):
        num = int(self._characters)
        if num != len(self._hit) + 1:
            raise ValueError(f"unexpected value in tag <Hsp_num> (found {num}, expected {len(self._hit) + 1})")
        self._characters = ""

    def _end_hsp_bit_score(self, name):
        self._hsp["bit-score"] = float(self._characters)
        self._characters = ""

    def _end_hsp_score(self, name):
        self._hsp["score"] = float(self._characters)
        self._characters = ""

    def _end_hsp_evalue(self, name):
        self._hsp["evalue"] = float(self._characters)
        self._characters = ""

    def _end_hsp_query_from(self, name):
        self._hsp["query-from"] = int(self._characters)
        self._characters = ""

    def _end_hsp_query_to(self, name):
        self._hsp["query-to"] = int(self._characters)
        self._characters = ""

    def _end_hsp_hit_from(self, name):
        self._hsp["hit-from"] = int(self._characters)
        self._characters = ""

    def _end_hsp_hit_to(self, name):
        self._hsp["hit-to"] = int(self._characters)
        self._characters = ""

    def _end_hsp_query_frame(self, name):
        query_frame = int(self._characters)
        if query_frame not in (-3, -2, -1, 0, +1, +2, +3):
            raise ValueError(f"unexpected value {query_frame} in tag <Hsp_query-frame>")
        self._hsp["query-frame"] = query_frame
        self._characters = ""

    def _end_hsp_hit_frame(self, name):
        hit_frame = int(self._characters)
        if hit_frame not in (-3, -2, -1, 0, +1, +2, +3):
            raise ValueError(f"unexpected value {fit_frame} in tag <Hsp_hit-frame>")
        self._hsp["hit-frame"] = hit_frame
        self._characters = ""

    def _end_hsp_identity(self, name):
        self._hsp["identity"] = int(self._characters)
        self._characters = ""

    def _end_hsp_positive(self, name):
        self._hsp["positive"] = int(self._characters)
        self._characters = ""

    def _end_hsp_gaps(self, name):
        self._hsp["gaps"] = int(self._characters)
        self._characters = ""

    def _end_hsp_align_len(self, name):
        self._hsp["align-len"] = int(self._characters)
        self._characters = ""

    def _end_hsp_qseq(self, name):
        self._hsp["qseq"] = self._characters
        self._characters = ""

    def _end_hsp_hseq(self, name):
        self._hsp["hseq"] = self._characters
        self._characters = ""

    def _end_hsp_midline(self, name):
        self._hsp["midline"] = self._characters
        self._characters = ""

    def _end_hsp(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        hsp = self._hsp
        del self._hsp
        align_len = hsp["align-len"]
        query = self._record.query
        query_id = query.id
        query_description = query.description
        query_length = len(query.seq)
        query_seq_aligned = hsp["qseq"]
        assert len(query_seq_aligned) == align_len
        target_seq_aligned = hsp["hseq"]
        assert len(target_seq_aligned) == align_len
        coordinates = Alignment.infer_coordinates([target_seq_aligned, query_seq_aligned])
        query_seq_data = query_seq_aligned.replace("-", "")
        query_frame =  hsp["query-frame"]
        if self.program == "blastx":
            query_start = hsp["query-from"] - 1
            query_end = hsp["query-to"]
            assert query_end - query_start == 3 * len(query_seq_data)
        else:
            if query_frame == +1 or query_frame == 0:
                query_start = hsp["query-from"] - 1
                query_end = hsp["query-to"]
            elif query_frame == -1:
                query_start = hsp["query-to"] - 1
                query_end = hsp["query-from"]
                query_seq_data = reverse_complement(query_seq_data)
                coordinates[1, :] = coordinates[1, ::-1]
            assert query_end - query_start == len(query_seq_data)
            query_seq_data = {query_start: query_seq_data}
            coordinates[1, :] += query_start
        query_seq = Seq(query_seq_data, query_length)
        query = SeqRecord(query_seq, query_id, description=query_description)
        query.annotations["start"] = query_start
        query.annotations["end"] = query_end
        query.annotations["frame"] = query_frame
        target = self._hit.target
        target_id = target.id
        target_name = target.name
        target_description = target.description
        target_length = len(target.seq)
        target_seq_data = target_seq_aligned.replace("-", "")
        target_frame =  hsp["hit-frame"]
        if target_frame == +1 or target_frame == 0:
            target_start = hsp["hit-from"] - 1
            target_end = hsp["hit-to"]
        elif target_frame == -1:
            target_start = hsp["hit-to"] - 1
            target_end = hsp["hit-from"]
            target_seq_data = reverse_complement(target_seq_data)
        assert target_end - target_start == len(target_seq_data)
        target_seq = Seq({target_start: target_seq_data}, target_length)
        target = SeqRecord(target_seq, target_id, target_name, description=target_description)
        target.annotations["start"] = target_start
        target.annotations["end"] = target_end
        target.annotations["frame"] = target_frame
        coordinates[0, :] += target_start
        if target_frame == -1:
            coordinates[0, :] = coordinates[0, ::-1]
        sequences = [target, query]
        alignment = Alignment(sequences, coordinates)
        alignment.score = hsp["score"]
        annotations = {}
        annotations["bit score"] = hsp["bit-score"]
        annotations["evalue"] = hsp["evalue"]
        annotations["identity"] = hsp["identity"]
        annotations["positive"] = hsp["positive"]
        annotations["gaps"] = hsp["gaps"]
        annotations["midline"] = hsp["midline"]
        alignment.annotations = annotations
        self._hit.append(alignment)

    def _end_iteration_stat(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_statistics(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_statistics_db_num(self, name):
        self._record.stat["db-num"] = int(self._characters)
        self._characters = ""

    def _end_statistics_db_len(self, name):
        self._record.stat["db-len"] = int(self._characters)
        self._characters = ""

    def _end_statistics_hsp_len(self, name):
        self._record.stat["hsp-len"] = int(self._characters)
        self._characters = ""

    def _end_statistics_eff_space(self, name):
        self._record.stat["eff-space"] = float(self._characters)
        self._characters = ""

    def _end_statistics_kappa(self, name):
        self._record.stat["kappa"] = float(self._characters)
        self._characters = ""

    def _end_statistics_lambda(self, name):
        self._record.stat["lambda"] = float(self._characters)
        self._characters = ""

    def _end_statistics_entropy(self, name):
        self._record.stat["entropy"] = float(self._characters)
        self._characters = ""

    def _xmlDeclHandler(self, version, encoding, standalone):
        parser = self._parser
        parser.ExternalEntityRefHandler = self._externalEntityRefHandler
        parser.StartElementHandler = self._startElementHandler
        parser.EndElementHandler = self._endElementHandler
        parser.CharacterDataHandler = self._characterDataHandler
        self._characters = ""

    def _externalEntityRefHandler(self, context, base, systemId, publicId):
        """Handle the DTD declaration.

        """
        assert context is None
        assert base is None
        if systemId != "NCBI_BlastOutput.dtd":
            raise ValueError("output from legacy BLAST program")
        assert publicId == "-//NCBI//NCBI BlastOutput/EN"
        return 1

    def _startElementHandler(self, name, attr):
        """Found XML start tag.

        No real need of attr, BLAST DTD doesn't use them

        Arguments:
         - name -- name of the tag
         - attr -- tag attributes

        """
        method = Records._start_methods.get(name)
        if method is None:
            raise ValueError("Failed to find method for %s" % name)
        method(self, name, attr)

    def _endElementHandler(self, name):
        """Found XML end tag.

        Arguments:
         - name -- tag name

        """
        method = Records._end_methods.get(name)
        if method is None:
            raise ValueError("Failed to find method for %s" % name)
        method(self, name)

    def _characterDataHandler(self, characters):
        """Found some text.

        Arguments:
         - characters -- characters read

        """
        self._characters += characters

    def __next__(self):
        stream = self._stream
        BLOCK = self.BLOCK
        try:
            cache = self._cache
        except AttributeError:
            raise StopIteration from None
        try:
            parser = self._parser
        except AttributeError:
            parser = None
        while True:
            try:
                record = cache.popleft()
            except IndexError:  # no record ready to be returned
                pass
            else:
                return record
            # Read in another block of data from the file.
            data = stream.read(BLOCK)
            if data == b"":
                del self._cache
                if parser is not None:
                    raise ValueError("premature end of XML file")
                raise StopIteration
            try:
                parser.Parse(data, False)
            except expat.ExpatError as e:
                if parser.StartElementHandler:
                    # We saw the initial <!xml declaration, so we can be sure
                    # that we are parsing XML data. Most likely, the XML file
                    # is corrupted.
                    raise CorruptedXMLError(e) from None
                else:
                    # We have not seen the initial <!xml declaration, so
                    # probably the input data is not in XML format.
                    raise NotXMLError(e) from None

    names = ["BlastOutput",
             "BlastOutput_program",
             "BlastOutput_version",
             "BlastOutput_reference",
             "BlastOutput_db",
             "BlastOutput_query-ID",
             "BlastOutput_query-def",
             "BlastOutput_query-len",
             "BlastOutput_param",
             "Parameters",
             "Parameters_matrix",
             "Parameters_expect",
             "Parameters_sc-match",
             "Parameters_sc-mismatch",
             "Parameters_gap-open",
             "Parameters_gap-extend",
             "Parameters_filter",
             "BlastOutput_iterations",
             "Iteration",
             "Iteration_iter-num",
             "Iteration_query-ID",
             "Iteration_query-def",
             "Iteration_query-len",
             "Iteration_hits",
             "Hit",
             "Hit_num",
             "Hit_id",
             "Hit_def",
             "Hit_accession",
             "Hit_len",
             "Hit_hsps",
             "Hsp",
             "Hsp_num",
             "Hsp_bit-score",
             "Hsp_score",
             "Hsp_evalue",
             "Hsp_query-from",
             "Hsp_query-to",
             "Hsp_hit-from",
             "Hsp_hit-to",
             "Hsp_query-frame",
             "Hsp_hit-frame",
             "Hsp_identity",
             "Hsp_positive",
             "Hsp_gaps",
             "Hsp_align-len",
             "Hsp_qseq",
             "Hsp_hseq",
             "Hsp_midline",
             "Iteration_stat",
             "Statistics",
             "Statistics_db-num",
             "Statistics_db-len",
             "Statistics_hsp-len",
             "Statistics_eff-space",
             "Statistics_kappa",
             "Statistics_lambda",
             "Statistics_entropy",
            ]
    for name in names:
        method_name = name.lower().replace("-", "_")
        method = eval("_start_" + method_name)
        _start_methods[name] = method
        method = eval("_end_" + method_name)
        _end_methods[name] = method
    del name
    del names
    del method
    del method_name


def parse(source):
    """Parse an XML file containing BLAST output and return a Bio.Blast.Records object.

    This is a generator function that returns BLAST records one by one.

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
    ...     
    ...
    >>> stream.close()

    """
    try:
        stream = open(source, "rb")
    except TypeError:  # not a path, assume we received a stream
        if source.read(0) != b"":
            raise StreamModeError(
                f"BLAST output files must be opened in binary mode."
            ) from None
        stream = source

    try:
        return Records(stream)
    finally:
        if stream != source:
            stream.close()


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
    >>> stream = open("Blast/wnts.xml", "rb")  # opened in binary mode
    >>> record = Blast.read(stream)
    >>> print(record)
    >>> stream.close()

    Use the Bio.Blast.parse function if you want to read a file containing
    BLAST output for more than one query.
    """
    records = parse(source)
    try:
        record = next(records)
    except StopIteration:
        raise ValueError("No BLAST output found.") from None
    try:
        next(records)
        raise ValueError("BLAST output for more than one query found.")
    except StopIteration:
        pass
    return record
