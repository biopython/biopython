# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# Copyright 2000 by Bertrand Frottier.  All rights reserved.
# Revisions 2005-2006 copyright Michiel de Hoon
# Revisions 2006-2009 copyright Peter Cock
# Revisions 2023 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code to parse BLAST XML output, and to parse the BLAST DTD file defining the XML.

The BLAST XML DTD file is available on the NCBI site at:
https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd
"""

import os.path
from xml.parsers import expat
from collections import deque

from Bio import StreamModeError
from Bio.Blast import Record
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment, Alignments
from Bio import Entrez


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


class DTDHandler:
    """Parser for the BLAST XML DTD file."""

    def __init__(self):
        """Initialize the parser and parse the BLAST XML DTD file."""
        parser = expat.ParserCreate()
        parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        parser.ExternalEntityRefHandler = self._externalEntityRefHandler
        self.parser = parser
        self.names = []
        self._externalEntityRefHandler(None, None, "NCBI_BlastOutput.dtd", None)

    def _elementDeclHandler(self, name, model):
        self.names.append(name)

    def _externalEntityRefHandler(self, context, base, systemId, publicId):
        assert context is None
        assert base is None
        directory = Entrez.__path__[0]
        path = os.path.join(directory, "DTDs", systemId)
        parser = self.parser.ExternalEntityParserCreate(None)
        parser.ElementDeclHandler = self._elementDeclHandler
        with open(path, "rb") as stream:
            parser.ParseFile(stream)
        return 1


class XMLHandler(deque):
    """Handler for BLAST XML data."""

    BLOCK = 2048  # default block size from expat

    _start_methods = {}
    _end_methods = {}

    def __init__(self, source):
        """Initialize the expat parser."""
        self.source = source
        try:
            self._stream = open(source, "rb")
        except TypeError:  # not a path, assume we received a stream
            if source.read(0) != b"":
                raise StreamModeError(
                    "BLAST output files must be opened in binary mode."
                ) from None
            self._stream = source

        parser = expat.ParserCreate()
        parser.XmlDeclHandler = self._xmlDeclHandler
        parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self._parser = parser
        self._stream = stream

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

    def read_header(self, records):
        """Read the BLAST XML file header and store as attributes on records."""
        self._records = records
        parser = self._parser
        stream = self._stream
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
                self._records
            except AttributeError:
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

    def _start_blastoutput_mbstat(self, name, attrs):
        # ignore for now
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_param(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        self._records.param = {}

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

    def _start_parameters_include(self, name, attrs):
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

    def _start_parameters_pattern(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters_entrez_query(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_iterations(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        del self._records

    def _start_blastoutput_query_len(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_query_seq(self, name, attrs):
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

    def _start_hit(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""
        self._alignment = Alignments()

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

    def _start_hsp_pattern_from(self, name, attrs):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_pattern_to(self, name, attrs):
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

    def _start_hsp_density(self, name, attrs):
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

    def _start_iteration_message(self, name, attrs):
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
        program = self._characters
        self._program = program
        self._records.program = program
        self._characters = ""

    def _end_blastoutput_version(self, name):
        self._records.version = self._characters
        self._characters = ""

    def _end_blastoutput_reference(self, name):
        self._records.reference = self._characters
        self._characters = ""

    def _end_blastoutput_db(self, name):
        self._records.db = self._characters
        self._characters = ""

    def _end_blastoutput_query_id(self, name):
        query_id = self._characters
        self._records.query = SeqRecord(None, query_id)
        self._characters = ""

    def _end_blastoutput_query_def(self, name):
        query_def = self._characters
        self._records.query.description = query_def
        self._characters = ""

    def _end_blastoutput_query_len(self, name):
        length = int(self._characters)
        self._records.query.seq = Seq(None, length=length)
        self._characters = ""

    def _end_blastoutput_query_seq(self, name):
        seq = Seq(self._characters)
        self._characters = ""
        assert len(seq) == len(self._records.query.seq)
        self._records.query.seq = seq

    def _end_blastoutput_mbstat(self, name):
        # ignore for now
        self._characters = ""

    def _end_blastoutput_param(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_parameters(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_parameters_matrix(self, name):
        self._records.param["matrix"] = self._characters
        self._characters = ""

    def _end_parameters_expect(self, name):
        self._records.param["expect"] = float(self._characters)
        self._characters = ""

    def _end_parameters_sc_match(self, name):
        self._records.param["sc-match"] = int(self._characters)
        self._characters = ""

    def _end_parameters_sc_mismatch(self, name):
        self._records.param["sc-mismatch"] = int(self._characters)
        self._characters = ""

    def _end_parameters_include(self, name):
        self._records.param["include"] = float(self._characters)
        self._characters = ""

    def _end_parameters_gap_open(self, name):
        self._records.param["gap-open"] = int(self._characters)
        self._characters = ""

    def _end_parameters_gap_extend(self, name):
        self._records.param["gap-extend"] = int(self._characters)
        self._characters = ""

    def _end_parameters_filter(self, name):
        self._records.param["filter"] = self._characters
        self._characters = ""

    def _end_parameters_pattern(self, name, attrs):
        self._records.param["pattern"] = self._characters
        self._characters = ""

    def _end_parameters_entrez_query(self, name, attrs):
        self._records.param["entrez-query"] = self._characters
        self._characters = ""

    def _end_blastoutput_iterations(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_iteration(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        self.append(self._record)
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

    def _end_hit(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        alignment = self._alignment
        del self._alignment
        self._record.append(alignment)

    def _end_hit_num(self, name):
        num = int(self._characters)
        if num != len(self._record) + 1:
            raise ValueError(
                f"unexpected value in tag <Hit_num> (found {num}, expected {len(self._record) + 1})"
            )
        self._characters = ""

    def _end_hit_id(self, name):
        hit_id = self._characters
        self._alignment.target = SeqRecord(None, hit_id)
        self._characters = ""

    def _end_hit_def(self, name):
        description = self._characters
        self._alignment.target.description = description
        self._characters = ""

    def _end_hit_accession(self, name):
        accession = self._characters
        self._alignment.target.name = accession
        self._characters = ""

    def _end_hit_len(self, name):
        length = int(self._characters)
        self._alignment.target.seq = Seq(None, length=length)
        self._characters = ""

    def _end_hit_hsps(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_hsp_num(self, name):
        num = int(self._characters)
        if num != len(self._alignment) + 1:
            raise ValueError(
                f"unexpected value in tag <Hsp_num> (found {num}, expected {len(self._alignment) + 1})"
            )
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

    def _end_hsp_pattern_from(self, name):
        # ignore for now
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_hsp_pattern_to(self, name):
        # ignore for now
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_hsp_query_frame(self, name):
        query_frame = int(self._characters)
        if self._program in ("blastx", "tblastx") and query_frame in (
            -3,
            -2,
            -1,
            1,
            2,
            3,
        ):
            pass
        elif self._program in ("blastp", "tblastn") and query_frame == 0:
            pass
        elif self._program == "blastn" and query_frame == 1:
            pass
        else:
            raise ValueError(
                f"unexpected value {query_frame} in tag <Hsp_query-frame> for program {self._program}"
            )
        self._hsp["query-frame"] = query_frame
        self._characters = ""

    def _end_hsp_hit_frame(self, name):
        hit_frame = int(self._characters)
        if self._program in ("blastp", "blastx") and hit_frame == 0:
            pass
        elif self._program in ("tblastn", "tblastx") and hit_frame in (
            -3,
            -2,
            -1,
            1,
            2,
            3,
        ):
            pass
        elif self._program == "blastn" and hit_frame in (-1, 1):
            pass
        else:
            raise ValueError(
                f"unexpected value {hit_frame} in tag <Hsp_hit-frame> for program {self._program}"
            )
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

    def _end_hsp_density(self, name):
        self._hsp["density"] = int(self._characters)
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
        coordinates = Alignment.infer_coordinates(
            [target_seq_aligned, query_seq_aligned]
        )
        query_seq_data = query_seq_aligned.replace("-", "")
        query_frame = hsp["query-frame"]
        if self._program in ("blastx", "tblastx"):
            query_start = hsp["query-from"] - 1
            query_end = hsp["query-to"]
            assert query_end - query_start == 3 * len(query_seq_data)
        else:
            query_start = hsp["query-from"] - 1
            query_end = hsp["query-to"]
            coordinates[1, :] += query_start
            assert query_end - query_start == len(query_seq_data)
            query_seq_data = {query_start: query_seq_data}
        query_seq = Seq(query_seq_data, query_length)
        query = SeqRecord(query_seq, query_id, description=query_description)
        query.annotations["start"] = query_start
        query.annotations["end"] = query_end
        query.annotations["frame"] = query_frame
        target = self._alignment.target
        target_id = target.id
        target_name = target.name
        target_description = target.description
        target_length = len(target.seq)
        target_seq_data = target_seq_aligned.replace("-", "")
        target_frame = hsp["hit-frame"]
        if self._program in ("tblastn", "tblastx"):
            target_start = hsp["hit-from"] - 1
            target_end = hsp["hit-to"]
            assert target_end - target_start == 3 * len(target_seq_data)
            target_seq = Seq(target_seq_data, target_length)
        else:
            if target_frame == +1 or target_frame == 0:
                target_start = hsp["hit-from"] - 1
                target_end = hsp["hit-to"]
                coordinates[0, :] += target_start
            elif target_frame == -1:
                target_start = hsp["hit-to"] - 1
                target_end = hsp["hit-from"]
                target_seq_data = reverse_complement(target_seq_data)
                coordinates[0, :] = target_end - coordinates[0, :]
            assert target_end - target_start == len(target_seq_data)
            target_seq = Seq({target_start: target_seq_data}, target_length)
        target = SeqRecord(
            target_seq, target_id, target_name, description=target_description
        )
        target.annotations["start"] = target_start
        target.annotations["end"] = target_end
        target.annotations["frame"] = target_frame
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
        self._alignment.append(alignment)

    def _end_iteration_stat(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_iteration_message(self, name):
        self._record.message = self._characters
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
        """Handle the DTD declaration."""
        assert context is None
        assert base is None
        if systemId not in (
            "NCBI_BlastOutput.dtd",
            "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd",
        ):
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
        method = XMLHandler._start_methods.get(name)
        if method is None:
            raise ValueError("Failed to find method for %s" % name)
        method(self, name, attr)

    def _endElementHandler(self, name):
        """Found XML end tag.

        Arguments:
         - name -- tag name

        """
        method = XMLHandler._end_methods.get(name)
        if method is None:
            raise ValueError("Failed to find method for %s" % name)
        method(self, name)

    def _characterDataHandler(self, characters):
        """Found some text.

        Arguments:
         - characters -- characters read

        """
        self._characters += characters

    def __iter__(self):
        return self

    def __next__(self):
        try:
            stream = self._stream
        except AttributeError:
            raise StopIteration from None
        try:
            parser = self._parser
        except AttributeError:
            parser = None
        BLOCK = self.BLOCK
        while True:
            try:
                record = self.popleft()
            except IndexError:  # no record ready to be returned
                pass
            else:
                return record
            # Read in another block of data from the file.
            data = stream.read(BLOCK)
            if data == b"":
                del self._stream
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

    for name in DTDHandler().names:
        method_name = name.lower().replace("-", "_")
        _start_methods[name] = eval("_start_" + method_name)
        _end_methods[name] = eval("_end_" + method_name)
    del name
    del method_name
