# Copyright 2000 by Bertrand Frottier.  All rights reserved.
# Revisions 2005-2006 copyright Michiel de Hoon
# Revisions 2006-2009 copyright Peter Cock
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code to work with the BLAST XML output.

The BLAST XML DTD file is on the NCBI FTP site at:
ftp://ftp.ncbi.nlm.nih.gov/blast/documents/xml/NCBI_BlastOutput.dtd
"""

# fmt: off
# flake8: noqa
# mypy: ignore-errors

from xml.parsers import expat
from collections import deque

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


class Record(dict):
    pass


class Hit(dict):
    pass


class DataHandler:

    def __init__(self):
        self.parser = expat.ParserCreate()
        self.parser.XmlDeclHandler = self.xmlDeclHandler
        self.parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)

    def read_header(self, stream):
        self.done = False
        self.header = {}
        BLOCK= 2048
        BLOCK = 1220  # default block size from expat
        while self.done is False:
            data = stream.read(BLOCK)
            try:
                self.parser.Parse(data, False)
            except expat.ExpatError as e:
                if self.parser.StartElementHandler:
                    # We saw the initial <!xml declaration, so we can be sure
                    # that we are parsing XML data. Most likely, the XML file
                    # is corrupted.
                    raise CorruptedXMLError(e) from None
                else:
                    # We have not seen the initial <!xml declaration, so
                    # probably the input data is not in XML format.
                    raise NotXMLError(e) from None
        del self.done
        header = self.header
        del self.header
        return header

    def _start_blastoutput(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_blastoutput_program(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_blastoutput_version(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_blastoutput_reference(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_blastoutput_db(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_blastoutput_query_id(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_blastoutput_query_def(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_blastoutput_param(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_parameters(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""
        self.parameters = {}

    def _start_parameters_matrix(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_parameters_expect(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_parameters_gap_open(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_parameters_gap_extend(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_parameters_filter(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_blastoutput_iterations(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""
        query_len = self.query_len
        query_id = self.query_id
        query_def = self.query_def
        del self.query_id
        del self.query_def
        del self.query_len
        sequence = Seq(None, length=query_len)
        query = SeqRecord(sequence, query_id, description=query_def)
        self.header["query"] = query
        self.done = True
        self.cache = deque()

    def _start_blastoutput_query_len(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_iteration(self, name, attrs):
        record = Record()
        self.record = record

    def _start_iteration_iter_num(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_iteration_query_id(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_iteration_query_def(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_iteration_query_len(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_iteration_hits(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""
        self.hits = []

    def _start_hit(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""
        self.hit = Hit()

    def _start_hit_num(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hit_id(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hit_def(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hit_hsps(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""
        self.hsps = []

    def _start_hit_len(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hit_accession(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""
        self.hsp = {}

    def _start_hsp_num(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_bit_score(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_score(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_evalue(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_query_from(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_query_to(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_hit_from(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_hit_to(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_query_frame(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_hit_frame(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_identity(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_positive(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_gaps(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_align_len(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_qseq(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_hseq(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_hsp_midline(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_iteration_stat(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_statistics(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""
        self.statistics = {}

    def _start_statistics_db_num(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_statistics_db_len(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_statistics_hsp_len(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_statistics_eff_space(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_statistics_kappa(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_statistics_lambda(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _start_statistics_entropy(self, name, attrs):
        assert self.characters.strip() == ""
        self.characters = ""

    def _end_blastoutput(self, name):
        assert self.characters.strip() == ""
        del self.characters

    def _end_blastoutput_program(self, name):
        self.header["program"] = self.characters
        self.characters = ""

    def _end_blastoutput_version(self, name):
        self.header["version"] = self.characters
        self.characters = ""

    def _end_blastoutput_reference(self, name):
        self.header["reference"] = self.characters
        self.characters = ""

    def _end_blastoutput_db(self, name):
        self.header["db"] = self.characters
        self.characters = ""

    def _end_blastoutput_query_id(self, name):
        self.query_id = self.characters
        self.characters = ""

    def _end_blastoutput_query_def(self, name):
        self.query_def = self.characters
        self.characters = ""

    def _end_blastoutput_query_len(self, name):
        self.query_len = int(self.characters)
        self.characters = ""

    def _end_blastoutput_param(self, name):
        assert self.characters.strip() == ""
        self.characters = ""

    def _end_parameters(self, name):
        assert self.characters.strip() == ""
        self.characters = ""
        self.header["param"] = self.parameters
        del self.parameters

    def _end_parameters_matrix(self, name):
        self.parameters["matrix"] = self.characters
        self.characters = ""

    def _end_parameters_expect(self, name):
        self.parameters["expect"] = float(self.characters)
        self.characters = ""

    def _end_parameters_gap_open(self, name):
        self.parameters["gap-open"] = int(self.characters)
        self.characters = ""

    def _end_parameters_gap_extend(self, name):
        self.parameters["gap-extend"] = int(self.characters)
        self.characters = ""

    def _end_parameters_filter(self, name):
        self.parameters["filter"] = self.characters
        self.characters = ""

    def _end_blastoutput_iterations(self, name):
        assert self.characters.strip() == ""
        self.characters = ""

    def _end_iteration(self, name):
        assert self.characters.strip() == ""
        self.characters = ""
        self.cache.append(self.record)
        del self.record

    def _end_iteration_iter_num(self, name):
        self.record["iter-num"] = int(self.characters)
        self.characters = ""

    def _end_iteration_query_id(self, name):
        self.record["query-ID"] = self.characters
        self.characters = ""

    def _end_iteration_query_def(self, name):
        self.record["query-def"] = self.characters
        self.characters = ""

    def _end_iteration_query_len(self, name):
        self.record["query-len"] = int(self.characters)
        self.characters = ""

    def _end_iteration_hits(self, name):
        assert self.characters.strip() == ""
        self.characters = ""
        hits = self.hits
        self.record.hits = hits
        del self.hits

    def _end_hit(self, name):
        assert self.characters.strip() == ""
        self.characters = ""
        hit = self.hit
        hit_id = self.hit_id
        hit_def = self.hit_def
        hit_len = self.hit_len
        hit_accession = self.hit_accession
        del self.hit_id
        del self.hit_def
        del self.hit_len
        del self.hit_accession
        sequence = Seq(None, length=hit_len)
        target = SeqRecord(sequence, hit_id, hit_accession)
        target.description = hit_def
        self.hit.target = target
        self.hits.append(hit)

    def _end_hit_num(self, name):
        num = int(self.characters)
        if num != len(self.hits) + 1:
            raise ValueError(f"unexpected value found in tag <Hit_num> (found f{num}, expected {len(self.hits) + 1})")
        self.characters = ""

    def _end_hit_id(self, name):
        self.hit_id = self.characters
        self.characters = ""

    def _end_hit_def(self, name):
        self.hit_def = self.characters
        self.characters = ""

    def _end_hit_hsps(self, name):
        assert self.characters.strip() == ""
        self.characters = ""
        hsps = self.hsps
        self.hit["hsps"] = hsps
        del self.hsps

    def _end_hit_len(self, name):
        self.hit_len = int(self.characters)
        self.characters = ""

    def _end_hit_accession(self, name):
        self.hit_accession = self.characters
        self.characters = ""

    def _end_hsp_num(self, name):
        self.hsp["num"] = int(self.characters)
        self.characters = ""

    def _end_hsp_bit_score(self, name):
        self.hsp["bit-score"] = float(self.characters)
        self.characters = ""

    def _end_hsp_score(self, name):
        self.hsp["score"] = float(self.characters)
        self.characters = ""

    def _end_hsp_evalue(self, name):
        self.hsp["evalue"] = float(self.characters)
        self.characters = ""

    def _end_hsp_query_from(self, name):
        self.hsp["query-from"] = int(self.characters)
        self.characters = ""

    def _end_hsp_query_to(self, name):
        self.hsp["query-to"] = int(self.characters)
        self.characters = ""

    def _end_hsp_hit_from(self, name):
        self.hsp["hit-from"] = int(self.characters)
        self.characters = ""

    def _end_hsp_hit_to(self, name):
        self.hsp["hit-to"] = int(self.characters)
        self.characters = ""

    def _end_hsp_query_frame(self, name):
        self.hsp["query-frame"] = int(self.characters)
        self.characters = ""

    def _end_hsp_hit_frame(self, name):
        self.hsp["hit-frame"] = int(self.characters)
        self.characters = ""

    def _end_hsp_identity(self, name):
        self.hsp["identity"] = int(self.characters)
        self.characters = ""

    def _end_hsp_positive(self, name):
        self.hsp["positive"] = int(self.characters)
        self.characters = ""

    def _end_hsp_gaps(self, name):
        self.hsp["gaps"] = int(self.characters)
        self.characters = ""

    def _end_hsp_align_len(self, name):
        self.hsp["align-len"] = int(self.characters)
        self.characters = ""

    def _end_hsp_qseq(self, name):
        self.hsp["qseq"] = self.characters
        self.characters = ""

    def _end_hsp_hseq(self, name):
        self.hsp["hseq"] = self.characters
        self.characters = ""

    def _end_hsp_midline(self, name):
        self.hsp["midline"] = self.characters
        self.characters = ""

    def _end_hsp(self, name):
        assert self.characters.strip() == ""
        self.characters = ""
        hsp = self.hsp
        self.hsps.append(hsp)
        del self.hsp

    def _end_iteration_stat(self, name):
        assert self.characters.strip() == ""
        self.characters = ""

    def _end_statistics(self, name):
        assert self.characters.strip() == ""
        self.characters = ""
        statistics = self.statistics
        self.record["stat"] = statistics
        del self.statistics

    def _end_statistics_db_num(self, name):
        self.statistics["db-num"] = int(self.characters)
        self.characters = ""

    def _end_statistics_db_len(self, name):
        self.statistics["db-len"] = int(self.characters)
        self.characters = ""

    def _end_statistics_hsp_len(self, name):
        self.statistics["hsp-len"] = int(self.characters)
        self.characters = ""

    def _end_statistics_eff_space(self, name):
        self.statistics["eff-space"] = float(self.characters)
        self.characters = ""

    def _end_statistics_kappa(self, name):
        self.statistics["kappa"] = float(self.characters)
        self.characters = ""

    def _end_statistics_lambda(self, name):
        self.statistics["lambda"] = float(self.characters)
        self.characters = ""

    def _end_statistics_entropy(self, name):
        self.statistics["entropy"] = float(self.characters)
        self.characters = ""

    def xmlDeclHandler(self, version, encoding, standalone):
        self.parser.ExternalEntityRefHandler = self.externalEntityRefHandler
        self.parser.StartElementHandler = self.startElementHandler
        self.parser.EndElementHandler = self.endElementHandler
        self.parser.CharacterDataHandler = self.characterDataHandler
        self.characters = ""

    def externalEntityRefHandler(self, context, base, systemId, publicId):
        """Handle the DTD declaration.

        """
        assert context is None
        assert base is None
        if systemId != "NCBI_BlastOutput.dtd":
            raise ValueError("output from legacy BLAST program")
        assert publicId == "-//NCBI//NCBI BlastOutput/EN"
        try:
            DataHandler._start_methods
            DataHandler._end_methods
        except AttributeError:
            start_methods = {}
            end_methods = {}
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
                method = getattr(DataHandler, "_start_" + method_name)
                start_methods[name] = method
                method = getattr(DataHandler, "_end_" + method_name)
                end_methods[name] = method
            DataHandler._start_methods = start_methods
            DataHandler._end_methods = end_methods

        return 1

    def startElementHandler(self, name, attr):
        """Found XML start tag.

        No real need of attr, BLAST DTD doesn't use them

        Arguments:
         - name -- name of the tag
         - attr -- tag attributes

        """
        method = DataHandler._start_methods.get(name)
        if method is None:
            raise ValueError("Failed to find method for %s" % name)
        method(self, name, attr)

    def endElementHandler(self, name):
        """Found XML end tag.

        Arguments:
         - name -- tag name

        """
        method = DataHandler._end_methods.get(name)
        if method is None:
            raise ValueError("Failed to find method for %s" % name)
        method(self, name)

    def characterDataHandler(self, characters):
        """Found some text.

        Arguments:
         - characters -- characters read

        """
        self.characters += characters

    def read_next_record(self, stream):
        BLOCK = 2048  # default block size from expat
        BLOCK = 1220  # default block size from expat
        try:
            cache = self.cache
        except AttributeError:
            return None
        while True:
            try:
                record = cache.popleft()
            except IndexError:  # no record in cache
                pass
            else:
                return record
            # Read in another block of data from the file.
            data = stream.read(BLOCK)
            if data == b"":
                assert len(cache) == 0
                del self.cache
                return None
            try:
                self.parser.Parse(data, False)
            except expat.ExpatError as e:
                if self.parser.StartElementHandler:
                    # We saw the initial <!xml declaration, so we can be sure
                    # that we are parsing XML data. Most likely, the XML file
                    # is corrupted.
                    raise CorruptedXMLError(e) from None
                else:
                    # We have not seen the initial <!xml declaration, so
                    # probably the input data is not in XML format.
                    raise NotXMLError(e) from None
