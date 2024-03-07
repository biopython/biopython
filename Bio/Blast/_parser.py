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
from collections import deque
from xml.parsers import expat
from typing import Dict, Callable

from Bio.Blast import Record, Hit, HSP
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.Align import Alignment
from Bio import Entrez


class DTDHandler:
    """Parser for the BLAST XML DTD file."""

    def __init__(self):
        """Initialize the parser and parse the BLAST XML DTD file."""
        parser = expat.ParserCreate()
        parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        parser.ExternalEntityRefHandler = self._externalEntityRefHandler
        self.parser = parser
        self._externalEntityRefHandler(None, None, "NCBI_BlastOutput.dtd", None)

    def _elementDeclHandler(self, name, model):
        method_name = name
        if name not in (
            "BlastOutput_query-ID",
            "BlastOutput_query-def",
            "BlastOutput_query-len",
        ):
            # Note that NCBI_BlastOutput.dtd defines both BlastOutput_query-ID
            # and Iteration_query-ID, while NCBI_BlastOutput2.xsd defines
            # query-id only, which corresponds to Iteration_query-ID.
            # Same for BlastOutput_query-def and BlastOutput_query-len.
            for prefix in ("BlastOutput_", "Iteration_", "Parameters_"):
                method_name = method_name.replace(prefix, "")
        method_name = method_name.lower().replace("-", "_")
        start_method = "_start_" + method_name
        end_method = "_end_" + method_name
        XMLHandler._start_methods[name] = getattr(XMLHandler, start_method)
        XMLHandler._end_methods[name] = getattr(XMLHandler, end_method)

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


class SchemaHandler:
    """XML Schema parser used to parse NCBI_BlastOutput2.xsd.

    The XML Schema for Blast XML2 is available from
    http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd
    """

    def __init__(self, parser):
        """Initialize the XML Schema parser."""
        self.parser = parser
        self.start_methods = {}
        self.end_methods = {}

    def _startElementHandler(self, name, attributes):
        """Found XML start tag.

        Arguments:
         - name       -- name of the tag
         - attributes -- tag attributes

        """
        namespace = "http://www.ncbi.nlm.nih.gov"
        if name == "http://www.w3.org/2001/XMLSchema include":
            filename = attributes["schemaLocation"]
            directory = Entrez.__path__[0]
            path = os.path.join(directory, "XSDs", filename)
            stream = open(path, "rb")
            parser = expat.ParserCreate(namespace_separator=" ")
            parser.StartElementHandler = self._startElementHandler
            parser.ParseFile(stream)
        elif name == "http://www.w3.org/2001/XMLSchema element":
            tag = attributes.get("name")
            if tag is None:
                return
            key = f"{namespace} {tag}"
            if tag == "BlastOutput2":
                self.start_methods[key] = XMLHandler._start_blastoutput
            elif tag in (
                "error",
                "Err",
                "code",
                "message",
                "subjects",
                "bl2seq",
                "iter-num",
                "query-masking",
            ):
                pass  # TBD
            else:
                if tag == "params":
                    tag = "param"
                elif tag == "search":
                    tag = "iterations"
                elif tag == "Search":
                    tag = "Iteration"
                elif tag == "query-title":
                    tag = "query-def"
                method_name = tag.lower().replace("-", "_")
                start_method = "_start_" + method_name
                end_method = "_end_" + method_name
                self.start_methods[key] = getattr(XMLHandler, start_method)
                self.end_methods[key] = getattr(XMLHandler, end_method)
        else:
            print("In SchemaHandler._startElementHandler for", name, attributes)


class XMLHandler:
    """Handler for BLAST XML data."""

    _start_methods: Dict[str, Callable] = {}
    _end_methods: Dict[str, Callable] = {}

    def __init__(self, parser):
        """Initialize the expat parser."""
        parser.XmlDeclHandler = self._xmlDeclHandler
        parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self._parser = parser

    def _startNamespaceDeclHandler(self, prefix, uri):
        if uri == "http://www.w3.org/2001/XMLSchema-instance":
            # This is an xml schema
            self.schema_namespace = uri
            self._parser.StartElementHandler = self._start_blastxml2

    def _endNamespaceDeclHandler(self, prefix):
        return

    def _start_blastxml2(self, name, attributes):
        """Process the XML schema (before processing the element)."""
        key = "%s schemaLocation" % self.schema_namespace
        assert name == "http://www.ncbi.nlm.nih.gov BlastXML2"
        domain, url = attributes[key].split()
        assert domain == "http://www.ncbi.nlm.nih.gov"
        filename = os.path.basename(url)
        directory = Entrez.__path__[0]
        path = os.path.join(directory, "XSDs", filename)
        stream = open(path, "rb")
        parser = expat.ParserCreate(namespace_separator=" ")
        handler = SchemaHandler(parser)
        parser.StartElementHandler = handler._startElementHandler
        with open(path, "rb") as stream:
            parser.ParseFile(stream)
        XMLHandler._start_methods = handler.start_methods
        parser = self._parser
        parser.StartElementHandler = self._startElementHandler
        parser.EndElementHandler = self._endElementHandler
        parser.CharacterDataHandler = self._characterDataHandler
        self._characters = ""

    def _start_blastoutput(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_program(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_version(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_reference(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_db(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_query_id(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_query_def(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_mbstat(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_param(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_parameters(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""
        self._records.param = {}

    def _start_matrix(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_expect(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_sc_match(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_sc_mismatch(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_include(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_gap_open(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_gap_extend(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_filter(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_cbs(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_db_gencode(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_gencode(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_bl2seq_mode(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_pattern(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_entrez_query(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_iterations(self, name, attributes):
        self._records._cache = deque()
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_blastoutput_query_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_seq(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_iteration(self, name, attributes):
        record = Record()
        self._record = record

    def _start_iter_num(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_id(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_def(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hits(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""
        self._alignments = Hit()

    def _start_hit_num(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_id(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_def(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_hsps(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_accession(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""
        self._hsp = {}

    def _start_hsp_num(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_bit_score(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_score(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_evalue(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_query_from(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_query_to(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_hit_from(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_hit_to(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_pattern_from(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_pattern_to(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_query_frame(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_hit_frame(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_identity(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_positive(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_gaps(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_align_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_density(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_qseq(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_hseq(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_midline(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_stat(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_message(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""
        self._stat = {}

    def _start_statistics_db_num(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_db_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_hsp_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_eff_space(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_kappa(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_lambda(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_statistics_entropy(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_report(self, name, attributes):
        return

    def _start_search_target(self, name, attributes):
        return

    def _start_target(self, name, attributes):
        return

    def _start_results(self, name, attributes):
        return

    def _end_blastoutput(self, name):
        assert self._characters.strip() == ""
        parser = self._parser
        parser.StartElementHandler = None
        parser.EndElementHandler = None
        parser.CharacterDataHandler = None
        del self._characters
        del self._records
        del self._parser

    def _end_blastxml2(self, name):
        return

    def _end_program(self, name):
        program = self._characters
        self._program = program
        self._records.program = program
        self._characters = ""

    def _end_version(self, name):
        self._records.version = self._characters
        self._characters = ""

    def _end_reference(self, name):
        self._records.reference = self._characters
        self._characters = ""

    def _end_db(self, name):
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

    def _end_query_seq(self, name):
        seq = Seq(self._characters)
        self._characters = ""
        assert len(seq) == len(self._records.query.seq)
        self._records.query.seq = seq

    def _end_mbstat(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        self._records.mbstat = self._stat
        del self._stat

    def _end_param(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_parameters(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_matrix(self, name):
        self._records.param["matrix"] = self._characters
        self._characters = ""

    def _end_expect(self, name):
        self._records.param["expect"] = float(self._characters)
        self._characters = ""

    def _end_sc_match(self, name):
        self._records.param["sc-match"] = int(self._characters)
        self._characters = ""

    def _end_sc_mismatch(self, name):
        self._records.param["sc-mismatch"] = int(self._characters)
        self._characters = ""

    def _end_include(self, name):
        self._records.param["include"] = float(self._characters)
        self._characters = ""

    def _end_gap_open(self, name):
        self._records.param["gap-open"] = int(self._characters)
        self._characters = ""

    def _end_gap_extend(self, name):
        self._records.param["gap-extend"] = int(self._characters)
        self._characters = ""

    def _end_filter(self, name):
        self._records.param["filter"] = self._characters
        self._characters = ""

    def _end_cbs(self, name):
        self._records.param["cbs"] = int(self._characters)
        self._characters = ""

    def _end_db_gencode(self, name):
        self._records.param["db-gencode"] = int(self._characters)
        self._characters = ""

    def _end_bl2seq_mode(self, name):
        self._records.param["bl2seq-mode"] = int(self._characters)
        self._characters = ""

    def _end_query_gencode(self, name):
        self._records.param["query-gencode"] = int(self._characters)
        self._characters = ""

    def _end_pattern(self, name):
        self._records.param["pattern"] = self._characters
        self._characters = ""

    def _end_entrez_query(self, name):
        self._records.param["entrez-query"] = self._characters
        self._characters = ""

    def _end_iterations(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_iteration(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        self._records._cache.append(self._record)
        del self._record

    def _end_iter_num(self, name):
        self._record.num = int(self._characters)
        self._characters = ""

    def _end_query_id(self, name):
        query_id = self._characters
        self._record.query = SeqRecord(None, query_id)
        self._characters = ""

    def _end_query_def(self, name):
        query_def = self._characters
        self._record.query.description = query_def
        self._characters = ""

    def _end_query_len(self, name):
        length = int(self._characters)
        self._record.query.seq = Seq(None, length=length)
        self._characters = ""

    def _end_hits(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_hit(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        hit = self._alignments
        del self._alignments
        self._record.append(hit)

    def _end_hit_num(self, name):
        self._alignments.num = int(self._characters)
        self._characters = ""

    def _end_hit_id(self, name):
        hit_id = self._characters
        self._alignments.target = SeqRecord(None, hit_id)
        self._characters = ""

    def _end_hit_def(self, name):
        description = self._characters
        self._alignments.target.description = description
        self._characters = ""

    def _end_hit_accession(self, name):
        accession = self._characters
        self._alignments.target.name = accession
        self._characters = ""

    def _end_hit_len(self, name):
        length = int(self._characters)
        self._alignments.target.seq = Seq(None, length=length)
        self._characters = ""

    def _end_hit_hsps(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_hsp_num(self, name):
        self._hsp["num"] = int(self._characters)
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
        self._hsp["pattern-from"] = int(self._characters)
        self._characters = ""

    def _end_hsp_pattern_to(self, name):
        self._hsp["pattern-to"] = int(self._characters)
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
        elif self._program == "blastp" and query_frame in (0, 1):
            pass
        elif self._program == "tblastn" and query_frame == 0:
            pass
        elif self._program in ("blastn", "megablast") and query_frame == 1:
            pass
        else:
            raise ValueError(
                f"unexpected value {query_frame} in tag <Hsp_query-frame> for program {self._program}"
            )
        self._hsp["query-frame"] = query_frame
        self._characters = ""

    def _end_hsp_hit_frame(self, name):
        hit_frame = int(self._characters)
        if self._program in "blastp" and hit_frame in (0, 1):
            pass
        elif self._program == "blastx" and hit_frame == 0:
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
        elif self._program in ("blastn", "megablast") and hit_frame in (-1, 1):
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
        if query is None:
            query = self._records.query
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
        query = SeqRecord(None, query_id, description=query_description)
        query_start = hsp["query-from"] - 1
        query_end = hsp["query-to"]
        if self._program in ("blastx", "tblastx"):
            assert query_end - query_start == 3 * len(query_seq_data)
            location = SimpleLocation(0, len(query_seq_data))
            coded_by = f"{query_id}:{hsp['query-from']}..{hsp['query-to']}"
            if query_frame > 0:
                assert query_start % 3 == query_frame - 1
            elif query_frame < 0:
                assert (query_length - query_end) % 3 == -query_frame - 1
                coded_by = f"complement({coded_by})"
            qualifiers = {"coded_by": coded_by}
            feature = SeqFeature(location, type="CDS", qualifiers=qualifiers)
            query.features.append(feature)
        else:
            coordinates[1, :] += query_start
            assert query_end - query_start == len(query_seq_data)
            query_seq_data = {query_start: query_seq_data}
        query.seq = Seq(query_seq_data, query_length)
        target = self._alignments.target
        target_id = target.id
        target_name = target.name
        target_description = target.description
        target_length = len(target.seq)
        target_seq_data = target_seq_aligned.replace("-", "")
        target_frame = hsp["hit-frame"]
        target = SeqRecord(None, target_id, target_name, description=target_description)
        if self._program in ("tblastn", "tblastx"):
            target_start = hsp["hit-from"] - 1
            target_end = hsp["hit-to"]
            assert target_end - target_start == 3 * len(target_seq_data)
            target_seq = Seq(target_seq_data, target_length)
            location = SimpleLocation(0, target_length)
            coded_by = f"{target_id}:{hsp['hit-from']}..{hsp['hit-to']}"
            if target_frame > 0:
                assert target_start % 3 == target_frame - 1
            elif query_frame < 0:
                assert (target_length - target_end) % 3 == -target_frame - 1
                coded_by = f"complement({coded_by})"
            qualifiers = {"coded_by": coded_by}
            feature = SeqFeature(location, type="CDS", qualifiers=qualifiers)
            target.features.append(feature)
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
            target_seq_data = {target_start: target_seq_data}
        target.seq = Seq(target_seq_data, target_length)
        sequences = [target, query]
        alignment = HSP(sequences, coordinates)
        alignment.num = hsp["num"]
        alignment.score = hsp["score"]
        annotations = {}
        annotations["bit score"] = hsp["bit-score"]
        annotations["evalue"] = hsp["evalue"]
        annotations["identity"] = hsp["identity"]
        annotations["positive"] = hsp["positive"]
        try:
            annotations["gaps"] = hsp["gaps"]
        except KeyError:  # missing in megablast
            pass
        annotations["midline"] = hsp["midline"]
        alignment.annotations = annotations
        self._alignments.append(alignment)

    def _end_stat(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        self._record.stat = self._stat
        del self._stat

    def _end_message(self, name):
        self._record.message = self._characters
        self._characters = ""

    def _end_statistics(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_statistics_db_num(self, name):
        self._stat["db-num"] = int(self._characters)
        self._characters = ""

    def _end_statistics_db_len(self, name):
        self._stat["db-len"] = int(self._characters)
        self._characters = ""

    def _end_statistics_hsp_len(self, name):
        self._stat["hsp-len"] = int(self._characters)
        self._characters = ""

    def _end_statistics_eff_space(self, name):
        self._stat["eff-space"] = float(self._characters)
        self._characters = ""

    def _end_statistics_kappa(self, name):
        self._stat["kappa"] = float(self._characters)
        self._characters = ""

    def _end_statistics_lambda(self, name):
        self._stat["lambda"] = float(self._characters)
        self._characters = ""

    def _end_statistics_entropy(self, name):
        self._stat["entropy"] = float(self._characters)
        self._characters = ""

    def _end_report(self, name):
        return

    def _end_search_target(self, name):
        return

    def _end_target(self, name):
        return

    def _end_results(self, name):
        return

    def _xmlDeclHandler(self, version, encoding, standalone):
        parser = self._parser
        parser.ExternalEntityRefHandler = self._externalEntityRefHandler
        parser.StartNamespaceDeclHandler = self._startNamespaceDeclHandler
        parser.EndNamespaceDeclHandler = self._endNamespaceDeclHandler
        parser.StartElementHandler = self._startElementHandler
        parser.EndElementHandler = self._endElementHandler
        parser.CharacterDataHandler = self._characterDataHandler
        self._characters = ""
        parser.XmlDeclHandler = None

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
        self._parser.ExternalEntityRefHandler = None
        return 1

    def _startElementHandler(self, name, attributes):
        """Found XML start tag.

        Arguments:
         - name       -- name of the tag
         - attributes -- tag attributes

        """
        method = XMLHandler._start_methods.get(name)
        if method is None:
            raise ValueError("Failed to find method for %s" % name)
        method(self, name, attributes)

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

    def __repr__(self):
        try:
            stream = self._stream
        except AttributeError:
            stream = None
        try:
            parser = self._parser
        except AttributeError:
            parser = None
        address = hex(id(self))
        if stream is None and parser is None:
            return f"<Bio.Blast._parser.XMLHandler object at {address} with no stream or parser>"
        elif stream is None:
            return f"<Bio.Blast._parser.XMLHandler object at {address} with parser {parser} and no stream>"
        elif parser is None:
            return f"<Bio.Blast._parser.XMLHandler object at {address} with stream {stream} and no parser>"
        else:
            return f"<Bio.Blast._parser.XMLHandler object at {address} with stream {stream} and parser {parser}>"


# Initialize XMLHandler by parsing the DTD
DTDHandler()
