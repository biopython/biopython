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

import html
import os.path
from collections import deque
from xml.parsers import expat

from Bio import Entrez
from Bio.Align import Alignment
from Bio.Blast import Hit
from Bio.Blast import HSP
from Bio.Blast import Record
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation
from Bio.SeqRecord import SeqRecord


class DTDHandler:
    """Parser for the BLAST XML DTD file."""

    def __init__(self):
        """Initialize the parser and parse the BLAST XML DTD file."""
        parser = expat.ParserCreate()
        parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        parser.ExternalEntityRefHandler = self._externalEntityRefHandler
        self.parser = parser
        self.start_methods = {}
        self.end_methods = {}
        XMLHandler._dtd_methods = self.start_methods, self.end_methods

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
            for prefix in (
                "BlastOutput_",
                "Iteration_",
                "Parameters_",
                "Statistics_",
                "Hit_",
                "Hsp_",
            ):
                method_name = method_name.replace(prefix, "")
        method_name = method_name.lower().replace("-", "_")
        start_method = "_start_" + method_name
        end_method = "_end_" + method_name
        self.start_methods[name] = getattr(XMLHandler, start_method)
        self.end_methods[name] = getattr(XMLHandler, end_method)

    def _externalEntityRefHandler(self, context, base, systemId, publicId):
        assert context is None
        assert base is None
        self.parseFile(systemId)
        return 1

    def parseFile(self, filename):
        """Parse a DTD file."""
        directory = Entrez.__path__[0]
        path = os.path.join(directory, "DTDs", filename)
        parser = self.parser.ExternalEntityParserCreate(None)
        parser.ElementDeclHandler = self._elementDeclHandler
        with open(path, "rb") as stream:
            parser.ParseFile(stream)


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
        XMLHandler._schema_methods = self.start_methods, self.end_methods

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
            parser = expat.ParserCreate(namespace_separator=" ")
            parser.StartElementHandler = self._startElementHandler
            with open(path, "rb") as stream:
                parser.ParseFile(stream)
        elif name == "http://www.w3.org/2001/XMLSchema element":
            tag = attributes.get("name")
            if tag is None:
                return
            key = f"{namespace} {tag}"
            if tag == "BlastOutput2":
                self.start_methods[key] = XMLHandler._start_blastoutput
                self.end_methods[key] = XMLHandler._end_blastoutput_xml2
            elif tag in (
                "error",
                "Err",
                "code",
                "message",
                "subjects",
                "bl2seq",
                "iter-num",
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
                elif tag == "title":
                    tag = "def"
                method_name = tag.lower().replace("-", "_")
                start_method = "_start_" + method_name
                end_method = "_end_" + method_name
                if tag == "eff-space":
                    end_method += "_xml2"
                self.start_methods[key] = getattr(XMLHandler, start_method)
                self.end_methods[key] = getattr(XMLHandler, end_method)


class _HSP_cache:
    __slots__ = (
        "num",
        "bit_score",
        "score",
        "evalue",
        "identity",
        "positive",
        "query_from",
        "query_to",
        "query_frame",
        "query_strand",
        "hit_from",
        "hit_to",
        "hit_frame",
        "hit_strand",
        "qseq",
        "hseq",
        "gaps",
        "align_len",
        "density",
        "midline",
    )


class XMLHandler:
    """Handler for BLAST XML data."""

    schema_namespace = "http://www.w3.org/2001/XMLSchema-instance"
    _schema_methods = None
    _dtd_methods = None

    def __init__(self, parser):
        """Initialize the expat parser."""
        parser.XmlDeclHandler = self._xmlDeclHandler
        parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self._parser = parser

    def _startNamespaceDeclHandler(self, prefix, uri):
        parser = self._parser
        if uri == XMLHandler.schema_namespace:
            # This is an xml schema
            parser.StartElementHandler = self._start_blastxml2

    def _endNamespaceDeclHandler(self, prefix):
        return

    def _start_blastxml2(self, name, attributes):
        """Process the XML schema (before processing the element)."""
        key = "%s schemaLocation" % XMLHandler.schema_namespace
        assert name == "http://www.ncbi.nlm.nih.gov BlastXML2"
        domain, url = attributes[key].split()
        assert domain == "http://www.ncbi.nlm.nih.gov"
        if XMLHandler._schema_methods is None:
            filename = os.path.basename(url)
            directory = Entrez.__path__[0]
            path = os.path.join(directory, "XSDs", filename)
            stream = open(path, "rb")
            parser = expat.ParserCreate(namespace_separator=" ")
            handler = SchemaHandler(parser)
            parser.StartElementHandler = handler._startElementHandler
            with open(path, "rb") as stream:
                parser.ParseFile(stream)
        self._start_methods, self._end_methods = XMLHandler._schema_methods
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

    def _start_query_masking(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_range(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_from(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_to(self, name, attributes):
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

    def _start_num(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_description(self, name, attributes):
        self._alignments.targets = []

    def _start_hitdescr(self, name, attributes):
        return

    def _start_id(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_def(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_taxid(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_sciname(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsps(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_accession(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""
        self._hsp = _HSP_cache()

    def _start_bit_score(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_score(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_evalue(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_from(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_to(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_strand(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_from(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_to(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_strand(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_pattern_from(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_pattern_to(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_query_frame(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hit_frame(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_identity(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_positive(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_gaps(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_align_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_density(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_qseq(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hseq(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_midline(self, name, attributes):
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

    def _start_db_num(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_db_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_hsp_len(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_eff_space(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_kappa(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_lambda(self, name, attributes):
        assert self._characters.strip() == ""
        self._characters = ""

    def _start_entropy(self, name, attributes):
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

    def _end_blastoutput_xml2(self, name):
        assert self._characters.strip() == ""

    def _end_blastxml2(self, name):
        self._end_blastoutput(name)

    def _end_program(self, name):
        program = self._characters
        self._program = program
        self._records.program = program
        self._characters = ""

    def _end_version(self, name):
        self._records.version = self._characters
        self._characters = ""

    def _end_reference(self, name):
        self._records.reference = html.unescape(self._characters)
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

    def _end_query_masking(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        location = self._location
        del self._location
        feature = SeqFeature(location, type="masking")
        self._record.query.features.append(feature)

    def _end_range(self, name):
        start = self._from - 1
        del self._from
        end = self._to
        del self._to
        self._location = SimpleLocation(start, end)

    def _end_from(self, name):
        self._from = int(self._characters)
        self._characters = ""

    def _end_to(self, name):
        self._to = int(self._characters)
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

    def _end_description(self, name):
        self._alignments.target = self._alignments.targets[0]

    def _end_hitdescr(self, name):
        self._alignments.targets.append(self._alignments.target)

    def _end_id(self, name):
        hit_id = self._characters
        self._alignments.target = SeqRecord(None, hit_id)
        self._characters = ""

    def _end_def(self, name):
        description = self._characters
        self._alignments.target.description = description
        self._characters = ""

    def _end_taxid(self, name):
        taxid = self._characters
        self._alignments.target.annotations["taxid"] = int(taxid)
        self._characters = ""

    def _end_sciname(self, name):
        sciname = self._characters
        self._alignments.target.annotations["sciname"] = sciname
        self._characters = ""

    def _end_accession(self, name):
        accession = self._characters
        self._alignments.target.name = accession
        self._characters = ""

    def _end_len(self, name):
        length = int(self._characters)
        seq = Seq(None, length=length)
        alignments = self._alignments
        try:
            targets = alignments.targets
        except AttributeError:
            # XML
            alignments.target.seq = seq
        else:
            # XML2
            for target in alignments.targets:
                target.seq = seq
        self._characters = ""

    def _end_hsps(self, name):
        assert self._characters.strip() == ""
        self._characters = ""

    def _end_num(self, name):
        try:
            element = self._hsp
        except AttributeError:
            element = self._alignments
        element.num = int(self._characters)
        self._characters = ""

    def _end_bit_score(self, name):
        self._hsp.bit_score = float(self._characters)
        self._characters = ""

    def _end_score(self, name):
        self._hsp.score = float(self._characters)
        self._characters = ""

    def _end_evalue(self, name):
        self._hsp.evalue = float(self._characters)
        self._characters = ""

    def _end_query_from(self, name):
        self._hsp.query_from = int(self._characters)
        self._characters = ""

    def _end_query_to(self, name):
        self._hsp.query_to = int(self._characters)
        self._characters = ""

    def _end_query_strand(self, name):
        query_strand = self._characters
        assert query_strand == "Plus"
        self._hsp.query_strand = query_strand
        self._characters = ""

    def _end_hit_from(self, name):
        self._hsp.hit_from = int(self._characters)
        self._characters = ""

    def _end_hit_to(self, name):
        self._hsp.hit_to = int(self._characters)
        self._characters = ""

    def _end_hit_strand(self, name):
        hit_strand = self._characters
        assert hit_strand in ("Plus", "Minus")
        self._hsp.hit_strand = hit_strand
        self._characters = ""

    def _end_pattern_from(self, name):
        self._hsp.pattern_from = int(self._characters)
        self._characters = ""

    def _end_pattern_to(self, name):
        self._hsp.pattern_to = int(self._characters)
        self._characters = ""

    def _end_query_frame(self, name):
        query_frame = int(self._characters)
        program = self._program
        if program in ("blastn", "megablast") and query_frame == 1:
            pass
        elif program in ("blastx", "tblastx") and query_frame in (-3, -2, -1, 1, 2, 3):
            pass
        elif program in ("blastp", "tblastn", "rpsblast") and query_frame == 0:
            pass
        else:
            raise ValueError(
                f"unexpected value {query_frame} in tag <Hsp_query-frame> for program {self._program}"
            )
        self._hsp.query_frame = query_frame
        self._characters = ""

    def _end_hit_frame(self, name):
        hit_frame = int(self._characters)
        program = self._program
        if program in ("blastn", "megablast") and hit_frame in (-1, 1):
            pass
        elif program in ("blastp", "blastx", "rpsblast") and hit_frame == 0:
            pass
        elif program in ("tblastn", "tblastx") and hit_frame in (
            -3,
            -2,
            -1,
            1,
            2,
            3,
        ):
            pass
        else:
            raise ValueError(
                f"unexpected value {hit_frame} in tag <Hsp_hit-frame> for program {self._program}"
            )
        self._hsp.hit_frame = hit_frame
        self._characters = ""

    def _end_identity(self, name):
        self._hsp.identity = int(self._characters)
        self._characters = ""

    def _end_positive(self, name):
        self._hsp.positive = int(self._characters)
        self._characters = ""

    def _end_gaps(self, name):
        self._hsp.gaps = int(self._characters)
        self._characters = ""

    def _end_align_len(self, name):
        self._hsp.align_len = int(self._characters)
        self._characters = ""

    def _end_density(self, name):
        self._hsp.density = int(self._characters)
        self._characters = ""

    def _end_qseq(self, name):
        self._hsp.qseq = self._characters
        self._characters = ""

    def _end_hseq(self, name):
        self._hsp.hseq = self._characters
        self._characters = ""

    def _end_midline(self, name):
        self._hsp.midline = self._characters
        self._characters = ""

    def _end_hsp(self, name):
        assert self._characters.strip() == ""
        self._characters = ""
        hsp = self._hsp
        del self._hsp
        program = self._program
        align_len = hsp.align_len
        query = self._record.query
        if query is None:
            query = self._records.query
        query_id = query.id
        query_description = query.description
        query_length = len(query.seq)
        query_seq_aligned = hsp.qseq.encode()
        assert len(query_seq_aligned) == align_len
        target_seq_aligned = hsp.hseq.encode()
        assert len(target_seq_aligned) == align_len
        (target_seq_data, query_seq_data), coordinates = (
            Alignment.parse_printed_alignment([target_seq_aligned, query_seq_aligned])
        )
        query = SeqRecord(None, query_id, description=query_description)
        query_start = hsp.query_from - 1
        query_end = hsp.query_to
        if program in ("blastx", "tblastx"):
            assert query_end - query_start == 3 * len(query_seq_data)
            location = SimpleLocation(0, len(query_seq_data))
            coded_by = f"{query_id}:{hsp.query_from}..{hsp.query_to}"
            query_frame = hsp.query_frame
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
            if program == "blastn":
                try:
                    query_strand = hsp.query_strand
                except AttributeError:
                    # v1 XML
                    pass
                else:
                    # v2 XML
                    assert query_strand == "Plus"
        query.seq = Seq(query_seq_data, query_length)
        target = self._alignments.target
        target_id = target.id
        target_name = target.name
        target_description = target.description
        target_length = len(target.seq)
        target = SeqRecord(None, target_id, target_name, description=target_description)
        if program in ("blastn", "megablast"):
            try:
                target_strand = hsp.hit_strand
            except AttributeError:
                # v1 XML
                target_frame = hsp.hit_frame
                if target_frame == +1:
                    target_strand = "Plus"
                elif target_frame == -1:
                    target_strand = "Minus"
            if target_strand == "Plus":
                target_start = hsp.hit_from - 1
                target_end = hsp.hit_to
                coordinates[0, :] += target_start
                assert target_end - target_start == len(target_seq_data)
                target_seq_data = {target_start: target_seq_data}
                target.seq = Seq(target_seq_data, target_length)
            elif target_strand == "Minus":
                target_start = hsp.hit_to - 1
                target_end = hsp.hit_from
                coordinates[0, :] = target_end - coordinates[0, :]
                assert target_end - target_start == len(target_seq_data)
                target_seq_data = {target_length - target_end: target_seq_data}
                seq = Seq(target_seq_data, target_length)
                target.seq = seq.reverse_complement()
        elif program in ("blastp", "blastx", "rpsblast"):
            target_start = hsp.hit_from - 1
            target_end = hsp.hit_to
            coordinates[0, :] += target_start
            assert target_end - target_start == len(target_seq_data)
            target_seq_data = {target_start: target_seq_data}
            target.seq = Seq(target_seq_data, target_length)
        elif program in ("tblastn", "tblastx"):
            target_start = hsp.hit_from - 1
            target_end = hsp.hit_to
            assert target_end - target_start == 3 * len(target_seq_data)
            location = SimpleLocation(0, target_length)
            coded_by = f"{target_id}:{hsp.hit_from}..{hsp.hit_to}"
            target_frame = hsp.hit_frame
            if target_frame >= 0:
                assert target_start % 3 == target_frame - 1
            elif target_frame < 0:
                assert (target_length - target_end) % 3 == -target_frame - 1
                coded_by = f"complement({coded_by})"
            qualifiers = {"coded_by": coded_by}
            feature = SeqFeature(location, type="CDS", qualifiers=qualifiers)
            target.features.append(feature)
            target.seq = Seq(target_seq_data, target_length)
        else:
            raise RuntimeError("Unexpected program name '%s'" % program)
        sequences = [target, query]
        alignment = HSP(sequences, coordinates)
        alignment.num = hsp.num
        alignment.score = hsp.score
        annotations = {}
        annotations["bit score"] = hsp.bit_score
        annotations["evalue"] = hsp.evalue
        annotations["identity"] = hsp.identity
        try:
            annotations["positive"] = hsp.positive
        except AttributeError:
            # missing in blastn for XML1
            pass
        try:
            annotations["gaps"] = hsp.gaps
        except AttributeError:
            # missing in legacy megablast
            pass
        annotations["midline"] = hsp.midline
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

    def _end_db_num(self, name):
        self._stat["db-num"] = int(self._characters)
        self._characters = ""

    def _end_db_len(self, name):
        self._stat["db-len"] = int(self._characters)
        self._characters = ""

    def _end_hsp_len(self, name):
        self._stat["hsp-len"] = int(self._characters)
        self._characters = ""

    def _end_eff_space(self, name):
        characters = self._characters
        if characters.isdigit():
            value = int(characters)
        else:
            value = float(characters)
        self._stat["eff-space"] = value
        self._characters = ""

    def _end_eff_space_xml2(self, name):
        self._stat["eff-space"] = int(self._characters)
        self._characters = ""

    def _end_kappa(self, name):
        self._stat["kappa"] = float(self._characters)
        self._characters = ""

    def _end_lambda(self, name):
        self._stat["lambda"] = float(self._characters)
        self._characters = ""

    def _end_entropy(self, name):
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
        if XMLHandler._dtd_methods is None:
            handler = DTDHandler()
            handler.parseFile("NCBI_BlastOutput.dtd")
            self._parser.ExternalEntityRefHandler = None
        self._start_methods, self._end_methods = XMLHandler._dtd_methods
        return 1

    def _startElementHandler(self, name, attributes):
        """Found XML start tag.

        Arguments:
         - name       -- name of the tag
         - attributes -- tag attributes

        """
        method = self._start_methods.get(name)
        if method is None:
            raise ValueError(
                "Failed to find method for %s (%s)" % (name, self._start_methods.keys())
            )
        method(self, name, attributes)

    def _endElementHandler(self, name):
        """Found XML end tag.

        Arguments:
         - name -- tag name

        """
        method = self._end_methods.get(name)
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
