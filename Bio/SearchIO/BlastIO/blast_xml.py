# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Bio.SearchIO parser for BLAST+ XML output formats."""
# for more info: http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd

import re
import warnings
from itertools import chain
from xml.etree import ElementTree
from xml.sax.saxutils import XMLGenerator, escape

from Bio import BiopythonParserWarning
from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment

__all__ = ("BlastXmlParser", "BlastXmlIndexer", "BlastXmlWriter")


# element - optional qresult attribute name mapping
_ELEM_QRESULT_OPT = {
    "Statistics_db-num": ("stat_db_num", int),
    "Statistics_db-len": ("stat_db_len", int),
    "Statistics_eff-space": ("stat_eff_space", float),
    "Statistics_hsp-len": ("stat_hsp_len", int),
    "Statistics_kappa": ("stat_kappa", float),
    "Statistics_lambda": ("stat_lambda", float),
    "Statistics_entropy": ("stat_entropy", float),
}
# element - hit attribute name mapping
_ELEM_HIT = {
    # 'Hit_def': ('description', str),   # not set by this dict
    "Hit_accession": ("accession", str),
    "Hit_len": ("seq_len", int),
}
# element - hsp attribute name mapping
_ELEM_HSP = {
    "Hsp_bit-score": ("bitscore", float),
    "Hsp_score": ("bitscore_raw", int),
    "Hsp_evalue": ("evalue", float),
    "Hsp_identity": ("ident_num", int),
    "Hsp_positive": ("pos_num", int),
    "Hsp_gaps": ("gap_num", int),
    "Hsp_density": ("density", float),
}
# element - fragment attribute name mapping
_ELEM_FRAG = {
    "Hsp_query-from": ("query_start", int),
    "Hsp_query-to": ("query_end", int),
    "Hsp_hit-from": ("hit_start", int),
    "Hsp_hit-to": ("hit_end", int),
    "Hsp_query-frame": ("query_frame", int),
    "Hsp_hit-frame": ("hit_frame", int),
    "Hsp_align-len": ("aln_span", int),
    "Hsp_pattern-from": ("pattern_start", int),
    "Hsp_pattern-to": ("pattern_end", int),
    "Hsp_hseq": ("hit", str),
    "Hsp_qseq": ("query", str),
}
# dictionary for mapping tag name and meta key name
_ELEM_META = {
    "BlastOutput_db": ("target", str),
    "BlastOutput_program": ("program", str),
    "BlastOutput_version": ("version", str),
    "BlastOutput_reference": ("reference", str),
    "Parameters_expect": ("param_evalue_threshold", float),
    "Parameters_entrez-query": ("param_entrez_query", str),
    "Parameters_filter": ("param_filter", str),
    "Parameters_gap-extend": ("param_gap_extend", int),
    "Parameters_gap-open": ("param_gap_open", int),
    "Parameters_include": ("param_include", str),
    "Parameters_matrix": ("param_matrix", str),
    "Parameters_pattern": ("param_pattern", str),
    "Parameters_sc-match": ("param_score_match", int),
    "Parameters_sc-mismatch": ("param_score_mismatch", int),
}
# these are fallback tags that store information on the first query
# outside the <Iteration> tag
# only used if query_{ID,def,len} is not found in <Iteration>
# (seen in legacy Blast <2.2.14)
_ELEM_QRESULT_FALLBACK = {
    "BlastOutput_query-ID": ("id", str),
    "BlastOutput_query-def": ("description", str),
    "BlastOutput_query-len": ("len", str),
}
# element-attribute maps, for writing
_WRITE_MAPS = {
    "preamble": (
        ("program", "program"),
        ("version", "version"),
        ("reference", "reference"),
        ("db", "target"),
        ("query-ID", "id"),
        ("query-def", "description"),
        ("query-len", "seq_len"),
        ("param", None),
    ),
    "param": (
        ("matrix", "param_matrix"),
        ("expect", "param_evalue_threshold"),
        ("sc-match", "param_score_match"),
        ("sc-mismatch", "param_score_mismatch"),
        ("gap-open", "param_gap_open"),
        ("gap-extend", "param_gap_extend"),
        ("filter", "param_filter"),
        ("pattern", "param_pattern"),
        ("entrez-query", "param_entrez_query"),
    ),
    "qresult": (
        ("query-ID", "id"),
        ("query-def", "description"),
        ("query-len", "seq_len"),
    ),
    "stat": (
        ("db-num", "stat_db_num"),
        ("db-len", "stat_db_len"),
        ("hsp-len", "stat_hsp_len"),
        ("eff-space", "stat_eff_space"),
        ("kappa", "stat_kappa"),
        ("lambda", "stat_lambda"),
        ("entropy", "stat_entropy"),
    ),
    "hit": (
        ("id", "id"),
        ("def", "description"),
        ("accession", "accession"),
        ("len", "seq_len"),
    ),
    "hsp": (
        ("bit-score", "bitscore"),
        ("score", "bitscore_raw"),
        ("evalue", "evalue"),
        ("query-from", "query_start"),
        ("query-to", "query_end"),
        ("hit-from", "hit_start"),
        ("hit-to", "hit_end"),
        ("pattern-from", "pattern_start"),
        ("pattern-to", "pattern_end"),
        ("query-frame", "query_frame"),
        ("hit-frame", "hit_frame"),
        ("identity", "ident_num"),
        ("positive", "pos_num"),
        ("gaps", "gap_num"),
        ("align-len", "aln_span"),
        ("density", "density"),
        ("qseq", "query"),
        ("hseq", "hit"),
        ("midline", None),
    ),
}
# optional elements, based on the DTD
_DTD_OPT = (
    "BlastOutput_query-seq",
    "BlastOutput_mbstat",
    "Iteration_query-def",
    "Iteration_query-len",
    "Iteration-hits",
    "Iteration_stat",
    "Iteration_message",
    "Parameters_matrix",
    "Parameters_include",
    "Parameters_sc-match",
    "Parameters_sc-mismatch",
    "Parameters_filter",
    "Parameters_pattern",
    "Parameters_entrez-query",
    "Hit_hsps",
    "Hsp_pattern-from",
    "Hsp_pattern-to",
    "Hsp_query-frame",
    "Hsp_hit-frame",
    "Hsp_identity",
    "Hsp_positive",
    "Hsp_gaps",
    "Hsp_align-len",
    "Hsp_density",
    "Hsp_midline",
)

# compile RE patterns
# for capturing BLAST version
_RE_VERSION = re.compile(r"\d+\.\d+\.\d+\+?")
# for splitting ID-description pairs
_RE_ID_DESC_PAIRS_PATTERN = re.compile(" +>")
# for splitting ID and description (must be used with maxsplit = 1)
_RE_ID_DESC_PATTERN = re.compile(" +")


def _extract_ids_and_descs(raw_id, raw_desc):
    """Extract IDs, descriptions, and raw ID from raw values (PRIVATE).

    Given values of the ``Hit_id`` and ``Hit_def`` elements, this function
    returns a tuple of three elements: all IDs, all descriptions, and the
    BLAST-generated ID. The BLAST-generated ID is set to ``None`` if no
    BLAST-generated IDs are present.

    """
    ids = []
    descs = []

    blast_gen_id = raw_id
    if raw_id.startswith("gnl|BL_ORD_ID|"):
        id_desc_line = raw_desc
    else:
        id_desc_line = raw_id + " " + raw_desc

    # create a list of lists, each list containing an ID and description
    # or just an ID, if description is not present
    id_desc_pairs = [
        re.split(_RE_ID_DESC_PATTERN, x, 1)
        for x in re.split(_RE_ID_DESC_PAIRS_PATTERN, id_desc_line)
    ]
    # make sure empty descriptions are added as empty strings
    # also, we return lists for compatibility reasons between Py2 and Py3
    for pair in id_desc_pairs:
        if len(pair) != 2:
            pair.append("")
        ids.append(pair[0])
        descs.append(pair[1])

    return (ids, descs, blast_gen_id)


class BlastXmlParser:
    """Parser for the BLAST XML format."""

    def __init__(self, handle, use_raw_query_ids=False, use_raw_hit_ids=False):
        """Initialize the class."""
        self.xml_iter = iter(ElementTree.iterparse(handle, events=("start", "end")))
        self._use_raw_query_ids = use_raw_query_ids
        self._use_raw_hit_ids = use_raw_hit_ids
        self._meta, self._fallback = self._parse_preamble()

    def __iter__(self):
        """Iterate over BlastXmlParser object yields query results."""
        yield from self._parse_qresult()

    def _parse_preamble(self):
        """Parse all tag data prior to the first query result (PRIVATE)."""
        # dictionary for containing all information prior to the first query
        meta = {}
        # dictionary for fallback information
        fallback = {}

        # parse the preamble part (anything prior to the first result)
        for event, elem in self.xml_iter:
            # get the tag values, cast appropriately, store into meta
            if event == "end" and elem.tag in _ELEM_META:
                attr_name, caster = _ELEM_META[elem.tag]

                if caster is not str:
                    meta[attr_name] = caster(elem.text)
                else:
                    meta[attr_name] = elem.text

                # delete element after we finish parsing it
                elem.clear()
                continue
            # capture fallback values
            # these are used only if the first <Iteration> does not have any
            # ID, ref, or len.
            elif event == "end" and elem.tag in _ELEM_QRESULT_FALLBACK:
                attr_name, caster = _ELEM_QRESULT_FALLBACK[elem.tag]

                if caster is not str:
                    fallback[attr_name] = caster(elem.text)
                else:
                    fallback[attr_name] = elem.text

                elem.clear()
                continue

            if event == "start" and elem.tag == "Iteration":
                break

        # we only want the version number, sans the program name or date
        if meta.get("version") is not None:
            meta["version"] = re.search(_RE_VERSION, meta["version"]).group(0)

        return meta, fallback

    def _parse_qresult(self):
        """Parse query results (PRIVATE)."""
        # parse the queries
        for event, qresult_elem in self.xml_iter:
            # </Iteration> marks the end of a single query
            # which means we can process it
            if event == "end" and qresult_elem.tag == "Iteration":

                # we'll use the following schema
                # <!ELEMENT Iteration (
                #        Iteration_iter-num,
                #        Iteration_query-ID?,
                #        Iteration_query-def?,
                #        Iteration_query-len?,
                #        Iteration_hits?,
                #        Iteration_stat?,
                #        Iteration_message?)>

                # assign query attributes with fallbacks
                query_id = qresult_elem.findtext("Iteration_query-ID")
                if query_id is None:
                    query_id = self._fallback["id"]

                query_desc = qresult_elem.findtext("Iteration_query-def")
                if query_desc is None:
                    query_desc = self._fallback["description"]

                query_len = qresult_elem.findtext("Iteration_query-len")
                if query_len is None:
                    query_len = self._fallback["len"]

                blast_query_id = query_id
                # handle blast searches against databases with Blast's IDs
                # 'Query_' marks the beginning of a BLAST+-generated ID,
                # 'lcl|' marks the beginning of a BLAST legacy-generated ID
                if not self._use_raw_query_ids and (
                    query_id.startswith("Query_") or query_id.startswith("lcl|")
                ):
                    # store the Blast-generated query ID
                    id_desc = query_desc.split(" ", 1)
                    query_id = id_desc[0]
                    try:
                        query_desc = id_desc[1]
                    except IndexError:
                        query_desc = ""

                hit_list, key_list = [], []
                for hit in self._parse_hit(
                    qresult_elem.find("Iteration_hits"), query_id
                ):
                    if hit:
                        # need to keep track of hit IDs, since there could be duplicates,
                        if hit.id in key_list:
                            warnings.warn(
                                "Renaming hit ID %r to a BLAST-generated ID "
                                "%r since the ID was already matched "
                                "by your query %r. Your BLAST database "
                                "may contain duplicate entries."
                                % (hit.id, hit.blast_id, query_id),
                                BiopythonParserWarning,
                            )
                            # fallback to Blast-generated IDs, if the ID is already present
                            # and restore the desc, too
                            hit.description = "%s %s" % (hit.id, hit.description)
                            hit.id = hit.blast_id
                            # and change the hit_id of the HSPs contained
                            for hsp in hit:
                                hsp.hit_id = hit.blast_id
                        else:
                            key_list.append(hit.id)

                        hit_list.append(hit)

                # create qresult and assign its attributes
                qresult = QueryResult(hit_list, query_id)
                qresult.description = query_desc
                qresult.seq_len = int(query_len)
                qresult.blast_id = blast_query_id
                for key, value in self._meta.items():
                    setattr(qresult, key, value)

                # statistics are stored in Iteration_stat's 'grandchildren' with the
                # following DTD
                # <!ELEMENT Statistics (
                #        Statistics_db-num,
                #        Statistics_db-len,
                #        Statistics_hsp-len,
                #        Statistics_eff-space,
                #        Statistics_kappa,
                #        Statistics_lambda,
                #        Statistics_entropy)>

                stat_iter_elem = qresult_elem.find("Iteration_stat")
                if stat_iter_elem is not None:
                    stat_elem = stat_iter_elem.find("Statistics")

                    for key, val_info in _ELEM_QRESULT_OPT.items():
                        value = stat_elem.findtext(key)
                        if value is not None:
                            caster = val_info[1]
                            # recast only if value is not intended to be str
                            if value is not None and caster is not str:
                                value = caster(value)
                            setattr(qresult, val_info[0], value)

                # delete element after we finish parsing it
                qresult_elem.clear()
                yield qresult

    def _parse_hit(self, root_hit_elem, query_id):
        """Yield a generator object that transforms Iteration_hits XML elements into Hit objects (PRIVATE).

        :param root_hit_elem: root element of the Iteration_hits tag.
        :type root_hit_elem: XML element tag
        :param query_id: QueryResult ID of this Hit
        :type query_id: string

        """
        # Hit level processing
        # Hits are stored in the Iteration_hits tag, with the following
        # DTD
        # <!ELEMENT Hit (
        #        Hit_num,
        #        Hit_id,
        #        Hit_def,
        #        Hit_accession,
        #        Hit_len,
        #        Hit_hsps?)>

        # feed the loop below an empty list so iteration still works
        if root_hit_elem is None:
            root_hit_elem = []

        for hit_elem in root_hit_elem:

            # BLAST sometimes mangles the sequence IDs and descriptions, so we need
            # to extract the actual values.
            raw_hit_id = hit_elem.findtext("Hit_id")
            raw_hit_desc = hit_elem.findtext("Hit_def")
            if not self._use_raw_hit_ids:
                ids, descs, blast_hit_id = _extract_ids_and_descs(
                    raw_hit_id, raw_hit_desc
                )
            else:
                ids, descs, blast_hit_id = [raw_hit_id], [raw_hit_desc], raw_hit_id

            hit_id, alt_hit_ids = ids[0], ids[1:]
            hit_desc, alt_hit_descs = descs[0], descs[1:]

            hsps = list(self._parse_hsp(hit_elem.find("Hit_hsps"), query_id, hit_id))

            hit = Hit(hsps)
            hit.description = hit_desc
            hit._id_alt = alt_hit_ids
            hit._description_alt = alt_hit_descs
            hit.blast_id = blast_hit_id

            for key, val_info in _ELEM_HIT.items():
                value = hit_elem.findtext(key)
                if value is not None:
                    caster = val_info[1]
                    # recast only if value is not intended to be str
                    if value is not None and caster is not str:
                        value = caster(value)
                    setattr(hit, val_info[0], value)

            # delete element after we finish parsing it
            hit_elem.clear()
            yield hit

    def _parse_hsp(self, root_hsp_frag_elem, query_id, hit_id):
        """Yield a generator object that transforms Hit_hsps XML elements into HSP objects (PRIVATE).

        :param root_hsp_frag_elem: the ``Hit_hsps`` tag
        :type root_hsp_frag_elem: XML element tag
        :param query_id: query ID
        :type query_id: string
        :param hit_id: hit ID
        :type hit_id: string

        """
        # Hit_hsps DTD:
        # <!ELEMENT Hsp (
        #        Hsp_num,
        #        Hsp_bit-score,
        #        Hsp_score,
        #        Hsp_evalue,
        #        Hsp_query-from,
        #        Hsp_query-to,
        #        Hsp_hit-from,
        #        Hsp_hit-to,
        #        Hsp_pattern-from?,
        #        Hsp_pattern-to?,
        #        Hsp_query-frame?,
        #        Hsp_hit-frame?,
        #        Hsp_identity?,
        #        Hsp_positive?,
        #        Hsp_gaps?,
        #        Hsp_align-len?,
        #        Hsp_density?,
        #        Hsp_qseq,
        #        Hsp_hseq,
        #        Hsp_midline?)>

        # if value is None, feed the loop below an empty list
        if root_hsp_frag_elem is None:
            root_hsp_frag_elem = []

        for hsp_frag_elem in root_hsp_frag_elem:
            coords = {}  # temporary container for coordinates
            frag = HSPFragment(hit_id, query_id)
            for key, val_info in _ELEM_FRAG.items():
                value = hsp_frag_elem.findtext(key)
                caster = val_info[1]

                # adjust 'from' and 'to' coordinates to 0-based ones
                if value is not None:
                    if key.endswith("-from") or key.endswith("-to"):
                        # store coordinates for further processing
                        coords[val_info[0]] = caster(value)
                        continue
                    # recast only if value is not intended to be str
                    elif caster is not str:
                        value = caster(value)
                    setattr(frag, val_info[0], value)

            # set the similarity characters into aln_annotation dict
            frag.aln_annotation["similarity"] = hsp_frag_elem.findtext("Hsp_midline")

            # process coordinates
            # since 'x-from' could be bigger than 'x-to', we need to figure
            # out which one is smaller/bigger since 'x_start' is always smaller
            # than 'x_end'
            for coord_type in ("query", "hit", "pattern"):
                start_type = coord_type + "_start"
                end_type = coord_type + "_end"
                try:
                    start = coords[start_type]
                    end = coords[end_type]
                except KeyError:
                    continue
                else:
                    # convert to python range and setattr
                    setattr(frag, start_type, min(start, end) - 1)
                    setattr(frag, end_type, max(start, end))

            # set molecule type, based on program
            prog = self._meta.get("program")
            if prog == "blastn":
                frag.molecule_type = "DNA"
            elif prog in ["blastp", "blastx", "tblastn", "tblastx"]:
                frag.molecule_type = "protein"

            hsp = HSP([frag])
            for key, val_info in _ELEM_HSP.items():
                value = hsp_frag_elem.findtext(key)
                caster = val_info[1]
                if value is not None:
                    if caster is not str:
                        value = caster(value)
                    setattr(hsp, val_info[0], value)
            # delete element after we finish parsing it
            hsp_frag_elem.clear()
            yield hsp


class BlastXmlIndexer(SearchIndexer):
    """Indexer class for BLAST XML output."""

    _parser = BlastXmlParser
    qstart_mark = b"<Iteration>"
    qend_mark = b"</Iteration>"
    block_size = 16384

    def __init__(self, filename, **kwargs):
        """Initialize the class."""
        SearchIndexer.__init__(self, filename)
        # TODO: better way to do this?
        iter_obj = self._parser(self._handle, **kwargs)
        self._meta, self._fallback = iter_obj._meta, iter_obj._fallback

    def __iter__(self):
        """Iterate over BlastXmlIndexer yields qstart_id, start_offset, block's length."""
        qstart_mark = self.qstart_mark
        qend_mark = self.qend_mark
        blast_id_mark = b"Query_"
        block_size = self.block_size
        handle = self._handle
        handle.seek(0)
        re_desc = re.compile(
            b"<Iteration_query-ID>(.*?)"
            br"</Iteration_query-ID>\s+?"
            b"<Iteration_query-def>"
            b"(.*?)</Iteration_query-def>"
        )
        re_desc_end = re.compile(b"</Iteration_query-def>")
        counter = 0

        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if not line:
                break
            if qstart_mark not in line:
                continue
            # The following requirements are to make supporting BGZF compressed
            # BLAST XML files simpler (avoids complex offset manipulations):
            assert line.count(qstart_mark) == 1, "XML without line breaks?"
            assert line.lstrip().startswith(qstart_mark), line
            if qend_mark in line:
                # Should cope with <Iteration>...</Iteration> on one long line
                block = line
            else:
                # Load the rest of this block up to and including </Iteration>
                block = [line]
                while line and qend_mark not in line:
                    line = handle.readline()
                    assert qstart_mark not in line, line
                    block.append(line)
                assert line.rstrip().endswith(qend_mark), line
                block = b"".join(block)
            assert block.count(qstart_mark) == 1, "XML without line breaks? %r" % block
            assert block.count(qend_mark) == 1, "XML without line breaks? %r" % block
            # Now we have a full <Iteration>...</Iteration> block, find the ID
            regx = re.search(re_desc, block)
            try:
                qstart_desc = regx.group(2)
                qstart_id = regx.group(1)
            except AttributeError:
                # use the fallback values
                assert re.search(re_desc_end, block)
                qstart_desc = self._fallback["description"].encode()
                qstart_id = self._fallback["id"].encode()
            if qstart_id.startswith(blast_id_mark):
                qstart_id = qstart_desc.split(b" ", 1)[0]
            yield qstart_id.decode(), start_offset, len(block)
            counter += 1

    def _parse(self, handle):
        """Overwrite SearchIndexer parse (PRIVATE).

        As we need to set the meta and fallback dictionaries to the parser.
        """
        generator = self._parser(handle, **self._kwargs)
        generator._meta = self._meta
        generator._fallback = self._fallback
        return next(iter(generator))

    def get_raw(self, offset):
        """Return the raw record from the file as a bytes string."""
        qend_mark = self.qend_mark
        handle = self._handle
        handle.seek(offset)

        qresult_raw = handle.readline()
        assert qresult_raw.lstrip().startswith(self.qstart_mark)
        while qend_mark not in qresult_raw:
            qresult_raw += handle.readline()
        assert qresult_raw.rstrip().endswith(qend_mark)
        assert qresult_raw.count(qend_mark) == 1
        # Note this will include any leading and trailing whitespace, in
        # general expecting "    <Iteration>\n...\n    </Iteration>\n"
        return qresult_raw


class _BlastXmlGenerator(XMLGenerator):
    """Event-based XML Generator."""

    def __init__(self, out, encoding="utf-8", indent=" ", increment=2):
        """Initialize the class."""
        XMLGenerator.__init__(self, out, encoding)
        # the indentation character
        self._indent = indent
        # nest level
        self._level = 0
        # how many indentation character should we increment per level
        self._increment = increment
        # container for names of tags with children
        self._parent_stack = []
        # determine writer method

    def startDocument(self):
        """Start the XML document."""
        self._write(
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" '
            '"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">\n'
        )

    def startElement(self, name, attrs=None, children=False):
        """Start an XML element.

        :param name: element name
        :type name: string
        :param attrs: element attributes
        :type attrs: dictionary {string: object}
        :param children: whether the element has children or not
        :type children: bool

        """
        if attrs is None:
            attrs = {}
        self.ignorableWhitespace(self._indent * self._level)
        XMLGenerator.startElement(self, name, attrs)

    def endElement(self, name):
        """End and XML element of the given name."""
        XMLGenerator.endElement(self, name)
        self._write("\n")

    def startParent(self, name, attrs=None):
        """Start an XML element which has children.

        :param name: element name
        :type name: string
        :param attrs: element attributes
        :type attrs: dictionary {string: object}

        """
        if attrs is None:
            attrs = {}
        self.startElement(name, attrs, children=True)
        self._level += self._increment
        self._write("\n")
        # append the element name, so we can end it later
        self._parent_stack.append(name)

    def endParent(self):
        """End an XML element with children."""
        # the element to end is the one on top of the stack
        name = self._parent_stack.pop()
        self._level -= self._increment
        self.ignorableWhitespace(self._indent * self._level)
        self.endElement(name)

    def startParents(self, *names):
        """Start XML elements without children."""
        for name in names:
            self.startParent(name)

    def endParents(self, num):
        """End XML elements, according to the given number."""
        for i in range(num):
            self.endParent()

    def simpleElement(self, name, content=None):
        """Create an XML element without children with the given content."""
        self.startElement(name, attrs={})
        if content:
            self.characters(content)
        self.endElement(name)

    def characters(self, content):
        """Replace quotes and apostrophe."""
        content = escape(str(content))
        for a, b in (('"', "&quot;"), ("'", "&apos;")):
            content = content.replace(a, b)
        self._write(content)


class BlastXmlWriter:
    """Stream-based BLAST+ XML Writer."""

    def __init__(self, handle, use_raw_query_ids=True, use_raw_hit_ids=True):
        """Initialize the class."""
        self.xml = _BlastXmlGenerator(handle, "utf-8")
        self._use_raw_query_ids = use_raw_query_ids
        self._use_raw_hit_ids = use_raw_hit_ids

    def write_file(self, qresults):
        """Write the XML contents to the output handle."""
        xml = self.xml
        self.qresult_counter, self.hit_counter, self.hsp_counter, self.frag_counter = (
            0,
            0,
            0,
            0,
        )

        # get the first qresult, since the preamble requires its attr values
        first_qresult = next(qresults)
        # start the XML document, set the root element, and create the preamble
        xml.startDocument()
        xml.startParent("BlastOutput")
        self._write_preamble(first_qresult)
        # and write the qresults
        xml.startParent("BlastOutput_iterations")
        self._write_qresults(chain([first_qresult], qresults))
        xml.endParents(2)
        xml.endDocument()

        return (
            self.qresult_counter,
            self.hit_counter,
            self.hsp_counter,
            self.frag_counter,
        )

    def _write_elem_block(self, block_name, map_name, obj, opt_dict=None):
        """Write sibling XML elements (PRIVATE).

        :param block_name: common element name prefix
        :type block_name: string
        :param map_name: name of mapping between element and attribute names
        :type map_name: string
        :param obj: object whose attribute value will be used
        :type obj: object
        :param opt_dict: custom element-attribute mapping
        :type opt_dict: dictionary {string: string}

        """
        if opt_dict is None:
            opt_dict = {}
        for elem, attr in _WRITE_MAPS[map_name]:
            elem = block_name + elem
            try:
                content = str(getattr(obj, attr))
            except AttributeError:
                # ensure attrs that is not present is optional
                if elem not in _DTD_OPT:
                    raise ValueError(
                        "Element %r (attribute %r) not found" % (elem, attr)
                    )
            else:
                # custom element-attribute mapping, for fallback values
                if elem in opt_dict:
                    content = opt_dict[elem]
                self.xml.simpleElement(elem, content)

    def _write_preamble(self, qresult):
        """Write the XML file preamble (PRIVATE)."""
        xml = self.xml

        for elem, attr in _WRITE_MAPS["preamble"]:
            elem = "BlastOutput_" + elem
            if elem == "BlastOutput_param":
                xml.startParent(elem)
                self._write_param(qresult)
                xml.endParent()
                continue
            try:
                content = str(getattr(qresult, attr))
            except AttributeError:
                if elem not in _DTD_OPT:
                    raise ValueError(
                        "Element %s (attribute %s) not found" % (elem, attr)
                    )
            else:
                if elem == "BlastOutput_version":
                    content = "%s %s" % (qresult.program.upper(), qresult.version)
                elif qresult.blast_id:
                    if elem == "BlastOutput_query-ID":
                        content = qresult.blast_id
                    elif elem == "BlastOutput_query-def":
                        content = " ".join([qresult.id, qresult.description]).strip()
                xml.simpleElement(elem, content)

    def _write_param(self, qresult):
        """Write the parameter block of the preamble (PRIVATE)."""
        xml = self.xml
        xml.startParent("Parameters")
        self._write_elem_block("Parameters_", "param", qresult)
        xml.endParent()

    def _write_qresults(self, qresults):
        """Write QueryResult objects into iteration elements (PRIVATE)."""
        xml = self.xml

        for num, qresult in enumerate(qresults):
            xml.startParent("Iteration")
            xml.simpleElement("Iteration_iter-num", str(num + 1))
            opt_dict = {}
            if self._use_raw_query_ids:
                query_id = qresult.blast_id
                query_desc = qresult.id + " " + qresult.description
            else:
                query_id = qresult.id
                query_desc = qresult.description

            opt_dict = {
                "Iteration_query-ID": query_id,
                "Iteration_query-def": query_desc,
            }
            self._write_elem_block("Iteration_", "qresult", qresult, opt_dict)
            # the Iteration_hits tag only has children if there are hits
            if qresult:
                xml.startParent("Iteration_hits")
                self._write_hits(qresult.hits)
                xml.endParent()
            # otherwise it's a simple element without any contents
            else:
                xml.simpleElement("Iteration_hits", "")

            xml.startParents("Iteration_stat", "Statistics")
            self._write_elem_block("Statistics_", "stat", qresult)
            xml.endParents(2)
            # there's a message if no hits is present
            if not qresult:
                xml.simpleElement("Iteration_message", "No hits found")
            self.qresult_counter += 1
            xml.endParent()

    def _write_hits(self, hits):
        """Write Hit objects (PRIVATE)."""
        xml = self.xml

        for num, hit in enumerate(hits):
            xml.startParent("Hit")
            xml.simpleElement("Hit_num", str(num + 1))
            # use custom hit_id and hit_def mapping if the hit has a
            # BLAST-generated ID
            opt_dict = {}

            if self._use_raw_hit_ids:
                hit_id = hit.blast_id
                hit_desc = " >".join(
                    [
                        "{} {}".format(x, y)
                        for x, y in zip(hit.id_all, hit.description_all)
                    ]
                )
            else:
                hit_id = hit.id
                hit_desc = hit.description + " >".join(
                    [
                        "{} {}".format(x, y)
                        for x, y in zip(hit.id_all[1:], hit.description_all[1:])
                    ]
                )

            opt_dict = {"Hit_id": hit_id, "Hit_def": hit_desc}
            self._write_elem_block("Hit_", "hit", hit, opt_dict)
            xml.startParent("Hit_hsps")
            self._write_hsps(hit.hsps)
            self.hit_counter += 1
            xml.endParents(2)

    def _write_hsps(self, hsps):
        """Write HSP objects (PRIVATE)."""
        xml = self.xml
        for num, hsp in enumerate(hsps):
            xml.startParent("Hsp")
            xml.simpleElement("Hsp_num", str(num + 1))
            for elem, attr in _WRITE_MAPS["hsp"]:
                elem = "Hsp_" + elem
                try:
                    content = self._adjust_output(hsp, elem, attr)
                # make sure any elements that is not present is optional
                # in the DTD
                except AttributeError:
                    if elem not in _DTD_OPT:
                        raise ValueError(
                            "Element %s (attribute %s) not found" % (elem, attr)
                        )
                else:
                    xml.simpleElement(elem, str(content))
            self.hsp_counter += 1
            self.frag_counter += len(hsp.fragments)
            xml.endParent()

    def _adjust_output(self, hsp, elem, attr):
        """Adjust output to mimic native BLAST+ XML as much as possible (PRIVATE)."""
        # adjust coordinates
        if attr in (
            "query_start",
            "query_end",
            "hit_start",
            "hit_end",
            "pattern_start",
            "pattern_end",
        ):
            content = getattr(hsp, attr) + 1
            if "_start" in attr:
                content = getattr(hsp, attr) + 1
            else:
                content = getattr(hsp, attr)

            # adjust for 'from' <--> 'to' flip if it's not a translated search
            # and frames are different
            # adapted from /src/algo/blast/format/blastxml_format.cpp#L216
            if hsp.query_frame != 0 and hsp.hit_frame < 0:
                if attr == "hit_start":
                    content = getattr(hsp, "hit_end")
                elif attr == "hit_end":
                    content = getattr(hsp, "hit_start") + 1

        # for seqrecord objects, we only need the sequence string
        elif elem in ("Hsp_hseq", "Hsp_qseq"):
            content = str(getattr(hsp, attr).seq)
        elif elem == "Hsp_midline":
            content = hsp.aln_annotation["similarity"]
        elif elem in ("Hsp_evalue", "Hsp_bit-score"):
            # adapted from src/algo/blast/format/blastxml_format.cpp#L138-140
            content = "%.*g" % (6, getattr(hsp, attr))
        else:
            content = getattr(hsp, attr)

        return content


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
