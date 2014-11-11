# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ XML output formats."""
# for more info: http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd

import sys
import re
import warnings
from itertools import chain
from xml.sax.saxutils import XMLGenerator, escape

from Bio import BiopythonParserWarning


# For speed try to use cElementTree rather than ElementTree
try:
    if (3, 0) <= sys.version_info[:2] <= (3, 1):
        # Workaround for bug in python 3.0 and 3.1,
        # see http://bugs.python.org/issue9257
        from xml.etree import ElementTree as ElementTree
    else:
        from xml.etree import cElementTree as ElementTree
except ImportError:
    from xml.etree import ElementTree as ElementTree


from Bio._py3k import _as_bytes, _bytes_to_string, unicode
_empty_bytes_string = _as_bytes("")

from Bio.Alphabet import generic_dna, generic_protein
from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment


__all__ = ['BlastXmlParser', 'BlastXmlIndexer', 'BlastXmlWriter']

__docformat__ = "restructuredtext en"


# element - optional qresult attribute name mapping
_ELEM_QRESULT_OPT = {
    'Statistics_db-num': ('stat_db_num', int),
    'Statistics_db-len': ('stat_db_len', int),
    'Statistics_eff-space': ('stat_eff_space', float),
    'Statistics_hsp-len': ('stat_hsp_len', int),
    'Statistics_kappa': ('stat_kappa', float),
    'Statistics_lambda': ('stat_lambda', float),
    'Statistics_entropy': ('stat_entropy', float),
}
# element - hit attribute name mapping
_ELEM_HIT = {
    # 'Hit_def': ('description', str),   # not set by this dict
    'Hit_accession': ('accession', str),
    'Hit_len': ('seq_len', int),
}
# element - hsp attribute name mapping
_ELEM_HSP = {
    'Hsp_bit-score': ('bitscore', float),
    'Hsp_score': ('bitscore_raw', int),
    'Hsp_evalue': ('evalue', float),
    'Hsp_identity': ('ident_num', int),
    'Hsp_positive': ('pos_num', int),
    'Hsp_gaps': ('gap_num', int),
    'Hsp_density': ('density', float),
}
# element - fragment attribute name mapping
_ELEM_FRAG = {
    'Hsp_query-from': ('query_start', int),
    'Hsp_query-to': ('query_end', int),
    'Hsp_hit-from': ('hit_start', int),
    'Hsp_hit-to': ('hit_end', int),
    'Hsp_query-frame': ('query_frame', int),
    'Hsp_hit-frame': ('hit_frame', int),
    'Hsp_align-len': ('aln_span', int),
    'Hsp_pattern-from': ('pattern_start', int),
    'Hsp_pattern-to': ('pattern_end', int),
    'Hsp_hseq': ('hit', str),
    'Hsp_qseq': ('query', str),
}
# dictionary for mapping tag name and meta key name
_ELEM_META = {
    'BlastOutput_db': ('target', str),
    'BlastOutput_program': ('program', str),
    'BlastOutput_version': ('version', str),
    'BlastOutput_reference': ('reference', str),
    'Parameters_expect': ('param_evalue_threshold', float),
    'Parameters_entrez-query': ('param_entrez_query', str),
    'Parameters_filter': ('param_filter', str),
    'Parameters_gap-extend': ('param_gap_extend', int),
    'Parameters_gap-open': ('param_gap_open', int),
    'Parameters_include': ('param_include', str),
    'Parameters_matrix': ('param_matrix', str),
    'Parameters_pattern': ('param_pattern', str),
    'Parameters_sc-match': ('param_score_match', int),
    'Parameters_sc-mismatch': ('param_score_mismatch', int),
}
# these are fallback tags that store information on the first query
# outside the <Iteration> tag
# only used if query_{ID,def,len} is not found in <Iteration>
# (seen in legacy Blast <2.2.14)
_ELEM_QRESULT_FALLBACK = {
    'BlastOutput_query-ID': ('id', str),
    'BlastOutput_query-def': ('description', str),
    'BlastOutput_query-len': ('len', str),
}
# element-attribute maps, for writing
_WRITE_MAPS = {
    'preamble': (
        ('program', 'program'),
        ('version', 'version'),
        ('reference', 'reference'),
        ('db', 'target'),
        ('query-ID', 'id'),
        ('query-def', 'description'),
        ('query-len', 'seq_len'),
        ('param', None),
    ),
    'param': (
        ('matrix', 'param_matrix'),
        ('expect', 'param_evalue_threshold'),
        ('sc-match', 'param_score_match'),
        ('sc-mismatch', 'param_score_mismatch'),
        ('gap-open', 'param_gap_open'),
        ('gap-extend', 'param_gap_extend'),
        ('filter', 'param_filter'),
        ('pattern', 'param_pattern'),
        ('entrez-query', 'param_entrez_query'),
    ),
    'qresult': (
        ('query-ID', 'id'),
        ('query-def', 'description'),
        ('query-len', 'seq_len'),
    ),
    'stat': (
        ('db-num', 'stat_db_num'),
        ('db-len', 'stat_db_len'),
        ('hsp-len', 'stat_hsp_len'),
        ('eff-space', 'stat_eff_space'),
        ('kappa', 'stat_kappa'),
        ('lambda', 'stat_lambda'),
        ('entropy', 'stat_entropy'),
    ),
    'hit': (
        ('id', 'id'),
        ('def', 'description'),
        ('accession', 'accession'),
        ('len', 'seq_len'),
    ),
    'hsp': (
        ('bit-score', 'bitscore'),
        ('score', 'bitscore_raw'),
        ('evalue', 'evalue'),
        ('query-from', 'query_start'),
        ('query-to', 'query_end'),
        ('hit-from', 'hit_start'),
        ('hit-to', 'hit_end'),
        ('pattern-from', 'pattern_start'),
        ('pattern-to', 'pattern_end'),
        ('query-frame', 'query_frame'),
        ('hit-frame', 'hit_frame'),
        ('identity', 'ident_num'),
        ('positive', 'pos_num'),
        ('gaps', 'gap_num'),
        ('align-len', 'aln_span'),
        ('density', 'density'),
        ('qseq', 'query'),
        ('hseq', 'hit'),
        ('midline', None),
    ),
}
# optional elements, based on the DTD
_DTD_OPT = (
    'BlastOutput_query-seq', 'BlastOutput_mbstat', 'Iteration_query-def',
    'Iteration_query-len', 'Iteration-hits', 'Iteration_stat',
    'Iteration_message', 'Parameters_matrix', 'Parameters_include',
    'Parameters_sc-match', 'Parameters_sc-mismatch', 'Parameters_filter',
    'Parameters_pattern', 'Parameters_entrez-query', 'Hit_hsps',
    'Hsp_pattern-from', 'Hsp_pattern-to', 'Hsp_query-frame', 'Hsp_hit-frame',
    'Hsp_identity', 'Hsp_positive', 'Hsp_gaps', 'Hsp_align-len', 'Hsp_density',
    'Hsp_midline',
)

# compile RE patterns
_RE_VERSION = re.compile(r'\d+\.\d+\.\d+\+?')


class BlastXmlParser(object):
    """Parser for the BLAST XML format"""

    def __init__(self, handle):
        self.xml_iter = iter(ElementTree.iterparse(handle, events=('start', 'end')))
        self._meta, self._fallback = self._parse_preamble()

    def __iter__(self):
        for qresult in self._parse_qresult():
            yield qresult

    def _parse_preamble(self):
        """Parses all tag data prior to the first query result."""
        # dictionary for containing all information prior to the first query
        meta = {}
        # dictionary for fallback information
        fallback = {}

        # parse the preamble part (anything prior to the first result)
        for event, elem in self.xml_iter:
            # get the tag values, cast appropriately, store into meta
            if event == 'end' and elem.tag in _ELEM_META:
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
            elif event == 'end' and elem.tag in _ELEM_QRESULT_FALLBACK:
                attr_name, caster = _ELEM_QRESULT_FALLBACK[elem.tag]

                if caster is not str:
                    fallback[attr_name] = caster(elem.text)
                else:
                    fallback[attr_name] = elem.text

                elem.clear()
                continue

            if event == 'start' and elem.tag == 'Iteration':
                break

        # we only want the version number, sans the program name or date
        if meta.get('version') is not None:
            meta['version'] = re.search(_RE_VERSION,
                    meta['version']).group(0)

        return meta, fallback

    def _parse_qresult(self):
        """Parses query results."""
        # parse the queries
        for event, qresult_elem in self.xml_iter:
            # </Iteration> marks the end of a single query
            # which means we can process it
            if event == 'end' and qresult_elem.tag == 'Iteration':

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
                query_id = qresult_elem.findtext('Iteration_query-ID')
                if query_id is None:
                    query_id = self._fallback['id']

                query_desc = qresult_elem.findtext('Iteration_query-def')
                if query_desc is None:
                    query_desc = self._fallback['description']

                query_len = qresult_elem.findtext('Iteration_query-len')
                if query_len is None:
                    query_len = self._fallback['len']

                # handle blast searches against databases with Blast's IDs
                # 'Query_' marks the beginning of a BLAST+-generated ID,
                # 'lcl|' marks the beginning of a BLAST legacy-generated ID
                if query_id.startswith('Query_') or query_id.startswith('lcl|'):
                    # store the Blast-generated query ID
                    blast_query_id = query_id
                    id_desc = query_desc.split(' ', 1)
                    query_id = id_desc[0]
                    try:
                        query_desc = id_desc[1]
                    except IndexError:
                        query_desc = ''
                else:
                    blast_query_id = ''

                hit_list, key_list = [], []
                for hit in self._parse_hit(qresult_elem.find('Iteration_hits'),
                        query_id):
                    if hit:
                        # need to keep track of hit IDs, since there could be duplicates,
                        if hit.id in key_list:
                            warnings.warn("Adding hit with BLAST-generated ID "
                                    "%r since hit ID %r is already present "
                                    "in query %r. Your BLAST database may contain "
                                    "duplicate entries." %
                                    (hit._blast_id, hit.id, query_id), BiopythonParserWarning)
                            # fallback to Blast-generated IDs, if the ID is already present
                            # and restore the desc, too
                            hit.description = '%s %s' % (hit.id, hit.description)
                            hit.id = hit._blast_id
                            # and change the hit_id of the HSPs contained
                            for hsp in hit:
                                hsp.hit_id = hit._blast_id
                        else:
                            key_list.append(hit.id)

                        hit_list.append(hit)

                # create qresult and assign its attributes
                qresult = QueryResult(hit_list, query_id)
                qresult.description = query_desc
                qresult.seq_len = int(query_len)
                qresult._blast_id = blast_query_id
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

                stat_iter_elem = qresult_elem.find('Iteration_stat')
                if stat_iter_elem is not None:
                    stat_elem = stat_iter_elem.find('Statistics')

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
        """Generator that transforms Iteration_hits XML elements into Hit objects.

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

            # create empty hit object
            hit_id = hit_elem.findtext('Hit_id')
            hit_desc = hit_elem.findtext('Hit_def')
            # handle blast searches against databases with Blast's IDs
            if hit_id.startswith('gnl|BL_ORD_ID|'):
                blast_hit_id = hit_id
                id_desc = hit_desc.split(' ', 1)
                hit_id = id_desc[0]
                try:
                    hit_desc = id_desc[1]
                except IndexError:
                    hit_desc = ''
            else:
                blast_hit_id = ''

            # combine primary ID and defline first before splitting
            full_id_desc = hit_id + ' ' + hit_desc
            id_descs = [(x.strip(), y.strip()) for x, y in
                    [a.split(' ', 1) for a in full_id_desc.split(' >')]]
            hit_id, hit_desc = id_descs[0]

            hsps = [hsp for hsp in
                    self._parse_hsp(hit_elem.find('Hit_hsps'),
                        query_id, hit_id)]

            hit = Hit(hsps)
            hit.description = hit_desc
            hit._id_alt = [x[0] for x in id_descs[1:]]
            hit._description_alt = [x[1] for x in id_descs[1:]]
            # blast_hit_id is only set if the hit ID is Blast-generated
            hit._blast_id = blast_hit_id

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
        """Iterator that transforms Hit_hsps XML elements into HSP objects.

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
                    if key.endswith('-from') or key.endswith('-to'):
                        # store coordinates for further processing
                        coords[val_info[0]] = caster(value)
                        continue
                    # recast only if value is not intended to be str
                    elif caster is not str:
                        value = caster(value)
                    setattr(frag, val_info[0], value)

            # set the similarity characters into aln_annotation dict
            frag.aln_annotation['similarity'] = \
                    hsp_frag_elem.findtext('Hsp_midline')

            # process coordinates
            # since 'x-from' could be bigger than 'x-to', we need to figure
            # out which one is smaller/bigger since 'x_start' is always smaller
            # than 'x_end'
            for coord_type in ('query', 'hit', 'pattern'):
                start_type = coord_type + '_start'
                end_type = coord_type + '_end'
                try:
                    start = coords[start_type]
                    end = coords[end_type]
                except KeyError:
                    continue
                else:
                    # convert to python range and setattr
                    setattr(frag, start_type, min(start, end) - 1)
                    setattr(frag, end_type, max(start, end))

            # set alphabet, based on program
            prog = self._meta.get('program')
            if prog == 'blastn':
                frag.alphabet = generic_dna
            elif prog in ['blastp', 'blastx', 'tblastn', 'tblastx']:
                frag.alphabet = generic_protein

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
    qstart_mark = _as_bytes('<Iteration>')
    qend_mark = _as_bytes('</Iteration>')
    block_size = 16384

    def __init__(self, filename):
        SearchIndexer.__init__(self, filename)
        # TODO: better way to do this?
        iter_obj = self._parser(self._handle)
        self._meta, self._fallback = iter_obj._meta, iter_obj._fallback

    def __iter__(self):
        qstart_mark = self.qstart_mark
        qend_mark = self.qend_mark
        blast_id_mark = _as_bytes('Query_')
        block_size = self.block_size
        handle = self._handle
        handle.seek(0)
        re_desc = re.compile(_as_bytes(r'<Iteration_query-ID>(.*?)'
                '</Iteration_query-ID>\s+?<Iteration_query-def>'
                '(.*?)</Iteration_query-def>'))
        re_desc_end = re.compile(_as_bytes(r'</Iteration_query-def>'))
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
                block = _empty_bytes_string.join(block)
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
                qstart_desc = _as_bytes(self._fallback['description'])
                qstart_id = _as_bytes(self._fallback['id'])
            if qstart_id.startswith(blast_id_mark):
                qstart_id = qstart_desc.split(_as_bytes(' '), 1)[0]
            yield _bytes_to_string(qstart_id), start_offset, len(block)
            counter += 1

    def _parse(self, handle):
        # overwrites SearchIndexer._parse, since we need to set the meta and
        # fallback dictionaries to the parser
        generator = self._parser(handle, **self._kwargs)
        generator._meta = self._meta
        generator._fallback = self._fallback
        return next(iter(generator))

    def get_raw(self, offset):
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

    def __init__(self, out, encoding='utf-8', indent=" ", increment=2):
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
        try:
            # this should work for all platforms except Jython
            self.write = self._write
        except AttributeError:
            # Jython uses self._out.write
            self.write = self._out.write

    def startDocument(self):
        """Starts the XML document."""
        self.write(u'<?xml version="1.0"?>\n'
                '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" '
                '"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">\n')

    def startElement(self, name, attrs={}, children=False):
        """Starts an XML element.

        :param name: element name
        :type name: string
        :param attrs: element attributes
        :type attrs: dictionary {string: object}
        :param children: whether the element has children or not
        :type children: bool

        """
        self.ignorableWhitespace(self._indent * self._level)
        XMLGenerator.startElement(self, name, attrs)

    def endElement(self, name):
        """Ends and XML element of the given name."""
        XMLGenerator.endElement(self, name)
        self.write(u'\n')

    def startParent(self, name, attrs={}):
        """Starts an XML element which has children.

        :param name: element name
        :type name: string
        :param attrs: element attributes
        :type attrs: dictionary {string: object}

        """
        self.startElement(name, attrs, children=True)
        self._level += self._increment
        self.write(u'\n')
        # append the element name, so we can end it later
        self._parent_stack.append(name)

    def endParent(self):
        """Ends an XML element with children."""
        # the element to end is the one on top of the stack
        name = self._parent_stack.pop()
        self._level -= self._increment
        self.ignorableWhitespace(self._indent * self._level)
        self.endElement(name)

    def startParents(self, *names):
        """Starts XML elements without children."""
        for name in names:
            self.startParent(name)

    def endParents(self, num):
        """Ends XML elements, according to the given number."""
        for i in range(num):
            self.endParent()

    def simpleElement(self, name, content=None):
        """Creates an XML element without children with the given content."""
        self.startElement(name, attrs={})
        if content:
            self.characters(content)
        self.endElement(name)

    def characters(self, content):
        content = escape(unicode(content))
        for a, b in ((u'"', u'&quot;'), (u"'", u'&apos;')):
            content = content.replace(a, b)
        self.write(content)


class BlastXmlWriter(object):
    """Stream-based BLAST+ XML Writer."""

    def __init__(self, handle):
        self.xml = _BlastXmlGenerator(handle, 'utf-8')

    def write_file(self, qresults):
        """Writes the XML contents to the output handle."""
        xml = self.xml
        self.qresult_counter, self.hit_counter, self.hsp_counter, \
            self.frag_counter = 0, 0, 0, 0

        # get the first qresult, since the preamble requires its attr values
        first_qresult = next(qresults)
        # start the XML document, set the root element, and create the preamble
        xml.startDocument()
        xml.startParent('BlastOutput')
        self._write_preamble(first_qresult)
        # and write the qresults
        xml.startParent('BlastOutput_iterations')
        self._write_qresults(chain([first_qresult], qresults))
        xml.endParents(2)
        xml.endDocument()

        return self.qresult_counter, self.hit_counter, self.hsp_counter, \
            self.frag_counter

    def _write_elem_block(self, block_name, map_name, obj, opt_dict={}):
        """Writes sibling XML elements.

        :param block_name: common element name prefix
        :type block_name: string
        :param map_name: name of mapping between element and attribute names
        :type map_name: string
        :param obj: object whose attribute value will be used
        :type obj: object
        :param opt_dict: custom element-attribute mapping
        :type opt_dict: dictionary {string: string}

        """
        for elem, attr in _WRITE_MAPS[map_name]:
            elem = block_name + elem
            try:
                content = str(getattr(obj, attr))
            except AttributeError:
                # ensure attrs that is not present is optional
                assert elem in _DTD_OPT, "Element %r (attribute %r) not " \
                    "found" % (elem, attr)
            else:
                # custom element-attribute mapping, for fallback values
                if elem in opt_dict:
                    content = opt_dict[elem]
                self.xml.simpleElement(elem, content)

    def _write_preamble(self, qresult):
        """Writes the XML file preamble."""
        xml = self.xml

        for elem, attr in _WRITE_MAPS['preamble']:
            elem = 'BlastOutput_' + elem
            if elem == 'BlastOutput_param':
                xml.startParent(elem)
                self._write_param(qresult)
                xml.endParent()
                continue
            try:
                content = str(getattr(qresult, attr))
            except AttributeError:
                assert elem in _DTD_OPT, "Element %s (attribute %s) not " \
                    "found" % (elem, attr)
            else:
                if elem == 'BlastOutput_version':
                    content = '%s %s' % (qresult.program.upper(),
                            qresult.version)
                elif qresult._blast_id:
                    if elem == 'BlastOutput_query-ID':
                        content = qresult._blast_id
                    elif elem == 'BlastOutput_query-def':
                        content = ' '.join([qresult.id,
                            qresult.description]).strip()
                xml.simpleElement(elem, content)

    def _write_param(self, qresult):
        """Writes the parameter block of the preamble."""
        xml = self.xml
        xml.startParent('Parameters')
        self._write_elem_block('Parameters_', 'param', qresult)
        xml.endParent()

    def _write_qresults(self, qresults):
        """Writes QueryResult objects into iteration elements."""
        xml = self.xml

        for num, qresult in enumerate(qresults):
            xml.startParent('Iteration')
            xml.simpleElement('Iteration_iter-num', str(num+1))
            opt_dict = {}
            # use custom Iteration_query-ID and Iteration_query-def mapping
            # if the query has a BLAST-generated ID
            if qresult._blast_id:
                opt_dict = {
                    'Iteration_query-ID': qresult._blast_id,
                    'Iteration_query-def': ' '.join([qresult.id,
                            qresult.description]).strip(),
                }
            self._write_elem_block('Iteration_', 'qresult', qresult, opt_dict)
            # the Iteration_hits tag only has children if there are hits
            if qresult:
                xml.startParent('Iteration_hits')
                self._write_hits(qresult.hits)
                xml.endParent()
            # otherwise it's a simple element without any contents
            else:
                xml.simpleElement('Iteration_hits', '')

            xml.startParents('Iteration_stat', 'Statistics')
            self._write_elem_block('Statistics_', 'stat', qresult)
            xml.endParents(2)
            # there's a message if no hits is present
            if not qresult:
                xml.simpleElement('Iteration_message', 'No hits found')
            self.qresult_counter += 1
            xml.endParent()

    def _write_hits(self, hits):
        """Writes Hit objects."""
        xml = self.xml

        for num, hit in enumerate(hits):
            xml.startParent('Hit')
            xml.simpleElement('Hit_num', str(num+1))
            # use custom hit_id and hit_def mapping if the hit has a
            # BLAST-generated ID
            opt_dict = {}
            if hit._blast_id:
                opt_dict = {
                    'Hit_id': hit._blast_id,
                    'Hit_def': ' '.join([hit.id, hit.description]).strip(),
                }
            self._write_elem_block('Hit_', 'hit', hit, opt_dict)
            xml.startParent('Hit_hsps')
            self._write_hsps(hit.hsps)
            self.hit_counter += 1
            xml.endParents(2)

    def _write_hsps(self, hsps):
        """Writes HSP objects."""
        xml = self.xml
        for num, hsp in enumerate(hsps):
            xml.startParent('Hsp')
            xml.simpleElement('Hsp_num', str(num+1))
            for elem, attr in _WRITE_MAPS['hsp']:
                elem = 'Hsp_' + elem
                try:
                    content = self._adjust_output(hsp, elem, attr)
                # make sure any elements that is not present is optional
                # in the DTD
                except AttributeError:
                    assert elem in _DTD_OPT, "Element %s (attribute %s) not found" \
                            % (elem, attr)
                else:
                    xml.simpleElement(elem, str(content))
            self.hsp_counter += 1
            self.frag_counter += len(hsp.fragments)
            xml.endParent()

    def _adjust_output(self, hsp, elem, attr):
        """Adjusts output to mimic native BLAST+ XML as much as possible."""

        # adjust coordinates
        if attr in ('query_start', 'query_end', 'hit_start', 'hit_end',
                'pattern_start', 'pattern_end'):
            content = getattr(hsp, attr) + 1
            if '_start' in attr:
                content = getattr(hsp, attr) + 1
            else:
                content = getattr(hsp, attr)

            # adjust for 'from' <--> 'to' flip if it's not a translated search
            # and frames are different
            # adapted from /src/algo/blast/format/blastxml_format.cpp#L216
            if hsp.query_frame != 0 and hsp.hit_frame < 0:
                if attr == 'hit_start':
                    content = getattr(hsp, 'hit_end')
                elif attr == 'hit_end':
                    content = getattr(hsp, 'hit_start') + 1

        # for seqrecord objects, we only need the sequence string
        elif elem in ('Hsp_hseq', 'Hsp_qseq'):
            content = str(getattr(hsp, attr).seq)
        elif elem == 'Hsp_midline':
            content = hsp.aln_annotation['similarity']
        elif elem in ('Hsp_evalue', 'Hsp_bit-score'):
            # adapted from src/algo/blast/format/blastxml_format.cpp#L138-140
            content = '%.*g' % (6, getattr(hsp, attr))
        else:
            content = getattr(hsp, attr)

        return content


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
