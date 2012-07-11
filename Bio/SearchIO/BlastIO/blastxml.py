# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ XML output formats."""
# for more info: http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd

import re
from itertools import chain
from xml.sax.saxutils import XMLGenerator, escape
try:
    from xml.etree import cElementTree as ET
except ImportError:
    from xml.etree import ElementTree as ET

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._objects import QueryResult, Hit, ContiguousHSP
from Bio.SearchIO._index import SearchIndexer


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
    'Hit_def': ('desc', str),
    'Hit_accession': ('acc', str),
    'Hit_len': ('seq_len', int),
}
# element - hsp attribute name mapping
_ELEM_HSP = {
    'Hsp_bit-score': ('bitscore', float),
    'Hsp_score': ('bitscore_raw', int),
    'Hsp_evalue': ('evalue', float),
    'Hsp_query-from': ('query_start', int),
    'Hsp_query-to': ('query_end', int),
    'Hsp_hit-from': ('hit_start', int),
    'Hsp_hit-to': ('hit_end', int),
    'Hsp_query-frame': ('query_frame', int),
    'Hsp_hit-frame': ('hit_frame', int),
    'Hsp_identity': ('ident_num', int),
    'Hsp_positive': ('pos_num', int),
    'Hsp_gaps': ('gap_num', int),
    'Hsp_align-len': ('ali_len', int),
    'Hsp_pattern-from': ('pattern_start', int),
    'Hsp_pattern-to': ('pattern_end', int),
    'Hsp_density': ('density', float),
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
    'BlastOutput_query-def': ('desc', str),
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
        ('query-def', 'desc'),
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
        ('query-def', 'desc'),
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
        ('def', 'desc'),
        ('accession', 'acc'),
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
        ('align-len', None),
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


class BlastXmlIterator(object):

    def __init__(self, handle):
        self.xml_iter = iter(ET.iterparse(handle, events=('start', 'end')))
        self._meta, self._fallback = self.parse_preamble()

    def __iter__(self):
        for qresult in self.parse_qresult():
            yield qresult

    def parse_preamble(self):
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
            meta['version'] = re.search(_RE_VERSION, \
                    meta['version']).group(0)

        return meta, fallback

    def parse_qresult(self):
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
                    query_desc = self._fallback['desc']

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

                # create qresult and assign its attributes
                qresult = QueryResult(query_id)
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

                for hit in self.parse_hit(qresult_elem.find('Iteration_hits'), \
                        query_id):
                    # only append the Hit object if we have HSPs
                    if hit:
                        if hit.id in qresult:
                            # fallback to Blast-generated IDs, if the ID is already present
                            # and restore the desc, too
                            hit.desc = '%s %s' % (hit.id, hit.desc)
                            hit.id = hit._blast_id
                            # and change the hit_id of the HSPs contained
                            for hsp in hit:
                                hsp.hit_id = hit._blast_id

                        qresult.append(hit)

                # delete element after we finish parsing it
                qresult_elem.clear()
                # set desc here so hsp.query.description is set as well
                qresult.desc = query_desc
                yield qresult

    def parse_hit(self, root_hit_elem, query_id):
        """Generator that transforms Iteration_hits XML elements into Hit objects.

        Arguments:
        root_hit_elem -- Element object of the Iteration_hits tag.
        query_id -- String of QueryResult ID of this Hit

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

            hit = Hit(hit_id, query_id)
            for key, val_info in _ELEM_HIT.items():
                value = hit_elem.findtext(key)
                if value is not None:
                    caster = val_info[1]
                    # recast only if value is not intended to be str
                    if value is not None and caster is not str:
                        value = caster(value)
                    setattr(hit, val_info[0], value)

            # blast_hit_id is only set if the hit ID is Blast-generated
            hit._blast_id = blast_hit_id

            for hsp in self.parse_hsp(hit_elem.find('Hit_hsps'), query_id, \
                    hit_id):
                hit.append(hsp)

            # delete element after we finish parsing it
            hit_elem.clear()
            # set desc here so hsp.hit.description is set as well
            hit.desc = hit_desc
            yield hit

    def parse_hsp(self, root_hsp_elem, query_id, hit_id):
        """Iterator that transforms Hit_hsps XML elements into HSP objects.

        Arguments:
        root_hsp_elem -- Element object of the Hit_hsps tag.
        hit_id -- String of hit ID to initialize all HSP objects.

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
        if root_hsp_elem is None:
            root_hsp_elem = []

        for hsp_elem in root_hsp_elem:
            hsp = ContiguousHSP(hit_id, query_id)
            # temporary container for coordinates
            coords = {}

            for key, val_info in _ELEM_HSP.items():
                value = hsp_elem.findtext(key)
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
                    setattr(hsp, val_info[0], value)

            # set the homology characters into alignment_annotation dict
            hsp.alignment_annotation['homology'] = \
                    hsp_elem.findtext('Hsp_midline')

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
                    setattr(hsp, start_type, min(start, end) - 1)
                    setattr(hsp, end_type, max(start, end))

            # delete element after we finish parsing it
            hsp_elem.clear()
            yield hsp


class BlastXmlIndexer(SearchIndexer):

    """Indexer class for BLAST+ XML output."""

    _parser = BlastXmlIterator
    qstart_mark = _as_bytes('<Iteration>')
    qend_mark = _as_bytes('</Iteration>')
    block_size = 16384

    def __init__(self, *args, **kwargs):
        SearchIndexer.__init__(self, *args, **kwargs)
        # TODO: better way to do this?
        iter_obj = self._parser(self._handle)
        self._meta, self._fallback = iter_obj._meta, iter_obj._fallback

    def get_offsets(self, string, sub, offset=0):
        """Retrieves the offsets of a given substring in a string."""
        idx = string.find(sub, offset)
        # yield offset as long as it's present in the string
        while idx >= 0:
            yield idx
            idx = string.find(sub, idx + 1)

    def __iter__(self):
        qstart_mark = self.qstart_mark
        qend_mark = self.qend_mark
        block_size = self.block_size
        handle = self._handle
        handle.seek(0)
        re_desc = re.compile(r'<Iteration_query-ID>(.*?)'
                '</Iteration_query-ID>\s+?<Iteration_query-def>'
                '(.*?)</Iteration_query-def>')
        re_desc_end = re.compile(r'</Iteration_query-def>')
        counter = 0

        while True:
            # read file content (every block size) and store into block
            block = handle.read(block_size)
            # for each iteration start mark
            for qstart_idx in self.get_offsets(block, qstart_mark):
                # get the id of the query (re.search gets the 1st result)
                regx = re.search(re_desc, block[qstart_idx:])
                try:
                    qstart_desc = regx.group(2)
                    qstart_id = regx.group(1)
                # handle cases where the iteration query ID tag lies in a block
                # split
                except AttributeError:
                    # if we've found the end of the Iteration_query-def element
                    # use the fallback values
                    if re.search(re_desc_end, block[qstart_idx:]):
                        qstart_desc = self._fallback['desc']
                        qstart_id = self._fallback['id']
                    # otherwise, extend the read block and retrieve the ID and
                    # desc from the extended read
                    else:
                        # extend the cached read
                        block_ext = block + handle.read(block_size)
                        regx = re.search(re_desc, block_ext[qstart_idx:])
                        qstart_desc = regx.group(2)
                        qstart_id = regx.group(1)
                        # set file pointer to the position before block_ext
                        handle.seek((counter + 1) * block_size)

                # now for getting the length
                # try finding it in the block
                qlen = block[qstart_idx:].find(qend_mark)
                # if not there, loop until query end is found
                block_ext = block
                while qlen < 0:
                    ext = handle.read(block_size)
                    if not ext:
                        raise ValueError("Query end for %r not found" % qstart_id)
                    block_ext += ext
                    qlen = block_ext[qstart_idx:].find(qend_mark)

                # adjust for read blocks
                qstart_idx = counter * block_size + qstart_idx
                qlen = len(qend_mark) + qlen
                # move pointer to original position
                handle.seek((counter + 1) * block_size)
                if qstart_id.startswith(_as_bytes('Query_')):
                    qstart_id = qstart_desc.split(' ', 1)[0]
                # yield key, offset, length
                yield _bytes_to_string(qstart_id), qstart_idx, qlen

            counter += 1

            if not block:
                break

    def get_raw(self, offset):
        qend_mark = self.qend_mark
        block_size = self.block_size
        handle = self._handle
        handle.seek(offset)
        counter = 0
        qresult_raw = ''

        while True:
            block = handle.read(block_size)

            # if we reach EOF without encountering any query end mark
            if not block:
                raise ValueError("Query end not found")

            qresult_raw += block
            qend_idx = qresult_raw.find(qend_mark)

            # if a match is found, return the raw qresult string
            if qend_idx > 0:
                return qresult_raw[:qend_idx + len(qend_mark)]
            # otherwise, increment the counter and go on to the next iteration
            counter += 1

    def get(self, offset):
        qresult = SearchIndexer.get(self, offset)
        # set qresult object attributes with values from preamble
        for key, value in self._meta.items():
            setattr(qresult, key, value)
        return qresult


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

    def startDocument(self):
        """Starts the XML document."""
        self._write('<?xml version="1.0"?>\n'
                '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" '
                '"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">\n')

    def startElement(self, name, attrs={}, children=False):
        """Starts an XML element.

        Arguments:
        name -- String of element name.
        attrs -- Dictionary of element attributes.
        children -- Boolean, whether the element has children or not.

        """
        self.ignorableWhitespace(self._indent * self._level)
        XMLGenerator.startElement(self, name, attrs)

    def endElement(self, name):
        """Ends and XML element of the given name."""
        XMLGenerator.endElement(self, name)
        self._write('\n')

    def startParent(self, name, attrs={}):
        """Starts an XML element which has children.

        Arguments:
        name -- String of element name.
        attrs -- Dictionary of element attributes.

        """
        self.startElement(name, attrs, children=True)
        self._level += self._increment
        self._write('\n')
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
        # apply sax's filter first, then ours
        content = escape(content)
        for a, b in (('"', '&quot;'), ("'", '&apos;')):
            content = content.replace(a, b)
        self._write(content)


class BlastXmlWriter(object):

    """Stream-based BLAST+ XML Writer."""

    def __init__(self, handle):
        self.xml = _BlastXmlGenerator(handle, 'utf-8')

    def write_file(self, qresults):
        """Writes the XML contents to the output handle."""
        xml = self.xml
        self.qresult_counter, self.hit_counter, self.hsp_counter = 0, 0, 0

        # get the first qresult, since the preamble requires its attr values
        first_qresult = qresults.next()
        # start the XML document, set the root element, and create the preamble
        xml.startDocument()
        xml.startParent('BlastOutput')
        self.write_preamble(first_qresult)
        # and write the qresults
        xml.startParent('BlastOutput_iterations')
        self.write_qresults(chain([first_qresult], qresults))
        xml.endParents(2)
        xml.endDocument()

        return self.qresult_counter, self.hit_counter, self.hsp_counter

    def write_elem_block(self, block_name, map_name, obj, opt_dict={}):
        """Writes sibling XML elements.

        Arguments:
        block_name -- String of common element name prefix.
        map_name -- Dictionary name to for mapping element and attribute names.
        obj -- Object whose attribute values will be used.

        """
        for elem, attr in _WRITE_MAPS[map_name]:
            elem = block_name + elem
            try:
                content = str(getattr(obj, attr))
            except AttributeError:
                # ensure attrs that is not present is optional
                assert elem in _DTD_OPT, "Element %r (attribute %r) not found" \
                        % (elem, attr)
            else:
                # custom element-attribute mapping, for fallback values
                if elem in opt_dict:
                    content = opt_dict[elem]
                self.xml.simpleElement(elem, content)

    def write_preamble(self, qresult):
        """Writes the XML file preamble."""
        xml = self.xml

        for elem, attr in _WRITE_MAPS['preamble']:
            elem = 'BlastOutput_' + elem
            if elem == 'BlastOutput_param':
                xml.startParent(elem)
                self.write_param(qresult)
                xml.endParent()
                continue
            try:
                content = str(getattr(qresult, attr))
            except AttributeError:
                assert elem in _DTD_OPT, "Element %s (attribute %s) not found" \
                        % (elem, attr)
            else:
                if elem == 'BlastOutput_version':
                    content = '%s %s' % (qresult.program.upper(), \
                            qresult.version)
                elif qresult._blast_id:
                    if elem == 'BlastOutput_query-ID':
                        content = qresult._blast_id
                    elif elem == 'BlastOutput_query-def':
                        content = ' '.join([qresult.id, qresult.desc]).strip()
                xml.simpleElement(elem, content)

    def write_param(self, qresult):
        """Writes the parameter block of the preamble."""
        xml = self.xml
        xml.startParent('Parameters')
        self.write_elem_block('Parameters_', 'param', qresult)
        xml.endParent()

    def write_qresults(self, qresults):
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
                    'Iteration_query-def': ' '.join([qresult.id, \
                            qresult.desc]).strip(),
                }
            self.write_elem_block('Iteration_', 'qresult', qresult, opt_dict)
            # the Iteration_hits tag only has children if there are hits
            if qresult:
                xml.startParent('Iteration_hits')
                self.write_hits(qresult.hits)
                xml.endParent()
            # otherwise it's a simple element without any contents
            else:
                xml.simpleElement('Iteration_hits', '')

            xml.startParents('Iteration_stat', 'Statistics')
            self.write_elem_block('Statistics_', 'stat', qresult)
            xml.endParents(2)
            # there's a message if no hits is present
            if not qresult:
                xml.simpleElement('Iteration_message', 'No hits found')
            self.qresult_counter += 1
            xml.endParent()

    def write_hits(self, hits):
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
                    'Hit_def': ' '.join([hit.id, hit.desc]).strip(),
                }
            self.write_elem_block('Hit_', 'hit', hit, opt_dict)
            xml.startParent('Hit_hsps')
            self.write_hsps(hit.hsps)
            self.hit_counter += 1
            xml.endParents(2)

    def write_hsps(self, hsps):
        """Writes HSP objects."""
        xml = self.xml
        for num, hsp in enumerate(hsps):
            xml.startParent('Hsp')
            xml.simpleElement('Hsp_num', str(num+1))
            for elem, attr in _WRITE_MAPS['hsp']:
                elem = 'Hsp_' + elem
                try:
                    content = self.adjust_output(hsp, elem, attr)
                # make sure any elements that is not present is optional
                # in the DTD
                except AttributeError:
                    assert elem in _DTD_OPT, "Element %s (attribute %s) not found" \
                            % (elem, attr)
                else:
                    xml.simpleElement(elem, str(content))
            self.hsp_counter += 1
            xml.endParent()

    def adjust_output(self, hsp, elem, attr):
        """Adjusts output to mimic native BLAST+ XML as much as possible."""

        # adjust coordinates
        if attr in ('query_start' ,'query_end' ,'hit_start', 'hit_end', \
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
            content = hsp.alignment_annotation['homology']
        elif elem == 'Hsp_align-len':
            content = len(hsp)
        elif elem in ('Hsp_evalue', 'Hsp_bit-score'):
            # adapted from src/algo/blast/format/blastxml_format.cpp#L138-140
            content = '%.*g' % (6, getattr(hsp, attr))
        else:
            content = getattr(hsp, attr)

        return content


def _test():
    """Run the Bio.SearchIO.BlastIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = 'Tests'

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod()
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
