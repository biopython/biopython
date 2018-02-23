"""Test."""
from Bio.Alphabet import generic_protein
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from xml.etree import cElementTree as ElementTree

# Namespace
NS = '{http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5}'

# element - hit attribute name mapping
_ELEM_HIT = {
    'name': ('id', str),
    'ac': ('accession', str),
    'desc': ('description', str),
    'library': ('target', str),
    'version': ('target_version', str),
}
# element - hsp attribute name mapping
_ELEM_HSP = {
    'score': ('bitscore', float),
    'evalue': ('evalue', float),
}
# element - fragment attribute name mapping
_ELEM_FRAG = {
    'start': ('query_start', int),
    'end': ('query_end', int),
    'hmm-start': ('hit_start', int),
    'hmm-end': ('hit_end', int),
}


class InterproXmlParser(object):
    """Parser for the InterProScan XML format."""

    def __init__(self, handle):
        """Initialize the class."""
        self.xml_iter = iter(ElementTree.iterparse(
            handle, events=('start', 'end')))
        self._meta = self._parse_header()

    def _parse_header(self):
        """Parse the header for the InterProScan version."""
        event, elem = next(self.xml_iter)
        meta = dict()
        meta['program'] = 'InterProScan'
        meta['version'] = elem.attrib['interproscan-version']
        return meta

    def _parse_qresult(self):
        """Parse query results."""
        for event, elem in self.xml_iter:
            if event == 'end' and elem.tag == NS + 'protein':
                # store the query sequence
                seq = elem.find(NS + 'sequence')
                query_seq = seq.text

                # store the query id and description
                xref = elem.find(NS + 'xref')
                query_id = xref.attrib['id']
                query_desc = xref.attrib['name']

                # parse each hit
                hit_list = [hit for hit in self._parse_hit(
                    elem.find(NS + 'matches'), query_id, query_seq)]

                # create qresult and assing attributes
                qresult = QueryResult(hit_list, query_id)
                setattr(qresult, 'description', query_desc)
                for key, value in self._meta.items():
                    setattr(qresult, key, value)
                yield qresult

    def _parse_hit(self, root_hit_elem, query_id, query_seq=None):
        """Parse hit."""
        # feed the loop below an empty list so iteration still works
        if root_hit_elem is None:
            root_hit_elem = []

        for hit_elem in root_hit_elem:
            # store the hit id
            signature = hit_elem.find(NS + 'signature')
            hit_id = signature.attrib['ac']

            # store alt_ids and alt_descs
            alt_ids, alt_descs = self._parse_alt_ids(
                signature.find(NS + 'entry'))

            # parse each hsp
            hsps = [hsp for hsp in self._parse_hsp(
                hit_elem.find(NS + 'locations'), query_id, hit_id, query_seq)]

            # create hit and assign attributes
            hit = Hit(hsps, hit_id)
            # setattr(hit, '_id_alt', alt_ids)
            for key, (attr, caster) in _ELEM_HIT.items():
                value = signature.attrib.get(key)
                if value is not None:
                    setattr(hit, attr, caster(value))
            yield hit

    def _parse_hsp(self, root_hsp_elem, query_id, hit_id, query_seq=None):
        """Parse hsp."""
        # feed the loop below an empty list so iteration still works
        if root_hsp_elem is None:
            root_hsp_elem = []

        for hsp_elem in root_hsp_elem:
            # create frag and assign attributes
            frag = HSPFragment(hit_id, query_id)
            setattr(frag, 'alphabet', generic_protein)
            if query_seq is not None:
                setattr(frag, 'query', query_seq)
            for key, (attr, caster) in _ELEM_FRAG.items():
                value = hsp_elem.attrib.get(key)
                if value is not None:
                    # start should be 0-based
                    if attr.endswith('start'):
                        value = caster(value) - 1
                    setattr(frag, attr, caster(value))

            # create hsp and assign attributes
            hsp = HSP([frag])
            setattr(hsp, 'query_id', query_id)
            setattr(hsp, 'hit_id', hit_id)
            for key, (attr, caster) in _ELEM_HSP.items():
                value = hsp_elem.attrib.get(key)
                if value is not None:
                    setattr(hsp, attr, caster(value))
            yield hsp

    def _parse_alt_ids(self, root_entry_elem):
        """Parse xrefs."""
        alt_ids, alt_descs = [], []
        # store entry id and description
        if root_entry_elem is not None:
            alt_ids.append(root_entry_elem.attrib['ac'])
            if (root_entry_elem.attrib['name'] is not None or
               root_entry_elem.attrib['desc'] is not None):
                alt_descs.append('%s %s %s' % (
                    root_entry_elem.attrib['ac'],
                    root_entry_elem.attrib.get('name'),
                    root_entry_elem.attrib.get('desc')))

        # store go-xrefs and pathway-refs id and description
        if root_entry_elem is not None:
            alt_elem = []
            alt_elem = alt_elem + root_entry_elem.findall(NS + 'go-xref')
            alt_elem = alt_elem + root_entry_elem.findall(NS + 'pathway-xref')

            for entry in alt_elem:
                if entry.attrib['id'] is not None:
                    alt_ids.append(entry.attrib['id'])
                    if (entry.attrib['name'] is not None or
                       entry.attrib['category'] is not None):
                        alt_descs.append('%s %s [%s]' % (
                            entry.attrib['id'],
                            entry.attrib.get('name'),
                            entry.attrib.get('db')))
        return alt_ids, alt_descs
