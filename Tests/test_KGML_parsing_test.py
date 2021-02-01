"""Test KGML dump to file."""

import unittest

from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.REST import kegg_get
from tempfile import TemporaryFile


class KGMLParsingTest(unittest.TestCase):
    """Download a kgml xml, parse into a pathway object, dump to xml and read it back."""

    def setUp(self):
        self.parsed_pathway = KGML_parser.read(kegg_get("ko00680", "kgml"))
        self.tmpfile_handle = TemporaryFile()

    def test_kgml_parse(self):
        self.tmpfile_handle.write(self.parsed_pathway.get_KGML().encode())
        self.tmpfile_handle.seek(0)
        KGML_parser.read(self.tmpfile_handle.read().decode())

    def tearDown(self):
        self.tmpfile_handle.close()
