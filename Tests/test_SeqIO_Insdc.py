# Copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from Bio._py3k import StringIO

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from seq_tests_common import compare_record


class TestEmbl(unittest.TestCase):
    def test_annotation1(self):
        """Check parsing of annotation from EMBL files (1)."""
        record = SeqIO.read("EMBL/TRBG361.embl", "embl")
        self.assertEqual(len(record), 1859)
        # Single keyword:
        self.assertEqual(record.annotations["keywords"], ["beta-glucosidase"])
        self.assertEqual(record.annotations["topology"], "linear")

    def test_annotation2(self):
        """Check parsing of annotation from EMBL files (2)."""
        record = SeqIO.read("EMBL/DD231055_edited.embl", "embl")
        self.assertEqual(len(record), 315)
        # Multiple keywords:
        self.assertEqual(record.annotations["keywords"],
                         ['JP 2005522996-A/12', 'test-data',
                          'lot and lots of keywords for this example',
                          'multi-line keywords'])
        self.assertEqual(record.annotations["topology"], "linear")

    def test_annotation3(self):
        """Check parsing of annotation from EMBL files (3)."""
        record = SeqIO.read("EMBL/AE017046.embl", "embl")
        self.assertEqual(len(record), 9609)
        # TODO: Should this be an empty list, or simply absent?
        self.assertEqual(record.annotations["keywords"], [""])
        self.assertEqual(record.annotations["topology"], "circular")

    def test_annotation4(self):
        """Check parsing of annotation from EMBL files (4)."""
        record = SeqIO.read("EMBL/location_wrap.embl", "embl")
        self.assertEqual(len(record), 120)
        self.assertNotIn("keywords", record.annotations)
        # The ID line has the topology as unspecified:
        self.assertNotIn("topology", record.annotations)

    def test_writing_empty_qualifiers(self):
        f = SeqFeature(FeatureLocation(5, 20, strand=+1),
                       type="region",
                       qualifiers={"empty": None,
                                   "zero": 0,
                                   "one": 1,
                                   "text": "blah"})
        record = SeqRecord(Seq("A" * 100, generic_dna), "dummy",
                           features=[f])
        gbk = record.format("gb")
        self.assertIn(' /empty\n', gbk)
        self.assertIn(' /zero=0\n', gbk)
        self.assertIn(' /one=1\n', gbk)
        self.assertIn(' /text="blah"\n', gbk)


class TestEmblRewrite(unittest.TestCase):

    def check_rewrite(self, filename):
        old = SeqIO.read(filename, "embl")

        # TODO - Check these properties:
        old.dbxrefs = []
        old.annotations['accessions'] = old.annotations['accessions'][:1]
        del old.annotations['references']

        buffer = StringIO()
        self.assertEqual(1, SeqIO.write(old, buffer, "embl"))
        buffer.seek(0)
        new = SeqIO.read(buffer, "embl")

        self.assertTrue(compare_record(old, new))

    def test_annotation1(self):
        """Check writing-and-parsing EMBL file (1)."""
        self.check_rewrite("EMBL/TRBG361.embl")

    def test_annotation2(self):
        """Check writing-and-parsing EMBL file (2)."""
        self.check_rewrite("EMBL/DD231055_edited.embl")

    def test_annotation3(self):
        """Check writing-and-parsing EMBL file (3)."""
        self.check_rewrite("EMBL/AE017046.embl")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
