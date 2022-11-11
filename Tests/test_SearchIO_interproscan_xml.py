# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO InterproscanIO parsers."""

import os
import sys
import unittest
import warnings

from Bio import BiopythonParserWarning
from Bio.SearchIO import parse


# test case files are in the Blast directory
TEST_DIR = "InterProScan"
FMT = "interproscan-xml"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class InterproscanXmlCases(unittest.TestCase):
    def test_xml_001(self):
        xml_file = get_file("test_001.xml")
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        self.assertEqual("5.26-65.0", qresult.version)

        # test parsed values of qresult
        self.assertEqual("AT5G23090.4", qresult.id)
        self.assertEqual(
            "pacid=19665592 transcript=AT5G23090.4 locus=AT5G23090 "
            "ID=AT5G23090.4.TAIR10 annot-version=TAIR10",
            qresult.description,
        )
        self.assertEqual(4, len(qresult))

        hit = qresult[0]
        self.assertEqual("PF00808", hit.id)
        self.assertEqual(
            "Histone-like transcription factor (CBF/NF-Y) and archaeal histone",
            hit.description,
        )
        self.assertEqual("PFAM", hit.attributes["Target"])
        self.assertEqual("31.0", hit.attributes["Target version"])
        self.assertEqual("hmmer3", hit.attributes["Hit type"])
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(76.7, hsp.bitscore)
        self.assertEqual(1.1e-21, hsp.evalue)
        self.assertEqual(13, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(65, hsp.hit_end)
        self.assertEqual(66, hsp.aln_span)
        self.assertEqual(
            "MDPMDIVGKSKEDASLPKATMTKIIKEMLPPDVRVARDAQDLLIECCVEFINLVSSESNDVCNKEDKRTIAPEH"
            "VLKALQVLGFGEYIEEVYAAYEQHKYETMDTQRSVKWNPGAQMTEEEAAAEQQRMFAEARARMNGGVSVPQPEH"
            "PETDQRSPQS",
            hsp.query.seq,
        )

        # parse last hit
        hit = qresult[-1]
        self.assertEqual("SSF47113", hit.id)
        self.assertEqual(1, len(hit))
        self.assertEqual("IPR:IPR009072", hit.dbxrefs[0])
        self.assertEqual("GO:0046982", hit.dbxrefs[1])

        hsp = hit.hsps[0]
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(141, hsp.query_end)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
