"""Tests for SearchIO InterproIO parsers."""

import os
import sys
import unittest
import warnings

from Bio import BiopythonParserWarning
from Bio import BiopythonExperimentalWarning


with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = 'InterProScan'
FMT = 'interpro-xml'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BlastnXmlCases(unittest.TestCase):

    def test_xml_001(self):
        xml_file = get_file('test_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = next(qresults)
        counter += 1

        self.assertEqual('2.2.12', qresult.version)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual(1, qresult.param_score_match)
        self.assertEqual(-3, qresult.param_score_mismatch)
        self.assertEqual(5, qresult.param_gap_open)
        self.assertEqual(2, qresult.param_gap_extend)

        # test parsed values of qresult
        self.assertEqual('gi|1348916|gb|G26684.1|G26684', qresult.id)
        self.assertEqual('human STS STS_D11570, sequence tagged site',
                qresult.description)
        self.assertEqual(285, qresult.seq_len)
        self.assertEqual(371021, qresult.stat_db_num)
        self.assertEqual(1233631384, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.710603, qresult.stat_kappa)
        self.assertEqual(1.37406, qresult.stat_lambda)
        self.assertEqual(1.30725, qresult.stat_entropy)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|9950606|gb|AE004854.1|', hit.id)
        self.assertEqual('Pseudomonas aeruginosa PAO1, section 415 of 529 of the complete genome', hit.description)
        self.assertEqual(11884, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit.hsps[0]
        self.assertEqual(38.1576, hsp.bitscore)
        self.assertEqual(19, hsp.bitscore_raw)
        self.assertEqual(1.0598, hsp.evalue)
        self.assertEqual(67, hsp.query_start)
        self.assertEqual(86, hsp.query_end)
        self.assertEqual(6011, hsp.hit_start)
        self.assertEqual(6030, hsp.hit_end)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(19, hsp.ident_num)
        self.assertEqual(19, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(19, hsp.aln_span)
        self.assertEqual('CAGGCCAGCGACTTCTGGG', str(hsp.query.seq))
        self.assertEqual('CAGGCCAGCGACTTCTGGG', str(hsp.hit.seq))
        self.assertEqual('|||||||||||||||||||', hsp.aln_annotation['similarity'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # parse last hit
        hit = qresult[-1]
        self.assertEqual('gi|15073988|emb|AL591786.1|SME591786', hit.id)
        self.assertEqual('Sinorhizobium meliloti 1021 complete chromosome; segment 5/12', hit.description)
        self.assertEqual(299350, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit.hsps[0]
        self.assertEqual(36.1753, hsp.bitscore)
        self.assertEqual(18, hsp.bitscore_raw)
        self.assertEqual(4.18768, hsp.evalue)
        self.assertEqual(203, hsp.query_start)
        self.assertEqual(224, hsp.query_end)
        self.assertEqual(83627, hsp.hit_start)
        self.assertEqual(83648, hsp.hit_end)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(20, hsp.ident_num)
        self.assertEqual(20, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(21, hsp.aln_span)
        self.assertEqual('TGAAAGGAAATNAAAATGGAA', str(hsp.query.seq))
        self.assertEqual('TGAAAGGAAATCAAAATGGAA', str(hsp.hit.seq))
        self.assertEqual('||||||||||| |||||||||', hsp.aln_annotation['similarity'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
