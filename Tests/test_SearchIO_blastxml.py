# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO BlastIO parsers."""


import os
import unittest

from Bio.SearchIO import parse, read

# test case files are in the Blast directory
TEST_DIR = 'Blast'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BlastXmlCases(unittest.TestCase):

    fmt = 'blast-xml'

    def test_blastn_multiple(self):
        xml_file = get_file('xbt012.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Zheng Zhang, Scott Schwartz, Lukas Wagner, and '
                'Webb Miller (2000), \"A greedy algorithm for aligning DNA '
                'sequences\", J Comput Biol 2000; 7(1-2):203-14.', \
                        qresult.meta['program_reference'])
        self.assertEqual(1.0, qresult.meta['param_score_match'])
        self.assertEqual(-2.0, qresult.meta['param_score_mismatch'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;m;', qresult.meta['param_filter'])
        self.assertEqual(0.0, qresult.meta['param_gap_open'])
        self.assertEqual(0.0, qresult.meta['param_gap_extend'])
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(7702306.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|356995852:1-490 Mus musculus POU domain, class 5, transcription factor 1 (Pou5f1), transcript variant 1, mRNA', qresult.description)
        self.assertEqual(490, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(32217922.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|0', hit.id)
        self.assertEqual('gi|356995852|ref|NM_013633.3| Mus musculus POU '
                'domain, class 5, transcription factor 1 (Pou5f1), '
                'transcript variant 1, mRNA', hit.description)
        self.assertEqual(1353, hit.seq_length)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(905.979, hsp.bitscore)
        self.assertEqual(490, hsp.bitscore_raw)
        self.assertEqual(0, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(490, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(490, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(490, hsp.identity_num)
        self.assertEqual(490, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, len(hsp))
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.query.seq.tostring())
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('hg19_dna range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 strand=+ repeatMasking=none', qresult.description)
        self.assertEqual(66, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(3545672.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|6', hit.id)
        self.assertEqual('gi|94721341|ref|NM_001040441.1| Homo sapiens zinc '
                'finger and BTB domain containing 8A (ZBTB8A), mRNA', \
                        hit.description)
        self.assertEqual(7333, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(115.613, hsp.bitscore)
        self.assertEqual(62, hsp.bitscore_raw)
        self.assertEqual(5.58273e-29, hsp.evalue)
        self.assertEqual(5, hsp.query_from)
        self.assertEqual(66, hsp.query_to)
        self.assertEqual(3677, hsp.hit_from)
        self.assertEqual(3738, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(62, hsp.identity_num)
        self.assertEqual(62, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, len(hsp))
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(98.9927, hsp.bitscore)
        self.assertEqual(53, hsp.bitscore_raw)
        self.assertEqual(5.62236e-24, hsp.evalue)
        self.assertEqual(6, hsp.query_from)
        self.assertEqual(58, hsp.query_to)
        self.assertEqual(2824, hsp.hit_from)
        self.assertEqual(2876, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(53, hsp.identity_num)
        self.assertEqual(53, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(53, len(hsp))
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('|||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|7', hit.id)
        self.assertEqual('gi|332865372|ref|XM_003318468.1| PREDICTED: Pan '
                'troglodytes zinc finger protein 273, transcript variant 1 '
                '(ZNF273), mRNA', hit.description)
        self.assertEqual(4430, hit.seq_length)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(111.919, hsp.bitscore)
        self.assertEqual(60, hsp.bitscore_raw)
        self.assertEqual(7.22171e-28, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(66, hsp.query_to)
        self.assertEqual(2800, hsp.hit_from)
        self.assertEqual(2735, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(64, hsp.identity_num)
        self.assertEqual(64, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCACGCCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_blastn_no_qresult(self):
        xml_file = get_file('xbt013.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Zheng Zhang, Scott Schwartz, Lukas Wagner, and '
                'Webb Miller (2000), \"A greedy algorithm for aligning DNA '
                'sequences\", J Comput Biol 2000; 7(1-2):203-14.', \
                        qresult.meta['program_reference'])
        self.assertEqual(1.0, qresult.meta['param_score_match'])
        self.assertEqual(-2.0, qresult.meta['param_score_mismatch'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;m;', qresult.meta['param_filter'])
        self.assertEqual(0.0, qresult.meta['param_gap_open'])
        self.assertEqual(0.0, qresult.meta['param_gap_extend'])
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(7702306.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))
        self.assertEqual([], qresult.hits)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_blastn_single_hsp(self):
        xml_file = get_file('xbt014.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Zheng Zhang, Scott Schwartz, Lukas Wagner, and '
                'Webb Miller (2000), \"A greedy algorithm for aligning DNA '
                'sequences\", J Comput Biol 2000; 7(1-2):203-14.', \
                        qresult.meta['program_reference'])
        self.assertEqual(1.0, qresult.meta['param_score_match'])
        self.assertEqual(-2.0, qresult.meta['param_score_mismatch'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;m;', qresult.meta['param_filter'])
        self.assertEqual(0.0, qresult.meta['param_gap_open'])
        self.assertEqual(0.0, qresult.meta['param_gap_extend'])
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('gi|356995852:1-490 Mus musculus POU domain, class 5, transcription '
                'factor 1 (Pou5f1), transcript variant 1, mRNA', \
                        qresult.description)
        self.assertEqual(490, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(32217922.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|0', hit.id)
        self.assertEqual('gi|356995852|ref|NM_013633.3| Mus musculus POU '
                'domain, class 5, transcription factor 1 (Pou5f1), '
                'transcript variant 1, mRNA', hit.description)
        self.assertEqual(1353, hit.seq_length)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(905.979, hsp.bitscore)
        self.assertEqual(490, hsp.bitscore_raw)
        self.assertEqual(0, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(490, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(490, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(490, hsp.identity_num)
        self.assertEqual(490, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, len(hsp))
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.query.seq.tostring())
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)


        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_blastn_multiple_hsps(self):
        xml_file = get_file('xbt015.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('hg19_dna range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none',
                        qresult.description)
        self.assertEqual(66, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(3545672.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|6', hit.id)
        self.assertEqual('gi|94721341|ref|NM_001040441.1| Homo sapiens zinc '
                'finger and BTB domain containing 8A (ZBTB8A), mRNA', \
                        hit.description)
        self.assertEqual(7333, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(115.613, hsp.bitscore)
        self.assertEqual(62, hsp.bitscore_raw)
        self.assertEqual(5.58273e-29, hsp.evalue)
        self.assertEqual(5, hsp.query_from)
        self.assertEqual(66, hsp.query_to)
        self.assertEqual(3677, hsp.hit_from)
        self.assertEqual(3738, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(62, hsp.identity_num)
        self.assertEqual(62, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, len(hsp))
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(98.9927, hsp.bitscore)
        self.assertEqual(53, hsp.bitscore_raw)
        self.assertEqual(5.62236e-24, hsp.evalue)
        self.assertEqual(6, hsp.query_from)
        self.assertEqual(58, hsp.query_to)
        self.assertEqual(2824, hsp.hit_from)
        self.assertEqual(2876, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(53, hsp.identity_num)
        self.assertEqual(53, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(53, len(hsp))
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('|||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|7', hit.id)
        self.assertEqual('gi|332865372|ref|XM_003318468.1| PREDICTED: Pan '
                'troglodytes zinc finger protein 273, transcript variant 1 '
                '(ZNF273), mRNA', hit.description)
        self.assertEqual(4430, hit.seq_length)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(111.919, hsp.bitscore)
        self.assertEqual(60, hsp.bitscore_raw)
        self.assertEqual(7.22171e-28, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(66, hsp.query_to)
        self.assertEqual(2800, hsp.hit_from)
        self.assertEqual(2735, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(64, hsp.identity_num)
        self.assertEqual(64, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCACGCCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_blastn_remote(self):
        xml_file = get_file('xbt016.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Zheng Zhang, Scott Schwartz, Lukas Wagner, and '
                'Webb Miller (2000), \"A greedy algorithm for aligning DNA '
                'sequences\", J Comput Biol 2000; 7(1-2):203-14.', \
                        qresult.meta['program_reference'])
        self.assertEqual(1.0, qresult.meta['param_score_match'])
        self.assertEqual(-2.0, qresult.meta['param_score_mismatch'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;m;', qresult.meta['param_filter'])
        self.assertEqual(0.0, qresult.meta['param_gap_open'])
        self.assertEqual(0.0, qresult.meta['param_gap_extend'])
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('refseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|356995852:1-490 Mus musculus POU domain, class 5, transcription '
                'factor 1 (Pou5f1), transcript variant 1, mRNA', \
                        qresult.description)
        self.assertEqual(490, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|356995852|ref|NM_013633.3|', hit.id)
        self.assertEqual('Mus musculus POU '
                'domain, class 5, transcription factor 1 (Pou5f1), '
                'transcript variant 1, mRNA', hit.description)
        self.assertEqual(1353, hit.seq_length)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(905.979, hsp.bitscore)
        self.assertEqual(490, hsp.bitscore_raw)
        self.assertEqual(0, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(490, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(490, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(490, hsp.identity_num)
        self.assertEqual(490, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, len(hsp))
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.query.seq.tostring())
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('hg19_dna range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none',
                        qresult.description)
        self.assertEqual(66, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|332237160|ref|XM_003267724.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys ATG14 autophagy '
                'related 14 homolog (S. cerevisiae) (ATG14), mRNA', \
                        hit.description)
        self.assertEqual(4771, hit.seq_length)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(115.613, hsp.bitscore)
        self.assertEqual(62, hsp.bitscore_raw)
        self.assertEqual(3.31905e-23, hsp.evalue)
        self.assertEqual(5, hsp.query_from)
        self.assertEqual(66, hsp.query_to)
        self.assertEqual(2865, hsp.hit_from)
        self.assertEqual(2926, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(62, hsp.identity_num)
        self.assertEqual(62, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, len(hsp))
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|332254616|ref|XM_003276378.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys S100P binding '
                'protein, transcript variant 2 (S100PBP), mRNA', \
                        hit.description)
        self.assertEqual(4345, hit.seq_length)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(111.919, hsp.bitscore)
        self.assertEqual(60, hsp.bitscore_raw)
        self.assertEqual(4.29347e-22, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(66, hsp.query_to)
        self.assertEqual(2857, hsp.hit_from)
        self.assertEqual(2792, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(64, hsp.identity_num)
        self.assertEqual(64, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCATGCCACTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_blastp_multiple(self):
        xml_file = get_file('xbt017.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('F', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('db/minirefseq_protein', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(32, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(156650.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis '
                'subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(361344, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|1', hit.id)
        self.assertEqual('gi|308175296|ref|YP_003922001.1| membrane bound '
                'lipoprotein [Bacillus amyloliquefaciens DSM 7]', \
                hit.description)
        self.assertEqual(100, hit.seq_length)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(139.428, hsp.bitscore)
        self.assertEqual(350, hsp.bitscore_raw)
        self.assertEqual(1.99275e-46, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(102, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.identity_num)
        self.assertEqual(81, hsp.positive_num)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(102, len(hsp))
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.query.seq.tostring())
        self.assertEqual('MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPVSFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN', hsp.hit.seq.tostring())
        self.assertEqual('MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(345626, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|5', hit.id)
        self.assertEqual('gi|11464971|ref|NP_062422.1| pleckstrin '
                '[Mus musculus]', hit.description)
        self.assertEqual(350, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(205.682, hsp.bitscore)
        self.assertEqual(522, hsp.bitscore_raw)
        self.assertEqual(2.24956e-69, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(4, hsp.hit_from)
        self.assertEqual(101, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(98, hsp.identity_num)
        self.assertEqual(98, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(43.5134, hsp.bitscore)
        self.assertEqual(101, hsp.bitscore_raw)
        self.assertEqual(2.90061e-09, hsp.evalue)
        self.assertEqual(3, hsp.query_from)
        self.assertEqual(96, hsp.query_to)
        self.assertEqual(246, hsp.hit_from)
        self.assertEqual(345, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(29, hsp.identity_num)
        self.assertEqual(48, hsp.positive_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|9', hit.id)
        self.assertEqual('gi|350596020|ref|XP_003360649.2| PREDICTED: '
                'pleckstrin-like [Sus scrofa]', hit.description)
        self.assertEqual(228, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.97058e-68, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(4, hsp.hit_from)
        self.assertEqual(101, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(94, hsp.identity_num)
        self.assertEqual(96, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_blastp_no_qresult(self):
        xml_file = get_file('xbt018.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('F', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('db/minirefseq_protein', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(32, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(156650.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_blastp_single_hsp(self):
        xml_file = get_file('xbt019.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('F', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('db/minirefseq_protein', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis '
                'subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(361344, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|1', hit.id)
        self.assertEqual('gi|308175296|ref|YP_003922001.1| membrane bound '
                'lipoprotein [Bacillus amyloliquefaciens DSM 7]', \
                hit.description)
        self.assertEqual(100, hit.seq_length)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(139.428, hsp.bitscore)
        self.assertEqual(350, hsp.bitscore_raw)
        self.assertEqual(1.99275e-46, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(102, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.identity_num)
        self.assertEqual(81, hsp.positive_num)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(102, len(hsp))
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.query.seq.tostring())
        self.assertEqual('MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPVSFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN', hsp.hit.seq.tostring())
        self.assertEqual('MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_blastp_multiple_hsps(self):
        xml_file = get_file('xbt020.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('F', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('db/minirefseq_protein', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(345626, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|5', hit.id)
        self.assertEqual('gi|11464971|ref|NP_062422.1| pleckstrin '
                '[Mus musculus]', hit.description)
        self.assertEqual(350, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(205.682, hsp.bitscore)
        self.assertEqual(522, hsp.bitscore_raw)
        self.assertEqual(2.24956e-69, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(4, hsp.hit_from)
        self.assertEqual(101, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(98, hsp.identity_num)
        self.assertEqual(98, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(43.5134, hsp.bitscore)
        self.assertEqual(101, hsp.bitscore_raw)
        self.assertEqual(2.90061e-09, hsp.evalue)
        self.assertEqual(3, hsp.query_from)
        self.assertEqual(96, hsp.query_to)
        self.assertEqual(246, hsp.hit_from)
        self.assertEqual(345, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(29, hsp.identity_num)
        self.assertEqual(48, hsp.positive_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|9', hit.id)
        self.assertEqual('gi|350596020|ref|XP_003360649.2| PREDICTED: '
                'pleckstrin-like [Sus scrofa]', hit.description)
        self.assertEqual(228, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.97058e-68, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(4, hsp.hit_from)
        self.assertEqual(101, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(94, hsp.identity_num)
        self.assertEqual(96, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_blastp_remote(self):
        xml_file = get_file('xbt021.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('F', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('refseq_protein', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(32, qresult.query_length)
        self.assertEqual(12504708, qresult.stat_db_num)
        self.assertEqual(4346765131, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis '
                'subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.query_length)
        self.assertEqual(12504708, qresult.stat_db_num)
        self.assertEqual(4346765131, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|16080617|ref|NP_391444.1|', hit.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168] >gi|221311516|ref|ZP_03593363.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168] >gi|221315843|ref|ZP_03597648.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. NCIB 3610] >gi|221320757|ref|ZP_03602051.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. JH642] >gi|221325043|ref|ZP_03606337.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. SMY] >gi|321313111|ref|YP_004205398.1| unnamed protein product [Bacillus subtilis BSn5]', hit.description)
        self.assertEqual(102, hit.seq_length)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(205.297, hsp.bitscore)
        self.assertEqual(521, hsp.bitscore_raw)
        self.assertEqual(1.43621e-66, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(102, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(102, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(102, hsp.identity_num)
        self.assertEqual(102, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(102, len(hsp))
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.query.seq.tostring())
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.hit.seq.tostring())
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.query_length)
        self.assertEqual(12504708, qresult.stat_db_num)
        self.assertEqual(4346765131, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|11464971|ref|NP_062422.1|', hit.id)
        self.assertEqual('pleckstrin [Mus musculus]', hit.description)
        self.assertEqual(350, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(205.682, hsp.bitscore)
        self.assertEqual(522, hsp.bitscore_raw)
        self.assertEqual(1.52643e-63, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(4, hsp.hit_from)
        self.assertEqual(101, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(98, hsp.identity_num)
        self.assertEqual(98, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(43.5134, hsp.bitscore)
        self.assertEqual(101, hsp.bitscore_raw)
        self.assertEqual(0.0019682, hsp.evalue)
        self.assertEqual(3, hsp.query_from)
        self.assertEqual(96, hsp.query_to)
        self.assertEqual(246, hsp.hit_from)
        self.assertEqual(345, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(29, hsp.identity_num)
        self.assertEqual(48, hsp.positive_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|350596020|ref|XP_003360649.2|', hit.id)
        self.assertEqual('PREDICTED: pleckstrin-like [Sus scrofa]', \
                hit.description)
        self.assertEqual(228, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.33713e-62, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(4, hsp.hit_from)
        self.assertEqual(101, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(94, hsp.identity_num)
        self.assertEqual(96, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_blastx_multiple(self):
        xml_file = get_file('xbt022.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), "Gapped BLAST '
                'and PSI-BLAST: a new generation of protein database search '
                'programs", Nucleic Acids Res. 25:3389-3402.',
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('db/minirefseq_protein', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|356995852:1-490 Mus musculus POU domain, class 5, transcription '
                'factor 1 (Pou5f1), transcript variant 1, mRNA', \
                        qresult.description)
        self.assertEqual(490, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(673486.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|10', hit.id)
        self.assertEqual('gi|125490392|ref|NP_038661.2| POU domain, class 5'
                ', transcription factor 1 isoform 1 [Mus musculus]', \
                        hit.description)
        self.assertEqual(352, hit.seq_length)
        self.assertEqual(3, len(hit))

        hsp = hit[0]
        self.assertEqual(192.2, hsp.bitscore)
        self.assertEqual(487, hsp.bitscore_raw)
        self.assertEqual(5.46624e-63, hsp.evalue)
        self.assertEqual(69, hsp.query_from)
        self.assertEqual(488, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(140, hsp.hit_to)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(140, hsp.identity_num)
        self.assertEqual(140, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(140, len(hsp))
        self.assertEqual('MAGHLAXXXXXXXXXXXXXXXXXLEPGWVDPRTWLSFQXXXXXXXXXXXSEVLGISPCPPAYEFCGGMAYCXXXXXXXXXXXXXXETLQPEGQAGARVESNSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', hsp.query.seq.tostring())
        self.assertEqual('MAGHLASDFAFSPPPGGGDGSAGLEPGWVDPRTWLSFQGPPGGPGIGPGSEVLGISPCPPAYEFCGGMAYCGPQVGLGLVPQVGVETLQPEGQAGARVESNSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', hsp.hit.seq.tostring())
        self.assertEqual('MAGHLASDFAFSPPPGGGDGSAGLEPGWVDPRTWLSFQGPPGGPGIGPGSEVLGISPCPPAYEFCGGMAYCGPQVGLGLVPQVGVETLQPEGQAGARVESNSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', hsp.homology)
        #self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('hg19_dna range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none',
                        qresult.description)
        self.assertEqual(485, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(662354, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|15', hit.id)
        self.assertEqual('gi|332258565|ref|XP_003278367.1| PREDICTED: '
                'UPF0764 protein C16orf89-like [Nomascus leucogenys]',
                        hit.description)
        self.assertEqual(132, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(121.709, hsp.bitscore)
        self.assertEqual(304, hsp.bitscore_raw)
        self.assertEqual(2.9522e-38, hsp.evalue)
        self.assertEqual(16, hsp.query_from)
        self.assertEqual(300, hsp.query_to)
        self.assertEqual(25, hsp.hit_from)
        self.assertEqual(119, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.identity_num)
        self.assertEqual(74, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(95, len(hsp))
        self.assertEqual('LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', hsp.query.seq.tostring())
        self.assertEqual('LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP', hsp.hit.seq.tostring())
        self.assertEqual('LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(51.6026, hsp.bitscore)
        self.assertEqual(122, hsp.bitscore_raw)
        self.assertEqual(2.73605e-12, hsp.evalue)
        self.assertEqual(244, hsp.query_from)
        self.assertEqual(459, hsp.query_to)
        self.assertEqual(32, hsp.hit_from)
        self.assertEqual(98, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(34, hsp.identity_num)
        self.assertEqual(41, hsp.positive_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|19', hit.id)
        self.assertEqual('gi|33188429|ref|NP_872601.1| histone demethylase '
                'UTY isoform 1 [Homo sapiens]', hit.description)
        self.assertEqual(1079, hit.seq_length)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(104.375, hsp.bitscore)
        self.assertEqual(259, hsp.bitscore_raw)
        self.assertEqual(6.31914e-29, hsp.evalue)
        self.assertEqual(19, hsp.query_from)
        self.assertEqual(291, hsp.query_to)
        self.assertEqual(989, hsp.hit_from)
        self.assertEqual(1079, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(59, hsp.identity_num)
        self.assertEqual(66, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, len(hsp))
        self.assertEqual('SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ', hsp.query.seq.tostring())
        self.assertEqual('SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQACLELLTSGGLLASASQSAGITGVSHHAR', hsp.hit.seq.tostring())
        self.assertEqual('SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_blastx_no_qresult(self):
        xml_file = get_file('xbt023.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), "Gapped BLAST '
                'and PSI-BLAST: a new generation of protein database search '
                'programs", Nucleic Acids Res. 25:3389-3402.',
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('db/minirefseq_protein', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_blastx_multiple_hsps(self):
        xml_file = get_file('xbt024.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), "Gapped BLAST '
                'and PSI-BLAST: a new generation of protein database search '
                'programs", Nucleic Acids Res. 25:3389-3402.',
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('db/minirefseq_protein', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('hg19_dna range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none',
                        qresult.description)
        self.assertEqual(485, qresult.query_length)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(662354, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|15', hit.id)
        self.assertEqual('gi|332258565|ref|XP_003278367.1| PREDICTED: '
                'UPF0764 protein C16orf89-like [Nomascus leucogenys]',
                        hit.description)
        self.assertEqual(132, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(121.709, hsp.bitscore)
        self.assertEqual(304, hsp.bitscore_raw)
        self.assertEqual(2.9522e-38, hsp.evalue)
        self.assertEqual(16, hsp.query_from)
        self.assertEqual(300, hsp.query_to)
        self.assertEqual(25, hsp.hit_from)
        self.assertEqual(119, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.identity_num)
        self.assertEqual(74, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(95, len(hsp))
        self.assertEqual('LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', hsp.query.seq.tostring())
        self.assertEqual('LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP', hsp.hit.seq.tostring())
        self.assertEqual('LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(51.6026, hsp.bitscore)
        self.assertEqual(122, hsp.bitscore_raw)
        self.assertEqual(2.73605e-12, hsp.evalue)
        self.assertEqual(244, hsp.query_from)
        self.assertEqual(459, hsp.query_to)
        self.assertEqual(32, hsp.hit_from)
        self.assertEqual(98, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(34, hsp.identity_num)
        self.assertEqual(41, hsp.positive_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|19', hit.id)
        self.assertEqual('gi|33188429|ref|NP_872601.1| histone demethylase '
                'UTY isoform 1 [Homo sapiens]', hit.description)
        self.assertEqual(1079, hit.seq_length)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(104.375, hsp.bitscore)
        self.assertEqual(259, hsp.bitscore_raw)
        self.assertEqual(6.31914e-29, hsp.evalue)
        self.assertEqual(19, hsp.query_from)
        self.assertEqual(291, hsp.query_to)
        self.assertEqual(989, hsp.hit_from)
        self.assertEqual(1079, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(59, hsp.identity_num)
        self.assertEqual(66, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, len(hsp))
        self.assertEqual('SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ', hsp.query.seq.tostring())
        self.assertEqual('SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQACLELLTSGGLLASASQSAGITGVSHHAR', hsp.hit.seq.tostring())
        self.assertEqual('SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_blastx_remote(self):
        xml_file = get_file('xbt025.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), "Gapped BLAST '
                'and PSI-BLAST: a new generation of protein database search '
                'programs", Nucleic Acids Res. 25:3389-3402.',
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('refseq_protein', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(12504708, qresult.stat_db_num)
        self.assertEqual(4346765131, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|356995852:1-490 Mus musculus POU domain, class 5, transcription '
                'factor 1 (Pou5f1), transcript variant 1, mRNA', \
                        qresult.description)
        self.assertEqual(490, qresult.query_length)
        self.assertEqual(12504708, qresult.stat_db_num)
        self.assertEqual(4346765131, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hit.id)
        self.assertEqual('POU domain, class 5, transcription factor 1 '
                'isoform 1 [Mus musculus]', hit.description)
        self.assertEqual(352, hit.seq_length)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(192.2, hsp.bitscore)
        self.assertEqual(487, hsp.bitscore_raw)
        self.assertEqual(3.7091e-57, hsp.evalue)
        self.assertEqual(69, hsp.query_from)
        self.assertEqual(488, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(140, hsp.hit_to)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(140, hsp.identity_num)
        self.assertEqual(140, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(140, len(hsp))
        self.assertEqual('MAGHLAXXXXXXXXXXXXXXXXXLEPGWVDPRTWLSFQXXXXXXXXXXXSEVLGISPCPPAYEFCGGMAYCXXXXXXXXXXXXXXETLQPEGQAGARVESNSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', hsp.query.seq.tostring())
        self.assertEqual('MAGHLASDFAFSPPPGGGDGSAGLEPGWVDPRTWLSFQGPPGGPGIGPGSEVLGISPCPPAYEFCGGMAYCGPQVGLGLVPQVGVETLQPEGQAGARVESNSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', hsp.hit.seq.tostring())
        self.assertEqual('MAGHLASDFAFSPPPGGGDGSAGLEPGWVDPRTWLSFQGPPGGPGIGPGSEVLGISPCPPAYEFCGGMAYCGPQVGLGLVPQVGVETLQPEGQAGARVESNSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('hg19_dna range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none',
                        qresult.description)
        self.assertEqual(485, qresult.query_length)
        self.assertEqual(12504708, qresult.stat_db_num)
        self.assertEqual(4346765131, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|332258565|ref|XP_003278367.1|', hit.id)
        self.assertEqual('PREDICTED: UPF0764 protein C16orf89-like [Nomascus '
                'leucogenys]', hit.description)
        self.assertEqual(132, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(121.709, hsp.bitscore)
        self.assertEqual(304, hsp.bitscore_raw)
        self.assertEqual(2.0032e-32, hsp.evalue)
        self.assertEqual(16, hsp.query_from)
        self.assertEqual(300, hsp.query_to)
        self.assertEqual(25, hsp.hit_from)
        self.assertEqual(119, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.identity_num)
        self.assertEqual(74, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(95, len(hsp))
        self.assertEqual('LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', hsp.query.seq.tostring())
        self.assertEqual('LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP', hsp.hit.seq.tostring())
        self.assertEqual('LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(51.6026, hsp.bitscore)
        self.assertEqual(122, hsp.bitscore_raw)
        self.assertEqual(1.85654e-06, hsp.evalue)
        self.assertEqual(244, hsp.query_from)
        self.assertEqual(459, hsp.query_to)
        self.assertEqual(32, hsp.hit_from)
        self.assertEqual(98, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(34, hsp.identity_num)
        self.assertEqual(41, hsp.positive_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|57113895|ref|NP_001009002.1|', hit.id)
        self.assertEqual('histone demethylase UTY [Pan troglodytes]', \
                hit.description)
        self.assertEqual(1079, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(104.76, hsp.bitscore)
        self.assertEqual(260, hsp.bitscore_raw)
        self.assertEqual(3.15863e-23, hsp.evalue)
        self.assertEqual(19, hsp.query_from)
        self.assertEqual(291, hsp.query_to)
        self.assertEqual(989, hsp.hit_from)
        self.assertEqual(1079, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(59, hsp.identity_num)
        self.assertEqual(66, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, len(hsp))
        self.assertEqual('SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ', hsp.query.seq.tostring())
        self.assertEqual('SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQAHLELLTSGGLLASASQSAGITGVSHHAR', hsp.hit.seq.tostring())
        self.assertEqual('SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_tblastn_multiple(self):
        xml_file = get_file('xbt026.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(32, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis '
                'subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(1217216.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(3, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|10', hit.id)
        self.assertEqual('gi|145479850|ref|XM_001425911.1| Paramecium '
                'tetraurelia hypothetical protein (GSPATT00004923001) '
                'partial mRNA', hit.description)
        self.assertEqual(4632, hit.seq_length)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(34.6538, hsp.bitscore)
        self.assertEqual(78, hsp.bitscore_raw)
        self.assertEqual(1.09472e-05, hsp.evalue)
        self.assertEqual(31, hsp.query_from)
        self.assertEqual(73, hsp.query_to)
        self.assertEqual(1744, hsp.hit_from)
        self.assertEqual(1872, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(15, hsp.identity_num)
        self.assertEqual(26, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(43, len(hsp))
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', hsp.query.seq.tostring())
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', hsp.hit.seq.tostring())
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(1130272.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|14', hit.id)
        self.assertEqual('gi|350596019|ref|XM_003360601.2| PREDICTED: Sus '
                'scrofa pleckstrin-like (LOC100626968), mRNA', hit.description)
        self.assertEqual(772, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.59038e-67, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(95, hsp.hit_from)
        self.assertEqual(388, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(94, hsp.identity_num)
        self.assertEqual(96, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(32.7278, hsp.bitscore)
        self.assertEqual(73, hsp.bitscore_raw)
        self.assertEqual(4.12155e-05, hsp.evalue)
        self.assertEqual(30, hsp.query_from)
        self.assertEqual(96, hsp.query_to)
        self.assertEqual(542, hsp.hit_from)
        self.assertEqual(754, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(21, hsp.identity_num)
        self.assertEqual(33, hsp.positive_num)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(71, len(hsp))
        self.assertEqual('IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('+ +Y       P G I L+G  +TS      GK  F+ +     T  +  +F QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|17', hit.id)
        self.assertEqual('gi|338714227|ref|XM_001492113.3| PREDICTED: Equus '
                'caballus pleckstrin-like (LOC100051039), mRNA', hit.description)
        self.assertEqual(1390, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(202.986, hsp.bitscore)
        self.assertEqual(515, hsp.bitscore_raw)
        self.assertEqual(1.5394e-66, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(173, hsp.hit_from)
        self.assertEqual(466, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(96, hsp.identity_num)
        self.assertEqual(97, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVK+GSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_tblastn_no_qresult(self):
        xml_file = get_file('xbt027.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(32, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_tblastn_single_hsp(self):
        xml_file = get_file('xbt028.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis '
                'subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(1217216.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(3, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|10', hit.id)
        self.assertEqual('gi|145479850|ref|XM_001425911.1| Paramecium '
                'tetraurelia hypothetical protein (GSPATT00004923001) '
                'partial mRNA', hit.description)
        self.assertEqual(4632, hit.seq_length)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(34.6538, hsp.bitscore)
        self.assertEqual(78, hsp.bitscore_raw)
        self.assertEqual(1.09472e-05, hsp.evalue)
        self.assertEqual(31, hsp.query_from)
        self.assertEqual(73, hsp.query_to)
        self.assertEqual(1744, hsp.hit_from)
        self.assertEqual(1872, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(15, hsp.identity_num)
        self.assertEqual(26, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(43, len(hsp))
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', hsp.query.seq.tostring())
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', hsp.hit.seq.tostring())
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_tblastn_multiple_hsps(self):
        xml_file = get_file('xbt029.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(1130272.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|14', hit.id)
        self.assertEqual('gi|350596019|ref|XM_003360601.2| PREDICTED: Sus '
                'scrofa pleckstrin-like (LOC100626968), mRNA', hit.description)
        self.assertEqual(772, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.59038e-67, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(95, hsp.hit_from)
        self.assertEqual(388, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(94, hsp.identity_num)
        self.assertEqual(96, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(32.7278, hsp.bitscore)
        self.assertEqual(73, hsp.bitscore_raw)
        self.assertEqual(4.12155e-05, hsp.evalue)
        self.assertEqual(30, hsp.query_from)
        self.assertEqual(96, hsp.query_to)
        self.assertEqual(542, hsp.hit_from)
        self.assertEqual(754, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(21, hsp.identity_num)
        self.assertEqual(33, hsp.positive_num)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(71, len(hsp))
        self.assertEqual('IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('+ +Y       P G I L+G  +TS      GK  F+ +     T  +  +F QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|17', hit.id)
        self.assertEqual('gi|338714227|ref|XM_001492113.3| PREDICTED: Equus '
                'caballus pleckstrin-like (LOC100051039), mRNA', hit.description)
        self.assertEqual(1390, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(202.986, hsp.bitscore)
        self.assertEqual(515, hsp.bitscore_raw)
        self.assertEqual(1.5394e-66, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(173, hsp.hit_from)
        self.assertEqual(466, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(96, hsp.identity_num)
        self.assertEqual(97, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVK+GSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_tblastn_remote(self):
        xml_file = get_file('xbt030.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('refseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(32, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis '
                'subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(3, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|145479850|ref|XM_001425911.1|', hit.id)
        self.assertEqual('Paramecium tetraurelia hypothetical protein '
                '(GSPATT00004923001) partial mRNA', hit.description)
        self.assertEqual(4632, hit.seq_length)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(34.6538, hsp.bitscore)
        self.assertEqual(78, hsp.bitscore_raw)
        self.assertEqual(0.746153, hsp.evalue)
        self.assertEqual(31, hsp.query_from)
        self.assertEqual(73, hsp.query_to)
        self.assertEqual(1744, hsp.hit_from)
        self.assertEqual(1872, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(15, hsp.identity_num)
        self.assertEqual(26, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(43, len(hsp))
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', hsp.query.seq.tostring())
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', hsp.hit.seq.tostring())
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|354480463|ref|XM_003502378.1|', hit.id)
        self.assertEqual('PREDICTED: Cricetulus griseus pleckstrin-like '
                '(LOC100773128), mRNA', hit.description)
        self.assertEqual(1119, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(205.297, hsp.bitscore)
        self.assertEqual(521, hsp.bitscore_raw)
        self.assertEqual(1.44426e-63, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(76, hsp.hit_from)
        self.assertEqual(369, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(98, hsp.identity_num)
        self.assertEqual(98, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(43.8986, hsp.bitscore)
        self.assertEqual(102, hsp.bitscore_raw)
        self.assertEqual(0.000535139, hsp.evalue)
        self.assertEqual(3, hsp.query_from)
        self.assertEqual(96, hsp.query_to)
        self.assertEqual(802, hsp.hit_from)
        self.assertEqual(1101, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(30, hsp.identity_num)
        self.assertEqual(50, hsp.positive_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDF-GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAIHLRGCVVTSVESNHDGKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +  GK+     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|338714227|ref|XM_001492113.3|', hit.id)
        self.assertEqual('PREDICTED: Equus caballus pleckstrin-like '
                '(LOC100051039), mRNA', hit.description)
        self.assertEqual(1390, hit.seq_length)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(202.986, hsp.bitscore)
        self.assertEqual(515, hsp.bitscore_raw)
        self.assertEqual(1.04924e-61, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(98, hsp.query_to)
        self.assertEqual(173, hsp.hit_from)
        self.assertEqual(466, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(96, hsp.identity_num)
        self.assertEqual(97, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVK+GSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_tblastx_multiple(self):
        xml_file = get_file('xbt031.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|347972582:1-35 Anopheles gambiae str. PEST AGAP011294-PA '
                '(DEFI_ANOGA) mRNA, complete cds', qresult.description)
        self.assertEqual(35, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(182144.0, qresult.stat_eff_space)
        self.assertEqual(0.133956144488482, qresult.stat_kappa)
        self.assertEqual(0.317605957635731, qresult.stat_lambda)
        self.assertEqual(0.401214524497119, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|18', hit.id)
        self.assertEqual('gi|347972582|ref|XM_309352.4| Anopheles gambiae '
                'str. PEST AGAP011294-PA (DEFI_ANOGA) mRNA, complete cds', \
                hit.description)
        self.assertEqual(309, hit.seq_length)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(34.0583, hsp.bitscore)
        self.assertEqual(68, hsp.bitscore_raw)
        self.assertEqual(2.03639e-05, hsp.evalue)
        self.assertEqual(3, hsp.query_from)
        self.assertEqual(35, hsp.query_to)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(35, hsp.hit_to)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(3, hsp.hit_frame)
        self.assertEqual(11, hsp.identity_num)
        self.assertEqual(11, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, len(hsp))
        self.assertEqual('EVCNDRLYHCR', hsp.query.seq.tostring())
        self.assertEqual('EVCNDRLYHCR', hsp.hit.seq.tostring())
        self.assertEqual('EVCNDRLYHCR', hsp.homology)
        #self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, '
                'complete cds', qresult.description)
        self.assertEqual(350, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(1975088.0, qresult.stat_eff_space)
        self.assertEqual(0.133956144488482, qresult.stat_kappa)
        self.assertEqual(0.317605957635731, qresult.stat_lambda)
        self.assertEqual(0.401214524497119, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|19', hit.id)
        self.assertEqual('gi|296147483|ref|NM_001183135.1| Saccharomyces '
                'cerevisiae S288c Mon2p (MON2) mRNA, complete cds', hit.description)
        self.assertEqual(4911, hit.seq_length)
        self.assertEqual(8, len(hit))

        hsp = hit[0]
        self.assertEqual(289.739, hsp.bitscore)
        self.assertEqual(626, hsp.bitscore_raw)
        self.assertEqual(2.37998e-81, hsp.evalue)
        self.assertEqual(2, hsp.query_from)
        self.assertEqual(349, hsp.query_to)
        self.assertEqual(2, hsp.hit_from)
        self.assertEqual(349, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(116, hsp.identity_num)
        self.assertEqual(116, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(116, len(hsp))
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.hit.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.query.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(18.9375, hsp.bitscore)
        self.assertEqual(35, hsp.bitscore_raw)
        self.assertEqual(7.86812, hsp.evalue)
        self.assertEqual(293, hsp.query_from)
        self.assertEqual(325, hsp.query_to)
        self.assertEqual(341, hsp.hit_from)
        self.assertEqual(373, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(-3, hsp.hit_frame)
        self.assertEqual(6, hsp.identity_num)
        self.assertEqual(9, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, len(hsp))
        self.assertEqual('KFWMPSLRLLI', hsp.query.seq.tostring())
        self.assertEqual('KKWVPPVKLLI', hsp.hit.seq.tostring())
        self.assertEqual('K W+P ++LLI', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|23', hit.id)
        self.assertEqual('gi|254579534|ref|XM_002495708.1| Zygosaccharomyces '
                'rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds', \
                hit.description)
        self.assertEqual(4866, hit.seq_length)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(141.279, hsp.bitscore)
        self.assertEqual(302, hsp.bitscore_raw)
        self.assertEqual(1.16777e-36, hsp.evalue)
        self.assertEqual(97, hsp.query_from)
        self.assertEqual(348, hsp.query_to)
        self.assertEqual(97, hsp.hit_from)
        self.assertEqual(348, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(57, hsp.identity_num)
        self.assertEqual(72, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(84, len(hsp))
        self.assertEqual('IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', hsp.query.seq.tostring())
        self.assertEqual('IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK', hsp.hit.seq.tostring())
        self.assertEqual('IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_tblastx_no_qresult(self):
        xml_file = get_file('xbt032.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_tblastx_multiple_hsps(self):
        xml_file = get_file('xbt033.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('db/minirefseq_rna', qresult.target)

        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, '
                'complete cds', qresult.description)
        self.assertEqual(350, qresult.query_length)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(1975088.0, qresult.stat_eff_space)
        self.assertEqual(0.133956144488482, qresult.stat_kappa)
        self.assertEqual(0.317605957635731, qresult.stat_lambda)
        self.assertEqual(0.401214524497119, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gnl|BL_ORD_ID|19', hit.id)
        self.assertEqual('gi|296147483|ref|NM_001183135.1| Saccharomyces '
                'cerevisiae S288c Mon2p (MON2) mRNA, complete cds', hit.description)
        self.assertEqual(4911, hit.seq_length)
        self.assertEqual(8, len(hit))

        hsp = hit[0]
        self.assertEqual(289.739, hsp.bitscore)
        self.assertEqual(626, hsp.bitscore_raw)
        self.assertEqual(2.37998e-81, hsp.evalue)
        self.assertEqual(2, hsp.query_from)
        self.assertEqual(349, hsp.query_to)
        self.assertEqual(2, hsp.hit_from)
        self.assertEqual(349, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(116, hsp.identity_num)
        self.assertEqual(116, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(116, len(hsp))
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.hit.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.query.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(18.9375, hsp.bitscore)
        self.assertEqual(35, hsp.bitscore_raw)
        self.assertEqual(7.86812, hsp.evalue)
        self.assertEqual(293, hsp.query_from)
        self.assertEqual(325, hsp.query_to)
        self.assertEqual(341, hsp.hit_from)
        self.assertEqual(373, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(-3, hsp.hit_frame)
        self.assertEqual(6, hsp.identity_num)
        self.assertEqual(9, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, len(hsp))
        self.assertEqual('KFWMPSLRLLI', hsp.query.seq.tostring())
        self.assertEqual('KKWVPPVKLLI', hsp.hit.seq.tostring())
        self.assertEqual('K W+P ++LLI', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|23', hit.id)
        self.assertEqual('gi|254579534|ref|XM_002495708.1| Zygosaccharomyces '
                'rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds', \
                hit.description)
        self.assertEqual(4866, hit.seq_length)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(141.279, hsp.bitscore)
        self.assertEqual(302, hsp.bitscore_raw)
        self.assertEqual(1.16777e-36, hsp.evalue)
        self.assertEqual(97, hsp.query_from)
        self.assertEqual(348, hsp.query_to)
        self.assertEqual(97, hsp.hit_from)
        self.assertEqual(348, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(57, hsp.identity_num)
        self.assertEqual(72, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(84, len(hsp))
        self.assertEqual('IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', hsp.query.seq.tostring())
        self.assertEqual('IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK', hsp.hit.seq.tostring())
        self.assertEqual('IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_tblastx_remote(self):
        xml_file = get_file('xbt034.xml')
        qresults = parse(xml_file, self.fmt)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('refseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('random_s00', qresult.description)
        self.assertEqual(128, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0, qresult.stat_kappa)
        self.assertEqual(0, qresult.stat_lambda)
        self.assertEqual(0, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_2', qresult.id)
        self.assertEqual('gi|347972582:1-35 Anopheles gambiae str. PEST AGAP011294-PA '
                '(DEFI_ANOGA) mRNA, complete cds', qresult.description)
        self.assertEqual(35, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0, qresult.stat_kappa)
        self.assertEqual(0, qresult.stat_lambda)
        self.assertEqual(0, qresult.stat_entropy)
        self.assertEqual(1, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|347972582|ref|XM_309352.4|', hit.id)
        self.assertEqual('Anopheles gambiae str. PEST AGAP011294-PA '
                '(DEFI_ANOGA) mRNA, complete cds', hit.description)
        self.assertEqual(309, hit.seq_length)
        self.assertEqual(3, len(hit))

        hsp = hit[0]
        self.assertEqual(34.0583, hsp.bitscore)
        self.assertEqual(68, hsp.bitscore_raw)
        self.assertEqual(1.38453, hsp.evalue)
        self.assertEqual(3, hsp.query_from)
        self.assertEqual(35, hsp.query_to)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(35, hsp.hit_to)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(3, hsp.hit_frame)
        self.assertEqual(11, hsp.identity_num)
        self.assertEqual(11, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, len(hsp))
        self.assertEqual('EVCNDRLYHCR', hsp.query.seq.tostring())
        self.assertEqual('EVCNDRLYHCR', hsp.hit.seq.tostring())
        self.assertEqual('EVCNDRLYHCR', hsp.homology)
        #self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_3', qresult.id)
        self.assertEqual('gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, '
                'complete cds', qresult.description)
        self.assertEqual(350, qresult.query_length)
        self.assertEqual(2922980, qresult.stat_db_num)
        self.assertEqual(4670255720, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0, qresult.stat_kappa)
        self.assertEqual(0, qresult.stat_lambda)
        self.assertEqual(0, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|296147483|ref|NM_001183135.1|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, '
                'complete cds >gi|116616412|gb|EF059095.1| Synthetic '
                'construct Saccharomyces cerevisiae clone FLH203015.01X '
                'MON2, complete sequence', hit.description)
        self.assertEqual(4911, hit.seq_length)
        self.assertEqual(7, len(hit))

        hsp = hit[0]
        self.assertEqual(289.739, hsp.bitscore)
        self.assertEqual(626, hsp.bitscore_raw)
        self.assertEqual(1.04512e-76, hsp.evalue)
        self.assertEqual(2, hsp.query_from)
        self.assertEqual(349, hsp.query_to)
        self.assertEqual(2, hsp.hit_from)
        self.assertEqual(349, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(116, hsp.identity_num)
        self.assertEqual(116, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(116, len(hsp))
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.hit.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.query.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.homology)

        hsp = hit[-1]
        self.assertEqual(36.3494, hsp.bitscore)
        self.assertEqual(73, hsp.bitscore_raw)
        self.assertEqual(9.00552e-54, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(42, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(42, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(14, hsp.identity_num)
        self.assertEqual(14, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(14, len(hsp))
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.query.seq.tostring())
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.hit.seq.tostring())
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|254579534|ref|XM_002495708.1|', hit.id)
        self.assertEqual('Zygosaccharomyces rouxii hypothetical protein '
                '(ZYRO0C02266g) mRNA, complete cds', hit.description)
        self.assertEqual(4866, hit.seq_length)
        self.assertEqual(4, len(hit))

        hsp = hit[0]
        self.assertEqual(141.279, hsp.bitscore)
        self.assertEqual(302, hsp.bitscore_raw)
        self.assertEqual(5.12803e-32, hsp.evalue)
        self.assertEqual(97, hsp.query_from)
        self.assertEqual(348, hsp.query_to)
        self.assertEqual(97, hsp.hit_from)
        self.assertEqual(348, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(57, hsp.identity_num)
        self.assertEqual(72, hsp.positive_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(84, len(hsp))
        self.assertEqual('IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', hsp.query.seq.tostring())
        self.assertEqual('IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK', hsp.hit.seq.tostring())
        self.assertEqual('IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)



if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
