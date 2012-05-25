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
FMT = 'blast-xml'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BlastnXmlCases(unittest.TestCase):

    def test_xbt002(self):
        xml_file = get_file('xbt002.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.12', qresult.meta['program_version'])
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual(1, qresult.meta['param_score_match'])
        self.assertEqual(-3, qresult.meta['param_score_mismatch'])
        self.assertEqual(5, qresult.meta['param_gap_open'])
        self.assertEqual(2, qresult.meta['param_gap_extend'])

        # test parsed values of qresult
        self.assertEqual('gi|1348916|gb|G26684.1|G26684', qresult.id)
        self.assertEqual('human STS STS_D11570, sequence tagged site', \
                qresult.desc)
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
        self.assertEqual('Pseudomonas aeruginosa PAO1, section 415 of 529 of the complete genome', hit.desc)
        self.assertEqual(11884, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(38.1576, hsp.bitscore)
        self.assertEqual(19, hsp.bitscore_raw)
        self.assertEqual(1.0598, hsp.evalue)
        self.assertEqual(68, hsp.query_from)
        self.assertEqual(86, hsp.query_to)
        self.assertEqual(6012, hsp.hit_from)
        self.assertEqual(6030, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(19, hsp.ident_num)
        self.assertEqual(19, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(19, len(hsp))
        self.assertEqual('CAGGCCAGCGACTTCTGGG', hsp.query.seq.tostring())
        self.assertEqual('CAGGCCAGCGACTTCTGGG', hsp.hit.seq.tostring())
        self.assertEqual('|||||||||||||||||||', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # parse last hit
        hit = qresult[-1]
        self.assertEqual('gi|15073988|emb|AL591786.1|SME591786', hit.id)
        self.assertEqual('Sinorhizobium meliloti 1021 complete chromosome; segment 5/12', hit.desc)
        self.assertEqual(299350, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(36.1753, hsp.bitscore)
        self.assertEqual(18, hsp.bitscore_raw)
        self.assertEqual(4.18768, hsp.evalue)
        self.assertEqual(204, hsp.query_from)
        self.assertEqual(224, hsp.query_to)
        self.assertEqual(83648, hsp.hit_from)
        self.assertEqual(83628, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(20, hsp.ident_num)
        self.assertEqual(20, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(21, len(hsp))
        self.assertEqual('TGAAAGGAAATNAAAATGGAA', hsp.query.seq.tostring())
        self.assertEqual('TGAAAGGAAATCAAAATGGAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||| |||||||||', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt012(self):
        xml_file = get_file('xbt012.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
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
        self.assertEqual('gi|356995852:1-490 Mus musculus POU domain, class 5, transcription factor 1 (Pou5f1), transcript variant 1, mRNA', qresult.desc)
        self.assertEqual(490, qresult.seq_len)
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
                'transcript variant 1, mRNA', hit.desc)
        self.assertEqual(1353, hit.seq_len)
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
        self.assertEqual(490, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
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
        self.assertEqual('hg19_dna range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 strand=+ repeatMasking=none', qresult.desc)
        self.assertEqual(66, qresult.seq_len)
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
                        hit.desc)
        self.assertEqual(7333, hit.seq_len)
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
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
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
        self.assertEqual(53, hsp.ident_num)
        self.assertEqual(53, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(53, len(hsp))
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('|||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|7', hit.id)
        self.assertEqual('gi|332865372|ref|XM_003318468.1| PREDICTED: Pan '
                'troglodytes zinc finger protein 273, transcript variant 1 '
                '(ZNF273), mRNA', hit.desc)
        self.assertEqual(4430, hit.seq_len)
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
        self.assertEqual(64, hsp.ident_num)
        self.assertEqual(64, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCACGCCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_xbt013(self):
        xml_file = get_file('xbt013.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
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

    def test_xbt014(self):
        xml_file = get_file('xbt014.xml')
        qresults = parse(xml_file, FMT)
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
                        qresult.desc)
        self.assertEqual(490, qresult.seq_len)
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
                'transcript variant 1, mRNA', hit.desc)
        self.assertEqual(1353, hit.seq_len)
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
        self.assertEqual(490, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, len(hsp))
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.query.seq.tostring())
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)


        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt015(self):
        xml_file = get_file('xbt015.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1
        self.assertEqual('Query_1', qresult.id)
        self.assertEqual('hg19_dna range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none',
                        qresult.desc)
        self.assertEqual(66, qresult.seq_len)
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
                        hit.desc)
        self.assertEqual(7333, hit.seq_len)
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
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
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
        self.assertEqual(53, hsp.ident_num)
        self.assertEqual(53, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(53, len(hsp))
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('|||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|7', hit.id)
        self.assertEqual('gi|332865372|ref|XM_003318468.1| PREDICTED: Pan '
                'troglodytes zinc finger protein 273, transcript variant 1 '
                '(ZNF273), mRNA', hit.desc)
        self.assertEqual(4430, hit.seq_len)
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
        self.assertEqual(64, hsp.ident_num)
        self.assertEqual(64, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCACGCCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt016(self):
        xml_file = get_file('xbt016.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
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
                        qresult.desc)
        self.assertEqual(490, qresult.seq_len)
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
                'transcript variant 1, mRNA', hit.desc)
        self.assertEqual(1353, hit.seq_len)
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
        self.assertEqual(490, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
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
                        qresult.desc)
        self.assertEqual(66, qresult.seq_len)
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
                        hit.desc)
        self.assertEqual(4771, hit.seq_len)
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
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, len(hsp))
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|332254616|ref|XM_003276378.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys S100P binding '
                'protein, transcript variant 2 (S100PBP), mRNA', \
                        hit.desc)
        self.assertEqual(4345, hit.seq_len)
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
        self.assertEqual(64, hsp.ident_num)
        self.assertEqual(64, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCATGCCACTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)


class BlastpXmlCases(unittest.TestCase):

    def test_xbt001(self):
        xml_file = get_file('xbt001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.12', qresult.meta['program_version'])
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('nr', qresult.target)
        
        self.assertEqual('gi|49176427|ref|NP_418280.3|', qresult.id)
        self.assertEqual('component of Sec-independent translocase [Escherichia coli K12]', \
                qresult.desc)
        self.assertEqual(103, qresult.seq_len)
        self.assertEqual(2934173, qresult.stat_db_num)
        self.assertEqual(1011751523, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(212, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|49176427|ref|NP_418280.3|', hit.id)
        self.assertEqual('component of Sec-independent translocase [Escherichia coli K12] >gi|26250604|ref|NP_756644.1| Sec-independent protein translocase protein tatA [Escherichia coli CFT073] >gi|30064867|ref|NP_839038.1| hypothetical protein S3840 [Shigella flexneri 2a str. 2457T] >gi|24115132|ref|NP_709642.1| hypothetical protein SF3914 [Shigella flexneri 2a str. 301] >gi|24054404|gb|AAN45349.1| orf, conserved hypothetical protein [Shigella flexneri 2a str. 301] >gi|2367310|gb|AAC76839.1| component of Sec-independent translocase [Escherichia coli K12] >gi|30043127|gb|AAP18849.1| hypothetical protein S3840 [Shigella flexneri 2a str. 2457T] >gi|26111035|gb|AAN83218.1| Sec-independent protein translocase protein tatA [Escherichia coli CFT073] >gi|3193217|gb|AAC19240.1| MttA1 [Escherichia coli] >gi|7444818|pir||E65188 hypothetical 11.3 kD protein in udp-rfaH intergenic region - Escherichia coli (strain K-12)', hit.desc)
        self.assertEqual(103, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(185.267, hsp.bitscore)
        self.assertEqual(469, hsp.bitscore_raw)
        self.assertEqual(4.20576e-46, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(103, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(103, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(103, hsp.ident_num)
        self.assertEqual(103, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(103, len(hsp))
        self.assertEqual('MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV', hsp.query.seq.tostring())
        self.assertEqual('MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV', hsp.hit.seq.tostring())
        self.assertEqual('MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # parse last hit
        hit = qresult[-1]
        self.assertEqual('gi|39593039|emb|CAE64508.1|', hit.id)
        self.assertEqual('Hypothetical protein CBG09238 [Caenorhabditis briggsae]', hit.desc)
        self.assertEqual(960, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(31.5722, hsp.bitscore)
        self.assertEqual(70, hsp.bitscore_raw)
        self.assertEqual(7.7721, hsp.evalue)
        self.assertEqual(55, hsp.query_from)
        self.assertEqual(102, hsp.query_to)
        self.assertEqual(410, hsp.hit_from)
        self.assertEqual(459, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(19, hsp.ident_num)
        self.assertEqual(33, hsp.pos_num)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(51, len(hsp))
        self.assertEqual('KAMSDDEPKQD---KTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ', hsp.query.seq.tostring())
        self.assertEqual('KKEADDKAKKDLEAKTKKEADEKAKKEADEKA-KKEAEAKTKEAEAKTKKE', hsp.hit.seq.tostring())
        self.assertEqual('K  +DD+ K+D   KT ++AD  AK  AD++A   + +AKT++A+   K++', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt006(self):
        xml_file = get_file('xbt006.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.18+', qresult.meta['program_version'])
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('nr', qresult.target)
        
        self.assertEqual('31493', qresult.id)
        self.assertEqual('unnamed protein product', qresult.desc)
        self.assertEqual(70, qresult.seq_len)
        self.assertEqual(15287, qresult.stat_db_num)
        self.assertEqual(7033566, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(10, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|151942244|gb|EDN60600.1|', hit.id)
        self.assertEqual('cytosolic iron-sulfur protein assembly protein [Saccharomyces cerevisiae YJM789]', hit.desc)
        self.assertEqual(330, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(33.113, hsp.bitscore)
        self.assertEqual(74, hsp.bitscore_raw)
        self.assertEqual(0.0185319, hsp.evalue)
        self.assertEqual(15, hsp.query_from)
        self.assertEqual(62, hsp.query_to)
        self.assertEqual(114, hsp.hit_from)
        self.assertEqual(163, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(16, hsp.ident_num)
        self.assertEqual(27, hsp.pos_num)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(50, len(hsp))
        self.assertEqual('AWNKDRTQIAICPNNHEVHIYE--KSGAKWNKVHELKEHNGQVTGIDWAP', hsp.query.seq.tostring())
        self.assertEqual('AWSNDGYYLATCSRDKSVWIWETDESGEEYECISVLQEHSQDVKHVIWHP', hsp.hit.seq.tostring())
        self.assertEqual('AW+ D   +A C  +  V I+E  +SG ++  +  L+EH+  V  + W P', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # parse last hit
        hit = qresult[-1]
        self.assertEqual('gi|151567870|pdb|2PM9|B', hit.id)
        self.assertEqual('Chain B, Crystal Structure Of Yeast Sec1331 VERTEX ELEMENT OF THE Copii Vesicular Coat', hit.desc)
        self.assertEqual(297, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(30.8018, hsp.bitscore)
        self.assertEqual(68, hsp.bitscore_raw)
        self.assertEqual(0.0919731, hsp.evalue)
        self.assertEqual(21, hsp.query_from)
        self.assertEqual(62, hsp.query_to)
        self.assertEqual(68, hsp.hit_from)
        self.assertEqual(109, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(11, hsp.ident_num)
        self.assertEqual(23, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(42, len(hsp))
        self.assertEqual('TQIAICPNNHEVHIYEKSGAKWNKVHELKEHNGQVTGIDWAP', hsp.query.seq.tostring())
        self.assertEqual('TILASCSYDGKVMIWKEENGRWSQIAVHAVHSASVNSVQWAP', hsp.hit.seq.tostring())
        self.assertEqual('T +A C  + +V I+++   +W+++     H+  V  + WAP', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt007(self):
        xml_file = get_file('xbt007.xml')
        qresults = parse(xml_file, FMT)
        counter = 0
        
        qresult = qresults.next()
        counter += 1
        
        # test meta variables, only for the first one
        self.assertEqual('2.2.18+', qresult.meta['program_version'])
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', 
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(0.01, qresult.meta['param_evalue_threshold'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('m L; R -d repeat/repeat_9606;', qresult.meta['param_filter'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('gpipe/9606/Previous/protein', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('gi|585505|sp|Q08386|MOPB_RHOCA', qresult.id)
        self.assertEqual('Molybdenum-pterin-binding protein mopB >gi|310278|gb|AAA71913.1| molybdenum-pterin-binding protein', qresult.desc)
        self.assertEqual(270, qresult.seq_len)
        self.assertEqual(27252, qresult.stat_db_num)
        self.assertEqual(13958303, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # second qresult
        qresult = qresults.next()
        counter += 1        

        self.assertEqual('gi|129628|sp|P07175.1|PARA_AGRTU', qresult.id)
        self.assertEqual('Protein parA', qresult.desc)
        self.assertEqual(222, qresult.seq_len)
        self.assertEqual(27252, qresult.stat_db_num)
        self.assertEqual(13958303, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(2, counter)

    def test_xbt008(self):
        xml_file = get_file('xbt008.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.18', qresult.meta['program_version'])
        self.assertEqual('~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.', 
                        qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(1e-05, qresult.meta['param_evalue_threshold'])
        self.assertEqual('F', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('/Users/pjcock/Downloads/Software/blast-2.2.18/data/nr', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('lcl|1_0', qresult.id)
        self.assertEqual('Fake', qresult.desc)
        self.assertEqual(9, qresult.seq_len)
        self.assertEqual(6589360, qresult.stat_db_num)
        self.assertEqual(2253133281, qresult.stat_db_len)
        self.assertEqual(2.02782e+10, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt010(self):
        xml_file = get_file('xbt010.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.22+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(1e-06, qresult.meta['param_evalue_threshold'])
        self.assertEqual('F', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('nr', qresult.target)
        
        self.assertEqual('1', qresult.id)
        self.assertEqual('gi|3298468|dbj|BAA31520.1| SAMIPF', qresult.desc)
        self.assertEqual(107, qresult.seq_len)
        self.assertEqual(8994603, qresult.stat_db_num)
        self.assertEqual(-1216159329, qresult.stat_db_len)
        self.assertEqual(76934807744, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(10, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|3298468|dbj|BAA31520.1|', hit.id)
        self.assertEqual('SAMIPF [Aster tripolium]', hit.desc)
        self.assertEqual(107, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(204.912011757068, hsp.bitscore)
        self.assertEqual(520, hsp.bitscore_raw)
        self.assertEqual(1.77242652875017e-51, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(107, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(107, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(107, hsp.ident_num)
        self.assertEqual(107, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(107, len(hsp))
        self.assertEqual('GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI', hsp.query.seq.tostring())
        self.assertEqual('GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI', hsp.hit.seq.tostring())
        self.assertEqual('GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # parse last hit
        hit = qresult[-1]
        self.assertEqual('gi|162809290|dbj|BAF95576.1|', hit.id)
        self.assertEqual('tonoplast intrinsic protein [Nicotiana tabacum]', hit.desc)
        self.assertEqual(251, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(177.948041442853, hsp.bitscore)
        self.assertEqual(450, hsp.bitscore_raw)
        self.assertEqual(2.0302699895292e-43, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(107, hsp.query_to)
        self.assertEqual(81, hsp.hit_from)
        self.assertEqual(187, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(91, hsp.ident_num)
        self.assertEqual(95, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(107, len(hsp))
        self.assertEqual('GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI', hsp.query.seq.tostring())
        self.assertEqual('GGHVNPAVTFGAFVGGNITLFRGILYIIAQLLGSTVACFLLEFATGGMSTGAFALSAGVSVWNAFVFEIVMTFGLVYTVYATAIDPKKGDLGVIAPIAIGFIVGANI', hsp.hit.seq.tostring())
        self.assertEqual('GGHVNPAVTFGAFVGGNITL RGI+YIIAQLLGSTVAC LL+F T  M+ G F+LSAGV V NA VFEIVMTFGLVYTVYATAIDPKKG LG IAPIAIGFIVGANI', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

    # def test_xbt011(self):
    # PSI-blast, handle later

    def test_xbt017(self):
        xml_file = get_file('xbt017.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
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
                'subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
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
                hit.desc)
        self.assertEqual(100, hit.seq_len)
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
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(81, hsp.pos_num)
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
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
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
                '[Mus musculus]', hit.desc)
        self.assertEqual(350, hit.seq_len)
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
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
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
        self.assertEqual(29, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|9', hit.id)
        self.assertEqual('gi|350596020|ref|XP_003360649.2| PREDICTED: '
                'pleckstrin-like [Sus scrofa]', hit.desc)
        self.assertEqual(228, hit.seq_len)
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
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_xbt018(self):
        xml_file = get_file('xbt018.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
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

    def test_xbt019(self):
        xml_file = get_file('xbt019.xml')
        qresults = parse(xml_file, FMT)
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
                'subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
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
                hit.desc)
        self.assertEqual(100, hit.seq_len)
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
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(81, hsp.pos_num)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(102, len(hsp))
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.query.seq.tostring())
        self.assertEqual('MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPVSFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN', hsp.hit.seq.tostring())
        self.assertEqual('MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt020(self):
        xml_file = get_file('xbt020.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
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
                '[Mus musculus]', hit.desc)
        self.assertEqual(350, hit.seq_len)
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
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
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
        self.assertEqual(29, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|9', hit.id)
        self.assertEqual('gi|350596020|ref|XP_003360649.2| PREDICTED: '
                'pleckstrin-like [Sus scrofa]', hit.desc)
        self.assertEqual(228, hit.seq_len)
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
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt021(self):
        xml_file = get_file('xbt021.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
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
                'subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual(12504708, qresult.stat_db_num)
        self.assertEqual(4346765131, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|16080617|ref|NP_391444.1|', hit.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168] >gi|221311516|ref|ZP_03593363.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168] >gi|221315843|ref|ZP_03597648.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. NCIB 3610] >gi|221320757|ref|ZP_03602051.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. JH642] >gi|221325043|ref|ZP_03606337.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. SMY] >gi|321313111|ref|YP_004205398.1| unnamed protein product [Bacillus subtilis BSn5]', hit.desc)
        self.assertEqual(102, hit.seq_len)
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
        self.assertEqual(102, hsp.ident_num)
        self.assertEqual(102, hsp.pos_num)
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
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual(12504708, qresult.stat_db_num)
        self.assertEqual(4346765131, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|11464971|ref|NP_062422.1|', hit.id)
        self.assertEqual('pleckstrin [Mus musculus]', hit.desc)
        self.assertEqual(350, hit.seq_len)
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
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
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
        self.assertEqual(29, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|350596020|ref|XP_003360649.2|', hit.id)
        self.assertEqual('PREDICTED: pleckstrin-like [Sus scrofa]', \
                hit.desc)
        self.assertEqual(228, hit.seq_len)
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
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)


class BlastxXmlCases(unittest.TestCase):

    def test_xbt003(self):
        xml_file = get_file('xbt003.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1
        
        self.assertEqual('2.2.12', qresult.meta['program_version'])
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(10.0, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('nr', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('gi|1347369|gb|G25137.1|G25137', qresult.id)
        self.assertEqual('human STS EST48004, sequence tagged site', qresult.desc)
        self.assertEqual(556, qresult.seq_len)
        self.assertEqual(2934173, qresult.stat_db_num)
        self.assertEqual(1011751523, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)

        # test parsed values of the first hit
        hit = qresult[0]
        self.assertEqual('gi|12654095|gb|AAH00859.1|', hit.id)
        self.assertEqual('Unknown (protein for IMAGE:3459481) [Homo sapiens]', \
                        hit.desc)
        self.assertEqual(319, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(247.284, hsp.bitscore)
        self.assertEqual(630, hsp.bitscore_raw)
        self.assertEqual(1.69599e-64, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(399, hsp.query_to)
        self.assertEqual(156, hsp.hit_from)
        self.assertEqual(288, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(122, hsp.ident_num)
        self.assertEqual(123, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(133, len(hsp))
        self.assertEqual('DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFXGPYTDFTPXTTEXQKLXEQALNTYPVNTXERWXXIAVAVPGRXKE', hsp.query.seq.tostring())
        self.assertEqual('DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFEGPYTDFTPWTTEEQKLLEQALKTYPVNTPERWEKIAEAVPGRTKK', hsp.hit.seq.tostring())
        self.assertEqual('DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERF GPYTDFTP TTE QKL EQAL TYPVNT ERW  IA AVPGR K+', hsp.homology)

        # test parsed values of last hit
        hit = qresult[-1]
        self.assertEqual('gi|72081091|ref|XP_800619.1|', hit.id)
        self.assertEqual('PREDICTED: hypothetical protein XP_795526 [Strongylocentrotus purpuratus]', hit.desc)
        self.assertEqual(337, hit.seq_len)

        hsp = hit[0]
        self.assertEqual(32.3426, hsp.bitscore)
        self.assertEqual(72, hsp.bitscore_raw)
        self.assertEqual(8.57476, hsp.evalue)
        self.assertEqual(40, hsp.query_from)
        self.assertEqual(231, hsp.query_to)
        self.assertEqual(106, hsp.hit_from)
        self.assertEqual(172, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(21, hsp.ident_num)
        self.assertEqual(37, hsp.pos_num)
        self.assertEqual(3, hsp.gap_num)
        self.assertEqual('AGTNSRWEVIANYMNI--HSSSGVKRT-AKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQ', hsp.query.seq.tostring())
        self.assertEqual('SSSNSSSKASASSSNVGASSSSGTKKSDSKSSNESSKSKRDKEDHKEGSINRSKDEKVSKEHRVVKE', hsp.hit.seq.tostring())
        self.assertEqual('+ +NS  +  A+  N+   SSSG K++ +K     +KS +  + H++  IN+   +K  KEH VV +', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt009(self):
        xml_file = get_file('xbt009.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1
        
        self.assertEqual('2.2.22+', qresult.meta['program_version'])
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(0.0001, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('nr', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('1', qresult.id)
        self.assertEqual('gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal RNA gene, partial sequence', qresult.desc)
        self.assertEqual(1002, qresult.seq_len)
        self.assertEqual(8994603, qresult.stat_db_num)
        self.assertEqual(-1216159329, qresult.stat_db_len)
        self.assertEqual(367397307882, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)

        # test parsed values of the first hit
        hit = qresult[0]
        self.assertEqual('gi|149390769|gb|ABR25402.1|', hit.id)
        self.assertEqual('unknown [Oryza sativa (indica cultivar-group)]', hit.desc)
        self.assertEqual(26, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(54.2989775733826, hsp.bitscore)
        self.assertEqual(129, hsp.bitscore_raw)
        self.assertEqual(1.83262460293058e-05, hsp.evalue)
        self.assertEqual(911, hsp.query_from)
        self.assertEqual(988, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(26, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(24, hsp.ident_num)
        self.assertEqual(25, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(26, len(hsp))
        self.assertEqual('HMLVSKIKPCMCKYEQIQTVKLRMAH', hsp.query.seq.tostring())
        self.assertEqual('HMLVSKIKPCMCKYELIRTVKLRMAH', hsp.hit.seq.tostring())
        self.assertEqual('HMLVSKIKPCMCKYE I+TVKLRMAH', hsp.homology)

    def test_xbt022(self):
        xml_file = get_file('xbt022.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
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
                        qresult.desc)
        self.assertEqual(490, qresult.seq_len)
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
                        hit.desc)
        self.assertEqual(352, hit.seq_len)
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
        self.assertEqual(140, hsp.ident_num)
        self.assertEqual(140, hsp.pos_num)
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
                        qresult.desc)
        self.assertEqual(485, qresult.seq_len)
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
                        hit.desc)
        self.assertEqual(132, hit.seq_len)
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
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(74, hsp.pos_num)
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
        self.assertEqual(34, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|19', hit.id)
        self.assertEqual('gi|33188429|ref|NP_872601.1| histone demethylase '
                'UTY isoform 1 [Homo sapiens]', hit.desc)
        self.assertEqual(1079, hit.seq_len)
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
        self.assertEqual(59, hsp.ident_num)
        self.assertEqual(66, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, len(hsp))
        self.assertEqual('SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ', hsp.query.seq.tostring())
        self.assertEqual('SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQACLELLTSGGLLASASQSAGITGVSHHAR', hsp.hit.seq.tostring())
        self.assertEqual('SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_xbt023(self):
        xml_file = get_file('xbt023.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
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

    def test_xbt024(self):
        xml_file = get_file('xbt024.xml')
        qresults = parse(xml_file, FMT)
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
                        qresult.desc)
        self.assertEqual(485, qresult.seq_len)
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
                        hit.desc)
        self.assertEqual(132, hit.seq_len)
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
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(74, hsp.pos_num)
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
        self.assertEqual(34, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|19', hit.id)
        self.assertEqual('gi|33188429|ref|NP_872601.1| histone demethylase '
                'UTY isoform 1 [Homo sapiens]', hit.desc)
        self.assertEqual(1079, hit.seq_len)
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
        self.assertEqual(59, hsp.ident_num)
        self.assertEqual(66, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, len(hsp))
        self.assertEqual('SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ', hsp.query.seq.tostring())
        self.assertEqual('SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQACLELLTSGGLLASASQSAGITGVSHHAR', hsp.hit.seq.tostring())
        self.assertEqual('SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt025(self):
        xml_file = get_file('xbt025.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
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
                        qresult.desc)
        self.assertEqual(490, qresult.seq_len)
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
                'isoform 1 [Mus musculus]', hit.desc)
        self.assertEqual(352, hit.seq_len)
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
        self.assertEqual(140, hsp.ident_num)
        self.assertEqual(140, hsp.pos_num)
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
                        qresult.desc)
        self.assertEqual(485, qresult.seq_len)
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
                'leucogenys]', hit.desc)
        self.assertEqual(132, hit.seq_len)
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
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(74, hsp.pos_num)
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
        self.assertEqual(34, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|57113895|ref|NP_001009002.1|', hit.id)
        self.assertEqual('histone demethylase UTY [Pan troglodytes]', \
                hit.desc)
        self.assertEqual(1079, hit.seq_len)
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
        self.assertEqual(59, hsp.ident_num)
        self.assertEqual(66, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, len(hsp))
        self.assertEqual('SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ', hsp.query.seq.tostring())
        self.assertEqual('SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQAHLELLTSGGLLASASQSAGITGVSHHAR', hsp.hit.seq.tostring())
        self.assertEqual('SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)


class TblastnXmlCases(unittest.TestCase):

    def test_xbt004(self):
        xml_file = get_file('xbt004.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.12', qresult.meta['program_version'])
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM62', qresult.meta['param_matrix'])
        self.assertEqual(0.001, qresult.meta['param_evalue_threshold'])
        self.assertEqual(11, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('nr', qresult.target)
        
        # test parsed values of qresult
        self.assertEqual('gi|729325|sp|P39483|DHG2_BACME', qresult.id)
        self.assertEqual('Glucose 1-dehydrogenase II (GLCDH-II)', \
                qresult.desc)
        self.assertEqual(261, qresult.seq_len)
        self.assertEqual(251887, qresult.stat_db_num)
        self.assertEqual(438542399, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(100, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|58264321|ref|XM_569317.1|', hit.id)
        self.assertEqual('Filobasidiella neoformans glucose 1-dehydrogenase, putative (CNB05760) mRNA, complete cds', hit.desc)
        self.assertEqual(904, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(148.288, hsp.bitscore)
        self.assertEqual(373, hsp.bitscore_raw)
        self.assertEqual(1.46834e-35, hsp.evalue)
        self.assertEqual(5, hsp.query_from)
        self.assertEqual(250, hsp.query_to)
        self.assertEqual(16, hsp.hit_from)
        self.assertEqual(762, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(84, hsp.ident_num)
        self.assertEqual(143, hsp.pos_num)
        self.assertEqual(9, hsp.gap_num)
        self.assertEqual(252, len(hsp))
        self.assertEqual('LKDKVVVVTGGSKGLGRAMAVRFGQEQSKVVVNYRSNXXXXXXXXXXXXXXXXGGQAIIVRGDVTKEEDVVNLVETAVKEFGSLDVMINNAGVENPVPSH---ELSLENWNQVIDTNLTGAFLGSREAIKYFVENDIKG-NVINMSSVHEMIPWPLFVHYAASKGGMKLMTETLALEYAPKGIRVNNIGPGAIDTPINAEKFADPEQRADVESMIPMGYIGKPEEIASVAAFLASSQASYVTGITLFADGGM', hsp.query.seq.tostring())
        self.assertEqual('LQGKVVAITGCSTGIGRAIAIGAAKNGANVVLHHLGDSTASDIAQVQEECKQAGAKTVVVPGDIAEAKTANEIVSAAVSSFSRIDVLISNAGI---CPFHSFLDLPHPLWKRVQDVNLNGSFYVVQAVANQMAKQEPKGGSIVAVSSISALMGGGEQCHYTPTKAGIKSLMESCAIALGPMGIRCNSVLPGTIETNINKEDLSNPEKRADQIRRVPLGRLGKPEDLVGPTLFFASDLSNYCTGASVLVDGGM', hsp.hit.seq.tostring())
        self.assertEqual('L+ KVV +TG S G+GRA+A+   +  + VV+++  +   +   + + E  +AG + ++V GD+ + +    +V  AV  F  +DV+I+NAG+    P H   +L    W +V D NL G+F   +       + + KG +++ +SS+  ++      HY  +K G+K + E+ A+   P GIR N++ PG I+T IN E  ++PE+RAD    +P+G +GKPE++     F AS  ++Y TG ++  DGGM', hsp.homology)

        # parse last hit
        hit = qresult[-1]
        self.assertEqual('gi|450259|gb|L27825.1|EMEVERA1AA', hit.id)
        self.assertEqual('Emericella nidulans (verA) gene, complete cds, ORF 1 gene, complete cds, and ORF 2 gene, 5\' end', hit.desc)
        self.assertEqual(4310, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(91.2781, hsp.bitscore)
        self.assertEqual(225, hsp.bitscore_raw)
        self.assertEqual(1.31998e-20, hsp.evalue)
        self.assertEqual(5, hsp.query_from)
        self.assertEqual(204, hsp.query_to)
        self.assertEqual(579, hsp.hit_from)
        self.assertEqual(1253, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(3, hsp.hit_frame)
        self.assertEqual(76, hsp.ident_num)
        self.assertEqual(113, hsp.pos_num)
        self.assertEqual(31, hsp.gap_num)
        self.assertEqual(228, len(hsp))
        self.assertEqual('LKDKVVVVTGGSKGLGRAMAVRFGQEQSKVVVNYRSNXXXXXXXXXXXXXXGGQAIIVRGDVTKEEDVVNLVETAVKEFGSLDVMINNAGV-----------ENPVPS-HELSLE-----NWNQVIDTNLTGAFLGSREAIKYFVENDIKGNVINMSSVHEMIPW-PLFVHYAASKGGMKLMTETLALEYAPKGIRVNNIGPGAIDTPI----------NAEKFADPE', hsp.query.seq.tostring())
        self.assertEqual('LDGKVALVTGAGRGIGAAIAVALGQPGAKVVVNYANSREAAEKVVDEIKSNAQSAISIQADVGDPDAVTKLMDQAVEHFGYLDIVSSNAGIVSFGHVKDVTPDVCVPSPYESPVEL*PQQEFDRVFRVNTRGQFFVAREAYRHLREG---GRIILTSSNTASVKGVPRHAVYSGSKGAIDTFVRCLAIDCGDKKITVNAVAPGAIKTDMFLSVSREYIPNGETFTDEQ', hsp.hit.seq.tostring())
        self.assertEqual('L  KV +VTG  +G+G A+AV  GQ  +KVVVNY ++ E A +V  EI+     AI ++ DV   + V  L++ AV+ FG LD++ +NAG+           +  VPS +E  +E      +++V   N  G F  +REA ++  E    G +I  SS    +   P    Y+ SKG +      LA++   K I VN + PGAI T +          N E F D +', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt026(self):
        xml_file = get_file('xbt026.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
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
                'subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
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
                'partial mRNA', hit.desc)
        self.assertEqual(4632, hit.seq_len)
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
        self.assertEqual(15, hsp.ident_num)
        self.assertEqual(26, hsp.pos_num)
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
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
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
                'scrofa pleckstrin-like (LOC100626968), mRNA', hit.desc)
        self.assertEqual(772, hit.seq_len)
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
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
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
        self.assertEqual(21, hsp.ident_num)
        self.assertEqual(33, hsp.pos_num)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(71, len(hsp))
        self.assertEqual('IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('+ +Y       P G I L+G  +TS      GK  F+ +     T  +  +F QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|17', hit.id)
        self.assertEqual('gi|338714227|ref|XM_001492113.3| PREDICTED: Equus '
                'caballus pleckstrin-like (LOC100051039), mRNA', hit.desc)
        self.assertEqual(1390, hit.seq_len)
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
        self.assertEqual(96, hsp.ident_num)
        self.assertEqual(97, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVK+GSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_xbt027(self):
        xml_file = get_file('xbt027.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt028(self):
        xml_file = get_file('xbt028.xml')
        qresults = parse(xml_file, FMT)
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
                'subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
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
                'partial mRNA', hit.desc)
        self.assertEqual(4632, hit.seq_len)
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
        self.assertEqual(15, hsp.ident_num)
        self.assertEqual(26, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(43, len(hsp))
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', hsp.query.seq.tostring())
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', hsp.hit.seq.tostring())
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.homology)
        self.assertRaises(IndexError, hit.__getitem__, 1)

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt029(self):
        xml_file = get_file('xbt029.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
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
                'scrofa pleckstrin-like (LOC100626968), mRNA', hit.desc)
        self.assertEqual(772, hit.seq_len)
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
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
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
        self.assertEqual(21, hsp.ident_num)
        self.assertEqual(33, hsp.pos_num)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(71, len(hsp))
        self.assertEqual('IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('+ +Y       P G I L+G  +TS      GK  F+ +     T  +  +F QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|17', hit.id)
        self.assertEqual('gi|338714227|ref|XM_001492113.3| PREDICTED: Equus '
                'caballus pleckstrin-like (LOC100051039), mRNA', hit.desc)
        self.assertEqual(1390, hit.seq_len)
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
        self.assertEqual(96, hsp.ident_num)
        self.assertEqual(97, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVK+GSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt030(self):
        xml_file = get_file('xbt030.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
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
                'subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
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
                '(GSPATT00004923001) partial mRNA', hit.desc)
        self.assertEqual(4632, hit.seq_len)
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
        self.assertEqual(15, hsp.ident_num)
        self.assertEqual(26, hsp.pos_num)
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
        self.assertEqual('gi|11464971:4-101 pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
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
                '(LOC100773128), mRNA', hit.desc)
        self.assertEqual(1119, hit.seq_len)
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
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
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
        self.assertEqual(30, hsp.ident_num)
        self.assertEqual(50, hsp.pos_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDF-GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAIHLRGCVVTSVESNHDGKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +  GK+     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|338714227|ref|XM_001492113.3|', hit.id)
        self.assertEqual('PREDICTED: Equus caballus pleckstrin-like '
                '(LOC100051039), mRNA', hit.desc)
        self.assertEqual(1390, hit.seq_len)
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
        self.assertEqual(96, hsp.ident_num)
        self.assertEqual(97, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVK+GSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)


class TblastxXmlCases(unittest.TestCase):

    def test_xbt005(self):
        xml_file = get_file('xbt005.xml')
        qresults = parse(xml_file, FMT)
        counter = 0
        
        # test the first qresult
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.12', qresult.meta['program_version'])
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.meta['program_reference'])
        self.assertEqual('BLOSUM80', qresult.meta['param_matrix'])
        self.assertEqual(1, qresult.meta['param_evalue_threshold'])
        self.assertEqual('L;', qresult.meta['param_filter'])
        self.assertEqual(10, qresult.meta['param_gap_open'])
        self.assertEqual(1, qresult.meta['param_gap_extend'])
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('nr', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('gi|1348853|gb|G26621.1|G26621', qresult.id)
        self.assertEqual('human STS STS_D12006, sequence tagged site', qresult.desc)
        self.assertEqual(615, qresult.seq_len)
        self.assertEqual(3533718, qresult.stat_db_num)
        # why is the value negative? is this a blast bug?
        self.assertEqual(-1496331888, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.177051, qresult.stat_kappa)
        self.assertEqual(0.342969, qresult.stat_lambda)
        self.assertEqual(0.656794, qresult.stat_entropy)
        self.assertEqual(10, len(qresult))        

        hit = qresult[0]
        self.assertEqual('gi|18072170|gb|AC010333.7|', hit.id)
        self.assertEqual('Homo sapiens chromosome 16 clone CTD-3037G24, complete sequence', hit.desc)
        self.assertEqual(159870, hit.seq_len)
        self.assertEqual(13, len(hit))

        hsp = hit[0]
        self.assertEqual(329.561, hsp.bitscore)
        self.assertEqual(661.0, hsp.bitscore_raw)
        self.assertEqual(5.29552e-90, hsp.evalue)
        self.assertEqual(2, hsp.query_from)
        self.assertEqual(355, hsp.query_to)
        self.assertEqual(44324, hsp.hit_from)
        self.assertEqual(44677, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(-3, hsp.hit_frame)
        self.assertEqual(117, hsp.ident_num)
        self.assertEqual(117, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(118, len(hsp))
        self.assertEqual('ECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD', hsp.query.seq.tostring())
        self.assertEqual('ECCFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD', hsp.hit.seq.tostring())
        self.assertEqual('EC FIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|4309961|gb|AC005993.2|AC005993', hit.id)
        self.assertEqual('Homo sapiens PAC clone RP6-114E22 from 14, complete sequence', hit.desc)
        self.assertEqual(143943, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(41.0922, hsp.bitscore)
        self.assertEqual(78, hsp.bitscore_raw)
        self.assertEqual(0.716571, hsp.evalue)
        self.assertEqual(167, hsp.query_from)
        self.assertEqual(250, hsp.query_to)
        self.assertEqual(43680, hsp.hit_from)
        self.assertEqual(43763, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(3, hsp.hit_frame)
        self.assertEqual(13, hsp.ident_num)
        self.assertEqual(19, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(28, len(hsp))
        self.assertEqual('PLTKAHRLFQTSIVFYVTCFTASSQQLL', hsp.query.seq.tostring())
        self.assertEqual('PLNKYHTIFQISLCFYLFCYNMAQKQLL', hsp.hit.seq.tostring())
        self.assertEqual('PL K H +FQ S+ FY+ C+  + +QLL', hsp.homology)        

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)
        
    def test_xbt031(self):
        xml_file = get_file('xbt031.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
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
                '(DEFI_ANOGA) mRNA, complete cds', qresult.desc)
        self.assertEqual(35, qresult.seq_len)
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
                hit.desc)
        self.assertEqual(309, hit.seq_len)
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
        self.assertEqual(11, hsp.ident_num)
        self.assertEqual(11, hsp.pos_num)
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
                'complete cds', qresult.desc)
        self.assertEqual(350, qresult.seq_len)
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
                'cerevisiae S288c Mon2p (MON2) mRNA, complete cds', hit.desc)
        self.assertEqual(4911, hit.seq_len)
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
        self.assertEqual(116, hsp.ident_num)
        self.assertEqual(116, hsp.pos_num)
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
        self.assertEqual(6, hsp.ident_num)
        self.assertEqual(9, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, len(hsp))
        self.assertEqual('KFWMPSLRLLI', hsp.query.seq.tostring())
        self.assertEqual('KKWVPPVKLLI', hsp.hit.seq.tostring())
        self.assertEqual('K W+P ++LLI', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|23', hit.id)
        self.assertEqual('gi|254579534|ref|XM_002495708.1| Zygosaccharomyces '
                'rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds', \
                hit.desc)
        self.assertEqual(4866, hit.seq_len)
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
        self.assertEqual(57, hsp.ident_num)
        self.assertEqual(72, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(84, len(hsp))
        self.assertEqual('IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', hsp.query.seq.tostring())
        self.assertEqual('IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK', hsp.hit.seq.tostring())
        self.assertEqual('IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_xbt032(self):
        xml_file = get_file('xbt032.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(24, qresult.stat_db_num)
        self.assertEqual(68522, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt033(self):
        xml_file = get_file('xbt033.xml')
        qresults = parse(xml_file, FMT)
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
                'complete cds', qresult.desc)
        self.assertEqual(350, qresult.seq_len)
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
                'cerevisiae S288c Mon2p (MON2) mRNA, complete cds', hit.desc)
        self.assertEqual(4911, hit.seq_len)
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
        self.assertEqual(116, hsp.ident_num)
        self.assertEqual(116, hsp.pos_num)
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
        self.assertEqual(6, hsp.ident_num)
        self.assertEqual(9, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, len(hsp))
        self.assertEqual('KFWMPSLRLLI', hsp.query.seq.tostring())
        self.assertEqual('KKWVPPVKLLI', hsp.hit.seq.tostring())
        self.assertEqual('K W+P ++LLI', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gnl|BL_ORD_ID|23', hit.id)
        self.assertEqual('gi|254579534|ref|XM_002495708.1| Zygosaccharomyces '
                'rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds', \
                hit.desc)
        self.assertEqual(4866, hit.seq_len)
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
        self.assertEqual(57, hsp.ident_num)
        self.assertEqual(72, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(84, len(hsp))
        self.assertEqual('IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', hsp.query.seq.tostring())
        self.assertEqual('IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK', hsp.hit.seq.tostring())
        self.assertEqual('IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK', hsp.homology)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xbt034(self):
        xml_file = get_file('xbt034.xml')
        qresults = parse(xml_file, FMT)
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
        self.assertEqual('random_s00', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
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
                '(DEFI_ANOGA) mRNA, complete cds', qresult.desc)
        self.assertEqual(35, qresult.seq_len)
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
                '(DEFI_ANOGA) mRNA, complete cds', hit.desc)
        self.assertEqual(309, hit.seq_len)
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
        self.assertEqual(11, hsp.ident_num)
        self.assertEqual(11, hsp.pos_num)
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
                'complete cds', qresult.desc)
        self.assertEqual(350, qresult.seq_len)
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
                'MON2, complete sequence', hit.desc)
        self.assertEqual(4911, hit.seq_len)
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
        self.assertEqual(116, hsp.ident_num)
        self.assertEqual(116, hsp.pos_num)
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
        self.assertEqual(14, hsp.ident_num)
        self.assertEqual(14, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(14, len(hsp))
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.query.seq.tostring())
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.hit.seq.tostring())
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.homology)

        hit = qresult[-1]
        self.assertEqual('gi|254579534|ref|XM_002495708.1|', hit.id)
        self.assertEqual('Zygosaccharomyces rouxii hypothetical protein '
                '(ZYRO0C02266g) mRNA, complete cds', hit.desc)
        self.assertEqual(4866, hit.seq_len)
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
        self.assertEqual(57, hsp.ident_num)
        self.assertEqual(72, hsp.pos_num)
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
