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

    def test_xml_2212L_blastn_001(self):
        xml_file = get_file('xml_2212L_blastn_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.12', qresult.version)
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual(1, qresult.param_score_match)
        self.assertEqual(-3, qresult.param_score_mismatch)
        self.assertEqual(5, qresult.param_gap_open)
        self.assertEqual(2, qresult.param_gap_extend)

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
        self.assertEqual(67, hsp.query_from)
        self.assertEqual(85, hsp.query_to)
        self.assertEqual(6011, hsp.hit_from)
        self.assertEqual(6029, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(19, hsp.ident_num)
        self.assertEqual(19, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(19, hsp.ali_len)
        self.assertEqual(19, len(hsp))
        self.assertEqual('CAGGCCAGCGACTTCTGGG', hsp.query.seq.tostring())
        self.assertEqual('CAGGCCAGCGACTTCTGGG', hsp.hit.seq.tostring())
        self.assertEqual('|||||||||||||||||||', hsp.alignment_annotation['homology'])
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
        self.assertEqual(203, hsp.query_from)
        self.assertEqual(223, hsp.query_to)
        self.assertEqual(83627, hsp.hit_from)
        self.assertEqual(83647, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(20, hsp.ident_num)
        self.assertEqual(20, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(21, hsp.ali_len)
        self.assertEqual(21, len(hsp))
        self.assertEqual('TGAAAGGAAATNAAAATGGAA', hsp.query.seq.tostring())
        self.assertEqual('TGAAAGGAAATCAAAATGGAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||| |||||||||', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_blastn_001(self):
        xml_file = get_file('xml_2226_blastn_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Zheng Zhang, Scott Schwartz, Lukas Wagner, and '
                'Webb Miller (2000), \"A greedy algorithm for aligning DNA '
                'sequences\", J Comput Biol 2000; 7(1-2):203-14.', \
                        qresult.reference)
        self.assertEqual(1.0, qresult.param_score_match)
        self.assertEqual(-2.0, qresult.param_score_mismatch)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;m;', qresult.param_filter)
        self.assertEqual(0.0, qresult.param_gap_open)
        self.assertEqual(0.0, qresult.param_gap_extend)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(7616765.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|356995852:1-490', qresult.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcription factor 1 (Pou5f1), transcript variant 1, mRNA', qresult.desc)
        self.assertEqual(490, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(31860807.0, qresult.stat_eff_space)
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
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(489, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(489, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(490, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, hsp.ali_len)
        self.assertEqual(490, len(hsp))
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.query.seq.tostring())
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 strand=+ repeatMasking=none', qresult.desc)
        self.assertEqual(66, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(3506256.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|94721341|ref|NM_001040441.1|', hit.id)
        self.assertEqual('Homo sapiens zinc '
                'finger and BTB domain containing 8A (ZBTB8A), mRNA', \
                        hit.desc)
        self.assertEqual(7333, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(115.613, hsp.bitscore)
        self.assertEqual(62, hsp.bitscore_raw)
        self.assertEqual(5.52066e-29, hsp.evalue)
        self.assertEqual(4, hsp.query_from)
        self.assertEqual(65, hsp.query_to)
        self.assertEqual(3676, hsp.hit_from)
        self.assertEqual(3737, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, hsp.ali_len)
        self.assertEqual(62, len(hsp))
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(98.9927, hsp.bitscore)
        self.assertEqual(53, hsp.bitscore_raw)
        self.assertEqual(5.55986e-24, hsp.evalue)
        self.assertEqual(5, hsp.query_from)
        self.assertEqual(57, hsp.query_to)
        self.assertEqual(2823, hsp.hit_from)
        self.assertEqual(2875, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(53, hsp.ident_num)
        self.assertEqual(53, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(53, hsp.ali_len)
        self.assertEqual(53, len(hsp))
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('|||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|332865372|ref|XM_003318468.1|', hit.id)
        self.assertEqual('PREDICTED: Pan '
                'troglodytes zinc finger protein 273, transcript variant 1 '
                '(ZNF273), mRNA', hit.desc)
        self.assertEqual(4430, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(111.919, hsp.bitscore)
        self.assertEqual(60, hsp.bitscore_raw)
        self.assertEqual(7.14143e-28, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(65, hsp.query_to)
        self.assertEqual(2734, hsp.hit_from)
        self.assertEqual(2799, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(64, hsp.ident_num)
        self.assertEqual(64, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, hsp.ali_len)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCACGCCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_xml_2226_blastn_002(self):
        xml_file = get_file('xml_2226_blastn_002.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Zheng Zhang, Scott Schwartz, Lukas Wagner, and '
                'Webb Miller (2000), \"A greedy algorithm for aligning DNA '
                'sequences\", J Comput Biol 2000; 7(1-2):203-14.', \
                        qresult.reference)
        self.assertEqual(1.0, qresult.param_score_match)
        self.assertEqual(-2.0, qresult.param_score_mismatch)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;m;', qresult.param_filter)
        self.assertEqual(0.0, qresult.param_gap_open)
        self.assertEqual(0.0, qresult.param_gap_extend)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(7616765.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))
        self.assertEqual([], qresult.hits)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_blastn_003(self):
        xml_file = get_file('xml_2226_blastn_003.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Zheng Zhang, Scott Schwartz, Lukas Wagner, and '
                'Webb Miller (2000), \"A greedy algorithm for aligning DNA '
                'sequences\", J Comput Biol 2000; 7(1-2):203-14.', \
                        qresult.reference)
        self.assertEqual(1.0, qresult.param_score_match)
        self.assertEqual(-2.0, qresult.param_score_mismatch)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;m;', qresult.param_filter)
        self.assertEqual(0.0, qresult.param_gap_open)
        self.assertEqual(0.0, qresult.param_gap_extend)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        self.assertEqual('gi|356995852:1-490', qresult.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcription '
                'factor 1 (Pou5f1), transcript variant 1, mRNA', \
                        qresult.desc)
        self.assertEqual(490, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(31860807.0, qresult.stat_eff_space)
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
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(489, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(489, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(490, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, hsp.ali_len)
        self.assertEqual(490, len(hsp))
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.query.seq.tostring())
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_blastn_004(self):
        xml_file = get_file('xml_2226_blastn_004.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none',
                        qresult.desc)
        self.assertEqual(66, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(3506256.0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|94721341|ref|NM_001040441.1|', hit.id)
        self.assertEqual('Homo sapiens zinc '
                'finger and BTB domain containing 8A (ZBTB8A), mRNA', \
                        hit.desc)
        self.assertEqual(7333, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(115.613, hsp.bitscore)
        self.assertEqual(62, hsp.bitscore_raw)
        self.assertEqual(5.52066e-29, hsp.evalue)
        self.assertEqual(4, hsp.query_from)
        self.assertEqual(65, hsp.query_to)
        self.assertEqual(3676, hsp.hit_from)
        self.assertEqual(3737, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, hsp.ali_len)
        self.assertEqual(62, len(hsp))
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(98.9927, hsp.bitscore)
        self.assertEqual(53, hsp.bitscore_raw)
        self.assertEqual(5.55986e-24, hsp.evalue)
        self.assertEqual(5, hsp.query_from)
        self.assertEqual(57, hsp.query_to)
        self.assertEqual(2823, hsp.hit_from)
        self.assertEqual(2875, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(53, hsp.ident_num)
        self.assertEqual(53, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(53, hsp.ali_len)
        self.assertEqual(53, len(hsp))
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('|||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|332865372|ref|XM_003318468.1|', hit.id)
        self.assertEqual('PREDICTED: Pan '
                'troglodytes zinc finger protein 273, transcript variant 1 '
                '(ZNF273), mRNA', hit.desc)
        self.assertEqual(4430, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(111.919, hsp.bitscore)
        self.assertEqual(60, hsp.bitscore_raw)
        self.assertEqual(7.14143e-28, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(65, hsp.query_to)
        self.assertEqual(2734, hsp.hit_from)
        self.assertEqual(2799, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(64, hsp.ident_num)
        self.assertEqual(64, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, hsp.ali_len)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCACGCCATTGCACTCCAGCCTGGGCAACAAGAGTGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_blastn_005(self):
        xml_file = get_file('xml_2226_blastn_005.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Zheng Zhang, Scott Schwartz, Lukas Wagner, and '
                'Webb Miller (2000), \"A greedy algorithm for aligning DNA '
                'sequences\", J Comput Biol 2000; 7(1-2):203-14.', \
                        qresult.reference)
        self.assertEqual(1.0, qresult.param_score_match)
        self.assertEqual(-2.0, qresult.param_score_mismatch)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;m;', qresult.param_filter)
        self.assertEqual(0.0, qresult.param_gap_open)
        self.assertEqual(0.0, qresult.param_gap_extend)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('refseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(2933984, qresult.stat_db_num)
        self.assertEqual(4726730735, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|356995852:1-490', qresult.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcription '
                'factor 1 (Pou5f1), transcript variant 1, mRNA', \
                        qresult.desc)
        self.assertEqual(490, qresult.seq_len)
        self.assertEqual(2933984, qresult.stat_db_num)
        self.assertEqual(4726730735, qresult.stat_db_len)
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
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(489, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(489, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(490, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, hsp.ali_len)
        self.assertEqual(490, len(hsp))
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.query.seq.tostring())
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGGCCCCCGGCTCGGGGTGCCCACCTTCCCCATGGCTGGACACCTGGCTTCAGACTTCGCCTTCTCACCCCCACCAGGTGGGGGTGATGGGTCAGCAGGGCTGGAGCCGGGCTGGGTGGATCCTCGAACCTGGCTAAGCTTCCAAGGGCCTCCAGGTGGGCCTGGAATCGGACCAGGCTCAGAGGTATTGGGGATCTCCCCATGTCCGCCCGCATACGAGTTCTGCGGAGGGATGGCATACTGTGGACCTCAGGTTGGACTGGGCCTAGTCCCCCAAGTTGGCGTGGAGACTTTGCAGCCTGAGGGCCAGGCAGGAGCACGAGTGGAAAGCAACTCAGAGGGAACCTCCTCTGAGCCCTGTGCCGACCGCCCCAATGCCGTGAAGTTGGAGAAGGTGGAACCAACTCCCGAGGAGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none', qresult.desc)
        self.assertEqual(66, qresult.seq_len)
        self.assertEqual(2933984, qresult.stat_db_num)
        self.assertEqual(4726730735, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.46, qresult.stat_kappa)
        self.assertEqual(1.28, qresult.stat_lambda)
        self.assertEqual(0.85, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|332237160|ref|XM_003267724.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys ATG14 autophagy '
                'related 14 homolog (S. cerevisiae) (ATG14), mRNA', hit.desc)
        self.assertEqual(4771, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(115.613, hsp.bitscore)
        self.assertEqual(62, hsp.bitscore_raw)
        self.assertEqual(3.35972e-23, hsp.evalue)
        self.assertEqual(4, hsp.query_from)
        self.assertEqual(65, hsp.query_to)
        self.assertEqual(2864, hsp.hit_from)
        self.assertEqual(2925, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, hsp.ali_len)
        self.assertEqual(62, len(hsp))
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|332254616|ref|XM_003276378.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys S100P binding '
                'protein, transcript variant 2 (S100PBP), mRNA', hit.desc)
        self.assertEqual(4345, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(111.919, hsp.bitscore)
        self.assertEqual(60, hsp.bitscore_raw)
        self.assertEqual(4.34607e-22, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(65, hsp.query_to)
        self.assertEqual(2791, hsp.hit_from)
        self.assertEqual(2856, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(64, hsp.ident_num)
        self.assertEqual(64, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(66, hsp.ali_len)
        self.assertEqual(66, len(hsp))
        self.assertEqual('TCAAGCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.query.seq.tostring())
        self.assertEqual('TCATGCCACTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', hsp.hit.seq.tostring())
        self.assertEqual('||| |||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)


class BlastpXmlCases(unittest.TestCase):

    def test_xml_2212L_blastp_001(self):
        xml_file = get_file('xml_2212L_blastp_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.12', qresult.version)
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
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
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(102, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(102, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(103, hsp.ident_num)
        self.assertEqual(103, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(103, hsp.ali_len)
        self.assertEqual(103, len(hsp))
        self.assertEqual('MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV', hsp.query.seq.tostring())
        self.assertEqual('MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV', hsp.hit.seq.tostring())
        self.assertEqual('MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV', hsp.alignment_annotation['homology'])
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
        self.assertEqual(54, hsp.query_from)
        self.assertEqual(101, hsp.query_to)
        self.assertEqual(409, hsp.hit_from)
        self.assertEqual(458, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(19, hsp.ident_num)
        self.assertEqual(33, hsp.pos_num)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(51, hsp.ali_len)
        self.assertEqual(51, len(hsp))
        self.assertEqual('KAMSDDEPKQD---KTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ', hsp.query.seq.tostring())
        self.assertEqual('KKEADDKAKKDLEAKTKKEADEKAKKEADEKA-KKEAEAKTKEAEAKTKKE', hsp.hit.seq.tostring())
        self.assertEqual('K  +DD+ K+D   KT ++AD  AK  AD++A   + +AKT++A+   K++', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2218_blastp_001(self):
        xml_file = get_file('xml_2218_blastp_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.18+', qresult.version)
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
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
        self.assertEqual(14, hsp.query_from)
        self.assertEqual(61, hsp.query_to)
        self.assertEqual(113, hsp.hit_from)
        self.assertEqual(162, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(16, hsp.ident_num)
        self.assertEqual(27, hsp.pos_num)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(50, hsp.ali_len)
        self.assertEqual(50, len(hsp))
        self.assertEqual('AWNKDRTQIAICPNNHEVHIYE--KSGAKWNKVHELKEHNGQVTGIDWAP', hsp.query.seq.tostring())
        self.assertEqual('AWSNDGYYLATCSRDKSVWIWETDESGEEYECISVLQEHSQDVKHVIWHP', hsp.hit.seq.tostring())
        self.assertEqual('AW+ D   +A C  +  V I+E  +SG ++  +  L+EH+  V  + W P', hsp.alignment_annotation['homology'])
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
        self.assertEqual(20, hsp.query_from)
        self.assertEqual(61, hsp.query_to)
        self.assertEqual(67, hsp.hit_from)
        self.assertEqual(108, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(11, hsp.ident_num)
        self.assertEqual(23, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(42, hsp.ali_len)
        self.assertEqual(42, len(hsp))
        self.assertEqual('TQIAICPNNHEVHIYEKSGAKWNKVHELKEHNGQVTGIDWAP', hsp.query.seq.tostring())
        self.assertEqual('TILASCSYDGKVMIWKEENGRWSQIAVHAVHSASVNSVQWAP', hsp.hit.seq.tostring())
        self.assertEqual('T +A C  + +V I+++   +W+++     H+  V  + WAP', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2218_blastp_002(self):
        xml_file = get_file('xml_2218_blastp_002.xml')
        qresults = parse(xml_file, FMT)
        counter = 0
        
        qresult = qresults.next()
        counter += 1
        
        # test meta variables, only for the first one
        self.assertEqual('2.2.18+', qresult.version)
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', 
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(0.01, qresult.param_evalue_threshold)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('m L; R -d repeat/repeat_9606;', qresult.param_filter)
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

    def test_xml_2218L_blastp_001(self):
        xml_file = get_file('xml_2218L_blastp_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.18', qresult.version)
        self.assertEqual('~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.', 
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(1e-05, qresult.param_evalue_threshold)
        self.assertEqual('F', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('/Users/pjcock/Downloads/Software/blast-2.2.18/data/nr', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('Fake', qresult.id)
        self.assertEqual('', qresult.desc)
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

    def test_xml_2222_blastp_001(self):
        xml_file = get_file('xml_2222_blastp_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.22+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(1e-06, qresult.param_evalue_threshold)
        self.assertEqual('F', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
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
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(106, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(106, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(107, hsp.ident_num)
        self.assertEqual(107, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(107, hsp.ali_len)
        self.assertEqual(107, len(hsp))
        self.assertEqual('GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI', hsp.query.seq.tostring())
        self.assertEqual('GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI', hsp.hit.seq.tostring())
        self.assertEqual('GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI', hsp.alignment_annotation['homology'])
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
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(106, hsp.query_to)
        self.assertEqual(80, hsp.hit_from)
        self.assertEqual(186, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(91, hsp.ident_num)
        self.assertEqual(95, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(107, hsp.ali_len)
        self.assertEqual(107, len(hsp))
        self.assertEqual('GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVGVTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI', hsp.query.seq.tostring())
        self.assertEqual('GGHVNPAVTFGAFVGGNITLFRGILYIIAQLLGSTVACFLLEFATGGMSTGAFALSAGVSVWNAFVFEIVMTFGLVYTVYATAIDPKKGDLGVIAPIAIGFIVGANI', hsp.hit.seq.tostring())
        self.assertEqual('GGHVNPAVTFGAFVGGNITL RGI+YIIAQLLGSTVAC LL+F T  M+ G F+LSAGV V NA VFEIVMTFGLVYTVYATAIDPKKG LG IAPIAIGFIVGANI', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

    # def test_xml_2218L_rpsblast_001(self):
    # PSI-blast, handle later

    def test_xml_2226_blastp_001(self):
        xml_file = get_file('xml_2226_blastp_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('F', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('db/minirefseq_prot', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
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
        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis '
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
        self.assertEqual('gi|308175296|ref|YP_003922001.1|', hit.id)
        self.assertEqual('membrane bound '
                'lipoprotein [Bacillus amyloliquefaciens DSM 7]', \
                hit.desc)
        self.assertEqual(100, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(139.428, hsp.bitscore)
        self.assertEqual(350, hsp.bitscore_raw)
        self.assertEqual(1.99275e-46, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(101, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(99, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(81, hsp.pos_num)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(102, hsp.ali_len)
        self.assertEqual(102, len(hsp))
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.query.seq.tostring())
        self.assertEqual('MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPVSFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN', hsp.hit.seq.tostring())
        self.assertEqual('MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(345626, qresult.stat_eff_space)
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
        self.assertEqual(2.24956e-69, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(43.5134, hsp.bitscore)
        self.assertEqual(101, hsp.bitscore_raw)
        self.assertEqual(2.90061e-09, hsp.evalue)
        self.assertEqual(2, hsp.query_from)
        self.assertEqual(95, hsp.query_to)
        self.assertEqual(245, hsp.hit_from)
        self.assertEqual(344, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(29, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, hsp.ali_len)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|350596020|ref|XP_003360649.2|', hit.id)
        self.assertEqual('PREDICTED: pleckstrin-like [Sus scrofa]', hit.desc)
        self.assertEqual(228, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.97058e-68, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_xml_2226_blastp_002(self):
        xml_file = get_file('xml_2226_blastp_002.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('F', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('db/minirefseq_prot', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
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

    def test_xml_2226_blastp_003(self):
        xml_file = get_file('xml_2226_blastp_003.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('F', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('db/minirefseq_prot', qresult.target)

        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis '
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
        self.assertEqual('gi|308175296|ref|YP_003922001.1|', hit.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus amyloliquefaciens DSM 7]', hit.desc)
        self.assertEqual(100, hit.seq_len)
        self.assertEqual(1, len(hit))
        self.assertRaises(IndexError, hit.__getitem__, 1)

        hsp = hit[0]
        self.assertEqual(139.428, hsp.bitscore)
        self.assertEqual(350, hsp.bitscore_raw)
        self.assertEqual(1.99275e-46, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(101, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(99, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(81, hsp.pos_num)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(102, hsp.ali_len)
        self.assertEqual(102, len(hsp))
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.query.seq.tostring())
        self.assertEqual('MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPVSFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN', hsp.hit.seq.tostring())
        self.assertEqual('MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_blastp_004(self):
        xml_file = get_file('xml_2226_blastp_004.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('F', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('db/minirefseq_prot', qresult.target)

        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(345626, qresult.stat_eff_space)
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
        self.assertEqual(2.24956e-69, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(43.5134, hsp.bitscore)
        self.assertEqual(101, hsp.bitscore_raw)
        self.assertEqual(2.90061e-09, hsp.evalue)
        self.assertEqual(2, hsp.query_from)
        self.assertEqual(95, hsp.query_to)
        self.assertEqual(245, hsp.hit_from)
        self.assertEqual(344, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(29, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, hsp.ali_len)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|350596020|ref|XP_003360649.2|', hit.id)
        self.assertEqual('PREDICTED: pleckstrin-like [Sus scrofa]', hit.desc)
        self.assertEqual(228, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.97058e-68, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_blastp_005(self):
        xml_file = get_file('xml_2226_blastp_005.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), \"Gapped BLAST and '
                'PSI-BLAST: a new generation of protein database search '
                'programs\", Nucleic Acids Res. 25:3389-3402.', \
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('F', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('refseq_protein', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(12646943, qresult.stat_db_num)
        self.assertEqual(4397139428, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual(12646943, qresult.stat_db_num)
        self.assertEqual(4397139428, qresult.stat_db_len)
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
        self.assertEqual(1.45285e-66, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(101, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(101, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(102, hsp.ident_num)
        self.assertEqual(102, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(102, hsp.ali_len)
        self.assertEqual(102, len(hsp))
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.query.seq.tostring())
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.hit.seq.tostring())
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual(12646943, qresult.stat_db_num)
        self.assertEqual(4397139428, qresult.stat_db_len)
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
        self.assertEqual(1.54412e-63, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(43.5134, hsp.bitscore)
        self.assertEqual(101, hsp.bitscore_raw)
        self.assertEqual(0.00199101, hsp.evalue)
        self.assertEqual(2, hsp.query_from)
        self.assertEqual(95, hsp.query_to)
        self.assertEqual(245, hsp.hit_from)
        self.assertEqual(344, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(29, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, hsp.ali_len)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|350596020|ref|XP_003360649.2|', hit.id)
        self.assertEqual('PREDICTED: pleckstrin-like [Sus scrofa]', hit.desc)
        self.assertEqual(228, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.35263e-62, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)


class BlastxXmlCases(unittest.TestCase):

    def test_xml_2212L_blastx_001(self):
        xml_file = get_file('xml_2212L_blastx_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1
        
        self.assertEqual('2.2.12', qresult.version)
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
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
        self.assertEqual('Unknown (protein for IMAGE:3459481) [Homo sapiens]', hit.desc)
        self.assertEqual(319, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(247.284, hsp.bitscore)
        self.assertEqual(630, hsp.bitscore_raw)
        self.assertEqual(1.69599e-64, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(398, hsp.query_to)
        self.assertEqual(155, hsp.hit_from)
        self.assertEqual(287, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(122, hsp.ident_num)
        self.assertEqual(123, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(133, hsp.ali_len)
        self.assertEqual(133, len(hsp))
        self.assertEqual('DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFXGPYTDFTPXTTEXQKLXEQALNTYPVNTXERWXXIAVAVPGRXKE', hsp.query.seq.tostring())
        self.assertEqual('DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFEGPYTDFTPWTTEEQKLLEQALKTYPVNTPERWEKIAEAVPGRTKK', hsp.hit.seq.tostring())
        self.assertEqual('DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERF GPYTDFTP TTE QKL EQAL TYPVNT ERW  IA AVPGR K+', hsp.alignment_annotation['homology'])

        # test parsed values of last hit
        hit = qresult[-1]
        self.assertEqual('gi|72081091|ref|XP_800619.1|', hit.id)
        self.assertEqual('PREDICTED: hypothetical protein XP_795526 [Strongylocentrotus purpuratus]', hit.desc)
        self.assertEqual(337, hit.seq_len)

        hsp = hit[0]
        self.assertEqual(32.3426, hsp.bitscore)
        self.assertEqual(72, hsp.bitscore_raw)
        self.assertEqual(8.57476, hsp.evalue)
        self.assertEqual(39, hsp.query_from)
        self.assertEqual(230, hsp.query_to)
        self.assertEqual(105, hsp.hit_from)
        self.assertEqual(171, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(21, hsp.ident_num)
        self.assertEqual(37, hsp.pos_num)
        self.assertEqual(3, hsp.gap_num)
        self.assertEqual('AGTNSRWEVIANYMNI--HSSSGVKRT-AKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQ', hsp.query.seq.tostring())
        self.assertEqual('SSSNSSSKASASSSNVGASSSSGTKKSDSKSSNESSKSKRDKEDHKEGSINRSKDEKVSKEHRVVKE', hsp.hit.seq.tostring())
        self.assertEqual('+ +NS  +  A+  N+   SSSG K++ +K     +KS +  + H++  IN+   +K  KEH VV +', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2222_blastx_001(self):
        xml_file = get_file('xml_2222_blastx_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1
        
        self.assertEqual('2.2.22+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(0.0001, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
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
        self.assertEqual(910, hsp.query_from)
        self.assertEqual(987, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(25, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(24, hsp.ident_num)
        self.assertEqual(25, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(26, hsp.ali_len)
        self.assertEqual(26, len(hsp))
        self.assertEqual('HMLVSKIKPCMCKYEQIQTVKLRMAH', hsp.query.seq.tostring())
        self.assertEqual('HMLVSKIKPCMCKYELIRTVKLRMAH', hsp.hit.seq.tostring())
        self.assertEqual('HMLVSKIKPCMCKYE I+TVKLRMAH', hsp.alignment_annotation['homology'])

    def test_xml_2226_blastx_001(self):
        xml_file = get_file('xml_2226_blastx_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), "Gapped BLAST '
                'and PSI-BLAST: a new generation of protein database search '
                'programs", Nucleic Acids Res. 25:3389-3402.',
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('db/minirefseq_prot', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none', qresult.desc)
        self.assertEqual(485, qresult.seq_len)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(662354, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|332258565|ref|XP_003278367.1|', hit.id)
        self.assertEqual('PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]', hit.desc)
        self.assertEqual(132, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(121.709, hsp.bitscore)
        self.assertEqual(304, hsp.bitscore_raw)
        self.assertEqual(2.9522e-38, hsp.evalue)
        self.assertEqual(15, hsp.query_from)
        self.assertEqual(299, hsp.query_to)
        self.assertEqual(24, hsp.hit_from)
        self.assertEqual(118, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(74, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(95, hsp.ali_len)
        self.assertEqual(95, len(hsp))
        self.assertEqual('LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', hsp.query.seq.tostring())
        self.assertEqual('LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP', hsp.hit.seq.tostring())
        self.assertEqual('LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(51.6026, hsp.bitscore)
        self.assertEqual(122, hsp.bitscore_raw)
        self.assertEqual(2.73605e-12, hsp.evalue)
        self.assertEqual(243, hsp.query_from)
        self.assertEqual(458, hsp.query_to)
        self.assertEqual(31, hsp.hit_from)
        self.assertEqual(97, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(34, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, hsp.ali_len)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|33188429|ref|NP_872601.1|', hit.id)
        self.assertEqual('histone demethylase UTY isoform 1 [Homo sapiens]', hit.desc)
        self.assertEqual(1079, hit.seq_len)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(104.375, hsp.bitscore)
        self.assertEqual(259, hsp.bitscore_raw)
        self.assertEqual(6.31914e-29, hsp.evalue)
        self.assertEqual(18, hsp.query_from)
        self.assertEqual(290, hsp.query_to)
        self.assertEqual(988, hsp.hit_from)
        self.assertEqual(1078, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(59, hsp.ident_num)
        self.assertEqual(66, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, hsp.ali_len)
        self.assertEqual(91, len(hsp))
        self.assertEqual('SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ', hsp.query.seq.tostring())
        self.assertEqual('SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQACLELLTSGGLLASASQSAGITGVSHHAR', hsp.hit.seq.tostring())
        self.assertEqual('SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(2, counter)

    def test_xml_2226_blastx_002(self):
        xml_file = get_file('xml_2226_blastx_002.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), "Gapped BLAST '
                'and PSI-BLAST: a new generation of protein database search '
                'programs", Nucleic Acids Res. 25:3389-3402.',
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('db/minirefseq_prot', qresult.target)

        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
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

    def test_xml_2226_blastx_003(self):
        xml_file = get_file('xml_2226_blastx_003.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), "Gapped BLAST '
                'and PSI-BLAST: a new generation of protein database search '
                'programs", Nucleic Acids Res. 25:3389-3402.',
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('db/minirefseq_prot', qresult.target)

        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 strand=+ repeatMasking=none', qresult.desc)
        self.assertEqual(485, qresult.seq_len)
        self.assertEqual(20, qresult.stat_db_num)
        self.assertEqual(6406, qresult.stat_db_len)
        self.assertEqual(662354, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|332258565|ref|XP_003278367.1|', hit.id)
        self.assertEqual('PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]', hit.desc)
        self.assertEqual(132, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(121.709, hsp.bitscore)
        self.assertEqual(304, hsp.bitscore_raw)
        self.assertEqual(2.9522e-38, hsp.evalue)
        self.assertEqual(15, hsp.query_from)
        self.assertEqual(299, hsp.query_to)
        self.assertEqual(24, hsp.hit_from)
        self.assertEqual(118, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(74, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(95, hsp.ali_len)
        self.assertEqual(95, len(hsp))
        self.assertEqual('LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', hsp.query.seq.tostring())
        self.assertEqual('LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP', hsp.hit.seq.tostring())
        self.assertEqual('LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(51.6026, hsp.bitscore)
        self.assertEqual(122, hsp.bitscore_raw)
        self.assertEqual(2.73605e-12, hsp.evalue)
        self.assertEqual(243, hsp.query_from)
        self.assertEqual(458, hsp.query_to)
        self.assertEqual(31, hsp.hit_from)
        self.assertEqual(97, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(34, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, hsp.ali_len)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|33188429|ref|NP_872601.1|', hit.id)
        self.assertEqual('histone demethylase UTY isoform 1 [Homo sapiens]', hit.desc)
        self.assertEqual(1079, hit.seq_len)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(104.375, hsp.bitscore)
        self.assertEqual(259, hsp.bitscore_raw)
        self.assertEqual(6.31914e-29, hsp.evalue)
        self.assertEqual(18, hsp.query_from)
        self.assertEqual(290, hsp.query_to)
        self.assertEqual(988, hsp.hit_from)
        self.assertEqual(1078, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(59, hsp.ident_num)
        self.assertEqual(66, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, hsp.ali_len)
        self.assertEqual(91, len(hsp))
        self.assertEqual('SFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQ', hsp.query.seq.tostring())
        self.assertEqual('SFQESLRAGMQWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQACLELLTSGGLLASASQSAGITGVSHHAR', hsp.hit.seq.tostring())
        self.assertEqual('SF    +AG+QW DL   QPPPPGFK FS LS P+SW+YRH+P C  NF   VETGF+HVGQA LE   SG L A ASQS GITGVSHHA+', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_blastx_004(self):
        xml_file = get_file('xml_2226_blastx_004.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb '
                'Miller, and David J. Lipman (1997), "Gapped BLAST '
                'and PSI-BLAST: a new generation of protein database search '
                'programs", Nucleic Acids Res. 25:3389-3402.',
                        qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('refseq_protein', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(12646943, qresult.stat_db_num)
        self.assertEqual(4397139428, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 '
                'strand=+ repeatMasking=none', qresult.desc)
        self.assertEqual(485, qresult.seq_len)
        self.assertEqual(12646943, qresult.stat_db_num)
        self.assertEqual(4397139428, qresult.stat_db_len)
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
        self.assertEqual(2.02642e-32, hsp.evalue)
        self.assertEqual(15, hsp.query_from)
        self.assertEqual(299, hsp.query_to)
        self.assertEqual(24, hsp.hit_from)
        self.assertEqual(118, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(74, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(95, hsp.ali_len)
        self.assertEqual(95, len(hsp))
        self.assertEqual('LRRSFALVAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', hsp.query.seq.tostring())
        self.assertEqual('LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLANFLFLVEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP', hsp.hit.seq.tostring())
        self.assertEqual('LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW+YRH+PP L NF+FLVE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(51.6026, hsp.bitscore)
        self.assertEqual(122, hsp.bitscore_raw)
        self.assertEqual(1.87805e-06, hsp.evalue)
        self.assertEqual(243, hsp.query_from)
        self.assertEqual(458, hsp.query_to)
        self.assertEqual(31, hsp.hit_from)
        self.assertEqual(97, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(34, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(72, hsp.ali_len)
        self.assertEqual(72, len(hsp))
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPPHPANXXXXXXXXFLRRSFALVAQAGVQWLDLG', hsp.query.seq.tostring())
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', hsp.hit.seq.tostring())
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|332815399|ref|XP_003309509.1|', hit.id)
        self.assertEqual('PREDICTED: histone demethylase UTY-like [Pan troglodytes]', hit.desc)
        self.assertEqual(101, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(97.0561, hsp.bitscore)
        self.assertEqual(240, hsp.bitscore_raw)
        self.assertEqual(2.76414e-23, hsp.evalue)
        self.assertEqual(6, hsp.query_from)
        self.assertEqual(278, hsp.query_to)
        self.assertEqual(9, hsp.hit_from)
        self.assertEqual(99, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(56, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(91, hsp.ali_len)
        self.assertEqual(91, len(hsp))
        self.assertEqual('VAQAGVQWLDLGXXXXXXPGFK*FSCLSHPSSWDYRHMPPCLINFVFLVETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQPLCE', hsp.query.seq.tostring())
        self.assertEqual('VPHAGVQWHNLSSLQPPPSRFKPFSYLSLLSSWDQRRPPPCLVTFVFLIETGFRHVGQAGLKLLTSGDPSASASQSAGIRGVSHCTWPECQ', hsp.hit.seq.tostring())
        self.assertEqual('V  AGVQW +L   QPPP  FK FS LS  SSWD R  PPCL+ FVFL+ETGF HVGQAGL+   SG+  A ASQS GI GVSH   P C+', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(2, counter)


class TblastnXmlCases(unittest.TestCase):

    def test_xml_2212L_tblastn_001(self):
        xml_file = get_file('xml_2212L_tblastn_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.12', qresult.version)
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(0.001, qresult.param_evalue_threshold)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('nr', qresult.target)
        
        # test parsed values of qresult
        self.assertEqual('gi|729325|sp|P39483|DHG2_BACME', qresult.id)
        self.assertEqual('Glucose 1-dehydrogenase II (GLCDH-II)', qresult.desc)
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
        self.assertEqual(4, hsp.query_from)
        self.assertEqual(249, hsp.query_to)
        self.assertEqual(15, hsp.hit_from)
        self.assertEqual(761, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(84, hsp.ident_num)
        self.assertEqual(143, hsp.pos_num)
        self.assertEqual(9, hsp.gap_num)
        self.assertEqual(252, hsp.ali_len)
        self.assertEqual(252, len(hsp))
        self.assertEqual('LKDKVVVVTGGSKGLGRAMAVRFGQEQSKVVVNYRSNXXXXXXXXXXXXXXXXGGQAIIVRGDVTKEEDVVNLVETAVKEFGSLDVMINNAGVENPVPSH---ELSLENWNQVIDTNLTGAFLGSREAIKYFVENDIKG-NVINMSSVHEMIPWPLFVHYAASKGGMKLMTETLALEYAPKGIRVNNIGPGAIDTPINAEKFADPEQRADVESMIPMGYIGKPEEIASVAAFLASSQASYVTGITLFADGGM', hsp.query.seq.tostring())
        self.assertEqual('LQGKVVAITGCSTGIGRAIAIGAAKNGANVVLHHLGDSTASDIAQVQEECKQAGAKTVVVPGDIAEAKTANEIVSAAVSSFSRIDVLISNAGI---CPFHSFLDLPHPLWKRVQDVNLNGSFYVVQAVANQMAKQEPKGGSIVAVSSISALMGGGEQCHYTPTKAGIKSLMESCAIALGPMGIRCNSVLPGTIETNINKEDLSNPEKRADQIRRVPLGRLGKPEDLVGPTLFFASDLSNYCTGASVLVDGGM', hsp.hit.seq.tostring())
        self.assertEqual('L+ KVV +TG S G+GRA+A+   +  + VV+++  +   +   + + E  +AG + ++V GD+ + +    +V  AV  F  +DV+I+NAG+    P H   +L    W +V D NL G+F   +       + + KG +++ +SS+  ++      HY  +K G+K + E+ A+   P GIR N++ PG I+T IN E  ++PE+RAD    +P+G +GKPE++     F AS  ++Y TG ++  DGGM', hsp.alignment_annotation['homology'])

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
        self.assertEqual(4, hsp.query_from)
        self.assertEqual(203, hsp.query_to)
        self.assertEqual(578, hsp.hit_from)
        self.assertEqual(1252, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(3, hsp.hit_frame)
        self.assertEqual(76, hsp.ident_num)
        self.assertEqual(113, hsp.pos_num)
        self.assertEqual(31, hsp.gap_num)
        self.assertEqual(228, hsp.ali_len)
        self.assertEqual(228, len(hsp))
        self.assertEqual('LKDKVVVVTGGSKGLGRAMAVRFGQEQSKVVVNYRSNXXXXXXXXXXXXXXGGQAIIVRGDVTKEEDVVNLVETAVKEFGSLDVMINNAGV-----------ENPVPS-HELSLE-----NWNQVIDTNLTGAFLGSREAIKYFVENDIKGNVINMSSVHEMIPW-PLFVHYAASKGGMKLMTETLALEYAPKGIRVNNIGPGAIDTPI----------NAEKFADPE', hsp.query.seq.tostring())
        self.assertEqual('LDGKVALVTGAGRGIGAAIAVALGQPGAKVVVNYANSREAAEKVVDEIKSNAQSAISIQADVGDPDAVTKLMDQAVEHFGYLDIVSSNAGIVSFGHVKDVTPDVCVPSPYESPVEL*PQQEFDRVFRVNTRGQFFVAREAYRHLREG---GRIILTSSNTASVKGVPRHAVYSGSKGAIDTFVRCLAIDCGDKKITVNAVAPGAIKTDMFLSVSREYIPNGETFTDEQ', hsp.hit.seq.tostring())
        self.assertEqual('L  KV +VTG  +G+G A+AV  GQ  +KVVVNY ++ E A +V  EI+     AI ++ DV   + V  L++ AV+ FG LD++ +NAG+           +  VPS +E  +E      +++V   N  G F  +REA ++  E    G +I  SS    +   P    Y+ SKG +      LA++   K I VN + PGAI T +          N E F D +', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_tblastn_001(self):
        xml_file = get_file('xml_2226_tblastn_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis '
                'subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(1205400.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(3, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|145479850|ref|XM_001425911.1|', hit.id)
        self.assertEqual('Paramecium tetraurelia hypothetical protein (GSPATT00004923001) partial mRNA', hit.desc)
        self.assertEqual(4632, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(34.6538, hsp.bitscore)
        self.assertEqual(78, hsp.bitscore_raw)
        self.assertEqual(1.08241e-05, hsp.evalue)
        self.assertEqual(30, hsp.query_from)
        self.assertEqual(72, hsp.query_to)
        self.assertEqual(1743, hsp.hit_from)
        self.assertEqual(1871, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(15, hsp.ident_num)
        self.assertEqual(26, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(43, hsp.ali_len)
        self.assertEqual(43, len(hsp))
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', hsp.query.seq.tostring())
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', hsp.hit.seq.tostring())
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(1119300.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|350596019|ref|XM_003360601.2|', hit.id)
        self.assertEqual('PREDICTED: Sus scrofa pleckstrin-like (LOC100626968), mRNA', hit.desc)
        self.assertEqual(772, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.57249e-67, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(94, hsp.hit_from)
        self.assertEqual(387, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(32.7278, hsp.bitscore)
        self.assertEqual(73, hsp.bitscore_raw)
        self.assertEqual(4.07518e-05, hsp.evalue)
        self.assertEqual(29, hsp.query_from)
        self.assertEqual(95, hsp.query_to)
        self.assertEqual(541, hsp.hit_from)
        self.assertEqual(753, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(21, hsp.ident_num)
        self.assertEqual(33, hsp.pos_num)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(71, hsp.ali_len)
        self.assertEqual(71, len(hsp))
        self.assertEqual('IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('+ +Y       P G I L+G  +TS      GK  F+ +     T  +  +F QAA  +ER  W++ I+ A', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|365982352|ref|XM_003667962.1|', hit.id)
        self.assertEqual('Naumovozyma dairenensis CBS 421 hypothetical protein (NDAI0A06120), mRNA', hit.desc)
        self.assertEqual(4932, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(19.631, hsp.bitscore)
        self.assertEqual(39, hsp.bitscore_raw)
        self.assertEqual(1.65923, hsp.evalue)
        self.assertEqual(11, hsp.query_from)
        self.assertEqual(53, hsp.query_to)
        self.assertEqual(3180, hsp.hit_from)
        self.assertEqual(3335, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(16, hsp.ident_num)
        self.assertEqual(23, hsp.pos_num)
        self.assertEqual(9, hsp.gap_num)
        self.assertEqual(52, hsp.ali_len)
        self.assertEqual(52, len(hsp))
        self.assertEqual('GSVFNTWKPMWVVLL---------EDGIEFYKKKSDNSPKGMIPLKGSTLTS', hsp.query.seq.tostring())
        self.assertEqual('GSCFPTWDLIFIEVLNPFLKEKLWEADNEEISKFVDLTLKGLVDLYPSHFTS', hsp.hit.seq.tostring())
        self.assertEqual('GS F TW  +++ +L         E   E   K  D + KG++ L  S  TS', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)

    def test_xml_2226_tblastn_002(self):
        xml_file = get_file('xml_2226_tblastn_002.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_tblastn_003(self):
        xml_file = get_file('xml_2226_tblastn_003.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(1205400.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(3, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|145479850|ref|XM_001425911.1|', hit.id)
        self.assertEqual('Paramecium tetraurelia hypothetical protein (GSPATT00004923001) partial mRNA', hit.desc)
        self.assertEqual(4632, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(34.6538, hsp.bitscore)
        self.assertEqual(78, hsp.bitscore_raw)
        self.assertEqual(1.08241e-05, hsp.evalue)
        self.assertEqual(30, hsp.query_from)
        self.assertEqual(72, hsp.query_to)
        self.assertEqual(1743, hsp.hit_from)
        self.assertEqual(1871, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(15, hsp.ident_num)
        self.assertEqual(26, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(43, hsp.ali_len)
        self.assertEqual(43, len(hsp))
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', hsp.query.seq.tostring())
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', hsp.hit.seq.tostring())
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.alignment_annotation['homology'])

        self.assertRaises(IndexError, hit.__getitem__, 1)

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_tblastn_004(self):
        xml_file = get_file('xml_2226_tblastn_004.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(1119300.0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|350596019|ref|XM_003360601.2|', hit.id)
        self.assertEqual('PREDICTED: Sus scrofa pleckstrin-like (LOC100626968), mRNA', hit.desc)
        self.assertEqual(772, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(199.519, hsp.bitscore)
        self.assertEqual(506, hsp.bitscore_raw)
        self.assertEqual(1.57249e-67, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(94, hsp.hit_from)
        self.assertEqual(387, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(32.7278, hsp.bitscore)
        self.assertEqual(73, hsp.bitscore_raw)
        self.assertEqual(4.07518e-05, hsp.evalue)
        self.assertEqual(29, hsp.query_from)
        self.assertEqual(95, hsp.query_to)
        self.assertEqual(541, hsp.hit_from)
        self.assertEqual(753, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(21, hsp.ident_num)
        self.assertEqual(33, hsp.pos_num)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(71, hsp.ali_len)
        self.assertEqual(71, len(hsp))
        self.assertEqual('IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('+ +Y       P G I L+G  +TS      GK  F+ +     T  +  +F QAA  +ER  W++ I+ A', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|365982352|ref|XM_003667962.1|', hit.id)
        self.assertEqual('Naumovozyma dairenensis CBS 421 hypothetical protein (NDAI0A06120), mRNA', hit.desc)
        self.assertEqual(4932, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(19.631, hsp.bitscore)
        self.assertEqual(39, hsp.bitscore_raw)
        self.assertEqual(1.65923, hsp.evalue)
        self.assertEqual(11, hsp.query_from)
        self.assertEqual(53, hsp.query_to)
        self.assertEqual(3180, hsp.hit_from)
        self.assertEqual(3335, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(16, hsp.ident_num)
        self.assertEqual(23, hsp.pos_num)
        self.assertEqual(9, hsp.gap_num)
        self.assertEqual(52, hsp.ali_len)
        self.assertEqual(52, len(hsp))
        self.assertEqual('GSVFNTWKPMWVVLL---------EDGIEFYKKKSDNSPKGMIPLKGSTLTS', hsp.query.seq.tostring())
        self.assertEqual('GSCFPTWDLIFIEVLNPFLKEKLWEADNEEISKFVDLTLKGLVDLYPSHFTS', hsp.hit.seq.tostring())
        self.assertEqual('GS F TW  +++ +L         E   E   K  D + KG++ L  S  TS', hsp.alignment_annotation['homology'])

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_tblastn_005(self):
        xml_file = get_file('xml_2226_tblastn_005.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('refseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(2933984, qresult.stat_db_num)
        self.assertEqual(4726730735, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]', qresult.desc)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual(2933984, qresult.stat_db_num)
        self.assertEqual(4726730735, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(3, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|145479850|ref|XM_001425911.1|', hit.id)
        self.assertEqual('Paramecium tetraurelia hypothetical protein (GSPATT00004923001) partial mRNA', hit.desc)
        self.assertEqual(4632, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(34.6538, hsp.bitscore)
        self.assertEqual(78, hsp.bitscore_raw)
        self.assertEqual(0.755176, hsp.evalue)
        self.assertEqual(30, hsp.query_from)
        self.assertEqual(72, hsp.query_to)
        self.assertEqual(1743, hsp.hit_from)
        self.assertEqual(1871, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(15, hsp.ident_num)
        self.assertEqual(26, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(43, hsp.ali_len)
        self.assertEqual(43, len(hsp))
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', hsp.query.seq.tostring())
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', hsp.hit.seq.tostring())
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.alignment_annotation['homology'])
        self.assertRaises(IndexError, hit.__getitem__, 1)

        # test parsed values of the third qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.desc)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual(2933984, qresult.stat_db_num)
        self.assertEqual(4726730735, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0.041, qresult.stat_kappa)
        self.assertEqual(0.267, qresult.stat_lambda)
        self.assertEqual(0.14, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|354480463|ref|XM_003502378.1|', hit.id)
        self.assertEqual('PREDICTED: Cricetulus griseus pleckstrin-like (LOC100773128), mRNA', hit.desc)
        self.assertEqual(1119, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(205.297, hsp.bitscore)
        self.assertEqual(521, hsp.bitscore_raw)
        self.assertEqual(1.46172e-63, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(75, hsp.hit_from)
        self.assertEqual(368, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(43.8986, hsp.bitscore)
        self.assertEqual(102, hsp.bitscore_raw)
        self.assertEqual(0.00054161, hsp.evalue)
        self.assertEqual(2, hsp.query_from)
        self.assertEqual(95, hsp.query_to)
        self.assertEqual(801, hsp.hit_from)
        self.assertEqual(1100, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(30, hsp.ident_num)
        self.assertEqual(50, hsp.pos_num)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(100, hsp.ali_len)
        self.assertEqual(100, len(hsp))
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDF-GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', hsp.query.seq.tostring())
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAIHLRGCVVTSVESNHDGKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA', hsp.hit.seq.tostring())
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +  GK+     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|390474391|ref|XM_002757683.2|', hit.id)
        self.assertEqual('PREDICTED: Callithrix jacchus pleckstrin (PLEK), mRNA', hit.desc)
        self.assertEqual(1402, hit.seq_len)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(202.986, hsp.bitscore)
        self.assertEqual(515, hsp.bitscore_raw)
        self.assertEqual(1.27031e-61, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual(160, hsp.hit_from)
        self.assertEqual(453, hsp.hit_to)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(96, hsp.ident_num)
        self.assertEqual(97, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ali_len)
        self.assertEqual(98, len(hsp))
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.query.seq.tostring())
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.hit.seq.tostring())
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(3, counter)


class TblastxXmlCases(unittest.TestCase):

    def test_xml_2212L_tblastx_001(self):
        xml_file = get_file('xml_2212L_tblastx_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0
        
        # test the first qresult
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.12', qresult.version)
        self.assertEqual(u'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch\xe4ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM80', qresult.param_matrix)
        self.assertEqual(1, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(10, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
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
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(354, hsp.query_to)
        self.assertEqual(44323, hsp.hit_from)
        self.assertEqual(44676, hsp.hit_to)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(-3, hsp.hit_frame)
        self.assertEqual(117, hsp.ident_num)
        self.assertEqual(117, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(118, hsp.ali_len)
        self.assertEqual(118, len(hsp))
        self.assertEqual('ECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD', hsp.query.seq.tostring())
        self.assertEqual('ECCFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD', hsp.hit.seq.tostring())
        self.assertEqual('EC FIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|4309961|gb|AC005993.2|AC005993', hit.id)
        self.assertEqual('Homo sapiens PAC clone RP6-114E22 from 14, complete sequence', hit.desc)
        self.assertEqual(143943, hit.seq_len)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(41.0922, hsp.bitscore)
        self.assertEqual(78, hsp.bitscore_raw)
        self.assertEqual(0.716571, hsp.evalue)
        self.assertEqual(166, hsp.query_from)
        self.assertEqual(249, hsp.query_to)
        self.assertEqual(43679, hsp.hit_from)
        self.assertEqual(43762, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(3, hsp.hit_frame)
        self.assertEqual(13, hsp.ident_num)
        self.assertEqual(19, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(28, hsp.ali_len)
        self.assertEqual(28, len(hsp))
        self.assertEqual('PLTKAHRLFQTSIVFYVTCFTASSQQLL', hsp.query.seq.tostring())
        self.assertEqual('PLNKYHTIFQISLCFYLFCYNMAQKQLL', hsp.hit.seq.tostring())
        self.assertEqual('PL K H +FQ S+ FY+ C+  + +QLL', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)
        
    def test_xml_2226_tblastx_001(self):
        xml_file = get_file('xml_2226_tblastx_001.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|296147483:1-350', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds', qresult.desc)
        self.assertEqual(350, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(1954618.0, qresult.stat_eff_space)
        self.assertEqual(0.133956144488482, qresult.stat_kappa)
        self.assertEqual(0.317605957635731, qresult.stat_lambda)
        self.assertEqual(0.401214524497119, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|296147483|ref|NM_001183135.1|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds', hit.desc)
        self.assertEqual(4911, hit.seq_len)
        self.assertEqual(8, len(hit))

        hsp = hit[0]
        self.assertEqual(289.739, hsp.bitscore)
        self.assertEqual(626, hsp.bitscore_raw)
        self.assertEqual(2.35531e-81, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(348, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(348, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(116, hsp.ident_num)
        self.assertEqual(116, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(116, hsp.ali_len)
        self.assertEqual(116, len(hsp))
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.hit.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.query.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(18.9375, hsp.bitscore)
        self.assertEqual(35, hsp.bitscore_raw)
        self.assertEqual(7.78658, hsp.evalue)
        self.assertEqual(292, hsp.query_from)
        self.assertEqual(324, hsp.query_to)
        self.assertEqual(340, hsp.hit_from)
        self.assertEqual(372, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(-3, hsp.hit_frame)
        self.assertEqual(6, hsp.ident_num)
        self.assertEqual(9, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, hsp.ali_len)
        self.assertEqual(11, len(hsp))
        self.assertEqual('KFWMPSLRLLI', hsp.query.seq.tostring())
        self.assertEqual('KKWVPPVKLLI', hsp.hit.seq.tostring())
        self.assertEqual('K W+P ++LLI', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|254579534|ref|XM_002495708.1|', hit.id)
        self.assertEqual('Zygosaccharomyces rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds', hit.desc)
        self.assertEqual(4866, hit.seq_len)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(141.279, hsp.bitscore)
        self.assertEqual(302, hsp.bitscore_raw)
        self.assertEqual(1.15566e-36, hsp.evalue)
        self.assertEqual(96, hsp.query_from)
        self.assertEqual(347, hsp.query_to)
        self.assertEqual(96, hsp.hit_from)
        self.assertEqual(347, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(57, hsp.ident_num)
        self.assertEqual(72, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(84, hsp.ali_len)
        self.assertEqual(84, len(hsp))
        self.assertEqual('IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', hsp.query.seq.tostring())
        self.assertEqual('IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK', hsp.hit.seq.tostring())
        self.assertEqual('IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(2, counter)

    def test_xml_2226_tblastx_002(self):
        xml_file = get_file('xml_2226_tblastx_002.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(-1, qresult.stat_kappa)
        self.assertEqual(-1, qresult.stat_lambda)
        self.assertEqual(-1, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_tblastx_003(self):
        xml_file = get_file('xml_2226_tblastx_003.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('db/minirefseq_mrna', qresult.target)

        self.assertEqual('gi|296147483:1-350', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds', qresult.desc)
        self.assertEqual(350, qresult.seq_len)
        self.assertEqual(23, qresult.stat_db_num)
        self.assertEqual(67750, qresult.stat_db_len)
        self.assertEqual(1954618.0, qresult.stat_eff_space)
        self.assertEqual(0.133956144488482, qresult.stat_kappa)
        self.assertEqual(0.317605957635731, qresult.stat_lambda)
        self.assertEqual(0.401214524497119, qresult.stat_entropy)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('gi|296147483|ref|NM_001183135.1|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds', hit.desc)
        self.assertEqual(4911, hit.seq_len)
        self.assertEqual(8, len(hit))

        hsp = hit[0]
        self.assertEqual(289.739, hsp.bitscore)
        self.assertEqual(626, hsp.bitscore_raw)
        self.assertEqual(2.35531e-81, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(348, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(348, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(116, hsp.ident_num)
        self.assertEqual(116, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(116, hsp.ali_len)
        self.assertEqual(116, len(hsp))
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.hit.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.query.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(18.9375, hsp.bitscore)
        self.assertEqual(35, hsp.bitscore_raw)
        self.assertEqual(7.78658, hsp.evalue)
        self.assertEqual(292, hsp.query_from)
        self.assertEqual(324, hsp.query_to)
        self.assertEqual(340, hsp.hit_from)
        self.assertEqual(372, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(-3, hsp.hit_frame)
        self.assertEqual(6, hsp.ident_num)
        self.assertEqual(9, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, hsp.ali_len)
        self.assertEqual(11, len(hsp))
        self.assertEqual('KFWMPSLRLLI', hsp.query.seq.tostring())
        self.assertEqual('KKWVPPVKLLI', hsp.hit.seq.tostring())
        self.assertEqual('K W+P ++LLI', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|254579534|ref|XM_002495708.1|', hit.id)
        self.assertEqual('Zygosaccharomyces rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds', hit.desc)
        self.assertEqual(4866, hit.seq_len)
        self.assertEqual(6, len(hit))

        hsp = hit[0]
        self.assertEqual(141.279, hsp.bitscore)
        self.assertEqual(302, hsp.bitscore_raw)
        self.assertEqual(1.15566e-36, hsp.evalue)
        self.assertEqual(96, hsp.query_from)
        self.assertEqual(347, hsp.query_to)
        self.assertEqual(96, hsp.hit_from)
        self.assertEqual(347, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(57, hsp.ident_num)
        self.assertEqual(72, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(84, hsp.ali_len)
        self.assertEqual(84, len(hsp))
        self.assertEqual('IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', hsp.query.seq.tostring())
        self.assertEqual('IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK', hsp.hit.seq.tostring())
        self.assertEqual('IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_xml_2226_tblastx_004(self):
        xml_file = get_file('xml_2226_tblastx_004.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test meta variables, only for the first one
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual('Stephen F. Altschul, Thomas L. Madden, Alejandro '
                'A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, '
                'and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a '
                'new generation of protein database search programs", '
                'Nucleic Acids Res. 25:3389-3402.', qresult.reference)
        self.assertEqual('BLOSUM62', qresult.param_matrix)
        self.assertEqual(10.0, qresult.param_evalue_threshold)
        self.assertEqual('L;', qresult.param_filter)
        self.assertEqual(11, qresult.param_gap_open)
        self.assertEqual(1, qresult.param_gap_extend)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('refseq_rna', qresult.target)

        # test parsed values of the first qresult
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.desc)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual(2933984, qresult.stat_db_num)
        self.assertEqual(4726730735, qresult.stat_db_len)
        self.assertEqual(0, qresult.stat_eff_space)
        self.assertEqual(0, qresult.stat_kappa)
        self.assertEqual(0, qresult.stat_lambda)
        self.assertEqual(0, qresult.stat_entropy)
        self.assertEqual(0, len(qresult))

        # test parsed values of the second qresult
        qresult = qresults.next()
        counter += 1
        self.assertEqual('gi|296147483:1-350', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds', qresult.desc)
        self.assertEqual(350, qresult.seq_len)
        self.assertEqual(2933984, qresult.stat_db_num)
        self.assertEqual(4726730735, qresult.stat_db_len)
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
        self.assertEqual(1.05874e-76, hsp.evalue)
        self.assertEqual(1, hsp.query_from)
        self.assertEqual(348, hsp.query_to)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(348, hsp.hit_to)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(116, hsp.ident_num)
        self.assertEqual(116, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(116, hsp.ali_len)
        self.assertEqual(116, len(hsp))
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.hit.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.query.seq.tostring())
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINPLRY*KEYIVLRS*KDIRISHSLSCWLANQGMLK*RPWQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.alignment_annotation['homology'])

        hsp = hit[-1]
        self.assertEqual(36.3494, hsp.bitscore)
        self.assertEqual(73, hsp.bitscore_raw)
        self.assertEqual(9.12288e-54, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(41, hsp.query_to)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(41, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(14, hsp.ident_num)
        self.assertEqual(14, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(14, hsp.ali_len)
        self.assertEqual(14, len(hsp))
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.query.seq.tostring())
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.hit.seq.tostring())
        self.assertEqual('MAMNTGGFDSMQRQ', hsp.alignment_annotation['homology'])

        hit = qresult[-1]
        self.assertEqual('gi|254579534|ref|XM_002495708.1|', hit.id)
        self.assertEqual('Zygosaccharomyces rouxii hypothetical protein '
                '(ZYRO0C02266g) mRNA, complete cds', hit.desc)
        self.assertEqual(4866, hit.seq_len)
        self.assertEqual(4, len(hit))

        hsp = hit[0]
        self.assertEqual(141.279, hsp.bitscore)
        self.assertEqual(302, hsp.bitscore_raw)
        self.assertEqual(5.19486e-32, hsp.evalue)
        self.assertEqual(96, hsp.query_from)
        self.assertEqual(347, hsp.query_to)
        self.assertEqual(96, hsp.hit_from)
        self.assertEqual(347, hsp.hit_to)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(57, hsp.ident_num)
        self.assertEqual(72, hsp.pos_num)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(84, hsp.ali_len)
        self.assertEqual(84, len(hsp))
        self.assertEqual('IRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRNAKMTTLAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', hsp.query.seq.tostring())
        self.assertEqual('IRNASDKSIEILKVVHSYEELSRHPDFIVPLVMSCASKNAKLTTISMQCFQKLATVPCIPVDKLSDVLDAFIEANQLAMDIKLK', hsp.hit.seq.tostring())
        self.assertEqual('IR+ASDKSIEILK VHS+EEL RHPDF +P V++C S+NAK+TT++MQC Q L+TVP IP  +LS++LDAFIEA  LAM+I+LK', hsp.alignment_annotation['homology'])

        # check if we've finished iteration over qresults
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(2, counter)

class BlastXmlSpecialCases(unittest.TestCase):

    def test_xml_2226_blastn_006(self):
        xml_file = get_file('xml_2226_blastn_006.xml')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test each qresult's attributes
        qresult = qresults.next()
        counter += 1

        # test the Hit IDs only, since this is a special case
        hit1 = qresult[0]
        hit2 = qresult[1]
        self.assertEqual('gi|347972582|ref|XM_309352.4|', hit1.id)
        self.assertEqual('Anopheles gambiae str. PEST AGAP011294-PA (DEFI_ANOGA) mRNA, complete cds', hit1.desc)
        self.assertEqual('gnl|BL_ORD_ID|17', hit2.id)
        self.assertEqual('gi|347972582|ref|XM_309352.4| Anopheles gambiae str. PEST AGAP011294-PA (DEFI_ANOGA) mRNA, complete cds', hit2.desc)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
