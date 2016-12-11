# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO ExonerateIO parsers."""

import os
import unittest

from Bio import BiopythonExperimentalWarning

import warnings


with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio.SearchIO import parse, read

# test case files are in the Blast directory
TEST_DIR = 'Exonerate'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class ExonerateSpcCases(unittest.TestCase):

    coord = ('start', 'end')
    coords = ('inter_ranges', )
    stype = ('hit_', 'query_')

    def check_vulgar_text(self, vulgar, text):
        """Compares coordinate parsing for vulgar and text formats."""
        vfile = get_file(vulgar)
        tfile = get_file(text)

        vqres = read(vfile, 'exonerate-vulgar')
        tqres = read(tfile, 'exonerate-text')

        # compare coordinates of vulgar and text formats
        # should be the same since the files are results of the same query
        # vs db search
        for vhit, thit in zip(vqres, tqres):
            for vhsp, thsp in zip(vhit.hsps, thit.hsps):
                self.assertEqual(vhsp.query_start, thsp.query_start)
                self.assertEqual(vhsp.hit_start, thsp.hit_start)
                self.assertEqual(vhsp.query_end, thsp.query_end)
                self.assertEqual(vhsp.hit_end, thsp.hit_end)
                self.assertEqual(vhsp.query_inter_ranges, thsp.query_inter_ranges)
                self.assertEqual(vhsp.hit_inter_ranges, thsp.hit_inter_ranges)
                self.assertEqual(vhsp.query_split_codons, thsp.query_split_codons)
                self.assertEqual(vhsp.hit_split_codons, thsp.hit_split_codons)
                self.assertEqual(vhsp.query_frame_all, thsp.query_frame_all)
                self.assertEqual(vhsp.hit_frame_all, thsp.hit_frame_all)

    def test_vulgar_text_similar_g2g(self):
        """Compares vulgar-text coordinate parsing for the genome2genome model."""
        self.check_vulgar_text('exn_22_o_vulgar.exn', 'exn_22_m_genome2genome.exn')

    def test_vulgar_text_similar_c2c(self):
        """Compares vulgar-text coordinate parsing for the coding2coding model."""
        self.check_vulgar_text('exn_22_o_vulgar_fshifts.exn', 'exn_22_m_coding2coding_fshifts.exn')

    def test_vulgar_text_similar_p2d(self):
        """Compares vulgar-text coordinate parsing for the protein2dna model."""
        self.check_vulgar_text('exn_22_o_vulgar_fshifts2.exn', 'exn_22_m_protein2dna_fshifts.exn')


class ExonerateTextCases(unittest.TestCase):

    fmt = 'exonerate-text'

    def test_exn_22_m_affine_local(self):
        """Test parsing exonerate output (exn_22_m_affine_local.exn)"""

        exn_file = get_file('exn_22_m_affine_local.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('affine:local:dna2dna', qresult.model)
        self.assertEqual(3, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(6150, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443688|ref|NC_001145.3|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XIII, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(359, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(83, hsp.query_start)
        self.assertEqual(253990, hsp.hit_start)
        self.assertEqual(552, hsp.query_end)
        self.assertEqual(254474, hsp.hit_end)
        self.assertEqual([(83, 552)], hsp.query_range_all)
        self.assertEqual([(253990, 254474)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAA-AGGAT', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||| |||| |||||   ||| | |  ||| |||  | |||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ACCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCAAGAT', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('AAGTTATGTGGAACA--TAGGCTCATGGAACGCTCCCAGT', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('| ||   | | ||||  ||   |||   || | ||| |||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('ATGT--GGGGAAACAATTACCTTCACCAAATGATCCAAGT', str(hsp.hit_all[0].seq)[-40:])

        # third hit
        hit = qresult[2]
        self.assertEqual('gi|330443715|ref|NC_001146.8|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XIV, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(219, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(60, hsp.query_start)
        self.assertEqual(454073, hsp.hit_start)
        self.assertEqual(517, hsp.query_end)
        self.assertEqual(454531, hsp.hit_end)
        self.assertEqual([(60, 517)], hsp.query_range_all)
        self.assertEqual([(454073, 454531)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTG-', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||| || || || | |||||   |   |||| ||    | ', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGATGATATATTA-GATGGGG-ATG-AAGATGAGCCAGA', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('G-TATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('| |  |||||   | | | | |   | ||| | |||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('GTTGAAGAAGCCAAACGGAAGAAAGACGAGGA-GAGAAAG', str(hsp.hit_all[0].seq)[-40:])

    def test_exn_22_m_cdna2genome(self):
        """Test parsing exonerate output (exn_22_m_cdna2genome.exn)"""

        exn_file = get_file('exn_22_m_cdna2genome.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('cdna2genome', qresult.model)
        self.assertEqual(2, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(6146, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(6146, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGA', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGA', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('CTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('CTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT', str(hsp.hit_all[0].seq)[-40:])

        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443688|ref|NC_001145.3|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XIII, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(518, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(85010, hsp.hit_start)
        self.assertEqual(516, hsp.query_end)
        self.assertEqual(667216, hsp.hit_end)
        self.assertEqual([(0, 65), (65, 225), (225, 320), (320, 346), (346, 516)], hsp.query_range_all)
        self.assertEqual([(85010, 85066), (253974, 254135), (350959, 351052), (473170, 473201), (667040, 667216)], hsp.hit_range_all)
        self.assertEqual([(65, 65), (225, 225), (320, 320), (346, 346)], hsp.query_inter_ranges)
        self.assertEqual([(85066, 253974), (254135, 350959), (351052, 473170), (473201, 667040)], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(5, len(hsp.query_all))
        self.assertEqual(5, len(hsp.hit_all))
        # first block
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||  ||  | ||||   | ||||||  |||| | | | ', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGTGAACCT-CTTCAAGACGGTCAG--AATA-A-TCAA', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTT', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||  |||| | | | ||||    ||||||||||||| | |', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('AG--AATA-A-TCAACAGG----ATGAAGAAGCAAAAGAT', str(hsp.hit_all[0].seq)[-40:])
        # last block
        self.assertEqual('TATTAGCCTTCC--TCGATGATCTGCA--A-GAACAACAG', str(hsp.query_all[-1].seq)[:40])
        self.assertEqual('|  |||| || |  ||||| | || ||  | ||| | |  ', hsp[-1].aln_annotation['similarity'][:40])
        self.assertEqual('TCATAGCGTTACGTTCGAT-ACCTTCACTACGAAGATCCA', str(hsp.hit_all[-1].seq)[:40])
        self.assertEqual('AAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAA', str(hsp.query_all[-1].seq)[-40:])
        self.assertEqual('   |||||||||||||     ||  ||| | ||  | |||', hsp[-1].aln_annotation['similarity'][-40:])
        self.assertEqual('TTCTATAGAAGTACAGTTATTCAAACAAAAAAAAAAAAAA', str(hsp.hit_all[-1].seq)[-40:])

    def test_exn_22_m_coding2coding(self):
        """Test parsing exonerate output (exn_22_m_coding2coding.exn)"""

        exn_file = get_file('exn_22_m_coding2coding.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('coding2coding', qresult.model)
        self.assertEqual(2, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(2151, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(1318047, hsp.hit_start)
        self.assertEqual(1228, hsp.query_end)
        self.assertEqual(1319274, hsp.hit_end)
        self.assertEqual([(1, 1228)], hsp.query_range_all)
        self.assertEqual([(1318047, 1319274)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA', str(hsp.hit_all[0].seq)[-40:])

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(2106, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443688|ref|NC_001145.3|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XIII, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(116, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(1065, hsp.query_start)
        self.assertEqual(255638, hsp.hit_start)
        self.assertEqual(1224, hsp.query_end)
        self.assertEqual(255794, hsp.hit_end)
        self.assertEqual([(1065, 1224)], hsp.query_range_all)
        self.assertEqual([(255638, 255794)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('TGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAAT', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||+!  .!.|||   !!:...||+:!::!:!  ||+||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('TGTTCGGAAATTTGGGATAGAATAACAACACATCCGAAAT', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('CAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTC', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('!!!.||+...|||:!:||+   |||+||  !!::!!.:!:', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('CAATGCAGAAGACGTTCAATTAGCTTTGAATAAGCATATG', str(hsp.hit_all[0].seq)[-40:])

    def test_exn_22_m_coding2genome(self):
        """Test parsing exonerate output (exn_22_m_coding2genome.exn)"""

        exn_file = get_file('exn_22_m_coding2genome.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('coding2genome', qresult.model)
        self.assertEqual(2, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(2151, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(1318047, hsp.hit_start)
        self.assertEqual(1228, hsp.query_end)
        self.assertEqual(1319274, hsp.hit_end)
        self.assertEqual([(1, 1228)], hsp.query_range_all)
        self.assertEqual([(1318047, 1319274)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA', str(hsp.hit_all[0].seq)[-40:])

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(2106, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443688|ref|NC_001145.3|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XIII, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(116, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(1065, hsp.query_start)
        self.assertEqual(255638, hsp.hit_start)
        self.assertEqual(1224, hsp.query_end)
        self.assertEqual(255794, hsp.hit_end)
        self.assertEqual([(1065, 1224)], hsp.query_range_all)
        self.assertEqual([(255638, 255794)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('TGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAAT', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||+!  .!.|||   !!:...||+:!::!:!  ||+||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('TGTTCGGAAATTTGGGATAGAATAACAACACATCCGAAAT', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('CAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTC', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('!!!.||+...|||:!:||+   |||+||  !!::!!.:!:', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('CAATGCAGAAGACGTTCAATTAGCTTTGAATAAGCATATG', str(hsp.hit_all[0].seq)[-40:])

    def test_exn_22_m_dna2protein(self):
        """Test parsing exonerate output (exn_22_m_dna2protein.exn)"""

        exn_file = get_file('exn_22_m_dna2protein.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('dna', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('ungapped:dna2protein', qresult.model)
        self.assertEqual(1, len(qresult))

        hit = qresult[0]
        self.assertEqual('protein', hit.id)
        self.assertEqual('', hit.description)
        self.assertEqual(1, len(hit.hsps))

        hsp = hit[0]
        self.assertEqual(105, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(0, hsp[0].hit_strand)
        self.assertEqual(0, hsp[0].query_start)
        self.assertEqual(93, hsp[0].query_end)
        self.assertEqual(313, hsp[0].hit_start)
        self.assertEqual(344, hsp[0].hit_end)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('NSPFXKGPLASVQNPVYHKQPLNPAPNAETH', str(hsp[0].query.seq)[:40])
        self.assertEqual(['|||', '...', ' !!', ' !!', '! !'], hsp[0].aln_annotation['similarity'][:5])
        self.assertEqual('NQSVPKRPAGSVQNPVYHNQPLNPAPSRDPH', str(hsp[0].hit.seq)[:40])

    def test_exn_22_m_est2genome(self):
        """Test parsing exonerate output (exn_22_m_est2genome.exn)"""

        exn_file = get_file('exn_22_m_est2genome.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('est2genome', qresult.model)
        self.assertEqual(2, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(6150, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443688|ref|NC_001145.3|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XIII, complete sequence', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(439, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(85010, hsp.hit_start)
        self.assertEqual(346, hsp.query_end)
        self.assertEqual(473201, hsp.hit_end)
        self.assertEqual([(0, 65), (65, 225), (225, 320), (320, 346)], hsp.query_range_all)
        self.assertEqual([(85010, 85066), (253974, 254135), (350959, 351052), (473170, 473201)], hsp.hit_range_all)
        self.assertEqual([(65, 65), (225, 225), (320, 320)], hsp.query_inter_ranges)
        self.assertEqual([(85066, 253974), (254135, 350959), (351052, 473170)], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(4, len(hsp.query_all))
        self.assertEqual(4, len(hsp.hit_all))
        # first block
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||  ||  | ||||   | ||||||  |||| | | | ', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGTGAACCT-CTTCAAGACGGTCAG--AATA-A-TCAA', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTT', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||  |||| | | | ||||    ||||||||||||| | |', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('AG--AATA-A-TCAACAGG----ATGAAGAAGCAAAAGAT', str(hsp.hit_all[0].seq)[-40:])
        # last block
        self.assertEqual('AGCTAAGAATTCTGATGATG-----AAAGAA', str(hsp.query_all[-1].seq))
        self.assertEqual('|   |||||||||||| |||     ||||||', hsp[-1].aln_annotation['similarity'])
        self.assertEqual('ATGGAAGAATTCTGATAATGCTGTAAAAGAA', str(hsp.hit_all[-1].seq))

        # second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(263, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(25, hsp.query_start)
        self.assertEqual(11338, hsp.hit_start)
        self.assertEqual(406, hsp.query_end)
        self.assertEqual(130198, hsp.hit_end)
        self.assertEqual([(25, 183), (183, 252), (252, 406)], hsp.query_range_all)
        self.assertEqual([(130038, 130198), (120612, 120681), (11338, 11487)], hsp.hit_range_all)
        self.assertEqual([(183, 183), (252, 252)], hsp.query_inter_ranges)
        self.assertEqual([(120681, 130038), (11487, 120612)], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(3, len(hsp.query_all))
        self.assertEqual(3, len(hsp.hit_all))
        # first block
        self.assertEqual('AGCAAATATATTTA-GCAGGTGACATGAAGAAGCAAATGT', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('| |||| |||   | ||||   | | || |||| | |  |', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ACCAAAGATAACAAGGCAG--AAAAAGAGGAAGAAGAAAT', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('AG-GACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('|| |||  ||||||  ||   |||  || ||   ||| ||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('AGTGAC--CCCAGAGGAGCCAAGCAAAAAGA---TTCGGA', str(hsp.hit_all[0].seq)[-40:])
        # last block
        self.assertEqual('AATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAA', str(hsp.query_all[-1].seq)[:40])
        self.assertEqual('|||||||  | ||| |    |||| | |  | | ||    ', hsp[-1].aln_annotation['similarity'][:40])
        self.assertEqual('AATAAGAGCAACACAG----TTTA-TCTTATATGTA----', str(hsp.hit_all[-1].seq)[:40])
        self.assertEqual('CTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAA', str(hsp.query_all[-1].seq)[-40:])
        self.assertEqual('|  | | || |  | || ||  ||||||||  ||  ||||', hsp[-1].aln_annotation['similarity'][-40:])
        self.assertEqual('CCACTAAAAAATTATAAGAGCCAAAACGAAGTAGATACAA', str(hsp.hit_all[-1].seq)[-40:])

    def test_exn_22_m_genome2genome(self):
        """Test parsing exonerate output (exn_22_m_genome2genome.exn)"""

        exn_file = get_file('exn_22_m_genome2genome.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('sacCer3_dna', qresult.id)
        self.assertEqual('range=chrIV:1319469-1319997 5\'pad=0 3\'pad=0 strand=+ repeatMasking=none', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('genome2genome', qresult.model)
        self.assertEqual(3, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(2641, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1319468, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(1319997, hsp.hit_end)
        self.assertEqual([(0, 529)], hsp.query_range_all)
        self.assertEqual([(1319468, 1319997)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATCCCTTATCTCTTTATCTTGTTGCCTGGTTCTCTTTTCC', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATCCCTTATCTCTTTATCTTGTTGCCTGGTTCTCTTTTCC', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('ACGGCAATACCTGGCATGTGATTGTCGGAAAGAACTTTGG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('ACGGCAATACCTGGCATGTGATTGTCGGAAAGAACTTTGG', str(hsp.hit_all[0].seq)[-40:])

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(2641, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1319468, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(1319997, hsp.hit_end)
        self.assertEqual([(0, 529)], hsp.query_range_all)
        self.assertEqual([(1319468, 1319997)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('CCAAAGTTCTTTCCGACAATCACATGCCAGGTATTGCCGT', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('CCAAAGTTCTTTCCGACAATCACATGCCAGGTATTGCCGT', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('GGAAAAGAGAACCAGGCAACAAGATAAAGAGATAAGGGAT', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('GGAAAAGAGAACCAGGCAACAAGATAAAGAGATAAGGGAT', str(hsp.hit_all[0].seq)[-40:])

        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443489|ref|NC_001135.5|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome III, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(267, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(162, hsp.query_start)
        self.assertEqual(23668, hsp.hit_start)
        self.assertEqual(491, hsp.query_end)
        self.assertEqual(115569, hsp.hit_end)
        self.assertEqual([(462, 491), (413, 462), (378, 413), (302, 378), (162, 302)], hsp.query_range_all)
        self.assertEqual([(23668, 23697), (32680, 32732), (42287, 42325), (97748, 97821), (115419, 115569)], hsp.hit_range_all)
        self.assertEqual([(462, 462), (413, 413), (378, 378), (302, 302)], hsp.query_inter_ranges)
        self.assertEqual([(23697, 32680), (32732, 42287), (42325, 97748), (97821, 115419)], hsp.hit_inter_ranges)
        self.assertEqual([(378, 379), (376, 378)], hsp.query_split_codons)
        self.assertEqual([(42324, 42325), (97748, 97750)], hsp.hit_split_codons)
        self.assertEqual(5, len(hsp.query_all))
        self.assertEqual(5, len(hsp.hit_all))
        # first block
        self.assertEqual('CCCTTTAAATGGAGATTACAAACTAGCGA', str(hsp.query_all[0].seq))
        self.assertEqual('||  | ||| | |||  ||||| |  | |', hsp[0].aln_annotation['similarity'])
        self.assertEqual('CCGCTGAAAGGAAGAGAACAAAGTTACAA', str(hsp.hit_all[0].seq))
        # last block
        self.assertEqual('TTTTCTTTACTAAC-TCGAGGAAGAGTGAGGTTTTCTTCC', str(hsp.query_all[-1].seq)[:40])
        self.assertEqual('| ||    || | | |  |||||| |||| | | |  |||', hsp[-1].aln_annotation['similarity'][:40])
        self.assertEqual('TCTTGAAGACCAGCATGTAGGAAG-GTGATGATATGCTCC', str(hsp.hit_all[-1].seq)[:40])
        self.assertEqual('TTTGTGTGTGTACATTTGAATATATATATTTAC-TAACAA', str(hsp.query_all[-1].seq)[-40:])
        self.assertEqual(' |||  ||| |   |||||||||||||   | | ||||||', hsp[-1].aln_annotation['similarity'][-40:])
        self.assertEqual('ATTGATTGTTTTGTTTTGAATATATATTGATGCTTAACAA', str(hsp.hit_all[-1].seq)[-40:])

        # third hit
        hit = qresult[2]
        self.assertEqual('gi|330443667|ref|NC_001143.9|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XI, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(267, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(78, hsp.query_start)
        self.assertEqual(71883, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(641760, hsp.hit_end)
        self.assertEqual([(449, 529), (319, 388), (198, 284), (161, 198), (78, 114)], hsp.query_range_all)
        self.assertEqual([(641682, 641760), (487327, 487387), (386123, 386207), (208639, 208677), (71883, 71917)], hsp.hit_range_all)
        self.assertEqual([(388, 449), (284, 319), (198, 198), (114, 161)], hsp.query_inter_ranges)
        self.assertEqual([(487387, 641682), (386207, 487327), (208677, 386123), (71917, 208639)], hsp.hit_inter_ranges)
        self.assertEqual([(198, 200), (197, 198)], hsp.query_split_codons)
        self.assertEqual([(386123, 386125), (208676, 208677)], hsp.hit_split_codons)
        self.assertEqual(5, len(hsp.query_all))
        self.assertEqual(5, len(hsp.hit_all))
        # first block
        self.assertEqual('ATCCCTTATCTCTTTATCTTGTTGCCTGGTTCTCTTTTCC', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||       |||  |||||   ||||  ', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATCCCTTATCTCTTCTAAAGATTGTGTGGTT---TTTT--', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('AAATGGAGATTACAA---ACTAGCGAA-ACTGCAGAAAAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('  ||     || |||    || || ||  || || | |||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('GCATATTTTTTCCAACCTTCTTGCCAATTCTTCA-ACAAG', str(hsp.hit_all[0].seq)[-40:])
        # last block
        self.assertEqual('TAAAGATGCTCTGGACAAGTACCAGTTGGAAAGAGA', str(hsp.query_all[-1].seq))
        self.assertEqual(' ||||||  |||  || | |  ||||||||||||||', hsp[-1].aln_annotation['similarity'])
        self.assertEqual('AAAAGATTTTCT--ACGACTTGCAGTTGGAAAGAGA', str(hsp.hit_all[-1].seq))

    def test_exn_22_m_ungapped(self):
        """Test parsing exonerate output (exn_22_m_ungapped.exn)"""

        exn_file = get_file('exn_22_m_ungapped.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('ungapped:dna2dna', qresult.model)
        self.assertEqual(2, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(6150, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443688|ref|NC_001145.3|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XIII, complete sequence', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(233, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(121, hsp.query_start)
        self.assertEqual(254031, hsp.hit_start)
        self.assertEqual(236, hsp.query_end)
        self.assertEqual(254146, hsp.hit_end)
        self.assertEqual([(121, 236)], hsp.query_range_all)
        self.assertEqual([(254031, 254146)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('TTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('| ||  ||||| |||||   | |||||||||||| ||| |', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('TGGATCCTGAAACTAAGCAGAAGAGGACTGCCCAAAATCG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('CCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('  || ||||||   ||| |  ||| |||| |     ||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('GGAAGATGAAGGAATTGGAGAAGAAGGTACAAAGTTTAGA', str(hsp.hit_all[0].seq)[-40:])

        # second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(151, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(1098, hsp.query_start)
        self.assertEqual(255671, hsp.hit_start)
        self.assertEqual(1166, hsp.query_end)
        self.assertEqual(255739, hsp.hit_end)
        self.assertEqual([(1098, 1166)], hsp.query_range_all)
        self.assertEqual([(255671, 255739)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('CCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('|| ||||| |||    | ||  | |||| ||||||   ||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('CCGAAATACTCAGATATTGATGTCGATGGTTTATGTTCCG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('ATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual(' ||||||   |||  ||||    |||||||||||| ||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('GTTTATGTTCCGAGCTAATGGCAAAGGCAAAATGTTCAGA', str(hsp.hit_all[0].seq)[-40:])

    def test_exn_22_m_ungapped_trans(self):
        """Test parsing exonerate output (exn_22_m_ungapped_trans.exn)"""

        exn_file = get_file('exn_22_m_ungapped_trans.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('ungapped:codon', qresult.model)
        self.assertEqual(1, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(3, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(2151, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(1318047, hsp.hit_start)
        self.assertEqual(1228, hsp.query_end)
        self.assertEqual(1319274, hsp.hit_end)
        self.assertEqual([(1, 1228)], hsp.query_range_all)
        self.assertEqual([(1318047, 1319274)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGC', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('GCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA', str(hsp.hit_all[0].seq)[-40:])

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(2106, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # first hit, third hsp
        hsp = qresult[0].hsps[2]
        self.assertEqual(2072, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGA', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGA', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('CTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('CTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT', str(hsp.hit_all[0].seq)[-40:])

    def test_exn_22_m_ner(self):
        """Test parsing exonerate output (exn_22_m_ner.exn)"""

        exn_file = get_file('exn_22_m_ner.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('NER:affine:local:dna2dna', qresult.model)
        self.assertEqual(2, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(6150, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all[:5])
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all[:5])
        self.assertEqual([], hsp.query_inter_ranges[:5])
        self.assertEqual([], hsp.hit_inter_ranges[:5])
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(440, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(509, hsp.query_start)
        self.assertEqual(183946, hsp.hit_start)
        self.assertEqual(1192, hsp.query_end)
        self.assertEqual(184603, hsp.hit_end)
        self.assertEqual([(509, 514), (537, 547), (567, 595), (607, 617), (636, 650)], hsp.query_range_all[:5])
        self.assertEqual([(183946, 183951), (183977, 183987), (184002, 184030), (184044, 184054), (184066, 184080)], hsp.hit_range_all[:5])
        self.assertEqual([(514, 537), (547, 567), (595, 607), (617, 636), (650, 667)], hsp.query_inter_ranges[:5])
        self.assertEqual([(183951, 183977), (183987, 184002), (184030, 184044), (184054, 184066), (184080, 184092)], hsp.hit_inter_ranges[:5])
        self.assertEqual(24, len(hsp.query_all))
        self.assertEqual(24, len(hsp.hit_all))
        # first block
        self.assertEqual('TGAGA', str(hsp.query_all[0].seq))
        self.assertEqual('|||||', hsp[0].aln_annotation['similarity'])
        self.assertEqual('TGAGA', str(hsp.hit_all[0].seq))
        # last block
        self.assertEqual('GACTGCAAAATAGTAGTCAAAGCTC', str(hsp.query_all[-1].seq))
        self.assertEqual('||| | ||||||||||||| | |||', hsp[-1].aln_annotation['similarity'])
        self.assertEqual('GACGGTAAAATAGTAGTCACACCTC', str(hsp.hit_all[-1].seq))

        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443681|ref|NC_001144.5|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome XII, complete sequence', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(502, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(110, hsp.query_start)
        self.assertEqual(297910, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(318994, hsp.hit_end)
        self.assertEqual([(110, 117), (148, 159), (169, 182), (184, 197), (227, 244)], hsp.query_range_all[:5])
        self.assertEqual([(297910, 297917), (297946, 297957), (297970, 297983), (297992, 298004), (298019, 298038)], hsp.hit_range_all[:5])
        self.assertEqual([(117, 148), (159, 169), (182, 184), (197, 227), (244, 255)], hsp.query_inter_ranges[:5])
        self.assertEqual([(297917, 297946), (297957, 297970), (297983, 297992), (298004, 298019), (298038, 298049)], hsp.hit_inter_ranges[:5])
        self.assertEqual(33, len(hsp.query_all))
        self.assertEqual(33, len(hsp.hit_all))
        # first block
        self.assertEqual('CAGAAAA', str(hsp.query_all[0].seq))
        self.assertEqual('| |||||', hsp[0].aln_annotation['similarity'])
        self.assertEqual('CTGAAAA', str(hsp.hit_all[0].seq))
        # last block
        self.assertEqual('TGGTTAGACAGCTCCTGTAG', str(hsp.query_all[-1].seq))
        self.assertEqual('|| |  ||||||||||||||', hsp[-1].aln_annotation['similarity'])
        self.assertEqual('TGATAGGACAGCTCCTGTAG', str(hsp.hit_all[-1].seq))

    def test_exn_22_q_multiple(self):
        """Test parsing exonerate output (exn_22_q_multiple.exn)"""

        exn_file = get_file('exn_22_q_multiple.exn')
        qresults = list(parse(exn_file, self.fmt))
        self.assertEqual(2, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|296142823|ref|NM_001178508.1|', qresult.id)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('est2genome', qresult.model)
        self.assertEqual(3, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|330443482|ref|NC_001134.8|', hit.id)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, first hit, first hsp
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(4485, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(560077, hsp.hit_start)
        self.assertEqual(897, hsp.query_end)
        self.assertEqual(560974, hsp.hit_end)
        self.assertEqual([(0, 897)], hsp.query_range_all)
        self.assertEqual([(560077, 560974)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('CAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('CAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA', str(hsp.hit_all[0].seq)[-40:])

        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|330443753|ref|NC_001148.4|', hit.id)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, second hit, first hsp
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(941, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(492033, hsp.hit_start)
        self.assertEqual(896, hsp.query_end)
        self.assertEqual(492933, hsp.hit_end)
        self.assertEqual([(2, 896)], hsp.query_range_all)
        self.assertEqual([(492033, 492933)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTC', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||  |  || | || |||  |  ||  | || ||  | ', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('GAGC--TCTATGAACAGATTTAAGCAG-TTAGAAAAGCTT', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('C-AGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('| || ||      ||| ||||| |  ||   ||| |  ||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('CAAGCAGGCTCTGCAT-CACCCTTGGTTTGCAGAGTACTA', str(hsp.hit_all[0].seq)[-40:])

        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, third hit, first hsp
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(651, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(34, hsp.query_start)
        self.assertEqual(267809, hsp.hit_start)
        self.assertEqual(721, hsp.query_end)
        self.assertEqual(300717, hsp.hit_end)
        self.assertEqual([(34, 691), (691, 721)], hsp.query_range_all)
        self.assertEqual([(267809, 268448), (300686, 300717)], hsp.hit_range_all)
        self.assertEqual([(691, 691)], hsp.query_inter_ranges)
        self.assertEqual([(268448, 300686)], hsp.hit_inter_ranges)
        self.assertEqual(2, len(hsp.query_all))
        self.assertEqual(2, len(hsp.hit_all))
        # first block
        self.assertEqual('AGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGC', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||| ||||| ||||| || |  ||||||||    | ', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('AGAAAGTTGGTGAGGGTACTTATGCGGTTGTTTA-CTTGG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('CGATCAGATTTTCAAG--ATATTCAGAGTATTGGGAACGC', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('|||||| |  |  |||  |  ||||| |  || || || |', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('CGATCAAA--TGGAAGTAACGTTCAGGGCCTTAGGGACAC', str(hsp.hit_all[0].seq)[-40:])
        # last block
        self.assertEqual('CGAATGAAGCTA-TATGGCCAGATATTGTCT', str(hsp.query_all[-1].seq))
        self.assertEqual('| ||   || ||  ||||||||||||| |||', hsp[-1].aln_annotation['similarity'])
        self.assertEqual('CAAACCGAGATAGAATGGCCAGATATTCTCT', str(hsp.hit_all[-1].seq))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('est2genome', qresult.model)
        self.assertEqual(2, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual(1, len(hit.hsps))
        # second qresult, first hit, first hsp
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(6150, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318045, hsp.hit_start)
        self.assertEqual(1230, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual([(0, 1230)], hsp.query_range_all)
        self.assertEqual([(1318045, 1319275)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG', str(hsp.hit_all[0].seq)[-40:])

        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|330443688|ref|NC_001145.3|', hit.id)
        self.assertEqual(2, len(hit.hsps))
        # second qresult, second hit, first hsp
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(439, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(85010, hsp.hit_start)
        self.assertEqual(346, hsp.query_end)
        self.assertEqual(473201, hsp.hit_end)
        self.assertEqual([(0, 65), (65, 225), (225, 320), (320, 346)], hsp.query_range_all)
        self.assertEqual([(85010, 85066), (253974, 254135), (350959, 351052), (473170, 473201)], hsp.hit_range_all)
        self.assertEqual([(65, 65), (225, 225), (320, 320)], hsp.query_inter_ranges)
        self.assertEqual([(85066, 253974), (254135, 350959), (351052, 473170)], hsp.hit_inter_ranges)
        self.assertEqual(4, len(hsp.query_all))
        self.assertEqual(4, len(hsp.hit_all))
        # first block
        self.assertEqual('ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAG', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('||||  ||  | ||||   | ||||||  |||| | | | ', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ATGGTGAACCT-CTTCAAGACGGTCAG--AATA-A-TCAA', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTT', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('||  |||| | | | ||||    ||||||||||||| | |', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('AG--AATA-A-TCAACAGG----ATGAAGAAGCAAAAGAT', str(hsp.hit_all[0].seq)[-40:])
        # last block
        self.assertEqual('AGCTAAGAATTCTGATGATG-----AAAGAA', str(hsp.query_all[-1].seq))
        self.assertEqual('|   |||||||||||| |||     ||||||', hsp[-1].aln_annotation['similarity'])
        self.assertEqual('ATGGAAGAATTCTGATAATGCTGTAAAAGAA', str(hsp.hit_all[-1].seq))

        # second qresult, second hit, second hsp
        # second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(263, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(25, hsp.query_start)
        self.assertEqual(11338, hsp.hit_start)
        self.assertEqual(406, hsp.query_end)
        self.assertEqual(130198, hsp.hit_end)
        self.assertEqual([(25, 183), (183, 252), (252, 406)], hsp.query_range_all)
        self.assertEqual([(130038, 130198), (120612, 120681), (11338, 11487)], hsp.hit_range_all)
        self.assertEqual([(183, 183), (252, 252)], hsp.query_inter_ranges)
        self.assertEqual([(120681, 130038), (11487, 120612)], hsp.hit_inter_ranges)
        self.assertEqual(3, len(hsp.query_all))
        self.assertEqual(3, len(hsp.hit_all))
        # first block
        self.assertEqual('AGCAAATATATTTA-GCAGGTGACATGAAGAAGCAAATGT', str(hsp.query_all[0].seq)[:40])
        self.assertEqual('| |||| |||   | ||||   | | || |||| | |  |', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ACCAAAGATAACAAGGCAG--AAAAAGAGGAAGAAGAAAT', str(hsp.hit_all[0].seq)[:40])
        self.assertEqual('AG-GACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGA', str(hsp.query_all[0].seq)[-40:])
        self.assertEqual('|| |||  ||||||  ||   |||  || ||   ||| ||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('AGTGAC--CCCAGAGGAGCCAAGCAAAAAGA---TTCGGA', str(hsp.hit_all[0].seq)[-40:])
        # last block
        self.assertEqual('AATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAA', str(hsp.query_all[-1].seq)[:40])
        self.assertEqual('|||||||  | ||| |    |||| | |  | | ||    ', hsp[-1].aln_annotation['similarity'][:40])
        self.assertEqual('AATAAGAGCAACACAG----TTTA-TCTTATATGTA----', str(hsp.hit_all[-1].seq)[:40])
        self.assertEqual('CTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAA', str(hsp.query_all[-1].seq)[-40:])
        self.assertEqual('|  | | || |  | || ||  ||||||||  ||  ||||', hsp[-1].aln_annotation['similarity'][-40:])
        self.assertEqual('CCACTAAAAAATTATAAGAGCCAAAACGAAGTAGATACAA', str(hsp.hit_all[-1].seq)[-40:])

    def test_exn_22_m_coding2coding_fshifts(self):
        """Test parsing exonerate output (exn_22_m_coding2coding_fshifts.exn)"""

        exn_file = get_file('exn_22_m_coding2coding_fshifts.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('coding2coding', qresult.model)
        self.assertEqual(1, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|296143771|ref|NM_001180731.1|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', hit.description)
        self.assertEqual(2, len(hit))
        # first hit, first hsp
        hsp = qresult[0][0]
        self.assertEqual(213, hsp.score)
        self.assertEqual([1, 1, 1, 1], hsp.query_strand_all)
        self.assertEqual([1, 1, 1, 1], hsp.hit_strand_all)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(465, hsp.hit_start)
        self.assertEqual(160, hsp.query_end)
        self.assertEqual(630, hsp.hit_end)
        self.assertEqual([(0, 93), (94, 127), (127, 139), (139, 160)], hsp.query_range_all)
        self.assertEqual([(465, 558), (558, 591), (593, 605), (609, 630)], hsp.hit_range_all)
        self.assertEqual([(93, 94), (127, 127), (139, 139)], hsp.query_inter_ranges)
        self.assertEqual([(558, 558), (591, 593), (605, 609)], hsp.hit_inter_ranges)
        self.assertEqual([1, 2, 2, 2], hsp.query_frame_all)
        self.assertEqual([1, 1, 3, 1], hsp.hit_frame_all)
        self.assertEqual(4, len(hsp.query_all))
        self.assertEqual(4, len(hsp.hit_all))
        self.assertEqual(4, len(hsp.aln_annotation_all))
        # first block
        self.assertEqual('ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAG', str(hsp[0].query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAG', str(hsp[0].hit.seq)[:40])
        self.assertEqual('TATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACC', str(hsp[0].query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[0].aln_annotation['similarity'][-40:])
        self.assertEqual('TATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACC', str(hsp[0].hit.seq)[-40:])
        # last block
        self.assertEqual('GACGAAAGTATTAATGGTAGT', str(hsp[-1].query.seq))
        self.assertEqual('|||||||||||||||||||||', hsp[-1].aln_annotation['similarity'])
        self.assertEqual('GACGAAAGTATTAATGGTAGT', str(hsp[-1].hit.seq))

        # first hit, second hsp
        hsp = qresult[0][1]
        self.assertEqual(201, hsp.score)
        self.assertEqual([-1, -1], hsp.query_strand_all)
        self.assertEqual([-1, -1], hsp.hit_strand_all)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(466, hsp.hit_start)
        self.assertEqual(158, hsp.query_end)
        self.assertEqual(628, hsp.hit_end)
        self.assertEqual([(95, 158), (1, 94)], hsp.query_range_all)
        self.assertEqual([(559, 628), (466, 559)], hsp.hit_range_all)
        self.assertEqual([(94, 95)], hsp.query_inter_ranges)
        self.assertEqual([(559, 559)], hsp.hit_inter_ranges)
        self.assertEqual([-3, -2], hsp.query_frame_all)
        self.assertEqual([-2, -2], hsp.hit_frame_all)
        self.assertEqual(2, len(hsp.query_all))
        self.assertEqual(2, len(hsp.hit_all))
        self.assertEqual(2, len(hsp.aln_annotation_all))
        # first block
        self.assertEqual('TACCATTAATACTTTCGTCATGGT<-><->AACGGCATGT', str(hsp[0].query.seq)[:40])
        self.assertEqual('||||||||||||||||||||+ !       ...  !:!!|', hsp[0].aln_annotation['similarity'][:40])
        self.assertEqual('TACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGT', str(hsp[0].hit.seq)[:40])
        # last block
        self.assertEqual('TGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACAT', str(hsp[-1].query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp[-1].aln_annotation['similarity'][:40])
        self.assertEqual('TGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACAT', str(hsp[-1].hit.seq)[:40])

    def test_exn_22_m_protein2dna_fshifts(self):
        """Test parsing exonerate output (exn_22_m_protein2dna_fshifts.exn)"""

        exn_file = get_file('exn_22_m_protein2dna_fshifts.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('sp|P24813|YAP2_YEAST', qresult.id)
        self.assertEqual('AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('protein2dna:local', qresult.model)
        self.assertEqual(1, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|296143771|ref|NM_001180731.1|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds', hit.description)
        self.assertEqual(2, len(hit))
        # first hit, first hsp
        hsp = qresult[0][0]
        self.assertEqual(367, hsp.score)
        self.assertEqual([0, 0], hsp.query_strand_all)
        self.assertEqual([1, 1], hsp.hit_strand_all)
        self.assertEqual(330, hsp.query_start)
        self.assertEqual(216, hsp.hit_start)
        self.assertEqual(409, hsp.query_end)
        self.assertEqual(455, hsp.hit_end)
        self.assertEqual([(330, 373), (373, 409)], hsp.query_range_all)
        self.assertEqual([(216, 345), (347, 455)], hsp.hit_range_all)
        self.assertEqual([(373, 373)], hsp.query_inter_ranges)
        self.assertEqual([(345, 347)], hsp.hit_inter_ranges)
        self.assertEqual([0, 0], hsp.query_frame_all)
        self.assertEqual([1, 3], hsp.hit_frame_all)
        self.assertEqual(2, len(hsp.query_all))
        self.assertEqual(2, len(hsp.hit_all))
        self.assertEqual(2, len(hsp.aln_annotation_all))
        # first block
        self.assertEqual('HTKTIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYS', str(hsp[0].query.seq)[:40])
        self.assertEqual('HTKTIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYS', str(hsp[0].hit.seq)[:40])
        self.assertEqual('TIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYSSLD', str(hsp[0].query.seq)[-40:])
        self.assertEqual('TIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYSSLD', str(hsp[0].hit.seq)[-40:])
        # last block
        self.assertEqual('IDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL', str(hsp[-1].query.seq))
        self.assertEqual('IDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL', str(hsp[-1].hit.seq))

        # first hit, second hsp
        hsp = qresult[0][1]
        self.assertEqual(322, hsp.score)
        self.assertEqual([0], hsp.query_strand_all)
        self.assertEqual([1], hsp.hit_strand_all)
        self.assertEqual(6, hsp.query_start)
        self.assertEqual(16, hsp.hit_start)
        self.assertEqual(70, hsp.query_end)
        self.assertEqual(208, hsp.hit_end)
        self.assertEqual([(6, 70)], hsp.query_range_all)
        self.assertEqual([(16, 208)], hsp.hit_range_all)
        self.assertEqual([], hsp.query_inter_ranges)
        self.assertEqual([], hsp.hit_inter_ranges)
        self.assertEqual([0], hsp.query_frame_all)
        self.assertEqual([2], hsp.hit_frame_all)
        self.assertEqual(1, len(hsp.query_all))
        self.assertEqual(1, len(hsp.hit_all))
        self.assertEqual(1, len(hsp.aln_annotation_all))
        self.assertEqual('KGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAK', str(hsp[0].query.seq)[:40])
        self.assertEqual('KGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAK', str(hsp[0].hit.seq)[:40])
        self.assertEqual('RKVGRPGRKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKS', str(hsp[0].query.seq)[-40:])
        self.assertEqual('RKVGRPGRKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKS', str(hsp[0].hit.seq)[-40:])

    def test_exn_22_m_protein2genome(self):
        """Test parsing exonerate output (exn_22_m_protein2genome.exn)"""

        exn_file = get_file('exn_22_m_protein2genome.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('sp|P24813|YAP2_YEAST', qresult.id)
        self.assertEqual('AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('protein2genome:local', qresult.model)
        self.assertEqual(3, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome IV, complete sequence', hit.description)
        self.assertEqual(1, len(hit))
        # first hit, first hsp
        hsp = qresult[0][0]
        self.assertEqual(2105, hsp.score)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(-1, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1318048, hsp.hit_start)
        self.assertEqual(409, hsp.query_end)
        self.assertEqual(1319275, hsp.hit_end)
        self.assertEqual('MGNILRKGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKR', str(hsp[0].query.seq)[:40])
        self.assertEqual('MGNILRKGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKR', str(hsp[0].hit.seq)[:40])
        self.assertEqual('SSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL', str(hsp[0].query.seq)[-40:])
        self.assertEqual('SSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL', str(hsp[0].hit.seq)[-40:])

        # last hit
        hit = qresult[-1]
        self.assertEqual('gi|330443590|ref|NC_001140.6|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c chromosome VIII, complete sequence', hit.description)
        self.assertEqual(1, len(hit))
        # last hit, first hsp
        hsp = qresult[-1][0]
        self.assertEqual(122, hsp.score)
        self.assertEqual([0, 0], hsp.query_strand_all)
        self.assertEqual([-1, -1], hsp.hit_strand_all)
        self.assertEqual('RKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKSLQERX', str(hsp[0].query.seq))
        self.assertEqual('NENVPDDSKAKKKAQNRAAQKAFRERKEARMKELQDKX', str(hsp[0].hit.seq))
        self.assertEqual('!.!', hsp.aln_annotation_all[0]['similarity'][0])
        self.assertEqual(':!', hsp.aln_annotation_all[0]['similarity'][-1])
        self.assertEqual('AAT', hsp.aln_annotation_all[0]['hit_annotation'][0])
        self.assertEqual('TT', hsp.aln_annotation_all[0]['hit_annotation'][-1])
        self.assertEqual('XELLEQKDAQNKTTTDFLLCSLKSLLSEITKYRAKNSDDERILAFLDDLQE', str(hsp[-1].query.seq))
        self.assertEqual('XNKILNRDPQFMSNSSFHQCVSLDSINTIEKDEEKNSDDDAGLQAATDARE', str(hsp[-1].hit.seq))
        self.assertEqual('!', hsp.aln_annotation_all[-1]['similarity'][0])
        self.assertEqual('|||', hsp.aln_annotation_all[-1]['similarity'][-1])
        self.assertEqual('A', hsp.aln_annotation_all[-1]['hit_annotation'][0])
        self.assertEqual('GAA', hsp.aln_annotation_all[-1]['hit_annotation'][-1])

        self.assertEqual([(37, 74), (75, 125)], hsp.query_range_all)
        self.assertEqual([(84533, 84646), (68450, 68601)], hsp.hit_range_all)
        self.assertEqual([(74, 75)], hsp.query_inter_ranges)
        self.assertEqual([(68601, 84533)], hsp.hit_inter_ranges)
        self.assertEqual([0, 0], hsp.query_frame_all)
        self.assertEqual([-3, -3], hsp.hit_frame_all)
        self.assertEqual(2, len(hsp.query_all))
        self.assertEqual(2, len(hsp.hit_all))
        self.assertEqual(2, len(hsp.aln_annotation_all))

    def test_exn_22_q_none(self):
        """Test parsing exonerate output (exn_22_q_none.exn)"""
        exn_file = get_file('exn_22_q_none.exn')
        qresults = parse(exn_file, 'exonerate-text')
        self.assertRaises(StopIteration, next, qresults)


class ExonerateVulgarCases(unittest.TestCase):

    fmt = 'exonerate-vulgar'

    def test_exn_22_o_vulgar(self):
        """Test parsing exonerate output (exn_22_o_vulgar.exn)"""

        exn_file = get_file('exn_22_o_vulgar.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('sacCer3_dna', qresult.id)
        self.assertEqual('<unknown description>', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual(3, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(2641, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1319468, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(1319997, hsp.hit_end)
        self.assertEqual([(0, 529)], hsp.query_range_all[:5])
        self.assertEqual([(1319468, 1319997)], hsp.hit_range_all[:5])
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(' M 26 26 C 3 3 M 500 500', hsp.vulgar_comp)
        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(2641, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1319468, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(1319997, hsp.hit_end)
        self.assertEqual([(0, 529)], hsp.query_range_all[:5])
        self.assertEqual([(1319468, 1319997)], hsp.hit_range_all[:5])
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(' M 90 90 C 3 3 M 436 436', hsp.vulgar_comp)
        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443489|ref|NC_001135.5|', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(267, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(162, hsp.query_start)
        self.assertEqual(23668, hsp.hit_start)
        self.assertEqual(491, hsp.query_end)
        self.assertEqual(115569, hsp.hit_end)
        self.assertEqual([(462, 491), (413, 462), (378, 413), (302, 378), (162, 302)], hsp.query_range_all[:5])
        self.assertEqual([(23668, 23697), (32680, 32732), (42287, 42325), (97748, 97821), (115419, 115569)], hsp.hit_range_all[:5])
        self.assertEqual([(378, 379), (376, 378)], hsp.query_split_codons)
        self.assertEqual([(42324, 42325), (97748, 97750)], hsp.hit_split_codons)
        self.assertEqual(' M 29 29 5 0 2 I 0 8979 3 0 2 M 32 32 G 0 2 M 2 2 G 0 1 M 15 15 5 0 2 I 0 9551 3 0 2 M 3 3 G 1 0 M 5 5 G 0 2 M 3 3 G 0 1 M 4 4 G 0 1 M 18 18 S 1 1 5 0 2 I 0 55419 3 0 2 S 2 2 C 3 3 M 22 22 G 3 0 M 46 46 5 0 2 I 0 17594 3 0 2 M 14 14 G 0 1 M 9 9 G 1 0 M 15 15 G 0 3 M 17 17 G 0 3 M 1 1 G 0 1 M 13 13 G 0 1 M 6 6 G 1 0 M 12 12 G 0 2 M 45 45 G 0 1 M 6 6', hsp.vulgar_comp)
        # third hit
        hit = qresult[2]
        self.assertEqual('gi|330443667|ref|NC_001143.9|', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(267, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(78, hsp.query_start)
        self.assertEqual(71883, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(641760, hsp.hit_end)
        self.assertEqual([(449, 529), (319, 388), (198, 284), (161, 198), (78, 114)], hsp.query_range_all[:5])
        self.assertEqual([(641682, 641760), (487327, 487387), (386123, 386207), (208639, 208677), (71883, 71917)], hsp.hit_range_all[:5])
        self.assertEqual([(198, 200), (197, 198)], hsp.query_split_codons)
        self.assertEqual([(386123, 386125), (208676, 208677)], hsp.hit_split_codons)
        self.assertEqual(' M 31 31 G 3 0 M 4 4 G 2 0 M 19 19 G 0 3 M 9 9 G 0 1 M 6 6 G 1 0 M 5 5 5 2 2 I 0 154244 I 57 0 I 0 47 3 2 2 M 25 25 G 5 0 M 4 4 G 1 0 M 3 3 G 3 0 M 4 4 G 1 0 M 9 9 G 0 1 M 14 14 5 2 2 I 0 101116 I 31 0 3 2 2 M 23 23 G 0 1 M 15 15 G 1 0 M 9 9 G 1 0 M 2 2 G 1 0 M 14 14 C 18 18 S 2 2 5 0 2 I 0 177442 3 0 2 S 1 1 C 12 12 M 2 2 G 0 1 M 22 22 5 2 2 I 0 136697 I 7 0 I 0 6 I 1 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 2 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 3 0 I 0 1 I 2 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 2 0 I 0 2 I 2 0 I 0 2 I 17 0 3 2 2 M 12 12 G 2 0 M 22 22', hsp.vulgar_comp)

    def test_exn_22_o_vulgar_fshifts(self):
        """Test parsing exonerate output (exn_22_o_vulgar_fshifts.exn)"""

        exn_file = get_file('exn_22_o_vulgar_fshifts.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('gi|296143771|ref|NM_001180731.1|', qresult.id)
        self.assertEqual('<unknown description>', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual(1, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|296143771|ref|NM_001180731.1|', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertEqual(2, len(hit))
        # first hit, first hsp
        hsp = qresult[0][0]
        self.assertEqual(213, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(465, hsp.hit_start)
        self.assertEqual(160, hsp.query_end)
        self.assertEqual(630, hsp.hit_end)
        self.assertEqual([(0, 93), (94, 127), (127, 139), (139, 160)], hsp.query_range_all[:5])
        self.assertEqual([(465, 558), (558, 591), (593, 605), (609, 630)], hsp.hit_range_all[:5])
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(' C 93 93 F 1 0 C 33 33 F 0 2 C 12 12 F 0 4 C 21 21', hsp.vulgar_comp)
        # first hit, second hsp
        hsp = qresult[0][1]
        self.assertEqual(201, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(466, hsp.hit_start)
        self.assertEqual(158, hsp.query_end)
        self.assertEqual(628, hsp.hit_end)
        self.assertEqual([(95, 158), (1, 94)], hsp.query_range_all[:5])
        self.assertEqual([(559, 628), (466, 559)], hsp.hit_range_all[:5])
        self.assertEqual([], hsp.query_split_codons)
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(' C 24 24 G 0 6 C 39 39 F 1 0 C 93 93', hsp.vulgar_comp)


class ExonerateCigarCases(unittest.TestCase):

    fmt = 'exonerate-cigar'

    def test_exn_22_o_vulgar_cigar(self):
        """Test parsing exonerate output (exn_22_o_vulgar_cigar.exn)"""

        exn_file = get_file('exn_22_o_vulgar_cigar.exn')
        qresult = read(exn_file, self.fmt)

        # check common attributes
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('sacCer3_dna', qresult.id)
        self.assertEqual('<unknown description>', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual(3, len(qresult))
        # first hit
        hit = qresult[0]
        self.assertEqual('gi|330443520|ref|NC_001136.10|', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertEqual(2, len(hit.hsps))
        # first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(2641, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1319468, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(1319997, hsp.hit_end)
        self.assertEqual('  M 26 M 3 M 500', hsp.cigar_comp)
        # first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(2641, hsp.score)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(1319468, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(1319997, hsp.hit_end)
        self.assertEqual('  M 90 M 3 M 436', hsp.cigar_comp)
        # second hit
        hit = qresult[1]
        self.assertEqual('gi|330443489|ref|NC_001135.5|', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(267, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(1, hsp[0].hit_strand)
        self.assertEqual(162, hsp.query_start)
        self.assertEqual(23668, hsp.hit_start)
        self.assertEqual(491, hsp.query_end)
        self.assertEqual(115569, hsp.hit_end)
        self.assertEqual('  M 29 D 8983 M 32 D 2 M 2 D 1 M 15 D 9555 M 3 I 1 M 5 D 2 M 3 D 1 M 4 D 1 M 18 M 1 D 55423 M 5 M 22 I 3 M 46 D 17598 M 14 D 1 M 9 I 1 M 15 D 3 M 17 D 3 M 1 D 1 M 13 D 1 M 6 I 1 M 12 D 2 M 45 D 1 M 6', hsp.cigar_comp)
        # third hit
        hit = qresult[2]
        self.assertEqual('gi|330443667|ref|NC_001143.9|', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertEqual(1, len(hit.hsps))
        # third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(267, hsp.score)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(78, hsp.query_start)
        self.assertEqual(71883, hsp.hit_start)
        self.assertEqual(529, hsp.query_end)
        self.assertEqual(641760, hsp.hit_end)
        self.assertEqual('  M 31 I 3 M 4 I 2 M 19 D 3 M 9 D 1 M 6 I 1 M 7 D 154244 I 57 D 47 M 27 I 5 M 4 I 1 M 3 I 3 M 4 I 1 M 9 D 1 M 16 D 101116 I 31 M 25 D 1 M 15 I 1 M 9 I 1 M 2 I 1 M 14 M 20 D 177446 M 13 M 2 D 1 M 24 D 136697 I 7 D 6 I 1 D 1 I 1 D 1 I 1 D 1 I 1 D 1 I 2 D 1 I 1 D 1 I 1 D 1 I 3 D 1 I 2 D 1 I 1 D 1 I 1 D 1 I 2 D 2 I 2 D 2 I 17 M 14 I 2 M 22', hsp.cigar_comp)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
