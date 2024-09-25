# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for SearchIO InfernalIO infernal-text parser."""

import os
import unittest
import itertools

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "Infernal"
FMT = "infernal-text"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


def next_result(qresults, counter):
    """Iterate over the results and counter."""
    return next(qresults), next(counter)


class CmscanCases(unittest.TestCase):
    """Test parsing cmsearch output."""

    def test_cmscan(self):
        """Test parsing infernal-text, cmscan, multiple queries"""
        tab_file = get_file("cmscan_115_IRES_5S_U2_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "ENA|BK006935|BK006935.2")
        self.assertEqual(qresult.seq_len, 230218)
        self.assertEqual(qresult.description, "<unknown description>")
        self.assertEqual(qresult.program, "cmscan")
        self.assertEqual(qresult.version, "1.1.5")
        self.assertEqual(qresult.target, "IRES_5S_U2.cm")
        # first hit
        hit = qresult[0]
        self.assertEqual(len(hit), 2)
        self.assertEqual(hit.id, "U2")
        self.assertEqual(hit.description, "U2 spliceosomal RNA")
        self.assertEqual(hit.query_id, "ENA|BK006935|BK006935.2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 2)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.44)
        self.assertEqual(hsp.evalue, 0.91)
        self.assertEqual(hsp.bitscore, 13.5)
        self.assertEqual(hsp.bias, 0.0)
        self.assertFalse(hsp.is_included)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.query_endtype, "[]")
        self.assertEqual(hsp.hit_start, 52929)
        self.assertEqual(hsp.hit_end, 53083)
        self.assertEqual(hsp.hit_endtype, "..")
        self.assertEqual(hsp.avg_acc, 0.80)
        # first fragment
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 46)
        self.assertEqual(frag.hit_start, 52929)
        self.assertEqual(frag.hit_end, 52974)
        self.assertEqual(frag.hit_strand, 0)
        # second fragment
        frag = hsp[1]
        self.assertEqual(frag.query_start, 84)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 52977)
        self.assertEqual(frag.hit_end, 53083)
        self.assertEqual(frag.hit_strand, 0)
        # second hit
        hsp = hit[1]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 196389)
        self.assertEqual(hsp.hit_end, 196571)
        # second qresult
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "ENA|BK006936|BK006936.2")
        # first hit
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "U2")
        self.assertEqual(hit.description, "U2 spliceosomal RNA")
        hsp = hit[0]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 681747)
        self.assertEqual(hsp.hit_end, 681858)
        # third (last) qresult
        qresult, count = next_result(qresults, counter)
        self.assertEqual(qresult.id, "ENA|BK006937|BK006937.2")

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(3, count)


class CmsearchCases(unittest.TestCase):
    """Test parsing cmsearch output."""

    def test_cmsearch_1q_0m(self):
        """Test parsing infernal-text, cmsearch, one query, no hits (IRES_Yeast)"""
        text_file = get_file("cmsearch_114_IRES_Yeast.txt")
        qresults = parse(text_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 0)
        self.assertEqual(qresult.id, "IRES_HCV")
        self.assertEqual(qresult.seq_len, 352)
        self.assertEqual(qresult.accession, "RF00061")
        self.assertEqual(
            qresult.description, "Hepatitis C virus internal ribosome entry site"
        )
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_1q_1m_1h(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, multiple fragments (U2_Yeast)"""
        tab_file = get_file("cmsearch_114_U2_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.",
        )
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 2)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.evalue, 5.9e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.query_endtype, "[]")
        self.assertEqual(hsp.hit_start, 681747)
        self.assertEqual(hsp.hit_end, 681858)
        self.assertEqual(hsp.hit_endtype, "..")
        self.assertEqual(hsp.avg_acc, 0.91)
        # first fragment
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 112)
        self.assertEqual(frag.hit_start, 681763)
        self.assertEqual(frag.hit_end, 681858)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(
            frag.query.seq,
            "AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAGUuUAAuAuCUGauAuggcccccAuugggggccaau-uauaUUAaauuaAUUUUUggaacua",
        )
        self.assertEqual(
            frag.hit.seq,
            "AUC---UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGA-CCUCAAUGAGGCUCAUUaCCUUUUAAUUUG-------------",
        )
        self.assertEqual(
            frag.aln_annotation["PP"],
            "***...************************************************************999.****************86555555555443.............",
        )
        self.assertEqual(
            frag.aln_annotation["NC"],
            "                                                                     v           v                               ",
        )
        self.assertEqual(
            frag.aln_annotation["CS"],
            "::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>>>>>>,,,.,,,,,,,,,,,,,,,,,,,,,,,,,",
        )
        self.assertEqual(
            frag.aln_annotation["similarity"],
            "AU+   UCU+:GCCUUUUGGC:+AGAUCAAGUGUAGUAUCUGUUCUU:UCAGU+UAA+A+CUGA:AUG: CC:CA+UG:GG+:CA+U   U+UUAA+UU              ",
        )
        # second fragment
        frag = hsp[1]
        self.assertEqual(frag.query_start, 186)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681754)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(frag.query.seq, "Acccuuu")
        self.assertEqual(frag.hit.seq, "ACAUUUU")
        self.assertEqual(frag.aln_annotation["PP"], "*******")
        self.assertEqual(frag.aln_annotation["NC"], "       ")
        self.assertEqual(frag.aln_annotation["CS"], ":::::::")
        self.assertEqual(frag.aln_annotation["similarity"], "AC +UUU")

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_1q_mm_1h(self):
        """Test parsing infernal-text, cmsearch, one queries, multiple hits, one hsp, multiple fragments (U2_Yeast_full)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_full.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 5)
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        # skip first hit (equivalent to test_cmsearch_1q_1m)
        # second hit (3 hsp, reverse strand)
        hit = qresult[1]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006948|BK006948.2")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome XV, complete sequence.",
        )
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 3)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.39)
        self.assertEqual(hsp.evalue, 0.49)
        self.assertEqual(hsp.bitscore, 19.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertFalse(hsp.is_included)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.query_endtype, "[]")
        self.assertEqual(hsp.hit_start, 737324)
        self.assertEqual(hsp.hit_end, 737498)
        self.assertEqual(hsp.hit_endtype, "..")
        self.assertEqual(hsp.avg_acc, 0.96)
        # first fragment
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 52)
        self.assertEqual(frag.hit_start, 737448)
        self.assertEqual(frag.hit_end, 737498)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(
            frag.query.seq, "AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAG"
        )
        self.assertEqual(
            frag.hit.seq, "AUCCCAUAUUUGCCAUC-GGCAUAUAUUAAGUAUAUUAGCAGUUCUAAUUAC"
        )
        self.assertEqual(
            frag.aln_annotation["PP"],
            "**************999.*******************************996",
        )
        self.assertEqual(
            frag.aln_annotation["NC"],
            "                                                    ",
        )
        self.assertEqual(
            frag.aln_annotation["CS"],
            "::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<",
        )
        self.assertEqual(
            frag.aln_annotation["similarity"],
            "AU+CC U U+ GCC U  GGC +A AU AAGU UA UA C GUUCU:A::A ",
        )
        # second fragment
        frag = hsp[1]
        self.assertEqual(frag.query_start, 60)
        self.assertEqual(frag.query_end, 112)
        self.assertEqual(frag.hit_start, 737334)
        self.assertEqual(frag.hit_end, 737360)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(
            frag.query.seq, "CUGauAuggcccccAuugggggccaauuauaUUAaauuaAUUUUUggaacua"
        )
        self.assertEqual(
            frag.hit.seq, "GUAGUUGGAAGGAUACUAUCCUUUAU--------------------------"
        )
        self.assertEqual(
            frag.aln_annotation["PP"],
            "69999999999999999999999987..........................",
        )
        self.assertEqual(
            frag.aln_annotation["NC"],
            "                                                    ",
        )
        self.assertEqual(
            frag.aln_annotation["CS"],
            ">>>>>>,<<<<<<<___>>>>>>>,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
        )
        self.assertEqual(
            frag.aln_annotation["similarity"],
            " U::U:  ::::::A U:::::: A+                          ",
        )
        # third fragment
        frag = hsp[2]
        self.assertEqual(frag.query_start, 186)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 737324)
        self.assertEqual(frag.hit_end, 737331)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(frag.query.seq, "Acccuuu")
        self.assertEqual(frag.hit.seq, "AUCCCCU")
        self.assertEqual(frag.aln_annotation["PP"], "*******")
        self.assertEqual(frag.aln_annotation["NC"], "       ")
        self.assertEqual(frag.aln_annotation["CS"], ":::::::")
        self.assertEqual(frag.aln_annotation["similarity"], "A CC++U")
        # third hit (2 hsp, forward strand)
        hit = qresult[2]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006947|BK006947.3")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV, complete sequence.",
        )
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 2)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.39)
        self.assertEqual(hsp.evalue, 5.7)
        self.assertEqual(hsp.bitscore, 15.3)
        self.assertEqual(hsp.bias, 0.0)
        self.assertFalse(hsp.is_included)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.query_endtype, "[]")
        self.assertEqual(hsp.hit_start, 266059)
        self.assertEqual(hsp.hit_end, 266208)
        self.assertEqual(hsp.hit_endtype, "..")
        self.assertEqual(hsp.avg_acc, 0.91)
        # first fragment
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 36)
        self.assertEqual(frag.hit_start, 266059)
        self.assertEqual(frag.hit_end, 266097)
        self.assertEqual(frag.hit_strand, 0)
        self.assertEqual(frag.query.seq, "AUacCU-UCu-cgGCcUUUUgGCuaaGAUCAA-GUGUAG")
        self.assertEqual(frag.hit.seq, "AUGUUGaUCUaUCGUCAAUUGACCCAGAUGAUaGUGUAG")
        self.assertEqual(
            frag.aln_annotation["PP"], "*****9****999****************9988999987"
        )
        self.assertEqual(
            frag.aln_annotation["NC"], "            v          v               "
        )
        self.assertEqual(
            frag.aln_annotation["CS"], "::::::.<<<.-<<<<____>>>>->>>,,,,.,,,,,,"
        )
        self.assertEqual(
            frag.aln_annotation["similarity"], "AU     UCU + G C  UUG C  AGAU A  GUGUAG"
        )
        # second fragment
        frag = hsp[1]
        self.assertEqual(frag.query_start, 84)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 266098)
        self.assertEqual(frag.hit_end, 266208)
        self.assertEqual(frag.hit_strand, 0)
        self.assertEqual(
            frag.query.seq,
            "aauuauaUUAaauuaAUUUUUggaacuaGugggggcauuu-uggGCUUGCccau--ugcccccaCacggguugaccuggcaUUGCAcUaccgccagguucagcccAcccuuu",
        )
        self.assertEqual(
            frag.hit.seq,
            "-GUUAUAUAGUUUUGAUAUUUUGGCGAAAAGUUGAGAAUAuUGCGCUUGCGUAUauAUUCCAUUUGAGGUGGCACUAGAGCUCGCAUUAU-UACCAGUAGUGGCAGGAUUGC",
        )
        self.assertEqual(
            frag.aln_annotation["PP"],
            ".3377778888888888888888888888888888888888**************99999******************************.99*******************",
        )
        self.assertEqual(
            frag.aln_annotation["NC"],
            "                                v                           v      v  v     v v             v v      v  v       ",
        )
        self.assertEqual(
            frag.aln_annotation["CS"],
            ",,,,,,,,,,,,,,,,,,,,,,,,,,,,<<<<<<<<----.<<<<<__>>>>>-..->>>>>>>>,,<<<<<<-<<<<<<___________>>>>>>-->>>>>>:::::::",
        )
        self.assertEqual(
            frag.aln_annotation["similarity"],
            "  UUAUAU   +UU AU UUU G   +A:: : G::A+U+ U :GCUUGC: AU  +::C : ::   G:  :AC: G   U GCA UA+   C :GU+:  :C    +U +",
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_1q_mm_1h_shuffled(self):
        """Test parsing infernal-text, cmsearch, one queries, multiple non-consecutive hits, one hsp, multiple fragments (U2_Yeast_full_shuffled)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_full_shuffled.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 2)
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "BK006936_7-8.fasta")
        # first hit (3 hsps at rank 1,3 and 4)
        hit = qresult[0]
        self.assertEqual(len(hit), 3)
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(hit.description, "")
        self.assertEqual(hit.query_id, "U2")
        # first hsp (rank 1)
        hsp = hit[0]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 681747)
        self.assertEqual(hsp.hit_end, 681858)
        # second hsp (rank 3)
        hsp = hit[1]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 1370418)
        self.assertEqual(hsp.hit_end, 1370563)
        # last hsp (rank 4)
        hsp = hit[2]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 1079243)
        self.assertEqual(hsp.hit_end, 1079392)
        # second hit
        hit = qresult[1]
        self.assertEqual(len(hit), 3)
        self.assertEqual(hit.id, "ENA|BK006948|BK006948.2")
        self.assertEqual(hit.description, "")
        self.assertEqual(hit.query_id, "U2")
        # first hsp (rank 2)
        hsp = hit[0]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 737324)
        self.assertEqual(hsp.hit_end, 737498)
        # second hsp (rank 5)
        hsp = hit[1]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 425490)
        self.assertEqual(hsp.hit_end, 425693)
        # last hsp (rank 6)
        hsp = hit[2]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 1073786)
        self.assertEqual(hsp.hit_end, 1073950)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_1q_mm_mh(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, multiple hsp, one fragment (5S_Yeast)"""
        tab_file = get_file("cmsearch_114_5S_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "5S_rRNA")
        self.assertEqual(qresult.seq_len, 119)
        self.assertEqual(qresult.accession, "RF00001")
        self.assertEqual(qresult.description, "5S ribosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(len(hit), 6)
        self.assertEqual(hit.id, "ENA|BK006945|BK006945.2")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.",
        )
        self.assertEqual(hit.query_id, "5S_rRNA")
        # first hit
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.52)
        self.assertEqual(hsp.evalue, 1.6e-18)
        self.assertEqual(hsp.bitscore, 88.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 119)
        self.assertEqual(hsp.query_endtype, "[]")
        self.assertEqual(hsp.hit_start, 459676)
        self.assertEqual(hsp.hit_end, 459796)
        self.assertEqual(hsp.hit_endtype, "..")
        self.assertEqual(hsp.avg_acc, 0.99)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 119)
        self.assertEqual(frag.hit_start, 459676)
        self.assertEqual(frag.hit_end, 459796)
        self.assertEqual(frag.hit_strand, 0)
        self.assertEqual(
            frag.query.seq,
            "gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc-gAAguUAAGcgcgcUugggCcagggUA-GUAcuagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc",
        )
        self.assertEqual(
            frag.hit.seq,
            "GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAGUGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAAUCU",
        )
        self.assertEqual(
            frag.aln_annotation["PP"],
            "***********************************************99***********************8756***************************9*****************",
        )
        self.assertEqual(
            frag.aln_annotation["NC"],
            "                                                                               vv                  vv                    ",
        )
        self.assertEqual(
            frag.aln_annotation["CS"],
            "(((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-<<-----<<____>>----->>->-->>->>>))))))))):",
        )
        self.assertEqual(
            frag.aln_annotation["similarity"],
            "G::UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA::C+",
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_1q_1m_1h_noali(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, noali (U2_Yeast_noali)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_noali.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(
            hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome II,"
        )
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.evalue, 5.9e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.hit_start, 681747)
        self.assertEqual(hsp.hit_end, 681858)
        # first fragment
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681858)
        self.assertEqual(frag.hit_strand, -1)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_1q_1m_mh_noali(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, multiple hsp, noali (5S_Yeast_noali)"""
        tab_file = get_file("cmsearch_114_5S_Yeast_noali.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "5S_rRNA")
        self.assertEqual(qresult.seq_len, 119)
        self.assertEqual(qresult.accession, "RF00001")
        self.assertEqual(qresult.description, "5S ribosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(len(hit), 6)
        self.assertEqual(hit.id, "ENA|BK006945|BK006945.2")
        self.assertEqual(
            hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII"
        )
        self.assertEqual(hit.query_id, "5S_rRNA")
        # first hsp
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.52)
        self.assertEqual(hsp.evalue, 1.6e-18)
        self.assertEqual(hsp.bitscore, 88.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.hit_start, 459676)
        self.assertEqual(hsp.hit_end, 459796)
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 459676)
        self.assertEqual(frag.hit_end, 459796)
        self.assertEqual(frag.hit_strand, 0)
        # last hsp
        hsp = hit[-1]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.53)
        self.assertEqual(hsp.evalue, 4.4e-17)
        self.assertEqual(hsp.bitscore, 83.2)
        self.assertEqual(hsp.bias, 0.0)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.hit_start, 485697)
        self.assertEqual(hsp.hit_end, 485817)
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 485697)
        self.assertEqual(frag.hit_end, 485817)
        self.assertEqual(frag.hit_strand, 0)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_1q_1m_mh_noali_inc(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, multiple hsp, noali, inclusion threshold (U2_Yeast_full_noali)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_full_noali.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 5)
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        # first hit
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(
            hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome II,"
        )
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.evalue, 5.9e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.hit_start, 681747)
        self.assertEqual(hsp.hit_end, 681858)
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681858)
        self.assertEqual(frag.hit_strand, -1)
        # second hit
        hit = qresult[1]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006948|BK006948.2")
        self.assertEqual(
            hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome XV,"
        )
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.39)
        self.assertEqual(hsp.evalue, 0.49)
        self.assertEqual(hsp.bitscore, 19.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertFalse(hsp.is_included)
        self.assertEqual(hsp.hit_start, 737324)
        self.assertEqual(hsp.hit_end, 737498)
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 737324)
        self.assertEqual(frag.hit_end, 737498)
        self.assertEqual(frag.hit_strand, -1)
        # last hit
        hit = qresult[-1]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006939|BK006939.2")
        self.assertEqual(
            hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome V,"
        )
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.41)
        self.assertEqual(hsp.evalue, 7.1)
        self.assertEqual(hsp.bitscore, 14.9)
        self.assertEqual(hsp.bias, 0.0)
        self.assertFalse(hsp.is_included)
        self.assertEqual(hsp.hit_start, 190882)
        self.assertEqual(hsp.hit_end, 191043)
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 190882)
        self.assertEqual(frag.hit_end, 191043)
        self.assertEqual(frag.hit_strand, 0)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_1q_1m_1h_hmmonly(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, one fragments, hmmonly (U2_Yeast_hmmonly)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_hmmonly.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.",
        )
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "hmm")
        self.assertEqual(hsp.truncated, "-")
        self.assertEqual(hsp.gc, 0.35)
        self.assertEqual(hsp.evalue, 1.5e-19)
        self.assertEqual(hsp.bitscore, 73.1)
        self.assertEqual(hsp.bias, 2.7)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.query_start, 7)
        self.assertEqual(hsp.query_end, 100)
        self.assertEqual(hsp.query_endtype, "..")
        self.assertEqual(hsp.hit_start, 681762)
        self.assertEqual(hsp.hit_end, 681855)
        self.assertEqual(hsp.hit_endtype, "..")
        self.assertEqual(hsp.avg_acc, 0.74)
        # first fragment
        frag = hsp[0]
        self.assertEqual(frag.query_start, 7)
        self.assertEqual(frag.query_end, 100)
        self.assertEqual(frag.hit_start, 681762)
        self.assertEqual(frag.hit_end, 681855)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(
            frag.query.seq,
            "UCucgGCcUUUUGGCUaaGAUCAAGUGUAGUAUCUGUUCUUuuCAGUuUAAuAuCUGauAuugucucuAuugggggccaau-uauaUUAaauuaA",
        )
        self.assertEqual(
            frag.hit.seq,
            "UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGACCUCAAU-GAGGCUCAUUaCCUUUUAAUUUGU",
        )
        self.assertEqual(
            frag.aln_annotation["PP"],
            "999*********************************************************7777666666.777776666544444444433221",
        )
        self.assertEqual(
            frag.aln_annotation["CS"],
            "<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>>>>>>,,,.,,,,,,,,,,,,,",
        )
        self.assertEqual(
            frag.aln_annotation["similarity"],
            "UCu+ GCcUUUUGGCU+aGAUCAAGUGUAGUAUCUGUUCUUuuCAGU+UAA+A+CUGa+Au+ +cuc Au g+gg  ca+u   u+UUAa+uu  ",
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, count)

    def test_cmsearch_mq(self):
        """Test parsing infernal-text, cmsearch, multiple queries"""
        tab_file = get_file("cmsearch_114_IRES_5S_U2_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # First qresult (empty)
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 0)
        # Second qresult (5S, multiple hits)
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 3)
        self.assertEqual(qresult.id, "5S_rRNA")
        self.assertEqual(qresult.seq_len, 119)
        self.assertEqual(qresult.accession, "RF00001")
        self.assertEqual(qresult.description, "5S ribosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        # first hit
        hit = qresult[0]
        self.assertEqual(len(hit), 6)
        self.assertEqual(hit.id, "ENA|BK006945|BK006945.2")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.",
        )
        self.assertEqual(hit.query_id, "5S_rRNA")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.52)
        self.assertEqual(hsp.evalue, 1.6e-18)
        self.assertEqual(hsp.bitscore, 88.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 119)
        self.assertEqual(hsp.query_endtype, "[]")
        self.assertEqual(hsp.hit_start, 459676)
        self.assertEqual(hsp.hit_end, 459796)
        self.assertEqual(hsp.hit_endtype, "..")
        self.assertEqual(hsp.avg_acc, 0.99)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 119)
        self.assertEqual(frag.hit_start, 459676)
        self.assertEqual(frag.hit_end, 459796)
        self.assertEqual(frag.hit_strand, 0)
        self.assertEqual(
            frag.query.seq,
            "gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc-gAAguUAAGcgcgcUugggCcagggUA-GUAcuagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc",
        )
        self.assertEqual(
            frag.hit.seq,
            "GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAGUGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAAUCU",
        )
        self.assertEqual(
            frag.aln_annotation["PP"],
            "***********************************************99***********************8756***************************9*****************",
        )
        self.assertEqual(
            frag.aln_annotation["NC"],
            "                                                                               vv                  vv                    ",
        )
        self.assertEqual(
            frag.aln_annotation["CS"],
            "(((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-<<-----<<____>>----->>->-->>->>>))))))))):",
        )
        self.assertEqual(
            frag.aln_annotation["similarity"],
            "G::UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA::C+",
        )
        # last hit
        hsp = hit[-1]
        hit = qresult[-1]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006947|BK006947.3")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV, complete sequence.",
        )
        self.assertEqual(hit.query_id, "5S_rRNA")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.41)
        self.assertEqual(hsp.evalue, 6.6)
        self.assertEqual(hsp.bitscore, 16.7)
        self.assertEqual(hsp.bias, 0.3)
        self.assertFalse(hsp.is_included)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 119)
        self.assertEqual(hsp.query_endtype, "[]")
        self.assertEqual(hsp.hit_start, 6968)
        self.assertEqual(hsp.hit_end, 7085)
        self.assertEqual(hsp.hit_endtype, "..")
        self.assertEqual(hsp.avg_acc, 0.91)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 119)
        self.assertEqual(frag.hit_start, 6968)
        self.assertEqual(frag.hit_end, 7085)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(
            frag.query.seq,
            "gccuGcggcCAUAccagc-gcg-aAagcACcgGa-uCCCAUCcGaACuCcgAAguUAAGcgcgcUugggCcagggUAGUAcuagGaUGgGuGAcCuCcUGggAAgaccagGu-gccgCaggcc",
        )
        self.assertEqual(
            frag.hit.seq,
            "GAGAUGGUAUAUACUGUAgCAUcCGUGUACGUAUgACCGAUCAGA--AUACAAGUGAAGGUGAGUAUGGCAUGUG--GUAGUGGGAUUAGAG-UGGUAGGGUAAGUAUAUGUgUAUUAUUUAC",
        )
        self.assertEqual(
            frag.aln_annotation["PP"],
            "**************976325541459999****989999999999..89**********9999999*********..**********99988.689999************************",
        )
        self.assertEqual(
            frag.aln_annotation["NC"],
            "v             v  v                 v        v                 v   v                                                      v ",
        )
        self.assertEqual(
            frag.aln_annotation["CS"],
            "(((((((((,,,,<<-<<.<<<.---<<--<<<<.<<______>>-->>>>-->>---->>>>>-->><<<-<<----<-<<-----<<____>>----->>->-->>->>>.))))))))):",
        )
        self.assertEqual(
            frag.aln_annotation["similarity"],
            " : :: ::: AUAC +   ::     G:AC::::  CC AUC+G   ::::AA:U AAG ::  U+ GGC:  :G  GUA U+GGAU :G G :     GG AAG+: A:GU ::: :: : C",
        )
        # third qresult (U2, multiple hits)
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 5)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(3, count)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
