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


class CmsearchCases(unittest.TestCase):
    """Test parsing cmsearch output."""


    def test_cmsearch_1q_0m(self):
        """Test parsing infernal-text, cmsearch, single query, no hits"""
        text_file = get_file("IRES_Yeast.txt")
        qresults = parse(text_file, FMT)
        counter = itertools.count(start=1)

        qresult, counter  = next_result(qresults, counter)
        self.assertEqual(0, len(qresult))
        self.assertEqual(qresult.id, "IRES_HCV")
        self.assertEqual(qresult.seq_len, 352)
        self.assertEqual(qresult.accession, "RF00061")
        self.assertEqual(qresult.description, "Hepatitis C virus internal ribosome entry site")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)


    def test_cmsearch_1q_1m_1h_1f(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, two fragments"""
        tab_file = get_file("U2_Yeast-threshold.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, counter  = next_result(qresults, counter)
        self.assertEqual(1, len(qresult))
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.")
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(2, len(hsp))
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.evalue, 5.9e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertEqual(hsp.is_included, True)
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
        self.assertEqual(frag.query.seq, 
            'AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAGUuUAAuAuCUGauAuggcccccAuugggggccaau-uauaUUAaauuaAUUUUUggaacua'
            )
        self.assertEqual(frag.hit.seq,
            'AUC---UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGA-CCUCAAUGAGGCUCAUUaCCUUUUAAUUUG-------------'
        )
        self.assertEqual(frag.aln_annotation['PP'],
            '***...************************************************************999.****************86555555555443.............'
        )
        self.assertEqual(frag.aln_annotation['NC'], 
            '                                                                     v           v                               '
        )
        self.assertEqual(frag.aln_annotation['CS'],
            '::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>>>>>>,,,.,,,,,,,,,,,,,,,,,,,,,,,,,'
        )
        self.assertEqual(frag.aln_annotation['similarity'],
            'AU+   UCU+:GCCUUUUGGC:+AGAUCAAGUGUAGUAUCUGUUCUU:UCAGU+UAA+A+CUGA:AUG: CC:CA+UG:GG+:CA+U   U+UUAA+UU              '
        )
        # second fragment
        frag = hsp[1]
        self.assertEqual(frag.query_start, 186)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681754)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(frag.query.seq, 'Acccuuu')
        self.assertEqual(frag.hit.seq, 'ACAUUUU')
        self.assertEqual(frag.aln_annotation['PP'], '*******')
        self.assertEqual(frag.aln_annotation['NC'], '       ')
        self.assertEqual(frag.aln_annotation['CS'], ':::::::')
        self.assertEqual(frag.aln_annotation['similarity'], 'AC +UUU')

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)


    def test_cmsearch_1q_mm_1h_mf(self):
        """Test parsing infernal-text, cmsearch, one queries, multiple hits, one hsp, multiple fragments"""
        tab_file = get_file("U2_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, counter  = next_result(qresults, counter)
        self.assertEqual(5, len(qresult))
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
        self.assertEqual(1, len(hit))
        self.assertEqual(hit.id, "ENA|BK006948|BK006948.2")
        self.assertEqual(hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome XV, complete sequence.")
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(3, len(hsp))
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.39)
        self.assertEqual(hsp.evalue, 0.49)
        self.assertEqual(hsp.bitscore, 19.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertEqual(hsp.is_included, False)
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
        self.assertEqual(frag.query.seq, 'AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAG')
        self.assertEqual(frag.hit.seq, 'AUCCCAUAUUUGCCAUC-GGCAUAUAUUAAGUAUAUUAGCAGUUCUAAUUAC')
        self.assertEqual(frag.aln_annotation['PP'], '**************999.*******************************996')
        self.assertEqual(frag.aln_annotation['NC'], '                                                    ')
        self.assertEqual(frag.aln_annotation['CS'], '::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<')
        self.assertEqual(frag.aln_annotation['similarity'], 'AU+CC U U+ GCC U  GGC +A AU AAGU UA UA C GUUCU:A::A ')
        # second fragment
        frag = hsp[1]
        self.assertEqual(frag.query_start, 60)
        self.assertEqual(frag.query_end, 112)
        self.assertEqual(frag.hit_start, 737334)
        self.assertEqual(frag.hit_end, 737360)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(frag.query.seq, 'CUGauAuggcccccAuugggggccaauuauaUUAaauuaAUUUUUggaacua')
        self.assertEqual(frag.hit.seq, 'GUAGUUGGAAGGAUACUAUCCUUUAU--------------------------')
        self.assertEqual(frag.aln_annotation['PP'], '69999999999999999999999987..........................')
        self.assertEqual(frag.aln_annotation['NC'], '                                                    ')
        self.assertEqual(frag.aln_annotation['CS'], '>>>>>>,<<<<<<<___>>>>>>>,,,,,,,,,,,,,,,,,,,,,,,,,,,,')
        self.assertEqual(frag.aln_annotation['similarity'], ' U::U:  ::::::A U:::::: A+                          ')
        # third fragment
        frag = hsp[2]
        self.assertEqual(frag.query_start, 186)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 737324)
        self.assertEqual(frag.hit_end, 737331)
        self.assertEqual(frag.hit_strand, -1)
        self.assertEqual(frag.query.seq, 'Acccuuu')
        self.assertEqual(frag.hit.seq, 'AUCCCCU')
        self.assertEqual(frag.aln_annotation['PP'], '*******')
        self.assertEqual(frag.aln_annotation['NC'], '       ')
        self.assertEqual(frag.aln_annotation['CS'], ':::::::')
        self.assertEqual(frag.aln_annotation['similarity'], 'A CC++U')
        # third hit (2 hsp, forward strand)
        hit = qresult[2]
        self.assertEqual(1, len(hit))
        self.assertEqual(hit.id, "ENA|BK006947|BK006947.3")
        self.assertEqual(hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV, complete sequence.")
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(2, len(hsp))
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.39)
        self.assertEqual(hsp.evalue, 5.7)
        self.assertEqual(hsp.bitscore, 15.3)
        self.assertEqual(hsp.bias, 0.0)
        self.assertEqual(hsp.is_included, False)
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
        self.assertEqual(frag.query.seq, 'AUacCU-UCu-cgGCcUUUUgGCuaaGAUCAA-GUGUAG')
        self.assertEqual(frag.hit.seq, 'AUGUUGaUCUaUCGUCAAUUGACCCAGAUGAUaGUGUAG')
        self.assertEqual(frag.aln_annotation['PP'], '*****9****999****************9988999987')
        self.assertEqual(frag.aln_annotation['NC'], '            v          v               ')
        self.assertEqual(frag.aln_annotation['CS'], '::::::.<<<.-<<<<____>>>>->>>,,,,.,,,,,,')
        self.assertEqual(frag.aln_annotation['similarity'], 'AU     UCU + G C  UUG C  AGAU A  GUGUAG')
        # second fragment
        frag = hsp[1]
        self.assertEqual(frag.query_start, 84)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 266098)
        self.assertEqual(frag.hit_end, 266208)
        self.assertEqual(frag.hit_strand, 0)
        self.assertEqual(frag.query.seq, 'aauuauaUUAaauuaAUUUUUggaacuaGugggggcauuu-uggGCUUGCccau--ugcccccaCacggguugaccuggcaUUGCAcUaccgccagguucagcccAcccuuu')
        self.assertEqual(frag.hit.seq, '-GUUAUAUAGUUUUGAUAUUUUGGCGAAAAGUUGAGAAUAuUGCGCUUGCGUAUauAUUCCAUUUGAGGUGGCACUAGAGCUCGCAUUAU-UACCAGUAGUGGCAGGAUUGC')
        self.assertEqual(frag.aln_annotation['PP'], '.3377778888888888888888888888888888888888**************99999******************************.99*******************')
        self.assertEqual(frag.aln_annotation['NC'], '                                v                           v      v  v     v v             v v      v  v       ')
        self.assertEqual(frag.aln_annotation['CS'], ',,,,,,,,,,,,,,,,,,,,,,,,,,,,<<<<<<<<----.<<<<<__>>>>>-..->>>>>>>>,,<<<<<<-<<<<<<___________>>>>>>-->>>>>>:::::::')
        self.assertEqual(frag.aln_annotation['similarity'], '  UUAUAU   +UU AU UUU G   +A:: : G::A+U+ U :GCUUGC: AU  +::C : ::   G:  :AC: G   U GCA UA+   C :GU+:  :C    +U +')

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)


    def test_cmsearch_1q_mm_mh_1f(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, multiple hsp, one fragment"""
        tab_file = get_file("5S_Yeast.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, counter  = next_result(qresults, counter)
        self.assertEqual(1, len(qresult))
        self.assertEqual(qresult.id, "5S_rRNA")
        self.assertEqual(qresult.seq_len, 119)
        self.assertEqual(qresult.accession, "RF00001")
        self.assertEqual(qresult.description, "5S ribosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(6, len(hit))
        self.assertEqual(hit.id, "ENA|BK006945|BK006945.2")
        self.assertEqual(hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.")
        self.assertEqual(hit.query_id, "5S_rRNA")
        # first hit
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.52)
        self.assertEqual(hsp.evalue, 1.6e-18)
        self.assertEqual(hsp.bitscore, 88.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertEqual(hsp.is_included, True)
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
        self.assertEqual(frag.query.seq, 'gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc-gAAguUAAGcgcgcUugggCcagggUA-GUAcuagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc')
        self.assertEqual(frag.hit.seq, 'GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAGUGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAAUCU')
        self.assertEqual(frag.aln_annotation['PP'], '***********************************************99***********************8756***************************9*****************')
        self.assertEqual(frag.aln_annotation['NC'], '                                                                               vv                  vv                    ')
        self.assertEqual(frag.aln_annotation['CS'], '(((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-<<-----<<____>>----->>->-->>->>>))))))))):')
        self.assertEqual(frag.aln_annotation['similarity'], 'G::UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA::C+')

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)


    def test_cmsearch_1q_1m_1h_noali(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, noali"""
        tab_file = get_file("U2_Yeast-noali.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, counter  = next_result(qresults, counter)
        self.assertEqual(1, len(qresult))
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.4")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome II,")
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.evalue, 5.9e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertEqual(hsp.is_included, True)
        self.assertEqual(hsp.hit_start, 681747)
        self.assertEqual(hsp.hit_end, 681858)
        # first fragment
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681858)
        self.assertEqual(frag.hit_strand, -1)


    def test_cmsearch_1q_1m_mh_noali(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, multiple hsp, noali"""
        tab_file = get_file("5S_Yeast-noali.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)
        
        qresult, counter  = next_result(qresults, counter)
        self.assertEqual(1, len(qresult))
        self.assertEqual(qresult.id, "5S_rRNA")
        self.assertEqual(qresult.seq_len, 119)
        self.assertEqual(qresult.accession, "RF00001")
        self.assertEqual(qresult.description, "5S ribosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.5")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(6, len(hit))
        self.assertEqual(hit.id, "ENA|BK006945|BK006945.2")
        self.assertEqual(hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII")
        self.assertEqual(hit.query_id, "5S_rRNA")
        # first hsp
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.52)
        self.assertEqual(hsp.evalue, 1.6e-18)
        self.assertEqual(hsp.bitscore, 88.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertEqual(hsp.is_included, True)
        self.assertEqual(hsp.hit_start, 459676)
        self.assertEqual(hsp.hit_end, 459796)
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 459676)
        self.assertEqual(frag.hit_end, 459796)
        self.assertEqual(frag.hit_strand, 0)
        # last hsp
        hsp = hit[-1]
        self.assertEqual(1, len(hsp))
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.gc, 0.53)
        self.assertEqual(hsp.evalue, 4.4e-17)
        self.assertEqual(hsp.bitscore, 83.2)
        self.assertEqual(hsp.bias, 0.0)
        self.assertEqual(hsp.is_included, True)
        self.assertEqual(hsp.hit_start, 485697)
        self.assertEqual(hsp.hit_end, 485817)
        frag = hsp[0]
        self.assertEqual(frag.hit_start, 485697)
        self.assertEqual(frag.hit_end, 485817)
        self.assertEqual(frag.hit_strand, 0)


    def test_cmsearch_1q_1m_1h_1f_hmmonly(self):
        """Test parsing infernal-text, cmsearch, one queries, one hit, one hsp, one fragments, hmmonly"""
        tab_file = get_file("U2_Yeast-hmmonly.txt")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, counter  = next_result(qresults, counter)
        self.assertEqual(1, len(qresult))
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.seq_len, 193)
        self.assertEqual(qresult.accession, "RF00004")
        self.assertEqual(qresult.description, "U2 spliceosomal RNA")
        self.assertEqual(qresult.program, "cmsearch")
        self.assertEqual(qresult.version, "1.1.5")
        self.assertEqual(qresult.target, "GCA_000146045.2.fasta")
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(hit.description, "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.")
        self.assertEqual(hit.query_id, "U2")
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(hsp.model, "hmm")
        self.assertEqual(hsp.truncated, "-")
        self.assertEqual(hsp.gc, 0.35)
        self.assertEqual(hsp.evalue, 1.5e-19)
        self.assertEqual(hsp.bitscore, 73.1)
        self.assertEqual(hsp.bias, 2.7)
        self.assertEqual(hsp.is_included, True)
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
        self.assertEqual(frag.query.seq, 'UCucgGCcUUUUGGCUaaGAUCAAGUGUAGUAUCUGUUCUUuuCAGUuUAAuAuCUGauAuugucucuAuugggggccaau-uauaUUAaauuaA')
        self.assertEqual(frag.hit.seq,'UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGACCUCAAU-GAGGCUCAUUaCCUUUUAAUUUGU')
        self.assertEqual(frag.aln_annotation['PP'],'999*********************************************************7777666666.777776666544444444433221')
        self.assertEqual(frag.aln_annotation['CS'],'<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>>>>>>,,,.,,,,,,,,,,,,,')
        self.assertEqual(frag.aln_annotation['similarity'],'UCu+ GCcUUUUGGCU+aGAUCAAGUGUAGUAUCUGUUCUUuuCAGU+UAA+A+CUGa+Au+ +cuc Au g+gg  ca+u   u+UUAa+uu  ')

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)





if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
