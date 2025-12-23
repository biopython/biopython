# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Tests for SearchIO InfernalIO infernal-text indexing"""

import os
import unittest

from search_tests_common import CheckIndex
from search_tests_common import CheckRaw


class InfernalTabRawCases(CheckRaw):
    fmt = "infernal-text"

    def test_infernal_text_1q(self):
        """Test infernal-text raw string retrieval, cmsearch, one query (U2_Yeast)."""
        filename = os.path.join("Infernal", "cmsearch_114_U2_Yeast.txt")
        raw = """# cmsearch :: search CM(s) against a sequence database
# INFERNAL 1.1.4 (Dec 2020)
# Copyright (C) 2020 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query CM file:                         RF00004.cm
# target sequence database:              GCA_000146045.2.fasta
# tabular output of hits:                U2_Yeast-threshold.tbl
# sequence reporting threshold:          score >= 46
# number of worker threads:              56
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       U2  [CLEN=193]
Accession:   RF00004
Description: U2 spliceosomal RNA
Hit scores:
 rank     E-value  score  bias  sequence                 start    end   mdl trunc   gc  description
 ----   --------- ------ -----  ----------------------- ------ ------   --- ----- ----  -----------
  (1) !   5.9e-20   98.7   0.1  ENA|BK006936|BK006936.2 681858 681747 -  cm    no 0.33  TPA_inf: Saccharomyces cerevisiae S288C chromosome II,


Hit alignments:
>> ENA|BK006936|BK006936.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (1) !   5.9e-20   98.7   0.1  cm        1      193 []      681858      681747 - .. 0.91    no 0.33

                                                                                                      v           NC
                                 ::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>> CS
                       U2      1 AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAGUuUAAuAuCUGauAuggcccccAuuggg 80    
                                 AU+   UCU+:GCCUUUUGGC:+AGAUCAAGUGUAGUAUCUGUUCUU:UCAGU+UAA+A+CUGA:AUG: CC:CA+UG:G
  ENA|BK006936|BK006936.2 681858 AUC---UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGA-CCUCAAUGAG 681783
                                 ***...************************************************************999.********** PP

                                  v                                                   NC
                                 >>>>,,,.,,,,,,,,,,,,,,,,,,,,,,,,,~~~~~~~~~~~~::::::: CS
                       U2     81 ggccaau.uauaUUAaauuaAUUUUUggaacua*[34]**[40]*Acccuuu 193   
                                 G+:CA+U   U+UUAA+UU                          AC +UUU
  ENA|BK006936|BK006936.2 681782 GCUCAUUaCCUUUUAAUUUG-------------*[ 6]**[ 3]*ACAUUUU 681747
                                 ******86555555555443................7.....9..******* PP



Internal CM pipeline statistics summary:
----------------------------------------
Query model(s):                                                  1  (193 consensus positions)
Target sequences:                                               16  (24142652 residues searched)
Target sequences re-searched for truncated hits:                16  (15424 residues re-searched)
Windows   passing  local HMM SSV           filter:           72732  (0.7978); expected (0.35)
Windows   passing  local HMM Viterbi       filter:           24175  (0.3233); expected (0.15)
Windows   passing  local HMM Viterbi  bias filter:            6671  (0.09669); expected (0.15)
Windows   passing  local HMM Forward       filter:            2037  (0.03161); expected (0.003)
Windows   passing  local HMM Forward  bias filter:            1133  (0.01757); expected (0.003)
Windows   passing glocal HMM Forward       filter:             596  (0.01251); expected (0.003)
Windows   passing glocal HMM Forward  bias filter:             438  (0.009175); expected (0.003)
Envelopes passing glocal HMM envelope defn filter:             460  (0.00429); expected (0.003)
Envelopes passing  local CM  CYK           filter:              38  (0.000201); expected (0.0001)
Total CM hits reported:                                          1  (4.636e-06); includes 0 truncated hit(s)

# CPU time: 64.67u 2.16s 00:01:06.83 Elapsed: 00:00:03.09
//
"""
        self.check_raw(filename, "U2", raw)

    def test_infernal_text_mq_first(self):
        """Test infernal-text raw string retrieval, cmsearch, multiple queries, first (IRES_5S_U2_Yeast)."""
        filename = os.path.join("Infernal", "cmsearch_114_IRES_5S_U2_Yeast.txt")
        raw = """# cmsearch :: search CM(s) against a sequence database
# INFERNAL 1.1.4 (Dec 2020)
# Copyright (C) 2020 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query CM file:                         IRES_5S_U2.cm
# target sequence database:              GCA_000146045.2.fasta
# tabular output of hits:                IRES_5S_U2_Yeast.tbl
# number of worker threads:              56
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       IRES_HCV  [CLEN=352]
Accession:   RF00061
Description: Hepatitis C virus internal ribosome entry site
Hit scores:
 rank     E-value  score  bias  sequence  start    end   mdl trunc   gc  description
 ----   --------- ------ -----  -------- ------ ------   --- ----- ----  -----------

   [No hits detected that satisfy reporting thresholds]


Hit alignments:

   [No hits detected that satisfy reporting thresholds]


Internal CM pipeline statistics summary:
----------------------------------------
Query model(s):                                                  1  (352 consensus positions)
Target sequences:                                               16  (24142652 residues searched)
Target sequences re-searched for truncated hits:                16  (28160 residues re-searched)
Windows   passing  local HMM SSV           filter:            6432  (0.1374); expected (0.35)
Windows   passing  local HMM Viterbi       filter:            1631  (0.03555); expected (0.15)
Windows   passing  local HMM Viterbi  bias filter:            1607  (0.03502); expected (0.15)
Windows   passing  local HMM Forward       filter:               9  (0.0001992); expected (0.003)
Windows   passing  local HMM Forward  bias filter:               9  (0.0001992); expected (0.003)
Windows   passing glocal HMM Forward       filter:               1  (2.16e-05); expected (0.003)
Windows   passing glocal HMM Forward  bias filter:               1  (2.16e-05); expected (0.003)
Envelopes passing glocal HMM envelope defn filter:               1  (1.705e-05); expected (0.003)
Envelopes passing  local CM  CYK           filter:               0  (0); expected (0.0001)
Total CM hits reported:                                          0  (0); includes 0 truncated hit(s)

# CPU time: 3.98u 0.19s 00:00:04.17 Elapsed: 00:00:00.87
//
"""
        self.check_raw(filename, "IRES_HCV", raw)

    def test_infernal_text_mq_middle(self):
        """Test infernal-text raw string retrieval, cmsearch, multiple queries, middle (IRES_5S_U2_Yeast)."""
        filename = os.path.join("Infernal", "cmsearch_114_IRES_5S_U2_Yeast.txt")
        raw = """# cmsearch :: search CM(s) against a sequence database
# INFERNAL 1.1.4 (Dec 2020)
# Copyright (C) 2020 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query CM file:                         IRES_5S_U2.cm
# target sequence database:              GCA_000146045.2.fasta
# tabular output of hits:                IRES_5S_U2_Yeast.tbl
# number of worker threads:              56
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       5S_rRNA  [CLEN=119]
Accession:   RF00001
Description: 5S ribosomal RNA
Hit scores:
 rank     E-value  score  bias  sequence                 start    end   mdl trunc   gc  description
 ----   --------- ------ -----  ----------------------- ------ ------   --- ----- ----  -----------
  (1) !   1.6e-18   88.8   0.0  ENA|BK006945|BK006945.2 459676 459796 +  cm    no 0.52  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII
  (2) !   1.6e-18   88.8   0.0  ENA|BK006945|BK006945.2 489349 489469 +  cm    no 0.52  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII
  (3) !   4.4e-17   83.2   0.0  ENA|BK006945|BK006945.2 468813 468933 +  cm    no 0.53  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII
  (4) !   4.4e-17   83.2   0.0  ENA|BK006945|BK006945.2 472465 472585 +  cm    no 0.53  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII
  (5) !   4.4e-17   83.2   0.0  ENA|BK006945|BK006945.2 482045 482165 +  cm    no 0.53  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII
  (6) !   4.4e-17   83.2   0.0  ENA|BK006945|BK006945.2 485697 485817 +  cm    no 0.53  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII
 ------ inclusion threshold ------
  (7) ?      0.56   20.9   0.0  ENA|BK006943|BK006943.2 357031 357144 +  cm    no 0.46  TPA_inf: Saccharomyces cerevisiae S288C chromosome X, 
  (8) ?       6.6   16.7   0.3  ENA|BK006947|BK006947.3   7085   6968 -  cm    no 0.41  TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV


Hit alignments:
>> ENA|BK006945|BK006945.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (1) !   1.6e-18   88.8   0.0  cm        1      119 []      459676      459796 + .. 0.99    no 0.52

                                                                                                                v NC
                                 (((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-< CS
                  5S_rRNA      1 gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc.gAAguUAAGcgcgcUugggCcagggUA.GUAc 78    
                                 G::UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA 
  ENA|BK006945|BK006945.2 459676 GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAG 459755
                                 ***********************************************99***********************8756**** PP

                                 v                  vv                     NC
                                 <-----<<____>>----->>->-->>->>>))))))))): CS
                  5S_rRNA     79 uagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc 119   
                                  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA::C+
  ENA|BK006945|BK006945.2 459756 UGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAAUCU 459796
                                 ***********************9***************** PP

>> ENA|BK006945|BK006945.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (2) !   1.6e-18   88.8   0.0  cm        1      119 []      489349      489469 + .. 0.99    no 0.52

                                                                                                                v NC
                                 (((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-< CS
                  5S_rRNA      1 gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc.gAAguUAAGcgcgcUugggCcagggUA.GUAc 78    
                                 G::UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA 
  ENA|BK006945|BK006945.2 489349 GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAG 489428
                                 ***********************************************99***********************8756**** PP

                                 v                  vv                     NC
                                 <-----<<____>>----->>->-->>->>>))))))))): CS
                  5S_rRNA     79 uagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc 119   
                                  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA::C+
  ENA|BK006945|BK006945.2 489429 UGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAAUCU 489469
                                 ***********************9***************** PP

>> ENA|BK006945|BK006945.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (3) !   4.4e-17   83.2   0.0  cm        1      119 []      468813      468933 + .. 0.99    no 0.53

                                                                                                                v NC
                                 (((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-< CS
                  5S_rRNA      1 gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc.gAAguUAAGcgcgcUugggCcagggUA.GUAc 78    
                                  : UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA 
  ENA|BK006945|BK006945.2 468813 GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAG 468892
                                 ***********************************************99***********************8756**** PP

                                 v                  vv                     NC
                                 <-----<<____>>----->>->-->>->>>))))))))): CS
                  5S_rRNA     79 uagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc 119   
                                  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA :  
  ENA|BK006945|BK006945.2 468893 UGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAGUUG 468933
                                 ***********************9***************** PP

>> ENA|BK006945|BK006945.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (4) !   4.4e-17   83.2   0.0  cm        1      119 []      472465      472585 + .. 0.99    no 0.53

                                                                                                                v NC
                                 (((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-< CS
                  5S_rRNA      1 gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc.gAAguUAAGcgcgcUugggCcagggUA.GUAc 78    
                                  : UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA 
  ENA|BK006945|BK006945.2 472465 GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAG 472544
                                 ***********************************************99***********************8756**** PP

                                 v                  vv                     NC
                                 <-----<<____>>----->>->-->>->>>))))))))): CS
                  5S_rRNA     79 uagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc 119   
                                  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA :  
  ENA|BK006945|BK006945.2 472545 UGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAGUUG 472585
                                 ***********************9***************** PP

>> ENA|BK006945|BK006945.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (5) !   4.4e-17   83.2   0.0  cm        1      119 []      482045      482165 + .. 0.99    no 0.53

                                                                                                                v NC
                                 (((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-< CS
                  5S_rRNA      1 gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc.gAAguUAAGcgcgcUugggCcagggUA.GUAc 78    
                                  : UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA 
  ENA|BK006945|BK006945.2 482045 GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAG 482124
                                 ***********************************************99***********************8756**** PP

                                 v                  vv                     NC
                                 <-----<<____>>----->>->-->>->>>))))))))): CS
                  5S_rRNA     79 uagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc 119   
                                  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA :  
  ENA|BK006945|BK006945.2 482125 UGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAGUUG 482165
                                 ***********************9***************** PP

>> ENA|BK006945|BK006945.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (6) !   4.4e-17   83.2   0.0  cm        1      119 []      485697      485817 + .. 0.99    no 0.53

                                                                                                                v NC
                                 (((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>.>-->>---->>>>>-->><<<-<<---.-<-< CS
                  5S_rRNA      1 gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCc.gAAguUAAGcgcgcUugggCcagggUA.GUAc 78    
                                  : UGC:GCCAUA:C :C::GAAAGCACCG :UCCC+UCCGA C: C G AGUUAAGC::G: +G:GCC G:    GUA 
  ENA|BK006945|BK006945.2 485697 GGUUGCGGCCAUAUCUACCAGAAAGCACCGUUUCCCGUCCGAUCAACuGUAGUUAAGCUGGUAAGAGCCUGACCGaGUAG 485776
                                 ***********************************************99***********************8756**** PP

                                 v                  vv                     NC
                                 <-----<<____>>----->>->-->>->>>))))))))): CS
                  5S_rRNA     79 uagGaUGgGuGAcCuCcUGggAAgaccagGugccgCaggcc 119   
                                  +  +UGGGUGACC+   G  AA  :CAGGUGC:GCA :  
  ENA|BK006945|BK006945.2 485777 UGUAGUGGGUGACCAUACGCGAAACUCAGGUGCUGCAGUUG 485817
                                 ***********************9***************** PP

>> ENA|BK006943|BK006943.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome X, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (7) ?      0.56   20.9   0.0  cm        1      119 []      357031      357144 + .. 0.86    no 0.46

                                         v        vv      v      vv      vv         v       vv     vv             NC
                                 (((((((((,,,,<<-<<<<<---<<--<<<<<<______>>-->>>>-.->>---->>>>>-->><<<-<<----<-<< CS
                  5S_rRNA      1 gccuGcggcCAUAccagcgcgaAagcACcgGauCCCAUCcGaACuCcgA.AguUAAGcgcgcUugggCcagggUAGUAcu 79    
                                 ::: GC:G CA A:: G  :G+AA: AC::  ++ CAU     C  ::A A :UAAGC:  CU+::  C G:G  GUACU
  ENA|BK006943|BK006943.2 357031 CAGGGCUGGCAGAGGUGUCGGGAAAAACAAGGAU-CAUAU--CCUUUUAcAAUUAAGCCAUCUACCACCUGAG--GUACU 357105
                                 ****************************888743.44433..555555579**********************..***** PP

                                       v    v                vvv           NC
                                 -----<<____>>----->>->-->>->>>)))))).))): CS
                  5S_rRNA     80 agGaUGgGuGAcCuCcUGggAAgaccagGugccgCa.ggcc 119   
                                 A  + G      CU C GGGAA+A:C+G   C:GC  :::+
  ENA|BK006943|BK006943.2 357106 AAAG-G-AAAGGCUACCGGGAAUAUCUGAAACAGCUgCUGU 357144
                                 9994.3.33334778889****************9879*** PP

>> ENA|BK006947|BK006947.3  TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (8) ?       6.6   16.7   0.3  cm        1      119 []        7085        6968 - .. 0.91    no 0.41

                               v             v  v                 v        v                 v   v                  NC
                               (((((((((,,,,<<-<<.<<<.---<<--<<<<.<<______>>-->>>>-->>---->>>>>-->><<<-<<----<-<<-- CS
                  5S_rRNA    1 gccuGcggcCAUAccagc.gcg.aAagcACcgGa.uCCCAUCcGaACuCcgAAguUAAGcgcgcUugggCcagggUAGUAcuag 81  
                                : :: ::: AUAC +   ::     G:AC::::  CC AUC+G   ::::AA:U AAG ::  U+ GGC:  :G  GUA U+G
  ENA|BK006947|BK006947.3 7085 GAGAUGGUAUAUACUGUAgCAUcCGUGUACGUAUgACCGAUCAGA--AUACAAGUGAAGGUGAGUAUGGCAUGUG--GUAGUGG 7006
                               **************976325541459999****989999999999..89**********9999999*********..******* PP

                                                                    v  NC
                               ---<<____>>----->>->-->>->>>.))))))))): CS
                  5S_rRNA   82 GaUGgGuGAcCuCcUGggAAgaccagGu.gccgCaggcc 119 
                               GAU :G G :     GG AAG+: A:GU ::: :: : C
  ENA|BK006947|BK006947.3 7005 GAUUAGAG-UGGUAGGGUAAGUAUAUGUgUAUUAUUUAC 6968
                               ***99988.689999************************ PP



Internal CM pipeline statistics summary:
----------------------------------------
Query model(s):                                                  1  (119 consensus positions)
Target sequences:                                               16  (24142652 residues searched)
Target sequences re-searched for truncated hits:                16  (12416 residues re-searched)
Windows   passing  local HMM SSV           filter:           24991  (0.2423); expected (0.35)
Windows   passing  local HMM Viterbi       filter:            8464  (0.08504); expected (0.15)
Windows   passing  local HMM Viterbi  bias filter:            8432  (0.08473); expected (0.15)
Windows   passing  local HMM Forward       filter:             135  (0.001502); expected (0.003)
Windows   passing  local HMM Forward  bias filter:             134  (0.001493); expected (0.003)
Windows   passing glocal HMM Forward       filter:              65  (0.0007404); expected (0.003)
Windows   passing glocal HMM Forward  bias filter:              65  (0.0007404); expected (0.003)
Envelopes passing glocal HMM envelope defn filter:              61  (0.0003446); expected (0.003)
Envelopes passing  local CM  CYK           filter:              15  (6.52e-05); expected (0.0001)
Total CM hits reported:                                          8  (3.966e-05); includes 0 truncated hit(s)

# CPU time: 4.17u 0.21s 00:00:04.38 Elapsed: 00:00:00.40
//
"""
        self.check_raw(filename, "5S_rRNA", raw)

    def test_infernal_text_mq_last(self):
        """Test infernal-text raw string retrieval, cmsearch, multiple queries, last (IRES_5S_U2_Yeast)."""
        filename = os.path.join("Infernal", "cmsearch_114_IRES_5S_U2_Yeast.txt")
        raw = """# cmsearch :: search CM(s) against a sequence database
# INFERNAL 1.1.4 (Dec 2020)
# Copyright (C) 2020 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query CM file:                         IRES_5S_U2.cm
# target sequence database:              GCA_000146045.2.fasta
# tabular output of hits:                IRES_5S_U2_Yeast.tbl
# number of worker threads:              56
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       U2  [CLEN=193]
Accession:   RF00004
Description: U2 spliceosomal RNA
Hit scores:
 rank     E-value  score  bias  sequence                 start    end   mdl trunc   gc  description
 ----   --------- ------ -----  ----------------------- ------ ------   --- ----- ----  -----------
  (1) !   5.9e-20   98.7   0.1  ENA|BK006936|BK006936.2 681858 681747 -  cm    no 0.33  TPA_inf: Saccharomyces cerevisiae S288C chromosome II,
 ------ inclusion threshold ------
  (2) ?      0.49   19.8   0.0  ENA|BK006948|BK006948.2 737498 737324 -  cm    no 0.39  TPA_inf: Saccharomyces cerevisiae S288C chromosome XV,
  (3) ?       5.7   15.3   0.0  ENA|BK006947|BK006947.3 266059 266208 +  cm    no 0.39  TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV
  (4) ?       6.6   15.1   0.4  ENA|BK006949|BK006949.2 443393 443253 -  cm    no 0.32  TPA_inf: Saccharomyces cerevisiae S288C chromosome XVI
  (5) ?       7.1   14.9   0.0  ENA|BK006939|BK006939.2 190882 191043 +  cm    no 0.41  TPA_inf: Saccharomyces cerevisiae S288C chromosome V, 


Hit alignments:
>> ENA|BK006936|BK006936.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (1) !   5.9e-20   98.7   0.1  cm        1      193 []      681858      681747 - .. 0.91    no 0.33

                                                                                                      v           NC
                                 ::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<________>>>>>>,<<<<<<<___>>> CS
                       U2      1 AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAGUuUAAuAuCUGauAuggcccccAuuggg 80    
                                 AU+   UCU+:GCCUUUUGGC:+AGAUCAAGUGUAGUAUCUGUUCUU:UCAGU+UAA+A+CUGA:AUG: CC:CA+UG:G
  ENA|BK006936|BK006936.2 681858 AUC---UCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUAACAACUGAAAUGA-CCUCAAUGAG 681783
                                 ***...************************************************************999.********** PP

                                  v                                                   NC
                                 >>>>,,,.,,,,,,,,,,,,,,,,,,,,,,,,,~~~~~~~~~~~~::::::: CS
                       U2     81 ggccaau.uauaUUAaauuaAUUUUUggaacua*[34]**[40]*Acccuuu 193   
                                 G+:CA+U   U+UUAA+UU                          AC +UUU
  ENA|BK006936|BK006936.2 681782 GCUCAUUaCCUUUUAAUUUG-------------*[ 6]**[ 3]*ACAUUUU 681747
                                 ******86555555555443................7.....9..******* PP

>> ENA|BK006948|BK006948.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome XV, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (2) ?      0.49   19.8   0.0  cm        1      193 []      737498      737324 - .. 0.96    no 0.39

                                                                                                                  NC
                                 ::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,<<<<<<~~~~~~>>>>>>,<<<<<<<___>>>>> CS
                       U2      1 AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCUUauCAG*[ 8]*CUGauAuggcccccAuuggggg 82    
                                 AU+CC U U+ GCC U  GGC +A AU AAGU UA UA C GUUCU:A::A        U::U:  ::::::A U:::::
  ENA|BK006948|BK006948.2 737498 AUCCCAUAUUUGCCAUC-GGCAUAUAUUAAGUAUAUUAGCAGUUCUAAUUAC*[88]*GUAGUUGGAAGGAUACUAUCCU 737338
                                 **************999.*******************************996...*..6999999999999999999999 PP

                                                                                   NC
                                 >>,,,,,,,,,,,,,,,,,,,,,,,,,,,,~~~~~~~~~~~~::::::: CS
                       U2     83 ccaauuauaUUAaauuaAUUUUUggaacua*[34]**[40]*Acccuuu 193   
                                 : A+                                      A CC++U
  ENA|BK006948|BK006948.2 737337 UUAU--------------------------*[ 2]**[ 1]*AUCCCCU 737324
                                 9987.............................6.....9..******* PP

>> ENA|BK006947|BK006947.3  TPA_inf: Saccharomyces cerevisiae S288C chromosome XIV, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (3) ?       5.7   15.3   0.0  cm        1      193 []      266059      266208 + .. 0.91    no 0.39

                                             v          v                                                     v   NC
                                 ::::::.<<<.-<<<<____>>>>->>>,,,,.,,,,,,~~~~~~,,,,,,,,,,,,,,,,,,,,,,,,,,,,<<<<<<< CS
                       U2      1 AUacCU.UCu.cgGCcUUUUgGCuaaGAUCAA.GUGUAG*[48]*aauuauaUUAaauuaAUUUUUggaacuaGuggggg 119   
                                 AU     UCU + G C  UUG C  AGAU A  GUGUAG        UUAUAU   +UU AU UUU G   +A:: : G:
  ENA|BK006947|BK006947.3 266059 AUGUUGaUCUaUCGUCAAUUGACCCAGAUGAUaGUGUAG*[ 1]*-GUUAUAUAGUUUUGAUAUUUUGGCGAAAAGUUGA 266132
                                 *****9****999****************9988999987...5...3377778888888888888888888888888888 PP

                                                          v      v  v     v v             v v      v  v        NC
                                 <----.<<<<<__>>>>>-..->>>>>>>>,,<<<<<<-<<<<<<___________>>>>>>-->>>>>>::::::: CS
                       U2    120 cauuu.uggGCUUGCccau..ugcccccaCacggguugaccuggcaUUGCAcUaccgccagguucagcccAcccuuu 193   
                                 :A+U+ U :GCUUGC: AU  +::C : ::   G:  :AC: G   U GCA UA+   C :GU+:  :C    +U +
  ENA|BK006947|BK006947.3 266133 GAAUAuUGCGCUUGCGUAUauAUUCCAUUUGAGGUGGCACUAGAGCUCGCAUUAU-UACCAGUAGUGGCAGGAUUGC 266208
                                 888888**************99999******************************.99******************* PP

>> ENA|BK006949|BK006949.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome XVI, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (4) ?       6.6   15.1   0.4  cm        1      193 []      443393      443253 - .. 0.71    no 0.32

                                                                                         v  v   v  v              NC
                                 ::::::<<<-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,~~~~~~,<<<<<<<___>>>>>>>,,,..,,,,, CS
                       U2      1 AUacCUUCucgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCU*[20]*uggcccccAuugggggccaau..uauaU 92    
                                   A  U C:+G  C    G  UA:G UCAAG  UAGUAU UGUUCU       ::: C  A+   G :::++U       
  ENA|BK006949|BK006949.2 443393 CCAUUUACAUGAACCCCAGUUUAUGUUCAAG--UAGUAUAUGUUCU*[ 4]*-UCAACGCAAACUGCUGAUUUcaCC--- 443322
                                 ***************999************6..7**********98...4...222222222222222222222222... PP

                                                         v           v  v         v      v   v    v               NC
                                 ,,,,,,,,,,,,,,,,,,,,<<<<<<<<----<<<<<__>>>>>-->>>>>>>>,,<<<<<<-<<<<<<___________ CS
                       U2     93 UAaauuaAUUUUUggaacuaGugggggcauuuuggGCUUGCccauugcccccaCacggguugaccuggcaUUGCAcUacc 172   
                                                     :U : :::+U U     UU+      ::: : A:A+ ::: G+C: : :A U CAC AC 
  ENA|BK006949|BK006949.2 443321 --------------------AUUGAAUAUUGU-----UUA------UAUAUGAUAUAUACCGUCAAAUUACUUCACGAC- 443274
                                 ....................222222222222.....221......3333345599***********************. PP

                                    v     v   v        NC
                                 >>>>>>-->>>>>>::::::: CS
                       U2    173 gccagguucagcccAcccuuu 193   
                                 : : :G++C :::  +  U+ 
  ENA|BK006949|BK006949.2 443273 AGUGUGAACUGUGAUAAAUCA 443253
                                 ********************* PP

>> ENA|BK006939|BK006939.2  TPA_inf: Saccharomyces cerevisiae S288C chromosome V, complete sequence.
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
  (5) ?       7.1   14.9   0.0  cm        1      193 []      190882      191043 + .. 0.92    no 0.41

                                             v        v                                                           NC
                                 ::::::<<<.-<<<<____>>>>->>>,,,,,,,,,,,,,,,,,,,,~~~~~~~~~~~~,,,,,,,,,,,,,,,,,,,,, CS
                       U2      1 AUacCUUCu.cgGCcUUUUgGCuaaGAUCAAGUGUAGUAUCUGUUCU*[20]**[18]*aauuauaUUAaauuaAUUUUU 105   
                                 AU CC U : +: C: UUU:G : : AUCA                             AA ++U+UUAAAU AA UUUU
  ENA|BK006939|BK006939.2 190882 AUUCCAUGAuUUCCUGUUUAGCU-UCAUCA-----------------*[ 4]**[ 4]*AACAUUUUUAAAUGAAAUUUU 190939
                                 *******9955678899999976.799997....................9.....9..77699**************** PP

                                            vv v                             v vv              v                  NC
                                 ,,,,,,,<<<<<<<<----.<<<<<_._>>>>>-.........->>>>>>>>,,<<<<<<-<<<<<<__________... CS
                       U2    106 ggaacuaGugggggcauuu.uggGCU.UGCccau.........ugcccccaCacggguugaccuggcaUUGCAcUac... 171   
                                  +  + +GU::  : AUU  :G:GCU UGC:C:          + :  ::ACA+ :: :: : ::::  U CAC AC   
  ENA|BK006939|BK006939.2 190940 AAUGUCUGUUUCCUUAUUGaAGAGCUuUGCUCUGgauuuuccaACAUUAAACAUGCCGCCGAGGCCUCCUCCACCACcac 191019
                                 ***********9999998899999964799999999999999999999**************************999888 PP

                                        v                 NC
                                 .._>>>>>>-->>>>>>::::::: CS
                       U2    172 ..cgccagguucagcccAcccuuu 193   
                                   +:::: :UU:: ::  + CU+U
  ENA|BK006939|BK006939.2 191020 caUUGGCAUUUGGUGGUGAACUAU 191043
                                 988********************* PP



Internal CM pipeline statistics summary:
----------------------------------------
Query model(s):                                                  1  (193 consensus positions)
Target sequences:                                               16  (24142652 residues searched)
Target sequences re-searched for truncated hits:                16  (15424 residues re-searched)
Windows   passing  local HMM SSV           filter:           72732  (0.7978); expected (0.35)
Windows   passing  local HMM Viterbi       filter:           24175  (0.3233); expected (0.15)
Windows   passing  local HMM Viterbi  bias filter:            6671  (0.09669); expected (0.15)
Windows   passing  local HMM Forward       filter:            2037  (0.03161); expected (0.003)
Windows   passing  local HMM Forward  bias filter:            1133  (0.01757); expected (0.003)
Windows   passing glocal HMM Forward       filter:             596  (0.01251); expected (0.003)
Windows   passing glocal HMM Forward  bias filter:             438  (0.009175); expected (0.003)
Envelopes passing glocal HMM envelope defn filter:             460  (0.00429); expected (0.003)
Envelopes passing  local CM  CYK           filter:              38  (0.000201); expected (0.0001)
Total CM hits reported:                                          5  (3.063e-05); includes 0 truncated hit(s)

# CPU time: 67.17u 2.23s 00:01:09.40 Elapsed: 00:00:03.24
//
"""
        self.check_raw(filename, "U2", raw)


class Hmmer3TextIndexCases(CheckIndex):
    fmt = "infernal-text"

    def test_infernal_text_1q_0m(self):
        """Test infernal-text indexing, cmsearch, one queries, no hits"""
        filename = os.path.join("Infernal", "cmsearch_114_IRES_Yeast.txt")
        self.check_index(filename, self.fmt)

    def test_infernal_text_1q_mm(self):
        """Test infernal-text indexing, cmsearch, one queries, multiple hits"""
        filename = os.path.join("Infernal", "cmsearch_114_5S_Yeast.txt")
        self.check_index(filename, self.fmt)

    def test_infernal_text_mq(self):
        """Test infernal-text indexing, cmsearch, multiple queries"""
        filename = os.path.join("Infernal", "cmsearch_114_IRES_5S_U2_Yeast.txt")
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
