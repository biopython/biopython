# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO hmmer3-text indexing."""

import unittest

from search_tests_common import CheckRaw, CheckIndex


class Hmmer3TextRawCases(CheckRaw):

    fmt = 'hmmer3-text'

    def test_hmmer3text_30_multiple_first(self):
        """Test hmmer3-text raw string retrieval, HMMER 3.0, multiple queries, first (text_30_hmmscan_001.out)"""
        filename = 'Hmmer/text_30_hmmscan_001.out'
        raw = """# hmmscan :: search sequence(s) against a profile database
# HMMER 3.0 (March 2010); http://hmmer.org/
# Copyright (C) 2010 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             mult.fasta
# target HMM database:             /home/bow/db/hmmer/Pfam-A.hmm
# output directed to file:         hmmer_cases/text_hmmscan_mult.out
# per-seq hits tabular output:     hmmer_cases/tab_hmmscan_mult.out
# per-dom hits tabular output:     hmmer_cases/domtab_hmmscan_mult.out
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       random_s00  [L=32]
Scores for complete sequence (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Model    Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------

   [No hits detected that satisfy reporting thresholds]


Domain annotation for each model (and alignments):

   [No targets detected that satisfy reporting thresholds]


Internal pipeline statistics summary:
-------------------------------------
Query sequence(s):                         1  (32 residues)
Target model(s):                       13672  (2396357 nodes)
Passed MSV filter:                       338  (0.0247221); expected 273.4 (0.02)
Passed bias filter:                       87  (0.00636337); expected 273.4 (0.02)
Passed Vit filter:                        23  (0.00168227); expected 13.7 (0.001)
Passed Fwd filter:                        14  (0.00102399); expected 0.1 (1e-05)
Initial search space (Z):              13672  [actual number of targets]
Domain search space  (domZ):               0  [number of targets reported over threshold]
# CPU time: 0.20u 0.12s 00:00:00.32 Elapsed: 00:00:00.19
# Mc/sec: 403.60
//
"""
        self.check_raw(filename, "random_s00", raw)

    def test_hmmer3text_30_multiple_middle(self):
        """Test hmmer3-text raw string retrieval, HMMER 3.0, multiple queries, middle (text_30_hmmscan_001.out)"""
        filename = 'Hmmer/text_30_hmmscan_001.out'
        raw = """# hmmscan :: search sequence(s) against a profile database
# HMMER 3.0 (March 2010); http://hmmer.org/
# Copyright (C) 2010 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             mult.fasta
# target HMM database:             /home/bow/db/hmmer/Pfam-A.hmm
# output directed to file:         hmmer_cases/text_hmmscan_mult.out
# per-seq hits tabular output:     hmmer_cases/tab_hmmscan_mult.out
# per-dom hits tabular output:     hmmer_cases/domtab_hmmscan_mult.out
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       gi|4885477|ref|NP_005359.1|  [L=154]
Description: myoglobin [Homo sapiens]
Scores for complete sequence (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Model    Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
      6e-21   74.6   0.3    9.2e-21   74.0   0.2    1.3  1  Globin   Globin


Domain annotation for each model (and alignments):
>> Globin  Globin
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   74.0   0.2   6.7e-25   9.2e-21       1     107 [.       7     112 ..       7     113 .. 0.97

  Alignments for each domain:
  == domain 1    score: 74.0 bits;  conditional E-value: 6.7e-25
                                  HHHHHHHHHHHHCHHHHHHHHHHHHHHHHHSGGGGGGGCCCTTTT.HHHHHTSCHHHHHHHHHHHHHHHHHHCTTSHHHHHH CS
                       Globin   1 qkalvkaswekvkanaeeigaeilkrlfkaypdtkklFkkfgdls.aedlksspkfkahakkvlaaldeavknldnddnlka 81 
                                  +++lv   w+kv+a+++ +g+e+l rlfk +p+t ++F kf+ l+  +++k s+++k+h+++vl al+ ++k+   ++ ++a
  gi|4885477|ref|NP_005359.1|   7 EWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKsEDEMKASEDLKKHGATVLTALGGILKK---KGHHEA 85 
                                  5789*********************************************************************...6899** PP

                                  HHHHHHHHHHTT-.--HHHHCCHHHHH CS
                       Globin  82 alkklgarHakrg.vdpanfklfgeal 107
                                  ++k l+++Ha+++ ++ ++ + ++e++
  gi|4885477|ref|NP_005359.1|  86 EIKPLAQSHATKHkIPVKYLEFISECI 112
                                  *********************999998 PP



Internal pipeline statistics summary:
-------------------------------------
Query sequence(s):                         1  (154 residues)
Target model(s):                       13672  (2396357 nodes)
Passed MSV filter:                       458  (0.0334991); expected 273.4 (0.02)
Passed bias filter:                      404  (0.0295494); expected 273.4 (0.02)
Passed Vit filter:                        31  (0.00226741); expected 13.7 (0.001)
Passed Fwd filter:                         1  (7.31422e-05); expected 0.1 (1e-05)
Initial search space (Z):              13672  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.33u 0.16s 00:00:00.49 Elapsed: 00:00:00.21
# Mc/sec: 1757.33
//
"""  # noqa : W291
        self.check_raw(filename, "gi|4885477|ref|NP_005359.1|", raw)

    def test_hmmer3text_30_multiple_last(self):
        """Test hmmer3-text raw string retrieval, HMMER 3.0, multiple queries, last (text_30_hmmscan_001.out)"""
        filename = 'Hmmer/text_30_hmmscan_001.out'
        raw = """# hmmscan :: search sequence(s) against a profile database
# HMMER 3.0 (March 2010); http://hmmer.org/
# Copyright (C) 2010 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             mult.fasta
# target HMM database:             /home/bow/db/hmmer/Pfam-A.hmm
# output directed to file:         hmmer_cases/text_hmmscan_mult.out
# per-seq hits tabular output:     hmmer_cases/tab_hmmscan_mult.out
# per-dom hits tabular output:     hmmer_cases/domtab_hmmscan_mult.out
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       gi|125490392|ref|NP_038661.2|  [L=352]
Description: POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]
Scores for complete sequence (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Model       Description
    ------- ------ -----    ------- ------ -----   ---- --  --------    -----------
      7e-37  124.8   0.5    1.4e-36  123.9   0.3    1.5  1  Pou         Pou domain - N-terminal to homeobox domain
    2.1e-18   65.5   1.1    4.1e-18   64.6   0.7    1.5  1  Homeobox    Homeobox domain
  ------ inclusion threshold ------
      0.012   15.6   0.0       0.16   12.0   0.0    2.2  2  HTH_31      Helix-turn-helix domain
      0.039   13.5   0.0      0.095   12.3   0.0    1.6  1  Homeobox_KN Homeobox KN domain
       0.14   10.5   0.1       0.26    9.6   0.1    1.4  1  DUF521      Protein of unknown function (DUF521)


Domain annotation for each model (and alignments):
>> Pou  Pou domain - N-terminal to homeobox domain
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  123.9   0.3     5e-40   1.4e-36       3      75 .]     133     205 ..     131     205 .. 0.97

  Alignments for each domain:
  == domain 1    score: 123.9 bits;  conditional E-value: 5e-40
                            Pou   3 eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae 75 
                                    ++ ++ele+fak +kq+ri+Lg+tqadvg +lg+l+Gk+fsqttIcrFEalqLslknmckL+pllekW+eea+
  gi|125490392|ref|NP_038661.2| 133 KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD 205
                                    67899******************************************************************96 PP

>> Homeobox  Homeobox domain
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   64.6   0.7   1.5e-21   4.1e-18       1      57 []     224     280 ..     224     280 .. 0.98

  Alignments for each domain:
  == domain 1    score: 64.6 bits;  conditional E-value: 1.5e-21
                                    SS--SS--HHHHHHHHHHCCTSSS--HHHHHHHHHH----HHHHHHHHHHHHHHHHH CS
                       Homeobox   1 rrkRttftkeqleeLeelFeknrypsaeereeLAkklgLterqVkvWFqNrRakekk 57 
                                    +rkRt++++     Le +F k+++ps ++++++A++lgL++++V+vWF+NrR+k k+
  gi|125490392|ref|NP_038661.2| 224 KRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQKGKR 280
                                    79****************************************************997 PP

>> HTH_31  Helix-turn-helix domain
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   12.0   0.0   5.7e-05      0.16       1      35 [.     141     181 ..     141     184 .. 0.96
   2 ?    0.8   0.0      0.19   5.2e+02      39      62 ..     245     268 ..     243     270 .. 0.86

  Alignments for each domain:
  == domain 1    score: 12.0 bits;  conditional E-value: 5.7e-05
                         HTH_31   1 aLGarLralReraGLtqeevAerlg......vSastlsrlE 35 
                                    +++ +L++ R + G tq++v+  lg      +S++t++r E
  gi|125490392|ref|NP_038661.2| 141 QFAKLLKQKRITLGYTQADVGLTLGvlfgkvFSQTTICRFE 181
                                    6999***********************************99 PP

  == domain 2    score: 0.8 bits;  conditional E-value: 0.19
                         HTH_31  39 rgrpsaavlaalaralgldpaera 62 
                                    ++ ps+++++ +a+ lgl+ + ++
  gi|125490392|ref|NP_038661.2| 245 CPKPSLQQITHIANQLGLEKDVVR 268
                                    678**************9988765 PP

>> Homeobox_KN  Homeobox KN domain
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   12.3   0.0   3.5e-05     0.095       7      39 ..     244     276 ..     241     277 .. 0.91

  Alignments for each domain:
  == domain 1    score: 12.3 bits;  conditional E-value: 3.5e-05
                    Homeobox_KN   7 hnPYPskevkeelakqTglsrkqidnWFiNaRr 39 
                                    + P Ps +++  +a+q gl  + +  WF N R 
  gi|125490392|ref|NP_038661.2| 244 KCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQ 276
                                    56779*************************996 PP

>> DUF521  Protein of unknown function (DUF521)
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?    9.6   0.1   9.4e-05      0.26     273     334 ..     221     280 ..     197     294 .. 0.77

  Alignments for each domain:
  == domain 1    score: 9.6 bits;  conditional E-value: 9.4e-05
                         DUF521 273 adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk 334
                                    + +++ + +++++   +++ ++l cP  sl++++++a++l  +k    v+++ +  r+  ++
  gi|125490392|ref|NP_038661.2| 221 QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR 280
                                    345666667778888899************************99..9999999988876554 PP



Internal pipeline statistics summary:
-------------------------------------
Query sequence(s):                         1  (352 residues)
Target model(s):                       13672  (2396357 nodes)
Passed MSV filter:                       603  (0.0441047); expected 273.4 (0.02)
Passed bias filter:                      465  (0.0340111); expected 273.4 (0.02)
Passed Vit filter:                        44  (0.00321826); expected 13.7 (0.001)
Passed Fwd filter:                         5  (0.000365711); expected 0.1 (1e-05)
Initial search space (Z):              13672  [actual number of targets]
Domain search space  (domZ):               5  [number of targets reported over threshold]
# CPU time: 0.51u 0.15s 00:00:00.66 Elapsed: 00:00:00.23
# Mc/sec: 3667.47
//
"""  # noqa : W291
        self.check_raw(filename, "gi|125490392|ref|NP_038661.2|", raw)

    def test_hmmer3text_30_single(self):
        """Test hmmer3-text raw string retrieval, HMMER 3.0, single query (text_30_hmmscan_003.out)"""
        filename = 'Hmmer/text_30_hmmscan_003.out'
        raw = """# hmmscan :: search sequence(s) against a profile database
# HMMER 3.0 (March 2010); http://hmmer.org/
# Copyright (C) 2010 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             s01.fasta
# target HMM database:             /home/bow/db/hmmer/Pfam-A.hmm
# output directed to file:         hmmer_cases/text_hmmscan_s01.out
# per-seq hits tabular output:     hmmer_cases/tab_hmmscan_s01.out
# per-dom hits tabular output:     hmmer_cases/domtab_hmmscan_s01.out
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       gi|4885477|ref|NP_005359.1|  [L=154]
Description: myoglobin [Homo sapiens]
Scores for complete sequence (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Model    Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
      6e-21   74.6   0.3    9.2e-21   74.0   0.2    1.3  1  Globin   Globin


Domain annotation for each model (and alignments):
>> Globin  Globin
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   74.0   0.2   6.7e-25   9.2e-21       1     107 [.       7     112 ..       7     113 .. 0.97

  Alignments for each domain:
  == domain 1    score: 74.0 bits;  conditional E-value: 6.7e-25
                                  HHHHHHHHHHHHCHHHHHHHHHHHHHHHHHSGGGGGGGCCCTTTT.HHHHHTSCHHHHHHHHHHHHHHHHHHCTTSHHHHHH CS
                       Globin   1 qkalvkaswekvkanaeeigaeilkrlfkaypdtkklFkkfgdls.aedlksspkfkahakkvlaaldeavknldnddnlka 81 
                                  +++lv   w+kv+a+++ +g+e+l rlfk +p+t ++F kf+ l+  +++k s+++k+h+++vl al+ ++k+   ++ ++a
  gi|4885477|ref|NP_005359.1|   7 EWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKsEDEMKASEDLKKHGATVLTALGGILKK---KGHHEA 85 
                                  5789*********************************************************************...6899** PP

                                  HHHHHHHHHHTT-.--HHHHCCHHHHH CS
                       Globin  82 alkklgarHakrg.vdpanfklfgeal 107
                                  ++k l+++Ha+++ ++ ++ + ++e++
  gi|4885477|ref|NP_005359.1|  86 EIKPLAQSHATKHkIPVKYLEFISECI 112
                                  *********************999998 PP



Internal pipeline statistics summary:
-------------------------------------
Query sequence(s):                         1  (154 residues)
Target model(s):                       13672  (2396357 nodes)
Passed MSV filter:                       458  (0.0334991); expected 273.4 (0.02)
Passed bias filter:                      404  (0.0295494); expected 273.4 (0.02)
Passed Vit filter:                        31  (0.00226741); expected 13.7 (0.001)
Passed Fwd filter:                         1  (7.31422e-05); expected 0.1 (1e-05)
Initial search space (Z):              13672  [actual number of targets]
Domain search space  (domZ):               1  [number of targets reported over threshold]
# CPU time: 0.28u 0.17s 00:00:00.45 Elapsed: 00:00:00.21
# Mc/sec: 1757.33
//
"""  # noqa : W291
        self.check_raw(filename, "gi|4885477|ref|NP_005359.1|", raw)


class Hmmer3TextIndexCases(CheckIndex):

    fmt = 'hmmer3-text'

    def test_hmmertext_text_30_hmmscan_001(self):
        """Test hmmer3-text indexing, HMMER 3.0, multiple queries"""
        filename = 'Hmmer/text_30_hmmscan_001.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmscan_002(self):
        """Test hmmer3-text indexing, HMMER 3.0, single query, no hits"""
        filename = 'Hmmer/text_30_hmmscan_002.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmscan_006(self):
        """Test hmmer3-text indexing, HMMER 3.0, single query, multiple hits"""
        filename = 'Hmmer/text_30_hmmscan_006.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmscan_007(self):
        """Test hmmer3-text indexing, HMMER 3.0, single query, no alignments"""
        filename = 'Hmmer/text_30_hmmscan_007.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmscan_008(self):
        """Test hmmer3-text indexing, HMMER 3.0, single query, no alignment width"""
        filename = 'Hmmer/text_30_hmmscan_008.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmsearch_005(self):
        """Test hmmer3-text indexing, HMMER 3.0, multiple queries"""
        filename = 'Hmmer/text_30_hmmsearch_005.out'
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
