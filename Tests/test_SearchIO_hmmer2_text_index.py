# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO hmmer2-text indexing."""

import os
import unittest

from search_tests_common import CheckRaw, CheckIndex


class Hmmer2TextRawCases(CheckRaw):

    fmt = 'hmmer2-text'

    def test_hmmer2text_22_single_hmmsearch(self):
        """Test hmmer2-text raw string retrieval, single query, hmmsearch"""
        filename = os.path.join('Hmmer', 'text_22_hmmsearch_001.out')
        raw = """hmmsearch - search a sequence database with a profile HMM
HMMER 2.2g (August 2001)
Copyright (C) 1992-2001 HHMI/Washington University School of Medicine
Freely distributed under the GNU General Public License (GPL)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
HMM file:                   Peptidase_C1.hmm [Peptidase_C1]
Sequence database:          cysprot1b.fa
per-sequence score cutoff:  [none]
per-domain score cutoff:    [none]
per-sequence Eval cutoff:   <= 10        
per-domain Eval cutoff:     [none]
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query HMM:   Peptidase_C1
Accession:   PF00112
Description: Papain family cysteine protease
  [HMM has been calibrated; E-values are empirical estimates]

Scores for complete sequences (score includes all domains):
Sequence   Description                                  Score    E-value  N 
--------   -----------                                  -----    ------- ---
CATL_RAT                                                449.4     2e-135   1
CATL_HUMAN                                              444.5   6.1e-134   1
CATH_RAT                                                381.8   4.8e-115   1
PAPA_CARPA                                              337.7     9e-102   1

Parsed for domains:
Sequence   Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
--------   ------- ----- -----    ----- -----      -----  -------
CATL_RAT     1/1     114   332 ..     1   337 []   449.4   2e-135
CATL_HUMAN   1/1     114   332 ..     1   337 []   444.5 6.1e-134
CATH_RAT     1/1     114   330 ..     1   337 []   381.8 4.8e-115
PAPA_CARPA   1/1     134   343 ..     1   337 []   337.7   9e-102

Alignments of top-scoring domains:
CATL_RAT: domain 1 of 1, from 114 to 332: score 449.4, E = 2e-135
                   *->lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgryciktgt
                      +P+++DWRe kg  VtpVK+QG qCGSCWAFSa g lEg+ ++kt  
    CATL_RAT   114    IPKTVDWRE-KG-CVTPVKNQG-QCGSCWAFSASGCLEGQMFLKT-- 155  

                   kawggklvsLSEQqLvDCdgedygnngesCGyGCnGGGlmdnAfeYikke
                       gkl+sLSEQ+LvDC++ d gn+      GCnG Glmd Af+Yik+ 
    CATL_RAT   156 ----GKLISLSEQNLVDCSH-DQGNQ------GCNG-GLMDFAFQYIKE- 192  

                   qIsnNgGlvtEsdYekgCkPYtdfPCgkdggndtyypCpgkaydpndTgt
                       NgGl++E++Y     PY+    +kd                   g+
    CATL_RAT   193 ----NGGLDSEESY-----PYE----AKD-------------------GS 210  

                   CkynckknskypktyakikgygdvpynvsTydEealqkalaknGPvsVai
                   Cky+  + ++     a+++g++d+p++     E+al+ka+a++GP+sVa+
    CATL_RAT   211 CKYR-AEYAV-----ANDTGFVDIPQQ-----EKALMKAVATVGPISVAM 249  

                   dasedskgDFqlYksGendvgyGvYkhtsageCggtpfteLdHAVliVGY
                   das+ s    q+Y+sG       +Y+++    C+++   +LdH+Vl+VGY
    CATL_RAT   250 DASHPS---LQFYSSG-------IYYEP---NCSSK---DLDHGVLVVGY 283  

                   GteneggtfdetssskksesgiqvssgsngssgSSgssgapiedkgkdYW
                   G e+                                      ++++ +YW
    CATL_RAT   284 GYEG-T------------------------------------DSNKDKYW 296  

                   IVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi<-*
                   +VKNSWG++WG++GY++ia+++n    n+CG+a+ asypi   
    CATL_RAT   297 LVKNSWGKEWGMDGYIKIAKDRN----NHCGLATAASYPI    332  

CATL_HUMAN: domain 1 of 1, from 114 to 332: score 444.5, E = 6.1e-134
                   *->lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgryciktgt
                      +P s+DWRe kg +VtpVK+QG qCGSCWAFSa+galEg+ ++kt  
  CATL_HUMAN   114    APRSVDWRE-KG-YVTPVKNQG-QCGSCWAFSATGALEGQMFRKT-- 155  

                   kawggklvsLSEQqLvDCdgedygnngesCGyGCnGGGlmdnAfeYikke
                       g l+sLSEQ+LvDC+g + gn+      GCnG Glmd+Af+Y+++ 
  CATL_HUMAN   156 ----GRLISLSEQNLVDCSG-PQGNE------GCNG-GLMDYAFQYVQD- 192  

                   qIsnNgGlvtEsdYekgCkPYtdfPCgkdggndtyypCpgkaydpndTgt
                       NgGl++E++Y     PY+    +++                    +
  CATL_HUMAN   193 ----NGGLDSEESY-----PYE----ATE-------------------ES 210  

                   CkynckknskypktyakikgygdvpynvsTydEealqkalaknGPvsVai
                   Ckyn +k s+     a+++g++d+p +     E+al+ka+a++GP+sVai
  CATL_HUMAN   211 CKYN-PKYSV-----ANDTGFVDIPKQ-----EKALMKAVATVGPISVAI 249  

                   dasedskgDFqlYksGendvgyGvYkhtsageCggtpfteLdHAVliVGY
                   da++ s   F +Yk G       +Y ++   +C+++   + dH+Vl+VGY
  CATL_HUMAN   250 DAGHES---FLFYKEG-------IYFEP---DCSSE---DMDHGVLVVGY 283  

                   GteneggtfdetssskksesgiqvssgsngssgSSgssgapiedkgkdYW
                   G e+                                      e+++ +YW
  CATL_HUMAN   284 GFES-T------------------------------------ESDNNKYW 296  

                   IVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi<-*
                   +VKNSWG++WG+ GY+++a+++     n+CGIas asyp+   
  CATL_HUMAN   297 LVKNSWGEEWGMGGYVKMAKDRR----NHCGIASAASYPT    332  

CATH_RAT: domain 1 of 1, from 114 to 330: score 381.8, E = 4.8e-115
                   *->lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgryciktgt
                       P s+DWR+ kg  V+pVK+QG  CGSCW FS++galE++ +i++  
    CATH_RAT   114    YPSSMDWRK-KGNVVSPVKNQG-ACGSCWTFSTTGALESAVAIAS-- 156  

                   kawggklvsLSEQqLvDCdgedygnngesCGyGCnGGGlmdnAfeYikke
                       gk   L EQqLvDC   +++n+      GC+G Gl+++AfeYi++ 
    CATH_RAT   157 ----GKMMTLAEQQLVDCAQ-NFNNH------GCQG-GLPSQAFEYILY- 193  

                   qIsnNgGlvtEsdYekgCkPYtdfPCgkdggndtyypCpgkaydpndTgt
                       N+G++ E++Y     PY     gk+                   g+
    CATH_RAT   194 ----NKGIMGEDSY-----PYI----GKN-------------------GQ 211  

                   CkynckknskypktyakikgygdvpynvsTydEealqkalaknGPvsVai
                   Ck+n +++++     a++k+ ++++ n    dE+a+ +a+a + Pvs a+
    CATH_RAT   212 CKFN-PEKAV-----AFVKNVVNITLN----DEAAMVEAVALYNPVSFAF 251  

                   dasedskgDFqlYksGendvgyGvYkhtsageCggtpfteLdHAVliVGY
                   +++e    DF++YksG       vY++    +C +tp + ++HAVl+VGY
    CATH_RAT   252 EVTE----DFMMYKSG-------VYSSN---SCHKTP-DKVNHAVLAVGY 286  

                   GteneggtfdetssskksesgiqvssgsngssgSSgssgapiedkgkdYW
                   G +n g                                          YW
    CATH_RAT   287 GEQN-GLL----------------------------------------YW 295  

                   IVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi<-*
                   IVKNSWG++WG nGYf i+Rgkn     +CG+a +asypi   
    CATH_RAT   296 IVKNSWGSNWGNNGYFLIERGKN-----MCGLAACASYPI    330  

PAPA_CARPA: domain 1 of 1, from 134 to 343: score 337.7, E = 9e-102
                   *->lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgryciktgt
                      +Pe +DWR+ kg aVtpVK+QG +CGSCWAFSav ++Eg+++i+t  
  PAPA_CARPA   134    IPEYVDWRQ-KG-AVTPVKNQG-SCGSCWAFSAVVTIEGIIKIRT-- 175  

                   kawggklvsLSEQqLvDCdgedygnngesCGyGCnGGGlmdnAfeYikke
                       g+l  +SEQ+L+DCd+    ++      GCnG G+++ A++ + + 
  PAPA_CARPA   176 ----GNLNEYSEQELLDCDR---RSY------GCNG-GYPWSALQLVAQ- 210  

                   qIsnNgGlvtEsdYekgCkPYtdfPCgkdggndtyypCpgkaydpndTgt
                         G++    Y     PY+    g++                     
  PAPA_CARPA   211 -----YGIHYRNTY-----PYE----GVQ-------------------RY 227  

                   CkynckknskypktyakikgygdvpynvsTydEealqkalaknGPvsVai
                   C+++ +k+ +    +ak +g ++v+++    +E al + +a+ +PvsV  
  PAPA_CARPA   228 CRSR-EKGPY----AAKTDGVRQVQPY----NEGALLYSIAN-QPVSVVL 267  

                   dasedskgDFqlYksGendvgyGvYkhtsageCggtpfteLdHAVliVGY
                   +a +    DFqlY++G       ++++    +Cg+     +dHAV++VGY
  PAPA_CARPA   268 EAAGK---DFQLYRGG-------IFVG----PCGN----KVDHAVAAVGY 299  

                   GteneggtfdetssskksesgiqvssgsngssgSSgssgapiedkgkdYW
                   G                                              +Y 
  PAPA_CARPA   300 G---------------------------------------------PNYI 304  

                   IVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi<-*
                   ++KNSWGt WGEnGY+ri+Rg+++s ++ CG+ ++  yp+   
  PAPA_CARPA   305 LIKNSWGTGWGENGYIRIKRGTGNS-YGVCGLYTSSFYPV    343  


Histogram of all scores:
score    obs    exp  (one = represents 1 sequences)
-----    ---    ---
> 337      4      -|====                                                       


% Statistical details of theoretical EVD fit:
              mu =  -195.8384
          lambda =     0.1423
chi-sq statistic =     0.0000
  P(chi-square)  =          0

Total sequences searched: 4

Whole sequence top hits:
tophits_s report:
     Total hits:           4
     Satisfying E cutoff:  4
     Total memory:         16K

Domain top hits:
tophits_s report:
     Total hits:           4
     Satisfying E cutoff:  4
     Total memory:         20K
"""  # noqa : W291
        self.check_raw(filename, 'Peptidase_C1', raw)

    def test_hmmer2text_22_single_hmmpfam(self):
        """Test hmmer2-text raw string retrieval, single query, hmmpfam"""
        filename = os.path.join('Hmmer', 'text_22_hmmpfam_001.out')
        raw = """hmmpfam - search one or more sequences against HMM database
HMMER 2.2g (August 2001)
Copyright (C) 1992-2001 HHMI/Washington University School of Medicine
Freely distributed under the GNU General Public License (GPL)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
HMM file:                 Pfam
Sequence file:            L77119.faa
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query sequence: gi|1522636|gb|AAC37060.1|
Accession:      [none]
Description:    M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]

Scores for sequence family classification (score includes all domains):
Model       Description                                 Score    E-value  N 
--------    -----------                                 -----    ------- ---
Methylase_M Type I restriction modification system, M  -105.2     0.0022   1

Parsed for domains:
Model       Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
--------    ------- ----- -----    ----- -----      -----  -------
Methylase_M   1/1     280   481 ..     1   279 []  -105.2   0.0022

Alignments of top-scoring domains:
Methylase_M: domain 1 of 1, from 280 to 481: score -105.2, E = 0.0022
                   *->lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerri
                       ++EL+++  av+   R              L+F K++ dk      
  gi|1522636   280    NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------ 306  

                   eieerktdtesepsldyakledqyeqlededlekedfyqkkGvFilPsql
                   +i+         p +   + +++y   ++   ++ ++y ++      + l
  gi|1522636   307 GIV---------PRDLLRRTYEDY---KKSNVLI-NYYDAY-L----KPL 338  

                   FwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdl
                   F++++   e ++  ++++ + +    ++      + +       Glf ++
  gi|1522636   339 FYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSN 374  

                   dfnsnkLgskaqarnetLtelidlfselelgtPmHNG.dfeelgikDlfG
                   ++  ++ +s+   +ne ++e+i+ +++ +++     G++ +el   D++G
  gi|1522636   375 NV--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILG 421  

                   DaYEYLLgkFAeneGKsGGeFYTPqeVSkLiaeiLtigqpsegdfsIYDP
                    +YE L+   Ae   K+ G +YTP e++  ia+ + i+  ++        
  gi|1522636   422 YVYEKLINILAEKGQKGLGAYYTPDEITSYIAKNT-IEPIVVE------- 463  

                   AcGSGSLllqaskflgehdgkrnaisyYGQEsn<-*
                             +++ ++    k+n+i +    s+   
  gi|1522636   464 ---------RFKEIIK--NWKINDINF----ST    481  

//
"""  # noqa : W291
        self.check_raw(filename, 'gi|1522636|gb|AAC37060.1|', raw)

    def test_hmmer2text_22_multiple_first_hmmpfam(self):
        """Test hmmer2-text raw string retrieval, multiple queries, hmmpfam"""
        filename = os.path.join('Hmmer', 'text_24_hmmpfam_001.out')
        raw = """hmmpfam - search one or more sequences against HMM database
HMMER 2.4i (December 2006)
Copyright (C) 1992-2006 HHMI Janelia Farm
Freely distributed under the GNU General Public License (GPL)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
HMM file:                 /home/bow/db/hmmer/Pfam_fs
Sequence file:            mult.fasta
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query sequence: random_s00
Accession:      [none]
Description:    [none]

Scores for sequence family classification (score includes all domains):
Model           Description                             Score    E-value  N 
--------        -----------                             -----    ------- ---
	[no hits above thresholds]

Parsed for domains:
Model           Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
--------        ------- ----- -----    ----- -----      -----  -------
	[no hits above thresholds]

Alignments of top-scoring domains:
	[no hits above thresholds]
//
"""  # noqa : W291
        self.check_raw(filename, 'random_s00', raw)  # noqa : E101

    def test_hmmer2text_22_multiple_middle_hmmpfam(self):
        """Test hmmer2-text raw string retrieval, multiple queries, hmmpfam"""
        filename = os.path.join('Hmmer', 'text_24_hmmpfam_001.out')
        raw = """hmmpfam - search one or more sequences against HMM database
HMMER 2.4i (December 2006)
Copyright (C) 1992-2006 HHMI Janelia Farm
Freely distributed under the GNU General Public License (GPL)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
HMM file:                 /home/bow/db/hmmer/Pfam_fs
Sequence file:            mult.fasta
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query sequence: gi|4885477|ref|NP_005359.1|
Accession:      [none]
Description:    myoglobin [Homo sapiens]

Scores for sequence family classification (score includes all domains):
Model         Description                               Score    E-value  N 
--------      -----------                               -----    ------- ---
Globin        Globin                                    129.8    5.8e-37   1
tRNA-synt_1b  tRNA synthetases class I (W and Y)          1.8        2.2   1
Rotavirus_VP3 Rotavirus VP3 protein                      -1.2        7.9   1
DTHCT         DTHCT (NUC029) region                       1.5        9.8   1

Parsed for domains:
Model         Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
--------      ------- ----- -----    ----- -----      -----  -------
Globin          1/1       7   143 ..     1   148 []   129.8  5.8e-37
DTHCT           1/1     101   106 ..   101   106 .]     1.5      9.8
tRNA-synt_1b    1/1     120   137 ..   320   337 .]     1.8      2.2
Rotavirus_VP3   1/1     134   147 ..     1    15 [.    -1.2      7.9

Alignments of top-scoring domains:
Globin: domain 1 of 1, from 7 to 143: score 129.8, E = 5.8e-37
                   *->dkalvkasWgkvkgtdnreelGaealarlFkayPdtktyFpkfgdls
                      +++lv+  Wgkv++  +++ +G+e+l rlFk +P+t ++F kf+ l+
  gi|4885477     7    EWQLVLNVWGKVEA--DIPGHGQEVLIRLFKGHPETLEKFDKFKHLK 51   

                   sadaikgspkfkaHgkkVlaalgeavkhLgnddddgnlkaalkkLaarHa
                   s d++k s+++k+Hg++Vl alg ++k       +g + a +k La +Ha
  gi|4885477    52 SEDEMKASEDLKKHGATVLTALGGILKK------KGHHEAEIKPLAQSHA 95   

                   erghvdpanFkllgeallIvvLaahlggeveftpevkaAWdkaldvvada
                   ++++++ ++ + ++e+++ +vL+++ +g  +f +++++A++kal  +  +
  gi|4885477    96 TKHKIPVKYLEFISECII-QVLQSKHPG--DFGADAQGAMNKALELFRKD 142  

                   l<-*
                   +   
  gi|4885477   143 M    143  

DTHCT: domain 1 of 1, from 101 to 106: score 1.5, E = 9.8
                   *->pvKYLe<-*
                      pvKYLe   
  gi|4885477   101    PVKYLE    106  

tRNA-synt_1b: domain 1 of 1, from 120 to 137: score 1.8, E = 2.2
                   *->hggelKkaaaeavnalls<-*
                      h+g++  +a+ a+n++l+   
  gi|4885477   120    HPGDFGADAQGAMNKALE    137  

Rotavirus_VP3: domain 1 of 1, from 134 to 147: score -1.2, E = 7.9
                   *->MkVLaLFrrgvalnY<-*
                       k L LFr+++a+nY   
  gi|4885477   134    -KALELFRKDMASNY    147  

//
"""  # noqa : W291
        self.check_raw(filename, 'gi|4885477|ref|NP_005359.1|', raw)

    def test_hmmer2text_22_multiple_last_hmmpfam(self):
        """Test hmmer2-text raw string retrieval, multiple queries, hmmpfam"""
        filename = os.path.join('Hmmer', 'text_24_hmmpfam_001.out')
        raw = """hmmpfam - search one or more sequences against HMM database
HMMER 2.4i (December 2006)
Copyright (C) 1992-2006 HHMI Janelia Farm
Freely distributed under the GNU General Public License (GPL)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
HMM file:                 /home/bow/db/hmmer/Pfam_fs
Sequence file:            mult.fasta
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query sequence: gi|125490392|ref|NP_038661.2|
Accession:      [none]
Description:    POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]

Scores for sequence family classification (score includes all domains):
Model           Description                             Score    E-value  N 
--------        -----------                             -----    ------- ---
Pou             Pou domain - N-terminal to homeobox d   171.2    2.4e-48   1
Homeobox        Homeobox domain                          86.7    6.5e-23   1
HTH_3           Helix-turn-helix                          7.5       0.33   1
WSC             WSC domain                                2.1        1.5   1
ComC            COMC family                               3.4        2.3   1
CBM_1           Fungal cellulose binding domain           3.3        3.6   1
Peptidase_M29   Thermophilic metalloprotease (M29)       -2.1        4.4   1
DUF1690         Protein of Unknown function (DUF1690)     1.1        4.5   1
DUF137          Protein of unknown function DUF137        1.0        4.7   1
DASH_Duo1       DASH complex subunit Duo1                 2.3        5.7   1
TFIIB           Transcription factor TFIIB repeat         2.6        6.1   1
DUF1392         Protein of unknown function (DUF1392)     0.8        6.7   1

Parsed for domains:
Model           Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
--------        ------- ----- -----    ----- -----      -----  -------
WSC               1/1      59    67 ..    80    87 .]     2.1      1.5
CBM_1             1/1      61    73 ..     1    13 [.     3.3      3.6
Pou               1/1     131   205 ..     1    78 []   171.2  2.4e-48
ComC              1/1     132   145 ..     1    20 [.     3.4      2.3
Peptidase_M29     1/1     132   145 ..     1    14 [.    -2.1      4.4
DASH_Duo1         1/1     133   141 ..     1     9 [.     2.3      5.7
DUF1690           1/1     137   147 ..   151   161 .]     1.1      4.5
HTH_3             1/1     146   166 ..     1    25 [.     7.5     0.33
DUF137            1/1     197   213 ..   164   180 ..     1.0      4.7
Homeobox          1/1     224   280 ..     1    57 []    86.7  6.5e-23
DUF1392           1/1     248   269 ..   127   164 .]     0.8      6.7
TFIIB             1/1     253   268 ..     1    19 [.     2.6      6.1

Alignments of top-scoring domains:
WSC: domain 1 of 1, from 59 to 67: score 2.1, E = 1.5
                   *->p.pseiCGG<-*
                      p+++e CGG   
  gi|1254903    59    PpAYEFCGG    67   

CBM_1: domain 1 of 1, from 61 to 73: score 3.3, E = 3.6
                   *->vygQCGGigysGp<-*
                      +y+ CGG+ y Gp   
  gi|1254903    61    AYEFCGGMAYCGP    73   

Pou: domain 1 of 1, from 131 to 205: score 171.2, E = 2.4e-48
                   *->deatdleeLEkFAkeFKqRRIkLGyTQadVGlALgalygPGvnafSQ
                      d+++ ++eLE+FAk +Kq+RI+LGyTQadVGl+Lg l+g   ++fSQ
  gi|1254903   131    DMKALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFG---KVFSQ 174  

                   tTICRFEaLqLSfKNmcKLKPlLekWLeeAE<-*
                   tTICRFEaLqLS KNmcKL+PlLekW+eeA+   
  gi|1254903   175 TTICRFEALQLSLKNMCKLRPLLEKWVEEAD    205  

ComC: domain 1 of 1, from 132 to 145: score 3.4, E = 2.3
                   *->MKntvkkkllkkeLeqFkeL<-*
                      MK ++k      eLeqF  L   
  gi|1254903   132    MKALQK------ELEQFAKL    145  

Peptidase_M29: domain 1 of 1, from 132 to 145: score -2.1, E = 4.4
                   *->MdaFkkeLekyAeL<-*
                      M a +keLe +A L   
  gi|1254903   132    MKALQKELEQFAKL    145  

DASH_Duo1: domain 1 of 1, from 133 to 141: score 2.3, E = 5.7
                   *->aaLqkELeq<-*
                      +aLqkELeq   
  gi|1254903   133    KALQKELEQ    141  

DUF1690: domain 1 of 1, from 137 to 147: score 1.1, E = 4.5
                   *->eEvEqFKklvr<-*
                      +E EqF+kl++   
  gi|1254903   137    KELEQFAKLLK    147  

HTH_3: domain 1 of 1, from 146 to 166: score 7.5, E = 0.33
                   *->lkelRkkkelglsqeeLAeklGskv<-*
                      lk++R    lg++q+++   lG  v   
  gi|1254903   146    LKQKRI--TLGYTQADVGLTLG--V    166  

DUF137: domain 1 of 1, from 197 to 213: score 1.0, E = 4.7
                   *->LeeIVEnyDNkKnLkEv<-*
                      Le+ VE+ DN++nL+E+   
  gi|1254903   197    LEKWVEEADNNENLQEI    213  

Homeobox: domain 1 of 1, from 224 to 280: score 86.7, E = 6.5e-23
                   *->rrkRTtftpeQleeLEkeFqknrYPsreeReeLAkkLgLterqVkvW
                      +rkRT++ +  +  LE +F k+++Ps +++ ++A++LgL++++V+vW
  gi|1254903   224    KRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEKDVVRVW 270  

                   FQNRRaKwKk<-*
                   F+NRR+K K+   
  gi|1254903   271 FCNRRQKGKR    280  

DUF1392: domain 1 of 1, from 248 to 269: score 0.8, E = 6.7
                   *->PtLsqtttEGlCIFPrssqGnkmpnRfsLvrerDLVrV<-*
                      P+L q+t                +n++ L  e+D VrV   
  gi|1254903   248    PSLQQITH--------------IANQLGL--EKDVVRV    269  

TFIIB: domain 1 of 1, from 253 to 268: score 2.6, E = 6.1
                   *->ikrfadaLeLpeKkikVad<-*
                      i+++a +L+L +    V++   
  gi|1254903   253    ITHIANQLGLEK---DVVR    268  

//
"""  # noqa : W291
        self.check_raw(filename, 'gi|125490392|ref|NP_038661.2|', raw)


class Hmmer2TextIndexCases(CheckIndex):

    fmt = 'hmmer2-text'

    def test_hmmertext_text_21_hmmpfam_001(self):
        """Test hmmer2-text indexing, HMMER 2.1"""
        filename = os.path.join('Hmmer', 'text_21_hmmpfam_001.out')
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_22_hmmpfam_001(self):
        """Test hmmer2-text indexing, HMMER 2.2"""
        filename = os.path.join('Hmmer', 'text_22_hmmpfam_001.out')
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_23_hmmpfam_001(self):
        """Test hmmer2-text indexing, HMMER 2.3"""
        filename = os.path.join('Hmmer', 'text_23_hmmpfam_001.out')
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_24_hmmpfam_001(self):
        """Test hmmer2-text indexing, HMMER 2.4"""
        filename = os.path.join('Hmmer', 'text_24_hmmpfam_001.out')
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_22_hmmsearch_001(self):
        """Test hmmer2-text indexing, HMMER 2.2"""
        filename = os.path.join('Hmmer', 'text_22_hmmsearch_001.out')
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
