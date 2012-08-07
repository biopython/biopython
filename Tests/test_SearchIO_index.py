# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO indexing.

Since the SearchIO indexing core shares a lot of similarity with SeqIO
indexing, here we are only testing the format-specific parsers and not
the SearchIO indexing core itself.

"""
# For using with statement in Python 2.5 or Jython
from __future__ import with_statement

import os
import unittest

from Bio import SearchIO
from Bio._py3k import _as_bytes

from search_tests_common import compare_search_obj


class BlastXmlRawCases(unittest.TestCase):

    fmt = 'blast-xml'

    def test_blastxml_2226_multiple_first(self):
        """Test blast-xml raw string retrieval, BLAST 2.2.26+, multiple queries, first (xml_2226_blastp_001.xml)"""
        filename = 'Blast/xml_2226_blastp_001.xml'
        idx = SearchIO.index(filename, self.fmt)
        raw = """<Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>random_s00</Iteration_query-def>
      <Iteration_query-len>32</Iteration_query-len>
      <Iteration_hits></Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>20</Statistics_db-num>
          <Statistics_db-len>6406</Statistics_db-len>
          <Statistics_hsp-len>7</Statistics_hsp-len>
          <Statistics_eff-space>156650</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
      <Iteration_message>No hits found</Iteration_message>
    </Iteration>"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('random_s00'))

    def test_blastxml_2226_multiple_middle(self):
        """Test blast-xml raw string retrieval, BLAST 2.2.26+, multiple queries, middle (xml_2226_blastp_001.xml)"""
        filename = 'Blast/xml_2226_blastp_001.xml'
        idx = SearchIO.index(filename, self.fmt)
        raw = """<Iteration>
      <Iteration_iter-num>2</Iteration_iter-num>
      <Iteration_query-ID>Query_2</Iteration_query-ID>
      <Iteration_query-def>gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]</Iteration_query-def>
      <Iteration_query-len>102</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|1</Hit_id>
          <Hit_def>gi|308175296|ref|YP_003922001.1| membrane bound lipoprotein [Bacillus amyloliquefaciens DSM 7]</Hit_def>
          <Hit_accession>1</Hit_accession>
          <Hit_len>100</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>139.428</Hsp_bit-score>
              <Hsp_score>350</Hsp_score>
              <Hsp_evalue>1.99275e-46</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>102</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>69</Hsp_identity>
              <Hsp_positive>81</Hsp_positive>
              <Hsp_gaps>2</Hsp_gaps>
              <Hsp_align-len>102</Hsp_align-len>
              <Hsp_qseq>MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN</Hsp_qseq>
              <Hsp_hseq>MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPVSFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN</Hsp_hseq>
              <Hsp_midline>MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|2</Hit_id>
          <Hit_def>gi|375363999|ref|YP_005132038.1| lytA gene product [Bacillus amyloliquefaciens subsp. plantarum CAU B946]</Hit_def>
          <Hit_accession>2</Hit_accession>
          <Hit_len>105</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>88.9669</Hsp_bit-score>
              <Hsp_score>219</Hsp_score>
              <Hsp_evalue>6.94052e-27</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>101</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>104</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>48</Hsp_identity>
              <Hsp_positive>69</Hsp_positive>
              <Hsp_gaps>5</Hsp_gaps>
              <Hsp_align-len>105</Hsp_align-len>
              <Hsp_qseq>MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERA</Hsp_qseq>
              <Hsp_hseq>MKKTIAASFLILLFSVVLAACGTAEQSKKGSG-SSENQAQKETAYYVGMADTHTIEVKVDDQPVSFEFSDDFSDVLNKFSENDKVSITYFTNDKGQKEIKEIEKA</Hsp_hseq>
              <Hsp_midline>MKK IA  F ILL    L+ CG   Q  +G   S ++  + +   YVG+ADTHTIEV VD++PVS + +++ +  L+KF+  DKV+ITY  ND+GQ  +K+IE+A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>3</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|3</Hit_id>
          <Hit_def>gi|154687679|ref|YP_001422840.1| LytA [Bacillus amyloliquefaciens FZB42]</Hit_def>
          <Hit_accession>3</Hit_accession>
          <Hit_len>105</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>88.9669</Hsp_bit-score>
              <Hsp_score>219</Hsp_score>
              <Hsp_evalue>8.41012e-27</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>101</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>104</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>48</Hsp_identity>
              <Hsp_positive>69</Hsp_positive>
              <Hsp_gaps>5</Hsp_gaps>
              <Hsp_align-len>105</Hsp_align-len>
              <Hsp_qseq>MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERA</Hsp_qseq>
              <Hsp_hseq>MKKTIAASFLILLFSVVLAACGTADQSKKGSG-SSENQAQKETAYYVGMADTHTIEVKVDDQPVSFEFSDDFSDVLNKFSENDKVSITYFTNDKGQKEIKEIEKA</Hsp_hseq>
              <Hsp_midline>MKK IA  F ILL    L+ CG   Q  +G   S ++  + +   YVG+ADTHTIEV VD++PVS + +++ +  L+KF+  DKV+ITY  ND+GQ  +K+IE+A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>4</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|4</Hit_id>
          <Hit_def>gi|311070071|ref|YP_003974994.1| unnamed protein product [Bacillus atrophaeus 1942]</Hit_def>
          <Hit_accession>4</Hit_accession>
          <Hit_len>105</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>83.1889</Hsp_bit-score>
              <Hsp_score>204</Hsp_score>
              <Hsp_evalue>1.37847e-24</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>103</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>45</Hsp_identity>
              <Hsp_positive>66</Hsp_positive>
              <Hsp_gaps>5</Hsp_gaps>
              <Hsp_align-len>104</Hsp_align-len>
              <Hsp_qseq>MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIER</Hsp_qseq>
              <Hsp_hseq>MKKNVASSFLILLFSIILAACGTAEQSKEG-NGSSSSQVQNETAYYVGMADTHTIEVKIDDQPVSFEFTDDFSEILNEFEENDKVNISYLTNDKGQKELTEIEK</Hsp_hseq>
              <Hsp_midline>MKK +A  F ILL    L+ CG   Q  +G + S  S ++ +   YVG+ADTHTIEV +D++PVS + T++ +  L++F   DKV I+Y  ND+GQ  L +IE+</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>5</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|15</Hit_id>
          <Hit_def>gi|332258565|ref|XP_003278367.1| PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]</Hit_def>
          <Hit_accession>15</Hit_accession>
          <Hit_len>132</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>15.779</Hsp_bit-score>
              <Hsp_score>29</Hsp_score>
              <Hsp_evalue>7.12269</Hsp_evalue>
              <Hsp_query-from>60</Hsp_query-from>
              <Hsp_query-to>84</Hsp_query-to>
              <Hsp_hit-from>80</Hsp_hit-from>
              <Hsp_hit-to>104</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>7</Hsp_identity>
              <Hsp_positive>11</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>25</Hsp_align-len>
              <Hsp_qseq>VSLDITEESTSDLDKFNSGDKVTIT</Hsp_qseq>
              <Hsp_hseq>VEMGFLHVGQAGLELVTSGDPPTLT</Hsp_hseq>
              <Hsp_midline>V +       + L+   SGD  T+T</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>20</Statistics_db-num>
          <Statistics_db-len>6406</Statistics_db-len>
          <Statistics_hsp-len>38</Statistics_hsp-len>
          <Statistics_eff-space>361344</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|16080617|ref|NP_391444.1|'))

    def test_blastxml_2226_multiple_last(self):
        """Test blast-xml raw string retrieval, BLAST 2.2.26+, multiple queries, last (xml_2226_blastp_001.xml)"""
        filename = 'Blast/xml_2226_blastp_001.xml'
        idx = SearchIO.index(filename, self.fmt)
        raw = """<Iteration>
      <Iteration_iter-num>3</Iteration_iter-num>
      <Iteration_query-ID>Query_3</Iteration_query-ID>
      <Iteration_query-def>gi|11464971:4-101 pleckstrin [Mus musculus]</Iteration_query-def>
      <Iteration_query-len>98</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|5</Hit_id>
          <Hit_def>gi|11464971|ref|NP_062422.1| pleckstrin [Mus musculus]</Hit_def>
          <Hit_accession>5</Hit_accession>
          <Hit_len>350</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>205.682</Hsp_bit-score>
              <Hsp_score>522</Hsp_score>
              <Hsp_evalue>2.24956e-69</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>98</Hsp_identity>
              <Hsp_positive>98</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>43.5134</Hsp_bit-score>
              <Hsp_score>101</Hsp_score>
              <Hsp_evalue>2.90061e-09</Hsp_evalue>
              <Hsp_query-from>3</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>246</Hsp_hit-from>
              <Hsp_hit-to>345</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>29</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>6</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA</Hsp_hseq>
              <Hsp_midline>I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|6</Hit_id>
          <Hit_def>gi|354480464|ref|XP_003502426.1| PREDICTED: pleckstrin-like [Cricetulus griseus]</Hit_def>
          <Hit_accession>6</Hit_accession>
          <Hit_len>350</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>205.297</Hsp_bit-score>
              <Hsp_score>521</Hsp_score>
              <Hsp_evalue>3.2078e-69</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>98</Hsp_identity>
              <Hsp_positive>98</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>43.8986</Hsp_bit-score>
              <Hsp_score>102</Hsp_score>
              <Hsp_evalue>1.81272e-09</Hsp_evalue>
              <Hsp_query-from>3</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>246</Hsp_hit-from>
              <Hsp_hit-to>345</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>30</Hsp_identity>
              <Hsp_positive>50</Hsp_positive>
              <Hsp_gaps>6</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDF-GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAIHLRGCVVTSVESNHDGKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA</Hsp_hseq>
              <Hsp_midline>I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +  GK+     + +I T  +  ++ QAA  +ER  W++ I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>3</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|7</Hit_id>
          <Hit_def>gi|156616273|ref|NP_002655.2| pleckstrin [Homo sapiens]</Hit_def>
          <Hit_accession>7</Hit_accession>
          <Hit_len>350</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>204.142</Hsp_bit-score>
              <Hsp_score>518</Hsp_score>
              <Hsp_evalue>1.081e-68</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>97</Hsp_identity>
              <Hsp_positive>97</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>47.3654</Hsp_bit-score>
              <Hsp_score>111</Hsp_score>
              <Hsp_evalue>1.50729e-10</Hsp_evalue>
              <Hsp_query-from>3</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>246</Hsp_hit-from>
              <Hsp_hit-to>345</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>31</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>6</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMF----VLKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGAEDPLGAIHLRGCVVTSVESNSNGRKSEEENLFEIITADEVHYFLQAATPKERTEWIRAIQMA</Hsp_hseq>
              <Hsp_midline>I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +   R      + +I T  +  +F QAA  +ER  W+R I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>4</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|8</Hit_id>
          <Hit_def>gi|297667453|ref|XP_002811995.1| PREDICTED: pleckstrin-like [Pongo abelii]</Hit_def>
          <Hit_accession>8</Hit_accession>
          <Hit_len>350</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>204.142</Hsp_bit-score>
              <Hsp_score>518</Hsp_score>
              <Hsp_evalue>1.10449e-68</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>97</Hsp_identity>
              <Hsp_positive>97</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>45.4394</Hsp_bit-score>
              <Hsp_score>106</Hsp_score>
              <Hsp_evalue>6.1425e-10</Hsp_evalue>
              <Hsp_query-from>3</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>246</Hsp_hit-from>
              <Hsp_hit-to>345</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>30</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>6</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMF----VLKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGAEDPLGAIHLRGCVVTSVESNSNGRKSEEENLFEIITADEVHYFLQAATPKERTEWIKAIQVA</Hsp_hseq>
              <Hsp_midline>I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +   R      + +I T  +  +F QAA  +ER  W++ I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>5</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|9</Hit_id>
          <Hit_def>gi|350596020|ref|XP_003360649.2| PREDICTED: pleckstrin-like [Sus scrofa]</Hit_def>
          <Hit_accession>9</Hit_accession>
          <Hit_len>228</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>199.519</Hsp_bit-score>
              <Hsp_score>506</Hsp_score>
              <Hsp_evalue>1.97058e-68</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>94</Hsp_identity>
              <Hsp_positive>96</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>32.3426</Hsp_bit-score>
              <Hsp_score>72</Hsp_score>
              <Hsp_evalue>1.12281e-05</Hsp_evalue>
              <Hsp_query-from>30</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>153</Hsp_hit-from>
              <Hsp_hit-to>223</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>21</Hsp_identity>
              <Hsp_positive>32</Hsp_positive>
              <Hsp_gaps>4</Hsp_gaps>
              <Hsp_align-len>71</Hsp_align-len>
              <Hsp_qseq>IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFV---LKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA</Hsp_hseq>
              <Hsp_midline>+ +Y       P G I L+G  +TS      GK  F+       T  +  +F QAA  +ER  W++ I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>20</Statistics_db-num>
          <Statistics_db-len>6406</Statistics_db-len>
          <Statistics_hsp-len>37</Statistics_hsp-len>
          <Statistics_eff-space>345626</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|11464971:4-101'))

    def test_blastxml_2226_single(self):
        """Test blast-xml raw string retrieval, BLAST 2.2.26+, single query (xml_2226_blastp_004.xml)"""
        filename = 'Blast/xml_2226_blastp_004.xml'
        idx = SearchIO.index(filename, self.fmt)
        raw = """<Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>gi|11464971:4-101 pleckstrin [Mus musculus]</Iteration_query-def>
      <Iteration_query-len>98</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|5</Hit_id>
          <Hit_def>gi|11464971|ref|NP_062422.1| pleckstrin [Mus musculus]</Hit_def>
          <Hit_accession>5</Hit_accession>
          <Hit_len>350</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>205.682</Hsp_bit-score>
              <Hsp_score>522</Hsp_score>
              <Hsp_evalue>2.24956e-69</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>98</Hsp_identity>
              <Hsp_positive>98</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>43.5134</Hsp_bit-score>
              <Hsp_score>101</Hsp_score>
              <Hsp_evalue>2.90061e-09</Hsp_evalue>
              <Hsp_query-from>3</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>246</Hsp_hit-from>
              <Hsp_hit-to>345</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>29</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>6</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTS--PCQDFGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAVHLRGCVVTSVESSHDVKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA</Hsp_hseq>
              <Hsp_midline>I++G L+K+G     WK    +L ED   + +Y       P G + L+G  +TS     D  K     + +I T  +  ++ QAA  +ER  W++ I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|6</Hit_id>
          <Hit_def>gi|354480464|ref|XP_003502426.1| PREDICTED: pleckstrin-like [Cricetulus griseus]</Hit_def>
          <Hit_accession>6</Hit_accession>
          <Hit_len>350</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>205.297</Hsp_bit-score>
              <Hsp_score>521</Hsp_score>
              <Hsp_evalue>3.2078e-69</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>98</Hsp_identity>
              <Hsp_positive>98</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>43.8986</Hsp_bit-score>
              <Hsp_score>102</Hsp_score>
              <Hsp_evalue>1.81272e-09</Hsp_evalue>
              <Hsp_query-from>3</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>246</Hsp_hit-from>
              <Hsp_hit-to>345</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>30</Hsp_identity>
              <Hsp_positive>50</Hsp_positive>
              <Hsp_gaps>6</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDF-GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAIHLRGCVVTSVESNHDGKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA</Hsp_hseq>
              <Hsp_midline>I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +  GK+     + +I T  +  ++ QAA  +ER  W++ I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>3</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|7</Hit_id>
          <Hit_def>gi|156616273|ref|NP_002655.2| pleckstrin [Homo sapiens]</Hit_def>
          <Hit_accession>7</Hit_accession>
          <Hit_len>350</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>204.142</Hsp_bit-score>
              <Hsp_score>518</Hsp_score>
              <Hsp_evalue>1.081e-68</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>97</Hsp_identity>
              <Hsp_positive>97</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>47.3654</Hsp_bit-score>
              <Hsp_score>111</Hsp_score>
              <Hsp_evalue>1.50729e-10</Hsp_evalue>
              <Hsp_query-from>3</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>246</Hsp_hit-from>
              <Hsp_hit-to>345</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>31</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>6</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMF----VLKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGAEDPLGAIHLRGCVVTSVESNSNGRKSEEENLFEIITADEVHYFLQAATPKERTEWIRAIQMA</Hsp_hseq>
              <Hsp_midline>I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +   R      + +I T  +  +F QAA  +ER  W+R I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>4</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|8</Hit_id>
          <Hit_def>gi|297667453|ref|XP_002811995.1| PREDICTED: pleckstrin-like [Pongo abelii]</Hit_def>
          <Hit_accession>8</Hit_accession>
          <Hit_len>350</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>204.142</Hsp_bit-score>
              <Hsp_score>518</Hsp_score>
              <Hsp_evalue>1.10449e-68</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>97</Hsp_identity>
              <Hsp_positive>97</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>45.4394</Hsp_bit-score>
              <Hsp_score>106</Hsp_score>
              <Hsp_evalue>6.1425e-10</Hsp_evalue>
              <Hsp_query-from>3</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>246</Hsp_hit-from>
              <Hsp_hit-to>345</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>30</Hsp_identity>
              <Hsp_positive>48</Hsp_positive>
              <Hsp_gaps>6</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMF----VLKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGAEDPLGAIHLRGCVVTSVESNSNGRKSEEENLFEIITADEVHYFLQAATPKERTEWIKAIQVA</Hsp_hseq>
              <Hsp_midline>I++G L+K+G     WK    +L ED   + +Y       P G I L+G  +TS   +   R      + +I T  +  +F QAA  +ER  W++ I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>5</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|9</Hit_id>
          <Hit_def>gi|350596020|ref|XP_003360649.2| PREDICTED: pleckstrin-like [Sus scrofa]</Hit_def>
          <Hit_accession>9</Hit_accession>
          <Hit_len>228</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>199.519</Hsp_bit-score>
              <Hsp_score>506</Hsp_score>
              <Hsp_evalue>1.97058e-68</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>98</Hsp_query-to>
              <Hsp_hit-from>4</Hsp_hit-from>
              <Hsp_hit-to>101</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>94</Hsp_identity>
              <Hsp_positive>96</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>98</Hsp_align-len>
              <Hsp_qseq>KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK</Hsp_qseq>
              <Hsp_hseq>KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK</Hsp_hseq>
              <Hsp_midline>KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK</Hsp_midline>
            </Hsp>
            <Hsp>
              <Hsp_num>2</Hsp_num>
              <Hsp_bit-score>32.3426</Hsp_bit-score>
              <Hsp_score>72</Hsp_score>
              <Hsp_evalue>1.12281e-05</Hsp_evalue>
              <Hsp_query-from>30</Hsp_query-from>
              <Hsp_query-to>96</Hsp_query-to>
              <Hsp_hit-from>153</Hsp_hit-from>
              <Hsp_hit-to>223</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>21</Hsp_identity>
              <Hsp_positive>32</Hsp_positive>
              <Hsp_gaps>4</Hsp_gaps>
              <Hsp_align-len>71</Hsp_align-len>
              <Hsp_qseq>IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFV---LKITTTKQQDHFFQAAFLEERDAWVRDIKKA</Hsp_qseq>
              <Hsp_hseq>LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA</Hsp_hseq>
              <Hsp_midline>+ +Y       P G I L+G  +TS      GK  F+       T  +  +F QAA  +ER  W++ I+ A</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>20</Statistics_db-num>
          <Statistics_db-len>6406</Statistics_db-len>
          <Statistics_hsp-len>37</Statistics_hsp-len>
          <Statistics_eff-space>345626</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|11464971:4-101'))


class BlastTabRawCases(unittest.TestCase):

    fmt = 'blast-tab'

    def test_blasttab_2226_multiple_first(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, first (tab_2226_tblastn_001.txt)"""
        filename = 'Blast/tab_2226_tblastn_001.txt'
        idx = SearchIO.index(filename, self.fmt)
        raw = """gi|16080617|ref|NP_391444.1|	gi|145479850|ref|XM_001425911.1|	34.88	43	28	0	31	73	1744	1872	1e-05	34.7
gi|16080617|ref|NP_391444.1|	gi|72012412|ref|XM_777959.1|	33.90	59	31	1	44	94	1057	1233	1e-04	31.6
gi|16080617|ref|NP_391444.1|	gi|115975252|ref|XM_001180111.1|	33.90	59	31	1	44	94	1057	1233	1e-04	31.6
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|16080617|ref|NP_391444.1|'))

    def test_blasttab_2226_multiple_last(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, last (tab_2226_tblastn_001.txt)"""
        filename = 'Blast/tab_2226_tblastn_001.txt'
        idx = SearchIO.index(filename, self.fmt)
        raw = """gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	95.92	98	4	0	1	98	95	388	2e-67	 199
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	29.58	71	46	2	30	96	542	754	4e-05	32.7
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	97.96	98	2	0	1	98	78	371	2e-67	 202
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	30.00	100	64	2	3	96	804	1103	3e-09	45.1
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	97.96	98	2	0	1	98	161	454	4e-67	 202
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	30.00	100	64	2	3	96	866	1165	3e-09	45.1
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	97.96	98	2	0	1	98	173	466	2e-66	 202
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	31.00	100	63	2	3	96	899	1198	1e-09	46.6
gi|11464971:4-101	gi|365982352|ref|XM_003667962.1|	30.77	52	27	1	12	54	3181	3336	1.7	19.6
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|11464971:4-101'))

    def test_blasttab_2226_single(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, single query (tab_2226_tblastn_004.txt)"""
        filename = 'Blast/tab_2226_tblastn_004.txt'
        idx = SearchIO.index(filename, self.fmt)
        raw = """gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	95.92	98	4	0	1	98	95	388	2e-67	 199
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	29.58	71	46	2	30	96	542	754	4e-05	32.7
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	97.96	98	2	0	1	98	78	371	2e-67	 202
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	30.00	100	64	2	3	96	804	1103	3e-09	45.1
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	97.96	98	2	0	1	98	161	454	4e-67	 202
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	30.00	100	64	2	3	96	866	1165	3e-09	45.1
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	97.96	98	2	0	1	98	173	466	2e-66	 202
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	31.00	100	63	2	3	96	899	1198	1e-09	46.6
gi|11464971:4-101	gi|365982352|ref|XM_003667962.1|	30.77	52	27	1	12	54	3181	3336	1.7	19.6
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|11464971:4-101'))

    def test_blasttab_2226_multiple_first_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, first, commented (tab_2226_tblastn_005.txt)"""
        filename = 'Blast/tab_2226_tblastn_005.txt'
        idx = SearchIO.index(filename, self.fmt, comments=True)
        raw = """# TBLASTN 2.2.26+
# Query: random_s00
# Database: db/minirefseq_mrna
# 0 hits found
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('random_s00'))

    def test_blasttab_2226_multiple_middle_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, middle, commented (tab_2226_tblastn_005.txt)"""
        filename = 'Blast/tab_2226_tblastn_005.txt'
        idx = SearchIO.index(filename, self.fmt, comments=True)
        raw = """# TBLASTN 2.2.26+
# Query: gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]
# Database: db/minirefseq_mrna
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 3 hits found
gi|16080617|ref|NP_391444.1|	gi|145479850|ref|XM_001425911.1|	34.88	43	28	0	31	73	1744	1872	1e-05	34.7
gi|16080617|ref|NP_391444.1|	gi|72012412|ref|XM_777959.1|	33.90	59	31	1	44	94	1057	1233	1e-04	31.6
gi|16080617|ref|NP_391444.1|	gi|115975252|ref|XM_001180111.1|	33.90	59	31	1	44	94	1057	1233	1e-04	31.6
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|16080617|ref|NP_391444.1|'))

    def test_blasttab_2226_multiple_last_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, last, commented (tab_2226_tblastn_005.txt)"""
        filename = 'Blast/tab_2226_tblastn_005.txt'
        idx = SearchIO.index(filename, self.fmt, comments=True)
        raw = """# TBLASTN 2.2.26+
# Query: gi|11464971:4-101 pleckstrin [Mus musculus]
# Database: db/minirefseq_mrna
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 9 hits found
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	95.92	98	4	0	1	98	95	388	2e-67	 199
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	29.58	71	46	2	30	96	542	754	4e-05	32.7
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	97.96	98	2	0	1	98	78	371	2e-67	 202
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	30.00	100	64	2	3	96	804	1103	3e-09	45.1
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	97.96	98	2	0	1	98	161	454	4e-67	 202
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	30.00	100	64	2	3	96	866	1165	3e-09	45.1
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	97.96	98	2	0	1	98	173	466	2e-66	 202
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	31.00	100	63	2	3	96	899	1198	1e-09	46.6
gi|11464971:4-101	gi|365982352|ref|XM_003667962.1|	30.77	52	27	1	12	54	3181	3336	1.7	19.6
# BLAST processed 3 queries
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|11464971:4-101'))

    def test_blasttab_2226_single_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, single query, commented (tab_2226_tblastn_008.txt)"""
        filename = 'Blast/tab_2226_tblastn_008.txt'
        idx = SearchIO.index(filename, self.fmt, comments=True)
        raw = """# TBLASTN 2.2.26+
# Query: gi|11464971:4-101 pleckstrin [Mus musculus]
# Database: db/minirefseq_mrna
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 9 hits found
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	95.92	98	4	0	1	98	95	388	2e-67	 199
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	29.58	71	46	2	30	96	542	754	4e-05	32.7
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	97.96	98	2	0	1	98	78	371	2e-67	 202
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	30.00	100	64	2	3	96	804	1103	3e-09	45.1
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	97.96	98	2	0	1	98	161	454	4e-67	 202
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	30.00	100	64	2	3	96	866	1165	3e-09	45.1
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	97.96	98	2	0	1	98	173	466	2e-66	 202
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	31.00	100	63	2	3	96	899	1198	1e-09	46.6
gi|11464971:4-101	gi|365982352|ref|XM_003667962.1|	30.77	52	27	1	12	54	3181	3336	1.7	19.6
# BLAST processed 1 queries
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|11464971:4-101'))


class HmmerTextRawCases(unittest.TestCase):

    fmt = 'hmmer-text'

    def test_hmmertext_30_multiple_first(self):
        """Test hmmer-text raw string retrieval, HMMER 3.0, multiple queries, first (text_30_hmmscan_001.out)"""
        filename = 'Hmmer/text_30_hmmscan_001.out'
        idx = SearchIO.index(filename, self.fmt)
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
        self.assertEqual(_as_bytes(raw), idx.get_raw('random_s00'))

    def test_hmmertext_30_multiple_middle(self):
        """Test hmmer-text raw string retrieval, HMMER 3.0, multiple queries, middle (text_30_hmmscan_001.out)"""
        filename = 'Hmmer/text_30_hmmscan_001.out'
        idx = SearchIO.index(filename, self.fmt)
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
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|4885477|ref|NP_005359.1|'))

    def test_hmmertext_30_multiple_last(self):
        """Test hmmer-text raw string retrieval, HMMER 3.0, multiple queries, last (text_30_hmmscan_001.out)"""
        filename = 'Hmmer/text_30_hmmscan_001.out'
        idx = SearchIO.index(filename, self.fmt)
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
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|125490392|ref|NP_038661.2|'))

    def test_hmmertext_30_single(self):
        """Test hmmer-text raw string retrieval, HMMER 3.0, single query (text_30_hmmscan_003.out)"""
        filename = 'Hmmer/text_30_hmmscan_003.out'
        idx = SearchIO.index(filename, self.fmt)
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
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|4885477|ref|NP_005359.1|'))


class HmmerTabRawCases(unittest.TestCase):

    fmt = 'hmmer-tab'

    def test_hmmertab_30_multiple_first(self):
        """Test hmmer-tab raw string retrieval, HMMER 3.0, multiple queries, first (tab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_001.out')
        idx = SearchIO.index(filename, self.fmt)
        raw = """Globin               PF00042.17 gi|4885477|ref|NP_005359.1| -              6e-21   74.6   0.3   9.2e-21   74.0   0.2   1.3   1   0   0   1   1   1   1 Globin
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|4885477|ref|NP_005359.1|'))

    def test_hmmertab_30_multiple_middle(self):
        """Test hmmer-tab raw string retrieval, HMMER 3.0, multiple queries, middle (tab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_001.out')
        idx = SearchIO.index(filename, self.fmt)
        raw = """Ig_3                 PF13927.1  gi|126362951:116-221 -            1.4e-09   38.2   0.4   2.1e-09   37.6   0.3   1.3   1   0   0   1   1   1   1 Immunoglobulin domain
Ig_2                 PF13895.1  gi|126362951:116-221 -            3.5e-05   23.7   0.1   4.3e-05   23.4   0.1   1.1   1   0   0   1   1   1   1 Immunoglobulin domain
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|126362951:116-221'))

    def test_hmmertab_30_multiple_last(self):
        """Test hmmer-tab raw string retrieval, HMMER 3.0, multiple queries, last (tab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_001.out')
        idx = SearchIO.index(filename, self.fmt)
        raw = """Pou                  PF00157.12 gi|125490392|ref|NP_038661.2| -              7e-37  124.8   0.5   1.4e-36  123.9   0.3   1.5   1   0   0   1   1   1   1 Pou domain - N-terminal to homeobox domain
Homeobox             PF00046.24 gi|125490392|ref|NP_038661.2| -            2.1e-18   65.5   1.1   4.1e-18   64.6   0.7   1.5   1   0   0   1   1   1   1 Homeobox domain
HTH_31               PF13560.1  gi|125490392|ref|NP_038661.2| -              0.012   15.6   0.0      0.16   12.0   0.0   2.2   2   0   0   2   2   2   0 Helix-turn-helix domain
Homeobox_KN          PF05920.6  gi|125490392|ref|NP_038661.2| -              0.039   13.5   0.0     0.095   12.3   0.0   1.6   1   0   0   1   1   1   0 Homeobox KN domain
DUF521               PF04412.8  gi|125490392|ref|NP_038661.2| -               0.14   10.5   0.1      0.26    9.6   0.1   1.4   1   0   0   1   1   1   0 Protein of unknown function (DUF521)
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|125490392|ref|NP_038661.2|'))

    def test_hmmertab_30_single(self):
        """Test hmmer-tab raw string retrieval, HMMER 3.0, single query (tab_30_hmmscan_004.out)"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_004.out')
        idx = SearchIO.index(filename, self.fmt)
        raw = """Ig_3                 PF13927.1  gi|126362951:116-221 -            1.4e-09   38.2   0.4   2.1e-09   37.6   0.3   1.3   1   0   0   1   1   1   1 Immunoglobulin domain
Ig_2                 PF13895.1  gi|126362951:116-221 -            3.5e-05   23.7   0.1   4.3e-05   23.4   0.1   1.1   1   0   0   1   1   1   1 Immunoglobulin domain
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|126362951:116-221'))


class HmmerDomtabRawCases(unittest.TestCase):

    fmt = 'hmmscan-domtab'

    def test_hmmerdomtab_30_multiple_first(self):
        """Test hmmscan-domtab raw string retrieval, HMMER 3.0, multiple queries, first (domtab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_001.out')
        idx = SearchIO.index(filename, self.fmt)
        raw = """Globin               PF00042.17   108 gi|4885477|ref|NP_005359.1| -            154     6e-21   74.6   0.3   1   1   6.7e-25   9.2e-21   74.0   0.2     1   107     7   112     7   113 0.97 Globin
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|4885477|ref|NP_005359.1|'))

    def test_hmmerdomtab_30_multiple_middle(self):
        """Test hmmscan-domtab raw string retrieval, HMMER 3.0, multiple queries, middle (domtab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_001.out')
        idx = SearchIO.index(filename, self.fmt)
        raw = """Ig_3                 PF13927.1     75 gi|126362951:116-221 -            106   1.4e-09   38.2   0.4   1   1     3e-13   2.1e-09   37.6   0.3     1    73     9    84     9    88 0.94 Immunoglobulin domain
Ig_2                 PF13895.1     80 gi|126362951:116-221 -            106   3.5e-05   23.7   0.1   1   1   6.2e-09   4.3e-05   23.4   0.1     1    80     9   104     9   104 0.71 Immunoglobulin domain
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|126362951:116-221'))

    def test_hmmerdomtab_30_multiple_last(self):
        """Test hmmscan-domtab raw string retrieval, HMMER 3.0, multiple queries, last (domtab_30_hmmscan_001.out)"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_001.out')
        idx = SearchIO.index(filename, self.fmt)
        raw = """Pou                  PF00157.12    75 gi|125490392|ref|NP_038661.2| -            352     7e-37  124.8   0.5   1   1     5e-40   1.4e-36  123.9   0.3     3    75   133   205   131   205 0.97 Pou domain - N-terminal to homeobox domain
Homeobox             PF00046.24    57 gi|125490392|ref|NP_038661.2| -            352   2.1e-18   65.5   1.1   1   1   1.5e-21   4.1e-18   64.6   0.7     1    57   224   280   224   280 0.98 Homeobox domain
HTH_31               PF13560.1     64 gi|125490392|ref|NP_038661.2| -            352     0.012   15.6   0.0   1   2   5.7e-05      0.16   12.0   0.0     1    35   141   181   141   184 0.96 Helix-turn-helix domain
HTH_31               PF13560.1     64 gi|125490392|ref|NP_038661.2| -            352     0.012   15.6   0.0   2   2      0.19   5.2e+02    0.8   0.0    39    62   245   268   243   270 0.86 Helix-turn-helix domain
Homeobox_KN          PF05920.6     40 gi|125490392|ref|NP_038661.2| -            352     0.039   13.5   0.0   1   1   3.5e-05     0.095   12.3   0.0     7    39   244   276   241   277 0.91 Homeobox KN domain
DUF521               PF04412.8    400 gi|125490392|ref|NP_038661.2| -            352      0.14   10.5   0.1   1   1   9.4e-05      0.26    9.6   0.1   273   334   221   280   197   294 0.77 Protein of unknown function (DUF521)
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|125490392|ref|NP_038661.2|'))

    def test_hmmerdomtab_30_single(self):
        """Test hmmscan-domtab raw string retrieval, HMMER 3.0, single query (domtab_30_hmmscan_004.out)"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_004.out')
        idx = SearchIO.index(filename, self.fmt)
        raw = """Ig_3                 PF13927.1     75 gi|126362951:116-221 -            106   1.4e-09   38.2   0.4   1   1     3e-13   2.1e-09   37.6   0.3     1    73     9    84     9    88 0.94 Immunoglobulin domain
Ig_2                 PF13895.1     80 gi|126362951:116-221 -            106   3.5e-05   23.7   0.1   1   1   6.2e-09   4.3e-05   23.4   0.1     1    80     9   104     9   104 0.71 Immunoglobulin domain
"""
        self.assertEqual(_as_bytes(raw), idx.get_raw('gi|126362951:116-221'))


class SearchIndexCases(unittest.TestCase):

    def check_index(self, filename, format, **kwargs):
        # check if Python3 installation has sqlite3
        try:
            import sqlite3
        except ImportError:
            sqlite3 = None

        parsed = list(SearchIO.parse(filename, format, **kwargs))
        # compare values by index
        indexed = SearchIO.index(filename, format, **kwargs)
        self.assertEqual(len(parsed), len(indexed.keys()))
        # compare values by index_db, only if sqlite3 is present
        if sqlite3 is not None:
            db_indexed = SearchIO.index_db(':memory:', [filename], format, **kwargs)
            self.assertEqual(len(parsed), len(db_indexed.keys()))

        for qres in parsed:
            idx_qres = indexed[qres.id]
            # parsed and indexed qresult are different objects!
            self.assertNotEqual(id(qres), id(idx_qres))
            # but they should have the same attribute values
            self.assertTrue(compare_search_obj(qres, idx_qres))
            # sqlite3 comparison, only if it's present
            if sqlite3 is not None:
                dbidx_qres = db_indexed[qres.id]
                self.assertNotEqual(id(qres), id(dbidx_qres))
                self.assertTrue(compare_search_obj(qres, dbidx_qres))


class BlastXmlIndexCases(SearchIndexCases):

    fmt = 'blast-xml'

    def test_blastxml_2212L_blastp_001(self):
        """Test blast-xml indexing, BLAST 2.2.12"""
        filename = 'Blast/xml_2212L_blastp_001.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_2218_blastp_001(self):
        """Test blast-xml indexing, BLAST 2.2.18+"""
        filename = 'Blast/xml_2218_blastp_001.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_2222_blastx_001(self):
        """Test blast-xml indexing, BLAST 2.2.22+"""
        filename = 'Blast/xml_2222_blastx_001.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_2226_tblastn_001(self):
        """Test blast-xml indexing, BLAST 2.2.26+, multiple queries"""
        filename = 'Blast/xml_2226_tblastn_001.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_2226_tblastn_002(self):
        """Test blast-xml indexing, BlAST 2.2.26+, single query, no hits"""
        filename = 'Blast/xml_2226_tblastn_002.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_2226_tblastn_004(self):
        """Test blast-xml indexing, BLAST 2.2.26+, single query, multiple hits"""
        filename = 'Blast/xml_2226_tblastn_004.xml'
        self.check_index(filename, self.fmt)


class BlastTabIndexCases(SearchIndexCases):

    fmt = 'blast-tab'

    def test_blasttab_2226_tblastn_001(self):
        """Test blast-tab indexing, BLAST 2.2.26+, multiple queries"""
        filename = 'Blast/tab_2226_tblastn_001.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_2226_tblastn_002(self):
        """Test blast-tab indexing, BLAST 2.2.26+, single query, no hits"""
        filename = 'Blast/tab_2226_tblastn_002.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_2226_tblastn_004(self):
        """Test blast-tab indexing, BLAST 2.2.26+, single query, multiple hits"""
        filename = 'Blast/tab_2226_tblastn_004.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_2226_tblastn_005(self):
        """Test blast-tab indexing, BLAST 2.2.26+, multiple queries, commented"""
        filename = 'Blast/tab_2226_tblastn_005.txt'
        self.check_index(filename, self.fmt, comments=True)

    def test_blasttab_2226_tblastn_006(self):
        """Test blast-tab indexing, BLAST 2.2.26+, single query, no hits, commented"""
        filename = 'Blast/tab_2226_tblastn_006.txt'
        self.check_index(filename, self.fmt, comments=True)

    def test_blasttab_comment_sing(self):
        """Test blast-tab indexing, BLAST 2.2.26+, single query, multiple hits, commented"""
        filename = 'Blast/tab_2226_tblastn_008.txt'
        self.check_index(filename, self.fmt, comments=True)

    def test_blasttab_2226_tblastn_009(self):
        """Test blast-tab indexing, BLAST 2.2.26+, custom columns"""
        filename = 'Blast/tab_2226_tblastn_009.txt'
        self.check_index(filename, self.fmt, fields=['qseqid', 'sseqid'])

    def test_blasttab_2226_tblastn_010(self):
        """Test blast-tab indexing, BLAST 2.2.26+, custom columns, commented"""
        filename = 'Blast/tab_2226_tblastn_010.txt'
        self.check_index(filename, self.fmt, comments=True)

    def test_blasttab_2226_tblastn_011(self):
        """Test blast-tab indexing, BLAST 2.2.26+, all columns, commented"""
        filename = 'Blast/tab_2226_tblastn_011.txt'
        self.check_index(filename, self.fmt, comments=True)


class HmmerTextIndexCases(SearchIndexCases):

    fmt = 'hmmer-text'

    def test_hmmertext_text_30_hmmscan_001(self):
        """Test hmmer-text indexing, HMMER 3.0, multiple queries"""
        filename = 'Hmmer/text_30_hmmscan_001.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmscan_002(self):
        """Test hmmer-text indexing, HMMER 3.0, single query, no hits"""
        filename = 'Hmmer/text_30_hmmscan_002.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmscan_006(self):
        """Test hmmer-text indexing, HMMER 3.0, single query, multiple hits"""
        filename = 'Hmmer/text_30_hmmscan_006.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmscan_007(self):
        """Test hmmer-text indexing, HMMER 3.0, single query, no alignments"""
        filename = 'Hmmer/text_30_hmmscan_007.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmscan_008(self):
        """Test hmmer-text indexing, HMMER 3.0, single query, no alignment width"""
        filename = 'Hmmer/text_30_hmmscan_008.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_text_30_hmmsearch_005(self):
        """Test hmmer-text indexing, HMMER 3.0, multiple queries"""
        filename = 'Hmmer/text_30_hmmsearch_005.out'
        self.check_index(filename, self.fmt)


class HmmerTabIndexCases(SearchIndexCases):

    fmt = 'hmmer-tab'

    def test_hmmertab_30_hmmscan_001(self):
        """Test hmmer-tab indexing, HMMER 3.0, multiple queries"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_001.out')
        self.check_index(filename, self.fmt)

    def test_hmmertab_30_hmmscan_002(self):
        """Test hmmer-tab indexing, HMMER 3.0, single query, no hits"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_002.out')
        self.check_index(filename, self.fmt)

    def test_hmmertab_30_hmmscan_003(self):
        """Test hmmer-tab indexing, HMMER 3.0, single query, multiple hits"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_003.out')
        self.check_index(filename, self.fmt)

    def test_hmmertab_30_hmmscan_004(self):
        """Test hmmer-tab indexing, HMMER 3.0, single query, no alignments"""
        filename = os.path.join('Hmmer', 'tab_30_hmmscan_004.out')
        self.check_index(filename, self.fmt)


class HmmerDomtabIndexCases(SearchIndexCases):

    def test_hmmerdomtab_30_hmmscan_001(self):
        """Test hmmscan-domtab indexing, HMMER 3.0, multiple queries"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_001.out')
        self.check_index(filename, 'hmmscan-domtab')

    def test_hmmerdomtab_30_hmmscan_002(self):
        """Test hmmscan-domtab indexing, HMMER 3.0, single query, no hits"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_002.out')
        self.check_index(filename, 'hmmscan-domtab')

    def test_hmmerdomtab_30_hmmscan_003(self):
        """Test hmmscan-domtab indexing, HMMER 3.0, single query, multiple hits"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_003.out')
        self.check_index(filename, 'hmmscan-domtab')

    def test_hmmerdomtab_30_hmmscan_004(self):
        """Test hmmscan-domtab indexing, HMMER 3.0, single query, no alignments"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmscan_004.out')
        self.check_index(filename, 'hmmscan-domtab')

    def test_hmmerdomtab_30_hmmsearch_001(self):
        """Test hmmsearch-domtab indexing, HMMER 3.0, single query, no alignments"""
        filename = os.path.join('Hmmer', 'domtab_30_hmmsearch_001.out')
        self.check_index(filename, 'hmmsearch-domtab')


class FastaM10IndexCases(SearchIndexCases):

    fmt = 'fasta-m10'

    def test_output_002(self):
        """Test fasta-m10 indexing, fasta34, multiple queries"""
        filename = os.path.join('Fasta', 'output002.m10')
        self.check_index(filename, self.fmt)

    def test_output_001(self):
        """Test fasta-m10 indexing, fasta35, multiple queries"""
        filename = os.path.join('Fasta', 'output001.m10')
        self.check_index(filename, self.fmt)

    def test_output_005(self):
        """Test fasta-m10 indexing, ssearch35, multiple queries"""
        filename = os.path.join('Fasta', 'output005.m10')
        self.check_index(filename, self.fmt)

    def test_output_008(self):
        """Test fasta-m10 indexing, tfastx36, multiple queries"""
        filename = os.path.join('Fasta', 'output008.m10')
        self.check_index(filename, self.fmt)

    def test_output_009(self):
        """Test fasta-m10 indexing, fasta36, multiple queries"""
        filename = os.path.join('Fasta', 'output009.m10')
        self.check_index(filename, self.fmt)

    def test_output_010(self):
        """Test fasta-m10 indexing, fasta36, single query, no hits"""
        filename = os.path.join('Fasta', 'output010.m10')
        self.check_index(filename, self.fmt)

    def test_output_011(self):
        """Test fasta-m10 indexing, fasta36, single query, hits with single hsp"""
        filename = os.path.join('Fasta', 'output011.m10')
        self.check_index(filename, self.fmt)

    def test_output_012(self):
        """Test fasta-m10 indexing, fasta36, single query with multiple hsps"""
        filename = os.path.join('Fasta', 'output012.m10')
        self.check_index(filename, self.fmt)


class BlatPslIndexCases(SearchIndexCases):

    fmt = 'blat-psl'

    def test_psl_34_001(self):
        """Test blat-psl indexing, multiple queries"""
        filename = os.path.join('Blat', 'psl_34_001.psl')
        self.check_index(filename, self.fmt)

    def test_psl_34_002(self):
        """Test blat-psl indexing, single query, no hits"""
        filename = os.path.join('Blat', 'psl_34_002.psl')
        self.check_index(filename, self.fmt)

    def test_psl_34_003(self):
        """Test blat-psl indexing, single query, single hit"""
        filename = os.path.join('Blat', 'psl_34_003.psl')
        self.check_index(filename, self.fmt)

    def test_psl_34_004(self):
        """Test blat-psl indexing, single query, multiple hits with multiple hsps"""
        filename = os.path.join('Blat', 'psl_34_004.psl')
        self.check_index(filename, self.fmt)


    def test_psl_34_005(self):
        """Test blat-psl indexing, multiple queries, no header"""
        filename = os.path.join('Blat', 'psl_34_005.psl')
        self.check_index(filename, self.fmt)


class BlatPslxIndexCases(SearchIndexCases):

    fmt = 'blat-pslx'

    def test_pslx_34_001(self):
        """Test blat-pslx indexing, multiple queries"""
        filename = os.path.join('Blat', 'pslx_34_001.pslx')
        self.check_index(filename, self.fmt)

    def test_pslx_34_002(self):
        """Test blat-pslx indexing, single query, no hits"""
        filename = os.path.join('Blat', 'pslx_34_002.pslx')
        self.check_index(filename, self.fmt)

    def test_pslx_34_003(self):
        """Test blat-pslx indexing, single query, single hit"""
        filename = os.path.join('Blat', 'pslx_34_003.pslx')
        self.check_index(filename, self.fmt)

    def test_pslx_34_004(self):
        """Test blat-pslx indexing, single query, multiple hits with multiple hsps"""
        filename = os.path.join('Blat', 'pslx_34_004.pslx')
        self.check_index(filename, self.fmt)


    def test_pslx_34_005(self):
        """Test blat-pslx indexing, multiple queries, no header"""
        filename = os.path.join('Blat', 'pslx_34_005.pslx')
        self.check_index(filename, self.fmt)


class ExonerateTextIndexCases(SearchIndexCases):

    fmt = 'exonerate-text'

    def test_exn_22_m_est2genome(self):
        """Test exonerate-text indexing, single"""
        filename = os.path.join('Exonerate', 'exn_22_m_est2genome.exn')
        self.check_index(filename, self.fmt)

    def test_exn_22_q_multiple(self):
        """Test exonerate-text indexing, single"""
        filename = os.path.join('Exonerate', 'exn_22_q_multiple.exn')
        self.check_index(filename, self.fmt)


class ExonerateVulgarIndexCases(SearchIndexCases):

    fmt = 'exonerate-vulgar'

    def test_exn_22_m_est2genome(self):
        """Test exonerate-vulgar indexing, single"""
        filename = os.path.join('Exonerate', 'exn_22_o_vulgar.exn')
        self.check_index(filename, self.fmt)

    def test_exn_22_q_multiple(self):
        """Test exonerate-vulgar indexing, single"""
        filename = os.path.join('Exonerate', 'exn_22_q_multiple_vulgar.exn')
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
