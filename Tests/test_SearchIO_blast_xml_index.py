# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO blast-xml indexing."""

import unittest

from search_tests_common import CheckRaw, CheckIndex


class BlastXmlRawCases(CheckRaw):
    """Check BLAST XML get_raw method."""

    fmt = "blast-xml"

    def test_blastxml_2226_multiple_first(self):
        """Test blast-xml raw string retrieval, BLAST 2.2.26+, multiple queries, first (xml_2226_blastp_001.xml)."""
        filename = "Blast/xml_2226_blastp_001.xml"
        raw = """    <Iteration>
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
    </Iteration>
"""
        self.check_raw(filename, "random_s00", raw)

    def test_blastxml_2226_multiple_middle(self):
        """Test blast-xml raw string retrieval, BLAST 2.2.26+, multiple queries, middle (xml_2226_blastp_001.xml)."""
        filename = "Blast/xml_2226_blastp_001.xml"
        raw = """    <Iteration>
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
    </Iteration>
"""
        self.check_raw(filename, "gi|16080617|ref|NP_391444.1|", raw)

    def test_blastxml_2226_multiple_last(self):
        """Test blast-xml raw string retrieval, BLAST 2.2.26+, multiple queries, last (xml_2226_blastp_001.xml)."""
        filename = "Blast/xml_2226_blastp_001.xml"
        raw = """    <Iteration>
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
    </Iteration>
"""
        self.check_raw(filename, "gi|11464971:4-101", raw)

    def test_blastxml_2226_single(self):
        """Test blast-xml raw string retrieval, BLAST 2.2.26+, single query (xml_2226_blastp_004.xml)."""
        filename = "Blast/xml_2226_blastp_004.xml"
        raw = """    <Iteration>
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
    </Iteration>
"""
        self.check_raw(filename, "gi|11464971:4-101", raw)


class BlastXmlIndexCases(CheckIndex):

    fmt = "blast-xml"

    def test_blastxml_2212L_blastp_001(self):
        """Test blast-xml indexing, BLAST 2.2.12."""
        filename = "Blast/xml_2212L_blastp_001.xml"
        self.check_index(filename, self.fmt)

    def test_blastxml_2218_blastp_001(self):
        """Test blast-xml indexing, BLAST 2.2.18+."""
        filename = "Blast/xml_2218_blastp_001.xml"
        self.check_index(filename, self.fmt)

    def test_blastxml_2222_blastx_001(self):
        """Test blast-xml indexing, BLAST 2.2.22+."""
        filename = "Blast/xml_2222_blastx_001.xml"
        self.check_index(filename, self.fmt)

    def test_blastxml_2226_tblastn_001(self):
        """Test blast-xml indexing, BLAST 2.2.26+, multiple queries."""
        filename = "Blast/xml_2226_tblastn_001.xml"
        self.check_index(filename, self.fmt)

    def test_blastxml_2226_tblastn_002(self):
        """Test blast-xml indexing, BlAST 2.2.26+, single query, no hits."""
        filename = "Blast/xml_2226_tblastn_002.xml"
        self.check_index(filename, self.fmt)

    def test_blastxml_2226_tblastn_004(self):
        """Test blast-xml indexing, BLAST 2.2.26+, single query, multiple hits."""
        filename = "Blast/xml_2226_tblastn_004.xml"
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
