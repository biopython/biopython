BLAST TEST DATA FILES
=====================

This directory contains various data files for testing the
BLAST-related code in Biopython. All files are grouped by BLAST
release version, from the most recent first.


BLAST 2.2.28+
-------------
tab_2228_tblastn_001         single query, custom columns: evalue sallseqid qseqid, commented

tab_2228_tblastx_001         single query with all columns (including new ones added in this version), commented


BLAST 2.2.26+
-------------
tab_2226_tblastn_001         multiple queries
tab_2226_tblastn_002         single query, no hits
tab_2226_tblastn_003         single query, hits @ single hsp
tab_2226_tblastn_004         single query, multiple hsps per hit present
tab_2226_tblastn_005         multiple queries, commented
tab_2226_tblastn_006         single query, no hits, commented
tab_2226_tblastn_007         single query, hits @ single hsp, commented
tab_2226_tblastn_008         single query, multiple hsps per hit present, commented
tab_2226_tblastn_009         multiple queries with custom columns
tab_2226_tblastn_010         multiple queries with custom columns, commented
tab_2226_tblastn_011         multiple queries with all columns, commented
tab_2226_tblastn_012         multiple quries, remote search, standard columns, commented
tab_2226_tblastn_013         single query, custom columns: qsseq std sseq

xml_2226_blastn_001.xml      multiple queries
xml_2226_blastn_002.xml      single query, no hits
xml_2226_blastn_003.xml      single query, hits @ single hsp
xml_2226_blastn_004.xml      single query, multiple hsps per hit present
xml_2226_blastn_005.xml      multiple queries, remote search on NCBI's database
xml_2226_blastn_006.xml      single query, against a database with duplicate entries
xml_2226_blastp_001.xml      multiple queries
xml_2226_blastp_002.xml      single query, no hits
xml_2226_blastp_003.xml      single query, hits @ single hsp
xml_2226_blastp_004.xml      single query, multiple hsps per hit present
xml_2226_blastp_005.xml      multiple queries, remote search on NCBI's database
xml_2226_blastx_001.xml      multiple queries
xml_2226_blastx_002.xml      single query, no hits
xml_2226_blastx_003.xml      single query, multiple hsps per hit present
xml_2226_blastx_004.xml      multiple queries, remote search on NCBI's database
xml_2226_tblastn_001.xml     multiple queries
xml_2226_tblastn_002.xml     single query, no hits
xml_2226_tblastn_003.xml     single query, hits @ single hsp
xml_2226_tblastn_004.xml     single query, multiple hsps per hit present
xml_2226_tblastn_005.xml     multiple queries, remote search on NCBI's database
xml_2226_tblastx_001.xml     multiple queries
xml_2226_tblastx_002.xml     single query, no hits
xml_2226_tblastx_003.xml     single query, multiple hsps per hit present
xml_2226_tblastx_004.xml     multiple queries, remote search on NCBI's database


BLAST 2.2.22+
-------------
xml_2222_blastp_001.xml      n/a
xml_2222_blastx_001.xml      multiple queries


BLAST 2.2.18+
-------------
xml_2218_blastp_001.xml      n/a
xml_2218_blastp_002.xml      n/a


BLAST 2.2.18 (Legacy)
---------------------
xml_2218L_blastp_001.xml     n/a
xml_2218L_rpsblast_001.xml   n/a


BLAST 2.2.12 (Legacy)
---------------------
xml_2212L_blastn_001.xml     n/a
xml_2212L_blastp_001.xml     n/a
xml_2212L_blastx_001.xml     n/a
xml_2212L_tblastn_001.xml    n/a
xml_2212L_tblastx_001.xml    n/a


Mock files for mocked qblast searches
-------------------------------------
This are minimal xml/html files mimicking NCBI qblast responses
for our mocked online tests in test_NCBI_qlast_mock.py:

mock_actin.xml
mock_disco.xml
mock_orchid.xml
mock_pcr.xml
mock_short_empty.xml
mock_short_result.xml
mock_wait.html
