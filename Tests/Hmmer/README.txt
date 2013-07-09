HMMER TEST DATA FILES
=====================

This directory contains various data files for testing the
HMMER-related code in Biopython.


HMMER3 OUTPUT FILES (plain text)
-------------------------------
text_30_hmmscan_001.out      multiple queries
text_30_hmmscan_002.out      single query, no match
text_30_hmmscan_003.out      single query, one match, one hsp per match
text_30_hmmscan_004.out      single query, multiple matches, one hsp per match
text_30_hmmscan_005.out      single query, multiple matches, multiple hsps per match
text_30_hmmscan_006.out      single query, multiple matches, multiple hsps per match, inclusion threshold present
text_30_hmmscan_007.out      single query, one match, no alignment
text_30_hmmscan_008.out      single query, multiple matches, no alignment width
text_30_hmmscan_009.out      single query, alignment block(s) with large gaps (bug 3399 in Redmine)
text_30_hmmscan_010.out      
text_30_hmmsearch_001.out    single query, no match
text_30_hmmsearch_002.out    single query, multiple matches, multiple hsps per match
text_30_hmmsearch_003.out    single query, multiple matches, multiple hsps per match, no alignment
text_30_hmmsearch_004.out    single query, multiple matches, multiple hsps per match, no alignment width
text_30_hmmsearch_005.out    multiple queries

tab_30_hmmscan_001.out      multiple queries
tab_30_hmmscan_002.out      single query, no match
tab_30_hmmscan_003.out      single query, one match, one hsp per match
tab_30_hmmscan_004.out      single query, multiple matches, one hsp per match

domtab_30_hmmscan_001.out   multiple queries, hmm as hit
domtab_30_hmmscan_002.out   single query, no match, hmm as hit
domtab_30_hmmscan_003.out   single query, one match, one hsp per match, hmm as hit
domtab_30_hmmscan_004.out   single query, multiple matches, one hsp per match, hmm as hit
domtab_30_hmmsearch_001.out single query, multiple matches, hmm as query


HMMER2 OUTPUT FILES (plain text)
--------------------------------
text_21_hmmpfam_001.out     single query, two matches, bioperl's hmmpfam.out file
text_22_hmmpfam_001.out     single query, one match, bioperl's L77119.hmmer file
text_23_hmmpfam_001.out     single query, multiple matches, bioperl's hmmpfam_cs.out file
text_23_hmmpfam_002.out     single query, no match
text_23_hmmpfam_003.out     single query, one match, missing some consensus content
text_24_hmmpfam_001.out     multiple queries
text_20_hmmsearch_001.out   single query, multiple matches, bioperl's hmmsearch.out file
text_22_hmmsearch_001.out   single query, multiple matches, bioperl's cysprot1b.hmmsearch file
