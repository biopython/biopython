BLAST TEST DATA FILES
=====================

This directory contains various data files for testing the
BLAST-related code in Biopython.

BLAST OUTPUT FILES (plain text)
-------------------------------

ID    PROGRAM  VER.   QUERY      DATABASE  PARAMETERS     DESCRIPTION
bt001 blastall 2.0.10 sp:P21297  swissprot -p blastp
bt002 blastall 2.0.10 sp:Q02112  pdbaa     -p blastp      no hits
bt003 blastall 2.0.10 sp:P16153  swissprot -p blastp -v 0 no descriptions 
bt004 blastall 2.0.10 sp:P16153  swissprot -p blastp -b 0 no alignments
bt005 blastall 2.0.10 sp:P15394  swissprot -p blastp -m 6 master-slave
bt006 blastpgp 2.0.10 sp:P12921  swissprot -j 10          converges in 1 round
bt007 blastpgp 2.0.10 sp:P17216  swissprot -j 10          converges in 3 rounds
*bt008 blastpgp 2.0.10 sp:P33654  swissprot -j 10          doesn't converge
bt009 blastpgp 2.0.10 sp:Q02134  swissprot -j 10          no new sequences before CONVERGED! message
bt010 blastall 2.0.10 gi:1348916 sts       -p blastn
bt011 blastall 2.0.10 gi:1348400 sts       -p blastn -b 0 no alignment
bt012 blastall 2.0.10 gi:1347201 sts       -p blastn -v 0 no descriptions
bt013 blastall 2.0.10 gi:859351  sts       -p blastn -m 1 master-slave
bt014 blastall 2.0.10 gi:1347369 swissprot -p blastx
bt015 blastall 2.0.10 gi:1347782 swissprot -p blastx      no hits
bt016 blastall 2.0.10 sp:P39483  sts       -p tblastn
bt017 blastall 2.0.10 sp:P19888  sts       -p tblastn     no hits
bt018 blastall 2.0.10 gi:1348853 sts       -p tblastx
bt019 NCBI WWW 2.0.10 sp:P21297  swissprot blastp
bt020 NCBI WWW 2.0.10 sp:P54929  pdb       blastp         no hits
bt021 NCBI WWW 2.0.10 sp:P16153  swissprot blastp descriptions=0
bt022 NCBI WWW 2.0.10 sp:P16153  swissprot blastp alignments=0
bt023 NCBI WWW 2.0.10 sp:P15394  swissprot blastp flat master-slave without identities
bt024 NCBI WWW 2.0.10 gi:1348916 dbsts     blastn
bt025 NCBI WWW 2.0.10 gi:1348400 dbsts     blastn alignments=0
bt026 NCBI WWW 2.0.10 gi:1347201 dbsts     blastn descriptions=0
bt027 NCBI WWW 2.0.10 gi:859351  dbsts     blastn master-slave with identities
bt028 NCBI WWW 2.0.10 gi:1347369 swissprot blastx
bt029 NCBI WWW 2.0.10 gi:1347782 swissprot blastx         no hits
bt030 NCBI WWW 2.0.10 sp:P39483  dbsts     tblastn
bt031 NCBI WWW 2.0.10 sp:P19888  dbsts     tblastn        no hits
bt032 NCBI WWW 2.0.10 gi:1348853 dbsts     tblastx
bt033 NCBI WWW 2.0.10 gi:1593528 dbsts     blastx         >1 HSP in alignment
bt034 NCBI WWW 2.0.10 gi:1347369 swissprot blastx         no header info
bt035 NCBI WWW 2.0.10 sp:P23832  swissprot blastp         >1 HSP in alignment
bt036 NCBI WWW 2.0.10 sp:P29247  swissprot blastp no graphical overview
bt037 NCBI WWW 2.0.10 sp:Q06110  swissprot blastp no graphical overview descriptions=0
bt038 NCBI WWW 2.0.10 sp:Q06110  swissprot blastp descriptions=0 master-slave without identities
bt039 blastall 2.0.10 sp:Q08386  swissprot -p blastp      >1 HSP in alignment
bt040 blastall 2.0.10 sp:Q05362  swissprot -p blastp -v 0 -m 4 no descriptions, master-slave
bt041 blastall 2.0.11 sp:P21297  swissprot -p blastp
bt042 blastall 2.0.11 sp:Q02112  pdbaa     -p blastp      no hits
bt043 blastall 2.0.11 sp:P16153  swissprot -p blastp -v 0 no descriptions 
bt044 blastall 2.0.11 sp:P16153  swissprot -p blastp -b 0 no alignments
bt045 blastall 2.0.11 sp:P15394  swissprot -p blastp -m 6 master-slave
bt046 blastpgp 2.0.11 sp:P12921  swissprot -j 10          converges in 2 rounds
bt047 blastpgp 2.0.11 sp:Q02134  swissprot -j 10          no new sequences before CONVERGED! message
bt048 blastall 2.0.11 gi:1348916 sts       -p blastn
bt049 blastall 2.0.11 gi:1348400 sts       -p blastn -b 0 no alignment
bt050 blastall 2.0.11 gi:1347201 sts       -p blastn -v 0 no descriptions
bt051 blastall 2.0.11 gi:859351  sts       -p blastn -m 1 master-slave
bt052 blastall 2.0.11 gi:1347369 swissprot -p blastx
bt053 blastall 2.0.11 gi:1347782 swissprot -p blastx      no hits
bt054 blastall 2.0.11 sp:P39483  sts       -p tblastn
bt055 blastall 2.0.11 sp:P19888  sts       -p tblastn     no hits
bt056 blastall 2.0.11 gi:1348853 sts       -p tblastx
bt057 blastall 2.0.11 sp:Q08386  swissprot -p blastp      >1 HSP in alignment
bt058 blastall 2.0.11 sp:Q05362  swissprot -p blastp -v 0 -m 4 no descriptions, master-slave
bt059 blastpgp 2.0.11 sp:P07175  swissprot -m 6           multiple HSP (41943)
bt060 blastpgp 2.0.12                                     multiple db, Mike Poidinger
bt061 NCBI WWW 2.0.14 sp:P06563  swissprot blastp
bt062 blastall 2.0.14 gi:1348916 sts       -p blastn
bt063 blastpgp 2.0.14 ad-hoc                              from Mike Poidinger, no model sequences found
bt064 NCBIWWW  2.1.12 ad-hoc                              from Brad Chapman
bt065 blastn   2.1.2  ad-hoc                              no hits found
bt066 NCBIWWW  2.1.3  ad-hoc                              6/3/01 new format
bt067 blastall 2.2.1  ad-hoc              -p blastx       No 'length of query'
...
bt072 NCBIWWW  2.2.4  ad-hoc
bt073 NCBIWWW  2.2.4  ad-hoc   Confusing sequences line in the database
bt074 NCBIWWW  2.2.4  ad-hoc              -p blastp       nov 2002 new format
bt075 blastpgp 2.2.8  ad-hoc from Jer-Yee John Chuang, no Searching... lines 
bt076 blastall 2.2.15 ad-hoc              -p blastx       no hits, single query
bt077 blastall 2.2.20 ad-hoc     nr       -p blastx       single query
bt078 blastall 2.2.20 opuntia    nr       -p blastx       multiple queries, most with no hits
bt079 blastall 2.2.21 ad-hoc     nr       -p blastp       multiple queries
bt080 blastall 2.2.22 ad-hoc     nr       -p blastx       multiple queries
bt081 blastx   2.2.22+           nr                       the new C++ BLAST toolkit, multiple queries

* bt008 is not included because the file is too big (>15Mb).

BLAST OUTPUT FILES (XML)
------------------------

xbt001.xml - BLASTP 2.2.12
xbt002.xml - BLASTN 2.2.12
xbt003.xml - BLASTX 2.2.12
xbt004.xml - TBLASTN 2.2.12
xbt005.xml - TBLASTX 2.2.12
xbt006.xml - BLASTP 2.2.18+
xbt007.xml - BLASTP 2.2.18+
xbt008.xml - blastp 2.2.18
xbt009.xml - blastx 2.2.22+ (the new C++ tool blastx, not blastall - compare with bt081.txt file)
xbt010.xml - blastp 2.2.22+ (the new C++ tool blastp, not blastall)
xbt011.xml - RPSBLAST 2.2.18 (the old C tool blastpgp)
