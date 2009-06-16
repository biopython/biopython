# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from Bio import ParserSupport
from Bio.Blast import NCBIStandalone


all_tests = [
    'bt001', 'bt002', 'bt003', 'bt004', 'bt005',
    'bt006', 'bt007', 'bt009', 'bt010', 'bt011',
    'bt012', 'bt013', 'bt014', 'bt015', 'bt016',
    'bt017', 'bt018', 'bt039', 'bt040', 'bt041',
    'bt042', 'bt043', 'bt044', 'bt045', 'bt046',
    'bt047', 'bt048', 'bt049', 'bt050', 'bt051',
    'bt052', 'bt053', 'bt054', 'bt055', 'bt056',
    'bt057', 'bt058', 'bt059', 'bt060', 'bt062',
    'bt063', 'bt067', 'bt071'
    ]
# In order to keep the output file sizes reasonable, only generate
# a bunch of output for a few of the tests.
detailed_tests = [
    'bt001',   # 2.0.10 blastp
    'bt003',   # 2.0.10 blastp master-slave
    'bt006',   # 2.0.10 blastpgp
    'bt010',   # 2.0.10 blastn
    'bt014',   # 2.0.10 blastx
    'bt016',   # 2.0.10 tblastn
    'bt018',   # 2.0.10 tblastx
    'bt042',   # 2.0.11 blastp
    ]


### _Scanner

print "Running tests on _Scanner"

scanner = NCBIStandalone._Scanner()
for test in all_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    scanner.feed(open(datafile), ParserSupport.AbstractConsumer())

for test in detailed_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    scanner.feed(open(datafile), ParserSupport.TaggingConsumer())

### BlastParser

print "Running tests on BlastParser"

parser = NCBIStandalone.BlastParser()
pb_parser = NCBIStandalone.PSIBlastParser()
for test in all_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    try:
        # First, try parsing it with the normal parser.
        rec = parser.parse(open(datafile))
    except ValueError, x:
        # If it complains that the input is psiblast data, then
        # parse it with the psiblast parser.
        if 'PSI-BLAST data' in str(x):
            rec = pb_parser.parse(open(datafile))
        else:
            raise

### Blast Record
# - incredibly incomplete. Just testing what I'm adding -- Brad.

print "Running tests on Blast Records"

datafile = os.path.join("Blast", all_tests[4]) # bt005
rec = parser.parse(open(datafile))

print "\tTesting conversion of multiple alignments"
from Bio.Alphabet import IUPAC
generic_align = rec.multiple_alignment.to_generic(IUPAC.protein)
test_seq = generic_align.get_seq_by_num(0)
assert test_seq.alphabet == IUPAC.protein
assert test_seq.data[:60] == rec.multiple_alignment.alignment[0][2]


### Now some real tests

reference = 'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, \nJinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), \n"Gapped BLAST and PSI-BLAST: a new generation of protein database search\nprograms",  Nucleic Acids Res. 25:3389-3402.'

handle = open('Blast/bt001')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|120291|sp|P21297|FLBT_CAUCR FLBT PROTEIN"
assert record.query_letters==141
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==3
assert record.descriptions[0].title=="gi|120291|sp|P21297|FLBT_CAUCR FLBT PROTEIN"
assert record.descriptions[0].score==284
assert record.descriptions[0].e==7.0000000000000003e-77
assert record.descriptions[1].title=="gi|3024946|sp|Q58368|Y958_METJA HYPOTHETICAL PROTEIN MJ0958 PRE..."
assert record.descriptions[1].score==29
assert record.descriptions[1].e==3.400000
assert record.descriptions[2].title=="gi|3024745|sp|O26320|THSA_METTH PROBABLE THERMOSOME SUBUNIT A"
assert record.descriptions[2].score==29
assert record.descriptions[2].e==4.500000
assert len(record.alignments)==3
assert record.alignments[0].title==">gi|120291|sp|P21297|FLBT_CAUCR FLBT PROTEIN"
assert record.alignments[0].length==141
assert record.alignments[1].title==">gi|3024946|sp|Q58368|Y958_METJA HYPOTHETICAL PROTEIN MJ0958 PRECURSOR"
assert record.alignments[1].length==426
assert record.alignments[2].title==">gi|3024745|sp|O26320|THSA_METTH PROBABLE THERMOSOME SUBUNIT A"
assert record.alignments[2].length==542
assert record.alignments[0].hsps[0].score==718
assert record.alignments[0].hsps[0].bits==284
assert record.alignments[0].hsps[0].expect==7e-77
assert record.alignments[0].hsps[0].identities==(141, 141)
assert record.alignments[0].hsps[0].positives==(141, 141)
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==64
assert record.alignments[1].hsps[0].bits==29.3
assert record.alignments[1].hsps[0].expect==3.4
assert record.alignments[1].hsps[0].identities==(15, 47)
assert record.alignments[1].hsps[0].positives==(23, 47)
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==63
assert record.alignments[2].hsps[0].bits==29.0
assert record.alignments[2].hsps[0].expect==4.5
assert record.alignments[2].hsps[0].identities==(31, 107)
assert record.alignments[2].hsps[0].positives==(46, 107)
assert record.alignments[2].hsps[0].gaps==(9, 107)
assert len(record.alignments[2].hsps)==1
assert record.alignments[0].hsps[0].query=="MPLKLSLKPGEKFVLNGAVVQNGDRRGVLVLQNKASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFYEEFATRLNEFMGVVRNPVVLQDCIAISKHVMAREYYKALMLSRKLIEYEDERLGHVSSGVSAGGDAG"
assert record.alignments[0].hsps[0].match=="MPLKLSLKPGEKFVLNGAVVQNGDRRGVLVLQNKASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFYEEFATRLNEFMGVVRNPVVLQDCIAISKHVMAREYYKALMLSRKLIEYEDERLGHVSSGVSAGGDAG"
assert record.alignments[0].hsps[0].sbjct=="MPLKLSLKPGEKFVLNGAVVQNGDRRGVLVLQNKASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFYEEFATRLNEFMGVVRNPVVLQDCIAISKHVMAREYYKALMLSRKLIEYEDERLGHVSSGVSAGGDAG"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==141
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==141
assert record.alignments[1].hsps[0].query=="VLVLQNKASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFY"
assert record.alignments[1].hsps[0].match=="+LVL N  ++   K     D  TT   +IY P+ +    +  A+KFY"
assert record.alignments[1].hsps[0].sbjct=="ILVLINNTNITELKKFEDDDYYTTFQHYIYQPIFIFTTYDSKAKKFY"
assert record.alignments[1].hsps[0].query_start==28
assert record.alignments[1].hsps[0].query_end==74
assert record.alignments[1].hsps[0].sbjct_start==169
assert record.alignments[1].hsps[0].sbjct_end==215
assert record.alignments[2].hsps[0].query=="KASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFYEEFATRLNEFMGVVRNPVVLQDCIAISKHVMAREYYKALMLSRKLIEYEDERLGHVSSGVSAGGDA"
assert record.alignments[2].hsps[0].match=="+A V+ EK I   D+V T       P  +  +     +   EE    L + +GVV +   L+D     + V+A      + ++RKL EY D   G     VSA  DA"
assert record.alignments[2].hsps[0].sbjct=="EAGVIYEKKIF--DEVLTFIEECRDPKAISIILRGSTKHVAEEMERALEDAIGVVAS--TLED-----REVVAGGGAPEVEIARKLREYADTISGREQLAVSAFADA"
assert record.alignments[2].hsps[0].query_start==34
assert record.alignments[2].hsps[0].query_end==140
assert record.alignments[2].hsps[0].sbjct_start==339
assert record.alignments[2].hsps[0].sbjct_end==436
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.321, 0.138, 0.390]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==7927717
assert record.num_sequences==82258
assert record.num_extends==284596
assert record.num_good_extends==567
assert record.num_seqs_better_e==3
assert record.hsps_no_gap==2
assert record.hsps_prelim_gapped==1
assert record.hsps_gapped==3
assert record.query_length==141
assert record.database_length==29652561
assert record.effective_hsp_length==50
assert record.effective_query_length==91
assert record.effective_database_length==25539661
assert record.effective_search_space==2324109151
assert record.effective_search_space_used==2324109151
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.4)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.8)
assert record.blast_cutoff==(61,28.2)

handle = open('Blast/bt002')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|400206|sp|Q02112|LYTA_BACSU MEMBRANE-BOUND PROTEIN LYTA\nPRECURSOR"
assert record.query_letters==102
assert record.database=="data/pdbaa"
assert record.database_sequences==6999
assert record.database_letters==1461753
assert len(record.descriptions)==0
assert len(record.alignments)==0
assert record.database_name==['data/pdbaa']
assert record.num_letters_in_database==[1461753]
assert record.num_sequences_in_database==[6999]
assert record.posted_date==[('Oct 11, 1999 11:30 AM',)]
assert record.ka_params==[0.310, 0.132, 0.353]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==289135
assert record.num_sequences==6999
assert record.num_extends==10742
assert record.num_good_extends==10
assert record.num_seqs_better_e==0
assert record.hsps_no_gap==0
assert record.hsps_prelim_gapped==0
assert record.hsps_gapped==0
assert record.query_length==102
assert record.database_length==1461753
assert record.effective_hsp_length==45
assert record.effective_query_length==57
assert record.effective_database_length==1146798
assert record.effective_search_space==65367486
assert record.effective_search_space_used==65367486
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.1)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(42,21.7)
assert record.blast_cutoff==(47,22.7)

handle = open('Blast/bt003')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1718062|sp|P16153|UTXA_CLODI UTXA PROTEIN"
assert record.query_letters==166
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==0
assert len(record.alignments)==6
assert record.alignments[0].title==">gi|1718062|sp|P16153|UTXA_CLODI UTXA PROTEIN"
assert record.alignments[0].length==166
assert record.alignments[1].title==">gi|140528|sp|P24811|YQXH_BACSU HYPOTHETICAL 15.7 KD PROTEIN IN SPOIIIC-CWLA INTERGENIC REGION (ORF2)"
assert record.alignments[1].length==140
assert record.alignments[2].title==">gi|141088|sp|P26835|YNGD_CLOPE HYPOTHETICAL 14.9 KD PROTEIN IN NAGH 3'REGION (ORFD)"
assert record.alignments[2].length==132
assert record.alignments[3].title==">gi|6014830|sp|O78935|CYB_MARAM CYTOCHROME B"
assert record.alignments[3].length==379
assert record.alignments[4].title==">gi|1351589|sp|P47694|Y456_MYCGE HYPOTHETICAL PROTEIN MG456"
assert record.alignments[4].length==334
assert record.alignments[5].title==">gi|2496246|sp|Q57881|Y439_METJA HYPOTHETICAL ATP-BINDING PROTEIN MJ0439"
assert record.alignments[5].length==361
assert record.alignments[0].hsps[0].score==843
assert record.alignments[0].hsps[0].bits==332
assert record.alignments[0].hsps[0].expect==2e-91
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==90
assert record.alignments[1].hsps[0].bits==39.5
assert record.alignments[1].hsps[0].expect==0.004
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==88
assert record.alignments[2].hsps[0].bits==38.7
assert record.alignments[2].hsps[0].expect==0.007
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==64
assert record.alignments[3].hsps[0].bits==29.3
assert record.alignments[3].hsps[0].expect==4.6
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==63
assert record.alignments[4].hsps[0].bits==29.0
assert record.alignments[4].hsps[0].expect==6.0
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==62
assert record.alignments[5].hsps[0].bits==28.6
assert record.alignments[5].hsps[0].expect==7.8
assert len(record.alignments[5].hsps)==1
assert record.alignments[0].hsps[0].identities==(166, 166)
assert record.alignments[0].hsps[0].positives==(166, 166)
assert record.alignments[1].hsps[0].identities==(27, 130)
assert record.alignments[1].hsps[0].positives==(55, 130)
assert record.alignments[1].hsps[0].gaps==(19, 130)
assert record.alignments[2].hsps[0].identities==(24, 110)
assert record.alignments[2].hsps[0].positives==(52, 110)
assert record.alignments[2].hsps[0].gaps==(18, 110)
assert record.alignments[3].hsps[0].identities==(19, 57)
assert record.alignments[3].hsps[0].positives==(33, 57)
assert record.alignments[3].hsps[0].gaps==(2, 57)
assert record.alignments[4].hsps[0].identities==(16, 44)
assert record.alignments[4].hsps[0].positives==(24, 44)
assert record.alignments[4].hsps[0].gaps==(2, 44)
assert record.alignments[5].hsps[0].identities==(19, 56)
assert record.alignments[5].hsps[0].positives==(30, 56)
assert record.alignments[5].hsps[0].gaps==(12, 56)
assert record.alignments[0].hsps[0].query=="MHSSSPFYISNGNKIFFYINLGGVMNMTISFLSEHIFIKLVILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCLCGLPVPKRLKEKIAILLDAMTDEMNAKDEK"
assert record.alignments[0].hsps[0].match=="MHSSSPFYISNGNKIFFYINLGGVMNMTISFLSEHIFIKLVILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCLCGLPVPKRLKEKIAILLDAMTDEMNAKDEK"
assert record.alignments[0].hsps[0].sbjct=="MHSSSPFYISNGNKIFFYINLGGVMNMTISFLSEHIFIKLVILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCLCGLPVPKRLKEKIAILLDAMTDEMNAKDEK"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==166
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==166
assert record.alignments[1].hsps[0].query=="FIKLVILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILF-ILYESVSILKNMCLCGLPVPKRLKEKIAILLDAMTDEMNAKDE"
assert record.alignments[1].hsps[0].match=="++ L+++    D L G + A K +K  S     G +RK+     +   +V+D +   N                  G+  F ++LF I  E +SI +N+   G+ +P  + +++  + +      N  D+"
assert record.alignments[1].hsps[0].sbjct=="YLDLLLVLSIIDVLTGVIKAWKFKKLRSRSAWFGYVRKLLNFFAVILANVIDTVLNLN------------------GVLTFGTVLFYIANEGLSITENLAQIGVKIPSSITDRLQTIENEKEQSKNNADK"
assert record.alignments[1].hsps[0].query_start==37
assert record.alignments[1].hsps[0].query_end==165
assert record.alignments[1].hsps[0].sbjct_start==26
assert record.alignments[1].hsps[0].sbjct_end==137
assert record.alignments[2].hsps[0].query=="VILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVD-ILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCLCGLPVPKRLKEKI"
assert record.alignments[2].hsps[0].match=="+++ I  D L G +   KS++  S+ G+ G  +K  ++  +    ++D +L    ++F                     +  +I+ E +SIL+N    G+P+P++LK+ +"
assert record.alignments[2].hsps[0].sbjct=="LLVFIFLDYLTGVIKGCKSKELCSNIGLRGITKKGLILVVLLVAVMLDRLLDNGTWMFRT-----------------LIAYFYIMNEGISILENCAALGVPIPEKLKQAL"
assert record.alignments[2].hsps[0].query_start==41
assert record.alignments[2].hsps[0].query_end==149
assert record.alignments[2].hsps[0].sbjct_start==33
assert record.alignments[2].hsps[0].sbjct_end==125
assert record.alignments[3].hsps[0].query=="CIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCL"
assert record.alignments[3].hsps[0].match=="C+F+L V D+LT   ++   P +   F+ +  L    +F+IL IL  ++SI++N  L"
assert record.alignments[3].hsps[0].sbjct=="CLFWLLVADLLT-LTWIGGQPVEH-PFITIGQLASILYFAILLILMPAISIIENNLL"
assert record.alignments[3].hsps[0].query_start==80
assert record.alignments[3].hsps[0].query_end==136
assert record.alignments[3].hsps[0].sbjct_start==323
assert record.alignments[3].hsps[0].sbjct_end==377
assert record.alignments[4].hsps[0].query=="LTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKN"
assert record.alignments[4].hsps[0].match=="LTKFN  F+ P     FLR+  +G+   FS++ I +   S  +N"
assert record.alignments[4].hsps[0].sbjct=="LTKFNKFFLTPNKLNAFLRV--IGLCGLFSVIAISFGIYSYTRN"
assert record.alignments[4].hsps[0].query_start==90
assert record.alignments[4].hsps[0].query_end==133
assert record.alignments[4].hsps[0].sbjct_start==4
assert record.alignments[4].hsps[0].sbjct_end==45
assert record.alignments[5].hsps[0].query=="FLRLKHLGIS---EFFSILFILYES----VSILKNMC-----LCGLPVPKRLKEKI"
assert record.alignments[5].hsps[0].match=="++ L+ + IS   +F  +LF  YE     V I+K++      LCG+P PK   E+I"
assert record.alignments[5].hsps[0].sbjct=="YINLRGIFISKYKDFIEVLFEEYEEDRKPVEIIKSLIKDVPSLCGIPTPKNTLEEI"
assert record.alignments[5].hsps[0].query_start==106
assert record.alignments[5].hsps[0].query_end==149
assert record.alignments[5].hsps[0].sbjct_start==68
assert record.alignments[5].hsps[0].sbjct_end==123
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.331, 0.146, 0.428]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==8801581
assert record.num_sequences==82258
assert record.num_extends==320828
assert record.num_good_extends==892
assert record.num_seqs_better_e==6
assert record.hsps_no_gap==3
assert record.hsps_prelim_gapped==3
assert record.hsps_gapped==6
assert record.query_length==166
assert record.database_length==29652561
assert record.effective_hsp_length==46
assert record.effective_query_length==120
assert record.effective_database_length==25868693
assert record.effective_search_space==3104243160
assert record.effective_search_space_used==3104243160
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(15,7.2)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(40,21.9)
assert record.blast_cutoff==(62,28.6)

handle = open('Blast/bt004')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1718062|sp|P16153|UTXA_CLODI UTXA PROTEIN"
assert record.query_letters==166
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==6
assert record.descriptions[0].title=="gi|1718062|sp|P16153|UTXA_CLODI UTXA PROTEIN"
assert record.descriptions[0].score==332
assert record.descriptions[0].e==2.e-91
assert record.descriptions[1].title=="gi|140528|sp|P24811|YQXH_BACSU HYPOTHETICAL 15.7 KD PROTEIN IN ..."
assert record.descriptions[1].score==39
assert record.descriptions[1].e==0.004000
assert record.descriptions[2].title=="gi|141088|sp|P26835|YNGD_CLOPE HYPOTHETICAL 14.9 KD PROTEIN IN ..."
assert record.descriptions[2].score==39
assert record.descriptions[2].e==0.007000
assert record.descriptions[3].title=="gi|6014830|sp|O78935|CYB_MARAM CYTOCHROME B"
assert record.descriptions[3].score==29
assert record.descriptions[3].e==4.600000
assert record.descriptions[4].title=="gi|1351589|sp|P47694|Y456_MYCGE HYPOTHETICAL PROTEIN MG456"
assert record.descriptions[4].score==29
assert record.descriptions[4].e==6.000000
assert record.descriptions[5].title=="gi|2496246|sp|Q57881|Y439_METJA HYPOTHETICAL ATP-BINDING PROTEI..."
assert record.descriptions[5].score==29
assert record.descriptions[5].e==7.800000
assert len(record.alignments)==0
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.331, 0.146, 0.428]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==8801581
assert record.num_sequences==82258
assert record.num_extends==320828
assert record.num_good_extends==892
assert record.num_seqs_better_e==6
assert record.hsps_no_gap==3
assert record.hsps_prelim_gapped==3
assert record.hsps_gapped==6
assert record.query_length==166
assert record.database_length==29652561
assert record.effective_hsp_length==46
assert record.effective_query_length==120
assert record.effective_database_length==25868693
assert record.effective_search_space==3104243160
assert record.effective_search_space_used==3104243160
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(15,7.2)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(40,21.9)
assert record.blast_cutoff==(62,28.6)

handle = open('Blast/bt005')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|132349|sp|P15394|REPA_AGRTU REPLICATING PROTEIN"
assert record.query_letters==250
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==15
assert record.descriptions[0].title=="gi|132349|sp|P15394|REPA_AGRTU REPLICATING PROTEIN"
assert record.descriptions[0].score==514
assert record.descriptions[0].e==1.e-146
assert record.descriptions[1].title=="gi|123932|sp|P19922|HYIN_BRAJA INDOLEACETAMIDE HYDROLASE (IAH) ..."
assert record.descriptions[1].score==34
assert record.descriptions[1].e==0.290000
assert record.descriptions[2].title=="gi|137670|sp|P06422|VE2_HPV08 REGULATORY PROTEIN E2"
assert record.descriptions[2].score==32
assert record.descriptions[2].e==0.860000
assert record.descriptions[3].title=="gi|5921693|sp|Q05152|CCAB_RABIT VOLTAGE-DEPENDENT N-TYPE CALCIU..."
assert record.descriptions[3].score==32
assert record.descriptions[3].e==1.500000
assert record.descriptions[4].title=="gi|121952|sp|P02256|H1_PARAN HISTONE H1, GONADAL"
assert record.descriptions[4].score==31
assert record.descriptions[4].e==2.500000
assert record.descriptions[5].title=="gi|124141|sp|P08392|ICP4_HSV11 TRANS-ACTING TRANSCRIPTIONAL PRO..."
assert record.descriptions[5].score==31
assert record.descriptions[5].e==3.300000
assert record.descriptions[6].title=="gi|3915729|sp|P51592|HYDP_DROME HYPERPLASTIC DISCS PROTEIN (HYD..."
assert record.descriptions[6].score==31
assert record.descriptions[6].e==3.300000
assert record.descriptions[7].title=="gi|462182|sp|P33438|GLT_DROME GLUTACTIN PRECURSOR"
assert record.descriptions[7].score==31
assert record.descriptions[7].e==3.300000
assert record.descriptions[8].title=="gi|731294|sp|P39713|YAG1_YEAST HYPOTHETICAL ZINC-TYPE ALCOHOL D..."
assert record.descriptions[8].score==30
assert record.descriptions[8].e==4.300000
assert record.descriptions[9].title=="gi|1708851|sp|P55268|LMB2_HUMAN LAMININ BETA-2 CHAIN PRECURSOR ..."
assert record.descriptions[9].score==30
assert record.descriptions[9].e==4.300000
assert record.descriptions[10].title=="gi|2495137|sp|Q24704|H1_DROVI HISTONE H1"
assert record.descriptions[10].score==29
assert record.descriptions[10].e==7.500000
assert record.descriptions[11].title=="gi|1172635|sp|P46466|PRS4_ORYSA 26S PROTEASE REGULATORY SUBUNIT..."
assert record.descriptions[11].score==29
assert record.descriptions[11].e==9.800000
assert record.descriptions[12].title=="gi|6093970|sp|Q61085|RHOP_MOUSE GTP-RHO BINDING PROTEIN 1 (RHOP..."
assert record.descriptions[12].score==29
assert record.descriptions[12].e==9.800000
assert record.descriptions[13].title=="gi|547963|sp|Q01989|MYS9_DROME MYOSIN HEAVY CHAIN 95F (95F MHC)"
assert record.descriptions[13].score==29
assert record.descriptions[13].e==9.800000
assert record.descriptions[14].title=="gi|6226905|sp|Q59567|TOP1_MYCTU DNA TOPOISOMERASE I (OMEGA-PROT..."
assert record.descriptions[14].score==29
assert record.descriptions[14].e==9.800000
assert len(record.alignments)==0
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.317, 0.133, 0.395]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==14679054
assert record.num_sequences==82258
assert record.num_extends==599117
assert record.num_good_extends==1508
assert record.num_seqs_better_e==15
assert record.hsps_no_gap==4
assert record.hsps_prelim_gapped==11
assert record.hsps_gapped==15
assert record.query_length==250
assert record.database_length==29652561
assert record.effective_hsp_length==51
assert record.effective_query_length==199
assert record.effective_database_length==25457403
assert record.effective_search_space==5066023197
assert record.effective_search_space_used==5066023197
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(63,29.0)

handle = open('Blast/bt006')
record = pb_parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1174726|sp|P12921|TMRB_BACSU TUNICAMYCIN RESISTANCE PROTEIN"
assert record.query_letters==197
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.rounds)==1
assert len(record.rounds[0].new_seqs)==4
assert record.rounds[0].new_seqs[0].title=="gi|1174726|sp|P12921|TMRB_BACSU TUNICAMYCIN RESISTANCE PROTEIN"
assert record.rounds[0].new_seqs[0].score==402
assert record.rounds[0].new_seqs[0].e==9.9999999999999995e-113
assert record.rounds[0].new_seqs[1].title=="gi|1352934|sp|P47169|YJ9F_YEAST HYPOTHETICAL 161.2 KD PROTEIN I..."
assert record.rounds[0].new_seqs[1].score==30
assert record.rounds[0].new_seqs[1].e==3.300000
assert record.rounds[0].new_seqs[2].title=="gi|3915741|sp|P04407|KITH_HSV23 THYMIDINE KINASE"
assert record.rounds[0].new_seqs[2].score==29
assert record.rounds[0].new_seqs[2].e==7.400000
assert record.rounds[0].new_seqs[3].title=="gi|3334489|sp|P15398|RPA1_SCHPO DNA-DIRECTED RNA POLYMERASE I 1..."
assert record.rounds[0].new_seqs[3].score==29
assert record.rounds[0].new_seqs[3].e==7.400000
assert len(record.rounds[0].alignments)==4
assert record.rounds[0].alignments[0].title==">gi|1174726|sp|P12921|TMRB_BACSU TUNICAMYCIN RESISTANCE PROTEIN"
assert record.rounds[0].alignments[0].length==197
assert record.rounds[0].alignments[1].title==">gi|1352934|sp|P47169|YJ9F_YEAST HYPOTHETICAL 161.2 KD PROTEIN IN NMD5-HOM6 INTERGENIC REGION"
assert record.rounds[0].alignments[1].length==1442
assert record.rounds[0].alignments[2].title==">gi|3915741|sp|P04407|KITH_HSV23 THYMIDINE KINASE"
assert record.rounds[0].alignments[2].length==376
assert record.rounds[0].alignments[3].title==">gi|3334489|sp|P15398|RPA1_SCHPO DNA-DIRECTED RNA POLYMERASE I 190 KD POLYPEPTIDE"
assert record.rounds[0].alignments[3].length==1689
assert record.rounds[0].alignments[0].hsps[0].score==1021
assert record.rounds[0].alignments[0].hsps[0].bits==402
assert record.rounds[0].alignments[0].hsps[0].expect==1e-112
assert len(record.rounds[0].alignments[0].hsps)==1
assert record.rounds[0].alignments[1].hsps[0].score==66
assert record.rounds[0].alignments[1].hsps[0].bits==30.1
assert record.rounds[0].alignments[1].hsps[0].expect==3.3
assert len(record.rounds[0].alignments[1].hsps)==1
assert record.rounds[0].alignments[2].hsps[0].score==63
assert record.rounds[0].alignments[2].hsps[0].bits==29.0
assert record.rounds[0].alignments[2].hsps[0].expect==7.4
assert len(record.rounds[0].alignments[2].hsps)==1
assert record.rounds[0].alignments[3].hsps[0].score==63
assert record.rounds[0].alignments[3].hsps[0].bits==29.0
assert record.rounds[0].alignments[3].hsps[0].expect==7.4
assert len(record.rounds[0].alignments[3].hsps)==1
assert record.rounds[0].alignments[0].hsps[0].identities==(197, 197)
assert record.rounds[0].alignments[0].hsps[0].positives==(197, 197)
assert record.rounds[0].alignments[1].hsps[0].identities==(23, 70)
assert record.rounds[0].alignments[1].hsps[0].positives==(35, 70)
assert record.rounds[0].alignments[1].hsps[0].gaps==(10, 70)
assert record.rounds[0].alignments[2].hsps[0].identities==(15, 37)
assert record.rounds[0].alignments[2].hsps[0].positives==(22, 37)
assert record.rounds[0].alignments[2].hsps[0].gaps==(2, 37)
assert record.rounds[0].alignments[3].hsps[0].identities==(12, 38)
assert record.rounds[0].alignments[3].hsps[0].positives==(20, 38)
assert record.rounds[0].alignments[0].hsps[0].query=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[0].alignments[0].hsps[0].match=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[0].alignments[0].hsps[0].sbjct=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[0].alignments[0].hsps[0].query_start==1
assert record.rounds[0].alignments[0].hsps[0].query_end==197
assert record.rounds[0].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[0].hsps[0].sbjct_end==197
assert record.rounds[0].alignments[1].hsps[0].query=="TLMASKETLLKR--------LRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARA"
assert record.rounds[0].alignments[1].hsps[0].match=="TL+  K+  L R          TR + K S AA   D+ +EGLS P    +  +D  +  ++A+ +AARA"
assert record.rounds[0].alignments[1].hsps[0].sbjct=="TLLTRKDPSLSRNLKQSAGDALTRKQEKRSKAA--FDQLLEGLSGPPLHVYYASDGGNAANLAKRLAARA"
assert record.rounds[0].alignments[1].hsps[0].query_start==107
assert record.rounds[0].alignments[1].hsps[0].query_end==168
assert record.rounds[0].alignments[1].hsps[0].sbjct_start==637
assert record.rounds[0].alignments[1].hsps[0].sbjct_end==704
assert record.rounds[0].alignments[2].hsps[0].query=="IWINGAFGSGKTQTAFELHRRLNP--SYVYDPEKMGF"
assert record.rounds[0].alignments[2].hsps[0].match=="++I+G  G GKT T+ +L   L P  + VY PE M +"
assert record.rounds[0].alignments[2].hsps[0].sbjct=="VYIDGPHGVGKTTTSAQLMEALGPRDNIVYVPEPMTY"
assert record.rounds[0].alignments[2].hsps[0].query_start==3
assert record.rounds[0].alignments[2].hsps[0].query_end==37
assert record.rounds[0].alignments[2].hsps[0].sbjct_start==52
assert record.rounds[0].alignments[2].hsps[0].sbjct_end==88
assert record.rounds[0].alignments[3].hsps[0].query=="GILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASK"
assert record.rounds[0].alignments[3].hsps[0].match=="G +++P+   HP +F+++   LR      HHF L   K"
assert record.rounds[0].alignments[3].hsps[0].sbjct=="GHIVLPIPAYHPLFFSQMYNLLRSTCLYCHHFKLSKVK"
assert record.rounds[0].alignments[3].hsps[0].query_start==75
assert record.rounds[0].alignments[3].hsps[0].query_end==112
assert record.rounds[0].alignments[3].hsps[0].sbjct_start==78
assert record.rounds[0].alignments[3].hsps[0].sbjct_end==115
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.323, 0.138, 0.411]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==11332752
assert record.num_sequences==82258
assert record.num_extends==435116
assert record.num_good_extends==987
assert record.num_seqs_better_e==4
assert record.hsps_no_gap==2
assert record.hsps_prelim_gapped==2
assert record.hsps_gapped==4
assert record.query_length==197
assert record.database_length==29652561
assert record.effective_hsp_length==48
assert record.effective_query_length==149
assert record.effective_database_length==25704177
assert record.effective_search_space==3829922373
assert record.effective_search_space_used==3829922373
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.5)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,22.0)
assert record.blast_cutoff==(62,28.6)

handle = open('Blast/bt007')
record = pb_parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|126343|sp|P17216|LIVK_SALTY LEUCINE-SPECIFIC BINDING PROTEIN\nPRECURSOR (LS-BP) (L-BP)"
assert record.query_letters==369
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.rounds)==3
assert len(record.rounds[0].new_seqs)==14
assert record.rounds[0].new_seqs[0].title=="gi|126343|sp|P17216|LIVK_SALTY LEUCINE-SPECIFIC BINDING PROTEIN..."
assert record.rounds[0].new_seqs[0].score==743
assert record.rounds[0].new_seqs[0].e==0.000000
assert record.rounds[0].new_seqs[1].title=="gi|126349|sp|P04816|LIVK_ECOLI LEUCINE-SPECIFIC BINDING PROTEIN..."
assert record.rounds[0].new_seqs[1].score==679
assert record.rounds[0].new_seqs[1].e==0.000000
assert record.rounds[0].new_seqs[2].title=="gi|126348|sp|P02917|LIVJ_ECOLI LEU/ILE/VAL-BINDING PROTEIN PREC..."
assert record.rounds[0].new_seqs[2].score==581
assert record.rounds[0].new_seqs[2].e==1e-166
assert record.rounds[0].new_seqs[3].title=="gi|126342|sp|P17215|LIVJ_SALTY LEU/ILE/VAL/THR-BINDING PROTEIN ..."
assert record.rounds[0].new_seqs[3].score==577
assert record.rounds[0].new_seqs[3].e==9.9999999999999996e-165
assert record.rounds[0].new_seqs[4].title=="gi|126347|sp|P25399|LIVJ_CITFR LEU/ILE/VAL-BINDING PROTEIN PREC..."
assert record.rounds[0].new_seqs[4].score==570
assert record.rounds[0].new_seqs[4].e==9.9999999999999995e-163
assert record.rounds[0].new_seqs[5].title=="gi|115120|sp|P21175|BRAC_PSEAE LEUCINE-, ISOLEUCINE-, VALINE-, ..."
assert record.rounds[0].new_seqs[5].score==413
assert record.rounds[0].new_seqs[5].e==1.0000000000000001e-115
assert record.rounds[0].new_seqs[6].title=="gi|113709|sp|P27017|AMIC_PSEAE ALIPHATIC AMIDASE EXPRESSION-REG..."
assert record.rounds[0].new_seqs[6].score==40
assert record.rounds[0].new_seqs[6].e==0.006000
assert record.rounds[0].new_seqs[7].title=="gi|127751|sp|P02567|MYSD_CAEEL MYOSIN HEAVY CHAIN D (MHC D)"
assert record.rounds[0].new_seqs[7].score==32
assert record.rounds[0].new_seqs[7].e==1.400000
assert record.rounds[0].new_seqs[8].title=="gi|131068|sp|P23596|PRTD_ERWCH PROTEASES SECRETION ATP-BINDING ..."
assert record.rounds[0].new_seqs[8].score==32
assert record.rounds[0].new_seqs[8].e==1.800000
assert record.rounds[0].new_seqs[9].title=="gi|2495081|sp|Q09630|MGR1_CAEEL PROBABLE METABOTROPIC GLUTAMATE..."
assert record.rounds[0].new_seqs[9].score==32
assert record.rounds[0].new_seqs[9].e==1.800000
assert record.rounds[0].new_seqs[10].title=="gi|2506848|sp|P80040|MDH_CHLAU MALATE DEHYDROGENASE"
assert record.rounds[0].new_seqs[10].score==32
assert record.rounds[0].new_seqs[10].e==1.800000
assert record.rounds[0].new_seqs[11].title=="gi|3915245|sp|Q38394|VG37_BPK3 TAIL FIBER PROTEIN GP37 (RECEPTO..."
assert record.rounds[0].new_seqs[11].score==31
assert record.rounds[0].new_seqs[11].e==4.000000
assert record.rounds[0].new_seqs[12].title=="gi|1351210|sp|P47209|TCPE_CAEEL T-COMPLEX PROTEIN 1, EPSILON SU..."
assert record.rounds[0].new_seqs[12].score==31
assert record.rounds[0].new_seqs[12].e==4.000000
assert record.rounds[0].new_seqs[13].title=="gi|1351310|sp|P43496|TRXB_PENCH THIOREDOXIN REDUCTASE"
assert record.rounds[0].new_seqs[13].score==31
assert record.rounds[0].new_seqs[13].e==5.200000
assert len(record.rounds[0].alignments)==14
assert record.rounds[0].alignments[0].title==">gi|126343|sp|P17216|LIVK_SALTY LEUCINE-SPECIFIC BINDING PROTEIN PRECURSOR (LS-BP) (L-BP)"
assert record.rounds[0].alignments[0].length==369
assert record.rounds[0].alignments[1].title==">gi|126349|sp|P04816|LIVK_ECOLI LEUCINE-SPECIFIC BINDING PROTEIN PRECURSOR (LS-BP) (L-BP)"
assert record.rounds[0].alignments[1].length==369
assert record.rounds[0].alignments[2].title==">gi|126348|sp|P02917|LIVJ_ECOLI LEU/ILE/VAL-BINDING PROTEIN PRECURSOR (LIV-BP)"
assert record.rounds[0].alignments[2].length==367
assert record.rounds[0].alignments[3].title==">gi|126342|sp|P17215|LIVJ_SALTY LEU/ILE/VAL/THR-BINDING PROTEIN PRECURSOR (LIVT-BP)"
assert record.rounds[0].alignments[3].length==365
assert record.rounds[0].alignments[4].title==">gi|126347|sp|P25399|LIVJ_CITFR LEU/ILE/VAL-BINDING PROTEIN PRECURSOR (LIV-BP)"
assert record.rounds[0].alignments[4].length==367
assert record.rounds[0].alignments[5].title==">gi|115120|sp|P21175|BRAC_PSEAE LEUCINE-, ISOLEUCINE-, VALINE-, THREONINE-, AND ALANINE-BINDING PROTEIN PRECURSOR (LIVAT-BP) (LEU/ILE/VAL/THR/ALA-BINDING PROTEIN)"
assert record.rounds[0].alignments[5].length==373
assert record.rounds[0].alignments[6].title==">gi|113709|sp|P27017|AMIC_PSEAE ALIPHATIC AMIDASE EXPRESSION-REGULATING PROTEIN"
assert record.rounds[0].alignments[6].length==385
assert record.rounds[0].alignments[7].title==">gi|127751|sp|P02567|MYSD_CAEEL MYOSIN HEAVY CHAIN D (MHC D)"
assert record.rounds[0].alignments[7].length==1938
assert record.rounds[0].alignments[8].title==">gi|131068|sp|P23596|PRTD_ERWCH PROTEASES SECRETION ATP-BINDING PROTEIN PRTD"
assert record.rounds[0].alignments[8].length==575
assert record.rounds[0].alignments[9].title==">gi|2495081|sp|Q09630|MGR1_CAEEL PROBABLE METABOTROPIC GLUTAMATE RECEPTOR MGL-1"
assert record.rounds[0].alignments[9].length==999
assert record.rounds[0].alignments[10].title==">gi|2506848|sp|P80040|MDH_CHLAU MALATE DEHYDROGENASE"
assert record.rounds[0].alignments[10].length==309
assert record.rounds[0].alignments[11].title==">gi|3915245|sp|Q38394|VG37_BPK3 TAIL FIBER PROTEIN GP37 (RECEPTOR RECOGNIZING PROTEIN)"
assert record.rounds[0].alignments[11].length==1243
assert record.rounds[0].alignments[12].title==">gi|1351210|sp|P47209|TCPE_CAEEL T-COMPLEX PROTEIN 1, EPSILON SUBUNIT (TCP-1-EPSILON) (CCT-EPSILON)"
assert record.rounds[0].alignments[12].length==542
assert record.rounds[0].alignments[13].title==">gi|1351310|sp|P43496|TRXB_PENCH THIOREDOXIN REDUCTASE"
assert record.rounds[0].alignments[13].length==334
assert len(record.rounds[1].new_seqs)==18
assert record.rounds[1].new_seqs[0].title=="gi|113709|sp|P27017|AMIC_PSEAE ALIPHATIC AMIDASE EXPRESSION-REG..."
assert record.rounds[1].new_seqs[0].score==49
assert record.rounds[1].new_seqs[0].e==2e-05
assert record.rounds[1].new_seqs[1].title=="gi|3024131|sp|P91685|MGR_DROME METABOTROPIC GLUTAMATE RECEPTOR ..."
assert record.rounds[1].new_seqs[1].score==41
assert record.rounds[1].new_seqs[1].e==0.003
assert record.rounds[1].new_seqs[2].title=="gi|121300|sp|P23385|MGR1_RAT METABOTROPIC GLUTAMATE RECEPTOR 1 ..."
assert record.rounds[1].new_seqs[2].score==39
assert record.rounds[1].new_seqs[2].e==0.012
assert record.rounds[1].new_seqs[3].title=="gi|2495074|sp|Q13255|MGR1_HUMAN METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[1].new_seqs[3].score==39
assert record.rounds[1].new_seqs[3].e==0.012
assert record.rounds[1].new_seqs[4].title=="gi|1170947|sp|P31424|MGR5_RAT METABOTROPIC GLUTAMATE RECEPTOR 5..."
assert record.rounds[1].new_seqs[4].score==36
assert record.rounds[1].new_seqs[4].e==0.11
assert record.rounds[1].new_seqs[5].title=="gi|1709020|sp|P41594|MGR5_HUMAN METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[1].new_seqs[5].score==36
assert record.rounds[1].new_seqs[5].e==0.14
assert record.rounds[1].new_seqs[6].title=="gi|2495081|sp|Q09630|MGR1_CAEEL PROBABLE METABOTROPIC GLUTAMATE..."
assert record.rounds[1].new_seqs[6].score==35
assert record.rounds[1].new_seqs[6].e==0.18
assert record.rounds[1].new_seqs[7].title=="gi|6014748|sp|Q10900|CTPI_MYCTU PROBABLE CATION-TRANSPORTING AT..."
assert record.rounds[1].new_seqs[7].score==33
assert record.rounds[1].new_seqs[7].e==0.71
assert record.rounds[1].new_seqs[8].title=="gi|2833574|sp|Q58178|Y768_METJA HYPOTHETICAL PROTEIN MJ0768"
assert record.rounds[1].new_seqs[8].score==32
assert record.rounds[1].new_seqs[8].e==1.6
assert record.rounds[1].new_seqs[9].title=="gi|1169299|sp|P45513|DHAT_CITFR 1,3-PROPANEDIOL DEHYDROGENASE (..."
assert record.rounds[1].new_seqs[9].score==31
assert record.rounds[1].new_seqs[9].e==3.6
assert record.rounds[1].new_seqs[10].title=="gi|2495080|sp|P70579|MGR8_RAT METABOTROPIC GLUTAMATE RECEPTOR 8..."
assert record.rounds[1].new_seqs[10].score==31
assert record.rounds[1].new_seqs[10].e==3.6
assert record.rounds[1].new_seqs[11].title=="gi|2495079|sp|O00222|MGR8_HUMAN METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[1].new_seqs[11].score==31
assert record.rounds[1].new_seqs[11].e==3.6
assert record.rounds[1].new_seqs[12].title=="gi|1346533|sp|P47743|MGR8_MOUSE METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[1].new_seqs[12].score==30
assert record.rounds[1].new_seqs[12].e==6.2
assert record.rounds[1].new_seqs[13].title=="gi|113381|sp|P06758|ADH2_ZYMMO ALCOHOL DEHYDROGENASE II (ADH II)"
assert record.rounds[1].new_seqs[13].score==30
assert record.rounds[1].new_seqs[13].e==8.1
assert record.rounds[1].new_seqs[14].title=="gi|127751|sp|P02567|MYSD_CAEEL MYOSIN HEAVY CHAIN D (MHC D)"
assert record.rounds[1].new_seqs[14].score==30
assert record.rounds[1].new_seqs[14].e==8.1
assert record.rounds[1].new_seqs[15].title=="gi|1346321|sp|P23897|HSER_RAT HEAT-STABLE ENTEROTOXIN RECEPTOR ..."
assert record.rounds[1].new_seqs[15].score==30
assert record.rounds[1].new_seqs[15].e==8.1
assert record.rounds[1].new_seqs[16].title=="gi|1175674|sp|P45116|YCIH_HAEIN HYPOTHETICAL PROTEIN HI1225"
assert record.rounds[1].new_seqs[16].score==30
assert record.rounds[1].new_seqs[16].e==8.1
assert record.rounds[1].new_seqs[17].title=="gi|3025270|sp|P77269|YPHF_ECOLI ABC TRANSPORTER PERIPLASMIC BIN..."
assert record.rounds[1].new_seqs[17].score==30
assert record.rounds[1].new_seqs[17].e==8.1
assert len(record.rounds[1].alignments)==24
assert record.rounds[1].alignments[0].title==">gi|126349|sp|P04816|LIVK_ECOLI LEUCINE-SPECIFIC BINDING PROTEIN PRECURSOR (LS-BP) (L-BP)"
assert record.rounds[1].alignments[0].length==369
assert record.rounds[1].alignments[1].title==">gi|126342|sp|P17215|LIVJ_SALTY LEU/ILE/VAL/THR-BINDING PROTEIN PRECURSOR (LIVT-BP)"
assert record.rounds[1].alignments[1].length==365
assert record.rounds[1].alignments[2].title==">gi|126347|sp|P25399|LIVJ_CITFR LEU/ILE/VAL-BINDING PROTEIN PRECURSOR (LIV-BP)"
assert record.rounds[1].alignments[2].length==367
assert record.rounds[1].alignments[3].title==">gi|126343|sp|P17216|LIVK_SALTY LEUCINE-SPECIFIC BINDING PROTEIN PRECURSOR (LS-BP) (L-BP)"
assert record.rounds[1].alignments[3].length==369
assert record.rounds[1].alignments[4].title==">gi|126348|sp|P02917|LIVJ_ECOLI LEU/ILE/VAL-BINDING PROTEIN PRECURSOR (LIV-BP)"
assert record.rounds[1].alignments[4].length==367
assert record.rounds[1].alignments[5].title==">gi|115120|sp|P21175|BRAC_PSEAE LEUCINE-, ISOLEUCINE-, VALINE-, THREONINE-, AND ALANINE-BINDING PROTEIN PRECURSOR (LIVAT-BP) (LEU/ILE/VAL/THR/ALA-BINDING PROTEIN)"
assert record.rounds[1].alignments[5].length==373
assert record.rounds[1].alignments[6].title==">gi|113709|sp|P27017|AMIC_PSEAE ALIPHATIC AMIDASE EXPRESSION-REGULATING PROTEIN"
assert record.rounds[1].alignments[6].length==385
assert record.rounds[1].alignments[7].title==">gi|3024131|sp|P91685|MGR_DROME METABOTROPIC GLUTAMATE RECEPTOR PRECURSOR"
assert record.rounds[1].alignments[7].length==976
assert record.rounds[1].alignments[8].title==">gi|121300|sp|P23385|MGR1_RAT METABOTROPIC GLUTAMATE RECEPTOR 1 PRECURSOR"
assert record.rounds[1].alignments[8].length==1199
assert record.rounds[1].alignments[9].title==">gi|2495074|sp|Q13255|MGR1_HUMAN METABOTROPIC GLUTAMATE RECEPTOR 1 PRECURSOR"
assert record.rounds[1].alignments[9].length==1194
assert record.rounds[1].alignments[10].title==">gi|1170947|sp|P31424|MGR5_RAT METABOTROPIC GLUTAMATE RECEPTOR 5 PRECURSOR"
assert record.rounds[1].alignments[10].length==1203
assert record.rounds[1].alignments[11].title==">gi|1709020|sp|P41594|MGR5_HUMAN METABOTROPIC GLUTAMATE RECEPTOR 5 PRECURSOR"
assert record.rounds[1].alignments[11].length==1212
assert record.rounds[1].alignments[12].title==">gi|2495081|sp|Q09630|MGR1_CAEEL PROBABLE METABOTROPIC GLUTAMATE RECEPTOR MGL-1"
assert record.rounds[1].alignments[12].length==999
assert record.rounds[1].alignments[13].title==">gi|6014748|sp|Q10900|CTPI_MYCTU PROBABLE CATION-TRANSPORTING ATPASE I"
assert record.rounds[1].alignments[13].length==1632
assert record.rounds[1].alignments[14].title==">gi|2833574|sp|Q58178|Y768_METJA HYPOTHETICAL PROTEIN MJ0768"
assert record.rounds[1].alignments[14].length==249
assert record.rounds[1].alignments[15].title==">gi|1169299|sp|P45513|DHAT_CITFR 1,3-PROPANEDIOL DEHYDROGENASE (3-HYDROXYPROPIONALDEHYDE REDUCTASE) (1,3-PROPANEDIOL OXIDOREDUCTASE)"
assert record.rounds[1].alignments[15].length==387
assert record.rounds[1].alignments[16].title==">gi|2495080|sp|P70579|MGR8_RAT METABOTROPIC GLUTAMATE RECEPTOR 8 PRECURSOR"
assert record.rounds[1].alignments[16].length==908
assert record.rounds[1].alignments[17].title==">gi|2495079|sp|O00222|MGR8_HUMAN METABOTROPIC GLUTAMATE RECEPTOR 8 PRECURSOR"
assert record.rounds[1].alignments[17].length==908
assert record.rounds[1].alignments[18].title==">gi|1346533|sp|P47743|MGR8_MOUSE METABOTROPIC GLUTAMATE RECEPTOR 8 PRECURSOR"
assert record.rounds[1].alignments[18].length==908
assert record.rounds[1].alignments[19].title==">gi|113381|sp|P06758|ADH2_ZYMMO ALCOHOL DEHYDROGENASE II (ADH II)"
assert record.rounds[1].alignments[19].length==383
assert record.rounds[1].alignments[20].title==">gi|127751|sp|P02567|MYSD_CAEEL MYOSIN HEAVY CHAIN D (MHC D)"
assert record.rounds[1].alignments[20].length==1938
assert record.rounds[1].alignments[21].title==">gi|1346321|sp|P23897|HSER_RAT HEAT-STABLE ENTEROTOXIN RECEPTOR PRECURSOR (GC-C) (INTESTINAL GUANYLATE CYCLASE) (STA RECEPTOR)"
assert record.rounds[1].alignments[21].length==1072
assert record.rounds[1].alignments[22].title==">gi|1175674|sp|P45116|YCIH_HAEIN HYPOTHETICAL PROTEIN HI1225"
assert record.rounds[1].alignments[22].length==106
assert record.rounds[1].alignments[23].title==">gi|3025270|sp|P77269|YPHF_ECOLI ABC TRANSPORTER PERIPLASMIC BINDING PROTEIN YPHF PRECURSOR"
assert record.rounds[1].alignments[23].length==327
assert len(record.rounds[2].new_seqs)==16
assert record.rounds[2].new_seqs[0].title=="gi|3024134|sp|O15303|MGR6_HUMAN METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[2].new_seqs[0].score==42
assert record.rounds[2].new_seqs[0].e==0.002
assert record.rounds[2].new_seqs[1].title=="gi|2495077|sp|Q14833|MGR4_HUMAN METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[2].new_seqs[1].score==39
assert record.rounds[2].new_seqs[1].e==0.012
assert record.rounds[2].new_seqs[2].title=="gi|3024131|sp|P91685|MGR_DROME METABOTROPIC GLUTAMATE RECEPTOR ..."
assert record.rounds[2].new_seqs[2].score==38
assert record.rounds[2].new_seqs[2].e==0.021
assert record.rounds[2].new_seqs[3].title=="gi|400255|sp|P31423|MGR4_RAT METABOTROPIC GLUTAMATE RECEPTOR 4 ..."
assert record.rounds[2].new_seqs[3].score==38
assert record.rounds[2].new_seqs[3].e==0.036
assert record.rounds[2].new_seqs[4].title=="gi|547903|sp|P35349|MGR6_RAT METABOTROPIC GLUTAMATE RECEPTOR 6 ..."
assert record.rounds[2].new_seqs[4].score==37
assert record.rounds[2].new_seqs[4].e==0.047
assert record.rounds[2].new_seqs[5].title=="gi|2495080|sp|P70579|MGR8_RAT METABOTROPIC GLUTAMATE RECEPTOR 8..."
assert record.rounds[2].new_seqs[5].score==36
assert record.rounds[2].new_seqs[5].e==0.10
assert record.rounds[2].new_seqs[6].title=="gi|2495079|sp|O00222|MGR8_HUMAN METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[2].new_seqs[6].score==36
assert record.rounds[2].new_seqs[6].e==0.10
assert record.rounds[2].new_seqs[7].title=="gi|1346533|sp|P47743|MGR8_MOUSE METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[2].new_seqs[7].score==35
assert record.rounds[2].new_seqs[7].e==0.18
assert record.rounds[2].new_seqs[8].title=="gi|400402|sp|P13595|NCA1_MOUSE NEURAL CELL ADHESION MOLECULE, L..."
assert record.rounds[2].new_seqs[8].score==33
assert record.rounds[2].new_seqs[8].e==0.69
assert record.rounds[2].new_seqs[9].title=="gi|1170947|sp|P31424|MGR5_RAT METABOTROPIC GLUTAMATE RECEPTOR 5..."
assert record.rounds[2].new_seqs[9].score==33
assert record.rounds[2].new_seqs[9].e==0.69
assert record.rounds[2].new_seqs[10].title=="gi|1706242|sp|P51840|CYGE_RAT GUANYLYL CYCLASE GC-E PRECURSOR (..."
assert record.rounds[2].new_seqs[10].score==33
assert record.rounds[2].new_seqs[10].e==0.91
assert record.rounds[2].new_seqs[11].title=="gi|1709020|sp|P41594|MGR5_HUMAN METABOTROPIC GLUTAMATE RECEPTOR..."
assert record.rounds[2].new_seqs[11].score==33
assert record.rounds[2].new_seqs[11].e==0.91
assert record.rounds[2].new_seqs[12].title=="gi|549451|sp|Q06278|ADO_HUMAN ALDEHYDE OXIDASE"
assert record.rounds[2].new_seqs[12].score==32
assert record.rounds[2].new_seqs[12].e==2.0
assert record.rounds[2].new_seqs[13].title=="gi|6014748|sp|Q10900|CTPI_MYCTU PROBABLE CATION-TRANSPORTING AT..."
assert record.rounds[2].new_seqs[13].score==32
assert record.rounds[2].new_seqs[13].e==2.0
assert record.rounds[2].new_seqs[14].title=="gi|127861|sp|P13594|NCA2_MOUSE NEURAL CELL ADHESION MOLECULE, P..."
assert record.rounds[2].new_seqs[14].score==31
assert record.rounds[2].new_seqs[14].e==2.7
assert record.rounds[2].new_seqs[15].title=="gi|127971|sp||NCPR_SALTR_1 [Segment 1 of 3] NADPH-CYTOCHROME P4..."
assert record.rounds[2].new_seqs[15].score==31
assert record.rounds[2].new_seqs[15].e==3.5
assert len(record.rounds[2].alignments)==23
assert record.rounds[2].alignments[0].title==">gi|126342|sp|P17215|LIVJ_SALTY LEU/ILE/VAL/THR-BINDING PROTEIN PRECURSOR (LIVT-BP)"
assert record.rounds[2].alignments[0].length==365
assert record.rounds[2].alignments[1].title==">gi|126349|sp|P04816|LIVK_ECOLI LEUCINE-SPECIFIC BINDING PROTEIN PRECURSOR (LS-BP) (L-BP)"
assert record.rounds[2].alignments[1].length==369
assert record.rounds[2].alignments[2].title==">gi|126348|sp|P02917|LIVJ_ECOLI LEU/ILE/VAL-BINDING PROTEIN PRECURSOR (LIV-BP)"
assert record.rounds[2].alignments[2].length==367
assert record.rounds[2].alignments[3].title==">gi|126347|sp|P25399|LIVJ_CITFR LEU/ILE/VAL-BINDING PROTEIN PRECURSOR (LIV-BP)"
assert record.rounds[2].alignments[3].length==367
assert record.rounds[2].alignments[4].title==">gi|126343|sp|P17216|LIVK_SALTY LEUCINE-SPECIFIC BINDING PROTEIN PRECURSOR (LS-BP) (L-BP)"
assert record.rounds[2].alignments[4].length==369
assert record.rounds[2].alignments[5].title==">gi|115120|sp|P21175|BRAC_PSEAE LEUCINE-, ISOLEUCINE-, VALINE-, THREONINE-, AND ALANINE-BINDING PROTEIN PRECURSOR (LIVAT-BP) (LEU/ILE/VAL/THR/ALA-BINDING PROTEIN)"
assert record.rounds[2].alignments[5].length==373
assert record.rounds[2].alignments[6].title==">gi|113709|sp|P27017|AMIC_PSEAE ALIPHATIC AMIDASE EXPRESSION-REGULATING PROTEIN"
assert record.rounds[2].alignments[6].length==385
assert record.rounds[2].alignments[7].title==">gi|3024134|sp|O15303|MGR6_HUMAN METABOTROPIC GLUTAMATE RECEPTOR 6 PRECURSOR"
assert record.rounds[2].alignments[7].length==877
assert record.rounds[2].alignments[8].title==">gi|2495077|sp|Q14833|MGR4_HUMAN METABOTROPIC GLUTAMATE RECEPTOR 4 PRECURSOR"
assert record.rounds[2].alignments[8].length==912
assert record.rounds[2].alignments[9].title==">gi|3024131|sp|P91685|MGR_DROME METABOTROPIC GLUTAMATE RECEPTOR PRECURSOR"
assert record.rounds[2].alignments[9].length==976
assert record.rounds[2].alignments[10].title==">gi|400255|sp|P31423|MGR4_RAT METABOTROPIC GLUTAMATE RECEPTOR 4 PRECURSOR"
assert record.rounds[2].alignments[10].length==912
assert record.rounds[2].alignments[11].title==">gi|547903|sp|P35349|MGR6_RAT METABOTROPIC GLUTAMATE RECEPTOR 6 PRECURSOR"
assert record.rounds[2].alignments[11].length==871
assert record.rounds[2].alignments[12].title==">gi|2495080|sp|P70579|MGR8_RAT METABOTROPIC GLUTAMATE RECEPTOR 8 PRECURSOR"
assert record.rounds[2].alignments[12].length==908
assert record.rounds[2].alignments[13].title==">gi|2495079|sp|O00222|MGR8_HUMAN METABOTROPIC GLUTAMATE RECEPTOR 8 PRECURSOR"
assert record.rounds[2].alignments[13].length==908
assert record.rounds[2].alignments[14].title==">gi|1346533|sp|P47743|MGR8_MOUSE METABOTROPIC GLUTAMATE RECEPTOR 8 PRECURSOR"
assert record.rounds[2].alignments[14].length==908
assert record.rounds[2].alignments[15].title==">gi|400402|sp|P13595|NCA1_MOUSE NEURAL CELL ADHESION MOLECULE, LARGE ISOFORM PRECURSOR (N-CAM 180) (NCAM-180) [CONTAINS: N-CAM 140 (NCAM-140)]"
assert record.rounds[2].alignments[15].length==1115
assert record.rounds[2].alignments[16].title==">gi|1170947|sp|P31424|MGR5_RAT METABOTROPIC GLUTAMATE RECEPTOR 5 PRECURSOR"
assert record.rounds[2].alignments[16].length==1203
assert record.rounds[2].alignments[17].title==">gi|1706242|sp|P51840|CYGE_RAT GUANYLYL CYCLASE GC-E PRECURSOR (GUANYLATE CYCLASE 2E)"
assert record.rounds[2].alignments[17].length==1108
assert record.rounds[2].alignments[18].title==">gi|1709020|sp|P41594|MGR5_HUMAN METABOTROPIC GLUTAMATE RECEPTOR 5 PRECURSOR"
assert record.rounds[2].alignments[18].length==1212
assert record.rounds[2].alignments[19].title==">gi|549451|sp|Q06278|ADO_HUMAN ALDEHYDE OXIDASE"
assert record.rounds[2].alignments[19].length==1338
assert record.rounds[2].alignments[20].title==">gi|6014748|sp|Q10900|CTPI_MYCTU PROBABLE CATION-TRANSPORTING ATPASE I"
assert record.rounds[2].alignments[20].length==1632
assert record.rounds[2].alignments[21].title==">gi|127861|sp|P13594|NCA2_MOUSE NEURAL CELL ADHESION MOLECULE, PHOSPHATIDYLINOSITOL-LINKED ISOFORM PRECURSOR (N-CAM 120) (NCAM-120)"
assert record.rounds[2].alignments[21].length==725
assert record.rounds[2].alignments[22].title==">gi|127971|sp||NCPR_SALTR_1 [Segment 1 of 3] NADPH-CYTOCHROME P450 REDUCTASE (CPR)"
assert record.rounds[2].alignments[22].length==426
assert record.rounds[0].alignments[0].hsps[0].score==1897
assert record.rounds[0].alignments[0].hsps[0].bits==743
assert record.rounds[0].alignments[0].hsps[0].expect==0.0
assert len(record.rounds[0].alignments[0].hsps)==1
assert record.rounds[0].alignments[1].hsps[0].score==1733
assert record.rounds[0].alignments[1].hsps[0].bits==679
assert record.rounds[0].alignments[1].hsps[0].expect==0.0
assert len(record.rounds[0].alignments[1].hsps)==1
assert record.rounds[0].alignments[2].hsps[0].score==1482
assert record.rounds[0].alignments[2].hsps[0].bits==581
assert record.rounds[0].alignments[2].hsps[0].expect==1e-166
assert len(record.rounds[0].alignments[2].hsps)==1
assert record.rounds[0].alignments[3].hsps[0].score==1470
assert record.rounds[0].alignments[3].hsps[0].bits==577
assert record.rounds[0].alignments[3].hsps[0].expect==1e-164
assert len(record.rounds[0].alignments[3].hsps)==1
assert record.rounds[0].alignments[4].hsps[0].score==1453
assert record.rounds[0].alignments[4].hsps[0].bits==570
assert record.rounds[0].alignments[4].hsps[0].expect==1e-162
assert len(record.rounds[0].alignments[4].hsps)==1
assert record.rounds[0].alignments[5].hsps[0].score==1050
assert record.rounds[0].alignments[5].hsps[0].bits==413
assert record.rounds[0].alignments[5].hsps[0].expect==1e-115
assert len(record.rounds[0].alignments[5].hsps)==1
assert record.rounds[0].alignments[6].hsps[0].score==92
assert record.rounds[0].alignments[6].hsps[0].bits==40.2
assert record.rounds[0].alignments[6].hsps[0].expect==0.006
assert len(record.rounds[0].alignments[6].hsps)==1
assert record.rounds[0].alignments[7].hsps[0].score==72
assert record.rounds[0].alignments[7].hsps[0].bits==32.5
assert record.rounds[0].alignments[7].hsps[0].expect==1.4
assert len(record.rounds[0].alignments[7].hsps)==1
assert record.rounds[0].alignments[8].hsps[0].score==71
assert record.rounds[0].alignments[8].hsps[0].bits==32.1
assert record.rounds[0].alignments[8].hsps[0].expect==1.8
assert len(record.rounds[0].alignments[8].hsps)==1
assert record.rounds[0].alignments[9].hsps[0].score==71
assert record.rounds[0].alignments[9].hsps[0].bits==32.1
assert record.rounds[0].alignments[9].hsps[0].expect==1.8
assert len(record.rounds[0].alignments[9].hsps)==1
assert record.rounds[0].alignments[10].hsps[0].score==71
assert record.rounds[0].alignments[10].hsps[0].bits==32.1
assert record.rounds[0].alignments[10].hsps[0].expect==1.8
assert len(record.rounds[0].alignments[10].hsps)==1
assert record.rounds[0].alignments[11].hsps[0].score==68
assert record.rounds[0].alignments[11].hsps[0].bits==30.9
assert record.rounds[0].alignments[11].hsps[0].expect==4.0
assert len(record.rounds[0].alignments[11].hsps)==1
assert record.rounds[0].alignments[12].hsps[0].score==68
assert record.rounds[0].alignments[12].hsps[0].bits==30.9
assert record.rounds[0].alignments[12].hsps[0].expect==4.0
assert len(record.rounds[0].alignments[12].hsps)==1
assert record.rounds[0].alignments[13].hsps[0].score==67
assert record.rounds[0].alignments[13].hsps[0].bits==30.5
assert record.rounds[0].alignments[13].hsps[0].expect==5.2
assert record.rounds[1].alignments[0].hsps[0].score==1683
assert record.rounds[1].alignments[0].hsps[0].bits==660
assert record.rounds[1].alignments[0].hsps[0].expect==0.0
assert len(record.rounds[1].alignments[0].hsps)==1
assert record.rounds[1].alignments[1].hsps[0].score==1673
assert record.rounds[1].alignments[1].hsps[0].bits==656
assert record.rounds[1].alignments[1].hsps[0].expect==0.0
assert len(record.rounds[1].alignments[1].hsps)==1
assert record.rounds[1].alignments[2].hsps[0].score==1672
assert record.rounds[1].alignments[2].hsps[0].bits==655
assert record.rounds[1].alignments[2].hsps[0].expect==0.0
assert len(record.rounds[1].alignments[2].hsps)==1
assert record.rounds[1].alignments[3].hsps[0].score==1671
assert record.rounds[1].alignments[3].hsps[0].bits==655
assert record.rounds[1].alignments[3].hsps[0].expect==0.0
assert len(record.rounds[1].alignments[3].hsps)==1
assert record.rounds[1].alignments[4].hsps[0].score==1670
assert record.rounds[1].alignments[4].hsps[0].bits==655
assert record.rounds[1].alignments[4].hsps[0].expect==0.0
assert len(record.rounds[1].alignments[4].hsps)==1
assert record.rounds[1].alignments[5].hsps[0].score==1569
assert record.rounds[1].alignments[5].hsps[0].bits==615
assert record.rounds[1].alignments[5].hsps[0].expect==1e-176
assert len(record.rounds[1].alignments[5].hsps)==1
assert record.rounds[1].alignments[6].hsps[0].score==113
assert record.rounds[1].alignments[6].hsps[0].bits==48.6
assert record.rounds[1].alignments[6].hsps[0].expect==2e-05
assert len(record.rounds[1].alignments[6].hsps)==1
assert record.rounds[1].alignments[7].hsps[0].score==94
assert record.rounds[1].alignments[7].hsps[0].bits==41.2
assert record.rounds[1].alignments[7].hsps[0].expect==0.003
assert len(record.rounds[1].alignments[7].hsps)==1
assert record.rounds[1].alignments[8].hsps[0].score==89
assert record.rounds[1].alignments[8].hsps[0].bits==39.2
assert record.rounds[1].alignments[8].hsps[0].expect==0.012
assert len(record.rounds[1].alignments[8].hsps)==1
assert record.rounds[1].alignments[9].hsps[0].score==89
assert record.rounds[1].alignments[9].hsps[0].bits==39.2
assert record.rounds[1].alignments[9].hsps[0].expect==0.012
assert len(record.rounds[1].alignments[9].hsps)==1
assert record.rounds[1].alignments[10].hsps[0].score==81
assert record.rounds[1].alignments[10].hsps[0].bits==36.1
assert record.rounds[1].alignments[10].hsps[0].expect==0.11
assert len(record.rounds[1].alignments[10].hsps)==1
assert record.rounds[1].alignments[11].hsps[0].score==80
assert record.rounds[1].alignments[11].hsps[0].bits==35.7
assert record.rounds[1].alignments[11].hsps[0].expect==0.14
assert len(record.rounds[1].alignments[11].hsps)==1
assert record.rounds[1].alignments[12].hsps[0].score==79
assert record.rounds[1].alignments[12].hsps[0].bits==35.3
assert record.rounds[1].alignments[12].hsps[0].expect==0.18
assert len(record.rounds[1].alignments[12].hsps)==1
assert record.rounds[1].alignments[13].hsps[0].score==74
assert record.rounds[1].alignments[13].hsps[0].bits==33.4
assert record.rounds[1].alignments[13].hsps[0].expect==0.71
assert len(record.rounds[1].alignments[13].hsps)==1
assert record.rounds[1].alignments[14].hsps[0].score==71
assert record.rounds[1].alignments[14].hsps[0].bits==32.2
assert record.rounds[1].alignments[14].hsps[0].expect==1.6
assert len(record.rounds[1].alignments[14].hsps)==1
assert record.rounds[1].alignments[15].hsps[0].score==68
assert record.rounds[1].alignments[15].hsps[0].bits==31.1
assert record.rounds[1].alignments[15].hsps[0].expect==3.6
assert len(record.rounds[1].alignments[15].hsps)==1
assert record.rounds[1].alignments[16].hsps[0].score==68
assert record.rounds[1].alignments[16].hsps[0].bits==31.1
assert record.rounds[1].alignments[16].hsps[0].expect==3.6
assert len(record.rounds[1].alignments[16].hsps)==1
assert record.rounds[1].alignments[17].hsps[0].score==68
assert record.rounds[1].alignments[17].hsps[0].bits==31.1
assert record.rounds[1].alignments[17].hsps[0].expect==3.6
assert len(record.rounds[1].alignments[17].hsps)==1
assert record.rounds[1].alignments[18].hsps[0].score==66
assert record.rounds[1].alignments[18].hsps[0].bits==30.3
assert record.rounds[1].alignments[18].hsps[0].expect==6.2
assert len(record.rounds[1].alignments[18].hsps)==1
assert record.rounds[1].alignments[19].hsps[0].score==65
assert record.rounds[1].alignments[19].hsps[0].bits==29.9
assert record.rounds[1].alignments[19].hsps[0].expect==8.1
assert len(record.rounds[1].alignments[19].hsps)==1
assert record.rounds[1].alignments[20].hsps[0].score==65
assert record.rounds[1].alignments[20].hsps[0].bits==29.9
assert record.rounds[1].alignments[20].hsps[0].expect==8.1
assert len(record.rounds[1].alignments[20].hsps)==1
assert record.rounds[1].alignments[21].hsps[0].score==65
assert record.rounds[1].alignments[21].hsps[0].bits==29.9
assert record.rounds[1].alignments[21].hsps[0].expect==8.1
assert len(record.rounds[1].alignments[21].hsps)==1
assert record.rounds[1].alignments[22].hsps[0].score==65
assert record.rounds[1].alignments[22].hsps[0].bits==29.9
assert record.rounds[1].alignments[22].hsps[0].expect==8.1
assert len(record.rounds[1].alignments[22].hsps)==1
assert record.rounds[1].alignments[23].hsps[0].score==65
assert record.rounds[1].alignments[23].hsps[0].bits==29.9
assert record.rounds[1].alignments[23].hsps[0].expect==8.1
assert record.rounds[2].alignments[0].hsps[0].score==1478
assert record.rounds[2].alignments[0].hsps[0].bits==580
assert record.rounds[2].alignments[0].hsps[0].expect==1e-165
assert len(record.rounds[2].alignments[0].hsps)==1
assert record.rounds[2].alignments[1].hsps[0].score==1477
assert record.rounds[2].alignments[1].hsps[0].bits==579
assert record.rounds[2].alignments[1].hsps[0].expect==1e-165
assert len(record.rounds[2].alignments[1].hsps)==1
assert record.rounds[2].alignments[2].hsps[0].score==1475
assert record.rounds[2].alignments[2].hsps[0].bits==579
assert record.rounds[2].alignments[2].hsps[0].expect==1e-165
assert len(record.rounds[2].alignments[2].hsps)==1
assert record.rounds[2].alignments[3].hsps[0].score==1474
assert record.rounds[2].alignments[3].hsps[0].bits==578
assert record.rounds[2].alignments[3].hsps[0].expect==1e-165
assert len(record.rounds[2].alignments[3].hsps)==1
assert record.rounds[2].alignments[4].hsps[0].score==1472
assert record.rounds[2].alignments[4].hsps[0].bits==577
assert record.rounds[2].alignments[4].hsps[0].expect==1e-165
assert len(record.rounds[2].alignments[4].hsps)==1
assert record.rounds[2].alignments[5].hsps[0].score==1351
assert record.rounds[2].alignments[5].hsps[0].bits==530
assert record.rounds[2].alignments[5].hsps[0].expect==1e-150
assert len(record.rounds[2].alignments[5].hsps)==1
assert record.rounds[2].alignments[6].hsps[0].score==976
assert record.rounds[2].alignments[6].hsps[0].bits==384
assert record.rounds[2].alignments[6].hsps[0].expect==1e-106
assert len(record.rounds[2].alignments[6].hsps)==1
assert record.rounds[2].alignments[7].hsps[0].score==95
assert record.rounds[2].alignments[7].hsps[0].bits==41.6
assert record.rounds[2].alignments[7].hsps[0].expect==0.002
assert len(record.rounds[2].alignments[7].hsps)==1
assert record.rounds[2].alignments[8].hsps[0].score==89
assert record.rounds[2].alignments[8].hsps[0].bits==39.3
assert record.rounds[2].alignments[8].hsps[0].expect==0.012
assert len(record.rounds[2].alignments[8].hsps)==1
assert record.rounds[2].alignments[9].hsps[0].score==87
assert record.rounds[2].alignments[9].hsps[0].bits==38.5
assert record.rounds[2].alignments[9].hsps[0].expect==0.021
assert len(record.rounds[2].alignments[9].hsps)==1
assert record.rounds[2].alignments[10].hsps[0].score==85
assert record.rounds[2].alignments[10].hsps[0].bits==37.7
assert record.rounds[2].alignments[10].hsps[0].expect==0.036
assert len(record.rounds[2].alignments[10].hsps)==1
assert record.rounds[2].alignments[11].hsps[0].score==84
assert record.rounds[2].alignments[11].hsps[0].bits==37.3
assert record.rounds[2].alignments[11].hsps[0].expect==0.047
assert len(record.rounds[2].alignments[11].hsps)==1
assert record.rounds[2].alignments[12].hsps[0].score==81
assert record.rounds[2].alignments[12].hsps[0].bits==36.1
assert record.rounds[2].alignments[12].hsps[0].expect==0.10
assert len(record.rounds[2].alignments[12].hsps)==1
assert record.rounds[2].alignments[13].hsps[0].score==81
assert record.rounds[2].alignments[13].hsps[0].bits==36.1
assert record.rounds[2].alignments[13].hsps[0].expect==0.10
assert len(record.rounds[2].alignments[13].hsps)==1
assert record.rounds[2].alignments[14].hsps[0].score==79
assert record.rounds[2].alignments[14].hsps[0].bits==35.4
assert record.rounds[2].alignments[14].hsps[0].expect==0.18
assert len(record.rounds[2].alignments[14].hsps)==1
assert record.rounds[2].alignments[15].hsps[0].score==74
assert record.rounds[2].alignments[15].hsps[0].bits==33.4
assert record.rounds[2].alignments[15].hsps[0].expect==0.69
assert len(record.rounds[2].alignments[15].hsps)==1
assert record.rounds[2].alignments[16].hsps[0].score==74
assert record.rounds[2].alignments[16].hsps[0].bits==33.4
assert record.rounds[2].alignments[16].hsps[0].expect==0.69
assert len(record.rounds[2].alignments[16].hsps)==1
assert record.rounds[2].alignments[17].hsps[0].score==73
assert record.rounds[2].alignments[17].hsps[0].bits==33.0
assert record.rounds[2].alignments[17].hsps[0].expect==0.91
assert len(record.rounds[2].alignments[17].hsps)==1
assert record.rounds[2].alignments[18].hsps[0].score==73
assert record.rounds[2].alignments[18].hsps[0].bits==33.0
assert record.rounds[2].alignments[18].hsps[0].expect==0.91
assert len(record.rounds[2].alignments[18].hsps)==1
assert record.rounds[2].alignments[19].hsps[0].score==70
assert record.rounds[2].alignments[19].hsps[0].bits==31.9
assert record.rounds[2].alignments[19].hsps[0].expect==2.0
assert len(record.rounds[2].alignments[19].hsps)==1
assert record.rounds[2].alignments[20].hsps[0].score==70
assert record.rounds[2].alignments[20].hsps[0].bits==31.9
assert record.rounds[2].alignments[20].hsps[0].expect==2.0
assert len(record.rounds[2].alignments[20].hsps)==1
assert record.rounds[2].alignments[21].hsps[0].score==69
assert record.rounds[2].alignments[21].hsps[0].bits==31.5
assert record.rounds[2].alignments[21].hsps[0].expect==2.7
assert len(record.rounds[2].alignments[21].hsps)==1
assert record.rounds[2].alignments[22].hsps[0].score==68
assert record.rounds[2].alignments[22].hsps[0].bits==31.1
assert record.rounds[2].alignments[22].hsps[0].expect==3.5
assert len(record.rounds[2].alignments[22].hsps)==1
assert record.rounds[0].alignments[0].hsps[0].identities==(369, 369)
assert record.rounds[0].alignments[0].hsps[0].positives==(369, 369)
assert record.rounds[0].alignments[1].hsps[0].identities==(333, 369)
assert record.rounds[0].alignments[1].hsps[0].positives==(351, 369)
assert record.rounds[0].alignments[2].hsps[0].identities==(282, 369)
assert record.rounds[0].alignments[2].hsps[0].positives==(320, 369)
assert record.rounds[0].alignments[2].hsps[0].gaps==(2, 369)
assert record.rounds[0].alignments[3].hsps[0].identities==(280, 366)
assert record.rounds[0].alignments[3].hsps[0].positives==(316, 366)
assert record.rounds[0].alignments[3].hsps[0].gaps==(2, 366)
assert record.rounds[0].alignments[4].hsps[0].identities==(276, 369)
assert record.rounds[0].alignments[4].hsps[0].positives==(316, 369)
assert record.rounds[0].alignments[4].hsps[0].gaps==(2, 369)
assert record.rounds[0].alignments[5].hsps[0].identities==(193, 368)
assert record.rounds[0].alignments[5].hsps[0].positives==(268, 368)
assert record.rounds[0].alignments[6].hsps[0].identities==(60, 334)
assert record.rounds[0].alignments[6].hsps[0].positives==(130, 334)
assert record.rounds[0].alignments[6].hsps[0].gaps==(23, 334)
assert record.rounds[0].alignments[7].hsps[0].identities==(18, 67)
assert record.rounds[0].alignments[7].hsps[0].positives==(32, 67)
assert record.rounds[0].alignments[7].hsps[0].gaps==(2, 67)
assert record.rounds[0].alignments[8].hsps[0].identities==(28, 89)
assert record.rounds[0].alignments[8].hsps[0].positives==(42, 89)
assert record.rounds[0].alignments[8].hsps[0].gaps==(11, 89)
assert record.rounds[0].alignments[9].hsps[0].identities==(29, 125)
assert record.rounds[0].alignments[9].hsps[0].positives==(57, 125)
assert record.rounds[0].alignments[9].hsps[0].gaps==(6, 125)
assert record.rounds[0].alignments[10].hsps[0].identities==(24, 93)
assert record.rounds[0].alignments[10].hsps[0].positives==(44, 93)
assert record.rounds[0].alignments[10].hsps[0].gaps==(2, 93)
assert record.rounds[0].alignments[11].hsps[0].identities==(24, 65)
assert record.rounds[0].alignments[11].hsps[0].positives==(31, 65)
assert record.rounds[0].alignments[11].hsps[0].gaps==(9, 65)
assert record.rounds[0].alignments[12].hsps[0].identities==(18, 55)
assert record.rounds[0].alignments[12].hsps[0].positives==(31, 55)
assert record.rounds[0].alignments[12].hsps[0].gaps==(1, 55)
assert record.rounds[0].alignments[13].hsps[0].identities==(24, 90)
assert record.rounds[0].alignments[13].hsps[0].positives==(41, 90)
assert record.rounds[0].alignments[13].hsps[0].gaps==(10, 90)
assert record.rounds[1].alignments[0].hsps[0].identities==(333, 369)
assert record.rounds[1].alignments[0].hsps[0].positives==(351, 369)
assert record.rounds[1].alignments[1].hsps[0].identities==(280, 367)
assert record.rounds[1].alignments[1].hsps[0].positives==(316, 367)
assert record.rounds[1].alignments[1].hsps[0].gaps==(2, 367)
assert record.rounds[1].alignments[2].hsps[0].identities==(276, 369)
assert record.rounds[1].alignments[2].hsps[0].positives==(316, 369)
assert record.rounds[1].alignments[2].hsps[0].gaps==(2, 369)
assert record.rounds[1].alignments[3].hsps[0].identities==(369, 369)
assert record.rounds[1].alignments[3].hsps[0].positives==(369, 369)
assert record.rounds[1].alignments[4].hsps[0].identities==(282, 369)
assert record.rounds[1].alignments[4].hsps[0].positives==(320, 369)
assert record.rounds[1].alignments[4].hsps[0].gaps==(2, 369)
assert record.rounds[1].alignments[5].hsps[0].identities==(193, 368)
assert record.rounds[1].alignments[5].hsps[0].positives==(268, 368)
assert record.rounds[1].alignments[6].hsps[0].identities==(60, 314)
assert record.rounds[1].alignments[6].hsps[0].positives==(114, 314)
assert record.rounds[1].alignments[6].hsps[0].gaps==(12, 314)
assert record.rounds[1].alignments[7].hsps[0].identities==(36, 154)
assert record.rounds[1].alignments[7].hsps[0].positives==(66, 154)
assert record.rounds[1].alignments[7].hsps[0].gaps==(5, 154)
assert record.rounds[1].alignments[8].hsps[0].identities==(32, 121)
assert record.rounds[1].alignments[8].hsps[0].positives==(55, 121)
assert record.rounds[1].alignments[8].hsps[0].gaps==(4, 121)
assert record.rounds[1].alignments[9].hsps[0].identities==(32, 121)
assert record.rounds[1].alignments[9].hsps[0].positives==(55, 121)
assert record.rounds[1].alignments[9].hsps[0].gaps==(4, 121)
assert record.rounds[1].alignments[10].hsps[0].identities==(30, 121)
assert record.rounds[1].alignments[10].hsps[0].positives==(54, 121)
assert record.rounds[1].alignments[10].hsps[0].gaps==(4, 121)
assert record.rounds[1].alignments[11].hsps[0].identities==(30, 120)
assert record.rounds[1].alignments[11].hsps[0].positives==(53, 120)
assert record.rounds[1].alignments[11].hsps[0].gaps==(4, 120)
assert record.rounds[1].alignments[12].hsps[0].identities==(29, 125)
assert record.rounds[1].alignments[12].hsps[0].positives==(57, 125)
assert record.rounds[1].alignments[12].hsps[0].gaps==(6, 125)
assert record.rounds[1].alignments[13].hsps[0].identities==(20, 68)
assert record.rounds[1].alignments[13].hsps[0].positives==(26, 68)
assert record.rounds[1].alignments[13].hsps[0].gaps==(7, 68)
assert record.rounds[1].alignments[14].hsps[0].identities==(33, 171)
assert record.rounds[1].alignments[14].hsps[0].positives==(67, 171)
assert record.rounds[1].alignments[14].hsps[0].gaps==(10, 171)
assert record.rounds[1].alignments[15].hsps[0].identities==(22, 88)
assert record.rounds[1].alignments[15].hsps[0].positives==(38, 88)
assert record.rounds[1].alignments[15].hsps[0].gaps==(8, 88)
assert record.rounds[1].alignments[16].hsps[0].identities==(50, 196)
assert record.rounds[1].alignments[16].hsps[0].positives==(77, 196)
assert record.rounds[1].alignments[16].hsps[0].gaps==(17, 196)
assert record.rounds[1].alignments[17].hsps[0].identities==(50, 196)
assert record.rounds[1].alignments[17].hsps[0].positives==(77, 196)
assert record.rounds[1].alignments[17].hsps[0].gaps==(17, 196)
assert record.rounds[1].alignments[18].hsps[0].identities==(50, 196)
assert record.rounds[1].alignments[18].hsps[0].positives==(76, 196)
assert record.rounds[1].alignments[18].hsps[0].gaps==(17, 196)
assert record.rounds[1].alignments[19].hsps[0].identities==(19, 69)
assert record.rounds[1].alignments[19].hsps[0].positives==(28, 69)
assert record.rounds[1].alignments[19].hsps[0].gaps==(3, 69)
assert record.rounds[1].alignments[20].hsps[0].identities==(18, 72)
assert record.rounds[1].alignments[20].hsps[0].positives==(31, 72)
assert record.rounds[1].alignments[20].hsps[0].gaps==(2, 72)
assert record.rounds[1].alignments[21].hsps[0].identities==(16, 85)
assert record.rounds[1].alignments[21].hsps[0].positives==(35, 85)
assert record.rounds[1].alignments[21].hsps[0].gaps==(2, 85)
assert record.rounds[1].alignments[22].hsps[0].identities==(11, 57)
assert record.rounds[1].alignments[22].hsps[0].positives==(27, 57)
assert record.rounds[1].alignments[22].hsps[0].gaps==(3, 57)
assert record.rounds[1].alignments[23].hsps[0].identities==(19, 79)
assert record.rounds[1].alignments[23].hsps[0].positives==(31, 79)
assert record.rounds[1].alignments[23].hsps[0].gaps==(2, 79)
assert record.rounds[2].alignments[0].hsps[0].identities==(280, 367)
assert record.rounds[2].alignments[0].hsps[0].positives==(316, 367)
assert record.rounds[2].alignments[0].hsps[0].gaps==(2, 367)
assert record.rounds[2].alignments[1].hsps[0].identities==(333, 369)
assert record.rounds[2].alignments[1].hsps[0].positives==(351, 369)
assert record.rounds[2].alignments[2].hsps[0].identities==(282, 369)
assert record.rounds[2].alignments[2].hsps[0].positives==(320, 369)
assert record.rounds[2].alignments[2].hsps[0].gaps==(2, 369)
assert record.rounds[2].alignments[3].hsps[0].identities==(276, 369)
assert record.rounds[2].alignments[3].hsps[0].positives==(316, 369)
assert record.rounds[2].alignments[3].hsps[0].gaps==(2, 369)
assert record.rounds[2].alignments[4].hsps[0].identities==(369, 369)
assert record.rounds[2].alignments[4].hsps[0].positives==(369, 369)
assert record.rounds[2].alignments[5].hsps[0].identities==(193, 368)
assert record.rounds[2].alignments[5].hsps[0].positives==(268, 368)
assert record.rounds[2].alignments[6].hsps[0].identities==(60, 314)
assert record.rounds[2].alignments[6].hsps[0].positives==(114, 314)
assert record.rounds[2].alignments[6].hsps[0].gaps==(12, 314)
assert record.rounds[2].alignments[7].hsps[0].identities==(47, 266)
assert record.rounds[2].alignments[7].hsps[0].positives==(90, 266)
assert record.rounds[2].alignments[7].hsps[0].gaps==(20, 266)
assert record.rounds[2].alignments[8].hsps[0].identities==(43, 216)
assert record.rounds[2].alignments[8].hsps[0].positives==(78, 216)
assert record.rounds[2].alignments[8].hsps[0].gaps==(13, 216)
assert record.rounds[2].alignments[9].hsps[0].identities==(42, 182)
assert record.rounds[2].alignments[9].hsps[0].positives==(74, 182)
assert record.rounds[2].alignments[9].hsps[0].gaps==(8, 182)
assert record.rounds[2].alignments[10].hsps[0].identities==(42, 215)
assert record.rounds[2].alignments[10].hsps[0].positives==(76, 215)
assert record.rounds[2].alignments[10].hsps[0].gaps==(11, 215)
assert record.rounds[2].alignments[11].hsps[0].identities==(46, 266)
assert record.rounds[2].alignments[11].hsps[0].positives==(88, 266)
assert record.rounds[2].alignments[11].hsps[0].gaps==(20, 266)
assert record.rounds[2].alignments[12].hsps[0].identities==(46, 222)
assert record.rounds[2].alignments[12].hsps[0].positives==(82, 222)
assert record.rounds[2].alignments[12].hsps[0].gaps==(15, 222)
assert record.rounds[2].alignments[13].hsps[0].identities==(46, 222)
assert record.rounds[2].alignments[13].hsps[0].positives==(82, 222)
assert record.rounds[2].alignments[13].hsps[0].gaps==(15, 222)
assert record.rounds[2].alignments[14].hsps[0].identities==(46, 222)
assert record.rounds[2].alignments[14].hsps[0].positives==(81, 222)
assert record.rounds[2].alignments[14].hsps[0].gaps==(15, 222)
assert record.rounds[2].alignments[15].hsps[0].identities==(46, 221)
assert record.rounds[2].alignments[15].hsps[0].positives==(73, 221)
assert record.rounds[2].alignments[15].hsps[0].gaps==(24, 221)
assert record.rounds[2].alignments[16].hsps[0].identities==(31, 155)
assert record.rounds[2].alignments[16].hsps[0].positives==(58, 155)
assert record.rounds[2].alignments[16].hsps[0].gaps==(4, 155)
assert record.rounds[2].alignments[17].hsps[0].identities==(29, 193)
assert record.rounds[2].alignments[17].hsps[0].positives==(67, 193)
assert record.rounds[2].alignments[17].hsps[0].gaps==(12, 193)
assert record.rounds[2].alignments[18].hsps[0].identities==(31, 155)
assert record.rounds[2].alignments[18].hsps[0].positives==(57, 155)
assert record.rounds[2].alignments[18].hsps[0].gaps==(4, 155)
assert record.rounds[2].alignments[19].hsps[0].identities==(26, 105)
assert record.rounds[2].alignments[19].hsps[0].positives==(37, 105)
assert record.rounds[2].alignments[19].hsps[0].gaps==(3, 105)
assert record.rounds[2].alignments[20].hsps[0].identities==(16, 55)
assert record.rounds[2].alignments[20].hsps[0].positives==(22, 55)
assert record.rounds[2].alignments[20].hsps[0].gaps==(3, 55)
assert record.rounds[2].alignments[21].hsps[0].identities==(39, 188)
assert record.rounds[2].alignments[21].hsps[0].positives==(62, 188)
assert record.rounds[2].alignments[21].hsps[0].gaps==(23, 188)
assert record.rounds[2].alignments[22].hsps[0].identities==(28, 146)
assert record.rounds[2].alignments[22].hsps[0].positives==(53, 146)
assert record.rounds[2].alignments[22].hsps[0].gaps==(8, 146)
assert record.rounds[0].alignments[0].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[0].alignments[0].hsps[0].match=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[0].alignments[0].hsps[0].sbjct=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[0].alignments[0].hsps[0].query_start==1
assert record.rounds[0].alignments[0].hsps[0].query_end==369
assert record.rounds[0].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[0].hsps[0].sbjct_end==369
assert record.rounds[0].alignments[1].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[0].alignments[1].hsps[0].match=="MKR AKTIIAG++ALA+S  AMADDIKVA+VGAMSGP+AQWG MEFNGA QAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGI+YVIGHLCSSSTQPASDIYEDEGILMISPGAT PELTQRGYQ+IMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLK  NAN+VFFDGITAGEKDFSALIARL+KENIDFVYYGGYYPEMGQ++RQAR+ GLKTQFMGPEGVGNASLSNIAG AAEGMLVTMPKRYDQDPAN+ IV+ALKADKKDPSGPYVWITYAAVQSLATA+ R+ S  PL LVKDLKANGA+TVIGPL WDEKGDLKGF+FGVFQWHADGSST AK"
assert record.rounds[0].alignments[1].hsps[0].sbjct=="MKRNAKTIIAGMIALAISHTAMADDIKVAVVGAMSGPIAQWGIMEFNGAEQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIKYVIGHLCSSSTQPASDIYEDEGILMISPGATAPELTQRGYQHIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKAANANVVFFDGITAGEKDFSALIARLKKENIDFVYYGGYYPEMGQMLRQARSVGLKTQFMGPEGVGNASLSNIAGDAAEGMLVTMPKRYDQDPANQGIVDALKADKKDPSGPYVWITYAAVQSLATALERTGSDEPLALVKDLKANGANTVIGPLNWDEKGDLKGFDFGVFQWHADGSSTAAK"
assert record.rounds[0].alignments[1].hsps[0].query_start==1
assert record.rounds[0].alignments[1].hsps[0].query_end==369
assert record.rounds[0].alignments[1].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[1].hsps[0].sbjct_end==369
assert record.rounds[0].alignments[2].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[0].alignments[2].hsps[0].match=="M  K K ++AG++ALA S  A+A+DIKVA+VGAMSGPVAQ+GD EF GA QA+ DINAKGGIKG+KL   +YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGYQ I+RT GLDS QGPTAAKYILE VKPQRIAI+HDKQQYGEGLAR+VQDGLK+GNAN+VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQARA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   + +S    P ++ K LKAN  DTV+GPL WDEKGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[0].alignments[2].hsps[0].sbjct=="MNTKGKALLAGLIALAFSNMALAEDIKVAVVGAMSGPVAQYGDQEFTGAEQAVADINAKGGIKGNKLQIAKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYQLILRTTGLDSDQGPTAAKYILEKVKPQRIAIVHDKQQYGEGLARAVQDGLKKGNANVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQARAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNQSDD--PAEIAKYLKANSVDTVMGPLTWDEKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[0].alignments[2].hsps[0].query_start==1
assert record.rounds[0].alignments[2].hsps[0].query_end==369
assert record.rounds[0].alignments[2].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[2].hsps[0].sbjct_end==367
assert record.rounds[0].alignments[3].hsps[0].query=="KAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[0].alignments[3].hsps[0].match=="K KT++AG +AL++S  A ADDIKVA+VGAMSGPVAQ+GD EF GA QAI DINAKGGIKGDKLV V+YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGY+ ++RT GLDS QGPTAAKYILE VKPQRIAIIHDKQQYGEGLAR+VQDGLK+G  N+VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQ+RA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   +  S    P ++ K LK    DTV+GPL WDEKGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[0].alignments[3].hsps[0].sbjct=="KGKTLLAGCIALSLSHMAFADDIKVAVVGAMSGPVAQYGDQEFTGAEQAIADINAKGGIKGDKLVAVKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYKLVLRTTGLDSDQGPTAAKYILEKVKPQRIAIIHDKQQYGEGLARAVQDGLKKGGVNVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQSRAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNHSDD--PAEIAKYLKGATVDTVMGPLSWDEKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[0].alignments[3].hsps[0].query_start==4
assert record.rounds[0].alignments[3].hsps[0].query_end==369
assert record.rounds[0].alignments[3].hsps[0].sbjct_start==2
assert record.rounds[0].alignments[3].hsps[0].sbjct_end==365
assert record.rounds[0].alignments[4].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[0].alignments[4].hsps[0].match=="M  K K ++AG +AL++S  A A DIKVA+VGAMSGPVAQ+GD EF GA QAI DINAKGG+KGDKLV V+YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGY  ++RT GLDS QGPTAAKYI+E VKP+RIAI+HDKQQYGEGLARSVQD LK+ NA++VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQARA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   + +S    P ++ K LKAN  +TV+GPL WD KGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[0].alignments[4].hsps[0].sbjct=="MNMKGKALLAGCIALSLSNMAFAKDIKVAVVGAMSGPVAQYGDQEFTGAEQAIADINAKGGVKGDKLVMVKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYNLVLRTTGLDSDQGPTAAKYIVEKVKPKRIAIVHDKQQYGEGLARSVQDNLKKANADVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQARAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNQSDD--PAEIAKYLKANSVETVMGPLSWDAKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[0].alignments[4].hsps[0].query_start==1
assert record.rounds[0].alignments[4].hsps[0].query_end==369
assert record.rounds[0].alignments[4].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[4].hsps[0].sbjct_end==367
assert record.rounds[0].alignments[5].hsps[0].query=="KRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[0].alignments[5].hsps[0].match=="+R ++   A  +A   S    AD IK+A+ G ++GPVAQ+GDM+  GA  AI+ IN  GG+ G +L GV YDDACDPKQAVAVANK+VNDG+++V+GH+CSSSTQPA+DIYEDEG+LMI+P AT PE+T RGY+ I RT GLD+ QGP A K+I E  K + IA++HDKQQYGEG+A  V+  ++     +  F+G+ AG+KDF+ALI++L+K  + FVY+GGY+PEMG ++RQA+  GL  +FMGPEGVGN+ ++ IAG A+EGML T+P+ ++QDP NKA+++A KA  +DPSG +V   Y+AV  +A  + ++    P  + + L+AN  +T  G L +DEKGDLK F+F V++WH D + T  K"
assert record.rounds[0].alignments[5].hsps[0].sbjct=="QRLSRLFAAMAIAGFASYSMAADTIKIALAGPVTGPVAQYGDMQRAGALMAIEQINKAGGVNGAQLEGVIYDDACDPKQAVAVANKVVNDGVKFVVGHVCSSSTQPATDIYEDEGVLMITPSATAPEITSRGYKLIFRTIGLDNMQGPVAGKFIAERYKVKTIAVLHDKQQYGEGIATEVKKTVEDAGIKVAVFEGLNAGDKDFNALISKLKKAGVQFVYFGGYHPEMGLLLRQAKQAGLDARFMGPEGVGNSEITAIAGDASEGMLATLPRAFEQDPKNKALIDAFKAKNQDPSGIFVLPAYSAVTVIAKGIEKAGEADPEKVAEALRANTFETPTGNLGFDEKGDLKNFDFTVYEWHKDATRTEVK"
assert record.rounds[0].alignments[5].hsps[0].query_start==2
assert record.rounds[0].alignments[5].hsps[0].query_end==369
assert record.rounds[0].alignments[5].hsps[0].sbjct_start==6
assert record.rounds[0].alignments[5].hsps[0].sbjct_end==373
assert record.rounds[0].alignments[6].hsps[0].query=="IVGAMSGPVAQWGDMEFN---GARQAIKDINAKGGIKGDKLVGVEYDDACDP-KQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQY---IMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIV--FFDGITAGEKDFSALIARLQKENIDFVY---YGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQD-PANKAIVEALKADKKDPSGPYVWITYAAVQS--LATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGD"
assert record.rounds[0].alignments[6].hsps[0].match=="++G +        D+E +   GA  A++ +N +GG+ G  +  +  D   DP +  +   + I N G+++++G   S + +    + E    L+  P          G++Y   I+      +      A Y++     +R+  I     Y       ++   +Q    ++   +  +   + D    + R+ +   D V+    G    E+ + + +   +G +   +       A ++ +    AEG +V  P     D PA++A V+A      + +    W   A  Q+  L  A   + + R  D+ + L     D   GP++ + + +"
assert record.rounds[0].alignments[6].hsps[0].sbjct=="LIGLLFSETGVTADIERSHAYGALLAVEQLNREGGVGGRPIETLSQDPGGDPDRYRLCAEDFIRNRGVRFLVGCYMSHTRKAVMPVVERADALLCYP------TPYEGFEYSPNIVYGGPAPNQNSAPLAAYLIRHY-GERVVFIGSDYIYPRESNHVMRHLYRQHGGTVLEEIYIPLYPSDDDLQRAVERIYQARADVVFSTVVGTGTAELYRAIARRYGDGRRPP-IASLTTSEAEVAKMESDVAEGQVVVAPYFSSIDTPASRAFVQACHGFFPENATITAWAEAAYWQTLLLGRAAQAAGNWRVEDVQRHLYDIDIDAPQGPVRVERQNN"
assert record.rounds[0].alignments[6].hsps[0].query_start==30
assert record.rounds[0].alignments[6].hsps[0].query_end==348
assert record.rounds[0].alignments[6].hsps[0].sbjct_start==9
assert record.rounds[0].alignments[6].hsps[0].sbjct_end==334
assert record.rounds[0].alignments[7].hsps[0].query=="GPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATA"
assert record.rounds[0].alignments[7].hsps[0].match=="G +G  + ++ + AG     +L  + K  ++DP N  +V  +KA KK+     +W  Y   +  A A"
assert record.rounds[0].alignments[7].hsps[0].sbjct=="GKQGEAHLAMRHYAGTVRYNVLNWLEK--NKDPLNDTVVSVMKASKKNDLLVEIWQDYTTQEEAAAA"
assert record.rounds[0].alignments[7].hsps[0].query_start==247
assert record.rounds[0].alignments[7].hsps[0].query_end==313
assert record.rounds[0].alignments[7].hsps[0].sbjct_start==570
assert record.rounds[0].alignments[7].hsps[0].sbjct_end==634
assert record.rounds[0].alignments[8].hsps[0].query=="ASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKP------QRIAIIHDKQQYGEGLARSVQDGLKQGNA"
assert record.rounds[0].alignments[8].hsps[0].match=="A  +Y D  +L++     N  L   G Q +M+       +G T    +L T +P      Q+I I+H+ QQ   GLAR V   L+Q +A"
assert record.rounds[0].alignments[8].hsps[0].sbjct=="ARAMYGDPCLLILDE--PNASLDSEGDQALMQAIVALQKRGATV---VLITHRPALTTLAQKILILHEGQQQRMGLARDVLTELQQRSA"
assert record.rounds[0].alignments[8].hsps[0].query_start==108
assert record.rounds[0].alignments[8].hsps[0].query_end==190
assert record.rounds[0].alignments[8].hsps[0].sbjct_start==478
assert record.rounds[0].alignments[8].hsps[0].sbjct_end==561
assert record.rounds[0].alignments[9].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSA----LIARLQKE"
assert record.rounds[0].alignments[9].hsps[0].match=="V+G   SS +   +++     I  +SP +TN +L+ +  ++Y  RT   D  Q   A   I    K   +++++   +YGE  A + +   ++    I   + I   ++ F+     L+ +LQ E"
assert record.rounds[0].alignments[9].hsps[0].sbjct=="VVGGSYSSVSVQLANLLRLFRIAQVSPASTNADLSDKNRFEYFARTVPSDDYQA-MAMVEIAVKFKWSYVSLVYSADEYGELGADAFKKEARKKGICIALEERIQNKKESFTESINNLVQKLQPE"
assert record.rounds[0].alignments[9].hsps[0].query_start==96
assert record.rounds[0].alignments[9].hsps[0].query_end==215
assert record.rounds[0].alignments[9].hsps[0].sbjct_start==196
assert record.rounds[0].alignments[9].hsps[0].sbjct_end==319
assert record.rounds[0].alignments[10].hsps[0].query=="KDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLK-TQFMGPEGVGN-ASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKK"
assert record.rounds[0].alignments[10].hsps[0].match=="+ F A+ A +  E++  +  GG+  EM  + R +  +G+  ++F+ P+ +      +   GG    +L T    Y    A   +VEA+  DKK"
assert record.rounds[0].alignments[10].hsps[0].sbjct=="RTFIAMEAGVSVEDVQAMLMGGHGDEMVPLPRFSTISGIPVSEFIAPDRLAQIVERTRKGGGEIVNLLKTGSAYYAPAAATAQMVEAVLKDKK"
assert record.rounds[0].alignments[10].hsps[0].query_start==203
assert record.rounds[0].alignments[10].hsps[0].query_end==293
assert record.rounds[0].alignments[10].hsps[0].sbjct_start==153
assert record.rounds[0].alignments[10].hsps[0].sbjct_end==245
assert record.rounds[0].alignments[11].hsps[0].query=="QQYGEGLARS-----VQDGLKQGNANIVFFDGITAGEKDFSALIARL-QKENIDFVYYG---GYY"
assert record.rounds[0].alignments[11].hsps[0].match=="Q YG GL  S     +   LK GNA +     IT G  +F+ L   L +K N  F+ YG   G+Y"
assert record.rounds[0].alignments[11].hsps[0].sbjct=="QSYGHGLVISDNFVSISKPLKVGNAQLGTDGNITGGSGNFANLNTTLNRKVNSGFITYGATSGWY"
assert record.rounds[0].alignments[11].hsps[0].query_start==171
assert record.rounds[0].alignments[11].hsps[0].query_end==226
assert record.rounds[0].alignments[11].hsps[0].sbjct_start==716
assert record.rounds[0].alignments[11].hsps[0].sbjct_end==780
assert record.rounds[0].alignments[12].hsps[0].query=="AIVEALKADKKDPSGPYVWITYA-AVQSLATAMTRSASHRPLDLVKDLKANGADT"
assert record.rounds[0].alignments[12].hsps[0].match=="AI  A +AD+ D    Y +  +A A++S+  A+  ++   P+D + DLKA   +T"
assert record.rounds[0].alignments[12].hsps[0].sbjct=="AIQVAKEADRIDGIEQYAFRAFADALESIPMALAENSGLAPIDALSDLKAKQIET"
assert record.rounds[0].alignments[12].hsps[0].query_start==283
assert record.rounds[0].alignments[12].hsps[0].query_end==336
assert record.rounds[0].alignments[12].hsps[0].sbjct_start==429
assert record.rounds[0].alignments[12].hsps[0].sbjct_end==483
assert record.rounds[0].alignments[13].hsps[0].query=="KQAVAVANKIVN-DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGA--TNPE-------LTQRGYQYIMRTAGLDSSQGPTAAKYILET"
assert record.rounds[0].alignments[13].hsps[0].match=="K  ++ A ++V  +G+ Y +GH  +S         +DEG ++  PG   TN E       +  + Y+  + +AG        A K+I ET"
assert record.rounds[0].alignments[13].hsps[0].sbjct=="KDVLSNAEEVVEANGLFYAVGHDPASGLVKGQVELDDEGYIITKPGTSFTNVEGVFACGDVQDKRYRQAITSAGSGCVAALEAEKFIAET"
assert record.rounds[0].alignments[13].hsps[0].query_start==79
assert record.rounds[0].alignments[13].hsps[0].query_end==158
assert record.rounds[0].alignments[13].hsps[0].sbjct_start==235
assert record.rounds[0].alignments[13].hsps[0].sbjct_end==324
assert record.rounds[1].alignments[0].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[1].alignments[0].hsps[0].match=="MKR AKTIIAG++ALA+S  AMADDIKVA+VGAMSGP+AQWG MEFNGA QAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGI+YVIGHLCSSSTQPASDIYEDEGILMISPGAT PELTQRGYQ+IMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLK  NAN+VFFDGITAGEKDFSALIARL+KENIDFVYYGGYYPEMGQ++RQAR+ GLKTQFMGPEGVGNASLSNIAG AAEGMLVTMPKRYDQDPAN+ IV+ALKADKKDPSGPYVWITYAAVQSLATA+ R+ S  PL LVKDLKANGA+TVIGPL WDEKGDLKGF+FGVFQWHADGSST AK"
assert record.rounds[1].alignments[0].hsps[0].sbjct=="MKRNAKTIIAGMIALAISHTAMADDIKVAVVGAMSGPIAQWGIMEFNGAEQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIKYVIGHLCSSSTQPASDIYEDEGILMISPGATAPELTQRGYQHIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKAANANVVFFDGITAGEKDFSALIARLKKENIDFVYYGGYYPEMGQMLRQARSVGLKTQFMGPEGVGNASLSNIAGDAAEGMLVTMPKRYDQDPANQGIVDALKADKKDPSGPYVWITYAAVQSLATALERTGSDEPLALVKDLKANGANTVIGPLNWDEKGDLKGFDFGVFQWHADGSSTAAK"
assert record.rounds[1].alignments[0].hsps[0].query_start==1
assert record.rounds[1].alignments[0].hsps[0].query_end==369
assert record.rounds[1].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[0].hsps[0].sbjct_end==369
assert record.rounds[1].alignments[1].hsps[0].query=="RKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[1].alignments[1].hsps[0].match==" K KT++AG +AL++S  A ADDIKVA+VGAMSGPVAQ+GD EF GA QAI DINAKGGIKGDKLV V+YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGY+ ++RT GLDS QGPTAAKYILE VKPQRIAIIHDKQQYGEGLAR+VQDGLK+G  N+VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQ+RA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   +  S    P ++ K LK    DTV+GPL WDEKGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[1].alignments[1].hsps[0].sbjct=="MKGKTLLAGCIALSLSHMAFADDIKVAVVGAMSGPVAQYGDQEFTGAEQAIADINAKGGIKGDKLVAVKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYKLVLRTTGLDSDQGPTAAKYILEKVKPQRIAIIHDKQQYGEGLARAVQDGLKKGGVNVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQSRAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNHSD--DPAEIAKYLKGATVDTVMGPLSWDEKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[1].alignments[1].hsps[0].query_start==3
assert record.rounds[1].alignments[1].hsps[0].query_end==369
assert record.rounds[1].alignments[1].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[1].hsps[0].sbjct_end==365
assert record.rounds[1].alignments[2].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[1].alignments[2].hsps[0].match=="M  K K ++AG +AL++S  A A DIKVA+VGAMSGPVAQ+GD EF GA QAI DINAKGG+KGDKLV V+YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGY  ++RT GLDS QGPTAAKYI+E VKP+RIAI+HDKQQYGEGLARSVQD LK+ NA++VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQARA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   + +S    P ++ K LKAN  +TV+GPL WD KGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[1].alignments[2].hsps[0].sbjct=="MNMKGKALLAGCIALSLSNMAFAKDIKVAVVGAMSGPVAQYGDQEFTGAEQAIADINAKGGVKGDKLVMVKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYNLVLRTTGLDSDQGPTAAKYIVEKVKPKRIAIVHDKQQYGEGLARSVQDNLKKANADVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQARAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNQSD--DPAEIAKYLKANSVETVMGPLSWDAKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[1].alignments[2].hsps[0].query_start==1
assert record.rounds[1].alignments[2].hsps[0].query_end==369
assert record.rounds[1].alignments[2].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[2].hsps[0].sbjct_end==367
assert record.rounds[1].alignments[3].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[1].alignments[3].hsps[0].match=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[1].alignments[3].hsps[0].sbjct=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[1].alignments[3].hsps[0].query_start==1
assert record.rounds[1].alignments[3].hsps[0].query_end==369
assert record.rounds[1].alignments[3].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[3].hsps[0].sbjct_end==369
assert record.rounds[1].alignments[4].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[1].alignments[4].hsps[0].match=="M  K K ++AG++ALA S  A+A+DIKVA+VGAMSGPVAQ+GD EF GA QA+ DINAKGGIKG+KL   +YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGYQ I+RT GLDS QGPTAAKYILE VKPQRIAI+HDKQQYGEGLAR+VQDGLK+GNAN+VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQARA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   + +S    P ++ K LKAN  DTV+GPL WDEKGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[1].alignments[4].hsps[0].sbjct=="MNTKGKALLAGLIALAFSNMALAEDIKVAVVGAMSGPVAQYGDQEFTGAEQAVADINAKGGIKGNKLQIAKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYQLILRTTGLDSDQGPTAAKYILEKVKPQRIAIVHDKQQYGEGLARAVQDGLKKGNANVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQARAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNQSD--DPAEIAKYLKANSVDTVMGPLTWDEKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[1].alignments[4].hsps[0].query_start==1
assert record.rounds[1].alignments[4].hsps[0].query_end==369
assert record.rounds[1].alignments[4].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[4].hsps[0].sbjct_end==367
assert record.rounds[1].alignments[5].hsps[0].query=="KRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[1].alignments[5].hsps[0].match=="+R ++   A  +A   S    AD IK+A+ G ++GPVAQ+GDM+  GA  AI+ IN  GG+ G +L GV YDDACDPKQAVAVANK+VNDG+++V+GH+CSSSTQPA+DIYEDEG+LMI+P AT PE+T RGY+ I RT GLD+ QGP A K+I E  K + IA++HDKQQYGEG+A  V+  ++     +  F+G+ AG+KDF+ALI++L+K  + FVY+GGY+PEMG ++RQA+  GL  +FMGPEGVGN+ ++ IAG A+EGML T+P+ ++QDP NKA+++A KA  +DPSG +V   Y+AV  +A  + ++    P  + + L+AN  +T  G L +DEKGDLK F+F V++WH D + T  K"
assert record.rounds[1].alignments[5].hsps[0].sbjct=="QRLSRLFAAMAIAGFASYSMAADTIKIALAGPVTGPVAQYGDMQRAGALMAIEQINKAGGVNGAQLEGVIYDDACDPKQAVAVANKVVNDGVKFVVGHVCSSSTQPATDIYEDEGVLMITPSATAPEITSRGYKLIFRTIGLDNMQGPVAGKFIAERYKVKTIAVLHDKQQYGEGIATEVKKTVEDAGIKVAVFEGLNAGDKDFNALISKLKKAGVQFVYFGGYHPEMGLLLRQAKQAGLDARFMGPEGVGNSEITAIAGDASEGMLATLPRAFEQDPKNKALIDAFKAKNQDPSGIFVLPAYSAVTVIAKGIEKAGEADPEKVAEALRANTFETPTGNLGFDEKGDLKNFDFTVYEWHKDATRTEVK"
assert record.rounds[1].alignments[5].hsps[0].query_start==2
assert record.rounds[1].alignments[5].hsps[0].query_end==369
assert record.rounds[1].alignments[5].hsps[0].sbjct_start==6
assert record.rounds[1].alignments[5].hsps[0].sbjct_end==373
assert record.rounds[1].alignments[6].hsps[0].query=="SGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIV-NDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIV--FFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQ-ARANGLKTQ-FMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQD-PANKAIVEALKADKKDPSGPYVW--ITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGP"
assert record.rounds[1].alignments[6].hsps[0].match=="+G  A        GA  A++ +N +GG+ G  +  +  D   DP +    A   + N G+++++G   S + +    + E    L+  P  T  E  +     +      + +  P AA  I      +R+  I     Y       ++   +Q    ++   +  +   + D    + R+ +   D V+         ++ R  AR  G   +  +       A ++ +    AEG +V  P     D PA++A V+A      + +    W    Y     L  A   + + R  D+ + L     D   GP"
assert record.rounds[1].alignments[6].hsps[0].sbjct=="TGVTADIERSHAYGALLAVEQLNREGGVGGRPIETLSQDPGGDPDRYRLCAEDFIRNRGVRFLVGCYMSHTRKAVMPVVERADALLCYP--TPYEGFEYSPNIVYGGPAPNQNSAPLAAYLI--RHYGERVVFIGSDYIYPRESNHVMRHLYRQHGGTVLEEIYIPLYPSDDDLQRAVERIYQARADVVFSTVVGTGTAELYRAIARRYGDGRRPPIASLTTSEAEVAKMESDVAEGQVVVAPYFSSIDTPASRAFVQACHGFFPENATITAWAEAAYWQTLLLGRAAQAAGNWRVEDVQRHLYDIDIDAPQGP"
assert record.rounds[1].alignments[6].hsps[0].query_start==35
assert record.rounds[1].alignments[6].hsps[0].query_end==340
assert record.rounds[1].alignments[6].hsps[0].sbjct_start==17
assert record.rounds[1].alignments[6].hsps[0].sbjct_end==326
assert record.rounds[1].alignments[7].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGI--TAGEKDFSALIARL-QKENIDFVYYGGYYPEMGQIVRQARANGLKTQF"
assert record.rounds[1].alignments[7].hsps[0].match=="VIG   SS +   +++     I  +SP +T   L+ +  +    RT   D+ Q   A   IL+      ++ IH +  YGE    ++     + N  I   + +   A +K F ++I++L +K N   V       +  +I++ A+   L   F"
assert record.rounds[1].alignments[7].hsps[0].sbjct=="VIGGSYSSVSLQVANLLRLFHIPQVSPASTAKTLSDKTRFDLFARTVPPDTFQS-VALVDILKNFNWSYVSTIHSEGSYGEYGIEALHKEATERNVCIAVAEKVPSAADDKVFDSIISKLQKKPNARGVVLFTRAEDARRILQAAKRANLSQPF"
assert record.rounds[1].alignments[7].hsps[0].query_start==96
assert record.rounds[1].alignments[7].hsps[0].query_end==245
assert record.rounds[1].alignments[7].hsps[0].sbjct_start==152
assert record.rounds[1].alignments[7].hsps[0].sbjct_end==304
assert record.rounds[1].alignments[8].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGI--TAGEKDFSALIARLQ"
assert record.rounds[1].alignments[8].hsps[0].match=="VIG   SS      ++ +   I  I+  AT+ +L+ +  Y+Y +R    D+ Q   A   I++      ++ +H +  YGE    + ++   Q    I   D I   AGEK F  L+ +L+"
assert record.rounds[1].alignments[8].hsps[0].sbjct=="VIGPGSSSVAIQVQNLLQLFDIPQIAYSATSIDLSDKTLYKYFLRVVPSDTLQA-RAMLDIVKRYNWTYVSAVHTEGNYGESGMDAFKELAAQEGLCIAHSDKIYSNAGEKSFDRLLRKLR"
assert record.rounds[1].alignments[8].hsps[0].query_start==96
assert record.rounds[1].alignments[8].hsps[0].query_end==213
assert record.rounds[1].alignments[8].hsps[0].sbjct_start==159
assert record.rounds[1].alignments[8].hsps[0].sbjct_end==278
assert record.rounds[1].alignments[9].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGI--TAGEKDFSALIARLQ"
assert record.rounds[1].alignments[9].hsps[0].match=="VIG   SS      ++ +   I  I+  AT+ +L+ +  Y+Y +R    D+ Q   A   I++      ++ +H +  YGE    + ++   Q    I   D I   AGEK F  L+ +L+"
assert record.rounds[1].alignments[9].hsps[0].sbjct=="VIGPGSSSVAIQVQNLLQLFDIPQIAYSATSIDLSDKTLYKYFLRVVPSDTLQA-RAMLDIVKRYNWTYVSAVHTEGNYGESGMDAFKELAAQEGLCIAHSDKIYSNAGEKSFDRLLRKLR"
assert record.rounds[1].alignments[9].hsps[0].query_start==96
assert record.rounds[1].alignments[9].hsps[0].query_end==213
assert record.rounds[1].alignments[9].hsps[0].sbjct_start==159
assert record.rounds[1].alignments[9].hsps[0].sbjct_end==278
assert record.rounds[1].alignments[10].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGI--TAGEKDFSALIARLQ"
assert record.rounds[1].alignments[10].hsps[0].match=="VIG   SS      ++ +   I  I+  AT+ +L+ +  ++Y MR    D+ Q   A   I++      ++ +H +  YGE    + +D   +    I     I   AGE+ F  L+ +L+"
assert record.rounds[1].alignments[10].hsps[0].sbjct=="VIGPGSSSVAIQVQNLLQLFNIPQIAYSATSMDLSDKTLFKYFMRVVPSDAQQA-RAMVDIVKRYNWTYVSAVHTEGNYGESGMEAFKDMSAKEGICIAHSYKIYSNAGEQSFDKLLKKLR"
assert record.rounds[1].alignments[10].hsps[0].query_start==96
assert record.rounds[1].alignments[10].hsps[0].query_end==213
assert record.rounds[1].alignments[10].hsps[0].sbjct_start==145
assert record.rounds[1].alignments[10].hsps[0].sbjct_end==264
assert record.rounds[1].alignments[11].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGI--TAGEKDFSALIARL"
assert record.rounds[1].alignments[11].hsps[0].match=="VIG   SS      ++ +   I  I+  AT+ +L+ +  ++Y MR    D+ Q   A   I++      ++ +H +  YGE    + +D   +    I     I   AGE+ F  L+ +L"
assert record.rounds[1].alignments[11].hsps[0].sbjct=="VIGPGSSSVAIQVQNLLQLFNIPQIAYSATSMDLSDKTLFKYFMRVVPSDAQQA-RAMVDIVKRYNWTYVSAVHTEGNYGESGMEAFKDMSAKEGICIAHSYKIYSNAGEQSFDKLLKKL"
assert record.rounds[1].alignments[11].hsps[0].query_start==96
assert record.rounds[1].alignments[11].hsps[0].query_end==212
assert record.rounds[1].alignments[11].hsps[0].sbjct_start==146
assert record.rounds[1].alignments[11].hsps[0].sbjct_end==264
assert record.rounds[1].alignments[12].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDF----SALIARLQKE"
assert record.rounds[1].alignments[12].hsps[0].match=="V+G   SS +   +++     I  +SP +TN +L+ +  ++Y  RT   D  Q   A   I    K   +++++   +YGE  A + +   ++    I   + I   ++ F    + L+ +LQ E"
assert record.rounds[1].alignments[12].hsps[0].sbjct=="VVGGSYSSVSVQLANLLRLFRIAQVSPASTNADLSDKNRFEYFARTVPSDDYQA-MAMVEIAVKFKWSYVSLVYSADEYGELGADAFKKEARKKGICIALEERIQNKKESFTESINNLVQKLQPE"
assert record.rounds[1].alignments[12].hsps[0].query_start==96
assert record.rounds[1].alignments[12].hsps[0].query_end==215
assert record.rounds[1].alignments[12].hsps[0].sbjct_start==196
assert record.rounds[1].alignments[12].hsps[0].sbjct_end==319
assert record.rounds[1].alignments[13].hsps[0].query=="AGIVALAVSQGAMADDIKVAI---VGAMSGPVAQWGDMEFNGARQAIKDINAKGG----IKGDKLVGV"
assert record.rounds[1].alignments[13].hsps[0].match=="A   A   S   +   I  AI    G+ +GPV Q+ +   NG+  A       GG      G  L GV"
assert record.rounds[1].alignments[13].hsps[0].sbjct=="AAAAAGEASHVVVGGSIDAAIDTAKGSRAGPVEQYVNQAANGSLIAAASALVAGGGTEDAAGAILAGV"
assert record.rounds[1].alignments[13].hsps[0].query_start==10
assert record.rounds[1].alignments[13].hsps[0].query_end==70
assert record.rounds[1].alignments[13].hsps[0].sbjct_start==315
assert record.rounds[1].alignments[13].hsps[0].sbjct_end==382
assert record.rounds[1].alignments[14].hsps[0].query=="AIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGL---DSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDF"
assert record.rounds[1].alignments[14].hsps[0].match=="AI+     G  K  K+V V  D+A   K+ V V    +    ++  G +C++S    S+IY+   IL  +P  +  ++ +   +   +  G+   DS   P     +   +    I  + D++   +   R ++         +     +  GE D    +  ++  N+ F"
assert record.rounds[1].alignments[14].hsps[0].sbjct=="AIELAKKTG--KDPKVVQVILDEA---KEIVKVGKNFIITETKH--GFVCANSGVDESNIYKGIKILPKNPDESAEKIRKEIEKLTGKRVGVIISDSVGRPFRKGAVGIAIGVSGILALWDRKGEKDLFGRELKTTEVAIADELASMANVVMGEADEGIPVVIIRGANVPF"
assert record.rounds[1].alignments[14].hsps[0].query_start==52
assert record.rounds[1].alignments[14].hsps[0].query_end==219
assert record.rounds[1].alignments[14].hsps[0].sbjct_start==67
assert record.rounds[1].alignments[14].hsps[0].sbjct_end==230
assert record.rounds[1].alignments[15].hsps[0].query=="GPTAAKYILETVKP---QRIAIIHDK--QQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSAL--IARLQKENID-FVYYGGYYP"
assert record.rounds[1].alignments[15].hsps[0].match=="GP A   + E  K    ++  ++ DK  +   +G        L++   ++V FDG+    KD +    +   +KE+ D  V  GG  P"
assert record.rounds[1].alignments[15].hsps[0].sbjct=="GPNAISVVGERCKLLGGKKALLVTDKGLRAIKDGAVDKTLTHLREAGIDVVVFDGVEPNPKDTNVRDGLEVFRKEHCDIIVTVGGGSP"
assert record.rounds[1].alignments[15].hsps[0].query_start==148
assert record.rounds[1].alignments[15].hsps[0].query_end==227
assert record.rounds[1].alignments[15].hsps[0].sbjct_start==17
assert record.rounds[1].alignments[15].hsps[0].sbjct_end==104
assert record.rounds[1].alignments[16].hsps[0].query=="DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQR-GYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGE-GLARSVQDGLKQGNANIVFFDGITAGEK--DFSALIAR-LQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFM--GPEGVGNASLSNIA-----GGAAEGMLVTMPKR"
assert record.rounds[1].alignments[16].hsps[0].match=="D I  VIG   SS +   ++I     I  IS  +T PEL+    Y +  R    DS Q   A   I+  +    ++ +  +  YGE G+    Q   + G   I     I    +  +F  +I R L+  N   V       ++ +I+  A+       F+  G +  G    S IA        AEG +  +PKR"
assert record.rounds[1].alignments[16].hsps[0].sbjct=="DKISGVIGAAASSVSIMVANILRLFKIPQISYASTAPELSDNTRYDFFSRVVPPDSYQA-QAMVDIVTALGWNYVSTLASEGNYGESGVEAFTQISREIGGVCIAQSQKIPREPRPGEFEKIIKRLLETPNARAVIMFANEDDIRRILEAAKKLNQSGHFLWIGSDSWG----SKIAPVYQQEEIAEGAVTILPKR"
assert record.rounds[1].alignments[16].hsps[0].query_start==91
assert record.rounds[1].alignments[16].hsps[0].query_end==274
assert record.rounds[1].alignments[16].hsps[0].sbjct_start==145
assert record.rounds[1].alignments[16].hsps[0].sbjct_end==335
assert record.rounds[1].alignments[17].hsps[0].query=="DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQR-GYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGE-GLARSVQDGLKQGNANIVFFDGITAGEK--DFSALIAR-LQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFM--GPEGVGNASLSNIA-----GGAAEGMLVTMPKR"
assert record.rounds[1].alignments[17].hsps[0].match=="D I  VIG   SS +   ++I     I  IS  +T PEL+    Y +  R    DS Q   A   I+  +    ++ +  +  YGE G+    Q   + G   I     I    +  +F  +I R L+  N   V       ++ +I+  A+       F+  G +  G    S IA        AEG +  +PKR"
assert record.rounds[1].alignments[17].hsps[0].sbjct=="DKISGVIGAAASSVSIMVANILRLFKIPQISYASTAPELSDNTRYDFFSRVVPPDSYQA-QAMVDIVTALGWNYVSTLASEGNYGESGVEAFTQISREIGGVCIAQSQKIPREPRPGEFEKIIKRLLETPNARAVIMFANEDDIRRILEAAKKLNQSGHFLWIGSDSWG----SKIAPVYQQEEIAEGAVTILPKR"
assert record.rounds[1].alignments[17].hsps[0].query_start==91
assert record.rounds[1].alignments[17].hsps[0].query_end==274
assert record.rounds[1].alignments[17].hsps[0].sbjct_start==145
assert record.rounds[1].alignments[17].hsps[0].sbjct_end==335
assert record.rounds[1].alignments[18].hsps[0].query=="DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQR-GYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGE-GLARSVQDGLKQGNANIVFFDGITAGEK--DFSALIAR-LQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFM--GPEGVGNASLSNIA-----GGAAEGMLVTMPKR"
assert record.rounds[1].alignments[18].hsps[0].match=="D I  VIG   SS +   ++I     I  IS  +T PEL+    Y +  R    DS Q   A   I+  +    ++ +  +  YGE G+    Q   + G   I     I    +  +F  +I R L+  N   V       ++  I+  A+       F+  G +  G    S IA        AEG +  +PKR"
assert record.rounds[1].alignments[18].hsps[0].sbjct=="DKISGVIGAAASSVSIMVANILRLFKIPQISYASTAPELSDNTRYDFFSRVVPPDSYQA-QAMVDIVTALGWNYVSTLASEGNYGESGVEAFTQISREIGGVCIAQSQKIPREPRPGEFEKIIKRLLETPNARAVIMFANEDDIRGILEAAKKLNQSGHFLWIGSDSWG----SKIAPVYQQEEIAEGAVTILPKR"
assert record.rounds[1].alignments[18].hsps[0].query_start==91
assert record.rounds[1].alignments[18].hsps[0].query_end==274
assert record.rounds[1].alignments[18].hsps[0].sbjct_start==145
assert record.rounds[1].alignments[18].hsps[0].sbjct_end==335
assert record.rounds[1].alignments[19].hsps[0].query=="QRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSAL--IARLQKENIDFVY-YGGYYP"
assert record.rounds[1].alignments[19].hsps[0].match=="+   I+ D      G+ + V D LK    N   +DG+       + L  +  L+  N DFV   GG  P"
assert record.rounds[1].alignments[19].hsps[0].sbjct=="KNALIVSDAFMNKSGVVKQVADLLKAQGINSAVYDGVMPNPTVTAVLEGLKILKDNNSDFVISLGGGSP"
assert record.rounds[1].alignments[19].hsps[0].query_start==162
assert record.rounds[1].alignments[19].hsps[0].query_end==227
assert record.rounds[1].alignments[19].hsps[0].sbjct_start==32
assert record.rounds[1].alignments[19].hsps[0].sbjct_end==100
assert record.rounds[1].alignments[20].hsps[0].query=="GPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSA"
assert record.rounds[1].alignments[20].hsps[0].match=="G +G  + ++ + AG     +L  + K   +DP N  +V  +KA KK+     +W  Y   +  A A     "
assert record.rounds[1].alignments[20].hsps[0].sbjct=="GKQGEAHLAMRHYAGTVRYNVLNWLEKN--KDPLNDTVVSVMKASKKNDLLVEIWQDYTTQEEAAAAAKAGG"
assert record.rounds[1].alignments[20].hsps[0].query_start==247
assert record.rounds[1].alignments[20].hsps[0].query_end==318
assert record.rounds[1].alignments[20].hsps[0].sbjct_start==570
assert record.rounds[1].alignments[20].hsps[0].sbjct_end==639
assert record.rounds[1].alignments[21].hsps[0].query=="MLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKG"
assert record.rounds[1].alignments[21].hsps[0].match=="+++T+P   ++  AN ++     +++ D S  Y+  T      L T +    S       +  +      + GP+  D+ GD+  "
assert record.rounds[1].alignments[21].hsps[0].sbjct=="LVLTLPP--EKFIANASVSGRFPSERSDFSLAYLEGTLLFGHMLQTFLENGESVTTPKFARAFRNLTFQGLEGPVTLDDSGDIDN"
assert record.rounds[1].alignments[21].hsps[0].query_start==267
assert record.rounds[1].alignments[21].hsps[0].query_end==351
assert record.rounds[1].alignments[21].hsps[0].sbjct_start==294
assert record.rounds[1].alignments[21].hsps[0].sbjct_end==376
assert record.rounds[1].alignments[22].hsps[0].query=="TVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQK"
assert record.rounds[1].alignments[22].hsps[0].match==" +K ++ +++  K   G+G+ R  +    +  A +    G+   +++   L A L+K"
assert record.rounds[1].alignments[22].hsps[0].sbjct=="RIKEEKASVVRPK---GDGVVRIQKQTSGRKGAGVSVITGLDLSDEELKKLAAELKK"
assert record.rounds[1].alignments[22].hsps[0].query_start==158
assert record.rounds[1].alignments[22].hsps[0].query_end==214
assert record.rounds[1].alignments[22].hsps[0].sbjct_start==14
assert record.rounds[1].alignments[22].hsps[0].sbjct_end==67
assert record.rounds[1].alignments[23].hsps[0].query=="AIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGE-KDFSALIARLQKENIDFVYYGGYYPE-MGQIVRQARANGL"
assert record.rounds[1].alignments[23].hsps[0].match=="AI  D Q Y  G+ + VQD  K  +  +   +    G+    S  +  L   N+D +           + VR+A   G+"
assert record.rounds[1].alignments[23].hsps[0].sbjct=="AIYLDTQGYYAGVRQGVQDAAKDSSVQVQLIETNAQGDISKESTFVDTLVARNVDAIILSAVSENGSSRTVRRASEAGI"
assert record.rounds[1].alignments[23].hsps[0].query_start==165
assert record.rounds[1].alignments[23].hsps[0].query_end==241
assert record.rounds[1].alignments[23].hsps[0].sbjct_start==35
assert record.rounds[1].alignments[23].hsps[0].sbjct_end==113
assert record.rounds[2].alignments[0].hsps[0].query=="RKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[2].alignments[0].hsps[0].match==" K KT++AG +AL++S  A ADDIKVA+VGAMSGPVAQ+GD EF GA QAI DINAKGGIKGDKLV V+YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGY+ ++RT GLDS QGPTAAKYILE VKPQRIAIIHDKQQYGEGLAR+VQDGLK+G  N+VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQ+RA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   +  S    P ++ K LK    DTV+GPL WDEKGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[2].alignments[0].hsps[0].sbjct=="MKGKTLLAGCIALSLSHMAFADDIKVAVVGAMSGPVAQYGDQEFTGAEQAIADINAKGGIKGDKLVAVKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYKLVLRTTGLDSDQGPTAAKYILEKVKPQRIAIIHDKQQYGEGLARAVQDGLKKGGVNVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQSRAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNHSD--DPAEIAKYLKGATVDTVMGPLSWDEKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[2].alignments[0].hsps[0].query_start==3
assert record.rounds[2].alignments[0].hsps[0].query_end==369
assert record.rounds[2].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[2].alignments[0].hsps[0].sbjct_end==365
assert record.rounds[2].alignments[1].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[2].alignments[1].hsps[0].match=="MKR AKTIIAG++ALA+S  AMADDIKVA+VGAMSGP+AQWG MEFNGA QAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGI+YVIGHLCSSSTQPASDIYEDEGILMISPGAT PELTQRGYQ+IMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLK  NAN+VFFDGITAGEKDFSALIARL+KENIDFVYYGGYYPEMGQ++RQAR+ GLKTQFMGPEGVGNASLSNIAG AAEGMLVTMPKRYDQDPAN+ IV+ALKADKKDPSGPYVWITYAAVQSLATA+ R+ S  PL LVKDLKANGA+TVIGPL WDEKGDLKGF+FGVFQWHADGSST AK"
assert record.rounds[2].alignments[1].hsps[0].sbjct=="MKRNAKTIIAGMIALAISHTAMADDIKVAVVGAMSGPIAQWGIMEFNGAEQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIKYVIGHLCSSSTQPASDIYEDEGILMISPGATAPELTQRGYQHIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKAANANVVFFDGITAGEKDFSALIARLKKENIDFVYYGGYYPEMGQMLRQARSVGLKTQFMGPEGVGNASLSNIAGDAAEGMLVTMPKRYDQDPANQGIVDALKADKKDPSGPYVWITYAAVQSLATALERTGSDEPLALVKDLKANGANTVIGPLNWDEKGDLKGFDFGVFQWHADGSSTAAK"
assert record.rounds[2].alignments[1].hsps[0].query_start==1
assert record.rounds[2].alignments[1].hsps[0].query_end==369
assert record.rounds[2].alignments[1].hsps[0].sbjct_start==1
assert record.rounds[2].alignments[1].hsps[0].sbjct_end==369
assert record.rounds[2].alignments[2].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[2].alignments[2].hsps[0].match=="M  K K ++AG++ALA S  A+A+DIKVA+VGAMSGPVAQ+GD EF GA QA+ DINAKGGIKG+KL   +YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGYQ I+RT GLDS QGPTAAKYILE VKPQRIAI+HDKQQYGEGLAR+VQDGLK+GNAN+VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQARA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   + +S    P ++ K LKAN  DTV+GPL WDEKGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[2].alignments[2].hsps[0].sbjct=="MNTKGKALLAGLIALAFSNMALAEDIKVAVVGAMSGPVAQYGDQEFTGAEQAVADINAKGGIKGNKLQIAKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYQLILRTTGLDSDQGPTAAKYILEKVKPQRIAIVHDKQQYGEGLARAVQDGLKKGNANVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQARAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNQSD--DPAEIAKYLKANSVDTVMGPLTWDEKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[2].alignments[2].hsps[0].query_start==1
assert record.rounds[2].alignments[2].hsps[0].query_end==369
assert record.rounds[2].alignments[2].hsps[0].sbjct_start==1
assert record.rounds[2].alignments[2].hsps[0].sbjct_end==367
assert record.rounds[2].alignments[3].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[2].alignments[3].hsps[0].match=="M  K K ++AG +AL++S  A A DIKVA+VGAMSGPVAQ+GD EF GA QAI DINAKGG+KGDKLV V+YDDACDPKQAVAVANK+VNDGI+YVIGHLCSSSTQPASDIYEDEGILMI+P AT PELT RGY  ++RT GLDS QGPTAAKYI+E VKP+RIAI+HDKQQYGEGLARSVQD LK+ NA++VFFDGITAGEKDFS L+ARL+KENIDFVYYGGY+PEMGQI+RQARA GLKTQFMGPEGV N SLSNIAG +AEG+LVT PK YDQ PANK IV+A+KA K+DPSG +VW TYAA+QSL   + +S    P ++ K LKAN  +TV+GPL WD KGDLKGFEFGVF WHA+G++T AK"
assert record.rounds[2].alignments[3].hsps[0].sbjct=="MNMKGKALLAGCIALSLSNMAFAKDIKVAVVGAMSGPVAQYGDQEFTGAEQAIADINAKGGVKGDKLVMVKYDDACDPKQAVAVANKVVNDGIKYVIGHLCSSSTQPASDIYEDEGILMITPAATAPELTARGYNLVLRTTGLDSDQGPTAAKYIVEKVKPKRIAIVHDKQQYGEGLARSVQDNLKKANADVVFFDGITAGEKDFSTLVARLKKENIDFVYYGGYHPEMGQILRQARAAGLKTQFMGPEGVANVSLSNIAGESAEGLLVTKPKNYDQVPANKPIVDAIKAKKQDPSGAFVWTTYAALQSLQAGLNQSD--DPAEIAKYLKANSVETVMGPLSWDAKGDLKGFEFGVFDWHANGTATDAK"
assert record.rounds[2].alignments[3].hsps[0].query_start==1
assert record.rounds[2].alignments[3].hsps[0].query_end==369
assert record.rounds[2].alignments[3].hsps[0].sbjct_start==1
assert record.rounds[2].alignments[3].hsps[0].sbjct_end==367
assert record.rounds[2].alignments[4].hsps[0].query=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[2].alignments[4].hsps[0].match=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[2].alignments[4].hsps[0].sbjct=="MKRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[2].alignments[4].hsps[0].query_start==1
assert record.rounds[2].alignments[4].hsps[0].query_end==369
assert record.rounds[2].alignments[4].hsps[0].sbjct_start==1
assert record.rounds[2].alignments[4].hsps[0].sbjct_end==369
assert record.rounds[2].alignments[5].hsps[0].query=="KRKAKTIIAGIVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDLKGFEFGVFQWHADGSSTVAK"
assert record.rounds[2].alignments[5].hsps[0].match=="+R ++   A  +A   S    AD IK+A+ G ++GPVAQ+GDM+  GA  AI+ IN  GG+ G +L GV YDDACDPKQAVAVANK+VNDG+++V+GH+CSSSTQPA+DIYEDEG+LMI+P AT PE+T RGY+ I RT GLD+ QGP A K+I E  K + IA++HDKQQYGEG+A  V+  ++     +  F+G+ AG+KDF+ALI++L+K  + FVY+GGY+PEMG ++RQA+  GL  +FMGPEGVGN+ ++ IAG A+EGML T+P+ ++QDP NKA+++A KA  +DPSG +V   Y+AV  +A  + ++    P  + + L+AN  +T  G L +DEKGDLK F+F V++WH D + T  K"
assert record.rounds[2].alignments[5].hsps[0].sbjct=="QRLSRLFAAMAIAGFASYSMAADTIKIALAGPVTGPVAQYGDMQRAGALMAIEQINKAGGVNGAQLEGVIYDDACDPKQAVAVANKVVNDGVKFVVGHVCSSSTQPATDIYEDEGVLMITPSATAPEITSRGYKLIFRTIGLDNMQGPVAGKFIAERYKVKTIAVLHDKQQYGEGIATEVKKTVEDAGIKVAVFEGLNAGDKDFNALISKLKKAGVQFVYFGGYHPEMGLLLRQAKQAGLDARFMGPEGVGNSEITAIAGDASEGMLATLPRAFEQDPKNKALIDAFKAKNQDPSGIFVLPAYSAVTVIAKGIEKAGEADPEKVAEALRANTFETPTGNLGFDEKGDLKNFDFTVYEWHKDATRTEVK"
assert record.rounds[2].alignments[5].hsps[0].query_start==2
assert record.rounds[2].alignments[5].hsps[0].query_end==369
assert record.rounds[2].alignments[5].hsps[0].sbjct_start==6
assert record.rounds[2].alignments[5].hsps[0].sbjct_end==373
assert record.rounds[2].alignments[6].hsps[0].query=="SGPVAQWGDMEFNGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIV-NDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIV--FFDGITAGEKDFSALIARLQKENIDFVYYGGYYPEMGQIVRQ-ARANGLKTQ-FMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQD-PANKAIVEALKADKKDPSGPYVW--ITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGP"
assert record.rounds[2].alignments[6].hsps[0].match=="+G  A        GA  A++ +N +GG+ G  +  +  D   DP +    A   + N G+++++G   S + +    + E    L+  P  T  E  +     +      + +  P AA  I      +R+  I     Y       ++   +Q    ++   +  +   + D    + R+ +   D V+         ++ R  AR  G   +  +       A ++ +    AEG +V  P     D PA++A V+A      + +    W    Y     L  A   + + R  D+ + L     D   GP"
assert record.rounds[2].alignments[6].hsps[0].sbjct=="TGVTADIERSHAYGALLAVEQLNREGGVGGRPIETLSQDPGGDPDRYRLCAEDFIRNRGVRFLVGCYMSHTRKAVMPVVERADALLCYP--TPYEGFEYSPNIVYGGPAPNQNSAPLAAYLI--RHYGERVVFIGSDYIYPRESNHVMRHLYRQHGGTVLEEIYIPLYPSDDDLQRAVERIYQARADVVFSTVVGTGTAELYRAIARRYGDGRRPPIASLTTSEAEVAKMESDVAEGQVVVAPYFSSIDTPASRAFVQACHGFFPENATITAWAEAAYWQTLLLGRAAQAAGNWRVEDVQRHLYDIDIDAPQGP"
assert record.rounds[2].alignments[6].hsps[0].query_start==35
assert record.rounds[2].alignments[6].hsps[0].query_end==340
assert record.rounds[2].alignments[6].hsps[0].sbjct_start==17
assert record.rounds[2].alignments[6].hsps[0].sbjct_end==326
assert record.rounds[2].alignments[7].hsps[0].query=="NGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELT-QRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIV---FFDGITAGEKDFSALIARLQKE-NIDFVYYGGYYPEMGQIVRQARANGLKTQFM----GPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWIT"
assert record.rounds[2].alignments[7].hsps[0].match=="    QA+  + A    +GD            P    A   ++V      V+G   SS +   +++     I  IS  +T PEL+    Y +  R    DS Q   A   I+  +    ++ +  +  YGE    +     ++     +             +FS +I RL +  N   +       ++ +++  AR   L   F+       G   + + ++    A G +  +PKR   D       +       + +   +W  "
assert record.rounds[2].alignments[7].hsps[0].sbjct=="YALEQALSFVQALIRGRGDGDEVGVRCPGGVPPLRPAPPERVVA-----VVGASASSVSIMVANVLRLFAIPQISYASTAPELSDSTRYDFFSRVVPPDSYQA-QAMVDIVRALGWNYVSTLASEGNYGESGVEAFVQISREAGGVCIAQSIKIPREPKPGEFSKVIRRLMETPNARGIIIFANEDDIRRVLEAARQANLTGHFLWVGSDSWGAKTSPILSLE-DVAVGAITILPKRASID----GFDQYFMTRSLENNRRNIWFA"
assert record.rounds[2].alignments[7].hsps[0].query_start==47
assert record.rounds[2].alignments[7].hsps[0].query_end==303
assert record.rounds[2].alignments[7].hsps[0].sbjct_start==104
assert record.rounds[2].alignments[7].hsps[0].sbjct_end==358
assert record.rounds[2].alignments[8].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGE-GLARSVQDGLKQGNANIV--FFDGITAGEKDFSALIAR-LQKENIDFVYYGGYYPEMGQIVRQARANGLKTQF--MGPEGVGN--ASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWIT"
assert record.rounds[2].alignments[8].hsps[0].match=="VIG   SS +   ++I     I  IS  +T P+L+               +    A   I+  +K   ++ +  +  YGE G+   +Q   + G   I             +F  +I R L+  N   V       ++ +++  AR       F  MG +  G+  A + ++    AEG +  +PKR       +       +   D +   +W  "
assert record.rounds[2].alignments[8].hsps[0].sbjct=="VIGASGSSVSIMVANILRLFKIPQISYASTAPDLSDNSRYDFFSRVVPSDTYQAQAMVDIVRALKWNYVSTVASEGSYGESGVEAFIQKSREDGGVCIAQSVKIPREPKAGEFDKIIRRLLETSNARAVIIFANEDDIRRVLEAARRANQTGHFFWMGSDSWGSKIAPVLHLE-EVAEGAVTILPKRMSV----RGFDRYFSSRTLDNNRRNIWFA"
assert record.rounds[2].alignments[8].hsps[0].query_start==96
assert record.rounds[2].alignments[8].hsps[0].query_end==303
assert record.rounds[2].alignments[8].hsps[0].sbjct_start==153
assert record.rounds[2].alignments[8].hsps[0].sbjct_end==363
assert record.rounds[2].alignments[9].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFD--GITAGEKDFSALIARL-QKENIDFVYYGGYYPEMGQIVRQARANGLKTQF--MGPEGVG-NASLSNIAGGAAEGMLVT"
assert record.rounds[2].alignments[9].hsps[0].match=="VIG   SS +   +++     I  +SP +T   L+ +  +    RT   D+ Q   A   IL+      ++ IH +  YGE    ++     + N  I   +     A +K F ++I++L +K N   V       +  +I++ A+   L   F  +  +G G    L       AEG +  "
assert record.rounds[2].alignments[9].hsps[0].sbjct=="VIGGSYSSVSLQVANLLRLFHIPQVSPASTAKTLSDKTRFDLFARTVPPDTFQS-VALVDILKNFNWSYVSTIHSEGSYGEYGIEALHKEATERNVCIAVAEKVPSAADDKVFDSIISKLQKKPNARGVVLFTRAEDARRILQAAKRANLSQPFHWIASDGWGKQQKLLEGLEDIAEGAITV"
assert record.rounds[2].alignments[9].hsps[0].query_start==96
assert record.rounds[2].alignments[9].hsps[0].query_end==270
assert record.rounds[2].alignments[9].hsps[0].sbjct_start==152
assert record.rounds[2].alignments[9].hsps[0].sbjct_end==332
assert record.rounds[2].alignments[10].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGE-GLARSVQDGLKQGNANIV--FFDGITAGEKDFSALIAR-LQKENIDFVYYGGYYPEMGQIVRQARANGLKTQF--MGPEGVGNASLSNIA-GGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWIT"
assert record.rounds[2].alignments[10].hsps[0].match=="VIG   SS +   ++I     I  IS  +T P+L+               +    A   I+  +K   ++ +  +  YGE G+   +Q   + G   I             +F  +I R L+  N   +       ++ +++  AR       F  MG +  G+ S   +     AEG +  +PKR       +       +   D +   +W  "
assert record.rounds[2].alignments[10].hsps[0].sbjct=="VIGASGSSVSIMVANILRLFKIPQISYASTAPDLSDNSRYDFFSRVVPSDTYQAQAMVDIVRALKWNYVSTLASEGSYGESGVEAFIQKSRENGGVCIAQSVKIPREPKTGEFDKIIKRLLETSNARGIIIFANEDDIRRVLEAARRANQTGHFFWMGSDSWGSKSAPVLRLEEVAEGAVTILPKRMSV----RGFDRYFSSRTLDNNRRNIWFA"
assert record.rounds[2].alignments[10].hsps[0].query_start==96
assert record.rounds[2].alignments[10].hsps[0].query_end==303
assert record.rounds[2].alignments[10].hsps[0].sbjct_start==153
assert record.rounds[2].alignments[10].hsps[0].sbjct_end==363
assert record.rounds[2].alignments[11].hsps[0].query=="NGARQAIKDINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELT-QRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIV---FFDGITAGEKDFSALIARLQKE-NIDFVYYGGYYPEMGQIVRQARANGLKTQFM----GPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWIT"
assert record.rounds[2].alignments[11].hsps[0].match=="    QA+  + A    +GD            P    A   ++V      V+G   SS +   +++     I  IS  +T PEL+    Y +  R    DS Q   A   I+  +    ++ +  +  YGE    +     ++     +             +F  +I RL +  N   +       ++ +++   R   L   F+       G   + + N+    A G +  +PKR   D       +       + +   +W  "
assert record.rounds[2].alignments[11].hsps[0].sbjct=="YALEQALSFVQALIRGRGDGDEASVRCPGGVPPLRSAPPERVVA-----VVGASASSVSIMVANVLRLFAIPQISYASTAPELSDSTRYDFFSRVVPPDSYQA-QAMVDIVRALGWNYVSTLASEGNYGESGVEAFVQISREAGGVCIAQSIKIPREPKPGEFHKVIRRLMETPNARGIIIFANEDDIRRVLEATRQANLTGHFLWVGSDSWGSKISPILNLE-EEAVGAITILPKRASID----GFDQYFMTRSLENNRRNIWFA"
assert record.rounds[2].alignments[11].hsps[0].query_start==47
assert record.rounds[2].alignments[11].hsps[0].query_end==303
assert record.rounds[2].alignments[11].hsps[0].sbjct_start==98
assert record.rounds[2].alignments[11].hsps[0].sbjct_end==352
assert record.rounds[2].alignments[12].hsps[0].query=="DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFF---DGITAGEKDFSALIAR-LQKENIDFVYYGGYYPEMGQIVRQARANGLKTQF--MGPEGVGNASLSNI--AGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWIT"
assert record.rounds[2].alignments[12].hsps[0].match=="D I  VIG   SS +   ++I     I  IS  +T PEL+    Y +  R    DS Q   A   I+  +    ++ +  +  YGE    +     ++     +             +F  +I R L+  N   V       ++ +I+  A+       F  +G +  G + ++ +      AEG +  +PKR   D  ++       A+ +      VW  "
assert record.rounds[2].alignments[12].hsps[0].sbjct=="DKISGVIGAAASSVSIMVANILRLFKIPQISYASTAPELSDNTRYDFFSRVVPPDSYQA-QAMVDIVTALGWNYVSTLASEGNYGESGVEAFTQISREIGGVCIAQSQKIPREPRPGEFEKIIKRLLETPNARAVIMFANEDDIRRILEAAKKLNQSGHFLWIGSDSWG-SKIAPVYQQEEIAEGAVTILPKRASIDGFDRYFRSRTLANNRRN----VWFA"
assert record.rounds[2].alignments[12].hsps[0].query_start==91
assert record.rounds[2].alignments[12].hsps[0].query_end==303
assert record.rounds[2].alignments[12].hsps[0].sbjct_start==145
assert record.rounds[2].alignments[12].hsps[0].sbjct_end==360
assert record.rounds[2].alignments[13].hsps[0].query=="DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFF---DGITAGEKDFSALIAR-LQKENIDFVYYGGYYPEMGQIVRQARANGLKTQF--MGPEGVGNASLSNI--AGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWIT"
assert record.rounds[2].alignments[13].hsps[0].match=="D I  VIG   SS +   ++I     I  IS  +T PEL+    Y +  R    DS Q   A   I+  +    ++ +  +  YGE    +     ++     +             +F  +I R L+  N   V       ++ +I+  A+       F  +G +  G + ++ +      AEG +  +PKR   D  ++       A+ +      VW  "
assert record.rounds[2].alignments[13].hsps[0].sbjct=="DKISGVIGAAASSVSIMVANILRLFKIPQISYASTAPELSDNTRYDFFSRVVPPDSYQA-QAMVDIVTALGWNYVSTLASEGNYGESGVEAFTQISREIGGVCIAQSQKIPREPRPGEFEKIIKRLLETPNARAVIMFANEDDIRRILEAAKKLNQSGHFLWIGSDSWG-SKIAPVYQQEEIAEGAVTILPKRASIDGFDRYFRSRTLANNRRN----VWFA"
assert record.rounds[2].alignments[13].hsps[0].query_start==91
assert record.rounds[2].alignments[13].hsps[0].query_end==303
assert record.rounds[2].alignments[13].hsps[0].sbjct_start==145
assert record.rounds[2].alignments[13].hsps[0].sbjct_end==360
assert record.rounds[2].alignments[14].hsps[0].query=="DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRG-YQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFF---DGITAGEKDFSALIAR-LQKENIDFVYYGGYYPEMGQIVRQARANGLKTQF--MGPEGVGNASLSNI--AGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWIT"
assert record.rounds[2].alignments[14].hsps[0].match=="D I  VIG   SS +   ++I     I  IS  +T PEL+    Y +  R    DS Q   A   I+  +    ++ +  +  YGE    +     ++     +             +F  +I R L+  N   V       ++  I+  A+       F  +G +  G + ++ +      AEG +  +PKR   D  ++       A+ +      VW  "
assert record.rounds[2].alignments[14].hsps[0].sbjct=="DKISGVIGAAASSVSIMVANILRLFKIPQISYASTAPELSDNTRYDFFSRVVPPDSYQA-QAMVDIVTALGWNYVSTLASEGNYGESGVEAFTQISREIGGVCIAQSQKIPREPRPGEFEKIIKRLLETPNARAVIMFANEDDIRGILEAAKKLNQSGHFLWIGSDSWG-SKIAPVYQQEEIAEGAVTILPKRASIDGFDRYFRSRTLANNRRN----VWFA"
assert record.rounds[2].alignments[14].hsps[0].query_start==91
assert record.rounds[2].alignments[14].hsps[0].query_end==303
assert record.rounds[2].alignments[14].hsps[0].sbjct_start==145
assert record.rounds[2].alignments[14].hsps[0].sbjct_end==360
assert record.rounds[2].alignments[15].hsps[0].query=="DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDG----------LKQGNANIVFFDGITA-GEKDFS----ALIARLQKENIDFVYYGGYYPEMGQIVR-QARANGLKTQFMGPEGVGNASLSNIAGG-AAEGMLVTMPKRYDQDPANKAIVEALKADKKD"
assert record.rounds[2].alignments[15].hsps[0].match=="DG   V  H   SS    S  Y D G  M +   T   + Q      +        QGP A      T +  ++ I  +   Y        +DG          +K  N     +  +T   E DF       + R+ +E+++F+      P    I R +  ++  + QF  PE  G   +        + G      K YD   AN   +  +   K +"
assert record.rounds[2].alignments[15].hsps[0].sbjct=="DGHMVVRSHARVSSLTLKSIQYRDAGEYMCTASNT---IGQDSQSIDLEFQYAPKLQGPVAVY----TWEGNQVNITCEVFAYPSATISWFRDGQLLPSSNYSNIKIYNTPSASYLEVTPDSENDFGNYNCTAVNRIGQESLEFILVQADTPSSPSIDRVEPYSSTAQVQFDEPEATGGVPILKYKAEWKSLGEESWHFKWYDAKEANMEGIVTIMGLKPE"
assert record.rounds[2].alignments[15].hsps[0].query_start==91
assert record.rounds[2].alignments[15].hsps[0].query_end==294
assert record.rounds[2].alignments[15].hsps[0].sbjct_start==357
assert record.rounds[2].alignments[15].hsps[0].sbjct_end==570
assert record.rounds[2].alignments[16].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGIT--AGEKDFSALIARLQK--ENIDFVYYGGYYPEMGQIVRQARANGLKTQFM"
assert record.rounds[2].alignments[16].hsps[0].match=="VIG   SS      ++ +   I  I+  AT+ +L+ +             +Q   A   I++      ++ +H +  YGE    + +D   +    I     I   AGE+ F  L+ +L+        V        +  ++   R  GL  +F+"
assert record.rounds[2].alignments[16].hsps[0].sbjct=="VIGPGSSSVAIQVQNLLQLFNIPQIAYSATSMDLSDKTLFKYFMRVVPSDAQQARAMVDIVKRYNWTYVSAVHTEGNYGESGMEAFKDMSAKEGICIAHSYKIYSNAGEQSFDKLLKKLRSHLPKARVVACFCEGMTVRGLLMAMRRLGLAGEFL"
assert record.rounds[2].alignments[16].hsps[0].query_start==96
assert record.rounds[2].alignments[16].hsps[0].query_end==246
assert record.rounds[2].alignments[16].hsps[0].sbjct_start==145
assert record.rounds[2].alignments[16].hsps[0].sbjct_end==299
assert record.rounds[2].alignments[17].hsps[0].query=="IVALAVSQGAMADDIKVAIVGAMSGPVAQWGDMEFNGARQAIKDINAKGGIKGDKL--VGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGITAGE"
assert record.rounds[2].alignments[17].hsps[0].match=="++ L  S  A +   KV ++G  +             AR A   +N    + G     V +  +    P    AV++ +    +  ++G +  ++ +PA  + ++ G+ ++  G                T     +    A   +L+  +  R+A+I   Q       R++   L+     +     +   +"
assert record.rounds[2].alignments[17].hsps[0].sbjct=="LLLLLPSPSAFSAVFKVGVLGPWACDPIFARARPDLAARLATDRLNRDLALDGGPWFEVTLLPEPCLTPGSLGAVSSALTR--VSGLVGPVNPAACRPAELLAQEAGVALVPWGCPGTRAA--------GTTAPAVTPAADALYVLLKAFRWARVALITAPQDLWVEAGRALSTALRARGLPVALVTSMVPSD"
assert record.rounds[2].alignments[17].hsps[0].query_start==12
assert record.rounds[2].alignments[17].hsps[0].query_end==202
assert record.rounds[2].alignments[17].hsps[0].sbjct_start==43
assert record.rounds[2].alignments[17].hsps[0].sbjct_end==225
assert record.rounds[2].alignments[18].hsps[0].query=="VIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDGLKQGNANIVFFDGIT--AGEKDFSALIARLQK--ENIDFVYYGGYYPEMGQIVRQARANGLKTQFM"
assert record.rounds[2].alignments[18].hsps[0].match=="VIG   SS      ++ +   I  I+  AT+ +L+ +             +Q   A   I++      ++ +H +  YGE    + +D   +    I     I   AGE+ F  L+ +L         V        +  ++   R  GL  +F+"
assert record.rounds[2].alignments[18].hsps[0].sbjct=="VIGPGSSSVAIQVQNLLQLFNIPQIAYSATSMDLSDKTLFKYFMRVVPSDAQQARAMVDIVKRYNWTYVSAVHTEGNYGESGMEAFKDMSAKEGICIAHSYKIYSNAGEQSFDKLLKKLTSHLPKARVVACFCEGMTVRGLLMAMRRLGLAGEFL"
assert record.rounds[2].alignments[18].hsps[0].query_start==96
assert record.rounds[2].alignments[18].hsps[0].query_end==246
assert record.rounds[2].alignments[18].hsps[0].sbjct_start==146
assert record.rounds[2].alignments[18].hsps[0].sbjct_end==300
assert record.rounds[2].alignments[19].hsps[0].query=="DINAKGGIKGDKLVGVEYDDACDPKQAVAVANKIVNDGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETV"
assert record.rounds[2].alignments[19].hsps[0].match=="D+N + G +G       Y  AC   +   +     N     V+   CS +  PA DI + EG  +   G    E      Q I+ T G D  + P       E  "
assert record.rounds[2].alignments[19].hsps[0].sbjct=="DMNWEKG-EGQPFEYFVYGAACSEVEIDCLTGDHKNIRTDIVMDVGCSIN--PAIDIGQIEGAFIQGMGLYTIEELNYSPQGILHTRGPDQYKIPAICDMPTELH"
assert record.rounds[2].alignments[19].hsps[0].query_start==55
assert record.rounds[2].alignments[19].hsps[0].query_end==159
assert record.rounds[2].alignments[19].hsps[0].sbjct_start==1148
assert record.rounds[2].alignments[19].hsps[0].sbjct_end==1249
assert record.rounds[2].alignments[20].hsps[0].query=="AGIVALAVSQGAMADDIKVAI---VGAMSGPVAQWGDMEFNGARQAIKDINAKGG"
assert record.rounds[2].alignments[20].hsps[0].match=="A   A   S   +   I  AI    G+ +GPV Q+ +   NG+  A       GG"
assert record.rounds[2].alignments[20].hsps[0].sbjct=="AAAAAGEASHVVVGGSIDAAIDTAKGSRAGPVEQYVNQAANGSLIAAASALVAGG"
assert record.rounds[2].alignments[20].hsps[0].query_start==10
assert record.rounds[2].alignments[20].hsps[0].query_end==61
assert record.rounds[2].alignments[20].hsps[0].sbjct_start==315
assert record.rounds[2].alignments[20].hsps[0].sbjct_end==369
assert record.rounds[2].alignments[21].hsps[0].query=="DGIQYVIGHLCSSSTQPASDIYEDEGILMISPGATNPELTQRGYQYIMRTAGLDSSQGPTAAKYILETVKPQRIAIIHDKQQYGEGLARSVQDG----------LKQGNANIVFFDGITA-GEKDFS----ALIARLQKENIDFVYYGGYYPEMGQIVR-QARANGLKTQFMGPEGVGNASLSNIAGG"
assert record.rounds[2].alignments[21].hsps[0].match=="DG   V  H   SS    S  Y D G  M +   T   + Q      +        QGP A      T +  ++ I  +   Y        +DG          +K  N     +  +T   E DF       + R+ +E+++F+      P    I R +  ++  + QF  PE  G   +      "
assert record.rounds[2].alignments[21].hsps[0].sbjct=="DGHMVVRSHARVSSLTLKSIQYRDAGEYMCTASNT---IGQDSQSIDLEFQYAPKLQGPVAVY----TWEGNQVNITCEVFAYPSATISWFRDGQLLPSSNYSNIKIYNTPSASYLEVTPDSENDFGNYNCTAVNRIGQESLEFILVQADTPSSPSIDRVEPYSSTAQVQFDEPEATGGVPILKYKAE"
assert record.rounds[2].alignments[21].hsps[0].query_start==91
assert record.rounds[2].alignments[21].hsps[0].query_end==262
assert record.rounds[2].alignments[21].hsps[0].sbjct_start==357
assert record.rounds[2].alignments[21].hsps[0].sbjct_end==537
assert record.rounds[2].alignments[22].hsps[0].query=="IARLQKENIDFVYYGGYYPEMGQ-----IVRQARANGLKTQFMGPEGVGNASLSNIAGGAAEGMLVTMPKRYDQDPANKAIVEALKADKKDPSGPYVWITYAAVQSLATAMTRSASHRPLDLVKDLKANGADTVIGPLKWDEKGDL"
assert record.rounds[2].alignments[22].hsps[0].match=="I +++K   + V + G     G+     + + A   G+ +    PE    + LS +A       +  M    + DP + A  +      ++  G    + Y  V +L        +     + K L+  GA  V      D+ G+L"
assert record.rounds[2].alignments[22].hsps[0].sbjct=="IEKMKKTGRNIVVFYGSQTGTGEEFANRLSKDAHRYGMGSMAADPEEYDMSELSRLAEIGNSLAIFCMATYGEGDPTDNA--QDFYDWLQETDGQLSGVNY-PVFALGDKTYEHYNAMGAYVDKRLEELGAKRVFDLGMGDDDGNL"
assert record.rounds[2].alignments[22].hsps[0].query_start==209
assert record.rounds[2].alignments[22].hsps[0].query_end==349
assert record.rounds[2].alignments[22].hsps[0].sbjct_start==15
assert record.rounds[2].alignments[22].hsps[0].sbjct_end==157
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.319, 0.118, 0.300]
assert record.ka_params_gap==[0.270, 0.0413, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==51436041
assert record.num_sequences==82258
assert record.num_extends==1847150
assert record.num_good_extends==5653
assert record.num_seqs_better_e==61
assert record.hsps_no_gap==21
assert record.hsps_prelim_gapped==40
assert record.hsps_gapped==65
assert record.query_length==369
assert record.database_length==29652561
assert record.effective_hsp_length==53
assert record.effective_query_length==316
assert record.effective_database_length==25292887
assert record.effective_search_space==7992552292
assert record.effective_search_space_used==7992552292
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.4)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.9)
assert record.blast_cutoff==(65,29.9)

handle = open('Blast/bt009')
record = pb_parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|399896|sp|Q02134|HIS7_LACLA IMIDAZOLEGLYCEROL-PHOSPHATE\nDEHYDRATASE (IGPD)"
assert record.query_letters==200
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.rounds)==2
assert len(record.rounds[0].new_seqs)==30
assert record.rounds[0].new_seqs[0].title=="gi|399896|sp|Q02134|HIS7_LACLA IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[0].score==409
assert record.rounds[0].new_seqs[0].e==1.0000000000000001e-114
assert record.rounds[0].new_seqs[1].title=="gi|462273|sp|P34047|HIS7_ARATH IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[1].score==198
assert record.rounds[0].new_seqs[1].e==6e-51
assert record.rounds[0].new_seqs[2].title=="gi|2495230|sp|Q43072|HIS7_PEA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYD..."
assert record.rounds[0].new_seqs[2].score==196
assert record.rounds[0].new_seqs[2].e==4e-50
assert record.rounds[0].new_seqs[3].title=="gi|123157|sp|P18787|HIS7_AZOBR IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[3].score==185
assert record.rounds[0].new_seqs[3].e==5.0000000000000001e-47
assert record.rounds[0].new_seqs[4].title=="gi|462275|sp|P34048|HIS7_WHEAT IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[4].score==181
assert record.rounds[0].new_seqs[4].e==record.rounds[0].new_seqs[4].e
assert record.rounds[0].new_seqs[5].title=="gi|123161|sp|P16247|HIS7_STRCO IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[5].score==178
assert record.rounds[0].new_seqs[5].e==7e-45
assert record.rounds[0].new_seqs[6].title=="gi|462272|sp|Q05068|HIS7_ANASP IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[6].score==178
assert record.rounds[0].new_seqs[6].e==7e-45
assert record.rounds[0].new_seqs[7].title=="gi|123158|sp|P06987|HIS7_ECOLI HISTIDINE BIOSYNTHESIS BIFUNCTIO..."
assert record.rounds[0].new_seqs[7].score==175
assert record.rounds[0].new_seqs[7].e==7.9999999999999996e-44
assert record.rounds[0].new_seqs[8].title=="gi|1346293|sp|P48054|HIS7_SYNY3 IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[8].score==174
assert record.rounds[0].new_seqs[8].e==1.0000000000000001e-43
assert record.rounds[0].new_seqs[9].title=="gi|1170286|sp|P44327|HIS7_HAEIN HISTIDINE BIOSYNTHESIS BIFUNCTI..."
assert record.rounds[0].new_seqs[9].score==168
assert record.rounds[0].new_seqs[9].e==8.0000000000000003e-42
assert record.rounds[0].new_seqs[10].title=="gi|2495224|sp|O06590|HIS7_MYCTU IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[10].score==167
assert record.rounds[0].new_seqs[10].e==2e-41
assert record.rounds[0].new_seqs[11].title=="gi|123160|sp|P10368|HIS7_SALTY HISTIDINE BIOSYNTHESIS BIFUNCTIO..."
assert record.rounds[0].new_seqs[11].score==166
assert record.rounds[0].new_seqs[11].e==2e-41
assert record.rounds[0].new_seqs[12].title=="gi|2495226|sp|Q50504|HIS7_METTH PROBABLE IMIDAZOLEGLYCEROL-PHOS..."
assert record.rounds[0].new_seqs[12].score==153
assert record.rounds[0].new_seqs[12].e==3e-37
assert record.rounds[0].new_seqs[13].title=="gi|729718|sp|P40919|HIS7_CRYNE IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[13].score==152
assert record.rounds[0].new_seqs[13].e==7.0000000000000003e-37
assert record.rounds[0].new_seqs[14].title=="gi|3334215|sp|O33773|HIS7_SULSO PROBABLE IMIDAZOLEGLYCEROL-PHOS..."
assert record.rounds[0].new_seqs[14].score==151
assert record.rounds[0].new_seqs[14].e==9.0000000000000008e-37
assert record.rounds[0].new_seqs[15].title=="gi|123159|sp|P28624|HIS7_PHYPR IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[15].score==149
assert record.rounds[0].new_seqs[15].e==3.0000000000000002e-36
assert record.rounds[0].new_seqs[16].title=="gi|729719|sp|P40374|HIS7_SCHPO IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[16].score==136
assert record.rounds[0].new_seqs[16].e==3e-32
assert record.rounds[0].new_seqs[17].title=="gi|2495227|sp|P56090|HIS7_CANAL IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[17].score==128
assert record.rounds[0].new_seqs[17].e==8.9999999999999993e-30
assert record.rounds[0].new_seqs[18].title=="gi|2495225|sp|Q58109|HIS7_METJA PROBABLE IMIDAZOLEGLYCEROL-PHOS..."
assert record.rounds[0].new_seqs[18].score==126
assert record.rounds[0].new_seqs[18].e==3.9999999999999998e-29
assert record.rounds[0].new_seqs[19].title=="gi|399897|sp|Q02986|HIS7_SACKL IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[19].score==125
assert record.rounds[0].new_seqs[19].e==6.0000000000000005e-29
assert record.rounds[0].new_seqs[20].title=="gi|2495229|sp|Q92447|HIS7_PICPA IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[20].score==125
assert record.rounds[0].new_seqs[20].e==6.0000000000000005e-29
assert record.rounds[0].new_seqs[21].title=="gi|2495228|sp|Q12578|HIS7_CANGA IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[21].score==123
assert record.rounds[0].new_seqs[21].e==1.9999999999999999e-28
assert record.rounds[0].new_seqs[22].title=="gi|2506514|sp|P06633|HIS7_YEAST IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[22].score==122
assert record.rounds[0].new_seqs[22].e==3.9999999999999999e-28
assert record.rounds[0].new_seqs[23].title=="gi|462274|sp|P34041|HIS7_TRIHA IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[23].score==106
assert record.rounds[0].new_seqs[23].e==3e-23
assert record.rounds[0].new_seqs[24].title=="gi|1345641|sp|P49264|C7B1_THLAR CYTOCHROME P450 71B1 (CYPLXXIB1)"
assert record.rounds[0].new_seqs[24].score==35
assert record.rounds[0].new_seqs[24].e==0.130000
assert record.rounds[0].new_seqs[25].title=="gi|1731346|sp|Q10698|YY29_MYCTU PROBABLE DIPEPTIDASE CY49.29C"
assert record.rounds[0].new_seqs[25].score==32
assert record.rounds[0].new_seqs[25].e==1.100000
assert record.rounds[0].new_seqs[26].title=="gi|3287839|sp|Q01812|GLK4_RAT GLUTAMATE RECEPTOR, IONOTROPIC KA..."
assert record.rounds[0].new_seqs[26].score==30
assert record.rounds[0].new_seqs[26].e==3.300000
assert record.rounds[0].new_seqs[27].title=="gi|3123025|sp|Q94637|VIT6_OSCBR VITELLOGENIN 6 PRECURSOR"
assert record.rounds[0].new_seqs[27].score==29
assert record.rounds[0].new_seqs[27].e==5.600000
assert record.rounds[0].new_seqs[28].title=="gi|3287848|sp|Q16099|GLK4_HUMAN GLUTAMATE RECEPTOR, IONOTROPIC ..."
assert record.rounds[0].new_seqs[28].score==29
assert record.rounds[0].new_seqs[28].e==9.700000
assert record.rounds[0].new_seqs[29].title=="gi|1174406|sp|P36126|SP14_YEAST PHOSPHOLIPASE D1 (PLD 1) (CHOLI..."
assert record.rounds[0].new_seqs[29].score==29
assert record.rounds[0].new_seqs[29].e==9.700000
assert len(record.rounds[0].alignments)==30
assert record.rounds[0].alignments[0].title==">gi|399896|sp|Q02134|HIS7_LACLA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[0].length==200
assert record.rounds[0].alignments[1].title==">gi|462273|sp|P34047|HIS7_ARATH IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[1].length==270
assert record.rounds[0].alignments[2].title==">gi|2495230|sp|Q43072|HIS7_PEA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[2].length==281
assert record.rounds[0].alignments[3].title==">gi|123157|sp|P18787|HIS7_AZOBR IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[3].length==207
assert record.rounds[0].alignments[4].title==">gi|462275|sp|P34048|HIS7_WHEAT IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[4].length==195
assert record.rounds[0].alignments[5].title==">gi|123161|sp|P16247|HIS7_STRCO IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[5].length==197
assert record.rounds[0].alignments[6].title==">gi|462272|sp|Q05068|HIS7_ANASP IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[6].length==209
assert record.rounds[0].alignments[7].title==">gi|123158|sp|P06987|HIS7_ECOLI HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[0].alignments[7].length==355
assert record.rounds[0].alignments[8].title==">gi|1346293|sp|P48054|HIS7_SYNY3 IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[8].length==210
assert record.rounds[0].alignments[9].title==">gi|1170286|sp|P44327|HIS7_HAEIN HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[0].alignments[9].length==362
assert record.rounds[0].alignments[10].title==">gi|2495224|sp|O06590|HIS7_MYCTU IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[10].length==210
assert record.rounds[0].alignments[11].title==">gi|123160|sp|P10368|HIS7_SALTY HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[0].alignments[11].length==354
assert record.rounds[0].alignments[12].title==">gi|2495226|sp|Q50504|HIS7_METTH PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[12].length==194
assert record.rounds[0].alignments[13].title==">gi|729718|sp|P40919|HIS7_CRYNE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[13].length==202
assert record.rounds[0].alignments[14].title==">gi|3334215|sp|O33773|HIS7_SULSO PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[14].length==193
assert record.rounds[0].alignments[15].title==">gi|123159|sp|P28624|HIS7_PHYPR IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[15].length==452
assert record.rounds[0].alignments[16].title==">gi|729719|sp|P40374|HIS7_SCHPO IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[16].length==216
assert record.rounds[0].alignments[17].title==">gi|2495227|sp|P56090|HIS7_CANAL IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[17].length==223
assert record.rounds[0].alignments[18].title==">gi|2495225|sp|Q58109|HIS7_METJA PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[18].length==197
assert record.rounds[0].alignments[19].title==">gi|399897|sp|Q02986|HIS7_SACKL IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[19].length==232
assert record.rounds[0].alignments[20].title==">gi|2495229|sp|Q92447|HIS7_PICPA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[20].length==224
assert record.rounds[0].alignments[21].title==">gi|2495228|sp|Q12578|HIS7_CANGA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[21].length==210
assert record.rounds[0].alignments[22].title==">gi|2506514|sp|P06633|HIS7_YEAST IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[22].length==220
assert record.rounds[0].alignments[23].title==">gi|462274|sp|P34041|HIS7_TRIHA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[23].length==208
assert record.rounds[0].alignments[24].title==">gi|1345641|sp|P49264|C7B1_THLAR CYTOCHROME P450 71B1 (CYPLXXIB1)"
assert record.rounds[0].alignments[24].length==496
assert record.rounds[0].alignments[25].title==">gi|1731346|sp|Q10698|YY29_MYCTU PROBABLE DIPEPTIDASE CY49.29C"
assert record.rounds[0].alignments[25].length==375
assert record.rounds[0].alignments[26].title==">gi|3287839|sp|Q01812|GLK4_RAT GLUTAMATE RECEPTOR, IONOTROPIC KAINATE 4 PRECURSOR (GLUTAMATE RECEPTOR KA-1) (KA1)"
assert record.rounds[0].alignments[26].length==956
assert record.rounds[0].alignments[27].title==">gi|3123025|sp|Q94637|VIT6_OSCBR VITELLOGENIN 6 PRECURSOR"
assert record.rounds[0].alignments[27].length==1660
assert record.rounds[0].alignments[28].title==">gi|3287848|sp|Q16099|GLK4_HUMAN GLUTAMATE RECEPTOR, IONOTROPIC KAINATE 4 PRECURSOR (GLUTAMATE RECEPTOR KA-1) (KA1) (EXCITATORY AMINO ACID RECEPTOR 1) (EAA1)"
assert record.rounds[0].alignments[28].length==956
assert record.rounds[0].alignments[29].title==">gi|1174406|sp|P36126|SP14_YEAST PHOSPHOLIPASE D1 (PLD 1) (CHOLINE PHOSPHATASE 1) (PHOSPHATIDYLCHOLINE-HYDROLYZING PHOSPHOLIPASE D1) (MEIOSIS-SPECIFIC SPORULATION PROTEIN SPO14)"
assert record.rounds[0].alignments[29].length==1380
assert len(record.rounds[1].new_seqs)==0
assert len(record.rounds[1].alignments)==24
assert record.rounds[1].alignments[0].title==">gi|2495230|sp|Q43072|HIS7_PEA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[0].length==281
assert record.rounds[1].alignments[1].title==">gi|462273|sp|P34047|HIS7_ARATH IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[1].length==270
assert record.rounds[1].alignments[2].title==">gi|399896|sp|Q02134|HIS7_LACLA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[2].length==200
assert record.rounds[1].alignments[3].title==">gi|1346293|sp|P48054|HIS7_SYNY3 IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[3].length==210
assert record.rounds[1].alignments[4].title==">gi|462272|sp|Q05068|HIS7_ANASP IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[4].length==209
assert record.rounds[1].alignments[5].title==">gi|462275|sp|P34048|HIS7_WHEAT IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[5].length==195
assert record.rounds[1].alignments[6].title==">gi|123161|sp|P16247|HIS7_STRCO IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[6].length==197
assert record.rounds[1].alignments[7].title==">gi|2506514|sp|P06633|HIS7_YEAST IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[7].length==220
assert record.rounds[1].alignments[8].title==">gi|2495227|sp|P56090|HIS7_CANAL IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[8].length==223
assert record.rounds[1].alignments[9].title==">gi|399897|sp|Q02986|HIS7_SACKL IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[9].length==232
assert record.rounds[1].alignments[10].title==">gi|2495228|sp|Q12578|HIS7_CANGA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[10].length==210
assert record.rounds[1].alignments[11].title==">gi|123158|sp|P06987|HIS7_ECOLI HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[1].alignments[11].length==355
assert record.rounds[1].alignments[12].title==">gi|123157|sp|P18787|HIS7_AZOBR IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[12].length==207
assert record.rounds[1].alignments[13].title==">gi|729718|sp|P40919|HIS7_CRYNE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[13].length==202
assert record.rounds[1].alignments[14].title==">gi|2495229|sp|Q92447|HIS7_PICPA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[14].length==224
assert record.rounds[1].alignments[15].title==">gi|1170286|sp|P44327|HIS7_HAEIN HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[1].alignments[15].length==362
assert record.rounds[1].alignments[16].title==">gi|729719|sp|P40374|HIS7_SCHPO IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[16].length==216
assert record.rounds[1].alignments[17].title==">gi|123159|sp|P28624|HIS7_PHYPR IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[17].length==452
assert record.rounds[1].alignments[18].title==">gi|123160|sp|P10368|HIS7_SALTY HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[1].alignments[18].length==354
assert record.rounds[1].alignments[19].title==">gi|2495224|sp|O06590|HIS7_MYCTU IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[19].length==210
assert record.rounds[1].alignments[20].title==">gi|2495226|sp|Q50504|HIS7_METTH PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[20].length==194
assert record.rounds[1].alignments[21].title==">gi|2495225|sp|Q58109|HIS7_METJA PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[21].length==197
assert record.rounds[1].alignments[22].title==">gi|462274|sp|P34041|HIS7_TRIHA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[22].length==208
assert record.rounds[1].alignments[23].title==">gi|3334215|sp|O33773|HIS7_SULSO PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[23].length==193
assert record.rounds[0].alignments[0].hsps[0].score==1040
assert record.rounds[0].alignments[0].hsps[0].bits==409
assert record.rounds[0].alignments[0].hsps[0].expect==1e-114
assert len(record.rounds[0].alignments[0].hsps)==1
assert record.rounds[0].alignments[1].hsps[0].score==499
assert record.rounds[0].alignments[1].hsps[0].bits==198
assert record.rounds[0].alignments[1].hsps[0].expect==6e-51
assert len(record.rounds[0].alignments[1].hsps)==1
assert record.rounds[0].alignments[2].hsps[0].score==492
assert record.rounds[0].alignments[2].hsps[0].bits==196
assert record.rounds[0].alignments[2].hsps[0].expect==4e-50
assert len(record.rounds[0].alignments[2].hsps)==1
assert record.rounds[0].alignments[3].hsps[0].score==465
assert record.rounds[0].alignments[3].hsps[0].bits==185
assert record.rounds[0].alignments[3].hsps[0].expect==5e-47
assert len(record.rounds[0].alignments[3].hsps)==1
assert record.rounds[0].alignments[4].hsps[0].score==455
assert record.rounds[0].alignments[4].hsps[0].bits==181
assert record.rounds[0].alignments[4].hsps[0].expect==8e-46
assert len(record.rounds[0].alignments[4].hsps)==1
assert record.rounds[0].alignments[5].hsps[0].score==447
assert record.rounds[0].alignments[5].hsps[0].bits==178
assert record.rounds[0].alignments[5].hsps[0].expect==7e-45
assert len(record.rounds[0].alignments[5].hsps)==1
assert record.rounds[0].alignments[6].hsps[0].score==447
assert record.rounds[0].alignments[6].hsps[0].bits==178
assert record.rounds[0].alignments[6].hsps[0].expect==7e-45
assert len(record.rounds[0].alignments[6].hsps)==1
assert record.rounds[0].alignments[7].hsps[0].score==438
assert record.rounds[0].alignments[7].hsps[0].bits==175
assert record.rounds[0].alignments[7].hsps[0].expect==8e-44
assert len(record.rounds[0].alignments[7].hsps)==1
assert record.rounds[0].alignments[8].hsps[0].score==437
assert record.rounds[0].alignments[8].hsps[0].bits==174
assert record.rounds[0].alignments[8].hsps[0].expect==1e-43
assert len(record.rounds[0].alignments[8].hsps)==1
assert record.rounds[0].alignments[9].hsps[0].score==421
assert record.rounds[0].alignments[9].hsps[0].bits==168
assert record.rounds[0].alignments[9].hsps[0].expect==8e-42
assert len(record.rounds[0].alignments[9].hsps)==1
assert record.rounds[0].alignments[10].hsps[0].score==418
assert record.rounds[0].alignments[10].hsps[0].bits==167
assert record.rounds[0].alignments[10].hsps[0].expect==2e-41
assert len(record.rounds[0].alignments[10].hsps)==1
assert record.rounds[0].alignments[11].hsps[0].score==417
assert record.rounds[0].alignments[11].hsps[0].bits==166
assert record.rounds[0].alignments[11].hsps[0].expect==2e-41
assert len(record.rounds[0].alignments[11].hsps)==1
assert record.rounds[0].alignments[12].hsps[0].score==382
assert record.rounds[0].alignments[12].hsps[0].bits==153
assert record.rounds[0].alignments[12].hsps[0].expect==3e-37
assert len(record.rounds[0].alignments[12].hsps)==1
assert record.rounds[0].alignments[13].hsps[0].score==379
assert record.rounds[0].alignments[13].hsps[0].bits==152
assert record.rounds[0].alignments[13].hsps[0].expect==7e-37
assert len(record.rounds[0].alignments[13].hsps)==1
assert record.rounds[0].alignments[14].hsps[0].score==378
assert record.rounds[0].alignments[14].hsps[0].bits==151
assert record.rounds[0].alignments[14].hsps[0].expect==9e-37
assert len(record.rounds[0].alignments[14].hsps)==1
assert record.rounds[0].alignments[15].hsps[0].score==373
assert record.rounds[0].alignments[15].hsps[0].bits==149
assert record.rounds[0].alignments[15].hsps[0].expect==3e-36
assert len(record.rounds[0].alignments[15].hsps)==1
assert record.rounds[0].alignments[16].hsps[0].score==339
assert record.rounds[0].alignments[16].hsps[0].bits==136
assert record.rounds[0].alignments[16].hsps[0].expect==3e-32
assert len(record.rounds[0].alignments[16].hsps)==1
assert record.rounds[0].alignments[17].hsps[0].score==318
assert record.rounds[0].alignments[17].hsps[0].bits==128
assert record.rounds[0].alignments[17].hsps[0].expect==9e-30
assert len(record.rounds[0].alignments[17].hsps)==1
assert record.rounds[0].alignments[18].hsps[0].score==313
assert record.rounds[0].alignments[18].hsps[0].bits==126
assert record.rounds[0].alignments[18].hsps[0].expect==4e-29
assert len(record.rounds[0].alignments[18].hsps)==1
assert record.rounds[0].alignments[19].hsps[0].score==311
assert record.rounds[0].alignments[19].hsps[0].bits==125
assert record.rounds[0].alignments[19].hsps[0].expect==6e-29
assert len(record.rounds[0].alignments[19].hsps)==1
assert record.rounds[0].alignments[20].hsps[0].score==311
assert record.rounds[0].alignments[20].hsps[0].bits==125
assert record.rounds[0].alignments[20].hsps[0].expect==6e-29
assert len(record.rounds[0].alignments[20].hsps)==1
assert record.rounds[0].alignments[21].hsps[0].score==306
assert record.rounds[0].alignments[21].hsps[0].bits==123
assert record.rounds[0].alignments[21].hsps[0].expect==2e-28
assert len(record.rounds[0].alignments[21].hsps)==1
assert record.rounds[0].alignments[22].hsps[0].score==304
assert record.rounds[0].alignments[22].hsps[0].bits==122
assert record.rounds[0].alignments[22].hsps[0].expect==4e-28
assert len(record.rounds[0].alignments[22].hsps)==1
assert record.rounds[0].alignments[23].hsps[0].score==263
assert record.rounds[0].alignments[23].hsps[0].bits==106
assert record.rounds[0].alignments[23].hsps[0].expect==3e-23
assert len(record.rounds[0].alignments[23].hsps)==1
assert record.rounds[0].alignments[24].hsps[0].score==78
assert record.rounds[0].alignments[24].hsps[0].bits==34.8
assert record.rounds[0].alignments[24].hsps[0].expect==0.13
assert len(record.rounds[0].alignments[24].hsps)==1
assert record.rounds[0].alignments[25].hsps[0].score==70
assert record.rounds[0].alignments[25].hsps[0].bits==31.7
assert record.rounds[0].alignments[25].hsps[0].expect==1.1
assert len(record.rounds[0].alignments[25].hsps)==1
assert record.rounds[0].alignments[26].hsps[0].score==66
assert record.rounds[0].alignments[26].hsps[0].bits==30.1
assert record.rounds[0].alignments[26].hsps[0].expect==3.3
assert len(record.rounds[0].alignments[26].hsps)==1
assert record.rounds[0].alignments[27].hsps[0].score==64
assert record.rounds[0].alignments[27].hsps[0].bits==29.3
assert record.rounds[0].alignments[27].hsps[0].expect==5.6
assert len(record.rounds[0].alignments[27].hsps)==1
assert record.rounds[0].alignments[28].hsps[0].score==62
assert record.rounds[0].alignments[28].hsps[0].bits==28.6
assert record.rounds[0].alignments[28].hsps[0].expect==9.7
assert len(record.rounds[0].alignments[28].hsps)==1
assert record.rounds[0].alignments[29].hsps[0].score==62
assert record.rounds[0].alignments[29].hsps[0].bits==28.6
assert record.rounds[0].alignments[29].hsps[0].expect==9.7
assert record.rounds[1].alignments[0].hsps[0].score==820
assert record.rounds[1].alignments[0].hsps[0].bits==323
assert record.rounds[1].alignments[0].hsps[0].expect==1e-88
assert len(record.rounds[1].alignments[0].hsps)==1
assert record.rounds[1].alignments[1].hsps[0].score==817
assert record.rounds[1].alignments[1].hsps[0].bits==322
assert record.rounds[1].alignments[1].hsps[0].expect==3e-88
assert len(record.rounds[1].alignments[1].hsps)==1
assert record.rounds[1].alignments[2].hsps[0].score==808
assert record.rounds[1].alignments[2].hsps[0].bits==318
assert record.rounds[1].alignments[2].hsps[0].expect==4e-87
assert len(record.rounds[1].alignments[2].hsps)==1
assert record.rounds[1].alignments[3].hsps[0].score==798
assert record.rounds[1].alignments[3].hsps[0].bits==315
assert record.rounds[1].alignments[3].hsps[0].expect==5e-86
assert len(record.rounds[1].alignments[3].hsps)==1
assert record.rounds[1].alignments[4].hsps[0].score==795
assert record.rounds[1].alignments[4].hsps[0].bits==313
assert record.rounds[1].alignments[4].hsps[0].expect==1e-85
assert len(record.rounds[1].alignments[4].hsps)==1
assert record.rounds[1].alignments[5].hsps[0].score==793
assert record.rounds[1].alignments[5].hsps[0].bits==313
assert record.rounds[1].alignments[5].hsps[0].expect==2e-85
assert len(record.rounds[1].alignments[5].hsps)==1
assert record.rounds[1].alignments[6].hsps[0].score==776
assert record.rounds[1].alignments[6].hsps[0].bits==306
assert record.rounds[1].alignments[6].hsps[0].expect==2e-83
assert len(record.rounds[1].alignments[6].hsps)==1
assert record.rounds[1].alignments[7].hsps[0].score==772
assert record.rounds[1].alignments[7].hsps[0].bits==304
assert record.rounds[1].alignments[7].hsps[0].expect==6e-83
assert len(record.rounds[1].alignments[7].hsps)==1
assert record.rounds[1].alignments[8].hsps[0].score==771
assert record.rounds[1].alignments[8].hsps[0].bits==304
assert record.rounds[1].alignments[8].hsps[0].expect==8e-83
assert len(record.rounds[1].alignments[8].hsps)==1
assert record.rounds[1].alignments[9].hsps[0].score==770
assert record.rounds[1].alignments[9].hsps[0].bits==304
assert record.rounds[1].alignments[9].hsps[0].expect==1e-82
assert len(record.rounds[1].alignments[9].hsps)==1
assert record.rounds[1].alignments[10].hsps[0].score==767
assert record.rounds[1].alignments[10].hsps[0].bits==303
assert record.rounds[1].alignments[10].hsps[0].expect==2e-82
assert len(record.rounds[1].alignments[10].hsps)==1
assert record.rounds[1].alignments[11].hsps[0].score==765
assert record.rounds[1].alignments[11].hsps[0].bits==302
assert record.rounds[1].alignments[11].hsps[0].expect==4e-82
assert len(record.rounds[1].alignments[11].hsps)==1
assert record.rounds[1].alignments[12].hsps[0].score==762
assert record.rounds[1].alignments[12].hsps[0].bits==301
assert record.rounds[1].alignments[12].hsps[0].expect==9e-82
assert len(record.rounds[1].alignments[12].hsps)==1
assert record.rounds[1].alignments[13].hsps[0].score==759
assert record.rounds[1].alignments[13].hsps[0].bits==299
assert record.rounds[1].alignments[13].hsps[0].expect==2e-81
assert len(record.rounds[1].alignments[13].hsps)==1
assert record.rounds[1].alignments[14].hsps[0].score==756
assert record.rounds[1].alignments[14].hsps[0].bits==298
assert record.rounds[1].alignments[14].hsps[0].expect==5e-81
assert len(record.rounds[1].alignments[14].hsps)==1
assert record.rounds[1].alignments[15].hsps[0].score==741
assert record.rounds[1].alignments[15].hsps[0].bits==292
assert record.rounds[1].alignments[15].hsps[0].expect==3e-79
assert len(record.rounds[1].alignments[15].hsps)==1
assert record.rounds[1].alignments[16].hsps[0].score==734
assert record.rounds[1].alignments[16].hsps[0].bits==290
assert record.rounds[1].alignments[16].hsps[0].expect==2e-78
assert len(record.rounds[1].alignments[16].hsps)==1
assert record.rounds[1].alignments[17].hsps[0].score==734
assert record.rounds[1].alignments[17].hsps[0].bits==290
assert record.rounds[1].alignments[17].hsps[0].expect==2e-78
assert len(record.rounds[1].alignments[17].hsps)==1
assert record.rounds[1].alignments[18].hsps[0].score==726
assert record.rounds[1].alignments[18].hsps[0].bits==287
assert record.rounds[1].alignments[18].hsps[0].expect==1e-77
assert len(record.rounds[1].alignments[18].hsps)==1
assert record.rounds[1].alignments[19].hsps[0].score==716
assert record.rounds[1].alignments[19].hsps[0].bits==283
assert record.rounds[1].alignments[19].hsps[0].expect==2e-76
assert len(record.rounds[1].alignments[19].hsps)==1
assert record.rounds[1].alignments[20].hsps[0].score==695
assert record.rounds[1].alignments[20].hsps[0].bits==274
assert record.rounds[1].alignments[20].hsps[0].expect==6e-74
assert len(record.rounds[1].alignments[20].hsps)==1
assert record.rounds[1].alignments[21].hsps[0].score==685
assert record.rounds[1].alignments[21].hsps[0].bits==271
assert record.rounds[1].alignments[21].hsps[0].expect==1e-72
assert len(record.rounds[1].alignments[21].hsps)==1
assert record.rounds[1].alignments[22].hsps[0].score==680
assert record.rounds[1].alignments[22].hsps[0].bits==269
assert record.rounds[1].alignments[22].hsps[0].expect==4e-72
assert len(record.rounds[1].alignments[22].hsps)==1
assert record.rounds[1].alignments[23].hsps[0].score==662
assert record.rounds[1].alignments[23].hsps[0].bits==262
assert record.rounds[1].alignments[23].hsps[0].expect==5e-70
assert len(record.rounds[1].alignments[23].hsps)==1
assert record.rounds[0].alignments[0].hsps[0].identities==(200, 200)
assert record.rounds[0].alignments[0].hsps[0].positives==(200, 200)
assert record.rounds[0].alignments[1].hsps[0].identities==(99, 198)
assert record.rounds[0].alignments[1].hsps[0].positives==(135, 198)
assert record.rounds[0].alignments[1].hsps[0].gaps==(4, 198)
assert record.rounds[0].alignments[2].hsps[0].identities==(96, 199)
assert record.rounds[0].alignments[2].hsps[0].positives==(136, 199)
assert record.rounds[0].alignments[2].hsps[0].gaps==(4, 199)
assert record.rounds[0].alignments[3].hsps[0].identities==(91, 194)
assert record.rounds[0].alignments[3].hsps[0].positives==(126, 194)
assert record.rounds[0].alignments[3].hsps[0].gaps==(4, 194)
assert record.rounds[0].alignments[4].hsps[0].identities==(93, 194)
assert record.rounds[0].alignments[4].hsps[0].positives==(128, 194)
assert record.rounds[0].alignments[4].hsps[0].gaps==(4, 194)
assert record.rounds[0].alignments[5].hsps[0].identities==(89, 200)
assert record.rounds[0].alignments[5].hsps[0].positives==(124, 200)
assert record.rounds[0].alignments[5].hsps[0].gaps==(3, 200)
assert record.rounds[0].alignments[6].hsps[0].identities==(91, 198)
assert record.rounds[0].alignments[6].hsps[0].positives==(131, 198)
assert record.rounds[0].alignments[6].hsps[0].gaps==(4, 198)
assert record.rounds[0].alignments[7].hsps[0].identities==(91, 198)
assert record.rounds[0].alignments[7].hsps[0].positives==(130, 198)
assert record.rounds[0].alignments[7].hsps[0].gaps==(9, 198)
assert record.rounds[0].alignments[8].hsps[0].identities==(88, 198)
assert record.rounds[0].alignments[8].hsps[0].positives==(129, 198)
assert record.rounds[0].alignments[8].hsps[0].gaps==(4, 198)
assert record.rounds[0].alignments[9].hsps[0].identities==(89, 198)
assert record.rounds[0].alignments[9].hsps[0].positives==(127, 198)
assert record.rounds[0].alignments[9].hsps[0].gaps==(9, 198)
assert record.rounds[0].alignments[10].hsps[0].identities==(92, 207)
assert record.rounds[0].alignments[10].hsps[0].positives==(125, 207)
assert record.rounds[0].alignments[10].hsps[0].gaps==(14, 207)
assert record.rounds[0].alignments[11].hsps[0].identities==(89, 198)
assert record.rounds[0].alignments[11].hsps[0].positives==(129, 198)
assert record.rounds[0].alignments[11].hsps[0].gaps==(10, 198)
assert record.rounds[0].alignments[12].hsps[0].identities==(81, 198)
assert record.rounds[0].alignments[12].hsps[0].positives==(122, 198)
assert record.rounds[0].alignments[12].hsps[0].gaps==(8, 198)
assert record.rounds[0].alignments[13].hsps[0].identities==(83, 203)
assert record.rounds[0].alignments[13].hsps[0].positives==(120, 203)
assert record.rounds[0].alignments[13].hsps[0].gaps==(11, 203)
assert record.rounds[0].alignments[14].hsps[0].identities==(88, 201)
assert record.rounds[0].alignments[14].hsps[0].positives==(128, 201)
assert record.rounds[0].alignments[14].hsps[0].gaps==(9, 201)
assert record.rounds[0].alignments[15].hsps[0].identities==(86, 198)
assert record.rounds[0].alignments[15].hsps[0].positives==(120, 198)
assert record.rounds[0].alignments[15].hsps[0].gaps==(6, 198)
assert record.rounds[0].alignments[16].hsps[0].identities==(84, 221)
assert record.rounds[0].alignments[16].hsps[0].positives==(114, 221)
assert record.rounds[0].alignments[16].hsps[0].gaps==(29, 221)
assert record.rounds[0].alignments[17].hsps[0].identities==(81, 227)
assert record.rounds[0].alignments[17].hsps[0].positives==(119, 227)
assert record.rounds[0].alignments[17].hsps[0].gaps==(33, 227)
assert record.rounds[0].alignments[18].hsps[0].identities==(80, 196)
assert record.rounds[0].alignments[18].hsps[0].positives==(107, 196)
assert record.rounds[0].alignments[18].hsps[0].gaps==(9, 196)
assert record.rounds[0].alignments[19].hsps[0].identities==(81, 222)
assert record.rounds[0].alignments[19].hsps[0].positives==(119, 222)
assert record.rounds[0].alignments[19].hsps[0].gaps==(30, 222)
assert record.rounds[0].alignments[20].hsps[0].identities==(84, 223)
assert record.rounds[0].alignments[20].hsps[0].positives==(116, 223)
assert record.rounds[0].alignments[20].hsps[0].gaps==(31, 223)
assert record.rounds[0].alignments[21].hsps[0].identities==(78, 215)
assert record.rounds[0].alignments[21].hsps[0].positives==(116, 215)
assert record.rounds[0].alignments[21].hsps[0].gaps==(24, 215)
assert record.rounds[0].alignments[22].hsps[0].identities==(79, 218)
assert record.rounds[0].alignments[22].hsps[0].positives==(114, 218)
assert record.rounds[0].alignments[22].hsps[0].gaps==(30, 218)
assert record.rounds[0].alignments[23].hsps[0].identities==(68, 202)
assert record.rounds[0].alignments[23].hsps[0].positives==(102, 202)
assert record.rounds[0].alignments[23].hsps[0].gaps==(28, 202)
assert record.rounds[0].alignments[24].hsps[0].identities==(34, 134)
assert record.rounds[0].alignments[24].hsps[0].positives==(60, 134)
assert record.rounds[0].alignments[24].hsps[0].gaps==(11, 134)
assert record.rounds[0].alignments[25].hsps[0].identities==(16, 45)
assert record.rounds[0].alignments[25].hsps[0].positives==(21, 45)
assert record.rounds[0].alignments[25].hsps[0].gaps==(3, 45)
assert record.rounds[0].alignments[26].hsps[0].identities==(17, 48)
assert record.rounds[0].alignments[26].hsps[0].positives==(24, 48)
assert record.rounds[0].alignments[26].hsps[0].gaps==(3, 48)
assert record.rounds[0].alignments[27].hsps[0].identities==(25, 70)
assert record.rounds[0].alignments[27].hsps[0].positives==(32, 70)
assert record.rounds[0].alignments[27].hsps[0].gaps==(5, 70)
assert record.rounds[0].alignments[28].hsps[0].identities==(16, 48)
assert record.rounds[0].alignments[28].hsps[0].positives==(24, 48)
assert record.rounds[0].alignments[28].hsps[0].gaps==(3, 48)
assert record.rounds[0].alignments[29].hsps[0].identities==(20, 65)
assert record.rounds[0].alignments[29].hsps[0].positives==(31, 65)
assert record.rounds[0].alignments[29].hsps[0].gaps==(7, 65)
assert record.rounds[1].alignments[0].hsps[0].identities==(96, 199)
assert record.rounds[1].alignments[0].hsps[0].positives==(136, 199)
assert record.rounds[1].alignments[0].hsps[0].gaps==(4, 199)
assert record.rounds[1].alignments[1].hsps[0].identities==(99, 198)
assert record.rounds[1].alignments[1].hsps[0].positives==(135, 198)
assert record.rounds[1].alignments[1].hsps[0].gaps==(4, 198)
assert record.rounds[1].alignments[2].hsps[0].identities==(200, 200)
assert record.rounds[1].alignments[2].hsps[0].positives==(200, 200)
assert record.rounds[1].alignments[3].hsps[0].identities==(88, 198)
assert record.rounds[1].alignments[3].hsps[0].positives==(129, 198)
assert record.rounds[1].alignments[3].hsps[0].gaps==(4, 198)
assert record.rounds[1].alignments[4].hsps[0].identities==(91, 198)
assert record.rounds[1].alignments[4].hsps[0].positives==(131, 198)
assert record.rounds[1].alignments[4].hsps[0].gaps==(4, 198)
assert record.rounds[1].alignments[5].hsps[0].identities==(93, 196)
assert record.rounds[1].alignments[5].hsps[0].positives==(128, 196)
assert record.rounds[1].alignments[5].hsps[0].gaps==(4, 196)
assert record.rounds[1].alignments[6].hsps[0].identities==(89, 200)
assert record.rounds[1].alignments[6].hsps[0].positives==(124, 200)
assert record.rounds[1].alignments[6].hsps[0].gaps==(3, 200)
assert record.rounds[1].alignments[7].hsps[0].identities==(78, 220)
assert record.rounds[1].alignments[7].hsps[0].positives==(115, 220)
assert record.rounds[1].alignments[7].hsps[0].gaps==(30, 220)
assert record.rounds[1].alignments[8].hsps[0].identities==(81, 227)
assert record.rounds[1].alignments[8].hsps[0].positives==(119, 227)
assert record.rounds[1].alignments[8].hsps[0].gaps==(33, 227)
assert record.rounds[1].alignments[9].hsps[0].identities==(79, 222)
assert record.rounds[1].alignments[9].hsps[0].positives==(118, 222)
assert record.rounds[1].alignments[9].hsps[0].gaps==(30, 222)
assert record.rounds[1].alignments[10].hsps[0].identities==(76, 214)
assert record.rounds[1].alignments[10].hsps[0].positives==(113, 214)
assert record.rounds[1].alignments[10].hsps[0].gaps==(24, 214)
assert record.rounds[1].alignments[11].hsps[0].identities==(91, 198)
assert record.rounds[1].alignments[11].hsps[0].positives==(130, 198)
assert record.rounds[1].alignments[11].hsps[0].gaps==(9, 198)
assert record.rounds[1].alignments[12].hsps[0].identities==(91, 196)
assert record.rounds[1].alignments[12].hsps[0].positives==(127, 196)
assert record.rounds[1].alignments[12].hsps[0].gaps==(4, 196)
assert record.rounds[1].alignments[13].hsps[0].identities==(83, 203)
assert record.rounds[1].alignments[13].hsps[0].positives==(120, 203)
assert record.rounds[1].alignments[13].hsps[0].gaps==(11, 203)
assert record.rounds[1].alignments[14].hsps[0].identities==(82, 223)
assert record.rounds[1].alignments[14].hsps[0].positives==(115, 223)
assert record.rounds[1].alignments[14].hsps[0].gaps==(31, 223)
assert record.rounds[1].alignments[15].hsps[0].identities==(89, 198)
assert record.rounds[1].alignments[15].hsps[0].positives==(127, 198)
assert record.rounds[1].alignments[15].hsps[0].gaps==(9, 198)
assert record.rounds[1].alignments[16].hsps[0].identities==(83, 221)
assert record.rounds[1].alignments[16].hsps[0].positives==(114, 221)
assert record.rounds[1].alignments[16].hsps[0].gaps==(29, 221)
assert record.rounds[1].alignments[17].hsps[0].identities==(86, 199)
assert record.rounds[1].alignments[17].hsps[0].positives==(121, 199)
assert record.rounds[1].alignments[17].hsps[0].gaps==(8, 199)
assert record.rounds[1].alignments[18].hsps[0].identities==(89, 198)
assert record.rounds[1].alignments[18].hsps[0].positives==(129, 198)
assert record.rounds[1].alignments[18].hsps[0].gaps==(10, 198)
assert record.rounds[1].alignments[19].hsps[0].identities==(92, 207)
assert record.rounds[1].alignments[19].hsps[0].positives==(124, 207)
assert record.rounds[1].alignments[19].hsps[0].gaps==(14, 207)
assert record.rounds[1].alignments[20].hsps[0].identities==(81, 198)
assert record.rounds[1].alignments[20].hsps[0].positives==(122, 198)
assert record.rounds[1].alignments[20].hsps[0].gaps==(8, 198)
assert record.rounds[1].alignments[21].hsps[0].identities==(79, 196)
assert record.rounds[1].alignments[21].hsps[0].positives==(106, 196)
assert record.rounds[1].alignments[21].hsps[0].gaps==(9, 196)
assert record.rounds[1].alignments[22].hsps[0].identities==(68, 202)
assert record.rounds[1].alignments[22].hsps[0].positives==(102, 202)
assert record.rounds[1].alignments[22].hsps[0].gaps==(28, 202)
assert record.rounds[1].alignments[23].hsps[0].identities==(83, 200)
assert record.rounds[1].alignments[23].hsps[0].positives==(124, 200)
assert record.rounds[1].alignments[23].hsps[0].gaps==(7, 200)
assert record.rounds[0].alignments[0].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[0].hsps[0].match=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[0].hsps[0].sbjct=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[0].hsps[0].query_start==1
assert record.rounds[0].alignments[0].hsps[0].query_end==200
assert record.rounds[0].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[0].hsps[0].sbjct_end==200
assert record.rounds[0].alignments[1].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[1].hsps[0].match=="RI  + R TKET + + INLDGTG AD S+GI FLDHML  L  H  FD+ +   GD   V +D HH  ED+A+A+G  + + LG + GI R+G FT P+DEAL+   LD+SGRPYL ++ ++   Q++G YDT++ E FF++L   +G+TLH+ +  G+N+HHIIE  FK+ ARAL+QA   D  + G IPSSKGVL"
assert record.rounds[0].alignments[1].hsps[0].sbjct=="RIGEVKRVTKETNVSVKINLDGTGVADSSSGIPFLDHMLDQLASHGLFDVHVRATGD---VHIDDHHTNEDIALAIGTALLKALGERKGINRFGDFTAPLDEALIHVSLDLSGRPYLGYNLEIP-TQRVGTYDTQLVEHFFQSLVNTSGMTLHIRQLAGENSHHIIEATFKAFARALRQATETDPRRGGTIPSSKGVL"
assert record.rounds[0].alignments[1].hsps[0].query_start==3
assert record.rounds[0].alignments[1].hsps[0].query_end==200
assert record.rounds[0].alignments[1].hsps[0].sbjct_start==74
assert record.rounds[0].alignments[1].hsps[0].sbjct_end==267
assert record.rounds[0].alignments[2].hsps[0].query=="TRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[2].hsps[0].match=="TR+  + R TKET + + INLDG+G AD STGI FLDHML  L  H  FD+ +   GD   V +D HH  EDVA+A+G  + + LG++ GI R+G F+ P+DEAL+   LD+SGRP+L ++ D+   Q++G YDT++ E F +++   +G+TLH+ +  G+N+HHIIE  FK+ ARAL+QA   D  + G +PSSKGVL"
assert record.rounds[0].alignments[2].hsps[0].sbjct=="TRVGEVKRVTKETNVSVKINLDGSGVADSSTGIPFLDHMLDQLASHGLFDVHVKATGD---VHIDDHHTNEDVALAIGTALLQALGDRKGINRFGDFSAPLDEALIHVSLDLSGRPHLSYNLDIP-TQRVGTYDTQVVEHFLQSIVNTSGMTLHIRQLAGRNSHHIIEATFKAFARALRQATEYDPRRRGSVPSSKGVL"
assert record.rounds[0].alignments[2].hsps[0].query_start==2
assert record.rounds[0].alignments[2].hsps[0].query_end==200
assert record.rounds[0].alignments[2].hsps[0].sbjct_start==84
assert record.rounds[0].alignments[2].hsps[0].sbjct_end==278
assert record.rounds[0].alignments[3].hsps[0].query=="ITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[3].hsps[0].match=="I RNT ET+I +++NLDGTG  D+ TG+GFLDHML  L+ HS  DL +   GD   V +D HH  E   IA+G+ +++ +G++ GI+RYG   +PMDE L    LD S RPYL++    S   K+G  DTE+  E+F+A A  AG+TLH+   YG+N HHI+E  +K+ ARAL+  + ID  K   +PS+KG L"
assert record.rounds[0].alignments[3].hsps[0].sbjct=="IERNTTETRIRVAVNLDGTGVYDVKTGVGFLDHMLEQLSRHSLMDLSVAAEGD---VHIDAHHTTEHSGIAIGQAVAKAVGDRKGIQRYGHAYVPMDETLTRVALDFSNRPYLIWKVSFS-RDKIGDMDTELFREWFQAFAMAAGVTLHVECLYGENNHHIVESCYKALARALRAGIEIDPRKRDAVPSTKGTL"
assert record.rounds[0].alignments[3].hsps[0].query_start==7
assert record.rounds[0].alignments[3].hsps[0].query_end==200
assert record.rounds[0].alignments[3].hsps[0].sbjct_start==14
assert record.rounds[0].alignments[3].hsps[0].sbjct_end==203
assert record.rounds[0].alignments[4].hsps[0].query=="ITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[4].hsps[0].match=="+ R TKET + + INLDGTG A+ STGI FLDHML  L  H  FD+ +   GD     +D HH  ED+A+A+G  + + LG++ GI R+G FT P+DEA V   LD+SGRP+L     +   +++G YDT++ E FF++L   +G+TLH+ +  G N+HHIIE  FK+ ARAL+QA   D  + G +PSSKGVL"
assert record.rounds[0].alignments[4].hsps[0].sbjct=="VKRVTKETNVHVKINLDGTGVANSSTGIPFLDHMLDQLASHGLFDVYVKATGDTH---IDDHHSNEDIALAIGTALLQALGDRKGINRFGHFTAPLDEAAVEVILDLSGRPHLSCGLSIP-TERVGTYDTQLVEHFFQSLVNTSGMTLHIRQLAGNNSHHIIEATFKAFARALRQATEYDLRRQGTMPSSKGVL"
assert record.rounds[0].alignments[4].hsps[0].query_start==7
assert record.rounds[0].alignments[4].hsps[0].query_end==200
assert record.rounds[0].alignments[4].hsps[0].sbjct_start==3
assert record.rounds[0].alignments[4].hsps[0].sbjct_end==192
assert record.rounds[0].alignments[5].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[5].hsps[0].match=="M+R+  + R TKET + + I+LDGTG+ DI+TG+GF DHML  L  H  FDL +   GD   + +D HH IED A+ALG    + LG+K+GI R+G+ T+P+DE+L    +D+SGRPYLV     +    +G YD  MT     +    A + LH++  YG+N HHI+E  FK+ ARAL+ A   D    G +PS+KG L"
assert record.rounds[0].alignments[5].hsps[0].sbjct=="MSRVGRVERTTKETSVLVEIDLDGTGKTDIATGVGFYDHMLDQLGRHGLFDLTVKTDGD---LHIDSHHTIEDTALALGAAFRQALGDKVGIYRFGNCTVPLDESLAQVTVDLSGRPYLVHTEPENMAPMIGEYDVTMTRHILESFVAQAQVALHVHVPYGRNAHHIVECQFKALARALRYASERDPRAAGILPSTKGAL"
assert record.rounds[0].alignments[5].hsps[0].query_start==1
assert record.rounds[0].alignments[5].hsps[0].query_end==200
assert record.rounds[0].alignments[5].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[5].hsps[0].sbjct_end==197
assert record.rounds[0].alignments[6].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[6].hsps[0].match=="RI+ + R T ET +++++NLDGTG    +TGI FLDHML  ++ H   DL +   GD E   +D HH  EDV I LG+ +++ LG++ GI R+G+F  P+DEALV   LD SGRP+L +   +   +++G YDT++  EFF AL  ++ +TLH+ +  G N+HHIIE  FK+ ARA + A+ +D  + G IPSSKGVL"
assert record.rounds[0].alignments[6].hsps[0].sbjct=="RIASVHRITGETNVQVTVNLDGTGICKAATGIPFLDHMLHQISSHGLIDLDVQAKGDWE---IDDHHTNEDVGITLGQALAKALGDRKGIVRFGNFLAPLDEALVQVALDFSGRPHLSYGLQIP-TERVGTYDTQLVREFFVALVNHSQMTLHIRQLDGINSHHIIEATFKAFARAARMAIEVDPRRAGTIPSSKGVL"
assert record.rounds[0].alignments[6].hsps[0].query_start==3
assert record.rounds[0].alignments[6].hsps[0].query_end==200
assert record.rounds[0].alignments[6].hsps[0].sbjct_start==16
assert record.rounds[0].alignments[6].hsps[0].sbjct_end==209
assert record.rounds[0].alignments[7].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[7].hsps[0].match=="R +H+ RNTKETQI++ + LD  G + I+TG+GF DHML  +  H  F ++I   GD   + +D HH +ED  +ALG+ +   LG+K GI R+G F +PMDE L  C LDISGRP+L + A+ +  Q++G   TEM E FFR+L++  G+TLHL    G+N HH +E +FK+  R L+QA+ ++      +PSSKGVL"
assert record.rounds[0].alignments[7].hsps[0].sbjct=="RYAHVVRNTKETQIDVQVWLDREGGSKINTGVGFFDHMLDQIATHGGFRMEINVKGD---LYIDDHHTVEDTGLALGEALKIALGDKRGICRFG-FVLPMDECLARCALDISGRPHLEYKAEFT-YQRVGDLSTEMIEHFFRSLSYTMGVTLHLKTK-GKNDHHRVESLFKAFGRTLRQAIRVEGD---TLPSSKGVL"
assert record.rounds[0].alignments[7].hsps[0].query_start==3
assert record.rounds[0].alignments[7].hsps[0].query_end==200
assert record.rounds[0].alignments[7].hsps[0].sbjct_start==167
assert record.rounds[0].alignments[7].hsps[0].sbjct_end==355
assert record.rounds[0].alignments[8].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[8].hsps[0].match=="R + + R TKET + +S+NL G+G   ++TG+ FLDHML  +  H   DL++   GD E   +D HH  EDV I LG+ ++E LG++ GI R+G F  P+DEALV   LD SGRP+L +   +   +++G YDT++  EFF A+  ++ +TLH+ +  G N+HHIIE  FK+ ARA++ A+ +D  +   IPSSKGVL"
assert record.rounds[0].alignments[8].hsps[0].sbjct=="RAAAVHRVTKETDVRVSLNLMGSGLCHVATGVPFLDHMLHQIASHGLIDLEVNATGDIE---IDDHHTNEDVGITLGQALAEALGDRRGINRFGHFIAPLDEALVQVTLDFSGRPHLSYGLQIP-TERVGTYDTQLVREFFVAVVNHSQMTLHIRQLDGINSHHIIEATFKAFARAMRMAIEVDPRRADTIPSSKGVL"
assert record.rounds[0].alignments[8].hsps[0].query_start==3
assert record.rounds[0].alignments[8].hsps[0].query_end==200
assert record.rounds[0].alignments[8].hsps[0].sbjct_start==17
assert record.rounds[0].alignments[8].hsps[0].sbjct_end==210
assert record.rounds[0].alignments[9].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[9].hsps[0].match=="R + + R TKET I++ + LD  G  +I TG+GF DHML  +  H  F + +   GD   + +D HH +ED A+ALG+ + + +G+K GI R+G F +PMDE    C LD+SGRP++ F+A      K+G + TE+TE FF++LAF+   TLHLN   G N HH IE +FK+  R L+QA+ I+ +   E+PSSKGVL"
assert record.rounds[0].alignments[9].hsps[0].sbjct=="RFAEVIRQTKETDIKVQVWLDEAGVNEIKTGVGFFDHMLDQIATHGGFRMNVQCKGD---LWIDEHHTVEDTALALGQALKQAVGDKRGIARFG-FVLPMDECKAECALDLSGRPWIKFNACFK-RDKVGDFSTELTEHFFQSLAFSMLATLHLNV-TGNNDHHKIESLFKAFGRTLRQAIRIEGN---EMPSSKGVL"
assert record.rounds[0].alignments[9].hsps[0].query_start==3
assert record.rounds[0].alignments[9].hsps[0].query_end==200
assert record.rounds[0].alignments[9].hsps[0].sbjct_start==174
assert record.rounds[0].alignments[9].hsps[0].sbjct_end==362
assert record.rounds[0].alignments[10].hsps[0].query=="TRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVF--------HADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[10].hsps[0].match=="+R + I R T+E+ I + ++LDGTGQ  + TG+ F DHMLT L  H+ FDL +   GD E   ++ HH IED AIALG  + + LG+K GIRR+G   IPMDE L    +D+SGRPY V         H  ++G+     Y T +    F +LA NA I LH+   YG++ HHI E  +K+ ARAL+QAV  D  +V  +PS+KG L"
assert record.rounds[0].alignments[10].hsps[0].sbjct=="SRRARIERRTRESDIVIELDLDGTGQVAVDTGVPFYDHMLTALGSHASFDLTVRATGDVE---IEAHHTIEDTAIALGTALGQALGDKRGIRRFGDAFIPMDETLAHAAVDLSGRPYCVHTGEPDHLQHTTIAGSSV--PYHTVINRHVFESLAANARIALHVRVLYGRDPHHITEAQYKAVARALRQAVEPD-PRVSGVPSTKGAL"
assert record.rounds[0].alignments[10].hsps[0].query_start==2
assert record.rounds[0].alignments[10].hsps[0].query_end==200
assert record.rounds[0].alignments[10].hsps[0].sbjct_start==10
assert record.rounds[0].alignments[10].hsps[0].sbjct_end==210
assert record.rounds[0].alignments[11].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[11].hsps[0].match=="R +H+ RNTKETQI++S+ LD  G + I+TG+GF DHML  +  H  F ++I   GD   + +D HH +ED  +AL + +   L +K GI R+G F +PMDE L  C LDISGRP+L + A+ +  Q++G   TEM E FFR+L++  G+TLHL    G+N HH +E +FK+  R ++QA+ ++      +PSSKGVL"
assert record.rounds[0].alignments[11].hsps[0].sbjct=="RYAHVVRNTKETQIDVSVWLDREGNSKINTGVGFFDHMLDQIATHGGFRMEITVKGD---LYIDDHHTVEDTGLALREALKLALRDKRGICRFG-FVLPMDECL-ACALDISGRPHLEYKAEFT-YQRVGNLSTEMIEHFFRSLSYTMGVTLHLKTK-GKNDHHRVESLFKAFGRTVRQAIRVEGD---TLPSSKGVL"
assert record.rounds[0].alignments[11].hsps[0].query_start==3
assert record.rounds[0].alignments[11].hsps[0].query_end==200
assert record.rounds[0].alignments[11].hsps[0].sbjct_start==167
assert record.rounds[0].alignments[11].hsps[0].sbjct_end==354
assert record.rounds[0].alignments[12].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[12].hsps[0].match=="R S  TR T ET +++ + +DG+G++ ++TG+GFLDHML  +  H   DL++   GD E   +D HH +EDVA+ LG+ + E LG+K GIRR     +PMD+AL T  LD+SGRPY V   +   +  +G   ++    F  +LA +A + +H +   G+N HH  E +FK+ A A++ AV ++    GEIPS+KG L"
assert record.rounds[0].alignments[12].hsps[0].sbjct=="RRSMKTRETLETHVKVDLEIDGSGKSSVNTGLGFLDHMLESVARHGLLDLEVEARGDLE---VDDHHTVEDVALTLGEALREALGDKSGIRRMAHAMVPMDDALATVALDLSGRPYTVLELEFD-DAVIGDVKSQNIGHFIESLAVSAAMNIHASVR-GRNDHHKAEALFKALALAIRDAVRVEH---GEIPSTKGKL"
assert record.rounds[0].alignments[12].hsps[0].query_start==3
assert record.rounds[0].alignments[12].hsps[0].query_end==200
assert record.rounds[0].alignments[12].hsps[0].sbjct_start==5
assert record.rounds[0].alignments[12].hsps[0].sbjct_end==194
assert record.rounds[0].alignments[13].hsps[0].query=="RISHITRNTKETQIELSINLDGTG-----QADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[13].hsps[0].match=="RI+ + R T ET I  +I+LD        + ++STGIGFLDHM T L  H    L++   GD   + +D HH  ED A+ALG+   + LG + GI+RYG    P+DE+L    +DIS RPY + H   +  +K+G   TEM     ++ AF AG+TLH++   G+N HHI E  FK+ A A++ A+S   +   ++PS+KGVL"
assert record.rounds[0].alignments[13].hsps[0].sbjct=="RIASVERTTSETHISCTIDLDHIPGVTEQKINVSTGIGFLDHMFTALAKHGGMSLQLQCKGD---LHIDDHHTAEDCALALGEAFKKALGERKGIKRYGYAYAPLDESLSRAVIDISSRPYFMCHLPFT-REKVGDLSTEMVSHLLQSFAFAAGVTLHIDSIRGENNHHIAESAFKALALAIRMAIS--RTGGDDVPSTKGVL"
assert record.rounds[0].alignments[13].hsps[0].query_start==3
assert record.rounds[0].alignments[13].hsps[0].query_end==200
assert record.rounds[0].alignments[13].hsps[0].sbjct_start==4
assert record.rounds[0].alignments[13].hsps[0].sbjct_end==200
assert record.rounds[0].alignments[14].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHML-TLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[14].hsps[0].match=="M+R ++ITR TKET+IE+ +++D  G+  +ST I F +HML TLLT+ +     I+   D   +  D HH++EDVAI LG  I   LG+K GI+R+    IPMD+ALV   LDIS R     + +L  ++ +GG  TE    FF++ A+N+GITLH+++  G NTHHIIE  FK+   AL +A  I ++   EI S+KG++"
assert record.rounds[0].alignments[14].hsps[0].sbjct=="MSRSANITRETKETKIEVLLDIDRKGEVKVSTPIPFFNHMLITLLTYMNS--TAIVSATDK--LPYDDHHIVEDVAITLGLAIKTALGDKRGIKRFSHQIIPMDDALVLVSLDISNRGMAFVNLNLKRSE-IGGLATENVPHFFQSFAYNSGITLHISQLSGYNTHHIIEASFKALGLALYEATRIVDN---EIRSTKGII"
assert record.rounds[0].alignments[14].hsps[0].query_start==1
assert record.rounds[0].alignments[14].hsps[0].query_end==200
assert record.rounds[0].alignments[14].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[14].hsps[0].sbjct_end==193
assert record.rounds[0].alignments[15].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[15].hsps[0].match=="R ++I+R TKET I + ++LDGTG++ +S+GIGFLDHMLT L  HS FDL++   GD     +D HH  ED A+ LG+     LG++ GI R+GS  +P+DEAL    +DIS R +   +  L     +G   +EM   FF + A  A  TLH++   G+N HH  E  FK+ A AL+ AV  D +    +PS+KGVL"
assert record.rounds[0].alignments[15].hsps[0].sbjct=="REANISRVTKETSISVKLSLDGTGKSKVSSGIGFLDHMLTALAKHSRFDLELDCKGD---TWIDDHHTTEDCALTLGEAFDVALGDRAGIARFGSACVPLDEALSRAIVDISSRAHSEINLQLV-RPSVGELSSEMITHFFESFASAALXTLHVDVLRGRNDHHRAEASFKALAVALRTAVKHDAT--AGVPSTKGVL"
assert record.rounds[0].alignments[15].hsps[0].query_start==3
assert record.rounds[0].alignments[15].hsps[0].query_end==200
assert record.rounds[0].alignments[15].hsps[0].sbjct_start==260
assert record.rounds[0].alignments[15].hsps[0].sbjct_end==451
assert record.rounds[0].alignments[16].hsps[0].query=="RISHITRNTKETQIELSINLD-----------------------GTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[16].hsps[0].match=="R + + RNT ET+I ++I LD                       G     + TGIGFLDHM   L  H+ + L++   GD   + +D HH  ED AIALG    + +GN  G++R+G    P+DEAL    +D+SGRPY V    L   +K+G    EM      + +  AGITLH+   YG N HH  E  FKS A A++ A S+  S   E+PS+KGVL"
assert record.rounds[0].alignments[16].hsps[0].sbjct=="RRAFVERNTNETKISVAIALDKAPLPEESNFIDELITSKHANQKGEQVIQVDTGIGFLDHMYHALAKHAGWSLRLYSRGD---LIIDDHHTAEDTAIALGIAFKQAMGNFAGVKRFGHAYCPLDEALSRSVVDLSGRPYAVIDLGLK-REKVGELSCEMIPHLLYSFSVAAGITLHVTCLYGSNDHHRAESAFKSLAVAMRAATSLTGS--SEVPSTKGVL"
assert record.rounds[0].alignments[16].hsps[0].query_start==3
assert record.rounds[0].alignments[16].hsps[0].query_end==200
assert record.rounds[0].alignments[16].hsps[0].sbjct_start==2
assert record.rounds[0].alignments[16].hsps[0].sbjct_end==216
assert record.rounds[0].alignments[17].hsps[0].query=="MTRISHITRNTKETQIELSINLDG---------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[17].hsps[0].match=="M+R + I R T ET+I++++NLDG                           +   ++ TGIGFLDHM+  L  HS + L +   GD   + +D HH  EDV I+LG    + LG   G++R+G    P+DEAL    +D+S RP+ V    L   +K+G   TEM      + A   GIT+H++   G N HH  E  FK+ A A+K+A+S  ++   +IPS+KGVL"
assert record.rounds[0].alignments[17].hsps[0].sbjct=="MSREALINRITNETKIQIALNLDGGKLELKESIFPNQSIIIDEHHAKQVSGSQYINVQTGIGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTAEDVGISLGMAFKQALGQIKGVKRFGHGFAPLDEALSRAVVDLSNRPFAVIELGLK-REKIGDLSTEMIPHVLESFAGAVGITIHVDCLRGFNDHHRAESAFKALAIAIKEAIS--KTGKNDIPSTKGVL"
assert record.rounds[0].alignments[17].hsps[0].query_start==1
assert record.rounds[0].alignments[17].hsps[0].query_end==200
assert record.rounds[0].alignments[17].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[17].hsps[0].sbjct_end==221
assert record.rounds[0].alignments[18].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKG"
assert record.rounds[0].alignments[18].hsps[0].match=="RI  + R TKET I L IN+DGTG+  I TGI F DH+L     H  FDL +   GD E   +D HH +EDV I LG  +++    K  I R+G   IPMD+A  T  +D+SGR Y V + +    + +G   TE    FF ++A    + +H  E  G+N HH  E +FK+   AL  A  IDE K   + S+KG"
assert record.rounds[0].alignments[18].hsps[0].sbjct=="RIFEVMRETKETNIYLKINIDGTGKYKIDTGIPFFDHLLASFAKHGCFDLIVKARGDLE---IDDHHTVEDVGICLGLALNQI--EKRNIFRFGWAIIPMDDARATVAIDLSGRSYCVGNYE-PKREFVGDLATENINHFFESVASYGMLNIHY-EVIGKNEHHKAEALFKAFGVALDLATKIDERK--GVISTKG"
assert record.rounds[0].alignments[18].hsps[0].query_start==3
assert record.rounds[0].alignments[18].hsps[0].query_end==198
assert record.rounds[0].alignments[18].hsps[0].sbjct_start==7
assert record.rounds[0].alignments[18].hsps[0].sbjct_end==193
assert record.rounds[0].alignments[19].hsps[0].query=="RISHITRNTKETQIELSINLDG------------------TGQA------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[19].hsps[0].match=="R + I+R T ET+I+++I+L+G                    QA      DI TG+GFLDHM+  L  HS + L +   GD   + +D HH  ED  IALG+   E +G   G++R+G+   P+DEAL    +D+S RP+ V    L   + +G   TEM   F  + A  A ITLH++   G N HH  E  FK+ A A+++A+S   +   ++PS+KGVL"
assert record.rounds[0].alignments[19].hsps[0].sbjct=="RKAFISRITNETKIQIAISLNGGYIQIKDSILPAKKDDDVASQATQSQVIDIHTGVGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTTEDCGIALGQAFKEAMGAVRGVKRFGTGFAPLDEALSRAVVDLSNRPFAVIDLGLK-REMIGDLSTEMIPHFLESFAEAARITLHVDCLRGFNDHHRSESAFKALAVAIREAIS--SNGTNDVPSTKGVL"
assert record.rounds[0].alignments[19].hsps[0].query_start==3
assert record.rounds[0].alignments[19].hsps[0].query_end==200
assert record.rounds[0].alignments[19].hsps[0].sbjct_start==16
assert record.rounds[0].alignments[19].hsps[0].sbjct_end==231
assert record.rounds[0].alignments[20].hsps[0].query=="RISHITRNTKETQIELSINLDG-------------------TGQA------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[20].hsps[0].match=="R S I R T ET+I+++++LDG                     QA       ++TGIGFLDHML  L  H  + + I   GD   + +D HH  ED  IALG    E LG+  GI+R+GS   P+DEAL    +D+S RPY V    L   +K+G    EM      + A  A +T+H++   G N HH  E  FK+ A A+K+A+S   +   +IPS+KGVL"
assert record.rounds[0].alignments[20].hsps[0].sbjct=="RSSLIKRITNETKIQIALSLDGGPVSLAQSLFKDKDYSAEHAAQATSSQFISVNTGIGFLDHMLHALAKHGGWSVIIECVGD---LHIDDHHSAEDTGIALGMAFKEALGHVRGIKRFGSGFAPLDEALSRAVIDMSNRPYAVVDLGLK-REKIGDLSCEMIPHVLESFAQGAHVTMHVDCLRGFNDHHRAESAFKALAIAIKEAIS--SNGTDDIPSTKGVL"
assert record.rounds[0].alignments[20].hsps[0].query_start==3
assert record.rounds[0].alignments[20].hsps[0].query_end==200
assert record.rounds[0].alignments[20].hsps[0].sbjct_start==7
assert record.rounds[0].alignments[20].hsps[0].sbjct_end==223
assert record.rounds[0].alignments[21].hsps[0].query=="ISHITRNTKETQIELSINLDG-----------------TGQA-DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[21].hsps[0].match=="++ + R T+ET I+L+++LDG                  GQ   + TG+GFLDHMLT L  H  + L +   GD   + +D HH +ED  IALG+   E LG+  GI+R+G    P+DEAL    +D S RP+ V    L   +++G   TEM   F  + A    IT+H++   G N HH  E  FK+ A A+++A +   +   ++PS+KGVL"
assert record.rounds[0].alignments[21].hsps[0].sbjct=="MAFVKRVTQETNIQLALDLDGGSVSVRESILGKEYASGDGQTIHVHTGVGFLDHMLTALAKHGGWSLILECIGD---LHIDDHHTVEDCGIALGQAFKEALGSVRGIKRFGHGFAPLDEALSRAVVDFSNRPFAVVELGLK-RERIGQLSTEMIPHFLESFATEGRITMHVDCLRGTNDHHRSESGFKALAIAIREART--PTGRDDVPSTKGVL"
assert record.rounds[0].alignments[21].hsps[0].query_start==4
assert record.rounds[0].alignments[21].hsps[0].query_end==200
assert record.rounds[0].alignments[21].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[21].hsps[0].sbjct_end==209
assert record.rounds[0].alignments[22].hsps[0].query=="ITRNTKETQIELSINLDGTGQA------------------------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[22].hsps[0].match=="+ R T ET+I+++I+L G   A                        ++ TGIGFLDHM+  L  HS + L +   GD   + +D HH  ED  IALG+   E LG   G++R+GS   P+DEAL    +D+S RPY V    L   +K+G    EM   F  + A  + ITLH++   G+N HH  E  FK+ A A+++A S   +   ++PS+KGVL"
assert record.rounds[0].alignments[22].hsps[0].sbjct=="VKRITNETKIQIAISLKGGPLAIEHSIFPEKEAEAVAEQATQSQVINVHTGIGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTTEDCGIALGQAFKEALGAVRGVKRFGSGFAPLDEALSRAVVDLSNRPYAVVELGLQ-REKVGDLSCEMIPHFLESFAEASRITLHVDCLRGKNDHHRSESAFKALAVAIREATS--PNGTNDVPSTKGVL"
assert record.rounds[0].alignments[22].hsps[0].query_start==7
assert record.rounds[0].alignments[22].hsps[0].query_end==200
assert record.rounds[0].alignments[22].hsps[0].sbjct_start==8
assert record.rounds[0].alignments[22].hsps[0].sbjct_end==219
assert record.rounds[0].alignments[23].hsps[0].query=="RISHITRNTKETQIELSINLDG------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALK"
assert record.rounds[0].alignments[23].hsps[0].match=="R + ++R+T ET I++++++DG                        +    I+TGIGFLDHML  L  H+ + + +   GD   + +D HH  ED  IA+G   ++ LG   G+ R+G    P+DEAL    +D+S RPY V    L   +KLG    EM     ++ A  A ITLH++   G N HH  E  FK+ A A++"
assert record.rounds[0].alignments[23].hsps[0].sbjct=="RAAALSRDTNETSIQIALSIDGGELPQDADPRLLEASSAHASQTSKSQVISINTGIGFLDHMLHALAKHAGWSMALNCKGD---LHIDDHHTAEDCCIAVGTTFAKALGALTGVARFGYAYAPLDEALSRAVVDLSNRPYTVVDLGLK-REKLGELSCEMIPHCLQSFAQAARITLHVDCLRGDNDHHRAESAFKALAVAVR"
assert record.rounds[0].alignments[23].hsps[0].query_start==3
assert record.rounds[0].alignments[23].hsps[0].query_end==180
assert record.rounds[0].alignments[23].hsps[0].sbjct_start==8
assert record.rounds[0].alignments[23].hsps[0].sbjct_end==205
assert record.rounds[0].alignments[24].hsps[0].query=="GDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQK---LGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDES"
assert record.rounds[0].alignments[24].hsps[0].match=="G+   +G  PH  + +++   G  +S  LG+   +      T+   + L T DL+   RPY+ + A ++ N K      YD     +++R +     + L+  +   Q+  HI E    S  R  KQA S +E+"
assert record.rounds[0].alignments[24].hsps[0].sbjct=="GNLHQLGEKPHRAMVELSKTYGPLMSLKLGSVTTVVATSVETVR--DVLKTYDLECCSRPYMTYPARITYNLKDLVFSPYD-----KYWRQVRKLTVVELYTAKRV-QSFRHIREEEVASFVRFNKQAASSEET"
assert record.rounds[0].alignments[24].hsps[0].query_start==58
assert record.rounds[0].alignments[24].hsps[0].query_end==188
assert record.rounds[0].alignments[24].hsps[0].sbjct_start==40
assert record.rounds[0].alignments[24].hsps[0].sbjct_end==165
assert record.rounds[0].alignments[25].hsps[0].query=="HSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLG"
assert record.rounds[0].alignments[25].hsps[0].match=="HS+    I+G G H   G DPHH   D  +  G  +  D+G   G"
assert record.rounds[0].alignments[25].hsps[0].sbjct=="HSEVAFVIVGSGPH---GADPHHGYSDRELREGDIVVVDIGGTYG"
assert record.rounds[0].alignments[25].hsps[0].query_start==47
assert record.rounds[0].alignments[25].hsps[0].query_end==91
assert record.rounds[0].alignments[25].hsps[0].sbjct_start==195
assert record.rounds[0].alignments[25].hsps[0].sbjct_end==236
assert record.rounds[0].alignments[26].hsps[0].query=="PYLVF---HADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYG"
assert record.rounds[0].alignments[26].hsps[0].match=="PYL+    H D+ GN +  G+  +M +E    L FN  I L  +  YG"
assert record.rounds[0].alignments[26].hsps[0].sbjct=="PYLMLKGNHQDMEGNDRYEGFCVDMLKELAEILRFNYKIRLVGDGVYG"
assert record.rounds[0].alignments[26].hsps[0].query_start==117
assert record.rounds[0].alignments[26].hsps[0].query_end==161
assert record.rounds[0].alignments[26].hsps[0].sbjct_start==427
assert record.rounds[0].alignments[26].hsps[0].sbjct_end==474
assert record.rounds[0].alignments[27].hsps[0].query=="GKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLG---GYDTEMTEEFFRA"
assert record.rounds[0].alignments[27].hsps[0].match=="GK ISED     GI + G     +D   VT + D  G    +  + L  N++ G    YD E T EF RA"
assert record.rounds[0].alignments[27].hsps[0].sbjct=="GKEISEDKWEDFGISQRGEEKFFIDAEKVTVEFD--GFQAKIQMSSLYKNKQCGLCGHYDGEKTNEFRRA"
assert record.rounds[0].alignments[27].hsps[0].query_start==79
assert record.rounds[0].alignments[27].hsps[0].query_end==145
assert record.rounds[0].alignments[27].hsps[0].sbjct_start==1436
assert record.rounds[0].alignments[27].hsps[0].sbjct_end==1503
assert record.rounds[0].alignments[28].hsps[0].query=="PYLVF---HADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYG"
assert record.rounds[0].alignments[28].hsps[0].match=="PYL+    H ++ GN +  G+  +M +E    L FN  I L  +  YG"
assert record.rounds[0].alignments[28].hsps[0].sbjct=="PYLMLKGNHQEMEGNDRYEGFCVDMLKELAEILRFNYKIRLVGDGVYG"
assert record.rounds[0].alignments[28].hsps[0].query_start==117
assert record.rounds[0].alignments[28].hsps[0].query_end==161
assert record.rounds[0].alignments[28].hsps[0].sbjct_start==427
assert record.rounds[0].alignments[28].hsps[0].sbjct_end==474
assert record.rounds[0].alignments[29].hsps[0].query=="RYGSFTIPMDEAL-----VTCDLDISGRPYLVFHADLSGNQKL--GGYDTEMTEEFFRALAFNAG"
assert record.rounds[0].alignments[29].hsps[0].match=="R GS ++P++E          D DI G  Y   H D+  NQ+L  G  D  +++   +   FN+G"
assert record.rounds[0].alignments[29].hsps[0].sbjct=="RNGSNSLPLNEKSNEGESTNVDQDIEGDEYHRLHEDILKNQELDDGSLDDLLSQIIPKITNFNSG"
assert record.rounds[0].alignments[29].hsps[0].query_start==94
assert record.rounds[0].alignments[29].hsps[0].query_end==151
assert record.rounds[0].alignments[29].hsps[0].sbjct_start==1141
assert record.rounds[0].alignments[29].hsps[0].sbjct_end==1205
assert record.rounds[1].alignments[0].hsps[0].query=="TRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[0].hsps[0].match=="TR+  + R TKET + + INLDG+G AD STGI FLDHML  L  H  FD+ +   GD   V +D HH  EDVA+A+G  + + LG++ GI R+G F+ P+DEAL+   LD+SGRP+L ++ D+   Q++G YDT++ E F +++   +G+TLH+ +  G+N+HHIIE  FK+ ARAL+QA   D  + G +PSSKGVL"
assert record.rounds[1].alignments[0].hsps[0].sbjct=="TRVGEVKRVTKETNVSVKINLDGSGVADSSTGIPFLDHMLDQLASHGLFDVHVKATGD---VHIDDHHTNEDVALAIGTALLQALGDRKGINRFGDFSAPLDEALIHVSLDLSGRPHLSYNLDIP-TQRVGTYDTQVVEHFLQSIVNTSGMTLHIRQLAGRNSHHIIEATFKAFARALRQATEYDPRRRGSVPSSKGVL"
assert record.rounds[1].alignments[0].hsps[0].query_start==2
assert record.rounds[1].alignments[0].hsps[0].query_end==200
assert record.rounds[1].alignments[0].hsps[0].sbjct_start==84
assert record.rounds[1].alignments[0].hsps[0].sbjct_end==278
assert record.rounds[1].alignments[1].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[1].hsps[0].match=="RI  + R TKET + + INLDGTG AD S+GI FLDHML  L  H  FD+ +   GD   V +D HH  ED+A+A+G  + + LG + GI R+G FT P+DEAL+   LD+SGRPYL ++ ++   Q++G YDT++ E FF++L   +G+TLH+ +  G+N+HHIIE  FK+ ARAL+QA   D  + G IPSSKGVL"
assert record.rounds[1].alignments[1].hsps[0].sbjct=="RIGEVKRVTKETNVSVKINLDGTGVADSSSGIPFLDHMLDQLASHGLFDVHVRATGD---VHIDDHHTNEDIALAIGTALLKALGERKGINRFGDFTAPLDEALIHVSLDLSGRPYLGYNLEIP-TQRVGTYDTQLVEHFFQSLVNTSGMTLHIRQLAGENSHHIIEATFKAFARALRQATETDPRRGGTIPSSKGVL"
assert record.rounds[1].alignments[1].hsps[0].query_start==3
assert record.rounds[1].alignments[1].hsps[0].query_end==200
assert record.rounds[1].alignments[1].hsps[0].sbjct_start==74
assert record.rounds[1].alignments[1].hsps[0].sbjct_end==267
assert record.rounds[1].alignments[2].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[2].hsps[0].match=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[2].hsps[0].sbjct=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[2].hsps[0].query_start==1
assert record.rounds[1].alignments[2].hsps[0].query_end==200
assert record.rounds[1].alignments[2].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[2].hsps[0].sbjct_end==200
assert record.rounds[1].alignments[3].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[3].hsps[0].match=="R + + R TKET + +S+NL G+G   ++TG+ FLDHML  +  H   DL++   GD E   +D HH  EDV I LG+ ++E LG++ GI R+G F  P+DEALV   LD SGRP+L +   +   +++G YDT++  EFF A+  ++ +TLH+ +  G N+HHIIE  FK+ ARA++ A+ +D  +   IPSSKGVL"
assert record.rounds[1].alignments[3].hsps[0].sbjct=="RAAAVHRVTKETDVRVSLNLMGSGLCHVATGVPFLDHMLHQIASHGLIDLEVNATGDIE---IDDHHTNEDVGITLGQALAEALGDRRGINRFGHFIAPLDEALVQVTLDFSGRPHLSYGLQIP-TERVGTYDTQLVREFFVAVVNHSQMTLHIRQLDGINSHHIIEATFKAFARAMRMAIEVDPRRADTIPSSKGVL"
assert record.rounds[1].alignments[3].hsps[0].query_start==3
assert record.rounds[1].alignments[3].hsps[0].query_end==200
assert record.rounds[1].alignments[3].hsps[0].sbjct_start==17
assert record.rounds[1].alignments[3].hsps[0].sbjct_end==210
assert record.rounds[1].alignments[4].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[4].hsps[0].match=="RI+ + R T ET +++++NLDGTG    +TGI FLDHML  ++ H   DL +   GD E   +D HH  EDV I LG+ +++ LG++ GI R+G+F  P+DEALV   LD SGRP+L +   +   +++G YDT++  EFF AL  ++ +TLH+ +  G N+HHIIE  FK+ ARA + A+ +D  + G IPSSKGVL"
assert record.rounds[1].alignments[4].hsps[0].sbjct=="RIASVHRITGETNVQVTVNLDGTGICKAATGIPFLDHMLHQISSHGLIDLDVQAKGDWE---IDDHHTNEDVGITLGQALAKALGDRKGIVRFGNFLAPLDEALVQVALDFSGRPHLSYGLQIP-TERVGTYDTQLVREFFVALVNHSQMTLHIRQLDGINSHHIIEATFKAFARAARMAIEVDPRRAGTIPSSKGVL"
assert record.rounds[1].alignments[4].hsps[0].query_start==3
assert record.rounds[1].alignments[4].hsps[0].query_end==200
assert record.rounds[1].alignments[4].hsps[0].sbjct_start==16
assert record.rounds[1].alignments[4].hsps[0].sbjct_end==209
assert record.rounds[1].alignments[5].hsps[0].query=="SHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[5].hsps[0].match=="  + R TKET + + INLDGTG A+ STGI FLDHML  L  H  FD+ +   GD     +D HH  ED+A+A+G  + + LG++ GI R+G FT P+DEA V   LD+SGRP+L     +   +++G YDT++ E FF++L   +G+TLH+ +  G N+HHIIE  FK+ ARAL+QA   D  + G +PSSKGVL"
assert record.rounds[1].alignments[5].hsps[0].sbjct=="GEVKRVTKETNVHVKINLDGTGVANSSTGIPFLDHMLDQLASHGLFDVYVKATGDT---HIDDHHSNEDIALAIGTALLQALGDRKGINRFGHFTAPLDEAAVEVILDLSGRPHLSCGLSIP-TERVGTYDTQLVEHFFQSLVNTSGMTLHIRQLAGNNSHHIIEATFKAFARALRQATEYDLRRQGTMPSSKGVL"
assert record.rounds[1].alignments[5].hsps[0].query_start==5
assert record.rounds[1].alignments[5].hsps[0].query_end==200
assert record.rounds[1].alignments[5].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[5].hsps[0].sbjct_end==192
assert record.rounds[1].alignments[6].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[6].hsps[0].match=="M+R+  + R TKET + + I+LDGTG+ DI+TG+GF DHML  L  H  FDL +   GD   + +D HH IED A+ALG    + LG+K+GI R+G+ T+P+DE+L    +D+SGRPYLV     +    +G YD  MT     +    A + LH++  YG+N HHI+E  FK+ ARAL+ A   D    G +PS+KG L"
assert record.rounds[1].alignments[6].hsps[0].sbjct=="MSRVGRVERTTKETSVLVEIDLDGTGKTDIATGVGFYDHMLDQLGRHGLFDLTVKTDGD---LHIDSHHTIEDTALALGAAFRQALGDKVGIYRFGNCTVPLDESLAQVTVDLSGRPYLVHTEPENMAPMIGEYDVTMTRHILESFVAQAQVALHVHVPYGRNAHHIVECQFKALARALRYASERDPRAAGILPSTKGAL"
assert record.rounds[1].alignments[6].hsps[0].query_start==1
assert record.rounds[1].alignments[6].hsps[0].query_end==200
assert record.rounds[1].alignments[6].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[6].hsps[0].sbjct_end==197
assert record.rounds[1].alignments[7].hsps[0].query=="SHITRNTKETQIELSINLDG------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[7].hsps[0].match=="+ + R T ET+I+++I+L G                        +   ++ TGIGFLDHM+  L  HS + L +   GD   + +D HH  ED  IALG+   E LG   G++R+GS   P+DEAL    +D+S RPY V    L   +K+G    EM   F  + A  + ITLH++   G+N HH  E  FK+ A A+++A S   +   ++PS+KGVL"
assert record.rounds[1].alignments[7].hsps[0].sbjct=="ALVKRITNETKIQIAISLKGGPLAIEHSIFPEKEAEAVAEQATQSQVINVHTGIGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTTEDCGIALGQAFKEALGAVRGVKRFGSGFAPLDEALSRAVVDLSNRPYAVVELGL-QREKVGDLSCEMIPHFLESFAEASRITLHVDCLRGKNDHHRSESAFKALAVAIREATS--PNGTNDVPSTKGVL"
assert record.rounds[1].alignments[7].hsps[0].query_start==5
assert record.rounds[1].alignments[7].hsps[0].query_end==200
assert record.rounds[1].alignments[7].hsps[0].sbjct_start==6
assert record.rounds[1].alignments[7].hsps[0].sbjct_end==219
assert record.rounds[1].alignments[8].hsps[0].query=="MTRISHITRNTKETQIELSINLDG---------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[8].hsps[0].match=="M+R + I R T ET+I++++NLDG                           +   ++ TGIGFLDHM+  L  HS + L +   GD   + +D HH  EDV I+LG    + LG   G++R+G    P+DEAL    +D+S RP+ V    L   +K+G   TEM      + A   GIT+H++   G N HH  E  FK+ A A+K+A+S  ++   +IPS+KGVL"
assert record.rounds[1].alignments[8].hsps[0].sbjct=="MSREALINRITNETKIQIALNLDGGKLELKESIFPNQSIIIDEHHAKQVSGSQYINVQTGIGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTAEDVGISLGMAFKQALGQIKGVKRFGHGFAPLDEALSRAVVDLSNRPFAVIELGLK-REKIGDLSTEMIPHVLESFAGAVGITIHVDCLRGFNDHHRAESAFKALAIAIKEAIS--KTGKNDIPSTKGVL"
assert record.rounds[1].alignments[8].hsps[0].query_start==1
assert record.rounds[1].alignments[8].hsps[0].query_end==200
assert record.rounds[1].alignments[8].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[8].hsps[0].sbjct_end==221
assert record.rounds[1].alignments[9].hsps[0].query=="RISHITRNTKETQIELSINLDG------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[9].hsps[0].match=="R + I+R T ET+I+++I+L+G                        +   DI TG+GFLDHM+  L  HS + L +   GD   + +D HH  ED  IALG+   E +G   G++R+G+   P+DEAL    +D+S RP+ V    L   + +G   TEM   F  + A  A ITLH++   G N HH  E  FK+ A A+++A+S   +   ++PS+KGVL"
assert record.rounds[1].alignments[9].hsps[0].sbjct=="RKAFISRITNETKIQIAISLNGGYIQIKDSILPAKKDDDVASQATQSQVIDIHTGVGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTTEDCGIALGQAFKEAMGAVRGVKRFGTGFAPLDEALSRAVVDLSNRPFAVIDLGLK-REMIGDLSTEMIPHFLESFAEAARITLHVDCLRGFNDHHRSESAFKALAVAIREAIS--SNGTNDVPSTKGVL"
assert record.rounds[1].alignments[9].hsps[0].query_start==3
assert record.rounds[1].alignments[9].hsps[0].query_end==200
assert record.rounds[1].alignments[9].hsps[0].sbjct_start==16
assert record.rounds[1].alignments[9].hsps[0].sbjct_end==231
assert record.rounds[1].alignments[10].hsps[0].query=="SHITRNTKETQIELSINLDGTGQA------------------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[10].hsps[0].match=="+ + R T+ET I+L+++LDG   +                   + TG+GFLDHMLT L  H  + L +   GD   + +D HH +ED  IALG+   E LG+  GI+R+G    P+DEAL    +D S RP+ V    L   +++G   TEM   F  + A    IT+H++   G N HH  E  FK+ A A+++A     +   ++PS+KGVL"
assert record.rounds[1].alignments[10].hsps[0].sbjct=="AFVKRVTQETNIQLALDLDGGSVSVRESILGKEYASGDGQTIHVHTGVGFLDHMLTALAKHGGWSLILECIGD---LHIDDHHTVEDCGIALGQAFKEALGSVRGIKRFGHGFAPLDEALSRAVVDFSNRPFAVVELGLK-RERIGQLSTEMIPHFLESFATEGRITMHVDCLRGTNDHHRSESGFKALAIAIREA--RTPTGRDDVPSTKGVL"
assert record.rounds[1].alignments[10].hsps[0].query_start==5
assert record.rounds[1].alignments[10].hsps[0].query_end==200
assert record.rounds[1].alignments[10].hsps[0].sbjct_start==2
assert record.rounds[1].alignments[10].hsps[0].sbjct_end==209
assert record.rounds[1].alignments[11].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[11].hsps[0].match=="R +H+ RNTKETQI++ + LD  G + I+TG+GF DHML  +  H  F ++I   GD   + +D HH +ED  +ALG+ +   LG+K GI R+G F +PMDE L  C LDISGRP+L + A+ +  Q++G   TEM E FFR+L++  G+TLHL    G+N HH +E +FK+  R L+QA+ ++      +PSSKGVL"
assert record.rounds[1].alignments[11].hsps[0].sbjct=="RYAHVVRNTKETQIDVQVWLDREGGSKINTGVGFFDHMLDQIATHGGFRMEINVKGD---LYIDDHHTVEDTGLALGEALKIALGDKRGICRFG-FVLPMDECLARCALDISGRPHLEYKAEFT-YQRVGDLSTEMIEHFFRSLSYTMGVTLHLKT-KGKNDHHRVESLFKAFGRTLRQAIRVEG---DTLPSSKGVL"
assert record.rounds[1].alignments[11].hsps[0].query_start==3
assert record.rounds[1].alignments[11].hsps[0].query_end==200
assert record.rounds[1].alignments[11].hsps[0].sbjct_start==167
assert record.rounds[1].alignments[11].hsps[0].sbjct_end==355
assert record.rounds[1].alignments[12].hsps[0].query=="SHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[12].hsps[0].match=="+ I RNT ET+I +++NLDGTG  D+ TG+GFLDHML  L+ HS  DL +   GD   V +D HH  E   IA+G+ +++ +G++ GI+RYG   +PMDE L    LD S RPYL++    S   K+G  DTE+  E+F+A A  AG+TLH+   YG+N HHI+E  +K+ ARAL+  + ID  K   +PS+KG L"
assert record.rounds[1].alignments[12].hsps[0].sbjct=="ASIERNTTETRIRVAVNLDGTGVYDVKTGVGFLDHMLEQLSRHSLMDLSVAAEGD---VHIDAHHTTEHSGIAIGQAVAKAVGDRKGIQRYGHAYVPMDETLTRVALDFSNRPYLIWKVSFS-RDKIGDMDTELFREWFQAFAMAAGVTLHVECLYGENNHHIVESCYKALARALRAGIEIDPRKRDAVPSTKGTL"
assert record.rounds[1].alignments[12].hsps[0].query_start==5
assert record.rounds[1].alignments[12].hsps[0].query_end==200
assert record.rounds[1].alignments[12].hsps[0].sbjct_start==12
assert record.rounds[1].alignments[12].hsps[0].sbjct_end==203
assert record.rounds[1].alignments[13].hsps[0].query=="RISHITRNTKETQIELSINLDG-----TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[13].hsps[0].match=="RI+ + R T ET I  +I+LD        + ++STGIGFLDHM T L  H    L++   GD   + +D HH  ED A+ALG+   + LG + GI+RYG    P+DE+L    +DIS RPY + H   +  +K+G   TEM     ++ AF AG+TLH++   G+N HHI E  FK+ A A++ A+S   +   ++PS+KGVL"
assert record.rounds[1].alignments[13].hsps[0].sbjct=="RIASVERTTSETHISCTIDLDHIPGVTEQKINVSTGIGFLDHMFTALAKHGGMSLQLQCKGD---LHIDDHHTAEDCALALGEAFKKALGERKGIKRYGYAYAPLDESLSRAVIDISSRPYFMCHLPFT-REKVGDLSTEMVSHLLQSFAFAAGVTLHIDSIRGENNHHIAESAFKALALAIRMAIS--RTGGDDVPSTKGVL"
assert record.rounds[1].alignments[13].hsps[0].query_start==3
assert record.rounds[1].alignments[13].hsps[0].query_end==200
assert record.rounds[1].alignments[13].hsps[0].sbjct_start==4
assert record.rounds[1].alignments[13].hsps[0].sbjct_end==200
assert record.rounds[1].alignments[14].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQA-------------------------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[14].hsps[0].match=="R S I R T ET+I+++++LDG   +                          ++TGIGFLDHML  L  H  + + I   GD   + +D HH  ED  IALG    E LG+  GI+R+GS   P+DEAL    +D+S RPY V    L   +K+G    EM      + A  A +T+H++   G N HH  E  FK+ A A+K+A+S   +   +IPS+KGVL"
assert record.rounds[1].alignments[14].hsps[0].sbjct=="RSSLIKRITNETKIQIALSLDGGPVSLAQSLFKDKDYSAEHAAQATSSQFISVNTGIGFLDHMLHALAKHGGWSVIIECVGD---LHIDDHHSAEDTGIALGMAFKEALGHVRGIKRFGSGFAPLDEALSRAVIDMSNRPYAVVDLGLK-REKIGDLSCEMIPHVLESFAQGAHVTMHVDCLRGFNDHHRAESAFKALAIAIKEAIS--SNGTDDIPSTKGVL"
assert record.rounds[1].alignments[14].hsps[0].query_start==3
assert record.rounds[1].alignments[14].hsps[0].query_end==200
assert record.rounds[1].alignments[14].hsps[0].sbjct_start==7
assert record.rounds[1].alignments[14].hsps[0].sbjct_end==223
assert record.rounds[1].alignments[15].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[15].hsps[0].match=="R + + R TKET I++ + LD  G  +I TG+GF DHML  +  H  F + +   GD   + +D HH +ED A+ALG+ + + +G+K GI R+G F +PMDE    C LD+SGRP++ F+A      K+G + TE+TE FF++LAF+   TLHLN   G N HH IE +FK+  R L+QA+ I+ +   E+PSSKGVL"
assert record.rounds[1].alignments[15].hsps[0].sbjct=="RFAEVIRQTKETDIKVQVWLDEAGVNEIKTGVGFFDHMLDQIATHGGFRMNVQCKGD---LWIDEHHTVEDTALALGQALKQAVGDKRGIARFG-FVLPMDECKAECALDLSGRPWIKFNACFK-RDKVGDFSTELTEHFFQSLAFSMLATLHLNV-TGNNDHHKIESLFKAFGRTLRQAIRIEGN---EMPSSKGVL"
assert record.rounds[1].alignments[15].hsps[0].query_start==3
assert record.rounds[1].alignments[15].hsps[0].query_end==200
assert record.rounds[1].alignments[15].hsps[0].sbjct_start==174
assert record.rounds[1].alignments[15].hsps[0].sbjct_end==362
assert record.rounds[1].alignments[16].hsps[0].query=="RISHITRNTKETQIELSINLD-----------------------GTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[16].hsps[0].match=="R + + RNT ET+I ++I LD                       G     + TGIGFLDHM   L  H+ + L++   GD   + +D HH  ED AIALG    + +GN  G++R+G    P+DEAL    +D+SGRPY V    L   +K+G    EM      + +  AGITLH+   YG N HH  E  FKS A A++ A S+  +   E+PS+KGVL"
assert record.rounds[1].alignments[16].hsps[0].sbjct=="RRAFVERNTNETKISVAIALDKAPLPEESNFIDELITSKHANQKGEQVIQVDTGIGFLDHMYHALAKHAGWSLRLYSRGD---LIIDDHHTAEDTAIALGIAFKQAMGNFAGVKRFGHAYCPLDEALSRSVVDLSGRPYAVIDLGLK-REKVGELSCEMIPHLLYSFSVAAGITLHVTCLYGSNDHHRAESAFKSLAVAMRAATSL--TGSSEVPSTKGVL"
assert record.rounds[1].alignments[16].hsps[0].query_start==3
assert record.rounds[1].alignments[16].hsps[0].query_end==200
assert record.rounds[1].alignments[16].hsps[0].sbjct_start==2
assert record.rounds[1].alignments[16].hsps[0].sbjct_end==216
assert record.rounds[1].alignments[17].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGI-TLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[17].hsps[0].match=="R ++I+R TKET I + ++LDGTG++ +S+GIGFLDHMLT L  HS FDL++   GD     +D HH  ED A+ LG+     LG++ GI R+GS  +P+DEAL    +DIS R +   +  L     +G   +EM   FF + A  A + TLH++   G+N HH  E  FK+ A AL+ AV  D +    +PS+KGVL"
assert record.rounds[1].alignments[17].hsps[0].sbjct=="REANISRVTKETSISVKLSLDGTGKSKVSSGIGFLDHMLTALAKHSRFDLELDCKGDT---WIDDHHTTEDCALTLGEAFDVALGDRAGIARFGSACVPLDEALSRAIVDISSRAHSEINLQLV-RPSVGELSSEMITHFFESFASAA-LXTLHVDVLRGRNDHHRAEASFKALAVALRTAVKHDAT--AGVPSTKGVL"
assert record.rounds[1].alignments[17].hsps[0].query_start==3
assert record.rounds[1].alignments[17].hsps[0].query_end==200
assert record.rounds[1].alignments[17].hsps[0].sbjct_start==260
assert record.rounds[1].alignments[17].hsps[0].sbjct_end==451
assert record.rounds[1].alignments[18].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[18].hsps[0].match=="R +H+ RNTKETQI++S+ LD  G + I+TG+GF DHML  +  H  F ++I   GD   + +D HH +ED  +AL + +   L +K GI R+G F +PMDE L  C LDISGRP+L + A+ +  Q++G   TEM E FFR+L++  G+TLHL    G+N HH +E +FK+  R ++QA+ ++      +PSSKGVL"
assert record.rounds[1].alignments[18].hsps[0].sbjct=="RYAHVVRNTKETQIDVSVWLDREGNSKINTGVGFFDHMLDQIATHGGFRMEITVKGD---LYIDDHHTVEDTGLALREALKLALRDKRGICRFG-FVLPMDECLA-CALDISGRPHLEYKAEFT-YQRVGNLSTEMIEHFFRSLSYTMGVTLHLKT-KGKNDHHRVESLFKAFGRTVRQAIRVEG---DTLPSSKGVL"
assert record.rounds[1].alignments[18].hsps[0].query_start==3
assert record.rounds[1].alignments[18].hsps[0].query_end==200
assert record.rounds[1].alignments[18].hsps[0].sbjct_start==167
assert record.rounds[1].alignments[18].hsps[0].sbjct_end==354
assert record.rounds[1].alignments[19].hsps[0].query=="TRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVF--------HADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[19].hsps[0].match=="+R + I R T+E+ I + ++LDGTGQ  + TG+ F DHMLT L  H+ FDL +   GD E   ++ HH IED AIALG  + + LG+K GIRR+G   IPMDE L    +D+SGRPY V         H  ++G+     Y T +    F +LA NA I LH+   YG++ HHI E  +K+ ARAL+QAV  D   V  +PS+KG L"
assert record.rounds[1].alignments[19].hsps[0].sbjct=="SRRARIERRTRESDIVIELDLDGTGQVAVDTGVPFYDHMLTALGSHASFDLTVRATGDVE---IEAHHTIEDTAIALGTALGQALGDKRGIRRFGDAFIPMDETLAHAAVDLSGRPYCVHTGEPDHLQHTTIAGSSV--PYHTVINRHVFESLAANARIALHVRVLYGRDPHHITEAQYKAVARALRQAVEPDPR-VSGVPSTKGAL"
assert record.rounds[1].alignments[19].hsps[0].query_start==2
assert record.rounds[1].alignments[19].hsps[0].query_end==200
assert record.rounds[1].alignments[19].hsps[0].sbjct_start==10
assert record.rounds[1].alignments[19].hsps[0].sbjct_end==210
assert record.rounds[1].alignments[20].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[20].hsps[0].match=="R S  TR T ET +++ + +DG+G++ ++TG+GFLDHML  +  H   DL++   GD E   +D HH +EDVA+ LG+ + E LG+K GIRR     +PMD+AL T  LD+SGRPY V   +   +  +G   ++    F  +LA +A + +H +   G+N HH  E +FK+ A A++ AV ++    GEIPS+KG L"
assert record.rounds[1].alignments[20].hsps[0].sbjct=="RRSMKTRETLETHVKVDLEIDGSGKSSVNTGLGFLDHMLESVARHGLLDLEVEARGDLE---VDDHHTVEDVALTLGEALREALGDKSGIRRMAHAMVPMDDALATVALDLSGRPYTVLELEFD-DAVIGDVKSQNIGHFIESLAVSAAMNIHASV-RGRNDHHKAEALFKALALAIRDAVRVEH---GEIPSTKGKL"
assert record.rounds[1].alignments[20].hsps[0].query_start==3
assert record.rounds[1].alignments[20].hsps[0].query_end==200
assert record.rounds[1].alignments[20].hsps[0].sbjct_start==5
assert record.rounds[1].alignments[20].hsps[0].sbjct_end==194
assert record.rounds[1].alignments[21].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKG"
assert record.rounds[1].alignments[21].hsps[0].match=="RI  + R TKET I L IN+DGTG+  I TGI F DH+L     H  FDL +   GD E   +D HH +EDV I LG  +++    K  I R+G   IPMD+A  T  +D+SGR Y V + +    + +G   TE    FF ++A    + +H     G+N HH  E +FK+   AL  A  IDE K   + S+KG"
assert record.rounds[1].alignments[21].hsps[0].sbjct=="RIFEVMRETKETNIYLKINIDGTGKYKIDTGIPFFDHLLASFAKHGCFDLIVKARGDLE---IDDHHTVEDVGICLGLALNQI--EKRNIFRFGWAIIPMDDARATVAIDLSGRSYCVGNYE-PKREFVGDLATENINHFFESVASYGMLNIHYEVI-GKNEHHKAEALFKAFGVALDLATKIDERK--GVISTKG"
assert record.rounds[1].alignments[21].hsps[0].query_start==3
assert record.rounds[1].alignments[21].hsps[0].query_end==198
assert record.rounds[1].alignments[21].hsps[0].sbjct_start==7
assert record.rounds[1].alignments[21].hsps[0].sbjct_end==193
assert record.rounds[1].alignments[22].hsps[0].query=="RISHITRNTKETQIELSINLDG------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALK"
assert record.rounds[1].alignments[22].hsps[0].match=="R + ++R+T ET I++++++DG                        +    I+TGIGFLDHML  L  H+ + + +   GD   + +D HH  ED  IA+G   ++ LG   G+ R+G    P+DEAL    +D+S RPY V    L   +KLG    EM     ++ A  A ITLH++   G N HH  E  FK+ A A++"
assert record.rounds[1].alignments[22].hsps[0].sbjct=="RAAALSRDTNETSIQIALSIDGGELPQDADPRLLEASSAHASQTSKSQVISINTGIGFLDHMLHALAKHAGWSMALNCKGD---LHIDDHHTAEDCCIAVGTTFAKALGALTGVARFGYAYAPLDEALSRAVVDLSNRPYTVVDLGLK-REKLGELSCEMIPHCLQSFAQAARITLHVDCLRGDNDHHRAESAFKALAVAVR"
assert record.rounds[1].alignments[22].hsps[0].query_start==3
assert record.rounds[1].alignments[22].hsps[0].query_end==180
assert record.rounds[1].alignments[22].hsps[0].sbjct_start==8
assert record.rounds[1].alignments[22].hsps[0].sbjct_end==205
assert record.rounds[1].alignments[23].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[23].hsps[0].match=="M+R ++ITR TKET+IE+ +++D  G+  +ST I F +HML  L  + +    +      + +  D HH++EDVAI LG  I   LG+K GI+R+    IPMD+ALV   LDIS R     + +L  ++ +GG  TE    FF++ A+N+GITLH+++  G NTHHIIE  FK+   AL +A  I ++   EI S+KG++"
assert record.rounds[1].alignments[23].hsps[0].sbjct=="MSRSANITRETKETKIEVLLDIDRKGEVKVSTPIPFFNHMLITLLTYMNSTAIVSAT---DKLPYDDHHIVEDVAITLGLAIKTALGDKRGIKRFSHQIIPMDDALVLVSLDISNRGMAFVNLNLKRSE-IGGLATENVPHFFQSFAYNSGITLHISQLSGYNTHHIIEASFKALGLALYEATRIVDN---EIRSTKGII"
assert record.rounds[1].alignments[23].hsps[0].query_start==1
assert record.rounds[1].alignments[23].hsps[0].query_end==200
assert record.rounds[1].alignments[23].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[23].hsps[0].sbjct_end==193
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.317, 0.149, 0.479]
assert record.ka_params_gap==[0.270, 0.0524, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==23641501
assert record.num_sequences==82258
assert record.num_extends==1047546
assert record.num_good_extends==2481
assert record.num_seqs_better_e==54
assert record.hsps_no_gap==48
assert record.hsps_prelim_gapped==6
assert record.hsps_gapped==56
assert record.query_length==200
assert record.database_length==29652561
assert record.effective_hsp_length==50
assert record.effective_query_length==150
assert record.effective_database_length==25539661
assert record.effective_search_space==3830949150
assert record.effective_search_space_used==3830949150
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.5)
assert record.blast_cutoff==(63,28.8)

handle = open('Blast/bt010')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1348916|gb|G26684|G26684 human STS\nSTS_D11570.\x01gi|1375195|gb|G26945|G26945 human STS SHGC-32699."
assert record.query_letters==285
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==4
assert record.descriptions[0].title=="gi|1348916|gb|G26684|G26684 human STS STS_D11570. >gi|1375195|g..."
assert record.descriptions[0].score==517
assert record.descriptions[0].e==1e-146
assert record.descriptions[1].title=="gi|720683|gb|G03725|G03725 human STS WI-344."
assert record.descriptions[1].score==32
assert record.descriptions[1].e==1.600000
assert record.descriptions[2].title=="gi|4516686|dbj|AU026763.1|AU026763 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[2].score==32
assert record.descriptions[2].e==1.600000
assert record.descriptions[3].title=="gi|6120827|gb|G55508.1|G55508 SHGC-100856 Human Homo sapiens ST..."
assert record.descriptions[3].score==32
assert record.descriptions[3].e==1.600000
assert len(record.alignments)==4
assert record.alignments[0].title==">gi|1348916|gb|G26684|G26684 human STS STS_D11570. >gi|1375195|gb|G26945|G26945 human STS SHGC-32699."
assert record.alignments[0].length==285
assert record.alignments[1].title==">gi|720683|gb|G03725|G03725 human STS WI-344."
assert record.alignments[1].length==246
assert record.alignments[2].title==">gi|4516686|dbj|AU026763.1|AU026763 Rattus norvegicus, OTSUKA clone, OT33.16/752f07, microsatellite sequence, sequence tagged site"
assert record.alignments[2].length==307
assert record.alignments[3].title==">gi|6120827|gb|G55508.1|G55508 SHGC-100856 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[3].length==711
assert record.alignments[0].hsps[0].score==261
assert record.alignments[0].hsps[0].bits==517
assert record.alignments[0].hsps[0].expect==1e-146
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==16
assert record.alignments[1].hsps[0].bits==32.2
assert record.alignments[1].hsps[0].expect==1.6
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==16
assert record.alignments[2].hsps[0].bits==32.2
assert record.alignments[2].hsps[0].expect==1.6
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==16
assert record.alignments[3].hsps[0].bits==32.2
assert record.alignments[3].hsps[0].expect==1.6
assert len(record.alignments[3].hsps)==1
assert record.alignments[0].hsps[0].identities==(285, 285)
assert record.alignments[1].hsps[0].identities==(16, 16)
assert record.alignments[2].hsps[0].identities==(16, 16)
assert record.alignments[3].hsps[0].identities==(18, 19)
assert record.alignments[0].hsps[0].strand==("Plus", "Plus")
assert record.alignments[1].hsps[0].strand==("Plus", "Minus")
assert record.alignments[2].hsps[0].strand==("Plus", "Plus")
assert record.alignments[3].hsps[0].strand==("Plus", "Plus")
assert record.alignments[0].hsps[0].query=="gatccctacccttnccgttggtctctntcgctgactcgaggcacctaacatccattcacacccaacacaggccagcgacttctggggctcagccacagacatggtttgtnactnttgagcttctgttcctagagaatcctagaggcttgattggcccaggctgctgtntgtnctggaggcaaagaatccctacctcctaggggtgaaaggaaatnaaaatggaaagttcttgtagcgcaaggcctgacatgggtagctgctcaataaatgctagtntgttatttc"
assert record.alignments[0].hsps[0].match=="|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
assert record.alignments[0].hsps[0].sbjct=="gatccctacccttnccgttggtctctntcgctgactcgaggcacctaacatccattcacacccaacacaggccagcgacttctggggctcagccacagacatggtttgtnactnttgagcttctgttcctagagaatcctagaggcttgattggcccaggctgctgtntgtnctggaggcaaagaatccctacctcctaggggtgaaaggaaatnaaaatggaaagttcttgtagcgcaaggcctgacatgggtagctgctcaataaatgctagtntgttatttc"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==285
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==285
assert record.alignments[1].hsps[0].query=="ctcaataaatgctagt"
assert record.alignments[1].hsps[0].match=="||||||||||||||||"
assert record.alignments[1].hsps[0].sbjct=="ctcaataaatgctagt"
assert record.alignments[1].hsps[0].query_start==260
assert record.alignments[1].hsps[0].query_end==275
assert record.alignments[1].hsps[0].sbjct_start==178
assert record.alignments[1].hsps[0].sbjct_end==163
assert record.alignments[2].hsps[0].query=="ggaaagttcttgtagc"
assert record.alignments[2].hsps[0].match=="||||||||||||||||"
assert record.alignments[2].hsps[0].sbjct=="ggaaagttcttgtagc"
assert record.alignments[2].hsps[0].query_start==221
assert record.alignments[2].hsps[0].query_end==236
assert record.alignments[2].hsps[0].sbjct_start==32
assert record.alignments[2].hsps[0].sbjct_end==47
assert record.alignments[3].hsps[0].query=="gaaatnaaaatggaaagtt"
assert record.alignments[3].hsps[0].match=="||||| |||||||||||||"
assert record.alignments[3].hsps[0].sbjct=="gaaataaaaatggaaagtt"
assert record.alignments[3].hsps[0].query_start==210
assert record.alignments[3].hsps[0].query_end==228
assert record.alignments[3].hsps[0].sbjct_start==588
assert record.alignments[3].hsps[0].sbjct_end==606
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Nov 26, 1999  5:52 PM',)]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==3835
assert record.num_sequences==87792
assert record.num_extends==3835
assert record.num_good_extends==1155
assert record.num_seqs_better_e==4
assert record.query_length==285
assert record.database_length==31998854
assert record.effective_hsp_length==17
assert record.effective_query_length==268
assert record.effective_database_length==30506390
assert record.effective_search_space==8175712520
assert record.effective_search_space_used==8175712520
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(25,49.6)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(15,30.2)

handle = open('Blast/bt011')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1348400|gb|G26168|G26168 human STS\nEST47274.\x01gi|1375380|gb|G27130|G27130 human STS SHGC-32751."
assert record.query_letters==434
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==19
assert record.descriptions[0].title=="gi|1348400|gb|G26168|G26168 human STS EST47274. >gi|1375380|gb|..."
assert record.descriptions[0].score==718
assert record.descriptions[0].e==0.000000
assert record.descriptions[1].title=="gi|6121436|gb|G56117.1|G56117 SHGC-101583 Human Homo sapiens ST..."
assert record.descriptions[1].score==34
assert record.descriptions[1].e==0.650000
assert record.descriptions[2].title=="gi|4632200|dbj|AU047565.1|AU047565 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[2].score==34
assert record.descriptions[2].e==0.650000
assert record.descriptions[3].title=="gi|3249175|gb|G38401|G38401 SHGC-57345 Human Homo sapiens STS g..."
assert record.descriptions[3].score==34
assert record.descriptions[3].e==0.650000
assert record.descriptions[4].title=="gi|605557|gb|L31312|HUMUT937B Human STS UT937, 3' primer bind."
assert record.descriptions[4].score==34
assert record.descriptions[4].e==0.650000
assert record.descriptions[5].title=="gi|720383|gb|G03425|G03425 human STS WI-5868."
assert record.descriptions[5].score==34
assert record.descriptions[5].e==0.650000
assert record.descriptions[6].title=="gi|859592|gb|G06347|G06347 human STS WI-7005."
assert record.descriptions[6].score==32
assert record.descriptions[6].e==2.600000
assert record.descriptions[7].title=="gi|6121347|gb|G56178.1|G56178 SHGC-101470 Human Homo sapiens ST..."
assert record.descriptions[7].score==32
assert record.descriptions[7].e==2.600000
assert record.descriptions[8].title=="gi|1233216|emb|Z51916|HSA082WB5 H.sapiens (D1S2890) DNA segment..."
assert record.descriptions[8].score==32
assert record.descriptions[8].e==2.600000
assert record.descriptions[9].title=="gi|1232198|emb|Z50898|HS038XD8 H.sapiens (D18S1106) DNA segment..."
assert record.descriptions[9].score==32
assert record.descriptions[9].e==2.600000
assert record.descriptions[10].title=="gi|1348143|gb|G25911|G25911 human STS EST349382."
assert record.descriptions[10].score==32
assert record.descriptions[10].e==2.600000
assert record.descriptions[11].title=="gi|6122805|gb|G57486.1|G57486 SHGC-103345 Human Homo sapiens ST..."
assert record.descriptions[11].score==32
assert record.descriptions[11].e==2.600000
assert record.descriptions[12].title=="gi|1396897|gb|G28178|G28178 human STS SHGC-34170."
assert record.descriptions[12].score==32
assert record.descriptions[12].e==2.600000
assert record.descriptions[13].title=="gi|3893806|emb|AL034156|HSPE59A01 H.sapiens flow-sorted chromos..."
assert record.descriptions[13].score==32
assert record.descriptions[13].e==2.600000
assert record.descriptions[14].title=="gi|1161890|gb|G16001|G16001 human STS CHLC.GCT8D08.P11278 clone..."
assert record.descriptions[14].score==32
assert record.descriptions[14].e==2.600000
assert record.descriptions[15].title=="gi|1017612|gb|G11520|G11520 human STS SHGC-14676."
assert record.descriptions[15].score==32
assert record.descriptions[15].e==2.600000
assert record.descriptions[16].title=="gi|1130137|gb|G14398|G14398 human STS SHGC-9310 clone pG-5191."
assert record.descriptions[16].score==32
assert record.descriptions[16].e==2.600000
assert record.descriptions[17].title=="gi|5224295|gb|G52968.1|G52968 SHGC-86325 Human Homo sapiens STS..."
assert record.descriptions[17].score==32
assert record.descriptions[17].e==2.600000
assert record.descriptions[18].title=="gi|6123581|gb|G58262.1|G58262 SHGC-104352 Human Homo sapiens ST..."
assert record.descriptions[18].score==32
assert record.descriptions[18].e==2.600000
assert len(record.alignments)==0
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Nov 26, 1999  5:52 PM',)]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==8762
assert record.num_sequences==87792
assert record.num_extends==8762
assert record.num_good_extends==2655
assert record.num_seqs_better_e==27
assert record.query_length==434
assert record.database_length==31998854
assert record.effective_hsp_length==17
assert record.effective_query_length==417
assert record.effective_database_length==30506390
assert record.effective_search_space==12721164630
assert record.effective_search_space_used==12721164630
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(25,49.6)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(16,32.2)

handle = open('Blast/bt012')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1347201|gb|G24969|G24969 human STS\nEST21946.\x01gi|1375315|gb|G27065|G27065 human STS SHGC-31731."
assert record.query_letters==331
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==0
assert len(record.alignments)==9
assert record.alignments[0].title==">gi|1347201|gb|G24969|G24969 human STS EST21946. >gi|1375315|gb|G27065|G27065 human STS SHGC-31731."
assert record.alignments[0].length==331
assert record.alignments[1].title==">gi|2665277|emb|AL010115|HSPE77H4 H.sapiens flow-sorted chromosome 1 HindIII fragment, SC1pE77H4, sequence tagged site [Homo sapiens]"
assert record.alignments[1].length==554
assert record.alignments[2].title==">gi|6120731|gb|G55412.1|G55412 SHGC-100745 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[2].length==652
assert record.alignments[3].title==">gi|4493602|gb|G47248.1|G47248 Z17392_1 Zebrafish AB Danio rerio STS genomic clone Z17392 5', sequence tagged site"
assert record.alignments[3].length==454
assert record.alignments[4].title==">gi|1342139|gb|G21813|G21813 human STS WI-12408."
assert record.alignments[4].length==418
assert record.alignments[5].title==">gi|1235411|emb|Z53965|HSC009WH1 H.sapiens (D2S2321) DNA segment containing (CA) repeat; clone AFMc009wh1; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[5].length==382
assert record.alignments[6].title==">gi|4757440|gb|G49267.1|G49267 stbK343C1_96809 chromosome 22 genomic clone Homo sapiens STS genomic clone 343C1, sequence tagged site"
assert record.alignments[6].length==360
assert record.alignments[7].title==">gi|719782|gb|G02824|G02824 human STS WI-1312."
assert record.alignments[7].length==349
assert record.alignments[8].title==">gi|939357|gb|G08807|G08807 human STS CHLC.GATA70E11.P18111 clone GATA70E11"
assert record.alignments[8].length==643
assert record.alignments[0].hsps[0].score==331
assert record.alignments[0].hsps[0].bits==656
assert record.alignments[0].hsps[0].expect==0.0
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==17
assert record.alignments[1].hsps[0].bits==34.2
assert record.alignments[1].hsps[0].expect==0.49
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==16
assert record.alignments[2].hsps[0].bits==32.2
assert record.alignments[2].hsps[0].expect==1.9
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==16
assert record.alignments[3].hsps[0].bits==32.2
assert record.alignments[3].hsps[0].expect==1.9
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==16
assert record.alignments[4].hsps[0].bits==32.2
assert record.alignments[4].hsps[0].expect==1.9
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==16
assert record.alignments[5].hsps[0].bits==32.2
assert record.alignments[5].hsps[0].expect==1.9
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==16
assert record.alignments[6].hsps[0].bits==32.2
assert record.alignments[6].hsps[0].expect==1.9
assert len(record.alignments[6].hsps)==1
assert record.alignments[7].hsps[0].score==16
assert record.alignments[7].hsps[0].bits==32.2
assert record.alignments[7].hsps[0].expect==1.9
assert len(record.alignments[7].hsps)==1
assert record.alignments[8].hsps[0].score==16
assert record.alignments[8].hsps[0].bits==32.2
assert record.alignments[8].hsps[0].expect==1.9
assert len(record.alignments[8].hsps)==1
assert record.alignments[0].hsps[0].identities==(331, 331)
assert record.alignments[1].hsps[0].identities==(17, 17)
assert record.alignments[2].hsps[0].identities==(16, 16)
assert record.alignments[3].hsps[0].identities==(16, 16)
assert record.alignments[4].hsps[0].identities==(16, 16)
assert record.alignments[5].hsps[0].identities==(16, 16)
assert record.alignments[6].hsps[0].identities==(22, 24)
assert record.alignments[7].hsps[0].identities==(16, 16)
assert record.alignments[8].hsps[0].identities==(16, 16)
assert record.alignments[0].hsps[0].strand==("Plus", "Plus")
assert record.alignments[1].hsps[0].strand==("Plus", "Minus")
assert record.alignments[2].hsps[0].strand==("Plus", "Minus")
assert record.alignments[3].hsps[0].strand==("Plus", "Plus")
assert record.alignments[4].hsps[0].strand==("Plus", "Minus")
assert record.alignments[5].hsps[0].strand==("Plus", "Plus")
assert record.alignments[6].hsps[0].strand==("Plus", "Minus")
assert record.alignments[7].hsps[0].strand==("Plus", "Plus")
assert record.alignments[8].hsps[0].strand==("Plus", "Minus")
assert record.alignments[0].hsps[0].query=="cctccaccctctcatgagcaacaggatatgtgaaagtacttgcagccagaagcaaaaccacaatcctcgggtgctagatggagctccccaaggagcagagaggaaaaggcaggaggagagggccaggcagcagggatggagactaagtttggcccaaggctgcccgcaagcactgatgccatcatgccctctggtaggtgtctatttctgtctgaaccagaaatacaccaagctccacacatgggggctttgctggcttcgacatcactggttcaactatgtcactgctttgttatatttagtgctccagaacctcaggttccttcagatt"
assert record.alignments[0].hsps[0].match=="|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
assert record.alignments[0].hsps[0].sbjct=="cctccaccctctcatgagcaacaggatatgtgaaagtacttgcagccagaagcaaaaccacaatcctcgggtgctagatggagctccccaaggagcagagaggaaaaggcaggaggagagggccaggcagcagggatggagactaagtttggcccaaggctgcccgcaagcactgatgccatcatgccctctggtaggtgtctatttctgtctgaaccagaaatacaccaagctccacacatgggggctttgctggcttcgacatcactggttcaactatgtcactgctttgttatatttagtgctccagaacctcaggttccttcagatt"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==331
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==331
assert record.alignments[1].hsps[0].query=="ccaggcagcagggatgg"
assert record.alignments[1].hsps[0].match=="|||||||||||||||||"
assert record.alignments[1].hsps[0].sbjct=="ccaggcagcagggatgg"
assert record.alignments[1].hsps[0].query_start==123
assert record.alignments[1].hsps[0].query_end==139
assert record.alignments[1].hsps[0].sbjct_start==434
assert record.alignments[1].hsps[0].sbjct_end==418
assert record.alignments[2].hsps[0].query=="agaggaaaaggcagga"
assert record.alignments[2].hsps[0].match=="||||||||||||||||"
assert record.alignments[2].hsps[0].sbjct=="agaggaaaaggcagga"
assert record.alignments[2].hsps[0].query_start==99
assert record.alignments[2].hsps[0].query_end==114
assert record.alignments[2].hsps[0].sbjct_start==431
assert record.alignments[2].hsps[0].sbjct_end==416
assert record.alignments[3].hsps[0].query=="agaagcaaaaccacaa"
assert record.alignments[3].hsps[0].match=="||||||||||||||||"
assert record.alignments[3].hsps[0].sbjct=="agaagcaaaaccacaa"
assert record.alignments[3].hsps[0].query_start==48
assert record.alignments[3].hsps[0].query_end==63
assert record.alignments[3].hsps[0].sbjct_start==434
assert record.alignments[3].hsps[0].sbjct_end==449
assert record.alignments[4].hsps[0].query=="cagaagcaaaaccaca"
assert record.alignments[4].hsps[0].match=="||||||||||||||||"
assert record.alignments[4].hsps[0].sbjct=="cagaagcaaaaccaca"
assert record.alignments[4].hsps[0].query_start==47
assert record.alignments[4].hsps[0].query_end==62
assert record.alignments[4].hsps[0].sbjct_start==193
assert record.alignments[4].hsps[0].sbjct_end==178
assert record.alignments[5].hsps[0].query=="agagaggaaaaggcag"
assert record.alignments[5].hsps[0].match=="||||||||||||||||"
assert record.alignments[5].hsps[0].sbjct=="agagaggaaaaggcag"
assert record.alignments[5].hsps[0].query_start==97
assert record.alignments[5].hsps[0].query_end==112
assert record.alignments[5].hsps[0].sbjct_start==107
assert record.alignments[5].hsps[0].sbjct_end==122
assert record.alignments[6].hsps[0].query=="ggaggagagggccaggcagcaggg"
assert record.alignments[6].hsps[0].match=="||||||||||| | ||||||||||"
assert record.alignments[6].hsps[0].sbjct=="ggaggagaggggctggcagcaggg"
assert record.alignments[6].hsps[0].query_start==112
assert record.alignments[6].hsps[0].query_end==135
assert record.alignments[6].hsps[0].sbjct_start==287
assert record.alignments[6].hsps[0].sbjct_end==264
assert record.alignments[7].hsps[0].query=="tgctttgttatattta"
assert record.alignments[7].hsps[0].match=="||||||||||||||||"
assert record.alignments[7].hsps[0].sbjct=="tgctttgttatattta"
assert record.alignments[7].hsps[0].query_start==286
assert record.alignments[7].hsps[0].query_end==301
assert record.alignments[7].hsps[0].sbjct_start==111
assert record.alignments[7].hsps[0].sbjct_end==126
assert record.alignments[8].hsps[0].query=="gaggaaaaggcaggag"
assert record.alignments[8].hsps[0].match=="||||||||||||||||"
assert record.alignments[8].hsps[0].sbjct=="gaggaaaaggcaggag"
assert record.alignments[8].hsps[0].query_start==100
assert record.alignments[8].hsps[0].query_end==115
assert record.alignments[8].hsps[0].sbjct_start==482
assert record.alignments[8].hsps[0].sbjct_end==467
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Nov 26, 1999  5:52 PM',)]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==6844
assert record.num_sequences==87792
assert record.num_extends==6844
assert record.num_good_extends==1887
assert record.num_seqs_better_e==14
assert record.query_length==331
assert record.database_length==31998854
assert record.effective_hsp_length==17
assert record.effective_query_length==314
assert record.effective_database_length==30506390
assert record.effective_search_space==9579006460
assert record.effective_search_space_used==9579006460
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(25,49.6)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(15,30.2)

handle = open('Blast/bt013')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|859351|gb|G06106|G06106 human STS WI-6344."
assert record.query_letters==183
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==6
assert record.descriptions[0].title=="gi|859351|gb|G06106|G06106 human STS WI-6344."
assert record.descriptions[0].score==327
assert record.descriptions[0].e==1e-89
assert record.descriptions[1].title=="gi|1341350|gb|G21024|G21024 human STS WI-30979."
assert record.descriptions[1].score==32
assert record.descriptions[1].e==1.000000
assert record.descriptions[2].title=="gi|6126285|gb|G60966.1|G60966 SHGC-84377 Human Homo sapiens STS..."
assert record.descriptions[2].score==30
assert record.descriptions[2].e==4.100000
assert record.descriptions[3].title=="gi|1340656|gb|G20319|G20319 human STS A005L39."
assert record.descriptions[3].score==30
assert record.descriptions[3].e==4.100000
assert record.descriptions[4].title=="gi|5222421|gb|G51244.1|G51244 SHGC-80725 Human Homo sapiens STS..."
assert record.descriptions[4].score==30
assert record.descriptions[4].e==4.100000
assert record.descriptions[5].title=="gi|860526|gb|G07281|G07281 human STS WI-9430."
assert record.descriptions[5].score==30
assert record.descriptions[5].e==4.100000
assert len(record.alignments)==0
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Nov 26, 1999  5:52 PM',)]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==1902
assert record.num_sequences==87792
assert record.num_extends==1902
assert record.num_good_extends==481
assert record.num_seqs_better_e==8
assert record.query_length==183
assert record.database_length==31998854
assert record.effective_hsp_length==16
assert record.effective_query_length==167
assert record.effective_database_length==30594182
assert record.effective_search_space==5109228394
assert record.effective_search_space_used==5109228394
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(25,49.6)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(15,30.2)

handle = open('Blast/bt014')
record = parser.parse(handle)
assert record.application=="BLASTX"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1347369|gb|G25137|G25137 human STS EST48004."
assert record.query_letters==556
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==4
assert record.descriptions[0].title=="gi|1731448|sp|P54103|ZRF1_MOUSE ZUOTIN RELATED FACTOR-1"
assert record.descriptions[0].score==87
assert record.descriptions[0].e==3.0000000000000001e-17
assert record.descriptions[1].title=="gi|465911|sp|P34454|YMA9_CAEEL HYPOTHETICAL 31.6 KD PROTEIN F54..."
assert record.descriptions[1].score==42
assert record.descriptions[1].e==0.001000
assert record.descriptions[2].title=="gi|2494160|sp|Q61712|MTJ1_MOUSE DNAJ PROTEIN HOMOLOG MTJ1"
assert record.descriptions[2].score==37
assert record.descriptions[2].e==0.033000
assert record.descriptions[3].title=="gi|1730688|sp|P53745|YN8X_YEAST HYPOTHETICAL 68.1 KD PROTEIN IN..."
assert record.descriptions[3].score==29
assert record.descriptions[3].e==7.400000
assert len(record.alignments)==4
assert record.alignments[0].title==">gi|1731448|sp|P54103|ZRF1_MOUSE ZUOTIN RELATED FACTOR-1"
assert record.alignments[0].length==514
assert record.alignments[1].title==">gi|465911|sp|P34454|YMA9_CAEEL HYPOTHETICAL 31.6 KD PROTEIN F54F2.9 IN CHROMOSOME III"
assert record.alignments[1].length==275
assert record.alignments[2].title==">gi|2494160|sp|Q61712|MTJ1_MOUSE DNAJ PROTEIN HOMOLOG MTJ1"
assert record.alignments[2].length==552
assert record.alignments[3].title==">gi|1730688|sp|P53745|YN8X_YEAST HYPOTHETICAL 68.1 KD PROTEIN IN BIO3-FRE4 INTERGENIC REGION"
assert record.alignments[3].length==580
assert record.alignments[0].hsps[0].score==211
assert record.alignments[0].hsps[0].bits==86.6
assert record.alignments[0].hsps[0].expect==3e-17
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==96
assert record.alignments[1].hsps[0].bits==41.8
assert record.alignments[1].hsps[0].expect==0.001
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==83
assert record.alignments[2].hsps[0].bits==36.7
assert record.alignments[2].hsps[0].expect==0.033
assert record.alignments[2].hsps[1].score==69
assert record.alignments[2].hsps[1].bits==31.3
assert record.alignments[2].hsps[1].expect==1.5
assert len(record.alignments[2].hsps)==2
assert record.alignments[3].hsps[0].score==63
assert record.alignments[3].hsps[0].bits==29.0
assert record.alignments[3].hsps[0].expect==7.4
assert len(record.alignments[3].hsps)==1
assert record.alignments[0].hsps[0].identities==(41, 47)
assert record.alignments[0].hsps[0].positives==(44, 47)
assert record.alignments[1].hsps[0].identities==(30, 122)
assert record.alignments[1].hsps[0].positives==(54, 122)
assert record.alignments[2].hsps[0].identities==(17, 36)
assert record.alignments[2].hsps[0].positives==(19, 36)
assert record.alignments[2].hsps[1].identities==(18, 50)
assert record.alignments[2].hsps[1].positives==(26, 50)
assert record.alignments[3].hsps[0].identities==(27, 99)
assert record.alignments[3].hsps[0].positives==(41, 99)
assert record.alignments[0].hsps[0].frame==("+1", )
assert record.alignments[1].hsps[0].frame==("+1", )
assert record.alignments[2].hsps[0].frame==("+1", )
assert record.alignments[2].hsps[1].frame==("+1", )
assert record.alignments[3].hsps[0].frame==("+1", )
assert record.alignments[0].hsps[0].query=="DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAK"
assert record.alignments[0].hsps[0].match=="DLQLLIKAVNLFPAG NSRW+VIANYMNIHSSSGVKRTAKDVI + +"
assert record.alignments[0].hsps[0].sbjct=="DLQLLIKAVNLFPAGRNSRWDVIANYMNIHSSSGVKRTAKDVISEVR"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==141
assert record.alignments[0].hsps[0].sbjct_start==458
assert record.alignments[0].hsps[0].sbjct_end==504
assert record.alignments[1].hsps[0].query=="FPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFXGPYTDFTPXTTEXQKLXEQALNTYPVNTXERWXXIAVAVPGRXKE"
assert record.alignments[1].hsps[0].match=="+PAGT +RWE +   +N        R+A+DVI  A  ++++   +++D  K      ++   V  ++++                 +   QK  E AL  YP  T ERW  I+  +  + K+"
assert record.alignments[1].hsps[0].sbjct=="YPAGTPNRWEQMGRVLN--------RSAEDVIAMAGKMKQM---KQEDYTKLLMTTIQQSVPVEEKSED---------------DWSQAEQKAFETALQKYPKGTDERWERISEEIGSKTKK"
assert record.alignments[1].hsps[0].query_start==34
assert record.alignments[1].hsps[0].query_end==399
assert record.alignments[1].hsps[0].sbjct_start==159
assert record.alignments[1].hsps[0].sbjct_end==254
assert record.alignments[2].hsps[0].query=="TTEXQKLXEQALNTYPVNTXERWXXIAVAVPGRXKE"
assert record.alignments[2].hsps[0].match=="T   QKL E AL  YP    +RW  IA  VP + KE"
assert record.alignments[2].hsps[0].sbjct=="TQSQQKLLELALQQYPKGASDRWDKIAKCVPSKSKE"
assert record.alignments[2].hsps[0].query_start==292
assert record.alignments[2].hsps[0].query_end==399
assert record.alignments[2].hsps[0].sbjct_start==496
assert record.alignments[2].hsps[0].sbjct_end==531
assert record.alignments[2].hsps[1].query=="DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQ"
assert record.alignments[2].hsps[1].match=="DL  L +++  FP GT  RW+ IA+ +         R+  DV  KAK L+"
assert record.alignments[2].hsps[1].sbjct=="DLSQLTRSMVKFPGGTPGRWDKIAHELG--------RSVTDVTTKAKELK"
assert record.alignments[2].hsps[1].query_start==1
assert record.alignments[2].hsps[1].query_end==150
assert record.alignments[2].hsps[1].sbjct_start==332
assert record.alignments[2].hsps[1].sbjct_end==373
assert record.alignments[3].hsps[0].query=="SRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFXGPYTDFTPXTTEXQKLXEQALNTYPVNT"
assert record.alignments[3].hsps[0].match=="+RW+   +Y        V R+ KDV   ++SL  LD +QK     +A       H +          E    PY +FT   +      EQ+ N +PV+T"
assert record.alignments[3].hsps[0].sbjct=="NRWKSFISY--------VTRSRKDVKTVSRSLSNLDLYQKCSKEIRADQDISLLHSI----------ETKLFPYINFTALNS------EQSHNFWPVHT"
assert record.alignments[3].hsps[0].query_start==52
assert record.alignments[3].hsps[0].query_end==348
assert record.alignments[3].hsps[0].sbjct_start==75
assert record.alignments[3].hsps[0].sbjct_end==149
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.318, 0.135, 0.401]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==23174157
assert record.num_sequences==82258
assert record.num_extends==387821
assert record.num_good_extends==980
assert record.num_seqs_better_e==8
assert record.hsps_no_gap==3
assert record.hsps_prelim_gapped==1
assert record.hsps_gapped==7
assert record.query_length==185
assert record.database_length==29652561
assert record.effective_hsp_length==49
assert record.effective_query_length==135
assert record.effective_database_length==25621919
assert record.effective_search_space==3458959065
assert record.effective_search_space_used==3458959065
assert record.frameshift==('50,','0.1')
assert record.threshold==12
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(62,28.6)

handle = open('Blast/bt015')
record = parser.parse(handle)
assert record.application=="BLASTX"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1347782|gb|G25550|G25550 human STS\nEST47652.\x01gi|1592937|gb|G29386|G29386 human STS SHGC-32770"
assert record.query_letters==379
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==0
assert len(record.alignments)==0
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 15, 1999  2:55 PM',)]
assert record.ka_params==[0.318, 0.135, 0.401]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==14686001
assert record.num_sequences==82258
assert record.num_extends==235383
assert record.num_good_extends==396
assert record.num_seqs_better_e==0
assert record.hsps_no_gap==0
assert record.hsps_prelim_gapped==0
assert record.hsps_gapped==0
assert record.query_length==126
assert record.database_length==29652561
assert record.effective_hsp_length==48
assert record.effective_query_length==77
assert record.effective_database_length==25704177
assert record.effective_search_space==1979221629
assert record.effective_search_space_used==1979221629
assert record.frameshift==('50,','0.1')
assert record.threshold==12
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(60,27.8)

handle = open('Blast/bt016')
record = parser.parse(handle)
assert record.application=="TBLASTN"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|729325|sp|P39483|DHG2_BACME GLUCOSE 1-DEHYDROGENASE II\n(GLCDH-II)"
assert record.query_letters==261
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==7
assert record.descriptions[0].title=="gi|3820341|emb|AJ229891|KLAJ9891 Kluyveromyces lactis DNA fragm..."
assert record.descriptions[0].score==47
assert record.descriptions[0].e==0.000010
assert record.descriptions[1].title=="gi|1375419|gb|G27169|G27169 human STS SHGC-31983."
assert record.descriptions[1].score==43
assert record.descriptions[1].e==0.000100
assert record.descriptions[2].title=="gi|3819804|emb|AJ230012|KLAJ0012 Kluyveromyces lactis DNA fragm..."
assert record.descriptions[2].score==39
assert record.descriptions[2].e==0.002000
assert record.descriptions[3].title=="gi|1375215|gb|G26965|G26965 human STS SHGC-31083."
assert record.descriptions[3].score==31
assert record.descriptions[3].e==0.730000
assert record.descriptions[4].title=="gi|177714|gb|L09988|HUM4STS889 Human Chromosome 4 (clone p4-109..."
assert record.descriptions[4].score==29
assert record.descriptions[4].e==2.200000
assert record.descriptions[5].title=="gi|5714409|gb|AF106665.1|AF106665 Mus musculus chromosome 6 clo..."
assert record.descriptions[5].score==29
assert record.descriptions[5].e==2.200000
assert record.descriptions[6].title=="gi|1341648|gb|G21322|G21322 human STS WI-12250."
assert record.descriptions[6].score==29
assert record.descriptions[6].e==3.700000
assert len(record.alignments)==7
assert record.alignments[0].title==">gi|3820341|emb|AJ229891|KLAJ9891 Kluyveromyces lactis DNA fragment for sequence tagged site, clone okam5d07r [Kluyveromyces lactis]"
assert record.alignments[0].length==230
assert record.alignments[1].title==">gi|1375419|gb|G27169|G27169 human STS SHGC-31983."
assert record.alignments[1].length==594
assert record.alignments[2].title==">gi|3819804|emb|AJ230012|KLAJ0012 Kluyveromyces lactis DNA fragment for sequence tagged site, clone okam6d01d [Kluyveromyces lactis]"
assert record.alignments[2].length==199
assert record.alignments[3].title==">gi|1375215|gb|G26965|G26965 human STS SHGC-31083."
assert record.alignments[3].length==268
assert record.alignments[4].title==">gi|177714|gb|L09988|HUM4STS889 Human Chromosome 4 (clone p4-1095) STS4-889."
assert record.alignments[4].length==412
assert record.alignments[5].title==">gi|5714409|gb|AF106665.1|AF106665 Mus musculus chromosome 6 clone D6wum9 map between Nkrp1 and Prp strain C57BL/6J, sequence tagged site"
assert record.alignments[5].length==299
assert record.alignments[6].title==">gi|1341648|gb|G21322|G21322 human STS WI-12250."
assert record.alignments[6].length==586
assert record.alignments[0].hsps[0].score==109
assert record.alignments[0].hsps[0].bits==46.9
assert record.alignments[0].hsps[0].expect==1e-05
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==100
assert record.alignments[1].hsps[0].bits==43.4
assert record.alignments[1].hsps[0].expect==1e-04
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==90
assert record.alignments[2].hsps[0].bits==39.5
assert record.alignments[2].hsps[0].expect==0.002
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==68
assert record.alignments[3].hsps[0].bits==30.9
assert record.alignments[3].hsps[0].expect==0.73
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==64
assert record.alignments[4].hsps[0].bits==29.3
assert record.alignments[4].hsps[0].expect==2.2
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==64
assert record.alignments[5].hsps[0].bits==29.3
assert record.alignments[5].hsps[0].expect==2.2
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==62
assert record.alignments[6].hsps[0].bits==28.6
assert record.alignments[6].hsps[0].expect==3.7
assert len(record.alignments[6].hsps)==1
assert record.alignments[0].hsps[0].identities==(25, 72)
assert record.alignments[0].hsps[0].positives==(44, 72)
assert record.alignments[0].hsps[0].gaps==(3, 72)
assert record.alignments[1].hsps[0].identities==(21, 73)
assert record.alignments[1].hsps[0].positives==(34, 73)
assert record.alignments[2].hsps[0].identities==(18, 49)
assert record.alignments[2].hsps[0].positives==(26, 49)
assert record.alignments[3].hsps[0].identities==(12, 37)
assert record.alignments[3].hsps[0].positives==(19, 37)
assert record.alignments[4].hsps[0].identities==(14, 34)
assert record.alignments[4].hsps[0].positives==(22, 34)
assert record.alignments[5].hsps[0].identities==(17, 55)
assert record.alignments[5].hsps[0].positives==(32, 55)
assert record.alignments[5].hsps[0].gaps==(2, 55)
assert record.alignments[6].hsps[0].identities==(16, 39)
assert record.alignments[6].hsps[0].positives==(20, 39)
assert record.alignments[6].hsps[0].gaps==(1, 39)
assert record.alignments[0].hsps[0].frame==("+1", )
assert record.alignments[1].hsps[0].frame==("-1", )
assert record.alignments[2].hsps[0].frame==("-1", )
assert record.alignments[3].hsps[0].frame==("-1", )
assert record.alignments[4].hsps[0].frame==("-1", )
assert record.alignments[5].hsps[0].frame==("-2", )
assert record.alignments[6].hsps[0].frame==("+1", )
assert record.alignments[0].hsps[0].query=="NWNQVIDTNLTGAFLGSREAIKYFVEN---DIKGNVINMSSVHEMIPWPLFVHYAASKGGMKLMTETLALEYAPK"
assert record.alignments[0].hsps[0].match=="+W QVIDTN+ G F   + A+     +   D +  V+N+S+V+ ++  P    Y A+K  +  +T+++ALEYA +"
assert record.alignments[0].hsps[0].sbjct=="SWRQVIDTNINGTFYTLKYALPLMESSSSPDSEAAVVNLSAVNGLVGIPGISPYTATKHAVIGITQSVALEYAER"
assert record.alignments[0].hsps[0].query_start==108
assert record.alignments[0].hsps[0].query_end==179
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==225
assert record.alignments[1].hsps[0].query=="APKGIRVNNIGPGAIDTPINAEKFADPEQRADVESMIPMGYIGKPEEIASVAAFLASSQASYVTGITLFADGG"
assert record.alignments[1].hsps[0].match=="AP   RVN + P  + T +    + DP +   + +  P+G   + E +     FL S ++   TG TL  +GG"
assert record.alignments[1].hsps[0].sbjct=="APXRXRVNAVXPXVVMTSMGQATWXDPXKAXTMLNRXPLGXFAEVEHVVKAILFLLSDRSGMTTGSTLPVEGG"
assert record.alignments[1].hsps[0].query_start==177
assert record.alignments[1].hsps[0].query_end==249
assert record.alignments[1].hsps[0].sbjct_start==312
assert record.alignments[1].hsps[0].sbjct_end==94
assert record.alignments[2].hsps[0].query=="FADPEQRADVESMIPMGYIGKPEEIASVAAFLASSQASYVTGITLFADG"
assert record.alignments[2].hsps[0].match=="F D + +    S+IPMG  G P+E+     + AS  ++Y TG  L  DG"
assert record.alignments[2].hsps[0].sbjct=="FVDEDLKNKWHSLIPMGREGLPQELVGAYLYFASDASTYTTGSDLLVDG"
assert record.alignments[2].hsps[0].query_start==200
assert record.alignments[2].hsps[0].query_end==248
assert record.alignments[2].hsps[0].sbjct_start==157
assert record.alignments[2].hsps[0].sbjct_end==11
assert record.alignments[3].hsps[0].query=="PMGYIGKPEEIASVAAFLASSQASYVTGITLFADGGM"
assert record.alignments[3].hsps[0].match=="PMG  G PE++  V A      +  +TG ++   GG+"
assert record.alignments[3].hsps[0].sbjct=="PMGXXGDPEDVXDVXAXXXXEXSGXITGTSVEVTGGL"
assert record.alignments[3].hsps[0].query_start==214
assert record.alignments[3].hsps[0].query_end==250
assert record.alignments[3].hsps[0].sbjct_start==268
assert record.alignments[3].hsps[0].sbjct_end==158
assert record.alignments[4].hsps[0].query=="DKVVVVTGGSKGLGRAMAVRFGQEQSKVVVNYRS"
assert record.alignments[4].hsps[0].match=="DKV  V GGS+G+GRA+A    ++  ++ V  R+"
assert record.alignments[4].hsps[0].sbjct=="DKVCAVFGGSRGIGRAVAQLMARKGYRLAVIARN"
assert record.alignments[4].hsps[0].query_start==7
assert record.alignments[4].hsps[0].query_end==40
assert record.alignments[4].hsps[0].sbjct_start==316
assert record.alignments[4].hsps[0].sbjct_end==215
assert record.alignments[5].hsps[0].query=="NMSSVHEMIPWPLFVHYAASKGGMKLMTETL--ALEYAPKGIRVNNIGPGAIDTPIN"
assert record.alignments[5].hsps[0].match=="++S    +I +P F+    S  G  L+  +L  A+ + P GI V+++GP ++ T +N"
assert record.alignments[5].hsps[0].sbjct=="SLSPTQYLIMFPSFLPCPLSHPGPFLLPSSLVIAVFFLPNGIEVSSLGPFSLRTLLN"
assert record.alignments[5].hsps[0].query_start==142
assert record.alignments[5].hsps[0].query_end==196
assert record.alignments[5].hsps[0].sbjct_start==172
assert record.alignments[5].hsps[0].sbjct_end==2
assert record.alignments[6].hsps[0].query=="PVPSHELSLENWNQ-VIDTNLTGAFLGSREAIKYFVENDI"
assert record.alignments[6].hsps[0].match=="PVP  ELS  +W+Q  + T+ T  F  S     YF  N I"
assert record.alignments[6].hsps[0].sbjct=="PVPMQELSKVHWSQFFLTTSPTMTFFFSHYLANYFFRNSI"
assert record.alignments[6].hsps[0].query_start==98
assert record.alignments[6].hsps[0].query_end==136
assert record.alignments[6].hsps[0].sbjct_start==220
assert record.alignments[6].hsps[0].sbjct_end==339
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Nov 26, 1999  5:52 PM',)]
assert record.ka_params==[0.315, 0.134, 0.378]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==6165875
assert record.num_sequences==87792
assert record.num_extends==55665
assert record.num_good_extends==148
assert record.num_seqs_better_e==14
assert record.hsps_no_gap==5
assert record.hsps_prelim_gapped==2
assert record.hsps_gapped==7
assert record.query_length==261
assert record.database_length==10666284
assert record.effective_hsp_length==50
assert record.effective_query_length==211
assert record.effective_database_length==6276684
assert record.effective_search_space==1324380324
assert record.effective_search_space_used==1324380324
assert record.frameshift==('50,','0.1')
assert record.threshold==13
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(42,22.0)
assert record.blast_cutoff==(58,27.0)

handle = open('Blast/bt017')
record = parser.parse(handle)
assert record.application=="TBLASTN"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|127420|sp|P19888|MTBA_BACAR MODIFICATION METHYLASE BANI\n(CYTOSINE-SPECIFIC METHYLTRANSFERASE BANI) (M.BANI)"
assert record.query_letters==428
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==0
assert len(record.alignments)==0
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Nov 26, 1999  5:52 PM',)]
assert record.ka_params==[0.320, 0.140, 0.403]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==13588598
assert record.num_sequences==87792
assert record.num_extends==162273
assert record.num_good_extends==546
assert record.num_seqs_better_e==0
assert record.hsps_no_gap==0
assert record.hsps_prelim_gapped==0
assert record.hsps_gapped==0
assert record.query_length==428
assert record.database_length==10666284
assert record.effective_hsp_length==48
assert record.effective_query_length==380
assert record.effective_database_length==6452268
assert record.effective_search_space==2451861840
assert record.effective_search_space_used==2451861840
assert record.frameshift==('50,','0.1')
assert record.threshold==13
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.4)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.8)
assert record.blast_cutoff==(61,28.2)

handle = open('Blast/bt018')
record = parser.parse(handle)
assert record.application=="TBLASTX"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|1348853|gb|G26621|G26621 human STS\nSTS_D12006.\x01gi|1396339|gb|G27620|G27620 human STS SHGC-32705."
assert record.query_letters==615
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==19
assert record.descriptions[0].title=="gi|1348853|gb|G26621|G26621 human STS STS_D12006. >gi|1396339|g..."
assert record.descriptions[0].score==398
assert record.descriptions[0].e==1.0000000000000001e-111
assert record.descriptions[1].title=="gi|1348016|gb|G25784|G25784 human STS EST47998."
assert record.descriptions[1].score==302
assert record.descriptions[1].e==9.9999999999999996e-83
assert record.descriptions[2].title=="gi|3403105|gb|G41148|G41148 Z7324 Zebrafish AB Danio rerio STS ..."
assert record.descriptions[2].score==31
assert record.descriptions[2].e==0.600000
assert record.descriptions[3].title=="gi|1234967|emb|Z53521|HSB316XA9 H.sapiens (D22S1166) DNA segmen..."
assert record.descriptions[3].score==30
assert record.descriptions[3].e==0.830000
assert record.descriptions[4].title=="gi|4185670|gb|G42865|G42865 Xq3845 KWOK Homo sapiens STS genomi..."
assert record.descriptions[4].score==30
assert record.descriptions[4].e==1.100000
assert record.descriptions[5].title=="gi|4757419|gb|G49246.1|G49246 stbK116F5_30376 chromosome 22 gen..."
assert record.descriptions[5].score==29
assert record.descriptions[5].e==1.600000
assert record.descriptions[6].title=="gi|6120694|gb|G55375.1|G55375 SHGC-100697 Human Homo sapiens ST..."
assert record.descriptions[6].score==29
assert record.descriptions[6].e==3.000000
assert record.descriptions[7].title=="gi|5225124|gb|G53947.1|G53947 SHGC-85304 Human Homo sapiens STS..."
assert record.descriptions[7].score==29
assert record.descriptions[7].e==3.000000
assert record.descriptions[8].title=="gi|1311530|gb|L77996|HUMSWX945 Human chromosome X STS sWXD945, ..."
assert record.descriptions[8].score==28
assert record.descriptions[8].e==4.100000
assert record.descriptions[9].title=="gi|2734624|gb|G36957|G36957 SHGC-56642 Human Homo sapiens STS cDNA"
assert record.descriptions[9].score==28
assert record.descriptions[9].e==4.100000
assert record.descriptions[10].title=="gi|4631600|dbj|AU046965.1|AU046965 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[10].score==28
assert record.descriptions[10].e==4.100000
assert record.descriptions[11].title=="gi|4631518|dbj|AU046883.1|AU046883 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[11].score==28
assert record.descriptions[11].e==4.100000
assert record.descriptions[12].title=="gi|2226478|gb|G33174|G33174 human STS SHGC-6097 clone pG-2470"
assert record.descriptions[12].score==28
assert record.descriptions[12].e==4.100000
assert record.descriptions[13].title=="gi|938611|gb|G08061|G08061 human STS CHLC.GGAA7E02.P7438 clone ..."
assert record.descriptions[13].score==27
assert record.descriptions[13].e==7.700000
assert record.descriptions[14].title=="gi|307789|gb|L18105|HUMUT1736 Human STS UT1736."
assert record.descriptions[14].score==27
assert record.descriptions[14].e==7.700000
assert record.descriptions[15].title=="gi|859804|gb|G06559|G06559 human STS WI-7401."
assert record.descriptions[15].score==27
assert record.descriptions[15].e==7.700000
assert record.descriptions[16].title=="gi|4493143|gb|G46852.1|G46852 Z14841_1 Zebrafish AB Danio rerio..."
assert record.descriptions[16].score==27
assert record.descriptions[16].e==7.700000
assert record.descriptions[17].title=="gi|4492122|gb|G45831.1|G45831 Z4588_1 Zebrafish AB Danio rerio ..."
assert record.descriptions[17].score==27
assert record.descriptions[17].e==7.700000
assert record.descriptions[18].title=="gi|6121804|gb|G56635.1|G56635 SHGC-102032 Human Homo sapiens ST..."
assert record.descriptions[18].score==27
assert record.descriptions[18].e==7.700000
assert len(record.alignments)==19
assert record.alignments[0].title==">gi|1348853|gb|G26621|G26621 human STS STS_D12006. >gi|1396339|gb|G27620|G27620 human STS SHGC-32705."
assert record.alignments[0].length==615
assert record.alignments[1].title==">gi|1348016|gb|G25784|G25784 human STS EST47998."
assert record.alignments[1].length==617
assert record.alignments[2].title==">gi|3403105|gb|G41148|G41148 Z7324 Zebrafish AB Danio rerio STS genomic"
assert record.alignments[2].length==351
assert record.alignments[3].title==">gi|1234967|emb|Z53521|HSB316XA9 H.sapiens (D22S1166) DNA segment containing (CA) repeat; clone AFMb316xa9; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[3].length==345
assert record.alignments[4].title==">gi|4185670|gb|G42865|G42865 Xq3845 KWOK Homo sapiens STS genomic, sequence tagged site [Homo sapiens]"
assert record.alignments[4].length==1200
assert record.alignments[5].title==">gi|4757419|gb|G49246.1|G49246 stbK116F5_30376 chromosome 22 genomic clone Homo sapiens STS genomic clone 116F5, sequence tagged site"
assert record.alignments[5].length==375
assert record.alignments[6].title==">gi|6120694|gb|G55375.1|G55375 SHGC-100697 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[6].length==460
assert record.alignments[7].title==">gi|5225124|gb|G53947.1|G53947 SHGC-85304 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[7].length==444
assert record.alignments[8].title==">gi|1311530|gb|L77996|HUMSWX945 Human chromosome X STS sWXD945, single read."
assert record.alignments[8].length==196
assert record.alignments[9].title==">gi|2734624|gb|G36957|G36957 SHGC-56642 Human Homo sapiens STS cDNA"
assert record.alignments[9].length==466
assert record.alignments[10].title==">gi|4631600|dbj|AU046965.1|AU046965 Rattus norvegicus, OTSUKA clone, 108a02, microsatellite sequence, sequence tagged site"
assert record.alignments[10].length==330
assert record.alignments[11].title==">gi|4631518|dbj|AU046883.1|AU046883 Rattus norvegicus, OTSUKA clone, 085f03, microsatellite sequence, sequence tagged site"
assert record.alignments[11].length==351
assert record.alignments[12].title==">gi|2226478|gb|G33174|G33174 human STS SHGC-6097 clone pG-2470"
assert record.alignments[12].length==299
assert record.alignments[13].title==">gi|938611|gb|G08061|G08061 human STS CHLC.GGAA7E02.P7438 clone GGAA7E02"
assert record.alignments[13].length==338
assert record.alignments[14].title==">gi|307789|gb|L18105|HUMUT1736 Human STS UT1736."
assert record.alignments[14].length==355
assert record.alignments[15].title==">gi|859804|gb|G06559|G06559 human STS WI-7401."
assert record.alignments[15].length==3280
assert record.alignments[16].title==">gi|4493143|gb|G46852.1|G46852 Z14841_1 Zebrafish AB Danio rerio STS genomic clone Z14841 5', sequence tagged site"
assert record.alignments[16].length==291
assert record.alignments[17].title==">gi|4492122|gb|G45831.1|G45831 Z4588_1 Zebrafish AB Danio rerio STS genomic clone Z4588 5', sequence tagged site"
assert record.alignments[17].length==398
assert record.alignments[18].title==">gi|6121804|gb|G56635.1|G56635 SHGC-102032 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[18].length==541
assert record.alignments[0].hsps[0].score==796
assert record.alignments[0].hsps[0].bits==367
assert record.alignments[0].hsps[0].expect==1e-102
assert record.alignments[0].hsps[1].score==759
assert record.alignments[0].hsps[1].bits==350
assert record.alignments[0].hsps[1].expect==3e-97
assert record.alignments[0].hsps[2].score==387
assert record.alignments[0].hsps[2].bits==180
assert record.alignments[0].hsps[2].expect==9e-91
assert record.alignments[0].hsps[3].score==368
assert record.alignments[0].hsps[3].bits==171
assert record.alignments[0].hsps[3].expect==9e-91
assert record.alignments[0].hsps[4].score==864
assert record.alignments[0].hsps[4].bits==398
assert record.alignments[0].hsps[4].expect==1e-111
assert record.alignments[0].hsps[5].score==846
assert record.alignments[0].hsps[5].bits==390
assert record.alignments[0].hsps[5].expect==1e-109
assert record.alignments[0].hsps[6].score==684
assert record.alignments[0].hsps[6].bits==316
assert record.alignments[0].hsps[6].expect==7e-87
assert len(record.alignments[0].hsps)==7
assert record.alignments[1].hsps[0].score==366
assert record.alignments[1].hsps[0].bits==170
assert record.alignments[1].hsps[0].expect==3e-63
assert record.alignments[1].hsps[1].score==188
assert record.alignments[1].hsps[1].bits==89.0
assert record.alignments[1].hsps[1].expect==3e-63
assert record.alignments[1].hsps[2].score==590
assert record.alignments[1].hsps[2].bits==273
assert record.alignments[1].hsps[2].expect==7e-74
assert record.alignments[1].hsps[3].score==593
assert record.alignments[1].hsps[3].bits==274
assert record.alignments[1].hsps[3].expect==8e-76
assert record.alignments[1].hsps[4].score==53
assert record.alignments[1].hsps[4].bits==27.2
assert record.alignments[1].hsps[4].expect==8e-76
assert record.alignments[1].hsps[5].score==653
assert record.alignments[1].hsps[5].bits==302
assert record.alignments[1].hsps[5].expect==1e-82
assert record.alignments[1].hsps[6].score==598
assert record.alignments[1].hsps[6].bits==276
assert record.alignments[1].hsps[6].expect==5e-75
assert record.alignments[1].hsps[7].score==628
assert record.alignments[1].hsps[7].bits==290
assert record.alignments[1].hsps[7].expect==4e-79
assert len(record.alignments[1].hsps)==8
assert record.alignments[2].hsps[0].score==61
assert record.alignments[2].hsps[0].bits==30.8
assert record.alignments[2].hsps[0].expect==0.60
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==60
assert record.alignments[3].hsps[0].bits==30.4
assert record.alignments[3].hsps[0].expect==0.83
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==59
assert record.alignments[4].hsps[0].bits==29.9
assert record.alignments[4].hsps[0].expect==1.1
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==58
assert record.alignments[5].hsps[0].bits==29.5
assert record.alignments[5].hsps[0].expect==1.6
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==56
assert record.alignments[6].hsps[0].bits==28.6
assert record.alignments[6].hsps[0].expect==3.0
assert len(record.alignments[6].hsps)==1
assert record.alignments[7].hsps[0].score==56
assert record.alignments[7].hsps[0].bits==28.6
assert record.alignments[7].hsps[0].expect==3.0
assert len(record.alignments[7].hsps)==1
assert record.alignments[8].hsps[0].score==55
assert record.alignments[8].hsps[0].bits==28.1
assert record.alignments[8].hsps[0].expect==4.1
assert len(record.alignments[8].hsps)==1
assert record.alignments[9].hsps[0].score==55
assert record.alignments[9].hsps[0].bits==28.1
assert record.alignments[9].hsps[0].expect==4.1
assert len(record.alignments[9].hsps)==1
assert record.alignments[10].hsps[0].score==55
assert record.alignments[10].hsps[0].bits==28.1
assert record.alignments[10].hsps[0].expect==4.1
assert len(record.alignments[10].hsps)==1
assert record.alignments[11].hsps[0].score==55
assert record.alignments[11].hsps[0].bits==28.1
assert record.alignments[11].hsps[0].expect==4.1
assert len(record.alignments[11].hsps)==1
assert record.alignments[12].hsps[0].score==55
assert record.alignments[12].hsps[0].bits==28.1
assert record.alignments[12].hsps[0].expect==4.1
assert len(record.alignments[12].hsps)==1
assert record.alignments[13].hsps[0].score==53
assert record.alignments[13].hsps[0].bits==27.2
assert record.alignments[13].hsps[0].expect==7.7
assert len(record.alignments[13].hsps)==1
assert record.alignments[14].hsps[0].score==53
assert record.alignments[14].hsps[0].bits==27.2
assert record.alignments[14].hsps[0].expect==7.7
assert len(record.alignments[14].hsps)==1
assert record.alignments[15].hsps[0].score==53
assert record.alignments[15].hsps[0].bits==27.2
assert record.alignments[15].hsps[0].expect==7.7
assert len(record.alignments[15].hsps)==1
assert record.alignments[16].hsps[0].score==53
assert record.alignments[16].hsps[0].bits==27.2
assert record.alignments[16].hsps[0].expect==7.7
assert len(record.alignments[16].hsps)==1
assert record.alignments[17].hsps[0].score==53
assert record.alignments[17].hsps[0].bits==27.2
assert record.alignments[17].hsps[0].expect==7.7
assert len(record.alignments[17].hsps)==1
assert record.alignments[18].hsps[0].score==53
assert record.alignments[18].hsps[0].bits==27.2
assert record.alignments[18].hsps[0].expect==7.7
assert len(record.alignments[18].hsps)==1
assert record.alignments[0].hsps[0].identities==(192, 200)
assert record.alignments[0].hsps[0].positives==(192, 200)
assert record.alignments[0].hsps[1].identities==(195, 205)
assert record.alignments[0].hsps[1].positives==(195, 205)
assert record.alignments[0].hsps[2].identities==(74, 74)
assert record.alignments[0].hsps[2].positives==(74, 74)
assert record.alignments[0].hsps[3].identities==(114, 114)
assert record.alignments[0].hsps[3].positives==(114, 114)
assert record.alignments[0].hsps[4].identities==(205, 205)
assert record.alignments[0].hsps[4].positives==(205, 205)
assert record.alignments[0].hsps[5].identities==(196, 196)
assert record.alignments[0].hsps[5].positives==(196, 196)
assert record.alignments[0].hsps[6].identities==(146, 146)
assert record.alignments[0].hsps[6].positives==(146, 146)
assert record.alignments[1].hsps[0].identities==(71, 74)
assert record.alignments[1].hsps[0].positives==(71, 74)
assert record.alignments[1].hsps[1].identities==(42, 67)
assert record.alignments[1].hsps[1].positives==(43, 67)
assert record.alignments[1].hsps[2].identities==(121, 133)
assert record.alignments[1].hsps[2].positives==(121, 133)
assert record.alignments[1].hsps[3].identities==(112, 131)
assert record.alignments[1].hsps[3].positives==(112, 131)
assert record.alignments[1].hsps[4].identities==(9, 13)
assert record.alignments[1].hsps[4].positives==(10, 13)
assert record.alignments[1].hsps[5].identities==(128, 157)
assert record.alignments[1].hsps[5].positives==(132, 157)
assert record.alignments[1].hsps[6].identities==(122, 130)
assert record.alignments[1].hsps[6].positives==(122, 130)
assert record.alignments[1].hsps[7].identities==(119, 131)
assert record.alignments[1].hsps[7].positives==(120, 131)
assert record.alignments[2].hsps[0].identities==(11, 27)
assert record.alignments[2].hsps[0].positives==(18, 27)
assert record.alignments[3].hsps[0].identities==(10, 19)
assert record.alignments[3].hsps[0].positives==(13, 19)
assert record.alignments[4].hsps[0].identities==(10, 24)
assert record.alignments[4].hsps[0].positives==(14, 24)
assert record.alignments[5].hsps[0].identities==(15, 34)
assert record.alignments[5].hsps[0].positives==(17, 34)
assert record.alignments[6].hsps[0].identities==(9, 28)
assert record.alignments[6].hsps[0].positives==(16, 28)
assert record.alignments[7].hsps[0].identities==(10, 24)
assert record.alignments[7].hsps[0].positives==(15, 24)
assert record.alignments[8].hsps[0].identities==(9, 33)
assert record.alignments[8].hsps[0].positives==(19, 33)
assert record.alignments[9].hsps[0].identities==(10, 15)
assert record.alignments[9].hsps[0].positives==(11, 15)
assert record.alignments[10].hsps[0].identities==(12, 39)
assert record.alignments[10].hsps[0].positives==(19, 39)
assert record.alignments[11].hsps[0].identities==(8, 26)
assert record.alignments[11].hsps[0].positives==(15, 26)
assert record.alignments[12].hsps[0].identities==(13, 38)
assert record.alignments[12].hsps[0].positives==(19, 38)
assert record.alignments[13].hsps[0].identities==(12, 24)
assert record.alignments[13].hsps[0].positives==(13, 24)
assert record.alignments[14].hsps[0].identities==(12, 24)
assert record.alignments[14].hsps[0].positives==(13, 24)
assert record.alignments[15].hsps[0].identities==(9, 21)
assert record.alignments[15].hsps[0].positives==(13, 21)
assert record.alignments[16].hsps[0].identities==(8, 20)
assert record.alignments[16].hsps[0].positives==(13, 20)
assert record.alignments[17].hsps[0].identities==(7, 22)
assert record.alignments[17].hsps[0].positives==(13, 22)
assert record.alignments[18].hsps[0].identities==(9, 20)
assert record.alignments[18].hsps[0].positives==(13, 20)
assert record.alignments[0].hsps[0].frame==("+2", "+2")
assert record.alignments[0].hsps[1].frame==("+1", "+1")
assert record.alignments[0].hsps[2].frame==("+3", "+3")
assert record.alignments[0].hsps[3].frame==("+3", "+3")
assert record.alignments[0].hsps[4].frame==("-1", "-1")
assert record.alignments[0].hsps[5].frame==("-3", "-3")
assert record.alignments[0].hsps[6].frame==("-2", "-2")
assert record.alignments[1].hsps[0].frame==("+3", "+3")
assert record.alignments[1].hsps[1].frame==("+3", "+3")
assert record.alignments[1].hsps[2].frame==("+1", "+1")
assert record.alignments[1].hsps[3].frame==("+2", "+2")
assert record.alignments[1].hsps[4].frame==("+2", "+2")
assert record.alignments[1].hsps[5].frame==("-3", "-2")
assert record.alignments[1].hsps[6].frame==("-2", "-1")
assert record.alignments[1].hsps[7].frame==("-1", "-3")
assert record.alignments[2].hsps[0].frame==("+1", "+2")
assert record.alignments[3].hsps[0].frame==("-3", "+1")
assert record.alignments[4].hsps[0].frame==("-3", "-2")
assert record.alignments[5].hsps[0].frame==("-3", "-2")
assert record.alignments[6].hsps[0].frame==("+3", "+1")
assert record.alignments[7].hsps[0].frame==("+1", "-2")
assert record.alignments[8].hsps[0].frame==("+3", "-2")
assert record.alignments[9].hsps[0].frame==("-3", "-3")
assert record.alignments[10].hsps[0].frame==("+2", "-3")
assert record.alignments[11].hsps[0].frame==("+1", "+3")
assert record.alignments[12].hsps[0].frame==("+1", "+2")
assert record.alignments[13].hsps[0].frame==("-3", "+3")
assert record.alignments[14].hsps[0].frame==("-3", "-1")
assert record.alignments[15].hsps[0].frame==("-1", "-1")
assert record.alignments[16].hsps[0].frame==("-3", "+1")
assert record.alignments[17].hsps[0].frame==("+1", "+1")
assert record.alignments[18].hsps[0].frame==("-1", "+2")
assert record.alignments[0].hsps[0].query=="IRMPLHS*DSSFCPL*QEKWECMXXXXXXXXRPKRCLQPHPLNWP*LGLNALMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPLFWWTNNLGTPNPGRKNIQHYEXALX*S*NGFPKLVTHGPGXVXLXLXGLSSFQEFXSTVANPWGX*XXXXXXXXXFXXGXRXXXLXXXXGGCXXVXXVXXXWXXXF"
assert record.alignments[0].hsps[0].match=="IRMPLHS*DSSFCPL*QEKWECM        RPKRCLQPHPLNWP*LGLNALMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPLFWWTNNLGTPNPGRKNIQHYE AL *S*NGFPKLVTHGPG V L L GLSSFQEF STVANPWG *         F  G R   L    GGC  V  V   W   F"
assert record.alignments[0].hsps[0].sbjct=="IRMPLHS*DSSFCPL*QEKWECMQSSQKKQKRPKRCLQPHPLNWP*LGLNALMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPLFWWTNNLGTPNPGRKNIQHYEXALX*S*NGFPKLVTHGPGXVXLXLXGLSSFQEFXSTVANPWGX*XXXXXXXXXFXXGXRXXXLXXXXGGCXXVXXVXXXWXXXF"
assert record.alignments[0].hsps[0].query_start==2
assert record.alignments[0].hsps[0].query_end==601
assert record.alignments[0].hsps[0].sbjct_start==2
assert record.alignments[0].hsps[0].sbjct_end==601
assert record.alignments[0].hsps[1].query=="DQNAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKCTNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL*XSTXVILKWISKAGYTWPWXXXIXFXXAFFXXXXXXXXXXXXXXLXGXVX*XXXXSXGPXXXXVXXPXGWVXXGXFXFXXXXXXFXXXLG"
assert record.alignments[0].hsps[1].match=="DQNAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKCTNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL* ST VILKWISKAGYTWPW   I F  AFF              L G V *    S GP    V  P GWV  G F F      F   LG"
assert record.alignments[0].hsps[1].sbjct=="DQNAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKCTNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL*XSTXVILKWISKAGYTWPWXXXIXFXXAFFXSGVXVNGGKSXGXLXGXVX*XXXXSXGPXXXXVXXPXGWVXXGXFXFXXXXXXFXXXLG"
assert record.alignments[0].hsps[1].query_start==1
assert record.alignments[0].hsps[1].query_end==615
assert record.alignments[0].hsps[1].sbjct_start==1
assert record.alignments[0].hsps[1].sbjct_end==615
assert record.alignments[0].hsps[2].query=="SECPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*MH*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[0].hsps[2].match=="SECPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*MH*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[0].hsps[2].sbjct=="SECPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*MH*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[0].hsps[2].query_start==3
assert record.alignments[0].hsps[2].query_end==224
assert record.alignments[0].hsps[2].sbjct_start==3
assert record.alignments[0].hsps[2].sbjct_end==224
assert record.alignments[0].hsps[3].query=="YSHCSGGQIT*VLQILAGKIYSIMXQHSXNPKMDFQSWLHMALXXSY*XXXGFLXFRSXGQRWQIXGVXXWXGXLXXXXFXGAXGGXGXXSXXVGAXGXXXFXXXGXXXXXXXW"
assert record.alignments[0].hsps[3].match=="YSHCSGGQIT*VLQILAGKIYSIM QHS NPKMDFQSWLHMAL  SY*   GFL FRS GQRWQI GV  W G L    F GA GG G  S  VGA G   F   G       W"
assert record.alignments[0].hsps[3].sbjct=="YSHCSGGQIT*VLQILAGKIYSIMXQHSXNPKMDFQSWLHMALXXSY*XXXGFLXFRSXGQRWQIXGVXXWXGXLXXXXFXGAXGGXGXXSXXVGAXGXXXFXXXGXXXXXXXW"
assert record.alignments[0].hsps[3].query_start==273
assert record.alignments[0].hsps[3].query_end==614
assert record.alignments[0].hsps[3].sbjct_start==273
assert record.alignments[0].hsps[3].sbjct_end==614
assert record.alignments[0].hsps[4].query=="PKXXXKXXXXXXKTKXTXXHPPXWXXNXXXLWPXXXXXXLXNXXX*XPXGFATVDXNS*XEESPXKXNXTXPGPCVTSFGNPF*DYXSAXS*CCIFFLPGFGVPRLFVHQNNGCKGEN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCISAFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGILI"
assert record.alignments[0].hsps[4].match=="PK   K      KTK T  HPP W  N   LWP      L N   * P GFATVD NS* EESP K N T PGPCVTSFGNPF*DY SA S*CCIFFLPGFGVPRLFVHQNNGCKGEN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCISAFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGILI"
assert record.alignments[0].hsps[4].sbjct=="PKXXXKXXXXXXKTKXTXXHPPXWXXNXXXLWPXXXXXXLXNXXX*XPXGFATVDXNS*XEESPXKXNXTXPGPCVTSFGNPF*DYXSAXS*CCIFFLPGFGVPRLFVHQNNGCKGEN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCISAFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGILI"
assert record.alignments[0].hsps[4].query_start==615
assert record.alignments[0].hsps[4].query_end==1
assert record.alignments[0].hsps[4].sbjct_start==615
assert record.alignments[0].hsps[4].sbjct_end==1
assert record.alignments[0].hsps[5].query=="PXXXNXXNPXAPTXXXX*PXPPXAPXKXXXXXXPXXLXTPXICHR*PKLLKXRKPX*XQXDXXRAMCNQLWKSILGLXECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD"
assert record.alignments[0].hsps[5].match=="P   N  NP APT    *P PP AP K      P  L TP ICHR*PKLLK RKP * Q D  RAMCNQLWKSILGL EC FIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD"
assert record.alignments[0].hsps[5].sbjct=="PXXXNXXNPXAPTXXXX*PXPPXAPXKXXXXXXPXXLXTPXICHR*PKLLKXRKPX*XQXDXXRAMCNQLWKSILGLXECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD"
assert record.alignments[0].hsps[5].query_start==589
assert record.alignments[0].hsps[5].query_end==2
assert record.alignments[0].hsps[5].sbjct_start==589
assert record.alignments[0].hsps[5].sbjct_end==2
assert record.alignments[0].hsps[6].query=="EXKKAXXXSIXXXQGHV*PALEIHFRIT*VLXHNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALVHLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGAF*"
assert record.alignments[0].hsps[6].match=="E KKA   SI   QGHV*PALEIHFRIT*VL HNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALVHLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGAF*"
assert record.alignments[0].hsps[6].sbjct=="EXKKAXXXSIXXXQGHV*PALEIHFRIT*VLXHNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALVHLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGAF*"
assert record.alignments[0].hsps[6].query_start==440
assert record.alignments[0].hsps[6].query_end==3
assert record.alignments[0].hsps[6].sbjct_start==440
assert record.alignments[0].hsps[6].sbjct_end==3
assert record.alignments[1].hsps[0].query=="SECPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*MH*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[1].hsps[0].match=="S CPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*  *CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[1].hsps[0].sbjct=="SXCPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*XX*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[1].hsps[0].query_start==3
assert record.alignments[1].hsps[0].query_end==224
assert record.alignments[1].hsps[0].sbjct_start==3
assert record.alignments[1].hsps[0].sbjct_end==224
assert record.alignments[1].hsps[1].query=="YSHCSGGQIT*VLQILAGKIYSIMXQHSXNPKMDFQSWLHMALXXSY*XXXGFLXFRSXGQRWQIXG"
assert record.alignments[1].hsps[1].match=="YSHCSGGQIT*VLQILAGKIYSIM QHS   K  FQSWLHM     +     F   R  GQR Q  G"
assert record.alignments[1].hsps[1].sbjct=="YSHCSGGQIT*VLQILAGKIYSIMKQHSVILKWIFQSWLHMXCKVLFKFKRPFSFTRGLGQRXQTPG"
assert record.alignments[1].hsps[1].query_start==273
assert record.alignments[1].hsps[1].query_end==473
assert record.alignments[1].hsps[1].sbjct_start==273
assert record.alignments[1].hsps[1].sbjct_end==473
assert record.alignments[1].hsps[2].query=="DQNAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKCTNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL*XSTXVILKWISKAGYTWP"
assert record.alignments[1].hsps[2].match=="DQ APPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTK  NAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL* ST       SKAGYT P"
assert record.alignments[1].hsps[2].sbjct=="DQXAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKXXNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL*NSTQ*S*NGFSKAGYTXP"
assert record.alignments[1].hsps[2].query_start==1
assert record.alignments[1].hsps[2].query_end==399
assert record.alignments[1].hsps[2].sbjct_start==1
assert record.alignments[1].hsps[2].sbjct_end==399
assert record.alignments[1].hsps[3].query=="IRMPLHS*DSSFCPL*QEKWECMXXXXXXXXRPKRCLQPHPLNWP*LGLNALMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPLFWWTNNLGTPNPGRKNIQHYEXALX*S*NGFPKLVTH"
assert record.alignments[1].hsps[3].match=="IR PLHS*DSSFCPL*QEKWECM        RPKRCLQPHPLNWP*LGL  LMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQF PLQPLFWWTNNLGTPNPGRKNIQHYE AL      FPKLVTH"
assert record.alignments[1].hsps[3].sbjct=="IRXPLHS*DSSFCPL*QEKWECMQSSQKKQKRPKRCLQPHPLNWP*LGLXPLMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFFPLQPLFWWTNNLGTPNPGRKNIQHYETALSNPKMDFPKLVTH"
assert record.alignments[1].hsps[3].query_start==2
assert record.alignments[1].hsps[3].query_end==394
assert record.alignments[1].hsps[3].sbjct_start==2
assert record.alignments[1].hsps[3].sbjct_end==394
assert record.alignments[1].hsps[4].query=="FQEFXSTVANPWG"
assert record.alignments[1].hsps[4].match=="+Q F ST ANPWG"
assert record.alignments[1].hsps[4].sbjct=="YQGFRSTXANPWG"
assert record.alignments[1].hsps[4].query_start==437
assert record.alignments[1].hsps[4].query_end==475
assert record.alignments[1].hsps[4].sbjct_start==437
assert record.alignments[1].hsps[4].sbjct_end==475
assert record.alignments[1].hsps[5].query=="PXICHR*PKLLKXRKPX*XQXDXXRAMCNQLWKSILGLXECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD"
assert record.alignments[1].hsps[5].match=="P +C R*PK L   K         + MCNQLWK    + EC FIMLYIFPARIWST*VICPPEQWL*R+KLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*  *SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEG SD"
assert record.alignments[1].hsps[5].sbjct=="PGVCXR*PKPLVNEKGLLNLNKTLQXMCNQLWKIHFRITECCFIMLYIFPARIWST*VICPPEQWL*RKKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*XX*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGXSD"
assert record.alignments[1].hsps[5].query_start==472
assert record.alignments[1].hsps[5].query_end==2
assert record.alignments[1].hsps[5].sbjct_start==472
assert record.alignments[1].hsps[5].sbjct_end==2
assert record.alignments[1].hsps[6].query=="GHV*PALEIHFRIT*VLXHNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALVHLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGA"
assert record.alignments[1].hsps[6].match=="GHV*PALE  F   *VL HNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFAL  LVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGA"
assert record.alignments[1].hsps[6].sbjct=="GHV*PALENPF*DY*VLFHNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALXXLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGA"
assert record.alignments[1].hsps[6].query_start==398
assert record.alignments[1].hsps[6].query_end==9
assert record.alignments[1].hsps[6].sbjct_start==398
assert record.alignments[1].hsps[6].sbjct_end==9
assert record.alignments[1].hsps[7].query=="CVTSFGNPF*DYXSAXS*CCIFFLPGFGVPRLFVHQNNGCKGEN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCISAFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGILI"
assert record.alignments[1].hsps[7].match=="CVTSFG       SA S*CCIFFLPGFGVPRLFVHQNNGCKG+N*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCI  FSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRG LI"
assert record.alignments[1].hsps[7].sbjct=="CVTSFGKSILGLLSAVS*CCIFFLPGFGVPRLFVHQNNGCKGKN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCIXGFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGXLI"
assert record.alignments[1].hsps[7].query_start==393
assert record.alignments[1].hsps[7].query_end==1
assert record.alignments[1].hsps[7].sbjct_start==393
assert record.alignments[1].hsps[7].sbjct_end==1
assert record.alignments[2].hsps[0].query=="LCNLFYRIFSAASDHCSIFSFTAIVLV"
assert record.alignments[2].hsps[0].match=="L NL +R+F++  DHCS+  F   V++"
assert record.alignments[2].hsps[0].sbjct=="LTNLIFRLFTSK*DHCSLKQFNNKVVI"
assert record.alignments[2].hsps[0].query_start==211
assert record.alignments[2].hsps[0].query_end==291
assert record.alignments[2].hsps[0].sbjct_start==2
assert record.alignments[2].hsps[0].sbjct_end==82
assert record.alignments[3].hsps[0].query=="IWST*VICPPEQWL*RRKL"
assert record.alignments[3].hsps[0].match=="+W T  +CPP  W *RR+L"
assert record.alignments[3].hsps[0].sbjct=="LWGTAPLCPPVPWA*RRQL"
assert record.alignments[3].hsps[0].query_start==316
assert record.alignments[3].hsps[0].query_end==260
assert record.alignments[3].hsps[0].sbjct_start==226
assert record.alignments[3].hsps[0].sbjct_end==282
assert record.alignments[4].hsps[0].query=="PFWPFLFLLRALHTFPLFLSKWTK"
assert record.alignments[4].hsps[0].match=="P WP + LL  +H+FP  +  W K"
assert record.alignments[4].hsps[0].sbjct=="PEWPSINLLPGMHSFPRVIPCWEK"
assert record.alignments[4].hsps[0].query_start==106
assert record.alignments[4].hsps[0].query_end==35
assert record.alignments[4].hsps[0].sbjct_start==716
assert record.alignments[4].hsps[0].sbjct_end==645
assert record.alignments[5].hsps[0].query=="PFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHS"
assert record.alignments[5].hsps[0].match=="P W   FLLR LHT P    K+    V+  EG S"
assert record.alignments[5].hsps[0].sbjct=="PGWTLRFLLRGLHTPPKAAEKFPSQGVAEPEGFS"
assert record.alignments[5].hsps[0].query_start==106
assert record.alignments[5].hsps[0].query_end==5
assert record.alignments[5].hsps[0].sbjct_start==176
assert record.alignments[5].hsps[0].sbjct_end==75
assert record.alignments[6].hsps[0].query=="STHETLHFVHFDRKSGNVCRALKRNKKG"
assert record.alignments[6].hsps[0].match=="+TH    F+H  +K  N+C  +  N++G"
assert record.alignments[6].hsps[0].sbjct=="TTHSMTVFLHVKKKLSNICIYMPENREG"
assert record.alignments[6].hsps[0].query_start==15
assert record.alignments[6].hsps[0].query_end==98
assert record.alignments[6].hsps[0].sbjct_start==82
assert record.alignments[6].hsps[0].sbjct_end==165
assert record.alignments[7].hsps[0].query=="HCFLCNLFYRIFSAASDHCSIFSF"
assert record.alignments[7].hsps[0].match=="HC L ++F RIF    +  S+F+F"
assert record.alignments[7].hsps[0].sbjct=="HCHLSSVFCRIFFTLEESLSLFAF"
assert record.alignments[7].hsps[0].query_start==202
assert record.alignments[7].hsps[0].query_end==273
assert record.alignments[7].hsps[0].sbjct_start==443
assert record.alignments[7].hsps[0].sbjct_end==372
assert record.alignments[8].hsps[0].query=="ETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHI"
assert record.alignments[8].hsps[0].match=="E +H+V F +  G +  + ++ ++ QN   SH+"
assert record.alignments[8].hsps[0].sbjct=="EKIHYVLF*KNCGKIMTS*RQRQENQNNLLSHV"
assert record.alignments[8].hsps[0].query_start==24
assert record.alignments[8].hsps[0].query_end==122
assert record.alignments[8].hsps[0].sbjct_start==141
assert record.alignments[8].hsps[0].sbjct_end==43
assert record.alignments[9].hsps[0].query=="LKAPFWPFLFLLRAL"
assert record.alignments[9].hsps[0].match=="L  PFW FLF L+AL"
assert record.alignments[9].hsps[0].sbjct=="LSVPFWKFLFYLQAL"
assert record.alignments[9].hsps[0].query_start==115
assert record.alignments[9].hsps[0].query_end==71
assert record.alignments[9].hsps[0].sbjct_start==242
assert record.alignments[9].hsps[0].sbjct_end==198
assert record.alignments[10].hsps[0].query=="PLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPL"
assert record.alignments[10].hsps[0].match=="PL    RLF + +  ++ C T     L  +  F P++PL"
assert record.alignments[10].hsps[0].sbjct=="PLKTQGRLFYSKVSLFLKCRTTLCLFLNVSEGFXPIEPL"
assert record.alignments[10].hsps[0].query_start==167
assert record.alignments[10].hsps[0].query_end==283
assert record.alignments[10].hsps[0].sbjct_start==211
assert record.alignments[10].hsps[0].sbjct_end==95
assert record.alignments[11].hsps[0].query=="CFLCNLFYRIFSAASDHCSIFSFTAI"
assert record.alignments[11].hsps[0].match=="CF+C L  +I+  A +   +  FT++"
assert record.alignments[11].hsps[0].sbjct=="CFICKLVLKIYLRAEERTQLIEFTSL"
assert record.alignments[11].hsps[0].query_start==205
assert record.alignments[11].hsps[0].query_end==282
assert record.alignments[11].hsps[0].sbjct_start==249
assert record.alignments[11].hsps[0].sbjct_end==326
assert record.alignments[12].hsps[0].query=="KSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVL"
assert record.alignments[12].hsps[0].match=="+SA +I   HC    +   +FS ASD  S    T ++L"
assert record.alignments[12].hsps[0].sbjct=="QSAWIIGMSHCAWAVIHCLLFSTASDEKSAVDLTGVLL"
assert record.alignments[12].hsps[0].query_start==175
assert record.alignments[12].hsps[0].query_end==288
assert record.alignments[12].hsps[0].sbjct_start==119
assert record.alignments[12].hsps[0].sbjct_end==232
assert record.alignments[13].hsps[0].query=="WLKAPFWPFLFLLRALHTFPLFLS"
assert record.alignments[13].hsps[0].match=="WL  PF PFL  L +L   P  LS"
assert record.alignments[13].hsps[0].sbjct=="WLLFPFLPFLPFLPSLPFLPFLLS"
assert record.alignments[13].hsps[0].query_start==118
assert record.alignments[13].hsps[0].query_end==47
assert record.alignments[13].hsps[0].sbjct_start==153
assert record.alignments[13].hsps[0].sbjct_end==224
assert record.alignments[14].hsps[0].query=="WLKAPFWPFLFLLRALHTFPLFLS"
assert record.alignments[14].hsps[0].match=="WL  PF PFL  L +L   P  LS"
assert record.alignments[14].hsps[0].sbjct=="WLLFPFLPFLPFLPSLPFLPFLLS"
assert record.alignments[14].hsps[0].query_start==118
assert record.alignments[14].hsps[0].query_end==47
assert record.alignments[14].hsps[0].sbjct_start==235
assert record.alignments[14].hsps[0].sbjct_end==164
assert record.alignments[15].hsps[0].query=="CCIFFLPGFGVPRLFVHQNNG"
assert record.alignments[15].hsps[0].match=="C IF++P F    LF+H+  G"
assert record.alignments[15].hsps[0].sbjct=="CFIFYVPDFPWSNLFLHRGRG"
assert record.alignments[15].hsps[0].query_start==339
assert record.alignments[15].hsps[0].query_end==277
assert record.alignments[15].hsps[0].sbjct_start==1396
assert record.alignments[15].hsps[0].sbjct_end==1334
assert record.alignments[16].hsps[0].query=="NQLWKSILGLXECXFIMLYI"
assert record.alignments[16].hsps[0].match=="NQ+WK  L +  C F+ +Y+"
assert record.alignments[16].hsps[0].sbjct=="NQMWKHSLEVCMCVFVYIYV"
assert record.alignments[16].hsps[0].query_start==388
assert record.alignments[16].hsps[0].query_end==329
assert record.alignments[16].hsps[0].sbjct_start==37
assert record.alignments[16].hsps[0].sbjct_end==96
assert record.alignments[17].hsps[0].query=="CFLCNLFYRIFSAASDHCSIFS"
assert record.alignments[17].hsps[0].match=="CF   + +R+F+    HC+ F+"
assert record.alignments[17].hsps[0].sbjct=="CFALTVIWRVFAGCRPHCATFT"
assert record.alignments[17].hsps[0].query_start==205
assert record.alignments[17].hsps[0].query_end==270
assert record.alignments[17].hsps[0].sbjct_start==13
assert record.alignments[17].hsps[0].sbjct_end==78
assert record.alignments[18].hsps[0].query=="VKGFCISAFSPSQGQFRGCG"
assert record.alignments[18].hsps[0].match=="V GFC+++FSP +    G G"
assert record.alignments[18].hsps[0].sbjct=="VDGFCVTSFSPKKDTHPGSG"
assert record.alignments[18].hsps[0].query_start==174
assert record.alignments[18].hsps[0].query_end==115
assert record.alignments[18].hsps[0].sbjct_start==125
assert record.alignments[18].hsps[0].sbjct_end==184
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Nov 26, 1999  5:52 PM',)]
assert record.ka_params==[0.318, 0.135, 0.401]
assert record.matrix=='BLOSUM62'
assert record.num_hits==40473548
assert record.num_sequences==87792
assert record.num_extends==487631
assert record.num_good_extends==13175
assert record.num_seqs_better_e==38
assert record.query_length==205
assert record.database_length==10666284
assert record.effective_hsp_length==46
assert record.effective_query_length==158
assert record.effective_database_length==6627852
assert record.effective_search_space==1047200616
assert record.effective_search_space_used==1047200616
assert record.frameshift==('50,','0.1')
assert record.threshold==13
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(0,0.0)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(52,26.7)

handle = open('Blast/bt039')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|585505|sp|Q08386|MOPB_RHOCA MOLYBDENUM-PTERIN BINDING\nPROTEIN MOPB"
assert record.query_letters==270
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==13
assert record.descriptions[0].title=="gi|585505|sp|Q08386|MOPB_RHOCA MOLYBDENUM-PTERIN BINDING PROTEI..."
assert record.descriptions[0].score==467
assert record.descriptions[0].e==9.9999999999999999e-132
assert record.descriptions[1].title=="gi|585504|sp|Q08385|MOPA_RHOCA MOLYBDENUM-PTERIN BINDING PROTEI..."
assert record.descriptions[1].score==207
assert record.descriptions[1].e==2.0000000000000001e-53
assert record.descriptions[2].title=="gi|585492|sp|P37733|MODA_AZOVI MOLYBDENUM TRANSPORT PROTEIN MODA"
assert record.descriptions[2].score==145
assert record.descriptions[2].e==9.0000000000000002e-35
assert record.descriptions[3].title=="gi|1709070|sp|P46930|MODE_ECOLI MOLYBDENUM TRANSPORT PROTEIN MODE"
assert record.descriptions[3].score==87
assert record.descriptions[3].e==4.9999999999999999e-17
assert record.descriptions[4].title=="gi|1709071|sp|P45324|MODE_HAEIN MOLYBDENUM TRANSPORT PROTEIN MO..."
assert record.descriptions[4].score==54
assert record.descriptions[4].e==1.9999999999999999e-07
assert record.descriptions[5].title=="gi|585502|sp|P04952|MOP1_CLOPA MOLYBDENUM-PTERIN BINDING PROTEIN I"
assert record.descriptions[5].score==53
assert record.descriptions[5].e==5.9999999999999997e-07
assert record.descriptions[6].title=="gi|127241|sp|P08854|MOP2_CLOPA MOLYBDENUM-PTERIN BINDING PROTEI..."
assert record.descriptions[6].score==52
assert record.descriptions[6].e==0.000001
assert record.descriptions[7].title=="gi|585503|sp|P38366|MOP3_CLOPA MOLYBDENUM-PTERIN BINDING PROTEI..."
assert record.descriptions[7].score==51
assert record.descriptions[7].e==0.000003
assert record.descriptions[8].title=="gi|1170996|sp|P45183|MOP_HAEIN PROBABLE MOLYBDENUM-PTERIN BINDI..."
assert record.descriptions[8].score==46
assert record.descriptions[8].e==0.000050
assert record.descriptions[9].title=="gi|1709069|sp|P09833|MODC_ECOLI MOLYBDENUM TRANSPORT ATP-BINDIN..."
assert record.descriptions[9].score==38
assert record.descriptions[9].e==0.021000
assert record.descriptions[10].title=="gi|585500|sp|P37732|MODD_AZOVI MOLYBDENUM TRANSPORT ATP-BINDING..."
assert record.descriptions[10].score==33
assert record.descriptions[10].e==0.530000
assert record.descriptions[11].title=="gi|2507168|sp|P08838|PT1_BACSU PHOSPHOENOLPYRUVATE-PROTEIN PHOS..."
assert record.descriptions[11].score==30
assert record.descriptions[11].e==4.600000
assert record.descriptions[12].title=="gi|729786|sp|Q05355|HYDL_STRHA PUTATIVE POLYKETIDE HYDROXYLASE"
assert record.descriptions[12].score==29
assert record.descriptions[12].e==7.900000
assert len(record.alignments)==13
assert record.alignments[0].title==">gi|585505|sp|Q08386|MOPB_RHOCA MOLYBDENUM-PTERIN BINDING PROTEIN MOPB"
assert record.alignments[0].length==270
assert record.alignments[1].title==">gi|585504|sp|Q08385|MOPA_RHOCA MOLYBDENUM-PTERIN BINDING PROTEIN MOPA"
assert record.alignments[1].length==265
assert record.alignments[2].title==">gi|585492|sp|P37733|MODA_AZOVI MOLYBDENUM TRANSPORT PROTEIN MODA"
assert record.alignments[2].length==270
assert record.alignments[3].title==">gi|1709070|sp|P46930|MODE_ECOLI MOLYBDENUM TRANSPORT PROTEIN MODE"
assert record.alignments[3].length==262
assert record.alignments[4].title==">gi|1709071|sp|P45324|MODE_HAEIN MOLYBDENUM TRANSPORT PROTEIN MODE HOMOLOG"
assert record.alignments[4].length==255
assert record.alignments[5].title==">gi|585502|sp|P04952|MOP1_CLOPA MOLYBDENUM-PTERIN BINDING PROTEIN I"
assert record.alignments[5].length==68
assert record.alignments[6].title==">gi|127241|sp|P08854|MOP2_CLOPA MOLYBDENUM-PTERIN BINDING PROTEIN II"
assert record.alignments[6].length==68
assert record.alignments[7].title==">gi|585503|sp|P38366|MOP3_CLOPA MOLYBDENUM-PTERIN BINDING PROTEIN III"
assert record.alignments[7].length==68
assert record.alignments[8].title==">gi|1170996|sp|P45183|MOP_HAEIN PROBABLE MOLYBDENUM-PTERIN BINDING PROTEIN"
assert record.alignments[8].length==69
assert record.alignments[9].title==">gi|1709069|sp|P09833|MODC_ECOLI MOLYBDENUM TRANSPORT ATP-BINDING PROTEIN MODC"
assert record.alignments[9].length==352
assert record.alignments[10].title==">gi|585500|sp|P37732|MODD_AZOVI MOLYBDENUM TRANSPORT ATP-BINDING PROTEIN MODD"
assert record.alignments[10].length==380
assert record.alignments[11].title==">gi|2507168|sp|P08838|PT1_BACSU PHOSPHOENOLPYRUVATE-PROTEIN PHOSPHOTRANSFERASE (PHOSPHOTRANSFERASE SYSTEM, ENZYME I)"
assert record.alignments[11].length==570
assert record.alignments[12].title==">gi|729786|sp|Q05355|HYDL_STRHA PUTATIVE POLYKETIDE HYDROXYLASE"
assert record.alignments[12].length==555
assert record.alignments[0].hsps[0].score==1189
assert record.alignments[0].hsps[0].bits==467
assert record.alignments[0].hsps[0].expect==1e-131
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==521
assert record.alignments[1].hsps[0].bits==207
assert record.alignments[1].hsps[0].expect==2e-53
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==362
assert record.alignments[2].hsps[0].bits==145
assert record.alignments[2].hsps[0].expect==9e-35
assert record.alignments[2].hsps[1].score==98
assert record.alignments[2].hsps[1].bits==42.6
assert record.alignments[2].hsps[1].expect==8e-04
assert len(record.alignments[2].hsps)==2
assert record.alignments[3].hsps[0].score==211
assert record.alignments[3].hsps[0].bits==86.6
assert record.alignments[3].hsps[0].expect==5e-17
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==128
assert record.alignments[4].hsps[0].bits==54.3
assert record.alignments[4].hsps[0].expect==2e-07
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==125
assert record.alignments[5].hsps[0].bits==53.1
assert record.alignments[5].hsps[0].expect==6e-07
assert record.alignments[5].hsps[1].score==84
assert record.alignments[5].hsps[1].bits==37.1
assert record.alignments[5].hsps[1].expect==0.036
assert len(record.alignments[5].hsps)==2
assert record.alignments[6].hsps[0].score==123
assert record.alignments[6].hsps[0].bits==52.3
assert record.alignments[6].hsps[0].expect==1e-06
assert record.alignments[6].hsps[1].score==86
assert record.alignments[6].hsps[1].bits==37.9
assert record.alignments[6].hsps[1].expect==0.021
assert len(record.alignments[6].hsps)==2
assert record.alignments[7].hsps[0].score==119
assert record.alignments[7].hsps[0].bits==50.8
assert record.alignments[7].hsps[0].expect==3e-06
assert record.alignments[7].hsps[1].score==83
assert record.alignments[7].hsps[1].bits==36.7
assert record.alignments[7].hsps[1].expect==0.047
assert len(record.alignments[7].hsps)==2
assert record.alignments[8].hsps[0].score==108
assert record.alignments[8].hsps[0].bits==46.5
assert record.alignments[8].hsps[0].expect==5e-05
assert len(record.alignments[8].hsps)==1
assert record.alignments[9].hsps[0].score==86
assert record.alignments[9].hsps[0].bits==37.9
assert record.alignments[9].hsps[0].expect==0.021
assert len(record.alignments[9].hsps)==1
assert record.alignments[10].hsps[0].score==74
assert record.alignments[10].hsps[0].bits==33.2
assert record.alignments[10].hsps[0].expect==0.53
assert len(record.alignments[10].hsps)==1
assert record.alignments[11].hsps[0].score==66
assert record.alignments[11].hsps[0].bits==30.1
assert record.alignments[11].hsps[0].expect==4.6
assert len(record.alignments[11].hsps)==1
assert record.alignments[12].hsps[0].score==64
assert record.alignments[12].hsps[0].bits==29.3
assert record.alignments[12].hsps[0].expect==7.9
assert len(record.alignments[12].hsps)==1
assert record.alignments[0].hsps[0].identities==(247, 270)
assert record.alignments[0].hsps[0].positives==(247, 270)
assert record.alignments[1].hsps[0].identities==(123, 259)
assert record.alignments[1].hsps[0].positives==(155, 259)
assert record.alignments[1].hsps[0].gaps==(13, 259)
assert record.alignments[2].hsps[0].identities==(93, 253)
assert record.alignments[2].hsps[0].positives==(132, 253)
assert record.alignments[2].hsps[0].gaps==(8, 253)
assert record.alignments[2].hsps[1].identities==(33, 99)
assert record.alignments[2].hsps[1].positives==(47, 99)
assert record.alignments[2].hsps[1].gaps==(7, 99)
assert record.alignments[3].hsps[0].identities==(76, 247)
assert record.alignments[3].hsps[0].positives==(114, 247)
assert record.alignments[3].hsps[0].gaps==(17, 247)
assert record.alignments[4].hsps[0].identities==(46, 170)
assert record.alignments[4].hsps[0].positives==(76, 170)
assert record.alignments[4].hsps[0].gaps==(3, 170)
assert record.alignments[5].hsps[0].identities==(25, 64)
assert record.alignments[5].hsps[0].positives==(43, 64)
assert record.alignments[5].hsps[1].identities==(19, 63)
assert record.alignments[5].hsps[1].positives==(36, 63)
assert record.alignments[6].hsps[0].identities==(24, 64)
assert record.alignments[6].hsps[0].positives==(43, 64)
assert record.alignments[6].hsps[1].identities==(21, 63)
assert record.alignments[6].hsps[1].positives==(36, 63)
assert record.alignments[7].hsps[0].identities==(24, 64)
assert record.alignments[7].hsps[0].positives==(43, 64)
assert record.alignments[7].hsps[1].identities==(20, 63)
assert record.alignments[7].hsps[1].positives==(37, 63)
assert record.alignments[8].hsps[0].identities==(19, 67)
assert record.alignments[8].hsps[0].positives==(46, 67)
assert record.alignments[9].hsps[0].identities==(23, 62)
assert record.alignments[9].hsps[0].positives==(37, 62)
assert record.alignments[9].hsps[0].gaps==(1, 62)
assert record.alignments[10].hsps[0].identities==(41, 143)
assert record.alignments[10].hsps[0].positives==(62, 143)
assert record.alignments[10].hsps[0].gaps==(12, 143)
assert record.alignments[11].hsps[0].identities==(32, 141)
assert record.alignments[11].hsps[0].positives==(61, 141)
assert record.alignments[11].hsps[0].gaps==(6, 141)
assert record.alignments[12].hsps[0].identities==(21, 62)
assert record.alignments[12].hsps[0].positives==(29, 62)
assert record.alignments[0].hsps[0].query=="MAATKQGGGDDGRCARGVVLERTGARMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAYGAIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILALP"
assert record.alignments[0].hsps[0].match=="MAATKQGGGDDGRCARGVVLERTGARMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNN                       AGEKLIAAYGAIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILALP"
assert record.alignments[0].hsps[0].sbjct=="MAATKQGGGDDGRCARGVVLERTGARMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNLLAAPVVTAAPGGKAGGGAVLTPAGEKLIAAYGAIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILALP"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==270
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==270
assert record.alignments[1].hsps[0].query=="LERTGA-RMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAYGAIEAGVAKLL-------SSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILALP"
assert record.alignments[1].hsps[0].match=="L+R GA R+G +R+ LL AI R G+I+ AAREVGLSYK AWD V  +NN                       AG+ LIA +G +E  + K L       S+ EK+LN      L  L++RTS RN   C V  V    V A+V + L +G  LTAVIT RSA EM LAPG EV AL+K++FV+LA  G P R+S  NR+ G V  R D P+++E+ LDLG  K+ITA IT  SA+ L L PGV  TAL K+SHVILA+P"
assert record.alignments[1].hsps[0].sbjct=="LQRAGAPRVGGDRIRLLEAIARHGTIAGAAREVGLSYKTAWDAVGTLNNLFEQPLVEAAPGGRTGGNARVTEAGQALIAGFGLLEGALTKALGVLEGGVSAPEKALN-----TLWSLTMRTSNRNTLRCTVTRVTLGAVNAEVELALTDGHSLTAVITERSATEMGLAPGVEVFALIKASFVMLAAGGDPGRISACNRLTGIVAARTDGPVNTEIILDLGNCKSITAVITHTSADALGLAPGVPATALFKASHVILAMP"
assert record.alignments[1].hsps[0].query_start==20
assert record.alignments[1].hsps[0].query_end==270
assert record.alignments[1].hsps[0].sbjct_start==12
assert record.alignments[1].hsps[0].sbjct_end==265
assert record.alignments[2].hsps[0].query=="GARMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAYGAIEAGVAKLLSSFEKSLNLDPA-------EVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILAL"
assert record.alignments[2].hsps[0].match=="G  +   R+ LL AI R GSI+ AA+ V LSYKAAWD +  MNN                        G +++A Y A+E      L    + LN            ++  +S++TSARN +A  V  +    V  +VR+RL    ++ AVIT  SA  + LA G EV ALVKS+ V+L       +L+ RN++ G VI+  + P+++EVTL L  G+++T  +T DS + L L PGV   A  KSS VILA+"
assert record.alignments[2].hsps[0].sbjct=="GTALSDTRIRLLEAIEREGSINRAAKVVPLSYKAAWDAIDTMNNLAPEPLVVRVAGGRQGGGTQLTDYGRRIVAMYRALEIEYQSALDRLSERLNEVTGGDIQAFQRLMHSMSMKTSARNQFAGIVTGLRVGGVDYEVRIRLDAENEIAAVITKASAENLELAIGKEVFALVKSSSVMLT-TEPSLKLTARNQLWGEVIDIHEGPVNNEVTLALPSGRSVTCVVTADSCKALGLAPGVAACAFFKSSSVILAV"
assert record.alignments[2].hsps[0].query_start==24
assert record.alignments[2].hsps[0].query_end==269
assert record.alignments[2].hsps[0].sbjct_start==17
assert record.alignments[2].hsps[0].sbjct_end==268
assert record.alignments[2].hsps[1].query=="AIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAG"
assert record.alignments[2].hsps[1].match=="AI   V  L+ S    L  +P       SL+ +ARN    +V  +    V  +V + L  G+ +T V+TA S   + LAPG    A  KS+ V+LA  G"
assert record.alignments[2].hsps[1].sbjct=="AIGKEVFALVKSSSVMLTTEP-------SLKLTARNQLWGEVIDIHEGPVNNEVTLALPSGRSVTCVVTADSCKALGLAPGVAACAFFKSSSVILAVYG"
assert record.alignments[2].hsps[1].query_start==101
assert record.alignments[2].hsps[1].query_end==199
assert record.alignments[2].hsps[1].sbjct_start==179
assert record.alignments[2].hsps[1].sbjct_end==270
assert record.alignments[3].hsps[0].query=="RVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAY---GAIEAGVAKLLSSFEK-SLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQD-LTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVR----NRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILA"
assert record.alignments[3].hsps[0].match=="R++LL  I  +GSIS  A++ G+SYK+AWD +  MN                         G++LI  Y     I+     +LS  +   LN   A + R  SL+TSARN W   + +   DDV   V + L +G+  L   ITA+S A + L  G EVL L+K+ +V     G+ +  +V     N++ G +          EV + L  G+T+ AT+  +  E   L  G   TA   +  VI+A"
assert record.alignments[3].hsps[0].sbjct=="RISLLKHIALSGSISQGAKDAGISYKSAWDAINEMNQLSEHILVERATGGKGGGGAVLTRYGQRLIQLYDLLAQIQQKAFDVLSDDDALPLNSLLAAISR-FSLQTSARNQWFGTITARDHDDVQQHVDVLLADGKTRLKVAITAQSGARLGLDEGKEVLILLKAPWV-----GITQDEAVAQNADNQLPGIISHIERGAEQCEVLMALPDGQTLCATVPVN--EATSLQQGQNVTAYFNADSVIIA"
assert record.alignments[3].hsps[0].query_start==31
assert record.alignments[3].hsps[0].query_end==268
assert record.alignments[3].hsps[0].sbjct_start==21
assert record.alignments[3].hsps[0].sbjct_end==259
assert record.alignments[4].hsps[0].query=="ERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAYGAIEAGVAKLLSSF-EKSLNLDP-AEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRL-GEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLA"
assert record.alignments[4].hsps[0].match=="+RV LL  I + GSI+ AA+   +SYK+AWD ++AMN                          E+L+  Y  +E           ++S+ LD         SL++SARN +  +V      D    V + + G    L   IT +S+A ++L    EV+ + K+ +V ++"
assert record.alignments[4].hsps[0].sbjct=="KRVRLLKEIQQCGSINQAAKNAKVSYKSAWDHLEAMNKISPRPLLERNTGGKNGGGTALTTYAERLLQLYDLLERTQEHAFHILQDESVPLDSLLTATARFSLQSSARNQFFGRVAQQRIIDSRCVVDVNVQGLPTPLQVSITTKSSARLKLITEKEVMLMFKAPWVKIS"
assert record.alignments[4].hsps[0].query_start==30
assert record.alignments[4].hsps[0].query_end==196
assert record.alignments[4].hsps[0].sbjct_start==21
assert record.alignments[4].hsps[0].sbjct_end==190
assert record.alignments[5].hsps[0].query=="LSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVIL"
assert record.alignments[5].hsps[0].match=="+S RN+++G+V+      +++EV L++ GG  IT+ I+ DS E L +  G E TA+IKS+ V++"
assert record.alignments[5].hsps[0].sbjct=="ISARNQLKGKVVGLKKGVITAEVVLEIAGGNKITSIISLDSVEELGVKEGAELTAVIKSTDVMI"
assert record.alignments[5].hsps[0].query_start==204
assert record.alignments[5].hsps[0].query_end==267
assert record.alignments[5].hsps[0].sbjct_start==3
assert record.alignments[5].hsps[0].sbjct_end==66
assert record.alignments[5].hsps[1].query=="SARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLL"
assert record.alignments[5].hsps[1].match=="SARN    KV  +    + A+V + +  G  +T++I+  S  E+ +  G+E+ A++KS  V++"
assert record.alignments[5].hsps[1].sbjct=="SARNQLKGKVVGLKKGVITAEVVLEIAGGNKITSIISLDSVEELGVKEGAELTAVIKSTDVMI"
assert record.alignments[5].hsps[1].query_start==133
assert record.alignments[5].hsps[1].query_end==195
assert record.alignments[5].hsps[1].sbjct_start==4
assert record.alignments[5].hsps[1].sbjct_end==66
assert record.alignments[6].hsps[0].query=="LSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVIL"
assert record.alignments[6].hsps[0].match=="+S RN+++G+V+      +++EV L++ GG  IT+ I+ DS E L +  G E TA++KS+ V++"
assert record.alignments[6].hsps[0].sbjct=="ISARNQLKGKVVGLKKGVVTAEVVLEIAGGNKITSIISLDSVEELGVKEGAELTAVVKSTDVMI"
assert record.alignments[6].hsps[0].query_start==204
assert record.alignments[6].hsps[0].query_end==267
assert record.alignments[6].hsps[0].sbjct_start==3
assert record.alignments[6].hsps[0].sbjct_end==66
assert record.alignments[6].hsps[1].query=="SARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLL"
assert record.alignments[6].hsps[1].match=="SARN    KV  +    V A+V + +  G  +T++I+  S  E+ +  G+E+ A+VKS  V++"
assert record.alignments[6].hsps[1].sbjct=="SARNQLKGKVVGLKKGVVTAEVVLEIAGGNKITSIISLDSVEELGVKEGAELTAVVKSTDVMI"
assert record.alignments[6].hsps[1].query_start==133
assert record.alignments[6].hsps[1].query_end==195
assert record.alignments[6].hsps[1].sbjct_start==4
assert record.alignments[6].hsps[1].sbjct_end==66
assert record.alignments[7].hsps[0].query=="LSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVIL"
assert record.alignments[7].hsps[0].match=="+S RN+++G+V+      +++EV L++ GG  +T+ I+ DS E L +  G E TA+IKS+ V++"
assert record.alignments[7].hsps[0].sbjct=="ISARNQLKGKVVAVKKGLVTAEVVLEIAGGDKVTSIISLDSIEDLGVKEGTELTAVIKSTDVMI"
assert record.alignments[7].hsps[0].query_start==204
assert record.alignments[7].hsps[0].query_end==267
assert record.alignments[7].hsps[0].sbjct_start==3
assert record.alignments[7].hsps[0].sbjct_end==66
assert record.alignments[7].hsps[1].query=="SARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLL"
assert record.alignments[7].hsps[1].match=="SARN    KV +V    V A+V + +  G  +T++I+  S  ++ +  G+E+ A++KS  V++"
assert record.alignments[7].hsps[1].sbjct=="SARNQLKGKVVAVKKGLVTAEVVLEIAGGDKVTSIISLDSIEDLGVKEGTELTAVIKSTDVMI"
assert record.alignments[7].hsps[1].query_start==133
assert record.alignments[7].hsps[1].query_end==195
assert record.alignments[7].hsps[1].sbjct_start==4
assert record.alignments[7].hsps[1].sbjct_end==66
assert record.alignments[8].hsps[0].query=="RLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILAL"
assert record.alignments[8].hsps[0].match=="++S RN+++G+V+   +  +++ V +D+GGG  +++T++  + + L+L  G E  A+IK++ V++ +"
assert record.alignments[8].hsps[0].sbjct=="KISARNQLKGKVVSIENGSVNAIVHIDIGGGNVLSSTVSLAAVKELNLEVGKEAYAIIKATSVMVGV"
assert record.alignments[8].hsps[0].query_start==203
assert record.alignments[8].hsps[0].query_end==269
assert record.alignments[8].hsps[0].sbjct_start==2
assert record.alignments[8].hsps[0].sbjct_end==68
assert record.alignments[9].hsps[0].query=="PERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKS"
assert record.alignments[9].hsps[0].match=="P++ S+RN +R +V+   D     EV L++ GGKT+ A I+  + + L + PG+   A IKS"
assert record.alignments[9].hsps[0].sbjct=="PQQTSIRNVLRAKVVNSYDDNGQVEVELEV-GGKTLWARISPWARDELAIKPGLWLYAQIKS"
assert record.alignments[9].hsps[0].query_start==201
assert record.alignments[9].hsps[0].query_end==262
assert record.alignments[9].hsps[0].sbjct_start==287
assert record.alignments[9].hsps[0].sbjct_end==347
assert record.alignments[10].hsps[0].query=="EVLRGLSLRTSARNAWACKVWSVAA--DDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGG-GKTITATITRDSAEMLDLHPGVETTALIKS"
assert record.alignments[10].hsps[0].match=="+++  L L T+        + SV A  DD     R+    G    AV+ AR       APG  +   V +  V LA + + E  S+ N +   V E ++A   + V + L   G  + A ITR S + L + PG    A IK+"
assert record.alignments[10].hsps[0].sbjct=="DIMARLDLPTAFHEDAGVVIESVVAEHDDHYHLTRLAFPGG----AVLVARRPE----APGQRLRLRVHARDVSLANSRI-EDSSITNVLPATVREVVEADTPAHVLVRLEAEGTPLIARITRRSCDQLGIAPGRRMWAQIKA"
assert record.alignments[10].hsps[0].query_start==123
assert record.alignments[10].hsps[0].query_end==262
assert record.alignments[10].hsps[0].sbjct_start==242
assert record.alignments[10].hsps[0].sbjct_end==375
assert record.alignments[11].hsps[0].query=="AAYGAIEAGVAKLLSSFEKSLNLDP-AEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTI"
assert record.alignments[11].hsps[0].match=="AA G I+ GV  ++      + +DP AE ++    + +A  A   + W+   ++       + G   +L A I      +  L  G E + L ++ F+ +    +P      +  +  V+ER++       TLD+GG K +"
assert record.alignments[11].hsps[0].sbjct=="AATGTIQNGVTVIVDGINGDVIIDPSAETVKEYEEKHNAYLAQKAE-WAKLVNEPTVS---KDGHHVELAANIGTPDDVKGVLENGGEAVGLYRTEFLYMGRDQLPTEDEQFDAYK-TVLERMEGKSVVVRTLDIGGDKEL"
assert record.alignments[11].hsps[0].query_start==97
assert record.alignments[11].hsps[0].query_end==236
assert record.alignments[11].hsps[0].sbjct_start==207
assert record.alignments[11].hsps[0].sbjct_end==342
assert record.alignments[12].hsps[0].query=="AIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQ"
assert record.alignments[12].hsps[0].match=="A+E G     S+  +S   DPA V   +  R S  +      + VAAD   + VR +LG GQ"
assert record.alignments[12].hsps[0].sbjct=="AVELGGEIRFSTELQSFEQDPAGVTAVIKSRRSGEHTTVRADYLVAADGPRSPVREQLGIGQ"
assert record.alignments[12].hsps[0].query_start==101
assert record.alignments[12].hsps[0].query_end==162
assert record.alignments[12].hsps[0].sbjct_start==136
assert record.alignments[12].hsps[0].sbjct_end==197
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 29, 1999  4:39 PM',)]
assert record.ka_params==[0.316, 0.131, 0.361]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==12068104
assert record.num_sequences==82258
assert record.num_extends==396723
assert record.num_good_extends==1066
assert record.num_seqs_better_e==13
assert record.hsps_no_gap==10
assert record.hsps_prelim_gapped==3
assert record.hsps_gapped==18
assert record.query_length==270
assert record.database_length==29652561
assert record.effective_hsp_length==56
assert record.effective_query_length==214
assert record.effective_database_length==25046113
assert record.effective_search_space==5359868182
assert record.effective_search_space_used==5359868182
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.6)
assert record.blast_cutoff==(64,29.3)

handle = open('Blast/bt040')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.10'
assert record.date=="Aug-26-1999" 
assert record.reference==reference
assert record.query=="gi|730725|sp|Q05362|SCHB_STRHA SCHB PROTEIN"
assert record.query_letters==138
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==0
assert len(record.alignments)==0
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Nov 29, 1999  4:39 PM',)]
assert record.ka_params==[0.319, 0.139, 0.415]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==8952568
assert record.num_sequences==82258
assert record.num_extends==387403
assert record.num_good_extends==727
assert record.num_seqs_better_e==23
assert record.hsps_no_gap==13
assert record.hsps_prelim_gapped==10
assert record.hsps_gapped==23
assert record.query_length==138
assert record.database_length==29652561
assert record.effective_hsp_length==47
assert record.effective_query_length==91
assert record.effective_database_length==25786435
assert record.effective_search_space==2346565585
assert record.effective_search_space_used==2346565585
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.4)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(61,28.2)

handle = open('Blast/bt041')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000" 
assert record.reference==reference
assert record.query=="gi|120291|sp|P21297|FLBT_CAUCR FLBT PROTEIN"
assert record.query_letters==141
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==2
assert record.descriptions[0].title=="gi|120291|sp|P21297|FLBT_CAUCR FLBT PROTEIN"
assert record.descriptions[0].score==256
assert record.descriptions[0].e==2.0000000000000001e-68
assert record.descriptions[1].title=="gi|3024946|sp|Q58368|Y958_METJA HYPOTHETICAL PROTEIN MJ0958 PRE..."
assert record.descriptions[1].score==29
assert record.descriptions[1].e==3.400000
assert len(record.alignments)==2
assert record.alignments[0].title==">gi|120291|sp|P21297|FLBT_CAUCR FLBT PROTEIN"
assert record.alignments[0].length==141
assert record.alignments[1].title==">gi|3024946|sp|Q58368|Y958_METJA HYPOTHETICAL PROTEIN MJ0958 PRECURSOR"
assert record.alignments[1].length==426
assert record.alignments[0].hsps[0].score==646
assert record.alignments[0].hsps[0].bits==256
assert record.alignments[0].hsps[0].expect==2e-68
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==64
assert record.alignments[1].hsps[0].bits==29.3
assert record.alignments[1].hsps[0].expect==3.4
assert len(record.alignments[1].hsps)==1
assert record.alignments[0].hsps[0].identities==(127, 127)
assert record.alignments[0].hsps[0].positives==(127, 127)
assert record.alignments[1].hsps[0].identities==(15, 47)
assert record.alignments[1].hsps[0].positives==(23, 47)
assert record.alignments[0].hsps[0].query=="MPLKLSLKPGEKFVLNGAVVQNGDRRGVLVLQNKASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFYEEFATRLNEFMGVVRNPVVLQDCIAISKHVMAREYYKALMLSRKLIEYEDERL"
assert record.alignments[0].hsps[0].match=="MPLKLSLKPGEKFVLNGAVVQNGDRRGVLVLQNKASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFYEEFATRLNEFMGVVRNPVVLQDCIAISKHVMAREYYKALMLSRKLIEYEDERL"
assert record.alignments[0].hsps[0].sbjct=="MPLKLSLKPGEKFVLNGAVVQNGDRRGVLVLQNKASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFYEEFATRLNEFMGVVRNPVVLQDCIAISKHVMAREYYKALMLSRKLIEYEDERL"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==127
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==127
assert record.alignments[1].hsps[0].query=="VLVLQNKASVLREKDIMQPDQVTTPARHIYFPVMMMYLDEVGAEKFY"
assert record.alignments[1].hsps[0].match=="+LVL N  ++   K     D  TT   +IY P+ +    +  A+KFY"
assert record.alignments[1].hsps[0].sbjct=="ILVLINNTNITELKKFEDDDYYTTFQHYIYQPIFIFTTYDSKAKKFY"
assert record.alignments[1].hsps[0].query_start==28
assert record.alignments[1].hsps[0].query_end==74
assert record.alignments[1].hsps[0].sbjct_start==169
assert record.alignments[1].hsps[0].sbjct_end==215
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.ka_params==[0.323, 0.139, 0.392]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==6932673
assert record.num_sequences==82258
assert record.num_extends==246623
assert record.num_good_extends==486
assert record.num_seqs_better_e==2
assert record.hsps_no_gap==2
assert record.hsps_prelim_gapped==0
assert record.hsps_gapped==2
assert record.query_length==141
assert record.database_length==29652561
assert record.effective_hsp_length==50
assert record.effective_query_length==91
assert record.effective_database_length==25539661
assert record.effective_search_space==2324109151
assert record.effective_search_space_used==2324109151
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.5)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.9)
assert record.blast_cutoff==(61,28.2)

handle = open('Blast/bt042')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000" 
assert record.reference==reference
assert record.query=="gi|400206|sp|Q02112|LYTA_BACSU MEMBRANE-BOUND PROTEIN LYTA\nPRECURSOR"
assert record.query_letters==102
assert record.database=="data/pdbaa"
assert record.database_sequences==6999
assert record.database_letters==1461753
assert len(record.descriptions)==0
assert len(record.alignments)==0
assert record.database_name==['data/pdbaa']
assert record.num_letters_in_database==[1461753]
assert record.num_sequences_in_database==[6999]
assert record.posted_date==[('Feb 11, 2000  2:32 PM',)]
assert record.ka_params==[0.302, 0.126, 0.332]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==270485
assert record.num_sequences==6999
assert record.num_extends==9861
assert record.num_good_extends==7
assert record.num_seqs_better_e==0
assert record.hsps_no_gap==0
assert record.hsps_prelim_gapped==0
assert record.hsps_gapped==0
assert record.query_length==102
assert record.database_length==1461753
assert record.effective_hsp_length==47
assert record.effective_query_length==55
assert record.effective_database_length==1132800
assert record.effective_search_space==62304000
assert record.effective_search_space_used==62304000
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(17,7.4)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(43,21.7)
assert record.blast_cutoff==(47,22.7)

handle = open('Blast/bt043')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000" 
assert record.reference==reference
assert record.query=="gi|1718062|sp|P16153|UTXA_CLODI UTXA PROTEIN"
assert record.query_letters==166
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==0
assert len(record.alignments)==6
assert record.alignments[0].title==">gi|1718062|sp|P16153|UTXA_CLODI UTXA PROTEIN"
assert record.alignments[0].length==166
assert record.alignments[1].title==">gi|140528|sp|P24811|YQXH_BACSU HYPOTHETICAL 15.7 KD PROTEIN IN SPOIIIC-CWLA INTERGENIC REGION (ORF2)"
assert record.alignments[1].length==140
assert record.alignments[2].title==">gi|141088|sp|P26835|YNGD_CLOPE HYPOTHETICAL 14.9 KD PROTEIN IN NAGH 3'REGION (ORFD)"
assert record.alignments[2].length==132
assert record.alignments[3].title==">gi|6014830|sp|O78935|CYB_MARAM CYTOCHROME B"
assert record.alignments[3].length==379
assert record.alignments[4].title==">gi|1351589|sp|P47694|Y456_MYCGE HYPOTHETICAL PROTEIN MG456"
assert record.alignments[4].length==334
assert record.alignments[5].title==">gi|2496246|sp|Q57881|Y439_METJA HYPOTHETICAL ATP-BINDING PROTEIN MJ0439"
assert record.alignments[5].length==361
assert record.alignments[0].hsps[0].score==843
assert record.alignments[0].hsps[0].bits==332
assert record.alignments[0].hsps[0].expect==2e-91
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==90
assert record.alignments[1].hsps[0].bits==39.5
assert record.alignments[1].hsps[0].expect==0.004
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==88
assert record.alignments[2].hsps[0].bits==38.7
assert record.alignments[2].hsps[0].expect==0.007
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==64
assert record.alignments[3].hsps[0].bits==29.3
assert record.alignments[3].hsps[0].expect==4.6
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==63
assert record.alignments[4].hsps[0].bits==29.0
assert record.alignments[4].hsps[0].expect==6.0
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==62
assert record.alignments[5].hsps[0].bits==28.6
assert record.alignments[5].hsps[0].expect==7.8
assert len(record.alignments[5].hsps)==1
assert record.alignments[0].hsps[0].identities==(166, 166)
assert record.alignments[0].hsps[0].positives==(166, 166)
assert record.alignments[1].hsps[0].identities==(27, 130)
assert record.alignments[1].hsps[0].positives==(55, 130)
assert record.alignments[1].hsps[0].gaps==(19, 130)
assert record.alignments[2].hsps[0].identities==(24, 110)
assert record.alignments[2].hsps[0].positives==(52, 110)
assert record.alignments[2].hsps[0].gaps==(18, 110)
assert record.alignments[3].hsps[0].identities==(19, 57)
assert record.alignments[3].hsps[0].positives==(33, 57)
assert record.alignments[3].hsps[0].gaps==(2, 57)
assert record.alignments[4].hsps[0].identities==(16, 44)
assert record.alignments[4].hsps[0].positives==(24, 44)
assert record.alignments[4].hsps[0].gaps==(2, 44)
assert record.alignments[5].hsps[0].identities==(19, 56)
assert record.alignments[5].hsps[0].positives==(30, 56)
assert record.alignments[5].hsps[0].gaps==(12, 56)
assert record.alignments[0].hsps[0].query=="MHSSSPFYISNGNKIFFYINLGGVMNMTISFLSEHIFIKLVILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCLCGLPVPKRLKEKIAILLDAMTDEMNAKDEK"
assert record.alignments[0].hsps[0].match=="MHSSSPFYISNGNKIFFYINLGGVMNMTISFLSEHIFIKLVILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCLCGLPVPKRLKEKIAILLDAMTDEMNAKDEK"
assert record.alignments[0].hsps[0].sbjct=="MHSSSPFYISNGNKIFFYINLGGVMNMTISFLSEHIFIKLVILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCLCGLPVPKRLKEKIAILLDAMTDEMNAKDEK"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==166
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==166
assert record.alignments[1].hsps[0].query=="FIKLVILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILF-ILYESVSILKNMCLCGLPVPKRLKEKIAILLDAMTDEMNAKDE"
assert record.alignments[1].hsps[0].match=="++ L+++    D L G + A K +K  S     G +RK+     +   +V+D +   N                  G+  F ++LF I  E +SI +N+   G+ +P  + +++  + +      N  D+"
assert record.alignments[1].hsps[0].sbjct=="YLDLLLVLSIIDVLTGVIKAWKFKKLRSRSAWFGYVRKLLNFFAVILANVIDTVLNLN------------------GVLTFGTVLFYIANEGLSITENLAQIGVKIPSSITDRLQTIENEKEQSKNNADK"
assert record.alignments[1].hsps[0].query_start==37
assert record.alignments[1].hsps[0].query_end==165
assert record.alignments[1].hsps[0].sbjct_start==26
assert record.alignments[1].hsps[0].sbjct_end==137
assert record.alignments[2].hsps[0].query=="VILTISFDTLLGCLSAIKSRKFNSSFGIDGGIRKVAMIACIFFLSVVD-ILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCLCGLPVPKRLKEKI"
assert record.alignments[2].hsps[0].match=="+++ I  D L G +   KS++  S+ G+ G  +K  ++  +    ++D +L    ++F                     +  +I+ E +SIL+N    G+P+P++LK+ +"
assert record.alignments[2].hsps[0].sbjct=="LLVFIFLDYLTGVIKGCKSKELCSNIGLRGITKKGLILVVLLVAVMLDRLLDNGTWMFRT-----------------LIAYFYIMNEGISILENCAALGVPIPEKLKQAL"
assert record.alignments[2].hsps[0].query_start==41
assert record.alignments[2].hsps[0].query_end==149
assert record.alignments[2].hsps[0].sbjct_start==33
assert record.alignments[2].hsps[0].sbjct_end==125
assert record.alignments[3].hsps[0].query=="CIFFLSVVDILTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKNMCL"
assert record.alignments[3].hsps[0].match=="C+F+L V D+LT   ++   P +   F+ +  L    +F+IL IL  ++SI++N  L"
assert record.alignments[3].hsps[0].sbjct=="CLFWLLVADLLT-LTWIGGQPVEH-PFITIGQLASILYFAILLILMPAISIIENNLL"
assert record.alignments[3].hsps[0].query_start==80
assert record.alignments[3].hsps[0].query_end==136
assert record.alignments[3].hsps[0].sbjct_start==323
assert record.alignments[3].hsps[0].sbjct_end==377
assert record.alignments[4].hsps[0].query=="LTKFNFLFMLPQDCINFLRLKHLGISEFFSILFILYESVSILKN"
assert record.alignments[4].hsps[0].match=="LTKFN  F+ P     FLR+  +G+   FS++ I +   S  +N"
assert record.alignments[4].hsps[0].sbjct=="LTKFNKFFLTPNKLNAFLRV--IGLCGLFSVIAISFGIYSYTRN"
assert record.alignments[4].hsps[0].query_start==90
assert record.alignments[4].hsps[0].query_end==133
assert record.alignments[4].hsps[0].sbjct_start==4
assert record.alignments[4].hsps[0].sbjct_end==45
assert record.alignments[5].hsps[0].query=="FLRLKHLGIS---EFFSILFILYES----VSILKNMC-----LCGLPVPKRLKEKI"
assert record.alignments[5].hsps[0].match=="++ L+ + IS   +F  +LF  YE     V I+K++      LCG+P PK   E+I"
assert record.alignments[5].hsps[0].sbjct=="YINLRGIFISKYKDFIEVLFEEYEEDRKPVEIIKSLIKDVPSLCGIPTPKNTLEEI"
assert record.alignments[5].hsps[0].query_start==106
assert record.alignments[5].hsps[0].query_end==149
assert record.alignments[5].hsps[0].sbjct_start==68
assert record.alignments[5].hsps[0].sbjct_end==123
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.ka_params==[0.331, 0.146, 0.428]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==8801581
assert record.num_sequences==82258
assert record.num_extends==320828
assert record.num_good_extends==892
assert record.num_seqs_better_e==6
assert record.hsps_no_gap==3
assert record.hsps_prelim_gapped==3
assert record.hsps_gapped==6
assert record.query_length==166
assert record.database_length==29652561
assert record.effective_hsp_length==46
assert record.effective_query_length==120
assert record.effective_database_length==25868693
assert record.effective_search_space==3104243160
assert record.effective_search_space_used==3104243160
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(15,7.2)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(40,21.9)
assert record.blast_cutoff==(62,28.6)

handle = open('Blast/bt044')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000" 
assert record.reference==reference
assert record.query=="gi|1718062|sp|P16153|UTXA_CLODI UTXA PROTEIN"
assert record.query_letters==166
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==6
assert record.descriptions[0].title=="gi|1718062|sp|P16153|UTXA_CLODI UTXA PROTEIN"
assert record.descriptions[0].score==332
assert record.descriptions[0].e==record.descriptions[0].e
assert record.descriptions[1].title=="gi|140528|sp|P24811|YQXH_BACSU HYPOTHETICAL 15.7 KD PROTEIN IN ..."
assert record.descriptions[1].score==39
assert record.descriptions[1].e==0.004000
assert record.descriptions[2].title=="gi|141088|sp|P26835|YNGD_CLOPE HYPOTHETICAL 14.9 KD PROTEIN IN ..."
assert record.descriptions[2].score==39
assert record.descriptions[2].e==0.007000
assert record.descriptions[3].title=="gi|6014830|sp|O78935|CYB_MARAM CYTOCHROME B"
assert record.descriptions[3].score==29
assert record.descriptions[3].e==4.600000
assert record.descriptions[4].title=="gi|1351589|sp|P47694|Y456_MYCGE HYPOTHETICAL PROTEIN MG456"
assert record.descriptions[4].score==29
assert record.descriptions[4].e==6.000000
assert record.descriptions[5].title=="gi|2496246|sp|Q57881|Y439_METJA HYPOTHETICAL ATP-BINDING PROTEI..."
assert record.descriptions[5].score==29
assert record.descriptions[5].e==7.800000
assert len(record.alignments)==0
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.ka_params==[0.331, 0.146, 0.428]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==8801581
assert record.num_sequences==82258
assert record.num_extends==320828
assert record.num_good_extends==892
assert record.num_seqs_better_e==6
assert record.hsps_no_gap==3
assert record.hsps_prelim_gapped==3
assert record.hsps_gapped==6
assert record.query_length==166
assert record.database_length==29652561
assert record.effective_hsp_length==46
assert record.effective_query_length==120
assert record.effective_database_length==25868693
assert record.effective_search_space==3104243160
assert record.effective_search_space_used==3104243160
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(15,7.2)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(40,21.9)
assert record.blast_cutoff==(62,28.6)

handle = open('Blast/bt045')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000" 
assert record.reference==reference
assert record.query=="gi|132349|sp|P15394|REPA_AGRTU REPLICATING PROTEIN"
assert record.query_letters==250
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==15
assert record.descriptions[0].title=="gi|132349|sp|P15394|REPA_AGRTU REPLICATING PROTEIN"
assert record.descriptions[0].score==514
assert record.descriptions[0].e==1e-146
assert record.descriptions[1].title=="gi|123932|sp|P19922|HYIN_BRAJA INDOLEACETAMIDE HYDROLASE (IAH) ..."
assert record.descriptions[1].score==34
assert record.descriptions[1].e==0.290000
assert record.descriptions[2].title=="gi|137670|sp|P06422|VE2_HPV08 REGULATORY PROTEIN E2"
assert record.descriptions[2].score==32
assert record.descriptions[2].e==0.860000
assert record.descriptions[3].title=="gi|5921693|sp|Q05152|CCAB_RABIT VOLTAGE-DEPENDENT N-TYPE CALCIU..."
assert record.descriptions[3].score==32
assert record.descriptions[3].e==1.500000
assert record.descriptions[4].title=="gi|121952|sp|P02256|H1_PARAN HISTONE H1, GONADAL"
assert record.descriptions[4].score==31
assert record.descriptions[4].e==2.500000
assert record.descriptions[5].title=="gi|3915729|sp|P51592|HYDP_DROME HYPERPLASTIC DISCS PROTEIN (HYD..."
assert record.descriptions[5].score==31
assert record.descriptions[5].e==3.300000
assert record.descriptions[6].title=="gi|124141|sp|P08392|ICP4_HSV11 TRANS-ACTING TRANSCRIPTIONAL PRO..."
assert record.descriptions[6].score==31
assert record.descriptions[6].e==3.300000
assert record.descriptions[7].title=="gi|462182|sp|P33438|GLT_DROME GLUTACTIN PRECURSOR"
assert record.descriptions[7].score==31
assert record.descriptions[7].e==3.300000
assert record.descriptions[8].title=="gi|1708851|sp|P55268|LMB2_HUMAN LAMININ BETA-2 CHAIN PRECURSOR ..."
assert record.descriptions[8].score==30
assert record.descriptions[8].e==4.300000
assert record.descriptions[9].title=="gi|731294|sp|P39713|YAG1_YEAST HYPOTHETICAL ZINC-TYPE ALCOHOL D..."
assert record.descriptions[9].score==30
assert record.descriptions[9].e==4.300000
assert record.descriptions[10].title=="gi|2495137|sp|Q24704|H1_DROVI HISTONE H1"
assert record.descriptions[10].score==29
assert record.descriptions[10].e==7.500000
assert record.descriptions[11].title=="gi|6226905|sp|Q59567|TOP1_MYCTU DNA TOPOISOMERASE I (OMEGA-PROT..."
assert record.descriptions[11].score==29
assert record.descriptions[11].e==9.800000
assert record.descriptions[12].title=="gi|6093970|sp|Q61085|RHOP_MOUSE GTP-RHO BINDING PROTEIN 1 (RHOP..."
assert record.descriptions[12].score==29
assert record.descriptions[12].e==9.800000
assert record.descriptions[13].title=="gi|1172635|sp|P46466|PRS4_ORYSA 26S PROTEASE REGULATORY SUBUNIT..."
assert record.descriptions[13].score==29
assert record.descriptions[13].e==9.800000
assert record.descriptions[14].title=="gi|547963|sp|Q01989|MYS9_DROME MYOSIN HEAVY CHAIN 95F (95F MHC)"
assert record.descriptions[14].score==29
assert record.descriptions[14].e==9.800000
assert len(record.alignments)==0
assert record.database_name==['data/swissprot']
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.ka_params==[0.317, 0.133, 0.395]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==14679054
assert record.num_sequences==82258
assert record.num_extends==599117
assert record.num_good_extends==1508
assert record.num_seqs_better_e==15
assert record.hsps_no_gap==4
assert record.hsps_prelim_gapped==11
assert record.hsps_gapped==15
assert record.query_length==250
assert record.database_length==29652561
assert record.effective_hsp_length==51
assert record.effective_query_length==199
assert record.effective_database_length==25457403
assert record.effective_search_space==5066023197
assert record.effective_search_space_used==5066023197
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(63,29.0)

handle = open('Blast/bt046')
record = pb_parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000" 
assert record.reference==reference
assert record.query=="gi|1174726|sp|P12921|TMRB_BACSU TUNICAMYCIN RESISTANCE PROTEIN"
assert record.query_letters==197
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.rounds)==2
assert len(record.rounds[0].new_seqs)==4
assert record.rounds[0].new_seqs[0].title=="gi|1174726|sp|P12921|TMRB_BACSU TUNICAMYCIN RESISTANCE PROTEIN"
assert record.rounds[0].new_seqs[0].score==402
assert record.rounds[0].new_seqs[0].e==9.9999999999999995e-113
assert record.rounds[0].new_seqs[1].title=="gi|1352934|sp|P47169|YJ9F_YEAST HYPOTHETICAL 161.2 KD PROTEIN I..."
assert record.rounds[0].new_seqs[1].score==30
assert record.rounds[0].new_seqs[1].e==3.300000
assert record.rounds[0].new_seqs[2].title=="gi|3915741|sp|P04407|KITH_HSV23 THYMIDINE KINASE"
assert record.rounds[0].new_seqs[2].score==29
assert record.rounds[0].new_seqs[2].e==7.400000
assert record.rounds[0].new_seqs[3].title=="gi|3334489|sp|P15398|RPA1_SCHPO DNA-DIRECTED RNA POLYMERASE I 1..."
assert record.rounds[0].new_seqs[3].score==29
assert record.rounds[0].new_seqs[3].e==7.400000
assert len(record.rounds[0].alignments)==4
assert record.rounds[0].alignments[0].title==">gi|1174726|sp|P12921|TMRB_BACSU TUNICAMYCIN RESISTANCE PROTEIN"
assert record.rounds[0].alignments[0].length==197
assert record.rounds[0].alignments[1].title==">gi|1352934|sp|P47169|YJ9F_YEAST HYPOTHETICAL 161.2 KD PROTEIN IN NMD5-HOM6 INTERGENIC REGION"
assert record.rounds[0].alignments[1].length==1442
assert record.rounds[0].alignments[2].title==">gi|3915741|sp|P04407|KITH_HSV23 THYMIDINE KINASE"
assert record.rounds[0].alignments[2].length==376
assert record.rounds[0].alignments[3].title==">gi|3334489|sp|P15398|RPA1_SCHPO DNA-DIRECTED RNA POLYMERASE I 190 KD POLYPEPTIDE"
assert record.rounds[0].alignments[3].length==1689
assert len(record.rounds[1].new_seqs)==3
assert record.rounds[1].new_seqs[0].title=="gi|1352934|sp|P47169|YJ9F_YEAST HYPOTHETICAL 161.2 KD PROTEIN I..."
assert record.rounds[1].new_seqs[0].score==30
assert record.rounds[1].new_seqs[0].e==3.3
assert record.rounds[1].new_seqs[1].title=="gi|3334489|sp|P15398|RPA1_SCHPO DNA-DIRECTED RNA POLYMERASE I 1..."
assert record.rounds[1].new_seqs[1].score==30
assert record.rounds[1].new_seqs[1].e==3.3
assert record.rounds[1].new_seqs[2].title=="gi|3915741|sp|P04407|KITH_HSV23 THYMIDINE KINASE"
assert record.rounds[1].new_seqs[2].score==29
assert record.rounds[1].new_seqs[2].e==7.4
assert len(record.rounds[1].alignments)==4
assert record.rounds[1].alignments[0].title==">gi|1174726|sp|P12921|TMRB_BACSU TUNICAMYCIN RESISTANCE PROTEIN"
assert record.rounds[1].alignments[0].length==197
assert record.rounds[1].alignments[1].title==">gi|1352934|sp|P47169|YJ9F_YEAST HYPOTHETICAL 161.2 KD PROTEIN IN NMD5-HOM6 INTERGENIC REGION"
assert record.rounds[1].alignments[1].length==1442
assert record.rounds[1].alignments[2].title==">gi|3334489|sp|P15398|RPA1_SCHPO DNA-DIRECTED RNA POLYMERASE I 190 KD POLYPEPTIDE"
assert record.rounds[1].alignments[2].length==1689
assert record.rounds[1].alignments[3].title==">gi|3915741|sp|P04407|KITH_HSV23 THYMIDINE KINASE"
assert record.rounds[1].alignments[3].length==376
assert record.rounds[0].alignments[0].hsps[0].score==1021
assert record.rounds[0].alignments[0].hsps[0].bits==402
assert record.rounds[0].alignments[0].hsps[0].expect==1e-112
assert len(record.rounds[0].alignments[0].hsps)==1
assert record.rounds[0].alignments[1].hsps[0].score==66
assert record.rounds[0].alignments[1].hsps[0].bits==30.1
assert record.rounds[0].alignments[1].hsps[0].expect==3.3
assert len(record.rounds[0].alignments[1].hsps)==1
assert record.rounds[0].alignments[2].hsps[0].score==63
assert record.rounds[0].alignments[2].hsps[0].bits==29.0
assert record.rounds[0].alignments[2].hsps[0].expect==7.4
assert len(record.rounds[0].alignments[2].hsps)==1
assert record.rounds[0].alignments[3].hsps[0].score==63
assert record.rounds[0].alignments[3].hsps[0].bits==29.0
assert record.rounds[0].alignments[3].hsps[0].expect==7.4
assert record.rounds[1].alignments[0].hsps[0].score==1031
assert record.rounds[1].alignments[0].hsps[0].bits==406
assert record.rounds[1].alignments[0].hsps[0].expect==1e-113
assert len(record.rounds[1].alignments[0].hsps)==1
assert record.rounds[1].alignments[1].hsps[0].score==66
assert record.rounds[1].alignments[1].hsps[0].bits==30.1
assert record.rounds[1].alignments[1].hsps[0].expect==3.3
assert len(record.rounds[1].alignments[1].hsps)==1
assert record.rounds[1].alignments[2].hsps[0].score==66
assert record.rounds[1].alignments[2].hsps[0].bits==30.1
assert record.rounds[1].alignments[2].hsps[0].expect==3.3
assert len(record.rounds[1].alignments[2].hsps)==1
assert record.rounds[1].alignments[3].hsps[0].score==63
assert record.rounds[1].alignments[3].hsps[0].bits==28.9
assert record.rounds[1].alignments[3].hsps[0].expect==7.4
assert len(record.rounds[1].alignments[3].hsps)==1
assert record.rounds[0].alignments[0].hsps[0].identities==(197, 197)
assert record.rounds[0].alignments[0].hsps[0].positives==(197, 197)
assert record.rounds[0].alignments[1].hsps[0].identities==(23, 70)
assert record.rounds[0].alignments[1].hsps[0].positives==(35, 70)
assert record.rounds[0].alignments[1].hsps[0].gaps==(10, 70)
assert record.rounds[0].alignments[2].hsps[0].identities==(15, 37)
assert record.rounds[0].alignments[2].hsps[0].positives==(22, 37)
assert record.rounds[0].alignments[2].hsps[0].gaps==(2, 37)
assert record.rounds[0].alignments[3].hsps[0].identities==(12, 38)
assert record.rounds[0].alignments[3].hsps[0].positives==(20, 38)
assert record.rounds[1].alignments[0].hsps[0].identities==(197, 197)
assert record.rounds[1].alignments[0].hsps[0].positives==(197, 197)
assert record.rounds[1].alignments[1].hsps[0].identities==(23, 70)
assert record.rounds[1].alignments[1].hsps[0].positives==(35, 70)
assert record.rounds[1].alignments[1].hsps[0].gaps==(10, 70)
assert record.rounds[1].alignments[2].hsps[0].identities==(12, 38)
assert record.rounds[1].alignments[2].hsps[0].positives==(20, 38)
assert record.rounds[1].alignments[3].hsps[0].identities==(15, 37)
assert record.rounds[1].alignments[3].hsps[0].positives==(22, 37)
assert record.rounds[1].alignments[3].hsps[0].gaps==(2, 37)
assert record.rounds[0].alignments[0].hsps[0].query=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[0].alignments[0].hsps[0].match=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[0].alignments[0].hsps[0].sbjct=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[0].alignments[0].hsps[0].query_start==1
assert record.rounds[0].alignments[0].hsps[0].query_end==197
assert record.rounds[0].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[0].hsps[0].sbjct_end==197
assert record.rounds[0].alignments[1].hsps[0].query=="TLMASKETLLKR--------LRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARA"
assert record.rounds[0].alignments[1].hsps[0].match=="TL+  K+  L R          TR + K S AA   D+ +EGLS P    +  +D  +  ++A+ +AARA"
assert record.rounds[0].alignments[1].hsps[0].sbjct=="TLLTRKDPSLSRNLKQSAGDALTRKQEKRSKAA--FDQLLEGLSGPPLHVYYASDGGNAANLAKRLAARA"
assert record.rounds[0].alignments[1].hsps[0].query_start==107
assert record.rounds[0].alignments[1].hsps[0].query_end==168
assert record.rounds[0].alignments[1].hsps[0].sbjct_start==637
assert record.rounds[0].alignments[1].hsps[0].sbjct_end==704
assert record.rounds[0].alignments[2].hsps[0].query=="IWINGAFGSGKTQTAFELHRRLNP--SYVYDPEKMGF"
assert record.rounds[0].alignments[2].hsps[0].match=="++I+G  G GKT T+ +L   L P  + VY PE M +"
assert record.rounds[0].alignments[2].hsps[0].sbjct=="VYIDGPHGVGKTTTSAQLMEALGPRDNIVYVPEPMTY"
assert record.rounds[0].alignments[2].hsps[0].query_start==3
assert record.rounds[0].alignments[2].hsps[0].query_end==37
assert record.rounds[0].alignments[2].hsps[0].sbjct_start==52
assert record.rounds[0].alignments[2].hsps[0].sbjct_end==88
assert record.rounds[0].alignments[3].hsps[0].query=="GILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASK"
assert record.rounds[0].alignments[3].hsps[0].match=="G +++P+   HP +F+++   LR      HHF L   K"
assert record.rounds[0].alignments[3].hsps[0].sbjct=="GHIVLPIPAYHPLFFSQMYNLLRSTCLYCHHFKLSKVK"
assert record.rounds[0].alignments[3].hsps[0].query_start==75
assert record.rounds[0].alignments[3].hsps[0].query_end==112
assert record.rounds[0].alignments[3].hsps[0].sbjct_start==78
assert record.rounds[0].alignments[3].hsps[0].sbjct_end==115
assert record.rounds[1].alignments[0].hsps[0].query=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[1].alignments[0].hsps[0].match=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[1].alignments[0].hsps[0].sbjct=="MIIWINGAFGSGKTQTAFELHRRLNPSYVYDPEKMGFALRSMVPQEIAKDDFQSYPLWRAFNYSLLASLTDTYRGILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASKETLLKRLRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARAELPLDPDTRGSLRRFADRLMVKLNHIRIK"
assert record.rounds[1].alignments[0].hsps[0].query_start==1
assert record.rounds[1].alignments[0].hsps[0].query_end==197
assert record.rounds[1].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[0].hsps[0].sbjct_end==197
assert record.rounds[1].alignments[1].hsps[0].query=="TLMASKETLLKR--------LRTRAEGKNSWAAKQIDRCVEGLSSPIFEDHIQTDNLSIQDVAENIAARA"
assert record.rounds[1].alignments[1].hsps[0].match=="TL+  K+  L R          TR + K S AA   D+ +EGLS P    +  +D  +  ++A+ +AARA"
assert record.rounds[1].alignments[1].hsps[0].sbjct=="TLLTRKDPSLSRNLKQSAGDALTRKQEKRSKAA--FDQLLEGLSGPPLHVYYASDGGNAANLAKRLAARA"
assert record.rounds[1].alignments[1].hsps[0].query_start==107
assert record.rounds[1].alignments[1].hsps[0].query_end==168
assert record.rounds[1].alignments[1].hsps[0].sbjct_start==637
assert record.rounds[1].alignments[1].hsps[0].sbjct_end==704
assert record.rounds[1].alignments[2].hsps[0].query=="GILIVPMTIVHPEYFNEIIGRLRQEGRIVHHFTLMASK"
assert record.rounds[1].alignments[2].hsps[0].match=="G +++P+   HP +F+++   LR      HHF L   K"
assert record.rounds[1].alignments[2].hsps[0].sbjct=="GHIVLPIPAYHPLFFSQMYNLLRSTCLYCHHFKLSKVK"
assert record.rounds[1].alignments[2].hsps[0].query_start==75
assert record.rounds[1].alignments[2].hsps[0].query_end==112
assert record.rounds[1].alignments[2].hsps[0].sbjct_start==78
assert record.rounds[1].alignments[2].hsps[0].sbjct_end==115
assert record.rounds[1].alignments[3].hsps[0].query=="IWINGAFGSGKTQTAFELHRRLNP--SYVYDPEKMGF"
assert record.rounds[1].alignments[3].hsps[0].match=="++I+G  G GKT T+ +L   L P  + VY PE M +"
assert record.rounds[1].alignments[3].hsps[0].sbjct=="VYIDGPHGVGKTTTSAQLMEALGPRDNIVYVPEPMTY"
assert record.rounds[1].alignments[3].hsps[0].query_start==3
assert record.rounds[1].alignments[3].hsps[0].query_end==37
assert record.rounds[1].alignments[3].hsps[0].sbjct_start==52
assert record.rounds[1].alignments[3].hsps[0].sbjct_end==88
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.ka_params==[0.318, 0.134, 0.412]
assert record.ka_params_gap==[0.270, 0.0471, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==23669955
assert record.num_sequences==82258
assert record.num_extends==918938
assert record.num_good_extends==2130
assert record.num_seqs_better_e==8
assert record.hsps_no_gap==4
assert record.hsps_prelim_gapped==4
assert record.hsps_gapped==8
assert record.query_length==197
assert record.database_length==29652561
assert record.effective_hsp_length==48
assert record.effective_query_length==149
assert record.effective_database_length==25704177
assert record.effective_search_space==3829922373
assert record.effective_search_space_used==3829922373
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(62,28.6)

handle = open('Blast/bt047')
record = pb_parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000" 
assert record.reference==reference
assert record.query=="gi|399896|sp|Q02134|HIS7_LACLA IMIDAZOLEGLYCEROL-PHOSPHATE\nDEHYDRATASE (IGPD)"
assert record.query_letters==200
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.rounds)==2
assert len(record.rounds[0].new_seqs)==30
assert record.rounds[0].new_seqs[0].title=="gi|399896|sp|Q02134|HIS7_LACLA IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[0].score==409
assert record.rounds[0].new_seqs[0].e==1.0000000000000001e-114
assert record.rounds[0].new_seqs[1].title=="gi|462273|sp|P34047|HIS7_ARATH IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[1].score==198
assert record.rounds[0].new_seqs[1].e==6e-51
assert record.rounds[0].new_seqs[2].title=="gi|2495230|sp|Q43072|HIS7_PEA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYD..."
assert record.rounds[0].new_seqs[2].score==196
assert record.rounds[0].new_seqs[2].e==4e-50
assert record.rounds[0].new_seqs[3].title=="gi|123157|sp|P18787|HIS7_AZOBR IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[3].score==185
assert record.rounds[0].new_seqs[3].e==5.0000000000000001e-47
assert record.rounds[0].new_seqs[4].title=="gi|462275|sp|P34048|HIS7_WHEAT IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[4].score==181
assert record.rounds[0].new_seqs[4].e==8.0000000000000002e-46
assert record.rounds[0].new_seqs[5].title=="gi|123161|sp|P16247|HIS7_STRCO IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[5].score==178
assert record.rounds[0].new_seqs[5].e==7e-45
assert record.rounds[0].new_seqs[6].title=="gi|462272|sp|Q05068|HIS7_ANASP IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[6].score==178
assert record.rounds[0].new_seqs[6].e==7e-45
assert record.rounds[0].new_seqs[7].title=="gi|123158|sp|P06987|HIS7_ECOLI HISTIDINE BIOSYNTHESIS BIFUNCTIO..."
assert record.rounds[0].new_seqs[7].score==175
assert record.rounds[0].new_seqs[7].e==7.9999999999999996e-44
assert record.rounds[0].new_seqs[8].title=="gi|1346293|sp|P48054|HIS7_SYNY3 IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[8].score==174
assert record.rounds[0].new_seqs[8].e==1.0000000000000001e-43
assert record.rounds[0].new_seqs[9].title=="gi|1170286|sp|P44327|HIS7_HAEIN HISTIDINE BIOSYNTHESIS BIFUNCTI..."
assert record.rounds[0].new_seqs[9].score==168
assert record.rounds[0].new_seqs[9].e==8.0000000000000003e-42
assert record.rounds[0].new_seqs[10].title=="gi|2495224|sp|O06590|HIS7_MYCTU IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[10].score==167
assert record.rounds[0].new_seqs[10].e==2e-41
assert record.rounds[0].new_seqs[11].title=="gi|123160|sp|P10368|HIS7_SALTY HISTIDINE BIOSYNTHESIS BIFUNCTIO..."
assert record.rounds[0].new_seqs[11].score==166
assert record.rounds[0].new_seqs[11].e==2e-41
assert record.rounds[0].new_seqs[12].title=="gi|2495226|sp|Q50504|HIS7_METTH PROBABLE IMIDAZOLEGLYCEROL-PHOS..."
assert record.rounds[0].new_seqs[12].score==153
assert record.rounds[0].new_seqs[12].e==3e-37
assert record.rounds[0].new_seqs[13].title=="gi|729718|sp|P40919|HIS7_CRYNE IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[13].score==152
assert record.rounds[0].new_seqs[13].e==7.0000000000000003e-37
assert record.rounds[0].new_seqs[14].title=="gi|3334215|sp|O33773|HIS7_SULSO PROBABLE IMIDAZOLEGLYCEROL-PHOS..."
assert record.rounds[0].new_seqs[14].score==151
assert record.rounds[0].new_seqs[14].e==9.0000000000000008e-37
assert record.rounds[0].new_seqs[15].title=="gi|123159|sp|P28624|HIS7_PHYPR IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[15].score==149
assert record.rounds[0].new_seqs[15].e==3.0000000000000002e-36
assert record.rounds[0].new_seqs[16].title=="gi|729719|sp|P40374|HIS7_SCHPO IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[16].score==136
assert record.rounds[0].new_seqs[16].e==3e-32
assert record.rounds[0].new_seqs[17].title=="gi|2495227|sp|P56090|HIS7_CANAL IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[17].score==128
assert record.rounds[0].new_seqs[17].e==8.9999999999999993e-30
assert record.rounds[0].new_seqs[18].title=="gi|2495225|sp|Q58109|HIS7_METJA PROBABLE IMIDAZOLEGLYCEROL-PHOS..."
assert record.rounds[0].new_seqs[18].score==126
assert record.rounds[0].new_seqs[18].e==3.9999999999999998e-29
assert record.rounds[0].new_seqs[19].title=="gi|2495229|sp|Q92447|HIS7_PICPA IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[19].score==125
assert record.rounds[0].new_seqs[19].e==6.0000000000000005e-29
assert record.rounds[0].new_seqs[20].title=="gi|399897|sp|Q02986|HIS7_SACKL IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[20].score==125
assert record.rounds[0].new_seqs[20].e==6.0000000000000005e-29
assert record.rounds[0].new_seqs[21].title=="gi|2495228|sp|Q12578|HIS7_CANGA IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[21].score==123
assert record.rounds[0].new_seqs[21].e==1.9999999999999999e-28
assert record.rounds[0].new_seqs[22].title=="gi|2506514|sp|P06633|HIS7_YEAST IMIDAZOLEGLYCEROL-PHOSPHATE DEH..."
assert record.rounds[0].new_seqs[22].score==122
assert record.rounds[0].new_seqs[22].e==3.9999999999999999e-28
assert record.rounds[0].new_seqs[23].title=="gi|462274|sp|P34041|HIS7_TRIHA IMIDAZOLEGLYCEROL-PHOSPHATE DEHY..."
assert record.rounds[0].new_seqs[23].score==106
assert record.rounds[0].new_seqs[23].e==3e-23
assert record.rounds[0].new_seqs[24].title=="gi|1345641|sp|P49264|C7B1_THLAR CYTOCHROME P450 71B1 (CYPLXXIB1)"
assert record.rounds[0].new_seqs[24].score==35
assert record.rounds[0].new_seqs[24].e==0.130000
assert record.rounds[0].new_seqs[25].title=="gi|1731346|sp|Q10698|YY29_MYCTU PROBABLE DIPEPTIDASE CY49.29C"
assert record.rounds[0].new_seqs[25].score==32
assert record.rounds[0].new_seqs[25].e==1.100000
assert record.rounds[0].new_seqs[26].title=="gi|3287839|sp|Q01812|GLK4_RAT GLUTAMATE RECEPTOR, IONOTROPIC KA..."
assert record.rounds[0].new_seqs[26].score==30
assert record.rounds[0].new_seqs[26].e==3.300000
assert record.rounds[0].new_seqs[27].title=="gi|3123025|sp|Q94637|VIT6_OSCBR VITELLOGENIN 6 PRECURSOR"
assert record.rounds[0].new_seqs[27].score==29
assert record.rounds[0].new_seqs[27].e==5.600000
assert record.rounds[0].new_seqs[28].title=="gi|1174406|sp|P36126|SP14_YEAST PHOSPHOLIPASE D1 (PLD 1) (CHOLI..."
assert record.rounds[0].new_seqs[28].score==29
assert record.rounds[0].new_seqs[28].e==9.700000
assert record.rounds[0].new_seqs[29].title=="gi|3287848|sp|Q16099|GLK4_HUMAN GLUTAMATE RECEPTOR, IONOTROPIC ..."
assert record.rounds[0].new_seqs[29].score==29
assert record.rounds[0].new_seqs[29].e==9.700000
assert len(record.rounds[0].alignments)==30
assert record.rounds[0].alignments[0].title==">gi|399896|sp|Q02134|HIS7_LACLA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[0].length==200
assert record.rounds[0].alignments[1].title==">gi|462273|sp|P34047|HIS7_ARATH IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[1].length==270
assert record.rounds[0].alignments[2].title==">gi|2495230|sp|Q43072|HIS7_PEA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[2].length==281
assert record.rounds[0].alignments[3].title==">gi|123157|sp|P18787|HIS7_AZOBR IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[3].length==207
assert record.rounds[0].alignments[4].title==">gi|462275|sp|P34048|HIS7_WHEAT IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[4].length==195
assert record.rounds[0].alignments[5].title==">gi|123161|sp|P16247|HIS7_STRCO IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[5].length==197
assert record.rounds[0].alignments[6].title==">gi|462272|sp|Q05068|HIS7_ANASP IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[6].length==209
assert record.rounds[0].alignments[7].title==">gi|123158|sp|P06987|HIS7_ECOLI HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[0].alignments[7].length==355
assert record.rounds[0].alignments[8].title==">gi|1346293|sp|P48054|HIS7_SYNY3 IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[8].length==210
assert record.rounds[0].alignments[9].title==">gi|1170286|sp|P44327|HIS7_HAEIN HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[0].alignments[9].length==362
assert record.rounds[0].alignments[10].title==">gi|2495224|sp|O06590|HIS7_MYCTU IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[10].length==210
assert record.rounds[0].alignments[11].title==">gi|123160|sp|P10368|HIS7_SALTY HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[0].alignments[11].length==354
assert record.rounds[0].alignments[12].title==">gi|2495226|sp|Q50504|HIS7_METTH PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[12].length==194
assert record.rounds[0].alignments[13].title==">gi|729718|sp|P40919|HIS7_CRYNE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[13].length==202
assert record.rounds[0].alignments[14].title==">gi|3334215|sp|O33773|HIS7_SULSO PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[14].length==193
assert record.rounds[0].alignments[15].title==">gi|123159|sp|P28624|HIS7_PHYPR IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[15].length==452
assert record.rounds[0].alignments[16].title==">gi|729719|sp|P40374|HIS7_SCHPO IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[16].length==216
assert record.rounds[0].alignments[17].title==">gi|2495227|sp|P56090|HIS7_CANAL IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[17].length==223
assert record.rounds[0].alignments[18].title==">gi|2495225|sp|Q58109|HIS7_METJA PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[18].length==197
assert record.rounds[0].alignments[19].title==">gi|2495229|sp|Q92447|HIS7_PICPA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[19].length==224
assert record.rounds[0].alignments[20].title==">gi|399897|sp|Q02986|HIS7_SACKL IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[20].length==232
assert record.rounds[0].alignments[21].title==">gi|2495228|sp|Q12578|HIS7_CANGA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[21].length==210
assert record.rounds[0].alignments[22].title==">gi|2506514|sp|P06633|HIS7_YEAST IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[22].length==220
assert record.rounds[0].alignments[23].title==">gi|462274|sp|P34041|HIS7_TRIHA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[0].alignments[23].length==208
assert record.rounds[0].alignments[24].title==">gi|1345641|sp|P49264|C7B1_THLAR CYTOCHROME P450 71B1 (CYPLXXIB1)"
assert record.rounds[0].alignments[24].length==496
assert record.rounds[0].alignments[25].title==">gi|1731346|sp|Q10698|YY29_MYCTU PROBABLE DIPEPTIDASE CY49.29C"
assert record.rounds[0].alignments[25].length==375
assert record.rounds[0].alignments[26].title==">gi|3287839|sp|Q01812|GLK4_RAT GLUTAMATE RECEPTOR, IONOTROPIC KAINATE 4 PRECURSOR (GLUTAMATE RECEPTOR KA-1) (KA1)"
assert record.rounds[0].alignments[26].length==956
assert record.rounds[0].alignments[27].title==">gi|3123025|sp|Q94637|VIT6_OSCBR VITELLOGENIN 6 PRECURSOR"
assert record.rounds[0].alignments[27].length==1660
assert record.rounds[0].alignments[28].title==">gi|1174406|sp|P36126|SP14_YEAST PHOSPHOLIPASE D1 (PLD 1) (CHOLINE PHOSPHATASE 1) (PHOSPHATIDYLCHOLINE-HYDROLYZING PHOSPHOLIPASE D1) (MEIOSIS-SPECIFIC SPORULATION PROTEIN SPO14)"
assert record.rounds[0].alignments[28].length==1380
assert record.rounds[0].alignments[29].title==">gi|3287848|sp|Q16099|GLK4_HUMAN GLUTAMATE RECEPTOR, IONOTROPIC KAINATE 4 PRECURSOR (GLUTAMATE RECEPTOR KA-1) (KA1) (EXCITATORY AMINO ACID RECEPTOR 1) (EAA1)"
assert record.rounds[0].alignments[29].length==956
assert len(record.rounds[1].new_seqs)==2
assert record.rounds[1].new_seqs[0].title=="gi|2833252|sp|Q14571|IP3S_HUMAN INOSITOL 1,4,5-TRISPHOSPHATE-BI..."
assert record.rounds[1].new_seqs[0].score==30
assert record.rounds[1].new_seqs[0].e==3.7
assert record.rounds[1].new_seqs[1].title=="gi|266389|sp|P29995|IP3S_RAT INOSITOL 1,4,5-TRISPHOSPHATE-BINDI..."
assert record.rounds[1].new_seqs[1].score==29
assert record.rounds[1].new_seqs[1].e==8.2
assert len(record.rounds[1].alignments)==26
assert record.rounds[1].alignments[0].title==">gi|2495230|sp|Q43072|HIS7_PEA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[0].length==281
assert record.rounds[1].alignments[1].title==">gi|462273|sp|P34047|HIS7_ARATH IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[1].length==270
assert record.rounds[1].alignments[2].title==">gi|399896|sp|Q02134|HIS7_LACLA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[2].length==200
assert record.rounds[1].alignments[3].title==">gi|1346293|sp|P48054|HIS7_SYNY3 IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[3].length==210
assert record.rounds[1].alignments[4].title==">gi|462272|sp|Q05068|HIS7_ANASP IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[4].length==209
assert record.rounds[1].alignments[5].title==">gi|462275|sp|P34048|HIS7_WHEAT IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[5].length==195
assert record.rounds[1].alignments[6].title==">gi|123161|sp|P16247|HIS7_STRCO IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[6].length==197
assert record.rounds[1].alignments[7].title==">gi|2506514|sp|P06633|HIS7_YEAST IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[7].length==220
assert record.rounds[1].alignments[8].title==">gi|2495227|sp|P56090|HIS7_CANAL IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[8].length==223
assert record.rounds[1].alignments[9].title==">gi|399897|sp|Q02986|HIS7_SACKL IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[9].length==232
assert record.rounds[1].alignments[10].title==">gi|2495228|sp|Q12578|HIS7_CANGA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[10].length==210
assert record.rounds[1].alignments[11].title==">gi|123158|sp|P06987|HIS7_ECOLI HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[1].alignments[11].length==355
assert record.rounds[1].alignments[12].title==">gi|123157|sp|P18787|HIS7_AZOBR IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[12].length==207
assert record.rounds[1].alignments[13].title==">gi|729718|sp|P40919|HIS7_CRYNE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[13].length==202
assert record.rounds[1].alignments[14].title==">gi|2495229|sp|Q92447|HIS7_PICPA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[14].length==224
assert record.rounds[1].alignments[15].title==">gi|1170286|sp|P44327|HIS7_HAEIN HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[1].alignments[15].length==362
assert record.rounds[1].alignments[16].title==">gi|123159|sp|P28624|HIS7_PHYPR IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[16].length==452
assert record.rounds[1].alignments[17].title==">gi|729719|sp|P40374|HIS7_SCHPO IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[17].length==216
assert record.rounds[1].alignments[18].title==">gi|123160|sp|P10368|HIS7_SALTY HISTIDINE BIOSYNTHESIS BIFUNCTIONAL PROTEIN HISB [INCLUDES: IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD); HISTIDINOL-PHOSPHATASE ]"
assert record.rounds[1].alignments[18].length==354
assert record.rounds[1].alignments[19].title==">gi|2495224|sp|O06590|HIS7_MYCTU IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[19].length==210
assert record.rounds[1].alignments[20].title==">gi|2495226|sp|Q50504|HIS7_METTH PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[20].length==194
assert record.rounds[1].alignments[21].title==">gi|2495225|sp|Q58109|HIS7_METJA PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[21].length==197
assert record.rounds[1].alignments[22].title==">gi|462274|sp|P34041|HIS7_TRIHA IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[22].length==208
assert record.rounds[1].alignments[23].title==">gi|3334215|sp|O33773|HIS7_SULSO PROBABLE IMIDAZOLEGLYCEROL-PHOSPHATE DEHYDRATASE (IGPD)"
assert record.rounds[1].alignments[23].length==193
assert record.rounds[1].alignments[24].title==">gi|2833252|sp|Q14571|IP3S_HUMAN INOSITOL 1,4,5-TRISPHOSPHATE-BINDING PROTEIN TYPE 2 RECEPTOR (TYPE 2 INSP3 RECEPTOR) (TYPE 2 INOSITOL 1,4,5-TRISPHOSPHATE RECEPTOR)"
assert record.rounds[1].alignments[24].length==2701
assert record.rounds[1].alignments[25].title==">gi|266389|sp|P29995|IP3S_RAT INOSITOL 1,4,5-TRISPHOSPHATE-BINDING PROTEIN TYPE 2 RECEPTOR (TYPE 2 INSP3 RECEPTOR) (TYPE 2 INOSITOL 1,4,5-TRISPHOSPHATE RECEPTOR)"
assert record.rounds[1].alignments[25].length==2701
assert record.rounds[0].alignments[0].hsps[0].score==1040
assert record.rounds[0].alignments[0].hsps[0].bits==409
assert record.rounds[0].alignments[0].hsps[0].expect==1e-114
assert len(record.rounds[0].alignments[0].hsps)==1
assert record.rounds[0].alignments[1].hsps[0].score==499
assert record.rounds[0].alignments[1].hsps[0].bits==198
assert record.rounds[0].alignments[1].hsps[0].expect==6e-51
assert len(record.rounds[0].alignments[1].hsps)==1
assert record.rounds[0].alignments[2].hsps[0].score==492
assert record.rounds[0].alignments[2].hsps[0].bits==196
assert record.rounds[0].alignments[2].hsps[0].expect==4e-50
assert len(record.rounds[0].alignments[2].hsps)==1
assert record.rounds[0].alignments[3].hsps[0].score==465
assert record.rounds[0].alignments[3].hsps[0].bits==185
assert record.rounds[0].alignments[3].hsps[0].expect==5e-47
assert len(record.rounds[0].alignments[3].hsps)==1
assert record.rounds[0].alignments[4].hsps[0].score==455
assert record.rounds[0].alignments[4].hsps[0].bits==181
assert record.rounds[0].alignments[4].hsps[0].expect==8e-46
assert len(record.rounds[0].alignments[4].hsps)==1
assert record.rounds[0].alignments[5].hsps[0].score==447
assert record.rounds[0].alignments[5].hsps[0].bits==178
assert record.rounds[0].alignments[5].hsps[0].expect==7e-45
assert len(record.rounds[0].alignments[5].hsps)==1
assert record.rounds[0].alignments[6].hsps[0].score==447
assert record.rounds[0].alignments[6].hsps[0].bits==178
assert record.rounds[0].alignments[6].hsps[0].expect==7e-45
assert len(record.rounds[0].alignments[6].hsps)==1
assert record.rounds[0].alignments[7].hsps[0].score==438
assert record.rounds[0].alignments[7].hsps[0].bits==175
assert record.rounds[0].alignments[7].hsps[0].expect==8e-44
assert len(record.rounds[0].alignments[7].hsps)==1
assert record.rounds[0].alignments[8].hsps[0].score==437
assert record.rounds[0].alignments[8].hsps[0].bits==174
assert record.rounds[0].alignments[8].hsps[0].expect==1e-43
assert len(record.rounds[0].alignments[8].hsps)==1
assert record.rounds[0].alignments[9].hsps[0].score==421
assert record.rounds[0].alignments[9].hsps[0].bits==168
assert record.rounds[0].alignments[9].hsps[0].expect==8e-42
assert len(record.rounds[0].alignments[9].hsps)==1
assert record.rounds[0].alignments[10].hsps[0].score==418
assert record.rounds[0].alignments[10].hsps[0].bits==167
assert record.rounds[0].alignments[10].hsps[0].expect==2e-41
assert len(record.rounds[0].alignments[10].hsps)==1
assert record.rounds[0].alignments[11].hsps[0].score==417
assert record.rounds[0].alignments[11].hsps[0].bits==166
assert record.rounds[0].alignments[11].hsps[0].expect==2e-41
assert len(record.rounds[0].alignments[11].hsps)==1
assert record.rounds[0].alignments[12].hsps[0].score==382
assert record.rounds[0].alignments[12].hsps[0].bits==153
assert record.rounds[0].alignments[12].hsps[0].expect==3e-37
assert len(record.rounds[0].alignments[12].hsps)==1
assert record.rounds[0].alignments[13].hsps[0].score==379
assert record.rounds[0].alignments[13].hsps[0].bits==152
assert record.rounds[0].alignments[13].hsps[0].expect==7e-37
assert len(record.rounds[0].alignments[13].hsps)==1
assert record.rounds[0].alignments[14].hsps[0].score==378
assert record.rounds[0].alignments[14].hsps[0].bits==151
assert record.rounds[0].alignments[14].hsps[0].expect==9e-37
assert len(record.rounds[0].alignments[14].hsps)==1
assert record.rounds[0].alignments[15].hsps[0].score==373
assert record.rounds[0].alignments[15].hsps[0].bits==149
assert record.rounds[0].alignments[15].hsps[0].expect==3e-36
assert len(record.rounds[0].alignments[15].hsps)==1
assert record.rounds[0].alignments[16].hsps[0].score==339
assert record.rounds[0].alignments[16].hsps[0].bits==136
assert record.rounds[0].alignments[16].hsps[0].expect==3e-32
assert len(record.rounds[0].alignments[16].hsps)==1
assert record.rounds[0].alignments[17].hsps[0].score==318
assert record.rounds[0].alignments[17].hsps[0].bits==128
assert record.rounds[0].alignments[17].hsps[0].expect==9e-30
assert len(record.rounds[0].alignments[17].hsps)==1
assert record.rounds[0].alignments[18].hsps[0].score==313
assert record.rounds[0].alignments[18].hsps[0].bits==126
assert record.rounds[0].alignments[18].hsps[0].expect==4e-29
assert len(record.rounds[0].alignments[18].hsps)==1
assert record.rounds[0].alignments[19].hsps[0].score==311
assert record.rounds[0].alignments[19].hsps[0].bits==125
assert record.rounds[0].alignments[19].hsps[0].expect==6e-29
assert len(record.rounds[0].alignments[19].hsps)==1
assert record.rounds[0].alignments[20].hsps[0].score==311
assert record.rounds[0].alignments[20].hsps[0].bits==125
assert record.rounds[0].alignments[20].hsps[0].expect==6e-29
assert len(record.rounds[0].alignments[20].hsps)==1
assert record.rounds[0].alignments[21].hsps[0].score==306
assert record.rounds[0].alignments[21].hsps[0].bits==123
assert record.rounds[0].alignments[21].hsps[0].expect==2e-28
assert len(record.rounds[0].alignments[21].hsps)==1
assert record.rounds[0].alignments[22].hsps[0].score==304
assert record.rounds[0].alignments[22].hsps[0].bits==122
assert record.rounds[0].alignments[22].hsps[0].expect==4e-28
assert len(record.rounds[0].alignments[22].hsps)==1
assert record.rounds[0].alignments[23].hsps[0].score==263
assert record.rounds[0].alignments[23].hsps[0].bits==106
assert record.rounds[0].alignments[23].hsps[0].expect==3e-23
assert len(record.rounds[0].alignments[23].hsps)==1
assert record.rounds[0].alignments[24].hsps[0].score==78
assert record.rounds[0].alignments[24].hsps[0].bits==34.8
assert record.rounds[0].alignments[24].hsps[0].expect==0.13
assert len(record.rounds[0].alignments[24].hsps)==1
assert record.rounds[0].alignments[25].hsps[0].score==70
assert record.rounds[0].alignments[25].hsps[0].bits==31.7
assert record.rounds[0].alignments[25].hsps[0].expect==1.1
assert len(record.rounds[0].alignments[25].hsps)==1
assert record.rounds[0].alignments[26].hsps[0].score==66
assert record.rounds[0].alignments[26].hsps[0].bits==30.1
assert record.rounds[0].alignments[26].hsps[0].expect==3.3
assert len(record.rounds[0].alignments[26].hsps)==1
assert record.rounds[0].alignments[27].hsps[0].score==64
assert record.rounds[0].alignments[27].hsps[0].bits==29.3
assert record.rounds[0].alignments[27].hsps[0].expect==5.6
assert len(record.rounds[0].alignments[27].hsps)==1
assert record.rounds[0].alignments[28].hsps[0].score==62
assert record.rounds[0].alignments[28].hsps[0].bits==28.6
assert record.rounds[0].alignments[28].hsps[0].expect==9.7
assert len(record.rounds[0].alignments[28].hsps)==1
assert record.rounds[0].alignments[29].hsps[0].score==62
assert record.rounds[0].alignments[29].hsps[0].bits==28.6
assert record.rounds[0].alignments[29].hsps[0].expect==9.7
assert record.rounds[1].alignments[0].hsps[0].score==820
assert record.rounds[1].alignments[0].hsps[0].bits==323
assert record.rounds[1].alignments[0].hsps[0].expect==1e-88
assert len(record.rounds[1].alignments[0].hsps)==1
assert record.rounds[1].alignments[1].hsps[0].score==817
assert record.rounds[1].alignments[1].hsps[0].bits==322
assert record.rounds[1].alignments[1].hsps[0].expect==3e-88
assert len(record.rounds[1].alignments[1].hsps)==1
assert record.rounds[1].alignments[2].hsps[0].score==808
assert record.rounds[1].alignments[2].hsps[0].bits==318
assert record.rounds[1].alignments[2].hsps[0].expect==4e-87
assert len(record.rounds[1].alignments[2].hsps)==1
assert record.rounds[1].alignments[3].hsps[0].score==798
assert record.rounds[1].alignments[3].hsps[0].bits==315
assert record.rounds[1].alignments[3].hsps[0].expect==5e-86
assert len(record.rounds[1].alignments[3].hsps)==1
assert record.rounds[1].alignments[4].hsps[0].score==795
assert record.rounds[1].alignments[4].hsps[0].bits==313
assert record.rounds[1].alignments[4].hsps[0].expect==1e-85
assert len(record.rounds[1].alignments[4].hsps)==1
assert record.rounds[1].alignments[5].hsps[0].score==793
assert record.rounds[1].alignments[5].hsps[0].bits==313
assert record.rounds[1].alignments[5].hsps[0].expect==2e-85
assert len(record.rounds[1].alignments[5].hsps)==1
assert record.rounds[1].alignments[6].hsps[0].score==776
assert record.rounds[1].alignments[6].hsps[0].bits==306
assert record.rounds[1].alignments[6].hsps[0].expect==2e-83
assert len(record.rounds[1].alignments[6].hsps)==1
assert record.rounds[1].alignments[7].hsps[0].score==772
assert record.rounds[1].alignments[7].hsps[0].bits==304
assert record.rounds[1].alignments[7].hsps[0].expect==6e-83
assert len(record.rounds[1].alignments[7].hsps)==1
assert record.rounds[1].alignments[8].hsps[0].score==771
assert record.rounds[1].alignments[8].hsps[0].bits==304
assert record.rounds[1].alignments[8].hsps[0].expect==8e-83
assert len(record.rounds[1].alignments[8].hsps)==1
assert record.rounds[1].alignments[9].hsps[0].score==770
assert record.rounds[1].alignments[9].hsps[0].bits==304
assert record.rounds[1].alignments[9].hsps[0].expect==1e-82
assert len(record.rounds[1].alignments[9].hsps)==1
assert record.rounds[1].alignments[10].hsps[0].score==767
assert record.rounds[1].alignments[10].hsps[0].bits==303
assert record.rounds[1].alignments[10].hsps[0].expect==2e-82
assert len(record.rounds[1].alignments[10].hsps)==1
assert record.rounds[1].alignments[11].hsps[0].score==765
assert record.rounds[1].alignments[11].hsps[0].bits==302
assert record.rounds[1].alignments[11].hsps[0].expect==4e-82
assert len(record.rounds[1].alignments[11].hsps)==1
assert record.rounds[1].alignments[12].hsps[0].score==762
assert record.rounds[1].alignments[12].hsps[0].bits==301
assert record.rounds[1].alignments[12].hsps[0].expect==9e-82
assert len(record.rounds[1].alignments[12].hsps)==1
assert record.rounds[1].alignments[13].hsps[0].score==759
assert record.rounds[1].alignments[13].hsps[0].bits==299
assert record.rounds[1].alignments[13].hsps[0].expect==2e-81
assert len(record.rounds[1].alignments[13].hsps)==1
assert record.rounds[1].alignments[14].hsps[0].score==756
assert record.rounds[1].alignments[14].hsps[0].bits==298
assert record.rounds[1].alignments[14].hsps[0].expect==5e-81
assert len(record.rounds[1].alignments[14].hsps)==1
assert record.rounds[1].alignments[15].hsps[0].score==741
assert record.rounds[1].alignments[15].hsps[0].bits==292
assert record.rounds[1].alignments[15].hsps[0].expect==3e-79
assert len(record.rounds[1].alignments[15].hsps)==1
assert record.rounds[1].alignments[16].hsps[0].score==734
assert record.rounds[1].alignments[16].hsps[0].bits==290
assert record.rounds[1].alignments[16].hsps[0].expect==2e-78
assert len(record.rounds[1].alignments[16].hsps)==1
assert record.rounds[1].alignments[17].hsps[0].score==734
assert record.rounds[1].alignments[17].hsps[0].bits==290
assert record.rounds[1].alignments[17].hsps[0].expect==2e-78
assert len(record.rounds[1].alignments[17].hsps)==1
assert record.rounds[1].alignments[18].hsps[0].score==726
assert record.rounds[1].alignments[18].hsps[0].bits==287
assert record.rounds[1].alignments[18].hsps[0].expect==1e-77
assert len(record.rounds[1].alignments[18].hsps)==1
assert record.rounds[1].alignments[19].hsps[0].score==716
assert record.rounds[1].alignments[19].hsps[0].bits==283
assert record.rounds[1].alignments[19].hsps[0].expect==2e-76
assert len(record.rounds[1].alignments[19].hsps)==1
assert record.rounds[1].alignments[20].hsps[0].score==695
assert record.rounds[1].alignments[20].hsps[0].bits==274
assert record.rounds[1].alignments[20].hsps[0].expect==6e-74
assert len(record.rounds[1].alignments[20].hsps)==1
assert record.rounds[1].alignments[21].hsps[0].score==685
assert record.rounds[1].alignments[21].hsps[0].bits==271
assert record.rounds[1].alignments[21].hsps[0].expect==1e-72
assert len(record.rounds[1].alignments[21].hsps)==1
assert record.rounds[1].alignments[22].hsps[0].score==680
assert record.rounds[1].alignments[22].hsps[0].bits==269
assert record.rounds[1].alignments[22].hsps[0].expect==4e-72
assert len(record.rounds[1].alignments[22].hsps)==1
assert record.rounds[1].alignments[23].hsps[0].score==662
assert record.rounds[1].alignments[23].hsps[0].bits==262
assert record.rounds[1].alignments[23].hsps[0].expect==5e-70
assert len(record.rounds[1].alignments[23].hsps)==1
assert record.rounds[1].alignments[24].hsps[0].score==66
assert record.rounds[1].alignments[24].hsps[0].bits==30.0
assert record.rounds[1].alignments[24].hsps[0].expect==3.7
assert len(record.rounds[1].alignments[24].hsps)==1
assert record.rounds[1].alignments[25].hsps[0].score==63
assert record.rounds[1].alignments[25].hsps[0].bits==28.8
assert record.rounds[1].alignments[25].hsps[0].expect==8.2
assert len(record.rounds[1].alignments[25].hsps)==1
assert record.rounds[0].alignments[0].hsps[0].identities==(200, 200)
assert record.rounds[0].alignments[0].hsps[0].positives==(200, 200)
assert record.rounds[0].alignments[1].hsps[0].identities==(99, 198)
assert record.rounds[0].alignments[1].hsps[0].positives==(135, 198)
assert record.rounds[0].alignments[1].hsps[0].gaps==(4, 198)
assert record.rounds[0].alignments[2].hsps[0].identities==(96, 199)
assert record.rounds[0].alignments[2].hsps[0].positives==(136, 199)
assert record.rounds[0].alignments[2].hsps[0].gaps==(4, 199)
assert record.rounds[0].alignments[3].hsps[0].identities==(91, 194)
assert record.rounds[0].alignments[3].hsps[0].positives==(126, 194)
assert record.rounds[0].alignments[3].hsps[0].gaps==(4, 194)
assert record.rounds[0].alignments[4].hsps[0].identities==(93, 194)
assert record.rounds[0].alignments[4].hsps[0].positives==(128, 194)
assert record.rounds[0].alignments[4].hsps[0].gaps==(4, 194)
assert record.rounds[0].alignments[5].hsps[0].identities==(89, 200)
assert record.rounds[0].alignments[5].hsps[0].positives==(124, 200)
assert record.rounds[0].alignments[5].hsps[0].gaps==(3, 200)
assert record.rounds[0].alignments[6].hsps[0].identities==(91, 198)
assert record.rounds[0].alignments[6].hsps[0].positives==(131, 198)
assert record.rounds[0].alignments[6].hsps[0].gaps==(4, 198)
assert record.rounds[0].alignments[7].hsps[0].identities==(91, 198)
assert record.rounds[0].alignments[7].hsps[0].positives==(130, 198)
assert record.rounds[0].alignments[7].hsps[0].gaps==(9, 198)
assert record.rounds[0].alignments[8].hsps[0].identities==(88, 198)
assert record.rounds[0].alignments[8].hsps[0].positives==(129, 198)
assert record.rounds[0].alignments[8].hsps[0].gaps==(4, 198)
assert record.rounds[0].alignments[9].hsps[0].identities==(89, 198)
assert record.rounds[0].alignments[9].hsps[0].positives==(127, 198)
assert record.rounds[0].alignments[9].hsps[0].gaps==(9, 198)
assert record.rounds[0].alignments[10].hsps[0].identities==(92, 207)
assert record.rounds[0].alignments[10].hsps[0].positives==(125, 207)
assert record.rounds[0].alignments[10].hsps[0].gaps==(14, 207)
assert record.rounds[0].alignments[11].hsps[0].identities==(89, 198)
assert record.rounds[0].alignments[11].hsps[0].positives==(129, 198)
assert record.rounds[0].alignments[11].hsps[0].gaps==(10, 198)
assert record.rounds[0].alignments[12].hsps[0].identities==(81, 198)
assert record.rounds[0].alignments[12].hsps[0].positives==(122, 198)
assert record.rounds[0].alignments[12].hsps[0].gaps==(8, 198)
assert record.rounds[0].alignments[13].hsps[0].identities==(83, 203)
assert record.rounds[0].alignments[13].hsps[0].positives==(120, 203)
assert record.rounds[0].alignments[13].hsps[0].gaps==(11, 203)
assert record.rounds[0].alignments[14].hsps[0].identities==(88, 201)
assert record.rounds[0].alignments[14].hsps[0].positives==(128, 201)
assert record.rounds[0].alignments[14].hsps[0].gaps==(9, 201)
assert record.rounds[0].alignments[15].hsps[0].identities==(86, 198)
assert record.rounds[0].alignments[15].hsps[0].positives==(120, 198)
assert record.rounds[0].alignments[15].hsps[0].gaps==(6, 198)
assert record.rounds[0].alignments[16].hsps[0].identities==(84, 221)
assert record.rounds[0].alignments[16].hsps[0].positives==(114, 221)
assert record.rounds[0].alignments[16].hsps[0].gaps==(29, 221)
assert record.rounds[0].alignments[17].hsps[0].identities==(81, 227)
assert record.rounds[0].alignments[17].hsps[0].positives==(119, 227)
assert record.rounds[0].alignments[17].hsps[0].gaps==(33, 227)
assert record.rounds[0].alignments[18].hsps[0].identities==(80, 196)
assert record.rounds[0].alignments[18].hsps[0].positives==(107, 196)
assert record.rounds[0].alignments[18].hsps[0].gaps==(9, 196)
assert record.rounds[0].alignments[19].hsps[0].identities==(84, 223)
assert record.rounds[0].alignments[19].hsps[0].positives==(116, 223)
assert record.rounds[0].alignments[19].hsps[0].gaps==(31, 223)
assert record.rounds[0].alignments[20].hsps[0].identities==(81, 222)
assert record.rounds[0].alignments[20].hsps[0].positives==(119, 222)
assert record.rounds[0].alignments[20].hsps[0].gaps==(30, 222)
assert record.rounds[0].alignments[21].hsps[0].identities==(78, 215)
assert record.rounds[0].alignments[21].hsps[0].positives==(116, 215)
assert record.rounds[0].alignments[21].hsps[0].gaps==(24, 215)
assert record.rounds[0].alignments[22].hsps[0].identities==(79, 218)
assert record.rounds[0].alignments[22].hsps[0].positives==(114, 218)
assert record.rounds[0].alignments[22].hsps[0].gaps==(30, 218)
assert record.rounds[0].alignments[23].hsps[0].identities==(68, 202)
assert record.rounds[0].alignments[23].hsps[0].positives==(102, 202)
assert record.rounds[0].alignments[23].hsps[0].gaps==(28, 202)
assert record.rounds[0].alignments[24].hsps[0].identities==(34, 134)
assert record.rounds[0].alignments[24].hsps[0].positives==(60, 134)
assert record.rounds[0].alignments[24].hsps[0].gaps==(11, 134)
assert record.rounds[0].alignments[25].hsps[0].identities==(16, 45)
assert record.rounds[0].alignments[25].hsps[0].positives==(21, 45)
assert record.rounds[0].alignments[25].hsps[0].gaps==(3, 45)
assert record.rounds[0].alignments[26].hsps[0].identities==(17, 48)
assert record.rounds[0].alignments[26].hsps[0].positives==(24, 48)
assert record.rounds[0].alignments[26].hsps[0].gaps==(3, 48)
assert record.rounds[0].alignments[27].hsps[0].identities==(25, 70)
assert record.rounds[0].alignments[27].hsps[0].positives==(32, 70)
assert record.rounds[0].alignments[27].hsps[0].gaps==(5, 70)
assert record.rounds[0].alignments[28].hsps[0].identities==(20, 65)
assert record.rounds[0].alignments[28].hsps[0].positives==(31, 65)
assert record.rounds[0].alignments[28].hsps[0].gaps==(7, 65)
assert record.rounds[0].alignments[29].hsps[0].identities==(16, 48)
assert record.rounds[0].alignments[29].hsps[0].positives==(24, 48)
assert record.rounds[0].alignments[29].hsps[0].gaps==(3, 48)
assert record.rounds[1].alignments[0].hsps[0].identities==(96, 199)
assert record.rounds[1].alignments[0].hsps[0].positives==(136, 199)
assert record.rounds[1].alignments[0].hsps[0].gaps==(4, 199)
assert record.rounds[1].alignments[1].hsps[0].identities==(99, 198)
assert record.rounds[1].alignments[1].hsps[0].positives==(135, 198)
assert record.rounds[1].alignments[1].hsps[0].gaps==(4, 198)
assert record.rounds[1].alignments[2].hsps[0].identities==(200, 200)
assert record.rounds[1].alignments[2].hsps[0].positives==(200, 200)
assert record.rounds[1].alignments[3].hsps[0].identities==(88, 198)
assert record.rounds[1].alignments[3].hsps[0].positives==(129, 198)
assert record.rounds[1].alignments[3].hsps[0].gaps==(4, 198)
assert record.rounds[1].alignments[4].hsps[0].identities==(91, 198)
assert record.rounds[1].alignments[4].hsps[0].positives==(131, 198)
assert record.rounds[1].alignments[4].hsps[0].gaps==(4, 198)
assert record.rounds[1].alignments[5].hsps[0].identities==(93, 196)
assert record.rounds[1].alignments[5].hsps[0].positives==(128, 196)
assert record.rounds[1].alignments[5].hsps[0].gaps==(4, 196)
assert record.rounds[1].alignments[6].hsps[0].identities==(89, 200)
assert record.rounds[1].alignments[6].hsps[0].positives==(124, 200)
assert record.rounds[1].alignments[6].hsps[0].gaps==(3, 200)
assert record.rounds[1].alignments[7].hsps[0].identities==(78, 220)
assert record.rounds[1].alignments[7].hsps[0].positives==(115, 220)
assert record.rounds[1].alignments[7].hsps[0].gaps==(30, 220)
assert record.rounds[1].alignments[8].hsps[0].identities==(81, 227)
assert record.rounds[1].alignments[8].hsps[0].positives==(119, 227)
assert record.rounds[1].alignments[8].hsps[0].gaps==(33, 227)
assert record.rounds[1].alignments[9].hsps[0].identities==(79, 222)
assert record.rounds[1].alignments[9].hsps[0].positives==(118, 222)
assert record.rounds[1].alignments[9].hsps[0].gaps==(30, 222)
assert record.rounds[1].alignments[10].hsps[0].identities==(76, 214)
assert record.rounds[1].alignments[10].hsps[0].positives==(113, 214)
assert record.rounds[1].alignments[10].hsps[0].gaps==(24, 214)
assert record.rounds[1].alignments[11].hsps[0].identities==(91, 198)
assert record.rounds[1].alignments[11].hsps[0].positives==(130, 198)
assert record.rounds[1].alignments[11].hsps[0].gaps==(9, 198)
assert record.rounds[1].alignments[12].hsps[0].identities==(91, 196)
assert record.rounds[1].alignments[12].hsps[0].positives==(127, 196)
assert record.rounds[1].alignments[12].hsps[0].gaps==(4, 196)
assert record.rounds[1].alignments[13].hsps[0].identities==(83, 203)
assert record.rounds[1].alignments[13].hsps[0].positives==(120, 203)
assert record.rounds[1].alignments[13].hsps[0].gaps==(11, 203)
assert record.rounds[1].alignments[14].hsps[0].identities==(82, 223)
assert record.rounds[1].alignments[14].hsps[0].positives==(115, 223)
assert record.rounds[1].alignments[14].hsps[0].gaps==(31, 223)
assert record.rounds[1].alignments[15].hsps[0].identities==(89, 198)
assert record.rounds[1].alignments[15].hsps[0].positives==(127, 198)
assert record.rounds[1].alignments[15].hsps[0].gaps==(9, 198)
assert record.rounds[1].alignments[16].hsps[0].identities==(86, 199)
assert record.rounds[1].alignments[16].hsps[0].positives==(121, 199)
assert record.rounds[1].alignments[16].hsps[0].gaps==(8, 199)
assert record.rounds[1].alignments[17].hsps[0].identities==(83, 221)
assert record.rounds[1].alignments[17].hsps[0].positives==(114, 221)
assert record.rounds[1].alignments[17].hsps[0].gaps==(29, 221)
assert record.rounds[1].alignments[18].hsps[0].identities==(89, 198)
assert record.rounds[1].alignments[18].hsps[0].positives==(129, 198)
assert record.rounds[1].alignments[18].hsps[0].gaps==(10, 198)
assert record.rounds[1].alignments[19].hsps[0].identities==(92, 207)
assert record.rounds[1].alignments[19].hsps[0].positives==(124, 207)
assert record.rounds[1].alignments[19].hsps[0].gaps==(14, 207)
assert record.rounds[1].alignments[20].hsps[0].identities==(81, 198)
assert record.rounds[1].alignments[20].hsps[0].positives==(122, 198)
assert record.rounds[1].alignments[20].hsps[0].gaps==(8, 198)
assert record.rounds[1].alignments[21].hsps[0].identities==(79, 196)
assert record.rounds[1].alignments[21].hsps[0].positives==(106, 196)
assert record.rounds[1].alignments[21].hsps[0].gaps==(9, 196)
assert record.rounds[1].alignments[22].hsps[0].identities==(68, 202)
assert record.rounds[1].alignments[22].hsps[0].positives==(102, 202)
assert record.rounds[1].alignments[22].hsps[0].gaps==(28, 202)
assert record.rounds[1].alignments[23].hsps[0].identities==(83, 200)
assert record.rounds[1].alignments[23].hsps[0].positives==(124, 200)
assert record.rounds[1].alignments[23].hsps[0].gaps==(7, 200)
assert record.rounds[1].alignments[24].hsps[0].identities==(18, 120)
assert record.rounds[1].alignments[24].hsps[0].positives==(39, 120)
assert record.rounds[1].alignments[24].hsps[0].gaps==(19, 120)
assert record.rounds[1].alignments[25].hsps[0].identities==(18, 120)
assert record.rounds[1].alignments[25].hsps[0].positives==(37, 120)
assert record.rounds[1].alignments[25].hsps[0].gaps==(19, 120)
assert record.rounds[0].alignments[0].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[0].hsps[0].match=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[0].hsps[0].sbjct=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[0].hsps[0].query_start==1
assert record.rounds[0].alignments[0].hsps[0].query_end==200
assert record.rounds[0].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[0].hsps[0].sbjct_end==200
assert record.rounds[0].alignments[1].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[1].hsps[0].match=="RI  + R TKET + + INLDGTG AD S+GI FLDHML  L  H  FD+ +   GD   V +D HH  ED+A+A+G  + + LG + GI R+G FT P+DEAL+   LD+SGRPYL ++ ++   Q++G YDT++ E FF++L   +G+TLH+ +  G+N+HHIIE  FK+ ARAL+QA   D  + G IPSSKGVL"
assert record.rounds[0].alignments[1].hsps[0].sbjct=="RIGEVKRVTKETNVSVKINLDGTGVADSSSGIPFLDHMLDQLASHGLFDVHVRATGD---VHIDDHHTNEDIALAIGTALLKALGERKGINRFGDFTAPLDEALIHVSLDLSGRPYLGYNLEIP-TQRVGTYDTQLVEHFFQSLVNTSGMTLHIRQLAGENSHHIIEATFKAFARALRQATETDPRRGGTIPSSKGVL"
assert record.rounds[0].alignments[1].hsps[0].query_start==3
assert record.rounds[0].alignments[1].hsps[0].query_end==200
assert record.rounds[0].alignments[1].hsps[0].sbjct_start==74
assert record.rounds[0].alignments[1].hsps[0].sbjct_end==267
assert record.rounds[0].alignments[2].hsps[0].query=="TRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[2].hsps[0].match=="TR+  + R TKET + + INLDG+G AD STGI FLDHML  L  H  FD+ +   GD   V +D HH  EDVA+A+G  + + LG++ GI R+G F+ P+DEAL+   LD+SGRP+L ++ D+   Q++G YDT++ E F +++   +G+TLH+ +  G+N+HHIIE  FK+ ARAL+QA   D  + G +PSSKGVL"
assert record.rounds[0].alignments[2].hsps[0].sbjct=="TRVGEVKRVTKETNVSVKINLDGSGVADSSTGIPFLDHMLDQLASHGLFDVHVKATGD---VHIDDHHTNEDVALAIGTALLQALGDRKGINRFGDFSAPLDEALIHVSLDLSGRPHLSYNLDIP-TQRVGTYDTQVVEHFLQSIVNTSGMTLHIRQLAGRNSHHIIEATFKAFARALRQATEYDPRRRGSVPSSKGVL"
assert record.rounds[0].alignments[2].hsps[0].query_start==2
assert record.rounds[0].alignments[2].hsps[0].query_end==200
assert record.rounds[0].alignments[2].hsps[0].sbjct_start==84
assert record.rounds[0].alignments[2].hsps[0].sbjct_end==278
assert record.rounds[0].alignments[3].hsps[0].query=="ITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[3].hsps[0].match=="I RNT ET+I +++NLDGTG  D+ TG+GFLDHML  L+ HS  DL +   GD   V +D HH  E   IA+G+ +++ +G++ GI+RYG   +PMDE L    LD S RPYL++    S   K+G  DTE+  E+F+A A  AG+TLH+   YG+N HHI+E  +K+ ARAL+  + ID  K   +PS+KG L"
assert record.rounds[0].alignments[3].hsps[0].sbjct=="IERNTTETRIRVAVNLDGTGVYDVKTGVGFLDHMLEQLSRHSLMDLSVAAEGD---VHIDAHHTTEHSGIAIGQAVAKAVGDRKGIQRYGHAYVPMDETLTRVALDFSNRPYLIWKVSFS-RDKIGDMDTELFREWFQAFAMAAGVTLHVECLYGENNHHIVESCYKALARALRAGIEIDPRKRDAVPSTKGTL"
assert record.rounds[0].alignments[3].hsps[0].query_start==7
assert record.rounds[0].alignments[3].hsps[0].query_end==200
assert record.rounds[0].alignments[3].hsps[0].sbjct_start==14
assert record.rounds[0].alignments[3].hsps[0].sbjct_end==203
assert record.rounds[0].alignments[4].hsps[0].query=="ITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[4].hsps[0].match=="+ R TKET + + INLDGTG A+ STGI FLDHML  L  H  FD+ +   GD     +D HH  ED+A+A+G  + + LG++ GI R+G FT P+DEA V   LD+SGRP+L     +   +++G YDT++ E FF++L   +G+TLH+ +  G N+HHIIE  FK+ ARAL+QA   D  + G +PSSKGVL"
assert record.rounds[0].alignments[4].hsps[0].sbjct=="VKRVTKETNVHVKINLDGTGVANSSTGIPFLDHMLDQLASHGLFDVYVKATGDTH---IDDHHSNEDIALAIGTALLQALGDRKGINRFGHFTAPLDEAAVEVILDLSGRPHLSCGLSIP-TERVGTYDTQLVEHFFQSLVNTSGMTLHIRQLAGNNSHHIIEATFKAFARALRQATEYDLRRQGTMPSSKGVL"
assert record.rounds[0].alignments[4].hsps[0].query_start==7
assert record.rounds[0].alignments[4].hsps[0].query_end==200
assert record.rounds[0].alignments[4].hsps[0].sbjct_start==3
assert record.rounds[0].alignments[4].hsps[0].sbjct_end==192
assert record.rounds[0].alignments[5].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[5].hsps[0].match=="M+R+  + R TKET + + I+LDGTG+ DI+TG+GF DHML  L  H  FDL +   GD   + +D HH IED A+ALG    + LG+K+GI R+G+ T+P+DE+L    +D+SGRPYLV     +    +G YD  MT     +    A + LH++  YG+N HHI+E  FK+ ARAL+ A   D    G +PS+KG L"
assert record.rounds[0].alignments[5].hsps[0].sbjct=="MSRVGRVERTTKETSVLVEIDLDGTGKTDIATGVGFYDHMLDQLGRHGLFDLTVKTDGD---LHIDSHHTIEDTALALGAAFRQALGDKVGIYRFGNCTVPLDESLAQVTVDLSGRPYLVHTEPENMAPMIGEYDVTMTRHILESFVAQAQVALHVHVPYGRNAHHIVECQFKALARALRYASERDPRAAGILPSTKGAL"
assert record.rounds[0].alignments[5].hsps[0].query_start==1
assert record.rounds[0].alignments[5].hsps[0].query_end==200
assert record.rounds[0].alignments[5].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[5].hsps[0].sbjct_end==197
assert record.rounds[0].alignments[6].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[6].hsps[0].match=="RI+ + R T ET +++++NLDGTG    +TGI FLDHML  ++ H   DL +   GD E   +D HH  EDV I LG+ +++ LG++ GI R+G+F  P+DEALV   LD SGRP+L +   +   +++G YDT++  EFF AL  ++ +TLH+ +  G N+HHIIE  FK+ ARA + A+ +D  + G IPSSKGVL"
assert record.rounds[0].alignments[6].hsps[0].sbjct=="RIASVHRITGETNVQVTVNLDGTGICKAATGIPFLDHMLHQISSHGLIDLDVQAKGDWE---IDDHHTNEDVGITLGQALAKALGDRKGIVRFGNFLAPLDEALVQVALDFSGRPHLSYGLQIP-TERVGTYDTQLVREFFVALVNHSQMTLHIRQLDGINSHHIIEATFKAFARAARMAIEVDPRRAGTIPSSKGVL"
assert record.rounds[0].alignments[6].hsps[0].query_start==3
assert record.rounds[0].alignments[6].hsps[0].query_end==200
assert record.rounds[0].alignments[6].hsps[0].sbjct_start==16
assert record.rounds[0].alignments[6].hsps[0].sbjct_end==209
assert record.rounds[0].alignments[7].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[7].hsps[0].match=="R +H+ RNTKETQI++ + LD  G + I+TG+GF DHML  +  H  F ++I   GD   + +D HH +ED  +ALG+ +   LG+K GI R+G F +PMDE L  C LDISGRP+L + A+ +  Q++G   TEM E FFR+L++  G+TLHL    G+N HH +E +FK+  R L+QA+ ++      +PSSKGVL"
assert record.rounds[0].alignments[7].hsps[0].sbjct=="RYAHVVRNTKETQIDVQVWLDREGGSKINTGVGFFDHMLDQIATHGGFRMEINVKGD---LYIDDHHTVEDTGLALGEALKIALGDKRGICRFG-FVLPMDECLARCALDISGRPHLEYKAEFT-YQRVGDLSTEMIEHFFRSLSYTMGVTLHLKTK-GKNDHHRVESLFKAFGRTLRQAIRVEGD---TLPSSKGVL"
assert record.rounds[0].alignments[7].hsps[0].query_start==3
assert record.rounds[0].alignments[7].hsps[0].query_end==200
assert record.rounds[0].alignments[7].hsps[0].sbjct_start==167
assert record.rounds[0].alignments[7].hsps[0].sbjct_end==355
assert record.rounds[0].alignments[8].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[8].hsps[0].match=="R + + R TKET + +S+NL G+G   ++TG+ FLDHML  +  H   DL++   GD E   +D HH  EDV I LG+ ++E LG++ GI R+G F  P+DEALV   LD SGRP+L +   +   +++G YDT++  EFF A+  ++ +TLH+ +  G N+HHIIE  FK+ ARA++ A+ +D  +   IPSSKGVL"
assert record.rounds[0].alignments[8].hsps[0].sbjct=="RAAAVHRVTKETDVRVSLNLMGSGLCHVATGVPFLDHMLHQIASHGLIDLEVNATGDIE---IDDHHTNEDVGITLGQALAEALGDRRGINRFGHFIAPLDEALVQVTLDFSGRPHLSYGLQIP-TERVGTYDTQLVREFFVAVVNHSQMTLHIRQLDGINSHHIIEATFKAFARAMRMAIEVDPRRADTIPSSKGVL"
assert record.rounds[0].alignments[8].hsps[0].query_start==3
assert record.rounds[0].alignments[8].hsps[0].query_end==200
assert record.rounds[0].alignments[8].hsps[0].sbjct_start==17
assert record.rounds[0].alignments[8].hsps[0].sbjct_end==210
assert record.rounds[0].alignments[9].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[9].hsps[0].match=="R + + R TKET I++ + LD  G  +I TG+GF DHML  +  H  F + +   GD   + +D HH +ED A+ALG+ + + +G+K GI R+G F +PMDE    C LD+SGRP++ F+A      K+G + TE+TE FF++LAF+   TLHLN   G N HH IE +FK+  R L+QA+ I+ +   E+PSSKGVL"
assert record.rounds[0].alignments[9].hsps[0].sbjct=="RFAEVIRQTKETDIKVQVWLDEAGVNEIKTGVGFFDHMLDQIATHGGFRMNVQCKGD---LWIDEHHTVEDTALALGQALKQAVGDKRGIARFG-FVLPMDECKAECALDLSGRPWIKFNACFK-RDKVGDFSTELTEHFFQSLAFSMLATLHLNV-TGNNDHHKIESLFKAFGRTLRQAIRIEGN---EMPSSKGVL"
assert record.rounds[0].alignments[9].hsps[0].query_start==3
assert record.rounds[0].alignments[9].hsps[0].query_end==200
assert record.rounds[0].alignments[9].hsps[0].sbjct_start==174
assert record.rounds[0].alignments[9].hsps[0].sbjct_end==362
assert record.rounds[0].alignments[10].hsps[0].query=="TRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVF--------HADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[10].hsps[0].match=="+R + I R T+E+ I + ++LDGTGQ  + TG+ F DHMLT L  H+ FDL +   GD E   ++ HH IED AIALG  + + LG+K GIRR+G   IPMDE L    +D+SGRPY V         H  ++G+     Y T +    F +LA NA I LH+   YG++ HHI E  +K+ ARAL+QAV  D  +V  +PS+KG L"
assert record.rounds[0].alignments[10].hsps[0].sbjct=="SRRARIERRTRESDIVIELDLDGTGQVAVDTGVPFYDHMLTALGSHASFDLTVRATGDVE---IEAHHTIEDTAIALGTALGQALGDKRGIRRFGDAFIPMDETLAHAAVDLSGRPYCVHTGEPDHLQHTTIAGSSV--PYHTVINRHVFESLAANARIALHVRVLYGRDPHHITEAQYKAVARALRQAVEPD-PRVSGVPSTKGAL"
assert record.rounds[0].alignments[10].hsps[0].query_start==2
assert record.rounds[0].alignments[10].hsps[0].query_end==200
assert record.rounds[0].alignments[10].hsps[0].sbjct_start==10
assert record.rounds[0].alignments[10].hsps[0].sbjct_end==210
assert record.rounds[0].alignments[11].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[11].hsps[0].match=="R +H+ RNTKETQI++S+ LD  G + I+TG+GF DHML  +  H  F ++I   GD   + +D HH +ED  +AL + +   L +K GI R+G F +PMDE L  C LDISGRP+L + A+ +  Q++G   TEM E FFR+L++  G+TLHL    G+N HH +E +FK+  R ++QA+ ++      +PSSKGVL"
assert record.rounds[0].alignments[11].hsps[0].sbjct=="RYAHVVRNTKETQIDVSVWLDREGNSKINTGVGFFDHMLDQIATHGGFRMEITVKGD---LYIDDHHTVEDTGLALREALKLALRDKRGICRFG-FVLPMDECL-ACALDISGRPHLEYKAEFT-YQRVGNLSTEMIEHFFRSLSYTMGVTLHLKTK-GKNDHHRVESLFKAFGRTVRQAIRVEGD---TLPSSKGVL"
assert record.rounds[0].alignments[11].hsps[0].query_start==3
assert record.rounds[0].alignments[11].hsps[0].query_end==200
assert record.rounds[0].alignments[11].hsps[0].sbjct_start==167
assert record.rounds[0].alignments[11].hsps[0].sbjct_end==354
assert record.rounds[0].alignments[12].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[12].hsps[0].match=="R S  TR T ET +++ + +DG+G++ ++TG+GFLDHML  +  H   DL++   GD E   +D HH +EDVA+ LG+ + E LG+K GIRR     +PMD+AL T  LD+SGRPY V   +   +  +G   ++    F  +LA +A + +H +   G+N HH  E +FK+ A A++ AV ++    GEIPS+KG L"
assert record.rounds[0].alignments[12].hsps[0].sbjct=="RRSMKTRETLETHVKVDLEIDGSGKSSVNTGLGFLDHMLESVARHGLLDLEVEARGDLE---VDDHHTVEDVALTLGEALREALGDKSGIRRMAHAMVPMDDALATVALDLSGRPYTVLELEFD-DAVIGDVKSQNIGHFIESLAVSAAMNIHASVR-GRNDHHKAEALFKALALAIRDAVRVEH---GEIPSTKGKL"
assert record.rounds[0].alignments[12].hsps[0].query_start==3
assert record.rounds[0].alignments[12].hsps[0].query_end==200
assert record.rounds[0].alignments[12].hsps[0].sbjct_start==5
assert record.rounds[0].alignments[12].hsps[0].sbjct_end==194
assert record.rounds[0].alignments[13].hsps[0].query=="RISHITRNTKETQIELSINLDGTG-----QADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[13].hsps[0].match=="RI+ + R T ET I  +I+LD        + ++STGIGFLDHM T L  H    L++   GD   + +D HH  ED A+ALG+   + LG + GI+RYG    P+DE+L    +DIS RPY + H   +  +K+G   TEM     ++ AF AG+TLH++   G+N HHI E  FK+ A A++ A+S   +   ++PS+KGVL"
assert record.rounds[0].alignments[13].hsps[0].sbjct=="RIASVERTTSETHISCTIDLDHIPGVTEQKINVSTGIGFLDHMFTALAKHGGMSLQLQCKGD---LHIDDHHTAEDCALALGEAFKKALGERKGIKRYGYAYAPLDESLSRAVIDISSRPYFMCHLPFT-REKVGDLSTEMVSHLLQSFAFAAGVTLHIDSIRGENNHHIAESAFKALALAIRMAIS--RTGGDDVPSTKGVL"
assert record.rounds[0].alignments[13].hsps[0].query_start==3
assert record.rounds[0].alignments[13].hsps[0].query_end==200
assert record.rounds[0].alignments[13].hsps[0].sbjct_start==4
assert record.rounds[0].alignments[13].hsps[0].sbjct_end==200
assert record.rounds[0].alignments[14].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHML-TLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[14].hsps[0].match=="M+R ++ITR TKET+IE+ +++D  G+  +ST I F +HML TLLT+ +     I+   D   +  D HH++EDVAI LG  I   LG+K GI+R+    IPMD+ALV   LDIS R     + +L  ++ +GG  TE    FF++ A+N+GITLH+++  G NTHHIIE  FK+   AL +A  I ++   EI S+KG++"
assert record.rounds[0].alignments[14].hsps[0].sbjct=="MSRSANITRETKETKIEVLLDIDRKGEVKVSTPIPFFNHMLITLLTYMNS--TAIVSATDK--LPYDDHHIVEDVAITLGLAIKTALGDKRGIKRFSHQIIPMDDALVLVSLDISNRGMAFVNLNLKRSE-IGGLATENVPHFFQSFAYNSGITLHISQLSGYNTHHIIEASFKALGLALYEATRIVDN---EIRSTKGII"
assert record.rounds[0].alignments[14].hsps[0].query_start==1
assert record.rounds[0].alignments[14].hsps[0].query_end==200
assert record.rounds[0].alignments[14].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[14].hsps[0].sbjct_end==193
assert record.rounds[0].alignments[15].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[15].hsps[0].match=="R ++I+R TKET I + ++LDGTG++ +S+GIGFLDHMLT L  HS FDL++   GD     +D HH  ED A+ LG+     LG++ GI R+GS  +P+DEAL    +DIS R +   +  L     +G   +EM   FF + A  A  TLH++   G+N HH  E  FK+ A AL+ AV  D +    +PS+KGVL"
assert record.rounds[0].alignments[15].hsps[0].sbjct=="REANISRVTKETSISVKLSLDGTGKSKVSSGIGFLDHMLTALAKHSRFDLELDCKGD---TWIDDHHTTEDCALTLGEAFDVALGDRAGIARFGSACVPLDEALSRAIVDISSRAHSEINLQLV-RPSVGELSSEMITHFFESFASAALXTLHVDVLRGRNDHHRAEASFKALAVALRTAVKHDAT--AGVPSTKGVL"
assert record.rounds[0].alignments[15].hsps[0].query_start==3
assert record.rounds[0].alignments[15].hsps[0].query_end==200
assert record.rounds[0].alignments[15].hsps[0].sbjct_start==260
assert record.rounds[0].alignments[15].hsps[0].sbjct_end==451
assert record.rounds[0].alignments[16].hsps[0].query=="RISHITRNTKETQIELSINLD-----------------------GTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[16].hsps[0].match=="R + + RNT ET+I ++I LD                       G     + TGIGFLDHM   L  H+ + L++   GD   + +D HH  ED AIALG    + +GN  G++R+G    P+DEAL    +D+SGRPY V    L   +K+G    EM      + +  AGITLH+   YG N HH  E  FKS A A++ A S+  S   E+PS+KGVL"
assert record.rounds[0].alignments[16].hsps[0].sbjct=="RRAFVERNTNETKISVAIALDKAPLPEESNFIDELITSKHANQKGEQVIQVDTGIGFLDHMYHALAKHAGWSLRLYSRGD---LIIDDHHTAEDTAIALGIAFKQAMGNFAGVKRFGHAYCPLDEALSRSVVDLSGRPYAVIDLGLK-REKVGELSCEMIPHLLYSFSVAAGITLHVTCLYGSNDHHRAESAFKSLAVAMRAATSLTGS--SEVPSTKGVL"
assert record.rounds[0].alignments[16].hsps[0].query_start==3
assert record.rounds[0].alignments[16].hsps[0].query_end==200
assert record.rounds[0].alignments[16].hsps[0].sbjct_start==2
assert record.rounds[0].alignments[16].hsps[0].sbjct_end==216
assert record.rounds[0].alignments[17].hsps[0].query=="MTRISHITRNTKETQIELSINLDG---------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[17].hsps[0].match=="M+R + I R T ET+I++++NLDG                           +   ++ TGIGFLDHM+  L  HS + L +   GD   + +D HH  EDV I+LG    + LG   G++R+G    P+DEAL    +D+S RP+ V    L   +K+G   TEM      + A   GIT+H++   G N HH  E  FK+ A A+K+A+S  ++   +IPS+KGVL"
assert record.rounds[0].alignments[17].hsps[0].sbjct=="MSREALINRITNETKIQIALNLDGGKLELKESIFPNQSIIIDEHHAKQVSGSQYINVQTGIGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTAEDVGISLGMAFKQALGQIKGVKRFGHGFAPLDEALSRAVVDLSNRPFAVIELGLK-REKIGDLSTEMIPHVLESFAGAVGITIHVDCLRGFNDHHRAESAFKALAIAIKEAIS--KTGKNDIPSTKGVL"
assert record.rounds[0].alignments[17].hsps[0].query_start==1
assert record.rounds[0].alignments[17].hsps[0].query_end==200
assert record.rounds[0].alignments[17].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[17].hsps[0].sbjct_end==221
assert record.rounds[0].alignments[18].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKG"
assert record.rounds[0].alignments[18].hsps[0].match=="RI  + R TKET I L IN+DGTG+  I TGI F DH+L     H  FDL +   GD E   +D HH +EDV I LG  +++    K  I R+G   IPMD+A  T  +D+SGR Y V + +    + +G   TE    FF ++A    + +H  E  G+N HH  E +FK+   AL  A  IDE K   + S+KG"
assert record.rounds[0].alignments[18].hsps[0].sbjct=="RIFEVMRETKETNIYLKINIDGTGKYKIDTGIPFFDHLLASFAKHGCFDLIVKARGDLE---IDDHHTVEDVGICLGLALNQI--EKRNIFRFGWAIIPMDDARATVAIDLSGRSYCVGNYE-PKREFVGDLATENINHFFESVASYGMLNIHY-EVIGKNEHHKAEALFKAFGVALDLATKIDERK--GVISTKG"
assert record.rounds[0].alignments[18].hsps[0].query_start==3
assert record.rounds[0].alignments[18].hsps[0].query_end==198
assert record.rounds[0].alignments[18].hsps[0].sbjct_start==7
assert record.rounds[0].alignments[18].hsps[0].sbjct_end==193
assert record.rounds[0].alignments[19].hsps[0].query=="RISHITRNTKETQIELSINLDG-------------------TGQA------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[19].hsps[0].match=="R S I R T ET+I+++++LDG                     QA       ++TGIGFLDHML  L  H  + + I   GD   + +D HH  ED  IALG    E LG+  GI+R+GS   P+DEAL    +D+S RPY V    L   +K+G    EM      + A  A +T+H++   G N HH  E  FK+ A A+K+A+S   +   +IPS+KGVL"
assert record.rounds[0].alignments[19].hsps[0].sbjct=="RSSLIKRITNETKIQIALSLDGGPVSLAQSLFKDKDYSAEHAAQATSSQFISVNTGIGFLDHMLHALAKHGGWSVIIECVGD---LHIDDHHSAEDTGIALGMAFKEALGHVRGIKRFGSGFAPLDEALSRAVIDMSNRPYAVVDLGLK-REKIGDLSCEMIPHVLESFAQGAHVTMHVDCLRGFNDHHRAESAFKALAIAIKEAIS--SNGTDDIPSTKGVL"
assert record.rounds[0].alignments[19].hsps[0].query_start==3
assert record.rounds[0].alignments[19].hsps[0].query_end==200
assert record.rounds[0].alignments[19].hsps[0].sbjct_start==7
assert record.rounds[0].alignments[19].hsps[0].sbjct_end==223
assert record.rounds[0].alignments[20].hsps[0].query=="RISHITRNTKETQIELSINLDG------------------TGQA------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[20].hsps[0].match=="R + I+R T ET+I+++I+L+G                    QA      DI TG+GFLDHM+  L  HS + L +   GD   + +D HH  ED  IALG+   E +G   G++R+G+   P+DEAL    +D+S RP+ V    L   + +G   TEM   F  + A  A ITLH++   G N HH  E  FK+ A A+++A+S   +   ++PS+KGVL"
assert record.rounds[0].alignments[20].hsps[0].sbjct=="RKAFISRITNETKIQIAISLNGGYIQIKDSILPAKKDDDVASQATQSQVIDIHTGVGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTTEDCGIALGQAFKEAMGAVRGVKRFGTGFAPLDEALSRAVVDLSNRPFAVIDLGLK-REMIGDLSTEMIPHFLESFAEAARITLHVDCLRGFNDHHRSESAFKALAVAIREAIS--SNGTNDVPSTKGVL"
assert record.rounds[0].alignments[20].hsps[0].query_start==3
assert record.rounds[0].alignments[20].hsps[0].query_end==200
assert record.rounds[0].alignments[20].hsps[0].sbjct_start==16
assert record.rounds[0].alignments[20].hsps[0].sbjct_end==231
assert record.rounds[0].alignments[21].hsps[0].query=="ISHITRNTKETQIELSINLDG-----------------TGQA-DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[21].hsps[0].match=="++ + R T+ET I+L+++LDG                  GQ   + TG+GFLDHMLT L  H  + L +   GD   + +D HH +ED  IALG+   E LG+  GI+R+G    P+DEAL    +D S RP+ V    L   +++G   TEM   F  + A    IT+H++   G N HH  E  FK+ A A+++A +   +   ++PS+KGVL"
assert record.rounds[0].alignments[21].hsps[0].sbjct=="MAFVKRVTQETNIQLALDLDGGSVSVRESILGKEYASGDGQTIHVHTGVGFLDHMLTALAKHGGWSLILECIGD---LHIDDHHTVEDCGIALGQAFKEALGSVRGIKRFGHGFAPLDEALSRAVVDFSNRPFAVVELGLK-RERIGQLSTEMIPHFLESFATEGRITMHVDCLRGTNDHHRSESGFKALAIAIREART--PTGRDDVPSTKGVL"
assert record.rounds[0].alignments[21].hsps[0].query_start==4
assert record.rounds[0].alignments[21].hsps[0].query_end==200
assert record.rounds[0].alignments[21].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[21].hsps[0].sbjct_end==209
assert record.rounds[0].alignments[22].hsps[0].query=="ITRNTKETQIELSINLDGTGQA------------------------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[0].alignments[22].hsps[0].match=="+ R T ET+I+++I+L G   A                        ++ TGIGFLDHM+  L  HS + L +   GD   + +D HH  ED  IALG+   E LG   G++R+GS   P+DEAL    +D+S RPY V    L   +K+G    EM   F  + A  + ITLH++   G+N HH  E  FK+ A A+++A S   +   ++PS+KGVL"
assert record.rounds[0].alignments[22].hsps[0].sbjct=="VKRITNETKIQIAISLKGGPLAIEHSIFPEKEAEAVAEQATQSQVINVHTGIGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTTEDCGIALGQAFKEALGAVRGVKRFGSGFAPLDEALSRAVVDLSNRPYAVVELGLQ-REKVGDLSCEMIPHFLESFAEASRITLHVDCLRGKNDHHRSESAFKALAVAIREATS--PNGTNDVPSTKGVL"
assert record.rounds[0].alignments[22].hsps[0].query_start==7
assert record.rounds[0].alignments[22].hsps[0].query_end==200
assert record.rounds[0].alignments[22].hsps[0].sbjct_start==8
assert record.rounds[0].alignments[22].hsps[0].sbjct_end==219
assert record.rounds[0].alignments[23].hsps[0].query=="RISHITRNTKETQIELSINLDG------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALK"
assert record.rounds[0].alignments[23].hsps[0].match=="R + ++R+T ET I++++++DG                        +    I+TGIGFLDHML  L  H+ + + +   GD   + +D HH  ED  IA+G   ++ LG   G+ R+G    P+DEAL    +D+S RPY V    L   +KLG    EM     ++ A  A ITLH++   G N HH  E  FK+ A A++"
assert record.rounds[0].alignments[23].hsps[0].sbjct=="RAAALSRDTNETSIQIALSIDGGELPQDADPRLLEASSAHASQTSKSQVISINTGIGFLDHMLHALAKHAGWSMALNCKGD---LHIDDHHTAEDCCIAVGTTFAKALGALTGVARFGYAYAPLDEALSRAVVDLSNRPYTVVDLGLK-REKLGELSCEMIPHCLQSFAQAARITLHVDCLRGDNDHHRAESAFKALAVAVR"
assert record.rounds[0].alignments[23].hsps[0].query_start==3
assert record.rounds[0].alignments[23].hsps[0].query_end==180
assert record.rounds[0].alignments[23].hsps[0].sbjct_start==8
assert record.rounds[0].alignments[23].hsps[0].sbjct_end==205
assert record.rounds[0].alignments[24].hsps[0].query=="GDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQK---LGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDES"
assert record.rounds[0].alignments[24].hsps[0].match=="G+   +G  PH  + +++   G  +S  LG+   +      T+   + L T DL+   RPY+ + A ++ N K      YD     +++R +     + L+  +   Q+  HI E    S  R  KQA S +E+"
assert record.rounds[0].alignments[24].hsps[0].sbjct=="GNLHQLGEKPHRAMVELSKTYGPLMSLKLGSVTTVVATSVETVR--DVLKTYDLECCSRPYMTYPARITYNLKDLVFSPYD-----KYWRQVRKLTVVELYTAKRV-QSFRHIREEEVASFVRFNKQAASSEET"
assert record.rounds[0].alignments[24].hsps[0].query_start==58
assert record.rounds[0].alignments[24].hsps[0].query_end==188
assert record.rounds[0].alignments[24].hsps[0].sbjct_start==40
assert record.rounds[0].alignments[24].hsps[0].sbjct_end==165
assert record.rounds[0].alignments[25].hsps[0].query=="HSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLG"
assert record.rounds[0].alignments[25].hsps[0].match=="HS+    I+G G H   G DPHH   D  +  G  +  D+G   G"
assert record.rounds[0].alignments[25].hsps[0].sbjct=="HSEVAFVIVGSGPH---GADPHHGYSDRELREGDIVVVDIGGTYG"
assert record.rounds[0].alignments[25].hsps[0].query_start==47
assert record.rounds[0].alignments[25].hsps[0].query_end==91
assert record.rounds[0].alignments[25].hsps[0].sbjct_start==195
assert record.rounds[0].alignments[25].hsps[0].sbjct_end==236
assert record.rounds[0].alignments[26].hsps[0].query=="PYLVF---HADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYG"
assert record.rounds[0].alignments[26].hsps[0].match=="PYL+    H D+ GN +  G+  +M +E    L FN  I L  +  YG"
assert record.rounds[0].alignments[26].hsps[0].sbjct=="PYLMLKGNHQDMEGNDRYEGFCVDMLKELAEILRFNYKIRLVGDGVYG"
assert record.rounds[0].alignments[26].hsps[0].query_start==117
assert record.rounds[0].alignments[26].hsps[0].query_end==161
assert record.rounds[0].alignments[26].hsps[0].sbjct_start==427
assert record.rounds[0].alignments[26].hsps[0].sbjct_end==474
assert record.rounds[0].alignments[27].hsps[0].query=="GKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLG---GYDTEMTEEFFRA"
assert record.rounds[0].alignments[27].hsps[0].match=="GK ISED     GI + G     +D   VT + D  G    +  + L  N++ G    YD E T EF RA"
assert record.rounds[0].alignments[27].hsps[0].sbjct=="GKEISEDKWEDFGISQRGEEKFFIDAEKVTVEFD--GFQAKIQMSSLYKNKQCGLCGHYDGEKTNEFRRA"
assert record.rounds[0].alignments[27].hsps[0].query_start==79
assert record.rounds[0].alignments[27].hsps[0].query_end==145
assert record.rounds[0].alignments[27].hsps[0].sbjct_start==1436
assert record.rounds[0].alignments[27].hsps[0].sbjct_end==1503
assert record.rounds[0].alignments[28].hsps[0].query=="RYGSFTIPMDEAL-----VTCDLDISGRPYLVFHADLSGNQKL--GGYDTEMTEEFFRALAFNAG"
assert record.rounds[0].alignments[28].hsps[0].match=="R GS ++P++E          D DI G  Y   H D+  NQ+L  G  D  +++   +   FN+G"
assert record.rounds[0].alignments[28].hsps[0].sbjct=="RNGSNSLPLNEKSNEGESTNVDQDIEGDEYHRLHEDILKNQELDDGSLDDLLSQIIPKITNFNSG"
assert record.rounds[0].alignments[28].hsps[0].query_start==94
assert record.rounds[0].alignments[28].hsps[0].query_end==151
assert record.rounds[0].alignments[28].hsps[0].sbjct_start==1141
assert record.rounds[0].alignments[28].hsps[0].sbjct_end==1205
assert record.rounds[0].alignments[29].hsps[0].query=="PYLVF---HADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYG"
assert record.rounds[0].alignments[29].hsps[0].match=="PYL+    H ++ GN +  G+  +M +E    L FN  I L  +  YG"
assert record.rounds[0].alignments[29].hsps[0].sbjct=="PYLMLKGNHQEMEGNDRYEGFCVDMLKELAEILRFNYKIRLVGDGVYG"
assert record.rounds[0].alignments[29].hsps[0].query_start==117
assert record.rounds[0].alignments[29].hsps[0].query_end==161
assert record.rounds[0].alignments[29].hsps[0].sbjct_start==427
assert record.rounds[0].alignments[29].hsps[0].sbjct_end==474
assert record.rounds[1].alignments[0].hsps[0].query=="TRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[0].hsps[0].match=="TR+  + R TKET + + INLDG+G AD STGI FLDHML  L  H  FD+ +   GD   V +D HH  EDVA+A+G  + + LG++ GI R+G F+ P+DEAL+   LD+SGRP+L ++ D+   Q++G YDT++ E F +++   +G+TLH+ +  G+N+HHIIE  FK+ ARAL+QA   D  + G +PSSKGVL"
assert record.rounds[1].alignments[0].hsps[0].sbjct=="TRVGEVKRVTKETNVSVKINLDGSGVADSSTGIPFLDHMLDQLASHGLFDVHVKATGD---VHIDDHHTNEDVALAIGTALLQALGDRKGINRFGDFSAPLDEALIHVSLDLSGRPHLSYNLDIP-TQRVGTYDTQVVEHFLQSIVNTSGMTLHIRQLAGRNSHHIIEATFKAFARALRQATEYDPRRRGSVPSSKGVL"
assert record.rounds[1].alignments[0].hsps[0].query_start==2
assert record.rounds[1].alignments[0].hsps[0].query_end==200
assert record.rounds[1].alignments[0].hsps[0].sbjct_start==84
assert record.rounds[1].alignments[0].hsps[0].sbjct_end==278
assert record.rounds[1].alignments[1].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[1].hsps[0].match=="RI  + R TKET + + INLDGTG AD S+GI FLDHML  L  H  FD+ +   GD   V +D HH  ED+A+A+G  + + LG + GI R+G FT P+DEAL+   LD+SGRPYL ++ ++   Q++G YDT++ E FF++L   +G+TLH+ +  G+N+HHIIE  FK+ ARAL+QA   D  + G IPSSKGVL"
assert record.rounds[1].alignments[1].hsps[0].sbjct=="RIGEVKRVTKETNVSVKINLDGTGVADSSSGIPFLDHMLDQLASHGLFDVHVRATGD---VHIDDHHTNEDIALAIGTALLKALGERKGINRFGDFTAPLDEALIHVSLDLSGRPYLGYNLEIP-TQRVGTYDTQLVEHFFQSLVNTSGMTLHIRQLAGENSHHIIEATFKAFARALRQATETDPRRGGTIPSSKGVL"
assert record.rounds[1].alignments[1].hsps[0].query_start==3
assert record.rounds[1].alignments[1].hsps[0].query_end==200
assert record.rounds[1].alignments[1].hsps[0].sbjct_start==74
assert record.rounds[1].alignments[1].hsps[0].sbjct_end==267
assert record.rounds[1].alignments[2].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[2].hsps[0].match=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[2].hsps[0].sbjct=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[2].hsps[0].query_start==1
assert record.rounds[1].alignments[2].hsps[0].query_end==200
assert record.rounds[1].alignments[2].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[2].hsps[0].sbjct_end==200
assert record.rounds[1].alignments[3].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[3].hsps[0].match=="R + + R TKET + +S+NL G+G   ++TG+ FLDHML  +  H   DL++   GD E   +D HH  EDV I LG+ ++E LG++ GI R+G F  P+DEALV   LD SGRP+L +   +   +++G YDT++  EFF A+  ++ +TLH+ +  G N+HHIIE  FK+ ARA++ A+ +D  +   IPSSKGVL"
assert record.rounds[1].alignments[3].hsps[0].sbjct=="RAAAVHRVTKETDVRVSLNLMGSGLCHVATGVPFLDHMLHQIASHGLIDLEVNATGDIE---IDDHHTNEDVGITLGQALAEALGDRRGINRFGHFIAPLDEALVQVTLDFSGRPHLSYGLQIP-TERVGTYDTQLVREFFVAVVNHSQMTLHIRQLDGINSHHIIEATFKAFARAMRMAIEVDPRRADTIPSSKGVL"
assert record.rounds[1].alignments[3].hsps[0].query_start==3
assert record.rounds[1].alignments[3].hsps[0].query_end==200
assert record.rounds[1].alignments[3].hsps[0].sbjct_start==17
assert record.rounds[1].alignments[3].hsps[0].sbjct_end==210
assert record.rounds[1].alignments[4].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[4].hsps[0].match=="RI+ + R T ET +++++NLDGTG    +TGI FLDHML  ++ H   DL +   GD E   +D HH  EDV I LG+ +++ LG++ GI R+G+F  P+DEALV   LD SGRP+L +   +   +++G YDT++  EFF AL  ++ +TLH+ +  G N+HHIIE  FK+ ARA + A+ +D  + G IPSSKGVL"
assert record.rounds[1].alignments[4].hsps[0].sbjct=="RIASVHRITGETNVQVTVNLDGTGICKAATGIPFLDHMLHQISSHGLIDLDVQAKGDWE---IDDHHTNEDVGITLGQALAKALGDRKGIVRFGNFLAPLDEALVQVALDFSGRPHLSYGLQIP-TERVGTYDTQLVREFFVALVNHSQMTLHIRQLDGINSHHIIEATFKAFARAARMAIEVDPRRAGTIPSSKGVL"
assert record.rounds[1].alignments[4].hsps[0].query_start==3
assert record.rounds[1].alignments[4].hsps[0].query_end==200
assert record.rounds[1].alignments[4].hsps[0].sbjct_start==16
assert record.rounds[1].alignments[4].hsps[0].sbjct_end==209
assert record.rounds[1].alignments[5].hsps[0].query=="SHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[5].hsps[0].match=="  + R TKET + + INLDGTG A+ STGI FLDHML  L  H  FD+ +   GD     +D HH  ED+A+A+G  + + LG++ GI R+G FT P+DEA V   LD+SGRP+L     +   +++G YDT++ E FF++L   +G+TLH+ +  G N+HHIIE  FK+ ARAL+QA   D  + G +PSSKGVL"
assert record.rounds[1].alignments[5].hsps[0].sbjct=="GEVKRVTKETNVHVKINLDGTGVANSSTGIPFLDHMLDQLASHGLFDVYVKATGDT---HIDDHHSNEDIALAIGTALLQALGDRKGINRFGHFTAPLDEAAVEVILDLSGRPHLSCGLSIP-TERVGTYDTQLVEHFFQSLVNTSGMTLHIRQLAGNNSHHIIEATFKAFARALRQATEYDLRRQGTMPSSKGVL"
assert record.rounds[1].alignments[5].hsps[0].query_start==5
assert record.rounds[1].alignments[5].hsps[0].query_end==200
assert record.rounds[1].alignments[5].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[5].hsps[0].sbjct_end==192
assert record.rounds[1].alignments[6].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[6].hsps[0].match=="M+R+  + R TKET + + I+LDGTG+ DI+TG+GF DHML  L  H  FDL +   GD   + +D HH IED A+ALG    + LG+K+GI R+G+ T+P+DE+L    +D+SGRPYLV     +    +G YD  MT     +    A + LH++  YG+N HHI+E  FK+ ARAL+ A   D    G +PS+KG L"
assert record.rounds[1].alignments[6].hsps[0].sbjct=="MSRVGRVERTTKETSVLVEIDLDGTGKTDIATGVGFYDHMLDQLGRHGLFDLTVKTDGD---LHIDSHHTIEDTALALGAAFRQALGDKVGIYRFGNCTVPLDESLAQVTVDLSGRPYLVHTEPENMAPMIGEYDVTMTRHILESFVAQAQVALHVHVPYGRNAHHIVECQFKALARALRYASERDPRAAGILPSTKGAL"
assert record.rounds[1].alignments[6].hsps[0].query_start==1
assert record.rounds[1].alignments[6].hsps[0].query_end==200
assert record.rounds[1].alignments[6].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[6].hsps[0].sbjct_end==197
assert record.rounds[1].alignments[7].hsps[0].query=="SHITRNTKETQIELSINLDG------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[7].hsps[0].match=="+ + R T ET+I+++I+L G                        +   ++ TGIGFLDHM+  L  HS + L +   GD   + +D HH  ED  IALG+   E LG   G++R+GS   P+DEAL    +D+S RPY V    L   +K+G    EM   F  + A  + ITLH++   G+N HH  E  FK+ A A+++A S   +   ++PS+KGVL"
assert record.rounds[1].alignments[7].hsps[0].sbjct=="ALVKRITNETKIQIAISLKGGPLAIEHSIFPEKEAEAVAEQATQSQVINVHTGIGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTTEDCGIALGQAFKEALGAVRGVKRFGSGFAPLDEALSRAVVDLSNRPYAVVELGL-QREKVGDLSCEMIPHFLESFAEASRITLHVDCLRGKNDHHRSESAFKALAVAIREATS--PNGTNDVPSTKGVL"
assert record.rounds[1].alignments[7].hsps[0].query_start==5
assert record.rounds[1].alignments[7].hsps[0].query_end==200
assert record.rounds[1].alignments[7].hsps[0].sbjct_start==6
assert record.rounds[1].alignments[7].hsps[0].sbjct_end==219
assert record.rounds[1].alignments[8].hsps[0].query=="MTRISHITRNTKETQIELSINLDG---------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[8].hsps[0].match=="M+R + I R T ET+I++++NLDG                           +   ++ TGIGFLDHM+  L  HS + L +   GD   + +D HH  EDV I+LG    + LG   G++R+G    P+DEAL    +D+S RP+ V    L   +K+G   TEM      + A   GIT+H++   G N HH  E  FK+ A A+K+A+S  ++   +IPS+KGVL"
assert record.rounds[1].alignments[8].hsps[0].sbjct=="MSREALINRITNETKIQIALNLDGGKLELKESIFPNQSIIIDEHHAKQVSGSQYINVQTGIGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTAEDVGISLGMAFKQALGQIKGVKRFGHGFAPLDEALSRAVVDLSNRPFAVIELGLK-REKIGDLSTEMIPHVLESFAGAVGITIHVDCLRGFNDHHRAESAFKALAIAIKEAIS--KTGKNDIPSTKGVL"
assert record.rounds[1].alignments[8].hsps[0].query_start==1
assert record.rounds[1].alignments[8].hsps[0].query_end==200
assert record.rounds[1].alignments[8].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[8].hsps[0].sbjct_end==221
assert record.rounds[1].alignments[9].hsps[0].query=="RISHITRNTKETQIELSINLDG------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[9].hsps[0].match=="R + I+R T ET+I+++I+L+G                        +   DI TG+GFLDHM+  L  HS + L +   GD   + +D HH  ED  IALG+   E +G   G++R+G+   P+DEAL    +D+S RP+ V    L   + +G   TEM   F  + A  A ITLH++   G N HH  E  FK+ A A+++A+S   +   ++PS+KGVL"
assert record.rounds[1].alignments[9].hsps[0].sbjct=="RKAFISRITNETKIQIAISLNGGYIQIKDSILPAKKDDDVASQATQSQVIDIHTGVGFLDHMIHALAKHSGWSLIVECIGD---LHIDDHHTTEDCGIALGQAFKEAMGAVRGVKRFGTGFAPLDEALSRAVVDLSNRPFAVIDLGLK-REMIGDLSTEMIPHFLESFAEAARITLHVDCLRGFNDHHRSESAFKALAVAIREAIS--SNGTNDVPSTKGVL"
assert record.rounds[1].alignments[9].hsps[0].query_start==3
assert record.rounds[1].alignments[9].hsps[0].query_end==200
assert record.rounds[1].alignments[9].hsps[0].sbjct_start==16
assert record.rounds[1].alignments[9].hsps[0].sbjct_end==231
assert record.rounds[1].alignments[10].hsps[0].query=="SHITRNTKETQIELSINLDGTGQA------------------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[10].hsps[0].match=="+ + R T+ET I+L+++LDG   +                   + TG+GFLDHMLT L  H  + L +   GD   + +D HH +ED  IALG+   E LG+  GI+R+G    P+DEAL    +D S RP+ V    L   +++G   TEM   F  + A    IT+H++   G N HH  E  FK+ A A+++A     +   ++PS+KGVL"
assert record.rounds[1].alignments[10].hsps[0].sbjct=="AFVKRVTQETNIQLALDLDGGSVSVRESILGKEYASGDGQTIHVHTGVGFLDHMLTALAKHGGWSLILECIGD---LHIDDHHTVEDCGIALGQAFKEALGSVRGIKRFGHGFAPLDEALSRAVVDFSNRPFAVVELGLK-RERIGQLSTEMIPHFLESFATEGRITMHVDCLRGTNDHHRSESGFKALAIAIREA--RTPTGRDDVPSTKGVL"
assert record.rounds[1].alignments[10].hsps[0].query_start==5
assert record.rounds[1].alignments[10].hsps[0].query_end==200
assert record.rounds[1].alignments[10].hsps[0].sbjct_start==2
assert record.rounds[1].alignments[10].hsps[0].sbjct_end==209
assert record.rounds[1].alignments[11].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[11].hsps[0].match=="R +H+ RNTKETQI++ + LD  G + I+TG+GF DHML  +  H  F ++I   GD   + +D HH +ED  +ALG+ +   LG+K GI R+G F +PMDE L  C LDISGRP+L + A+ +  Q++G   TEM E FFR+L++  G+TLHL    G+N HH +E +FK+  R L+QA+ ++      +PSSKGVL"
assert record.rounds[1].alignments[11].hsps[0].sbjct=="RYAHVVRNTKETQIDVQVWLDREGGSKINTGVGFFDHMLDQIATHGGFRMEINVKGD---LYIDDHHTVEDTGLALGEALKIALGDKRGICRFG-FVLPMDECLARCALDISGRPHLEYKAEFT-YQRVGDLSTEMIEHFFRSLSYTMGVTLHLKT-KGKNDHHRVESLFKAFGRTLRQAIRVEG---DTLPSSKGVL"
assert record.rounds[1].alignments[11].hsps[0].query_start==3
assert record.rounds[1].alignments[11].hsps[0].query_end==200
assert record.rounds[1].alignments[11].hsps[0].sbjct_start==167
assert record.rounds[1].alignments[11].hsps[0].sbjct_end==355
assert record.rounds[1].alignments[12].hsps[0].query=="SHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[12].hsps[0].match=="+ I RNT ET+I +++NLDGTG  D+ TG+GFLDHML  L+ HS  DL +   GD   V +D HH  E   IA+G+ +++ +G++ GI+RYG   +PMDE L    LD S RPYL++    S   K+G  DTE+  E+F+A A  AG+TLH+   YG+N HHI+E  +K+ ARAL+  + ID  K   +PS+KG L"
assert record.rounds[1].alignments[12].hsps[0].sbjct=="ASIERNTTETRIRVAVNLDGTGVYDVKTGVGFLDHMLEQLSRHSLMDLSVAAEGD---VHIDAHHTTEHSGIAIGQAVAKAVGDRKGIQRYGHAYVPMDETLTRVALDFSNRPYLIWKVSFS-RDKIGDMDTELFREWFQAFAMAAGVTLHVECLYGENNHHIVESCYKALARALRAGIEIDPRKRDAVPSTKGTL"
assert record.rounds[1].alignments[12].hsps[0].query_start==5
assert record.rounds[1].alignments[12].hsps[0].query_end==200
assert record.rounds[1].alignments[12].hsps[0].sbjct_start==12
assert record.rounds[1].alignments[12].hsps[0].sbjct_end==203
assert record.rounds[1].alignments[13].hsps[0].query=="RISHITRNTKETQIELSINLDG-----TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[13].hsps[0].match=="RI+ + R T ET I  +I+LD        + ++STGIGFLDHM T L  H    L++   GD   + +D HH  ED A+ALG+   + LG + GI+RYG    P+DE+L    +DIS RPY + H   +  +K+G   TEM     ++ AF AG+TLH++   G+N HHI E  FK+ A A++ A+S   +   ++PS+KGVL"
assert record.rounds[1].alignments[13].hsps[0].sbjct=="RIASVERTTSETHISCTIDLDHIPGVTEQKINVSTGIGFLDHMFTALAKHGGMSLQLQCKGD---LHIDDHHTAEDCALALGEAFKKALGERKGIKRYGYAYAPLDESLSRAVIDISSRPYFMCHLPFT-REKVGDLSTEMVSHLLQSFAFAAGVTLHIDSIRGENNHHIAESAFKALALAIRMAIS--RTGGDDVPSTKGVL"
assert record.rounds[1].alignments[13].hsps[0].query_start==3
assert record.rounds[1].alignments[13].hsps[0].query_end==200
assert record.rounds[1].alignments[13].hsps[0].sbjct_start==4
assert record.rounds[1].alignments[13].hsps[0].sbjct_end==200
assert record.rounds[1].alignments[14].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQA-------------------------DISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[14].hsps[0].match=="R S I R T ET+I+++++LDG   +                          ++TGIGFLDHML  L  H  + + I   GD   + +D HH  ED  IALG    E LG+  GI+R+GS   P+DEAL    +D+S RPY V    L   +K+G    EM      + A  A +T+H++   G N HH  E  FK+ A A+K+A+S   +   +IPS+KGVL"
assert record.rounds[1].alignments[14].hsps[0].sbjct=="RSSLIKRITNETKIQIALSLDGGPVSLAQSLFKDKDYSAEHAAQATSSQFISVNTGIGFLDHMLHALAKHGGWSVIIECVGD---LHIDDHHSAEDTGIALGMAFKEALGHVRGIKRFGSGFAPLDEALSRAVIDMSNRPYAVVDLGLK-REKIGDLSCEMIPHVLESFAQGAHVTMHVDCLRGFNDHHRAESAFKALAIAIKEAIS--SNGTDDIPSTKGVL"
assert record.rounds[1].alignments[14].hsps[0].query_start==3
assert record.rounds[1].alignments[14].hsps[0].query_end==200
assert record.rounds[1].alignments[14].hsps[0].sbjct_start==7
assert record.rounds[1].alignments[14].hsps[0].sbjct_end==223
assert record.rounds[1].alignments[15].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[15].hsps[0].match=="R + + R TKET I++ + LD  G  +I TG+GF DHML  +  H  F + +   GD   + +D HH +ED A+ALG+ + + +G+K GI R+G F +PMDE    C LD+SGRP++ F+A      K+G + TE+TE FF++LAF+   TLHLN   G N HH IE +FK+  R L+QA+ I+ +   E+PSSKGVL"
assert record.rounds[1].alignments[15].hsps[0].sbjct=="RFAEVIRQTKETDIKVQVWLDEAGVNEIKTGVGFFDHMLDQIATHGGFRMNVQCKGD---LWIDEHHTVEDTALALGQALKQAVGDKRGIARFG-FVLPMDECKAECALDLSGRPWIKFNACFK-RDKVGDFSTELTEHFFQSLAFSMLATLHLNV-TGNNDHHKIESLFKAFGRTLRQAIRIEGN---EMPSSKGVL"
assert record.rounds[1].alignments[15].hsps[0].query_start==3
assert record.rounds[1].alignments[15].hsps[0].query_end==200
assert record.rounds[1].alignments[15].hsps[0].sbjct_start==174
assert record.rounds[1].alignments[15].hsps[0].sbjct_end==362
assert record.rounds[1].alignments[16].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGI-TLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[16].hsps[0].match=="R ++I+R TKET I + ++LDGTG++ +S+GIGFLDHMLT L  HS FDL++   GD     +D HH  ED A+ LG+     LG++ GI R+GS  +P+DEAL    +DIS R +   +  L     +G   +EM   FF + A  A + TLH++   G+N HH  E  FK+ A AL+ AV  D +    +PS+KGVL"
assert record.rounds[1].alignments[16].hsps[0].sbjct=="REANISRVTKETSISVKLSLDGTGKSKVSSGIGFLDHMLTALAKHSRFDLELDCKGDT---WIDDHHTTEDCALTLGEAFDVALGDRAGIARFGSACVPLDEALSRAIVDISSRAHSEINLQLV-RPSVGELSSEMITHFFESFASAA-LXTLHVDVLRGRNDHHRAEASFKALAVALRTAVKHDAT--AGVPSTKGVL"
assert record.rounds[1].alignments[16].hsps[0].query_start==3
assert record.rounds[1].alignments[16].hsps[0].query_end==200
assert record.rounds[1].alignments[16].hsps[0].sbjct_start==260
assert record.rounds[1].alignments[16].hsps[0].sbjct_end==451
assert record.rounds[1].alignments[17].hsps[0].query=="RISHITRNTKETQIELSINLD-----------------------GTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[17].hsps[0].match=="R + + RNT ET+I ++I LD                       G     + TGIGFLDHM   L  H+ + L++   GD   + +D HH  ED AIALG    + +GN  G++R+G    P+DEAL    +D+SGRPY V    L   +K+G    EM      + +  AGITLH+   YG N HH  E  FKS A A++ A S+  +   E+PS+KGVL"
assert record.rounds[1].alignments[17].hsps[0].sbjct=="RRAFVERNTNETKISVAIALDKAPLPEESNFIDELITSKHANQKGEQVIQVDTGIGFLDHMYHALAKHAGWSLRLYSRGD---LIIDDHHTAEDTAIALGIAFKQAMGNFAGVKRFGHAYCPLDEALSRSVVDLSGRPYAVIDLGLK-REKVGELSCEMIPHLLYSFSVAAGITLHVTCLYGSNDHHRAESAFKSLAVAMRAATSL--TGSSEVPSTKGVL"
assert record.rounds[1].alignments[17].hsps[0].query_start==3
assert record.rounds[1].alignments[17].hsps[0].query_end==200
assert record.rounds[1].alignments[17].hsps[0].sbjct_start==2
assert record.rounds[1].alignments[17].hsps[0].sbjct_end==216
assert record.rounds[1].alignments[18].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[18].hsps[0].match=="R +H+ RNTKETQI++S+ LD  G + I+TG+GF DHML  +  H  F ++I   GD   + +D HH +ED  +AL + +   L +K GI R+G F +PMDE L  C LDISGRP+L + A+ +  Q++G   TEM E FFR+L++  G+TLHL    G+N HH +E +FK+  R ++QA+ ++      +PSSKGVL"
assert record.rounds[1].alignments[18].hsps[0].sbjct=="RYAHVVRNTKETQIDVSVWLDREGNSKINTGVGFFDHMLDQIATHGGFRMEITVKGD---LYIDDHHTVEDTGLALREALKLALRDKRGICRFG-FVLPMDECLA-CALDISGRPHLEYKAEFT-YQRVGNLSTEMIEHFFRSLSYTMGVTLHLKT-KGKNDHHRVESLFKAFGRTVRQAIRVEG---DTLPSSKGVL"
assert record.rounds[1].alignments[18].hsps[0].query_start==3
assert record.rounds[1].alignments[18].hsps[0].query_end==200
assert record.rounds[1].alignments[18].hsps[0].sbjct_start==167
assert record.rounds[1].alignments[18].hsps[0].sbjct_end==354
assert record.rounds[1].alignments[19].hsps[0].query=="TRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVF--------HADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[19].hsps[0].match=="+R + I R T+E+ I + ++LDGTGQ  + TG+ F DHMLT L  H+ FDL +   GD E   ++ HH IED AIALG  + + LG+K GIRR+G   IPMDE L    +D+SGRPY V         H  ++G+     Y T +    F +LA NA I LH+   YG++ HHI E  +K+ ARAL+QAV  D   V  +PS+KG L"
assert record.rounds[1].alignments[19].hsps[0].sbjct=="SRRARIERRTRESDIVIELDLDGTGQVAVDTGVPFYDHMLTALGSHASFDLTVRATGDVE---IEAHHTIEDTAIALGTALGQALGDKRGIRRFGDAFIPMDETLAHAAVDLSGRPYCVHTGEPDHLQHTTIAGSSV--PYHTVINRHVFESLAANARIALHVRVLYGRDPHHITEAQYKAVARALRQAVEPDPR-VSGVPSTKGAL"
assert record.rounds[1].alignments[19].hsps[0].query_start==2
assert record.rounds[1].alignments[19].hsps[0].query_end==200
assert record.rounds[1].alignments[19].hsps[0].sbjct_start==10
assert record.rounds[1].alignments[19].hsps[0].sbjct_end==210
assert record.rounds[1].alignments[20].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[20].hsps[0].match=="R S  TR T ET +++ + +DG+G++ ++TG+GFLDHML  +  H   DL++   GD E   +D HH +EDVA+ LG+ + E LG+K GIRR     +PMD+AL T  LD+SGRPY V   +   +  +G   ++    F  +LA +A + +H +   G+N HH  E +FK+ A A++ AV ++    GEIPS+KG L"
assert record.rounds[1].alignments[20].hsps[0].sbjct=="RRSMKTRETLETHVKVDLEIDGSGKSSVNTGLGFLDHMLESVARHGLLDLEVEARGDLE---VDDHHTVEDVALTLGEALREALGDKSGIRRMAHAMVPMDDALATVALDLSGRPYTVLELEFD-DAVIGDVKSQNIGHFIESLAVSAAMNIHASV-RGRNDHHKAEALFKALALAIRDAVRVEH---GEIPSTKGKL"
assert record.rounds[1].alignments[20].hsps[0].query_start==3
assert record.rounds[1].alignments[20].hsps[0].query_end==200
assert record.rounds[1].alignments[20].hsps[0].sbjct_start==5
assert record.rounds[1].alignments[20].hsps[0].sbjct_end==194
assert record.rounds[1].alignments[21].hsps[0].query=="RISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKG"
assert record.rounds[1].alignments[21].hsps[0].match=="RI  + R TKET I L IN+DGTG+  I TGI F DH+L     H  FDL +   GD E   +D HH +EDV I LG  +++    K  I R+G   IPMD+A  T  +D+SGR Y V + +    + +G   TE    FF ++A    + +H     G+N HH  E +FK+   AL  A  IDE K   + S+KG"
assert record.rounds[1].alignments[21].hsps[0].sbjct=="RIFEVMRETKETNIYLKINIDGTGKYKIDTGIPFFDHLLASFAKHGCFDLIVKARGDLE---IDDHHTVEDVGICLGLALNQI--EKRNIFRFGWAIIPMDDARATVAIDLSGRSYCVGNYE-PKREFVGDLATENINHFFESVASYGMLNIHYEVI-GKNEHHKAEALFKAFGVALDLATKIDERK--GVISTKG"
assert record.rounds[1].alignments[21].hsps[0].query_start==3
assert record.rounds[1].alignments[21].hsps[0].query_end==198
assert record.rounds[1].alignments[21].hsps[0].sbjct_start==7
assert record.rounds[1].alignments[21].hsps[0].sbjct_end==193
assert record.rounds[1].alignments[22].hsps[0].query=="RISHITRNTKETQIELSINLDG------------------------TGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALK"
assert record.rounds[1].alignments[22].hsps[0].match=="R + ++R+T ET I++++++DG                        +    I+TGIGFLDHML  L  H+ + + +   GD   + +D HH  ED  IA+G   ++ LG   G+ R+G    P+DEAL    +D+S RPY V    L   +KLG    EM     ++ A  A ITLH++   G N HH  E  FK+ A A++"
assert record.rounds[1].alignments[22].hsps[0].sbjct=="RAAALSRDTNETSIQIALSIDGGELPQDADPRLLEASSAHASQTSKSQVISINTGIGFLDHMLHALAKHAGWSMALNCKGD---LHIDDHHTAEDCCIAVGTTFAKALGALTGVARFGYAYAPLDEALSRAVVDLSNRPYTVVDLGLK-REKLGELSCEMIPHCLQSFAQAARITLHVDCLRGDNDHHRAESAFKALAVAVR"
assert record.rounds[1].alignments[22].hsps[0].query_start==3
assert record.rounds[1].alignments[22].hsps[0].query_end==180
assert record.rounds[1].alignments[22].hsps[0].sbjct_start==8
assert record.rounds[1].alignments[22].hsps[0].sbjct_end==205
assert record.rounds[1].alignments[23].hsps[0].query=="MTRISHITRNTKETQIELSINLDGTGQADISTGIGFLDHMLTLLTFHSDFDLKIIGHGDHETVGMDPHHLIEDVAIALGKCISEDLGNKLGIRRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEEFFRALAFNAGITLHLNEHYGQNTHHIIEGMFKSTARALKQAVSIDESKVGEIPSSKGVL"
assert record.rounds[1].alignments[23].hsps[0].match=="M+R ++ITR TKET+IE+ +++D  G+  +ST I F +HML  L  + +    +      + +  D HH++EDVAI LG  I   LG+K GI+R+    IPMD+ALV   LDIS R     + +L  ++ +GG  TE    FF++ A+N+GITLH+++  G NTHHIIE  FK+   AL +A  I ++   EI S+KG++"
assert record.rounds[1].alignments[23].hsps[0].sbjct=="MSRSANITRETKETKIEVLLDIDRKGEVKVSTPIPFFNHMLITLLTYMNSTAIVSAT---DKLPYDDHHIVEDVAITLGLAIKTALGDKRGIKRFSHQIIPMDDALVLVSLDISNRGMAFVNLNLKRSE-IGGLATENVPHFFQSFAYNSGITLHISQLSGYNTHHIIEASFKALGLALYEATRIVDN---EIRSTKGII"
assert record.rounds[1].alignments[23].hsps[0].query_start==1
assert record.rounds[1].alignments[23].hsps[0].query_end==200
assert record.rounds[1].alignments[23].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[23].hsps[0].sbjct_end==193
assert record.rounds[1].alignments[24].hsps[0].query=="LTLLTFHSDFDLKIIGHGDHETVGMDPH-HLIEDVAIAL------GKCISEDLG----NKLGI--------RRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEE"
assert record.rounds[1].alignments[24].hsps[0].match=="L  +  +    + +       T+ M  H ++++  A+        G    E LG    +   I                PM +A  +  +D+   P L+F        + G + +++   "
assert record.rounds[1].alignments[24].hsps[0].sbjct=="LAEVAKNRGIAIPVDLDSQVNTLFMKSHSNMVQRAAMGWRLSARSGPRFKEALGGPAWDYRNIIEKLQDVVASLEHQFSPMMQAEFSVLVDVLYSPELLFPEGSDARIRCGAFMSKLINH"
assert record.rounds[1].alignments[24].hsps[0].query_start==41
assert record.rounds[1].alignments[24].hsps[0].query_end==141
assert record.rounds[1].alignments[24].hsps[0].sbjct_start==1540
assert record.rounds[1].alignments[24].hsps[0].sbjct_end==1659
assert record.rounds[1].alignments[25].hsps[0].query=="LTLLTFHSDFDLKIIGHGDHETVGMDPHH-LIEDVAIAL------GKCISEDLG----NKLGI--------RRYGSFTIPMDEALVTCDLDISGRPYLVFHADLSGNQKLGGYDTEMTEE"
assert record.rounds[1].alignments[25].hsps[0].match=="L  +  +    + +       T+ M  H   ++  A+        G    E LG    +   I                PM +A  +  +D+   P L+F        + G + +++   "
assert record.rounds[1].alignments[25].hsps[0].sbjct=="LAEVAKNRGIAIPVDLDSQVNTLFMKNHSSTVQRAAMGWRLSARSGPRFKEALGGPAWDYRNIIEKLQDVVASLEQQFSPMMQAEFSVLVDVLYSPELLFPEGSDARIRCGAFMSKLINH"
assert record.rounds[1].alignments[25].hsps[0].query_start==41
assert record.rounds[1].alignments[25].hsps[0].query_end==141
assert record.rounds[1].alignments[25].hsps[0].sbjct_start==1540
assert record.rounds[1].alignments[25].hsps[0].sbjct_end==1659
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.ka_params==[0.317, 0.149, 0.479]
assert record.ka_params_gap==[0.270, 0.0524, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==26734132
assert record.num_sequences==82258
assert record.num_extends==1229035
assert record.num_good_extends==2616
assert record.num_seqs_better_e==56
assert record.hsps_no_gap==48
assert record.hsps_prelim_gapped==8
assert record.hsps_gapped==58
assert record.query_length==200
assert record.database_length==29652561
assert record.effective_hsp_length==50
assert record.effective_query_length==150
assert record.effective_database_length==25539661
assert record.effective_search_space==3830949150
assert record.effective_search_space_used==3830949150
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.5)
assert record.blast_cutoff==(63,28.8)

handle = open('Blast/bt048')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000" 
assert record.reference==reference
assert record.query=="gi|1348916|gb|G26684|G26684 human STS\nSTS_D11570.\x01gi|1375195|gb|G26945|G26945 human STS SHGC-32699."
assert record.query_letters==285
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==23
assert record.descriptions[0].title=="gi|1348916|gb|G26684|G26684 human STS STS_D11570. >gi|1375195|g..."
assert record.descriptions[0].score==517
assert record.descriptions[0].e==1e-146
assert record.descriptions[1].title=="gi|4516686|dbj|AU026763.1|AU026763 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[1].score==32
assert record.descriptions[1].e==1.600000
assert record.descriptions[2].title=="gi|6120827|gb|G55508.1|G55508 SHGC-100856 Human Homo sapiens ST..."
assert record.descriptions[2].score==32
assert record.descriptions[2].e==1.600000
assert record.descriptions[3].title=="gi|720683|gb|G03725|G03725 human STS WI-344."
assert record.descriptions[3].score==32
assert record.descriptions[3].e==1.600000
assert record.descriptions[4].title=="gi|5690111|gb|G54226.1|G54226 B124N23/SP6 Human Chromosome 12 H..."
assert record.descriptions[4].score==30
assert record.descriptions[4].e==6.500000
assert record.descriptions[5].title=="gi|4493307|gb|G47007.1|G47007 Z15259_1 Zebrafish AB Danio rerio..."
assert record.descriptions[5].score==30
assert record.descriptions[5].e==6.500000
assert record.descriptions[6].title=="gi|4491799|gb|G45508.1|G45508 Z24506_1 Zebrafish AB Danio rerio..."
assert record.descriptions[6].score==30
assert record.descriptions[6].e==6.500000
assert record.descriptions[7].title=="gi|6121596|gb|G56277.1|G56277 SHGC-101791 Human Homo sapiens ST..."
assert record.descriptions[7].score==30
assert record.descriptions[7].e==6.500000
assert record.descriptions[8].title=="gi|5222417|gb|G51240.1|G51240 SHGC-80720 Human Homo sapiens STS..."
assert record.descriptions[8].score==30
assert record.descriptions[8].e==6.500000
assert record.descriptions[9].title=="gi|5221977|gb|G50800.1|G50800 SHGC-83850 Human Homo sapiens STS..."
assert record.descriptions[9].score==30
assert record.descriptions[9].e==6.500000
assert record.descriptions[10].title=="gi|5224501|gb|G53324.1|G53324 SHGC-82315 Human Homo sapiens STS..."
assert record.descriptions[10].score==30
assert record.descriptions[10].e==6.500000
assert record.descriptions[11].title=="gi|4529247|gb|G48587.1|G48587 SHGC-82546 Human Homo sapiens STS..."
assert record.descriptions[11].score==30
assert record.descriptions[11].e==6.500000
assert record.descriptions[12].title=="gi|3359917|gb|G40708|G40708 Z8947 Zebrafish AB Danio rerio STS ..."
assert record.descriptions[12].score==30
assert record.descriptions[12].e==6.500000
assert record.descriptions[13].title=="gi|3359244|gb|G40035|G40035 Z13538 Zebrafish AB Danio rerio STS..."
assert record.descriptions[13].score==30
assert record.descriptions[13].e==6.500000
assert record.descriptions[14].title=="gi|1347715|gb|G25483|G25483 human STS EST334642."
assert record.descriptions[14].score==30
assert record.descriptions[14].e==6.500000
assert record.descriptions[15].title=="gi|1244262|gb|G19475|G19475 human STS SHGC-18755."
assert record.descriptions[15].score==30
assert record.descriptions[15].e==6.500000
assert record.descriptions[16].title=="gi|1232611|emb|Z51311|HS302WC9 H.sapiens (D5S2069) DNA segment ..."
assert record.descriptions[16].score==30
assert record.descriptions[16].e==6.500000
assert record.descriptions[17].title=="gi|1223022|gb|G18565|G18565 BMS485 cow Bos taurus STS genomic, ..."
assert record.descriptions[17].score==30
assert record.descriptions[17].e==6.500000
assert record.descriptions[18].title=="gi|1161779|gb|G15890|G15890 human STS CHLC.UTR_01448_M84721.P56..."
assert record.descriptions[18].score==30
assert record.descriptions[18].e==6.500000
assert record.descriptions[19].title=="gi|858803|gb|G05558|G05558 human STS WI-7105."
assert record.descriptions[19].score==30
assert record.descriptions[19].e==6.500000
assert record.descriptions[20].title=="gi|1342455|gb|G22129|G22129 human STS WI-14200."
assert record.descriptions[20].score==30
assert record.descriptions[20].e==6.500000
assert record.descriptions[21].title=="gi|1347001|gb|G24769|G24769 human STS EST129834."
assert record.descriptions[21].score==30
assert record.descriptions[21].e==6.500000
assert record.descriptions[22].title=="gi|605469|gb|L31223|HUMUT821B Human STS UT821, 3' primer bind."
assert record.descriptions[22].score==30
assert record.descriptions[22].e==6.500000
assert len(record.alignments)==23
assert record.alignments[0].title==">gi|1348916|gb|G26684|G26684 human STS STS_D11570. >gi|1375195|gb|G26945|G26945 human STS SHGC-32699."
assert record.alignments[0].length==285
assert record.alignments[1].title==">gi|4516686|dbj|AU026763.1|AU026763 Rattus norvegicus, OTSUKA clone, OT33.16/752f07, microsatellite sequence, sequence tagged site"
assert record.alignments[1].length==307
assert record.alignments[2].title==">gi|6120827|gb|G55508.1|G55508 SHGC-100856 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[2].length==711
assert record.alignments[3].title==">gi|720683|gb|G03725|G03725 human STS WI-344."
assert record.alignments[3].length==246
assert record.alignments[4].title==">gi|5690111|gb|G54226.1|G54226 B124N23/SP6 Human Chromosome 12 Homo sapiens STS genomic clone RPCI-11-B124N23 SP6, sequence tagged site"
assert record.alignments[4].length==550
assert record.alignments[5].title==">gi|4493307|gb|G47007.1|G47007 Z15259_1 Zebrafish AB Danio rerio STS genomic clone Z15259 5', sequence tagged site"
assert record.alignments[5].length==442
assert record.alignments[6].title==">gi|4491799|gb|G45508.1|G45508 Z24506_1 Zebrafish AB Danio rerio STS genomic clone Z24506 5', sequence tagged site"
assert record.alignments[6].length==272
assert record.alignments[7].title==">gi|6121596|gb|G56277.1|G56277 SHGC-101791 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[7].length==641
assert record.alignments[8].title==">gi|5222417|gb|G51240.1|G51240 SHGC-80720 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[8].length==712
assert record.alignments[9].title==">gi|5221977|gb|G50800.1|G50800 SHGC-83850 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[9].length==422
assert record.alignments[10].title==">gi|5224501|gb|G53324.1|G53324 SHGC-82315 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[10].length==428
assert record.alignments[11].title==">gi|4529247|gb|G48587.1|G48587 SHGC-82546 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[11].length==694
assert record.alignments[12].title==">gi|3359917|gb|G40708|G40708 Z8947 Zebrafish AB Danio rerio STS genomic"
assert record.alignments[12].length==549
assert record.alignments[13].title==">gi|3359244|gb|G40035|G40035 Z13538 Zebrafish AB Danio rerio STS genomic"
assert record.alignments[13].length==536
assert record.alignments[14].title==">gi|1347715|gb|G25483|G25483 human STS EST334642."
assert record.alignments[14].length==407
assert record.alignments[15].title==">gi|1244262|gb|G19475|G19475 human STS SHGC-18755."
assert record.alignments[15].length==400
assert record.alignments[16].title==">gi|1232611|emb|Z51311|HS302WC9 H.sapiens (D5S2069) DNA segment containing (CA) repeat; clone AFM302wc9; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[16].length==374
assert record.alignments[17].title==">gi|1223022|gb|G18565|G18565 BMS485 cow Bos taurus STS genomic, sequence tagged site [Bos taurus]"
assert record.alignments[17].length==181
assert record.alignments[18].title==">gi|1161779|gb|G15890|G15890 human STS CHLC.UTR_01448_M84721.P56085 clone UTR_01448_M84721."
assert record.alignments[18].length==729
assert record.alignments[19].title==">gi|858803|gb|G05558|G05558 human STS WI-7105."
assert record.alignments[19].length==735
assert record.alignments[20].title==">gi|1342455|gb|G22129|G22129 human STS WI-14200."
assert record.alignments[20].length==373
assert record.alignments[21].title==">gi|1347001|gb|G24769|G24769 human STS EST129834."
assert record.alignments[21].length==306
assert record.alignments[22].title==">gi|605469|gb|L31223|HUMUT821B Human STS UT821, 3' primer bind."
assert record.alignments[22].length==127
assert record.alignments[0].hsps[0].score==261
assert record.alignments[0].hsps[0].bits==517
assert record.alignments[0].hsps[0].expect==1e-146
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==16
assert record.alignments[1].hsps[0].bits==32.2
assert record.alignments[1].hsps[0].expect==1.6
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==16
assert record.alignments[2].hsps[0].bits==32.2
assert record.alignments[2].hsps[0].expect==1.6
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==16
assert record.alignments[3].hsps[0].bits==32.2
assert record.alignments[3].hsps[0].expect==1.6
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==15
assert record.alignments[4].hsps[0].bits==30.2
assert record.alignments[4].hsps[0].expect==6.5
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==15
assert record.alignments[5].hsps[0].bits==30.2
assert record.alignments[5].hsps[0].expect==6.5
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==15
assert record.alignments[6].hsps[0].bits==30.2
assert record.alignments[6].hsps[0].expect==6.5
assert len(record.alignments[6].hsps)==1
assert record.alignments[7].hsps[0].score==15
assert record.alignments[7].hsps[0].bits==30.2
assert record.alignments[7].hsps[0].expect==6.5
assert len(record.alignments[7].hsps)==1
assert record.alignments[8].hsps[0].score==15
assert record.alignments[8].hsps[0].bits==30.2
assert record.alignments[8].hsps[0].expect==6.5
assert len(record.alignments[8].hsps)==1
assert record.alignments[9].hsps[0].score==15
assert record.alignments[9].hsps[0].bits==30.2
assert record.alignments[9].hsps[0].expect==6.5
assert len(record.alignments[9].hsps)==1
assert record.alignments[10].hsps[0].score==15
assert record.alignments[10].hsps[0].bits==30.2
assert record.alignments[10].hsps[0].expect==6.5
assert len(record.alignments[10].hsps)==1
assert record.alignments[11].hsps[0].score==15
assert record.alignments[11].hsps[0].bits==30.2
assert record.alignments[11].hsps[0].expect==6.5
assert len(record.alignments[11].hsps)==1
assert record.alignments[12].hsps[0].score==15
assert record.alignments[12].hsps[0].bits==30.2
assert record.alignments[12].hsps[0].expect==6.5
assert len(record.alignments[12].hsps)==1
assert record.alignments[13].hsps[0].score==15
assert record.alignments[13].hsps[0].bits==30.2
assert record.alignments[13].hsps[0].expect==6.5
assert len(record.alignments[13].hsps)==1
assert record.alignments[14].hsps[0].score==15
assert record.alignments[14].hsps[0].bits==30.2
assert record.alignments[14].hsps[0].expect==6.5
assert len(record.alignments[14].hsps)==1
assert record.alignments[15].hsps[0].score==15
assert record.alignments[15].hsps[0].bits==30.2
assert record.alignments[15].hsps[0].expect==6.5
assert len(record.alignments[15].hsps)==1
assert record.alignments[16].hsps[0].score==15
assert record.alignments[16].hsps[0].bits==30.2
assert record.alignments[16].hsps[0].expect==6.5
assert len(record.alignments[16].hsps)==1
assert record.alignments[17].hsps[0].score==15
assert record.alignments[17].hsps[0].bits==30.2
assert record.alignments[17].hsps[0].expect==6.5
assert len(record.alignments[17].hsps)==1
assert record.alignments[18].hsps[0].score==15
assert record.alignments[18].hsps[0].bits==30.2
assert record.alignments[18].hsps[0].expect==6.5
assert len(record.alignments[18].hsps)==1
assert record.alignments[19].hsps[0].score==15
assert record.alignments[19].hsps[0].bits==30.2
assert record.alignments[19].hsps[0].expect==6.5
assert len(record.alignments[19].hsps)==1
assert record.alignments[20].hsps[0].score==15
assert record.alignments[20].hsps[0].bits==30.2
assert record.alignments[20].hsps[0].expect==6.5
assert len(record.alignments[20].hsps)==1
assert record.alignments[21].hsps[0].score==15
assert record.alignments[21].hsps[0].bits==30.2
assert record.alignments[21].hsps[0].expect==6.5
assert len(record.alignments[21].hsps)==1
assert record.alignments[22].hsps[0].score==15
assert record.alignments[22].hsps[0].bits==30.2
assert record.alignments[22].hsps[0].expect==6.5
assert len(record.alignments[22].hsps)==1
assert record.alignments[0].hsps[0].identities==(285, 285)
assert record.alignments[1].hsps[0].identities==(16, 16)
assert record.alignments[2].hsps[0].identities==(18, 19)
assert record.alignments[3].hsps[0].identities==(16, 16)
assert record.alignments[4].hsps[0].identities==(15, 15)
assert record.alignments[5].hsps[0].identities==(15, 15)
assert record.alignments[6].hsps[0].identities==(15, 15)
assert record.alignments[7].hsps[0].identities==(15, 15)
assert record.alignments[8].hsps[0].identities==(17, 18)
assert record.alignments[9].hsps[0].identities==(15, 15)
assert record.alignments[10].hsps[0].identities==(15, 15)
assert record.alignments[11].hsps[0].identities==(17, 18)
assert record.alignments[12].hsps[0].identities==(15, 15)
assert record.alignments[13].hsps[0].identities==(15, 15)
assert record.alignments[14].hsps[0].identities==(15, 15)
assert record.alignments[15].hsps[0].identities==(15, 15)
assert record.alignments[16].hsps[0].identities==(15, 15)
assert record.alignments[17].hsps[0].identities==(15, 15)
assert record.alignments[18].hsps[0].identities==(18, 19)
assert record.alignments[19].hsps[0].identities==(15, 15)
assert record.alignments[20].hsps[0].identities==(15, 15)
assert record.alignments[21].hsps[0].identities==(15, 15)
assert record.alignments[22].hsps[0].identities==(15, 15)
assert record.alignments[0].hsps[0].strand==("Plus", "Plus")
assert record.alignments[1].hsps[0].strand==("Plus", "Plus")
assert record.alignments[2].hsps[0].strand==("Plus", "Plus")
assert record.alignments[3].hsps[0].strand==("Plus", "Minus")
assert record.alignments[4].hsps[0].strand==("Plus", "Plus")
assert record.alignments[5].hsps[0].strand==("Plus", "Minus")
assert record.alignments[6].hsps[0].strand==("Plus", "Minus")
assert record.alignments[7].hsps[0].strand==("Plus", "Minus")
assert record.alignments[8].hsps[0].strand==("Plus", "Plus")
assert record.alignments[9].hsps[0].strand==("Plus", "Minus")
assert record.alignments[10].hsps[0].strand==("Plus", "Minus")
assert record.alignments[11].hsps[0].strand==("Plus", "Plus")
assert record.alignments[12].hsps[0].strand==("Plus", "Minus")
assert record.alignments[13].hsps[0].strand==("Plus", "Minus")
assert record.alignments[14].hsps[0].strand==("Plus", "Plus")
assert record.alignments[15].hsps[0].strand==("Plus", "Minus")
assert record.alignments[16].hsps[0].strand==("Plus", "Plus")
assert record.alignments[17].hsps[0].strand==("Plus", "Plus")
assert record.alignments[18].hsps[0].strand==("Plus", "Plus")
assert record.alignments[19].hsps[0].strand==("Plus", "Minus")
assert record.alignments[20].hsps[0].strand==("Plus", "Minus")
assert record.alignments[21].hsps[0].strand==("Plus", "Plus")
assert record.alignments[22].hsps[0].strand==("Plus", "Minus")
assert record.alignments[0].hsps[0].query=="gatccctacccttnccgttggtctctntcgctgactcgaggcacctaacatccattcacacccaacacaggccagcgacttctggggctcagccacagacatggtttgtnactnttgagcttctgttcctagagaatcctagaggcttgattggcccaggctgctgtntgtnctggaggcaaagaatccctacctcctaggggtgaaaggaaatnaaaatggaaagttcttgtagcgcaaggcctgacatgggtagctgctcaataaatgctagtntgttatttc"
assert record.alignments[0].hsps[0].match=="|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
assert record.alignments[0].hsps[0].sbjct=="gatccctacccttnccgttggtctctntcgctgactcgaggcacctaacatccattcacacccaacacaggccagcgacttctggggctcagccacagacatggtttgtnactnttgagcttctgttcctagagaatcctagaggcttgattggcccaggctgctgtntgtnctggaggcaaagaatccctacctcctaggggtgaaaggaaatnaaaatggaaagttcttgtagcgcaaggcctgacatgggtagctgctcaataaatgctagtntgttatttc"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==285
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==285
assert record.alignments[1].hsps[0].query=="ggaaagttcttgtagc"
assert record.alignments[1].hsps[0].match=="||||||||||||||||"
assert record.alignments[1].hsps[0].sbjct=="ggaaagttcttgtagc"
assert record.alignments[1].hsps[0].query_start==221
assert record.alignments[1].hsps[0].query_end==236
assert record.alignments[1].hsps[0].sbjct_start==32
assert record.alignments[1].hsps[0].sbjct_end==47
assert record.alignments[2].hsps[0].query=="gaaatnaaaatggaaagtt"
assert record.alignments[2].hsps[0].match=="||||| |||||||||||||"
assert record.alignments[2].hsps[0].sbjct=="gaaataaaaatggaaagtt"
assert record.alignments[2].hsps[0].query_start==210
assert record.alignments[2].hsps[0].query_end==228
assert record.alignments[2].hsps[0].sbjct_start==588
assert record.alignments[2].hsps[0].sbjct_end==606
assert record.alignments[3].hsps[0].query=="ctcaataaatgctagt"
assert record.alignments[3].hsps[0].match=="||||||||||||||||"
assert record.alignments[3].hsps[0].sbjct=="ctcaataaatgctagt"
assert record.alignments[3].hsps[0].query_start==260
assert record.alignments[3].hsps[0].query_end==275
assert record.alignments[3].hsps[0].sbjct_start==178
assert record.alignments[3].hsps[0].sbjct_end==163
assert record.alignments[4].hsps[0].query=="aaaatggaaagttct"
assert record.alignments[4].hsps[0].match=="|||||||||||||||"
assert record.alignments[4].hsps[0].sbjct=="aaaatggaaagttct"
assert record.alignments[4].hsps[0].query_start==216
assert record.alignments[4].hsps[0].query_end==230
assert record.alignments[4].hsps[0].sbjct_start==330
assert record.alignments[4].hsps[0].sbjct_end==344
assert record.alignments[5].hsps[0].query=="ttctgttcctagaga"
assert record.alignments[5].hsps[0].match=="|||||||||||||||"
assert record.alignments[5].hsps[0].sbjct=="ttctgttcctagaga"
assert record.alignments[5].hsps[0].query_start==121
assert record.alignments[5].hsps[0].query_end==135
assert record.alignments[5].hsps[0].sbjct_start==384
assert record.alignments[5].hsps[0].sbjct_end==370
assert record.alignments[6].hsps[0].query=="ggaaagttcttgtag"
assert record.alignments[6].hsps[0].match=="|||||||||||||||"
assert record.alignments[6].hsps[0].sbjct=="ggaaagttcttgtag"
assert record.alignments[6].hsps[0].query_start==221
assert record.alignments[6].hsps[0].query_end==235
assert record.alignments[6].hsps[0].sbjct_start==138
assert record.alignments[6].hsps[0].sbjct_end==124
assert record.alignments[7].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[7].hsps[0].match=="|||||||||||||||"
assert record.alignments[7].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[7].hsps[0].query_start==258
assert record.alignments[7].hsps[0].query_end==272
assert record.alignments[7].hsps[0].sbjct_start==216
assert record.alignments[7].hsps[0].sbjct_end==202
assert record.alignments[8].hsps[0].query=="taaatgctagtntgttat"
assert record.alignments[8].hsps[0].match=="||||||||||| ||||||"
assert record.alignments[8].hsps[0].sbjct=="taaatgctagtttgttat"
assert record.alignments[8].hsps[0].query_start==265
assert record.alignments[8].hsps[0].query_end==282
assert record.alignments[8].hsps[0].sbjct_start==293
assert record.alignments[8].hsps[0].sbjct_end==310
assert record.alignments[9].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[9].hsps[0].match=="|||||||||||||||"
assert record.alignments[9].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[9].hsps[0].query_start==258
assert record.alignments[9].hsps[0].query_end==272
assert record.alignments[9].hsps[0].sbjct_start==32
assert record.alignments[9].hsps[0].sbjct_end==18
assert record.alignments[10].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[10].hsps[0].match=="|||||||||||||||"
assert record.alignments[10].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[10].hsps[0].query_start==258
assert record.alignments[10].hsps[0].query_end==272
assert record.alignments[10].hsps[0].sbjct_start==47
assert record.alignments[10].hsps[0].sbjct_end==33
assert record.alignments[11].hsps[0].query=="taaatgctagtntgttat"
assert record.alignments[11].hsps[0].match=="||||||||||| ||||||"
assert record.alignments[11].hsps[0].sbjct=="taaatgctagtttgttat"
assert record.alignments[11].hsps[0].query_start==265
assert record.alignments[11].hsps[0].query_end==282
assert record.alignments[11].hsps[0].sbjct_start==292
assert record.alignments[11].hsps[0].sbjct_end==309
assert record.alignments[12].hsps[0].query=="aacatccattcacac"
assert record.alignments[12].hsps[0].match=="|||||||||||||||"
assert record.alignments[12].hsps[0].sbjct=="aacatccattcacac"
assert record.alignments[12].hsps[0].query_start==47
assert record.alignments[12].hsps[0].query_end==61
assert record.alignments[12].hsps[0].sbjct_start==479
assert record.alignments[12].hsps[0].sbjct_end==465
assert record.alignments[13].hsps[0].query=="ttctgttcctagaga"
assert record.alignments[13].hsps[0].match=="|||||||||||||||"
assert record.alignments[13].hsps[0].sbjct=="ttctgttcctagaga"
assert record.alignments[13].hsps[0].query_start==121
assert record.alignments[13].hsps[0].query_end==135
assert record.alignments[13].hsps[0].sbjct_start==433
assert record.alignments[13].hsps[0].sbjct_end==419
assert record.alignments[14].hsps[0].query=="ctaacatccattcac"
assert record.alignments[14].hsps[0].match=="|||||||||||||||"
assert record.alignments[14].hsps[0].sbjct=="ctaacatccattcac"
assert record.alignments[14].hsps[0].query_start==45
assert record.alignments[14].hsps[0].query_end==59
assert record.alignments[14].hsps[0].sbjct_start==389
assert record.alignments[14].hsps[0].sbjct_end==403
assert record.alignments[15].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[15].hsps[0].match=="|||||||||||||||"
assert record.alignments[15].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[15].hsps[0].query_start==258
assert record.alignments[15].hsps[0].query_end==272
assert record.alignments[15].hsps[0].sbjct_start==324
assert record.alignments[15].hsps[0].sbjct_end==310
assert record.alignments[16].hsps[0].query=="acagacatggtttgt"
assert record.alignments[16].hsps[0].match=="|||||||||||||||"
assert record.alignments[16].hsps[0].sbjct=="acagacatggtttgt"
assert record.alignments[16].hsps[0].query_start==95
assert record.alignments[16].hsps[0].query_end==109
assert record.alignments[16].hsps[0].sbjct_start==246
assert record.alignments[16].hsps[0].sbjct_end==260
assert record.alignments[17].hsps[0].query=="ctcaataaatgctag"
assert record.alignments[17].hsps[0].match=="|||||||||||||||"
assert record.alignments[17].hsps[0].sbjct=="ctcaataaatgctag"
assert record.alignments[17].hsps[0].query_start==260
assert record.alignments[17].hsps[0].query_end==274
assert record.alignments[17].hsps[0].sbjct_start==145
assert record.alignments[17].hsps[0].sbjct_end==159
assert record.alignments[18].hsps[0].query=="gtagctgctcaataaatgc"
assert record.alignments[18].hsps[0].match=="|||| ||||||||||||||"
assert record.alignments[18].hsps[0].sbjct=="gtaggtgctcaataaatgc"
assert record.alignments[18].hsps[0].query_start==253
assert record.alignments[18].hsps[0].query_end==271
assert record.alignments[18].hsps[0].sbjct_start==698
assert record.alignments[18].hsps[0].sbjct_end==716
assert record.alignments[19].hsps[0].query=="gaaagttcttgtagc"
assert record.alignments[19].hsps[0].match=="|||||||||||||||"
assert record.alignments[19].hsps[0].sbjct=="gaaagttcttgtagc"
assert record.alignments[19].hsps[0].query_start==222
assert record.alignments[19].hsps[0].query_end==236
assert record.alignments[19].hsps[0].sbjct_start==543
assert record.alignments[19].hsps[0].sbjct_end==529
assert record.alignments[20].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[20].hsps[0].match=="|||||||||||||||"
assert record.alignments[20].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[20].hsps[0].query_start==258
assert record.alignments[20].hsps[0].query_end==272
assert record.alignments[20].hsps[0].sbjct_start==33
assert record.alignments[20].hsps[0].sbjct_end==19
assert record.alignments[21].hsps[0].query=="tggaaagttcttgta"
assert record.alignments[21].hsps[0].match=="|||||||||||||||"
assert record.alignments[21].hsps[0].sbjct=="tggaaagttcttgta"
assert record.alignments[21].hsps[0].query_start==220
assert record.alignments[21].hsps[0].query_end==234
assert record.alignments[21].hsps[0].sbjct_start==144
assert record.alignments[21].hsps[0].sbjct_end==158
assert record.alignments[22].hsps[0].query=="acagacatggtttgt"
assert record.alignments[22].hsps[0].match=="|||||||||||||||"
assert record.alignments[22].hsps[0].sbjct=="acagacatggtttgt"
assert record.alignments[22].hsps[0].query_start==95
assert record.alignments[22].hsps[0].query_end==109
assert record.alignments[22].hsps[0].sbjct_start==106
assert record.alignments[22].hsps[0].sbjct_end==92
assert record.database_name==['data/sts']
assert record.posted_date==[('Feb 11, 2000  2:37 PM',)]
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==3835
assert record.num_sequences==87792
assert record.num_extends==3835
assert record.num_good_extends==1155
assert record.num_seqs_better_e==24
assert record.query_length==285
assert record.database_length==31998854
assert record.effective_hsp_length==17
assert record.effective_query_length==268
assert record.effective_database_length==30506390
assert record.effective_search_space==8175712520
assert record.effective_search_space_used==8175712520
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(10,19.8)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(15,30.2)

handle = open('Blast/bt049')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|1348400|gb|G26168|G26168 human STS\nEST47274.\x01gi|1375380|gb|G27130|G27130 human STS SHGC-32751."
assert record.query_letters==434
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==19
assert record.descriptions[0].title=="gi|1348400|gb|G26168|G26168 human STS EST47274. >gi|1375380|gb|..."
assert record.descriptions[0].score==718
assert record.descriptions[0].e==0.000000
assert record.descriptions[1].title=="gi|4632200|dbj|AU047565.1|AU047565 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[1].score==34
assert record.descriptions[1].e==0.650000
assert record.descriptions[2].title=="gi|6121436|gb|G56117.1|G56117 SHGC-101583 Human Homo sapiens ST..."
assert record.descriptions[2].score==34
assert record.descriptions[2].e==0.650000
assert record.descriptions[3].title=="gi|3249175|gb|G38401|G38401 SHGC-57345 Human Homo sapiens STS g..."
assert record.descriptions[3].score==34
assert record.descriptions[3].e==0.650000
assert record.descriptions[4].title=="gi|720383|gb|G03425|G03425 human STS WI-5868."
assert record.descriptions[4].score==34
assert record.descriptions[4].e==0.650000
assert record.descriptions[5].title=="gi|605557|gb|L31312|HUMUT937B Human STS UT937, 3' primer bind."
assert record.descriptions[5].score==34
assert record.descriptions[5].e==0.650000
assert record.descriptions[6].title=="gi|6123581|gb|G58262.1|G58262 SHGC-104352 Human Homo sapiens ST..."
assert record.descriptions[6].score==32
assert record.descriptions[6].e==2.600000
assert record.descriptions[7].title=="gi|6122805|gb|G57486.1|G57486 SHGC-103345 Human Homo sapiens ST..."
assert record.descriptions[7].score==32
assert record.descriptions[7].e==2.600000
assert record.descriptions[8].title=="gi|6121347|gb|G56178.1|G56178 SHGC-101470 Human Homo sapiens ST..."
assert record.descriptions[8].score==32
assert record.descriptions[8].e==2.600000
assert record.descriptions[9].title=="gi|3893806|emb|AL034156|HSPE59A01 H.sapiens flow-sorted chromos..."
assert record.descriptions[9].score==32
assert record.descriptions[9].e==2.600000
assert record.descriptions[10].title=="gi|5224295|gb|G52968.1|G52968 SHGC-86325 Human Homo sapiens STS..."
assert record.descriptions[10].score==32
assert record.descriptions[10].e==2.600000
assert record.descriptions[11].title=="gi|1348143|gb|G25911|G25911 human STS EST349382."
assert record.descriptions[11].score==32
assert record.descriptions[11].e==2.600000
assert record.descriptions[12].title=="gi|1233216|emb|Z51916|HSA082WB5 H.sapiens (D1S2890) DNA segment..."
assert record.descriptions[12].score==32
assert record.descriptions[12].e==2.600000
assert record.descriptions[13].title=="gi|1232198|emb|Z50898|HS038XD8 H.sapiens (D18S1106) DNA segment..."
assert record.descriptions[13].score==32
assert record.descriptions[13].e==2.600000
assert record.descriptions[14].title=="gi|1161890|gb|G16001|G16001 human STS CHLC.GCT8D08.P11278 clone..."
assert record.descriptions[14].score==32
assert record.descriptions[14].e==2.600000
assert record.descriptions[15].title=="gi|1130137|gb|G14398|G14398 human STS SHGC-9310 clone pG-5191."
assert record.descriptions[15].score==32
assert record.descriptions[15].e==2.600000
assert record.descriptions[16].title=="gi|1017612|gb|G11520|G11520 human STS SHGC-14676."
assert record.descriptions[16].score==32
assert record.descriptions[16].e==2.600000
assert record.descriptions[17].title=="gi|1396897|gb|G28178|G28178 human STS SHGC-34170."
assert record.descriptions[17].score==32
assert record.descriptions[17].e==2.600000
assert record.descriptions[18].title=="gi|859592|gb|G06347|G06347 human STS WI-7005."
assert record.descriptions[18].score==32
assert record.descriptions[18].e==2.600000
assert len(record.alignments)==0
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Feb 11, 2000  2:37 PM',)]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==8762
assert record.num_sequences==87792
assert record.num_extends==8762
assert record.num_good_extends==2655
assert record.num_seqs_better_e==27
assert record.query_length==434
assert record.database_length==31998854
assert record.effective_hsp_length==17
assert record.effective_query_length==417
assert record.effective_database_length==30506390
assert record.effective_search_space==12721164630
assert record.effective_search_space_used==12721164630
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(10,19.8)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(16,32.2)

handle = open('Blast/bt050')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|1347201|gb|G24969|G24969 human STS\nEST21946.\x01gi|1375315|gb|G27065|G27065 human STS SHGC-31731."
assert record.query_letters==331
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==0
assert len(record.alignments)==45
assert record.alignments[0].title==">gi|1347201|gb|G24969|G24969 human STS EST21946. >gi|1375315|gb|G27065|G27065 human STS SHGC-31731."
assert record.alignments[0].length==331
assert record.alignments[1].title==">gi|2665277|emb|AL010115|HSPE77H4 H.sapiens flow-sorted chromosome 1 HindIII fragment, SC1pE77H4, sequence tagged site [Homo sapiens]"
assert record.alignments[1].length==554
assert record.alignments[2].title==">gi|4757440|gb|G49267.1|G49267 stbK343C1_96809 chromosome 22 genomic clone Homo sapiens STS genomic clone 343C1, sequence tagged site"
assert record.alignments[2].length==360
assert record.alignments[3].title==">gi|4493602|gb|G47248.1|G47248 Z17392_1 Zebrafish AB Danio rerio STS genomic clone Z17392 5', sequence tagged site"
assert record.alignments[3].length==454
assert record.alignments[4].title==">gi|6120731|gb|G55412.1|G55412 SHGC-100745 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[4].length==652
assert record.alignments[5].title==">gi|1235411|emb|Z53965|HSC009WH1 H.sapiens (D2S2321) DNA segment containing (CA) repeat; clone AFMc009wh1; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[5].length==382
assert record.alignments[6].title==">gi|939357|gb|G08807|G08807 human STS CHLC.GATA70E11.P18111 clone GATA70E11"
assert record.alignments[6].length==643
assert record.alignments[7].title==">gi|1342139|gb|G21813|G21813 human STS WI-12408."
assert record.alignments[7].length==418
assert record.alignments[8].title==">gi|719782|gb|G02824|G02824 human STS WI-1312."
assert record.alignments[8].length==349
assert record.alignments[9].title==">gi|5902652|gb|G54536.1|G54536 Xq4072 KWOK Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[9].length==997
assert record.alignments[10].title==">gi|5566455|gb|AF167528.1|AF167528 Mus musculus chromosome 10 STS D10Jhu41, sequence tagged site"
assert record.alignments[10].length==749
assert record.alignments[11].title==">gi|5566418|gb|AF167504.1|AF167504 Mus musculus chromosome 10 STS D10Jhu59, sequence tagged site"
assert record.alignments[11].length==550
assert record.alignments[12].title==">gi|5565757|gb|AF096569.1|AF096569 Rattus norvegicus clone D5Uwm59, sequence tagged site"
assert record.alignments[12].length==580
assert record.alignments[13].title==">gi|4757436|gb|G49263.1|G49263 stbK343C1_70109 chromosome 22 genomic clone Homo sapiens STS genomic clone 343C1, sequence tagged site"
assert record.alignments[13].length==359
assert record.alignments[14].title==">gi|4518571|dbj|AU028648.1|AU028648 Rattus norvegicus, OTSUKA clone, OT10.41/02694, microsatellite sequence, sequence tagged site"
assert record.alignments[14].length==465
assert record.alignments[15].title==">gi|6123014|gb|G57845.1|G57845 SHGC-103599 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[15].length==680
assert record.alignments[16].title==">gi|6124406|gb|G59237.1|G59237 SHGC-110189 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[16].length==505
assert record.alignments[17].title==">gi|6121931|gb|G56612.1|G56612 SHGC-102181 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[17].length==489
assert record.alignments[18].title==">gi|3599905|gb|G41846|G41846 Z1061 Zebrafish AB Danio rerio STS genomic, sequence tagged site [Danio rerio]"
assert record.alignments[18].length==410
assert record.alignments[19].title==">gi|5224504|gb|G53327.1|G53327 SHGC-82321 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[19].length==421
assert record.alignments[20].title==">gi|5221389|gb|G50212.1|G50212 SHGC-82917 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[20].length==297
assert record.alignments[21].title==">gi|3168757|gb|G38183|G38183 RPCI-6-164E4Sp6 Human Homo sapiens STS genomic"
assert record.alignments[21].length==395
assert record.alignments[22].title==">gi|3150191|emb|AL023622|HS863J11T H.sapiens STS from genomic clone 863J11, sequence tagged site [Homo sapiens]"
assert record.alignments[22].length==558
assert record.alignments[23].title==">gi|3123402|emb|AL023351|DM171A7S Drosophila melanogaster STS determined from European Mapping Project cosmid, sequence tagged site [Drosophila melanogaster]"
assert record.alignments[23].length==600
assert record.alignments[24].title==">gi|2665031|emb|AL009868|HSPE36B05 H.sapiens flow-sorted chromosome 1 HindIII fragment, SC1pE36B05, sequence tagged site [Homo sapiens]"
assert record.alignments[24].length==397
assert record.alignments[25].title==">gi|2641969|dbj|AB004264|AB004264 Mus spretus genomic DNA, RLGS spot, D5Rik122, sequence tagged site [Mus spretus]"
assert record.alignments[25].length==237
assert record.alignments[26].title==">gi|4529138|gb|G48478.1|G48478 SHGC-68947 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[26].length==415
assert record.alignments[27].title==">gi|2484048|gb|G36284|G36284   STS h14a2375 5"
assert record.alignments[27].length==305
assert record.alignments[28].title==">gi|1871265|gb|G31295|G31295 sy625c11-R Human (A.Gnirke) Homo sapiens STS genomic, sequence tagged site [Homo sapiens]"
assert record.alignments[28].length==279
assert record.alignments[29].title==">gi|1263801|emb|Z70926|DM63D12T D. melanogaster STS determined from European Mapping Project cosmid, sequence tagged site [Drosophila melanogaster]"
assert record.alignments[29].length==350
assert record.alignments[30].title==">gi|1244248|gb|G19461|G19461 human STS SHGC-11979."
assert record.alignments[30].length==400
assert record.alignments[31].title==">gi|1233253|emb|Z51953|HSA084XF1 H.sapiens (D19S932) DNA segment containing (CA) repeat; clone AFMa084xf1; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[31].length==276
assert record.alignments[32].title==">gi|1233152|emb|Z51852|HSA070WH1 H.sapiens (D9S1879) DNA segment containing (CA) repeat; clone AFMa070wh1; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[32].length==363
assert record.alignments[33].title==">gi|1215106|gb|G17680|G17680 human STS SHGC-3112 clone pG-956."
assert record.alignments[33].length==313
assert record.alignments[34].title==">gi|1594113|gb|G30562|G30562 human STS SHGC-37420"
assert record.alignments[34].length==447
assert record.alignments[35].title==">gi|1593022|gb|G29471|G29471 human STS SHGC-33749"
assert record.alignments[35].length==437
assert record.alignments[36].title==">gi|1052040|emb|Z67257|HSA090WF9 H.sapiens DNA segment containing (CA) repeat; clone AFMa090wf9; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[36].length==312
assert record.alignments[37].title==">gi|881858|gb|G07651|G07651 human STS SHGC-5830 clone pG-2175."
assert record.alignments[37].length==336
assert record.alignments[38].title==">gi|1344417|gb|G24091|G24091 human STS WI-12544."
assert record.alignments[38].length==452
assert record.alignments[39].title==">gi|860391|gb|G07146|G07146 human STS WI-9143."
assert record.alignments[39].length==1946
assert record.alignments[40].title==">gi|858959|gb|G05714|G05714 human STS WI-8961."
assert record.alignments[40].length==1337
assert record.alignments[41].title==">gi|1341506|gb|G21180|G21180 human STS WI-14638."
assert record.alignments[41].length==375
assert record.alignments[42].title==">gi|1341551|gb|G21225|G21225 human STS WI-11777."
assert record.alignments[42].length==348
assert record.alignments[43].title==">gi|465219|gb|L29850|HUMUT50 Human STS UT50."
assert record.alignments[43].length==382
assert record.alignments[44].title==">gi|860599|gb|G07354|G07354 human STS WI-9600."
assert record.alignments[44].length==226
assert record.alignments[0].hsps[0].score==331
assert record.alignments[0].hsps[0].bits==656
assert record.alignments[0].hsps[0].expect==0.0
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==17
assert record.alignments[1].hsps[0].bits==34.2
assert record.alignments[1].hsps[0].expect==0.49
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==16
assert record.alignments[2].hsps[0].bits==32.2
assert record.alignments[2].hsps[0].expect==1.9
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==16
assert record.alignments[3].hsps[0].bits==32.2
assert record.alignments[3].hsps[0].expect==1.9
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==16
assert record.alignments[4].hsps[0].bits==32.2
assert record.alignments[4].hsps[0].expect==1.9
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==16
assert record.alignments[5].hsps[0].bits==32.2
assert record.alignments[5].hsps[0].expect==1.9
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==16
assert record.alignments[6].hsps[0].bits==32.2
assert record.alignments[6].hsps[0].expect==1.9
assert len(record.alignments[6].hsps)==1
assert record.alignments[7].hsps[0].score==16
assert record.alignments[7].hsps[0].bits==32.2
assert record.alignments[7].hsps[0].expect==1.9
assert len(record.alignments[7].hsps)==1
assert record.alignments[8].hsps[0].score==16
assert record.alignments[8].hsps[0].bits==32.2
assert record.alignments[8].hsps[0].expect==1.9
assert len(record.alignments[8].hsps)==1
assert record.alignments[9].hsps[0].score==15
assert record.alignments[9].hsps[0].bits==30.2
assert record.alignments[9].hsps[0].expect==7.6
assert len(record.alignments[9].hsps)==1
assert record.alignments[10].hsps[0].score==15
assert record.alignments[10].hsps[0].bits==30.2
assert record.alignments[10].hsps[0].expect==7.6
assert len(record.alignments[10].hsps)==1
assert record.alignments[11].hsps[0].score==15
assert record.alignments[11].hsps[0].bits==30.2
assert record.alignments[11].hsps[0].expect==7.6
assert len(record.alignments[11].hsps)==1
assert record.alignments[12].hsps[0].score==15
assert record.alignments[12].hsps[0].bits==30.2
assert record.alignments[12].hsps[0].expect==7.6
assert len(record.alignments[12].hsps)==1
assert record.alignments[13].hsps[0].score==15
assert record.alignments[13].hsps[0].bits==30.2
assert record.alignments[13].hsps[0].expect==7.6
assert len(record.alignments[13].hsps)==1
assert record.alignments[14].hsps[0].score==15
assert record.alignments[14].hsps[0].bits==30.2
assert record.alignments[14].hsps[0].expect==7.6
assert len(record.alignments[14].hsps)==1
assert record.alignments[15].hsps[0].score==15
assert record.alignments[15].hsps[0].bits==30.2
assert record.alignments[15].hsps[0].expect==7.6
assert len(record.alignments[15].hsps)==1
assert record.alignments[16].hsps[0].score==15
assert record.alignments[16].hsps[0].bits==30.2
assert record.alignments[16].hsps[0].expect==7.6
assert len(record.alignments[16].hsps)==1
assert record.alignments[17].hsps[0].score==15
assert record.alignments[17].hsps[0].bits==30.2
assert record.alignments[17].hsps[0].expect==7.6
assert len(record.alignments[17].hsps)==1
assert record.alignments[18].hsps[0].score==15
assert record.alignments[18].hsps[0].bits==30.2
assert record.alignments[18].hsps[0].expect==7.6
assert len(record.alignments[18].hsps)==1
assert record.alignments[19].hsps[0].score==15
assert record.alignments[19].hsps[0].bits==30.2
assert record.alignments[19].hsps[0].expect==7.6
assert len(record.alignments[19].hsps)==1
assert record.alignments[20].hsps[0].score==15
assert record.alignments[20].hsps[0].bits==30.2
assert record.alignments[20].hsps[0].expect==7.6
assert len(record.alignments[20].hsps)==1
assert record.alignments[21].hsps[0].score==15
assert record.alignments[21].hsps[0].bits==30.2
assert record.alignments[21].hsps[0].expect==7.6
assert len(record.alignments[21].hsps)==1
assert record.alignments[22].hsps[0].score==15
assert record.alignments[22].hsps[0].bits==30.2
assert record.alignments[22].hsps[0].expect==7.6
assert len(record.alignments[22].hsps)==1
assert record.alignments[23].hsps[0].score==15
assert record.alignments[23].hsps[0].bits==30.2
assert record.alignments[23].hsps[0].expect==7.6
assert len(record.alignments[23].hsps)==1
assert record.alignments[24].hsps[0].score==15
assert record.alignments[24].hsps[0].bits==30.2
assert record.alignments[24].hsps[0].expect==7.6
assert len(record.alignments[24].hsps)==1
assert record.alignments[25].hsps[0].score==15
assert record.alignments[25].hsps[0].bits==30.2
assert record.alignments[25].hsps[0].expect==7.6
assert len(record.alignments[25].hsps)==1
assert record.alignments[26].hsps[0].score==15
assert record.alignments[26].hsps[0].bits==30.2
assert record.alignments[26].hsps[0].expect==7.6
assert len(record.alignments[26].hsps)==1
assert record.alignments[27].hsps[0].score==15
assert record.alignments[27].hsps[0].bits==30.2
assert record.alignments[27].hsps[0].expect==7.6
assert len(record.alignments[27].hsps)==1
assert record.alignments[28].hsps[0].score==15
assert record.alignments[28].hsps[0].bits==30.2
assert record.alignments[28].hsps[0].expect==7.6
assert len(record.alignments[28].hsps)==1
assert record.alignments[29].hsps[0].score==15
assert record.alignments[29].hsps[0].bits==30.2
assert record.alignments[29].hsps[0].expect==7.6
assert len(record.alignments[29].hsps)==1
assert record.alignments[30].hsps[0].score==15
assert record.alignments[30].hsps[0].bits==30.2
assert record.alignments[30].hsps[0].expect==7.6
assert len(record.alignments[30].hsps)==1
assert record.alignments[31].hsps[0].score==15
assert record.alignments[31].hsps[0].bits==30.2
assert record.alignments[31].hsps[0].expect==7.6
assert len(record.alignments[31].hsps)==1
assert record.alignments[32].hsps[0].score==15
assert record.alignments[32].hsps[0].bits==30.2
assert record.alignments[32].hsps[0].expect==7.6
assert len(record.alignments[32].hsps)==1
assert record.alignments[33].hsps[0].score==15
assert record.alignments[33].hsps[0].bits==30.2
assert record.alignments[33].hsps[0].expect==7.6
assert len(record.alignments[33].hsps)==1
assert record.alignments[34].hsps[0].score==15
assert record.alignments[34].hsps[0].bits==30.2
assert record.alignments[34].hsps[0].expect==7.6
assert len(record.alignments[34].hsps)==1
assert record.alignments[35].hsps[0].score==15
assert record.alignments[35].hsps[0].bits==30.2
assert record.alignments[35].hsps[0].expect==7.6
assert len(record.alignments[35].hsps)==1
assert record.alignments[36].hsps[0].score==15
assert record.alignments[36].hsps[0].bits==30.2
assert record.alignments[36].hsps[0].expect==7.6
assert len(record.alignments[36].hsps)==1
assert record.alignments[37].hsps[0].score==15
assert record.alignments[37].hsps[0].bits==30.2
assert record.alignments[37].hsps[0].expect==7.6
assert len(record.alignments[37].hsps)==1
assert record.alignments[38].hsps[0].score==15
assert record.alignments[38].hsps[0].bits==30.2
assert record.alignments[38].hsps[0].expect==7.6
assert len(record.alignments[38].hsps)==1
assert record.alignments[39].hsps[0].score==15
assert record.alignments[39].hsps[0].bits==30.2
assert record.alignments[39].hsps[0].expect==7.6
assert len(record.alignments[39].hsps)==1
assert record.alignments[40].hsps[0].score==15
assert record.alignments[40].hsps[0].bits==30.2
assert record.alignments[40].hsps[0].expect==7.6
assert len(record.alignments[40].hsps)==1
assert record.alignments[41].hsps[0].score==15
assert record.alignments[41].hsps[0].bits==30.2
assert record.alignments[41].hsps[0].expect==7.6
assert len(record.alignments[41].hsps)==1
assert record.alignments[42].hsps[0].score==15
assert record.alignments[42].hsps[0].bits==30.2
assert record.alignments[42].hsps[0].expect==7.6
assert len(record.alignments[42].hsps)==1
assert record.alignments[43].hsps[0].score==15
assert record.alignments[43].hsps[0].bits==30.2
assert record.alignments[43].hsps[0].expect==7.6
assert len(record.alignments[43].hsps)==1
assert record.alignments[44].hsps[0].score==15
assert record.alignments[44].hsps[0].bits==30.2
assert record.alignments[44].hsps[0].expect==7.6
assert len(record.alignments[44].hsps)==1
assert record.alignments[0].hsps[0].identities==(331, 331)
assert record.alignments[1].hsps[0].identities==(17, 17)
assert record.alignments[2].hsps[0].identities==(22, 24)
assert record.alignments[3].hsps[0].identities==(16, 16)
assert record.alignments[4].hsps[0].identities==(16, 16)
assert record.alignments[5].hsps[0].identities==(16, 16)
assert record.alignments[6].hsps[0].identities==(16, 16)
assert record.alignments[7].hsps[0].identities==(16, 16)
assert record.alignments[8].hsps[0].identities==(16, 16)
assert record.alignments[9].hsps[0].identities==(15, 15)
assert record.alignments[10].hsps[0].identities==(15, 15)
assert record.alignments[11].hsps[0].identities==(15, 15)
assert record.alignments[12].hsps[0].identities==(21, 23)
assert record.alignments[13].hsps[0].identities==(15, 15)
assert record.alignments[14].hsps[0].identities==(15, 15)
assert record.alignments[15].hsps[0].identities==(15, 15)
assert record.alignments[16].hsps[0].identities==(15, 15)
assert record.alignments[17].hsps[0].identities==(18, 19)
assert record.alignments[18].hsps[0].identities==(15, 15)
assert record.alignments[19].hsps[0].identities==(18, 19)
assert record.alignments[20].hsps[0].identities==(15, 15)
assert record.alignments[21].hsps[0].identities==(15, 15)
assert record.alignments[22].hsps[0].identities==(15, 15)
assert record.alignments[23].hsps[0].identities==(15, 15)
assert record.alignments[24].hsps[0].identities==(15, 15)
assert record.alignments[25].hsps[0].identities==(15, 15)
assert record.alignments[26].hsps[0].identities==(15, 15)
assert record.alignments[27].hsps[0].identities==(15, 15)
assert record.alignments[28].hsps[0].identities==(15, 15)
assert record.alignments[29].hsps[0].identities==(15, 15)
assert record.alignments[30].hsps[0].identities==(15, 15)
assert record.alignments[31].hsps[0].identities==(15, 15)
assert record.alignments[32].hsps[0].identities==(15, 15)
assert record.alignments[33].hsps[0].identities==(18, 19)
assert record.alignments[34].hsps[0].identities==(18, 19)
assert record.alignments[35].hsps[0].identities==(15, 15)
assert record.alignments[36].hsps[0].identities==(15, 15)
assert record.alignments[37].hsps[0].identities==(15, 15)
assert record.alignments[38].hsps[0].identities==(18, 19)
assert record.alignments[39].hsps[0].identities==(15, 15)
assert record.alignments[40].hsps[0].identities==(15, 15)
assert record.alignments[41].hsps[0].identities==(15, 15)
assert record.alignments[42].hsps[0].identities==(15, 15)
assert record.alignments[43].hsps[0].identities==(15, 15)
assert record.alignments[44].hsps[0].identities==(15, 15)
assert record.alignments[0].hsps[0].strand==("Plus", "Plus")
assert record.alignments[1].hsps[0].strand==("Plus", "Minus")
assert record.alignments[2].hsps[0].strand==("Plus", "Minus")
assert record.alignments[3].hsps[0].strand==("Plus", "Plus")
assert record.alignments[4].hsps[0].strand==("Plus", "Minus")
assert record.alignments[5].hsps[0].strand==("Plus", "Plus")
assert record.alignments[6].hsps[0].strand==("Plus", "Minus")
assert record.alignments[7].hsps[0].strand==("Plus", "Minus")
assert record.alignments[8].hsps[0].strand==("Plus", "Plus")
assert record.alignments[9].hsps[0].strand==("Plus", "Minus")
assert record.alignments[10].hsps[0].strand==("Plus", "Plus")
assert record.alignments[11].hsps[0].strand==("Plus", "Minus")
assert record.alignments[12].hsps[0].strand==("Plus", "Plus")
assert record.alignments[13].hsps[0].strand==("Plus", "Minus")
assert record.alignments[14].hsps[0].strand==("Plus", "Minus")
assert record.alignments[15].hsps[0].strand==("Plus", "Plus")
assert record.alignments[16].hsps[0].strand==("Plus", "Minus")
assert record.alignments[17].hsps[0].strand==("Plus", "Plus")
assert record.alignments[18].hsps[0].strand==("Plus", "Plus")
assert record.alignments[19].hsps[0].strand==("Plus", "Plus")
assert record.alignments[20].hsps[0].strand==("Plus", "Plus")
assert record.alignments[21].hsps[0].strand==("Plus", "Plus")
assert record.alignments[22].hsps[0].strand==("Plus", "Plus")
assert record.alignments[23].hsps[0].strand==("Plus", "Plus")
assert record.alignments[24].hsps[0].strand==("Plus", "Plus")
assert record.alignments[25].hsps[0].strand==("Plus", "Minus")
assert record.alignments[26].hsps[0].strand==("Plus", "Plus")
assert record.alignments[27].hsps[0].strand==("Plus", "Plus")
assert record.alignments[28].hsps[0].strand==("Plus", "Plus")
assert record.alignments[29].hsps[0].strand==("Plus", "Minus")
assert record.alignments[30].hsps[0].strand==("Plus", "Minus")
assert record.alignments[31].hsps[0].strand==("Plus", "Minus")
assert record.alignments[32].hsps[0].strand==("Plus", "Plus")
assert record.alignments[33].hsps[0].strand==("Plus", "Minus")
assert record.alignments[34].hsps[0].strand==("Plus", "Plus")
assert record.alignments[35].hsps[0].strand==("Plus", "Minus")
assert record.alignments[36].hsps[0].strand==("Plus", "Plus")
assert record.alignments[37].hsps[0].strand==("Plus", "Plus")
assert record.alignments[38].hsps[0].strand==("Plus", "Minus")
assert record.alignments[39].hsps[0].strand==("Plus", "Minus")
assert record.alignments[40].hsps[0].strand==("Plus", "Plus")
assert record.alignments[41].hsps[0].strand==("Plus", "Minus")
assert record.alignments[42].hsps[0].strand==("Plus", "Plus")
assert record.alignments[43].hsps[0].strand==("Plus", "Minus")
assert record.alignments[44].hsps[0].strand==("Plus", "Plus")
assert record.alignments[0].hsps[0].query=="cctccaccctctcatgagcaacaggatatgtgaaagtacttgcagccagaagcaaaaccacaatcctcgggtgctagatggagctccccaaggagcagagaggaaaaggcaggaggagagggccaggcagcagggatggagactaagtttggcccaaggctgcccgcaagcactgatgccatcatgccctctggtaggtgtctatttctgtctgaaccagaaatacaccaagctccacacatgggggctttgctggcttcgacatcactggttcaactatgtcactgctttgttatatttagtgctccagaacctcaggttccttcagatt"
assert record.alignments[0].hsps[0].match=="|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
assert record.alignments[0].hsps[0].sbjct=="cctccaccctctcatgagcaacaggatatgtgaaagtacttgcagccagaagcaaaaccacaatcctcgggtgctagatggagctccccaaggagcagagaggaaaaggcaggaggagagggccaggcagcagggatggagactaagtttggcccaaggctgcccgcaagcactgatgccatcatgccctctggtaggtgtctatttctgtctgaaccagaaatacaccaagctccacacatgggggctttgctggcttcgacatcactggttcaactatgtcactgctttgttatatttagtgctccagaacctcaggttccttcagatt"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==331
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==331
assert record.alignments[1].hsps[0].query=="ccaggcagcagggatgg"
assert record.alignments[1].hsps[0].match=="|||||||||||||||||"
assert record.alignments[1].hsps[0].sbjct=="ccaggcagcagggatgg"
assert record.alignments[1].hsps[0].query_start==123
assert record.alignments[1].hsps[0].query_end==139
assert record.alignments[1].hsps[0].sbjct_start==434
assert record.alignments[1].hsps[0].sbjct_end==418
assert record.alignments[2].hsps[0].query=="ggaggagagggccaggcagcaggg"
assert record.alignments[2].hsps[0].match=="||||||||||| | ||||||||||"
assert record.alignments[2].hsps[0].sbjct=="ggaggagaggggctggcagcaggg"
assert record.alignments[2].hsps[0].query_start==112
assert record.alignments[2].hsps[0].query_end==135
assert record.alignments[2].hsps[0].sbjct_start==287
assert record.alignments[2].hsps[0].sbjct_end==264
assert record.alignments[3].hsps[0].query=="agaagcaaaaccacaa"
assert record.alignments[3].hsps[0].match=="||||||||||||||||"
assert record.alignments[3].hsps[0].sbjct=="agaagcaaaaccacaa"
assert record.alignments[3].hsps[0].query_start==48
assert record.alignments[3].hsps[0].query_end==63
assert record.alignments[3].hsps[0].sbjct_start==434
assert record.alignments[3].hsps[0].sbjct_end==449
assert record.alignments[4].hsps[0].query=="agaggaaaaggcagga"
assert record.alignments[4].hsps[0].match=="||||||||||||||||"
assert record.alignments[4].hsps[0].sbjct=="agaggaaaaggcagga"
assert record.alignments[4].hsps[0].query_start==99
assert record.alignments[4].hsps[0].query_end==114
assert record.alignments[4].hsps[0].sbjct_start==431
assert record.alignments[4].hsps[0].sbjct_end==416
assert record.alignments[5].hsps[0].query=="agagaggaaaaggcag"
assert record.alignments[5].hsps[0].match=="||||||||||||||||"
assert record.alignments[5].hsps[0].sbjct=="agagaggaaaaggcag"
assert record.alignments[5].hsps[0].query_start==97
assert record.alignments[5].hsps[0].query_end==112
assert record.alignments[5].hsps[0].sbjct_start==107
assert record.alignments[5].hsps[0].sbjct_end==122
assert record.alignments[6].hsps[0].query=="gaggaaaaggcaggag"
assert record.alignments[6].hsps[0].match=="||||||||||||||||"
assert record.alignments[6].hsps[0].sbjct=="gaggaaaaggcaggag"
assert record.alignments[6].hsps[0].query_start==100
assert record.alignments[6].hsps[0].query_end==115
assert record.alignments[6].hsps[0].sbjct_start==482
assert record.alignments[6].hsps[0].sbjct_end==467
assert record.alignments[7].hsps[0].query=="cagaagcaaaaccaca"
assert record.alignments[7].hsps[0].match=="||||||||||||||||"
assert record.alignments[7].hsps[0].sbjct=="cagaagcaaaaccaca"
assert record.alignments[7].hsps[0].query_start==47
assert record.alignments[7].hsps[0].query_end==62
assert record.alignments[7].hsps[0].sbjct_start==193
assert record.alignments[7].hsps[0].sbjct_end==178
assert record.alignments[8].hsps[0].query=="tgctttgttatattta"
assert record.alignments[8].hsps[0].match=="||||||||||||||||"
assert record.alignments[8].hsps[0].sbjct=="tgctttgttatattta"
assert record.alignments[8].hsps[0].query_start==286
assert record.alignments[8].hsps[0].query_end==301
assert record.alignments[8].hsps[0].sbjct_start==111
assert record.alignments[8].hsps[0].sbjct_end==126
assert record.alignments[9].hsps[0].query=="cctccaccctctcat"
assert record.alignments[9].hsps[0].match=="|||||||||||||||"
assert record.alignments[9].hsps[0].sbjct=="cctccaccctctcat"
assert record.alignments[9].hsps[0].query_start==1
assert record.alignments[9].hsps[0].query_end==15
assert record.alignments[9].hsps[0].sbjct_start==618
assert record.alignments[9].hsps[0].sbjct_end==604
assert record.alignments[10].hsps[0].query=="agagggccaggcagc"
assert record.alignments[10].hsps[0].match=="|||||||||||||||"
assert record.alignments[10].hsps[0].sbjct=="agagggccaggcagc"
assert record.alignments[10].hsps[0].query_start==117
assert record.alignments[10].hsps[0].query_end==131
assert record.alignments[10].hsps[0].sbjct_start==487
assert record.alignments[10].hsps[0].sbjct_end==501
assert record.alignments[11].hsps[0].query=="gagctccccaaggag"
assert record.alignments[11].hsps[0].match=="|||||||||||||||"
assert record.alignments[11].hsps[0].sbjct=="gagctccccaaggag"
assert record.alignments[11].hsps[0].query_start==81
assert record.alignments[11].hsps[0].query_end==95
assert record.alignments[11].hsps[0].sbjct_start==278
assert record.alignments[11].hsps[0].sbjct_end==264
assert record.alignments[12].hsps[0].query=="gcagagaggaaaaggcaggagga"
assert record.alignments[12].hsps[0].match=="||||| |||| ||||||||||||"
assert record.alignments[12].hsps[0].sbjct=="gcagacaggagaaggcaggagga"
assert record.alignments[12].hsps[0].query_start==95
assert record.alignments[12].hsps[0].query_end==117
assert record.alignments[12].hsps[0].sbjct_start==84
assert record.alignments[12].hsps[0].sbjct_end==106
assert record.alignments[13].hsps[0].query=="aggcagcagggatgg"
assert record.alignments[13].hsps[0].match=="|||||||||||||||"
assert record.alignments[13].hsps[0].sbjct=="aggcagcagggatgg"
assert record.alignments[13].hsps[0].query_start==125
assert record.alignments[13].hsps[0].query_end==139
assert record.alignments[13].hsps[0].sbjct_start==58
assert record.alignments[13].hsps[0].sbjct_end==44
assert record.alignments[14].hsps[0].query=="atgccctctggtagg"
assert record.alignments[14].hsps[0].match=="|||||||||||||||"
assert record.alignments[14].hsps[0].sbjct=="atgccctctggtagg"
assert record.alignments[14].hsps[0].query_start==184
assert record.alignments[14].hsps[0].query_end==198
assert record.alignments[14].hsps[0].sbjct_start==114
assert record.alignments[14].hsps[0].sbjct_end==100
assert record.alignments[15].hsps[0].query=="ggatggagactaagt"
assert record.alignments[15].hsps[0].match=="|||||||||||||||"
assert record.alignments[15].hsps[0].sbjct=="ggatggagactaagt"
assert record.alignments[15].hsps[0].query_start==134
assert record.alignments[15].hsps[0].query_end==148
assert record.alignments[15].hsps[0].sbjct_start==451
assert record.alignments[15].hsps[0].sbjct_end==465
assert record.alignments[16].hsps[0].query=="cagcagggatggaga"
assert record.alignments[16].hsps[0].match=="|||||||||||||||"
assert record.alignments[16].hsps[0].sbjct=="cagcagggatggaga"
assert record.alignments[16].hsps[0].query_start==128
assert record.alignments[16].hsps[0].query_end==142
assert record.alignments[16].hsps[0].sbjct_start==163
assert record.alignments[16].hsps[0].sbjct_end==149
assert record.alignments[17].hsps[0].query=="tggagactaagtttggccc"
assert record.alignments[17].hsps[0].match=="||||||||||||| |||||"
assert record.alignments[17].hsps[0].sbjct=="tggagactaagttgggccc"
assert record.alignments[17].hsps[0].query_start==137
assert record.alignments[17].hsps[0].query_end==155
assert record.alignments[17].hsps[0].sbjct_start==302
assert record.alignments[17].hsps[0].sbjct_end==320
assert record.alignments[18].hsps[0].query=="tcactggttcaacta"
assert record.alignments[18].hsps[0].match=="|||||||||||||||"
assert record.alignments[18].hsps[0].sbjct=="tcactggttcaacta"
assert record.alignments[18].hsps[0].query_start==265
assert record.alignments[18].hsps[0].query_end==279
assert record.alignments[18].hsps[0].sbjct_start==185
assert record.alignments[18].hsps[0].sbjct_end==199
assert record.alignments[19].hsps[0].query=="aaggagcagagaggaaaag"
assert record.alignments[19].hsps[0].match=="|||||| ||||||||||||"
assert record.alignments[19].hsps[0].sbjct=="aaggagaagagaggaaaag"
assert record.alignments[19].hsps[0].query_start==90
assert record.alignments[19].hsps[0].query_end==108
assert record.alignments[19].hsps[0].sbjct_start==209
assert record.alignments[19].hsps[0].sbjct_end==227
assert record.alignments[20].hsps[0].query=="gagagggccaggcag"
assert record.alignments[20].hsps[0].match=="|||||||||||||||"
assert record.alignments[20].hsps[0].sbjct=="gagagggccaggcag"
assert record.alignments[20].hsps[0].query_start==116
assert record.alignments[20].hsps[0].query_end==130
assert record.alignments[20].hsps[0].sbjct_start==155
assert record.alignments[20].hsps[0].sbjct_end==169
assert record.alignments[21].hsps[0].query=="gcagccagaagcaaa"
assert record.alignments[21].hsps[0].match=="|||||||||||||||"
assert record.alignments[21].hsps[0].sbjct=="gcagccagaagcaaa"
assert record.alignments[21].hsps[0].query_start==42
assert record.alignments[21].hsps[0].query_end==56
assert record.alignments[21].hsps[0].sbjct_start==259
assert record.alignments[21].hsps[0].sbjct_end==273
assert record.alignments[22].hsps[0].query=="tgtctatttctgtct"
assert record.alignments[22].hsps[0].match=="|||||||||||||||"
assert record.alignments[22].hsps[0].sbjct=="tgtctatttctgtct"
assert record.alignments[22].hsps[0].query_start==199
assert record.alignments[22].hsps[0].query_end==213
assert record.alignments[22].hsps[0].sbjct_start==155
assert record.alignments[22].hsps[0].sbjct_end==169
assert record.alignments[23].hsps[0].query=="acaccaagctccaca"
assert record.alignments[23].hsps[0].match=="|||||||||||||||"
assert record.alignments[23].hsps[0].sbjct=="acaccaagctccaca"
assert record.alignments[23].hsps[0].query_start==225
assert record.alignments[23].hsps[0].query_end==239
assert record.alignments[23].hsps[0].sbjct_start==163
assert record.alignments[23].hsps[0].sbjct_end==177
assert record.alignments[24].hsps[0].query=="tgtctatttctgtct"
assert record.alignments[24].hsps[0].match=="|||||||||||||||"
assert record.alignments[24].hsps[0].sbjct=="tgtctatttctgtct"
assert record.alignments[24].hsps[0].query_start==199
assert record.alignments[24].hsps[0].query_end==213
assert record.alignments[24].hsps[0].sbjct_start==244
assert record.alignments[24].hsps[0].sbjct_end==258
assert record.alignments[25].hsps[0].query=="ggagagggccaggca"
assert record.alignments[25].hsps[0].match=="|||||||||||||||"
assert record.alignments[25].hsps[0].sbjct=="ggagagggccaggca"
assert record.alignments[25].hsps[0].query_start==115
assert record.alignments[25].hsps[0].query_end==129
assert record.alignments[25].hsps[0].sbjct_start==53
assert record.alignments[25].hsps[0].sbjct_end==39
assert record.alignments[26].hsps[0].query=="acatcactggttcaa"
assert record.alignments[26].hsps[0].match=="|||||||||||||||"
assert record.alignments[26].hsps[0].sbjct=="acatcactggttcaa"
assert record.alignments[26].hsps[0].query_start==262
assert record.alignments[26].hsps[0].query_end==276
assert record.alignments[26].hsps[0].sbjct_start==331
assert record.alignments[26].hsps[0].sbjct_end==345
assert record.alignments[27].hsps[0].query=="ggccaggcagcaggg"
assert record.alignments[27].hsps[0].match=="|||||||||||||||"
assert record.alignments[27].hsps[0].sbjct=="ggccaggcagcaggg"
assert record.alignments[27].hsps[0].query_start==121
assert record.alignments[27].hsps[0].query_end==135
assert record.alignments[27].hsps[0].sbjct_start==108
assert record.alignments[27].hsps[0].sbjct_end==122
assert record.alignments[28].hsps[0].query=="tttagtgctccagaa"
assert record.alignments[28].hsps[0].match=="|||||||||||||||"
assert record.alignments[28].hsps[0].sbjct=="tttagtgctccagaa"
assert record.alignments[28].hsps[0].query_start==298
assert record.alignments[28].hsps[0].query_end==312
assert record.alignments[28].hsps[0].sbjct_start==261
assert record.alignments[28].hsps[0].sbjct_end==275
assert record.alignments[29].hsps[0].query=="caggttccttcagat"
assert record.alignments[29].hsps[0].match=="|||||||||||||||"
assert record.alignments[29].hsps[0].sbjct=="caggttccttcagat"
assert record.alignments[29].hsps[0].query_start==316
assert record.alignments[29].hsps[0].query_end==330
assert record.alignments[29].hsps[0].sbjct_start==323
assert record.alignments[29].hsps[0].sbjct_end==309
assert record.alignments[30].hsps[0].query=="cagccagaagcaaaa"
assert record.alignments[30].hsps[0].match=="|||||||||||||||"
assert record.alignments[30].hsps[0].sbjct=="cagccagaagcaaaa"
assert record.alignments[30].hsps[0].query_start==43
assert record.alignments[30].hsps[0].query_end==57
assert record.alignments[30].hsps[0].sbjct_start==356
assert record.alignments[30].hsps[0].sbjct_end==342
assert record.alignments[31].hsps[0].query=="aacctcaggttcctt"
assert record.alignments[31].hsps[0].match=="|||||||||||||||"
assert record.alignments[31].hsps[0].sbjct=="aacctcaggttcctt"
assert record.alignments[31].hsps[0].query_start==311
assert record.alignments[31].hsps[0].query_end==325
assert record.alignments[31].hsps[0].sbjct_start==54
assert record.alignments[31].hsps[0].sbjct_end==40
assert record.alignments[32].hsps[0].query=="cctccaccctctcat"
assert record.alignments[32].hsps[0].match=="|||||||||||||||"
assert record.alignments[32].hsps[0].sbjct=="cctccaccctctcat"
assert record.alignments[32].hsps[0].query_start==1
assert record.alignments[32].hsps[0].query_end==15
assert record.alignments[32].hsps[0].sbjct_start==59
assert record.alignments[32].hsps[0].sbjct_end==73
assert record.alignments[33].hsps[0].query=="ccccaaggagcagagagga"
assert record.alignments[33].hsps[0].match=="|||| ||||||||||||||"
assert record.alignments[33].hsps[0].sbjct=="ccccgaggagcagagagga"
assert record.alignments[33].hsps[0].query_start==86
assert record.alignments[33].hsps[0].query_end==104
assert record.alignments[33].hsps[0].sbjct_start==213
assert record.alignments[33].hsps[0].sbjct_end==195
assert record.alignments[34].hsps[0].query=="ggaggagagggccaggcag"
assert record.alignments[34].hsps[0].match=="||||||| |||||||||||"
assert record.alignments[34].hsps[0].sbjct=="ggaggagggggccaggcag"
assert record.alignments[34].hsps[0].query_start==112
assert record.alignments[34].hsps[0].query_end==130
assert record.alignments[34].hsps[0].sbjct_start==234
assert record.alignments[34].hsps[0].sbjct_end==252
assert record.alignments[35].hsps[0].query=="ccaggcagcagggat"
assert record.alignments[35].hsps[0].match=="|||||||||||||||"
assert record.alignments[35].hsps[0].sbjct=="ccaggcagcagggat"
assert record.alignments[35].hsps[0].query_start==123
assert record.alignments[35].hsps[0].query_end==137
assert record.alignments[35].hsps[0].sbjct_start==320
assert record.alignments[35].hsps[0].sbjct_end==306
assert record.alignments[36].hsps[0].query=="gtaggtgtctatttc"
assert record.alignments[36].hsps[0].match=="|||||||||||||||"
assert record.alignments[36].hsps[0].sbjct=="gtaggtgtctatttc"
assert record.alignments[36].hsps[0].query_start==194
assert record.alignments[36].hsps[0].query_end==208
assert record.alignments[36].hsps[0].sbjct_start==201
assert record.alignments[36].hsps[0].sbjct_end==215
assert record.alignments[37].hsps[0].query=="ggccaggcagcaggg"
assert record.alignments[37].hsps[0].match=="|||||||||||||||"
assert record.alignments[37].hsps[0].sbjct=="ggccaggcagcaggg"
assert record.alignments[37].hsps[0].query_start==121
assert record.alignments[37].hsps[0].query_end==135
assert record.alignments[37].hsps[0].sbjct_start==86
assert record.alignments[37].hsps[0].sbjct_end==100
assert record.alignments[38].hsps[0].query=="aaaaggcaggaggagaggg"
assert record.alignments[38].hsps[0].match=="|||||||||||| ||||||"
assert record.alignments[38].hsps[0].sbjct=="aaaaggcaggagaagaggg"
assert record.alignments[38].hsps[0].query_start==104
assert record.alignments[38].hsps[0].query_end==122
assert record.alignments[38].hsps[0].sbjct_start==234
assert record.alignments[38].hsps[0].sbjct_end==216
assert record.alignments[39].hsps[0].query=="tgtcactgctttgtt"
assert record.alignments[39].hsps[0].match=="|||||||||||||||"
assert record.alignments[39].hsps[0].sbjct=="tgtcactgctttgtt"
assert record.alignments[39].hsps[0].query_start==280
assert record.alignments[39].hsps[0].query_end==294
assert record.alignments[39].hsps[0].sbjct_start==492
assert record.alignments[39].hsps[0].sbjct_end==478
assert record.alignments[40].hsps[0].query=="acacatgggggcttt"
assert record.alignments[40].hsps[0].match=="|||||||||||||||"
assert record.alignments[40].hsps[0].sbjct=="acacatgggggcttt"
assert record.alignments[40].hsps[0].query_start==237
assert record.alignments[40].hsps[0].query_end==251
assert record.alignments[40].hsps[0].sbjct_start==117
assert record.alignments[40].hsps[0].sbjct_end==131
assert record.alignments[41].hsps[0].query=="ctttgttatatttag"
assert record.alignments[41].hsps[0].match=="|||||||||||||||"
assert record.alignments[41].hsps[0].sbjct=="ctttgttatatttag"
assert record.alignments[41].hsps[0].query_start==288
assert record.alignments[41].hsps[0].query_end==302
assert record.alignments[41].hsps[0].sbjct_start==74
assert record.alignments[41].hsps[0].sbjct_end==60
assert record.alignments[42].hsps[0].query=="aggagcagagaggaa"
assert record.alignments[42].hsps[0].match=="|||||||||||||||"
assert record.alignments[42].hsps[0].sbjct=="aggagcagagaggaa"
assert record.alignments[42].hsps[0].query_start==91
assert record.alignments[42].hsps[0].query_end==105
assert record.alignments[42].hsps[0].sbjct_start==216
assert record.alignments[42].hsps[0].sbjct_end==230
assert record.alignments[43].hsps[0].query=="aggcagcagggatgg"
assert record.alignments[43].hsps[0].match=="|||||||||||||||"
assert record.alignments[43].hsps[0].sbjct=="aggcagcagggatgg"
assert record.alignments[43].hsps[0].query_start==125
assert record.alignments[43].hsps[0].query_end==139
assert record.alignments[43].hsps[0].sbjct_start==273
assert record.alignments[43].hsps[0].sbjct_end==259
assert record.alignments[44].hsps[0].query=="aggagagggccaggc"
assert record.alignments[44].hsps[0].match=="|||||||||||||||"
assert record.alignments[44].hsps[0].sbjct=="aggagagggccaggc"
assert record.alignments[44].hsps[0].query_start==114
assert record.alignments[44].hsps[0].query_end==128
assert record.alignments[44].hsps[0].sbjct_start==124
assert record.alignments[44].hsps[0].sbjct_end==138
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Feb 11, 2000  2:37 PM',)]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==6844
assert record.num_sequences==87792
assert record.num_extends==6844
assert record.num_good_extends==1887
assert record.num_seqs_better_e==51
assert record.query_length==331
assert record.database_length==31998854
assert record.effective_hsp_length==17
assert record.effective_query_length==314
assert record.effective_database_length==30506390
assert record.effective_search_space==9579006460
assert record.effective_search_space_used==9579006460
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(10,19.8)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(15,30.2)

handle = open('Blast/bt051')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|859351|gb|G06106|G06106 human STS WI-6344."
assert record.query_letters==183
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==6
assert record.descriptions[0].title=="gi|859351|gb|G06106|G06106 human STS WI-6344."
assert record.descriptions[0].score==327
assert record.descriptions[0].e==1e-89
assert record.descriptions[1].title=="gi|1341350|gb|G21024|G21024 human STS WI-30979."
assert record.descriptions[1].score==32
assert record.descriptions[1].e==1.000000
assert record.descriptions[2].title=="gi|6126285|gb|G60966.1|G60966 SHGC-84377 Human Homo sapiens STS..."
assert record.descriptions[2].score==30
assert record.descriptions[2].e==4.100000
assert record.descriptions[3].title=="gi|5222421|gb|G51244.1|G51244 SHGC-80725 Human Homo sapiens STS..."
assert record.descriptions[3].score==30
assert record.descriptions[3].e==4.100000
assert record.descriptions[4].title=="gi|1340656|gb|G20319|G20319 human STS A005L39."
assert record.descriptions[4].score==30
assert record.descriptions[4].e==4.100000
assert record.descriptions[5].title=="gi|860526|gb|G07281|G07281 human STS WI-9430."
assert record.descriptions[5].score==30
assert record.descriptions[5].e==4.100000
assert len(record.alignments)==0
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Feb 11, 2000  2:37 PM',)]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==1902
assert record.num_sequences==87792
assert record.num_extends==1902
assert record.num_good_extends==481
assert record.num_seqs_better_e==8
assert record.query_length==183
assert record.database_length==31998854
assert record.effective_hsp_length==16
assert record.effective_query_length==167
assert record.effective_database_length==30594182
assert record.effective_search_space==5109228394
assert record.effective_search_space_used==5109228394
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(10,19.8)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(15,30.2)

handle = open('Blast/bt052')
record = parser.parse(handle)
assert record.application=="BLASTX"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|1347369|gb|G25137|G25137 human STS EST48004."
assert record.query_letters==556
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==4
assert record.descriptions[0].title=="gi|1731448|sp|P54103|ZRF1_MOUSE ZUOTIN RELATED FACTOR-1"
assert record.descriptions[0].score==87
assert record.descriptions[0].e==3.0000000000000001e-17
assert record.descriptions[1].title=="gi|465911|sp|P34454|YMA9_CAEEL HYPOTHETICAL 31.6 KD PROTEIN F54..."
assert record.descriptions[1].score==42
assert record.descriptions[1].e==0.001000
assert record.descriptions[2].title=="gi|2494160|sp|Q61712|MTJ1_MOUSE DNAJ PROTEIN HOMOLOG MTJ1"
assert record.descriptions[2].score==37
assert record.descriptions[2].e==0.033000
assert record.descriptions[3].title=="gi|1730688|sp|P53745|YN8X_YEAST HYPOTHETICAL 68.1 KD PROTEIN IN..."
assert record.descriptions[3].score==29
assert record.descriptions[3].e==7.400000
assert len(record.alignments)==4
assert record.alignments[0].title==">gi|1731448|sp|P54103|ZRF1_MOUSE ZUOTIN RELATED FACTOR-1"
assert record.alignments[0].length==514
assert record.alignments[1].title==">gi|465911|sp|P34454|YMA9_CAEEL HYPOTHETICAL 31.6 KD PROTEIN F54F2.9 IN CHROMOSOME III"
assert record.alignments[1].length==275
assert record.alignments[2].title==">gi|2494160|sp|Q61712|MTJ1_MOUSE DNAJ PROTEIN HOMOLOG MTJ1"
assert record.alignments[2].length==552
assert record.alignments[3].title==">gi|1730688|sp|P53745|YN8X_YEAST HYPOTHETICAL 68.1 KD PROTEIN IN BIO3-FRE4 INTERGENIC REGION"
assert record.alignments[3].length==580
assert record.alignments[0].hsps[0].score==211
assert record.alignments[0].hsps[0].bits==86.6
assert record.alignments[0].hsps[0].expect==3e-17
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==96
assert record.alignments[1].hsps[0].bits==41.8
assert record.alignments[1].hsps[0].expect==0.001
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==83
assert record.alignments[2].hsps[0].bits==36.7
assert record.alignments[2].hsps[0].expect==0.033
assert record.alignments[2].hsps[1].score==69
assert record.alignments[2].hsps[1].bits==31.3
assert record.alignments[2].hsps[1].expect==1.5
assert len(record.alignments[2].hsps)==2
assert record.alignments[3].hsps[0].score==63
assert record.alignments[3].hsps[0].bits==29.0
assert record.alignments[3].hsps[0].expect==7.4
assert len(record.alignments[3].hsps)==1
assert record.alignments[0].hsps[0].identities==(41, 47)
assert record.alignments[0].hsps[0].positives==(44, 47)
assert record.alignments[1].hsps[0].identities==(30, 122)
assert record.alignments[1].hsps[0].positives==(54, 122)
assert record.alignments[2].hsps[0].identities==(17, 36)
assert record.alignments[2].hsps[0].positives==(19, 36)
assert record.alignments[2].hsps[1].identities==(18, 50)
assert record.alignments[2].hsps[1].positives==(26, 50)
assert record.alignments[3].hsps[0].identities==(27, 99)
assert record.alignments[3].hsps[0].positives==(41, 99)
assert record.alignments[0].hsps[0].frame==("+1", )
assert record.alignments[1].hsps[0].frame==("+1", )
assert record.alignments[2].hsps[0].frame==("+1", )
assert record.alignments[2].hsps[1].frame==("+1", )
assert record.alignments[3].hsps[0].frame==("+1", )
assert record.alignments[0].hsps[0].query=="DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAK"
assert record.alignments[0].hsps[0].match=="DLQLLIKAVNLFPAG NSRW+VIANYMNIHSSSGVKRTAKDVI + +"
assert record.alignments[0].hsps[0].sbjct=="DLQLLIKAVNLFPAGRNSRWDVIANYMNIHSSSGVKRTAKDVISEVR"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==141
assert record.alignments[0].hsps[0].sbjct_start==458
assert record.alignments[0].hsps[0].sbjct_end==504
assert record.alignments[1].hsps[0].query=="FPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFXGPYTDFTPXTTEXQKLXEQALNTYPVNTXERWXXIAVAVPGRXKE"
assert record.alignments[1].hsps[0].match=="+PAGT +RWE +   +N        R+A+DVI  A  ++++   +++D  K      ++   V  ++++                 +   QK  E AL  YP  T ERW  I+  +  + K+"
assert record.alignments[1].hsps[0].sbjct=="YPAGTPNRWEQMGRVLN--------RSAEDVIAMAGKMKQM---KQEDYTKLLMTTIQQSVPVEEKSED---------------DWSQAEQKAFETALQKYPKGTDERWERISEEIGSKTKK"
assert record.alignments[1].hsps[0].query_start==34
assert record.alignments[1].hsps[0].query_end==399
assert record.alignments[1].hsps[0].sbjct_start==159
assert record.alignments[1].hsps[0].sbjct_end==254
assert record.alignments[2].hsps[0].query=="TTEXQKLXEQALNTYPVNTXERWXXIAVAVPGRXKE"
assert record.alignments[2].hsps[0].match=="T   QKL E AL  YP    +RW  IA  VP + KE"
assert record.alignments[2].hsps[0].sbjct=="TQSQQKLLELALQQYPKGASDRWDKIAKCVPSKSKE"
assert record.alignments[2].hsps[0].query_start==292
assert record.alignments[2].hsps[0].query_end==399
assert record.alignments[2].hsps[0].sbjct_start==496
assert record.alignments[2].hsps[0].sbjct_end==531
assert record.alignments[2].hsps[1].query=="DLQLLIKAVNLFPAGTNSRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQ"
assert record.alignments[2].hsps[1].match=="DL  L +++  FP GT  RW+ IA+ +         R+  DV  KAK L+"
assert record.alignments[2].hsps[1].sbjct=="DLSQLTRSMVKFPGGTPGRWDKIAHELG--------RSVTDVTTKAKELK"
assert record.alignments[2].hsps[1].query_start==1
assert record.alignments[2].hsps[1].query_end==150
assert record.alignments[2].hsps[1].sbjct_start==332
assert record.alignments[2].hsps[1].sbjct_end==373
assert record.alignments[3].hsps[0].query=="SRWEVIANYMNIHSSSGVKRTAKDVIGKAKSLQKLDPHQKDDINKKAFDKFKKEHGVVPQADNATPSERFXGPYTDFTPXTTEXQKLXEQALNTYPVNT"
assert record.alignments[3].hsps[0].match=="+RW+   +Y        V R+ KDV   ++SL  LD +QK     +A       H +          E    PY +FT   +      EQ+ N +PV+T"
assert record.alignments[3].hsps[0].sbjct=="NRWKSFISY--------VTRSRKDVKTVSRSLSNLDLYQKCSKEIRADQDISLLHSI----------ETKLFPYINFTALNS------EQSHNFWPVHT"
assert record.alignments[3].hsps[0].query_start==52
assert record.alignments[3].hsps[0].query_end==348
assert record.alignments[3].hsps[0].sbjct_start==75
assert record.alignments[3].hsps[0].sbjct_end==149
assert record.database_name==['data/swissprot']
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.ka_params==[0.318, 0.135, 0.401]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==23174157
assert record.num_sequences==82258
assert record.num_extends==387821
assert record.num_good_extends==980
assert record.num_seqs_better_e==8
assert record.hsps_no_gap==3
assert record.hsps_prelim_gapped==1
assert record.hsps_gapped==7
assert record.query_length==185
assert record.database_length==29652561
assert record.effective_hsp_length==49
assert record.effective_query_length==135
assert record.effective_database_length==25621919
assert record.effective_search_space==3458959065
assert record.effective_search_space_used==3458959065
assert record.frameshift==('50,','0.1')
assert record.threshold==12
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(62,28.6)

handle = open('Blast/bt053')
record = parser.parse(handle)
assert record.application=="BLASTX"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|1347782|gb|G25550|G25550 human STS\nEST47652.\x01gi|1592937|gb|G29386|G29386 human STS SHGC-32770"
assert record.query_letters==379
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==0
assert len(record.alignments)==0
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.ka_params==[0.318, 0.135, 0.401]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==14686001
assert record.num_sequences==82258
assert record.num_extends==235383
assert record.num_good_extends==396
assert record.num_seqs_better_e==0
assert record.hsps_no_gap==0
assert record.hsps_prelim_gapped==0
assert record.hsps_gapped==0
assert record.query_length==126
assert record.database_length==29652561
assert record.effective_hsp_length==48
assert record.effective_query_length==77
assert record.effective_database_length==25704177
assert record.effective_search_space==1979221629
assert record.effective_search_space_used==1979221629
assert record.frameshift==('50,','0.1')
assert record.threshold==12
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(60,27.8)

handle = open('Blast/bt054')
record = parser.parse(handle)
assert record.application=="TBLASTN"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|729325|sp|P39483|DHG2_BACME GLUCOSE 1-DEHYDROGENASE II\n(GLCDH-II)"
assert record.query_letters==261
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==7
assert record.descriptions[0].title=="gi|3820341|emb|AJ229891|KLAJ9891 Kluyveromyces lactis DNA fragm..."
assert record.descriptions[0].score==47
assert record.descriptions[0].e==0.000010
assert record.descriptions[1].title=="gi|1375419|gb|G27169|G27169 human STS SHGC-31983."
assert record.descriptions[1].score==43
assert record.descriptions[1].e==0.000100
assert record.descriptions[2].title=="gi|3819804|emb|AJ230012|KLAJ0012 Kluyveromyces lactis DNA fragm..."
assert record.descriptions[2].score==39
assert record.descriptions[2].e==0.002000
assert record.descriptions[3].title=="gi|1375215|gb|G26965|G26965 human STS SHGC-31083."
assert record.descriptions[3].score==31
assert record.descriptions[3].e==0.730000
assert record.descriptions[4].title=="gi|5714409|gb|AF106665.1|AF106665 Mus musculus chromosome 6 clo..."
assert record.descriptions[4].score==29
assert record.descriptions[4].e==2.200000
assert record.descriptions[5].title=="gi|177714|gb|L09988|HUM4STS889 Human Chromosome 4 (clone p4-109..."
assert record.descriptions[5].score==29
assert record.descriptions[5].e==2.200000
assert record.descriptions[6].title=="gi|1341648|gb|G21322|G21322 human STS WI-12250."
assert record.descriptions[6].score==29
assert record.descriptions[6].e==3.700000
assert len(record.alignments)==7
assert record.alignments[0].title==">gi|3820341|emb|AJ229891|KLAJ9891 Kluyveromyces lactis DNA fragment for sequence tagged site, clone okam5d07r [Kluyveromyces lactis]"
assert record.alignments[0].length==230
assert record.alignments[1].title==">gi|1375419|gb|G27169|G27169 human STS SHGC-31983."
assert record.alignments[1].length==594
assert record.alignments[2].title==">gi|3819804|emb|AJ230012|KLAJ0012 Kluyveromyces lactis DNA fragment for sequence tagged site, clone okam6d01d [Kluyveromyces lactis]"
assert record.alignments[2].length==199
assert record.alignments[3].title==">gi|1375215|gb|G26965|G26965 human STS SHGC-31083."
assert record.alignments[3].length==268
assert record.alignments[4].title==">gi|5714409|gb|AF106665.1|AF106665 Mus musculus chromosome 6 clone D6wum9 map between Nkrp1 and Prp strain C57BL/6J, sequence tagged site"
assert record.alignments[4].length==299
assert record.alignments[5].title==">gi|177714|gb|L09988|HUM4STS889 Human Chromosome 4 (clone p4-1095) STS4-889."
assert record.alignments[5].length==412
assert record.alignments[6].title==">gi|1341648|gb|G21322|G21322 human STS WI-12250."
assert record.alignments[6].length==586
assert record.alignments[0].hsps[0].score==109
assert record.alignments[0].hsps[0].bits==46.9
assert record.alignments[0].hsps[0].expect==1e-05
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==100
assert record.alignments[1].hsps[0].bits==43.4
assert record.alignments[1].hsps[0].expect==1e-04
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==90
assert record.alignments[2].hsps[0].bits==39.5
assert record.alignments[2].hsps[0].expect==0.002
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==68
assert record.alignments[3].hsps[0].bits==30.9
assert record.alignments[3].hsps[0].expect==0.73
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==64
assert record.alignments[4].hsps[0].bits==29.3
assert record.alignments[4].hsps[0].expect==2.2
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==64
assert record.alignments[5].hsps[0].bits==29.3
assert record.alignments[5].hsps[0].expect==2.2
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==62
assert record.alignments[6].hsps[0].bits==28.6
assert record.alignments[6].hsps[0].expect==3.7
assert len(record.alignments[6].hsps)==1
assert record.alignments[0].hsps[0].identities==(25, 72)
assert record.alignments[0].hsps[0].positives==(44, 72)
assert record.alignments[0].hsps[0].gaps==(3, 72)
assert record.alignments[1].hsps[0].identities==(21, 73)
assert record.alignments[1].hsps[0].positives==(34, 73)
assert record.alignments[2].hsps[0].identities==(18, 49)
assert record.alignments[2].hsps[0].positives==(26, 49)
assert record.alignments[3].hsps[0].identities==(12, 37)
assert record.alignments[3].hsps[0].positives==(19, 37)
assert record.alignments[4].hsps[0].identities==(17, 55)
assert record.alignments[4].hsps[0].positives==(32, 55)
assert record.alignments[4].hsps[0].gaps==(2, 55)
assert record.alignments[5].hsps[0].identities==(14, 34)
assert record.alignments[5].hsps[0].positives==(22, 34)
assert record.alignments[6].hsps[0].identities==(16, 39)
assert record.alignments[6].hsps[0].positives==(20, 39)
assert record.alignments[6].hsps[0].gaps==(1, 39)
assert record.alignments[0].hsps[0].frame==("+1", )
assert record.alignments[1].hsps[0].frame==("-1", )
assert record.alignments[2].hsps[0].frame==("-1", )
assert record.alignments[3].hsps[0].frame==("-1", )
assert record.alignments[4].hsps[0].frame==("-2", )
assert record.alignments[5].hsps[0].frame==("-1", )
assert record.alignments[6].hsps[0].frame==("+1", )
assert record.alignments[0].hsps[0].query=="NWNQVIDTNLTGAFLGSREAIKYFVEN---DIKGNVINMSSVHEMIPWPLFVHYAASKGGMKLMTETLALEYAPK"
assert record.alignments[0].hsps[0].match=="+W QVIDTN+ G F   + A+     +   D +  V+N+S+V+ ++  P    Y A+K  +  +T+++ALEYA +"
assert record.alignments[0].hsps[0].sbjct=="SWRQVIDTNINGTFYTLKYALPLMESSSSPDSEAAVVNLSAVNGLVGIPGISPYTATKHAVIGITQSVALEYAER"
assert record.alignments[0].hsps[0].query_start==108
assert record.alignments[0].hsps[0].query_end==179
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==225
assert record.alignments[1].hsps[0].query=="APKGIRVNNIGPGAIDTPINAEKFADPEQRADVESMIPMGYIGKPEEIASVAAFLASSQASYVTGITLFADGG"
assert record.alignments[1].hsps[0].match=="AP   RVN + P  + T +    + DP +   + +  P+G   + E +     FL S ++   TG TL  +GG"
assert record.alignments[1].hsps[0].sbjct=="APXRXRVNAVXPXVVMTSMGQATWXDPXKAXTMLNRXPLGXFAEVEHVVKAILFLLSDRSGMTTGSTLPVEGG"
assert record.alignments[1].hsps[0].query_start==177
assert record.alignments[1].hsps[0].query_end==249
assert record.alignments[1].hsps[0].sbjct_start==312
assert record.alignments[1].hsps[0].sbjct_end==94
assert record.alignments[2].hsps[0].query=="FADPEQRADVESMIPMGYIGKPEEIASVAAFLASSQASYVTGITLFADG"
assert record.alignments[2].hsps[0].match=="F D + +    S+IPMG  G P+E+     + AS  ++Y TG  L  DG"
assert record.alignments[2].hsps[0].sbjct=="FVDEDLKNKWHSLIPMGREGLPQELVGAYLYFASDASTYTTGSDLLVDG"
assert record.alignments[2].hsps[0].query_start==200
assert record.alignments[2].hsps[0].query_end==248
assert record.alignments[2].hsps[0].sbjct_start==157
assert record.alignments[2].hsps[0].sbjct_end==11
assert record.alignments[3].hsps[0].query=="PMGYIGKPEEIASVAAFLASSQASYVTGITLFADGGM"
assert record.alignments[3].hsps[0].match=="PMG  G PE++  V A      +  +TG ++   GG+"
assert record.alignments[3].hsps[0].sbjct=="PMGXXGDPEDVXDVXAXXXXEXSGXITGTSVEVTGGL"
assert record.alignments[3].hsps[0].query_start==214
assert record.alignments[3].hsps[0].query_end==250
assert record.alignments[3].hsps[0].sbjct_start==268
assert record.alignments[3].hsps[0].sbjct_end==158
assert record.alignments[4].hsps[0].query=="NMSSVHEMIPWPLFVHYAASKGGMKLMTETL--ALEYAPKGIRVNNIGPGAIDTPIN"
assert record.alignments[4].hsps[0].match=="++S    +I +P F+    S  G  L+  +L  A+ + P GI V+++GP ++ T +N"
assert record.alignments[4].hsps[0].sbjct=="SLSPTQYLIMFPSFLPCPLSHPGPFLLPSSLVIAVFFLPNGIEVSSLGPFSLRTLLN"
assert record.alignments[4].hsps[0].query_start==142
assert record.alignments[4].hsps[0].query_end==196
assert record.alignments[4].hsps[0].sbjct_start==172
assert record.alignments[4].hsps[0].sbjct_end==2
assert record.alignments[5].hsps[0].query=="DKVVVVTGGSKGLGRAMAVRFGQEQSKVVVNYRS"
assert record.alignments[5].hsps[0].match=="DKV  V GGS+G+GRA+A    ++  ++ V  R+"
assert record.alignments[5].hsps[0].sbjct=="DKVCAVFGGSRGIGRAVAQLMARKGYRLAVIARN"
assert record.alignments[5].hsps[0].query_start==7
assert record.alignments[5].hsps[0].query_end==40
assert record.alignments[5].hsps[0].sbjct_start==316
assert record.alignments[5].hsps[0].sbjct_end==215
assert record.alignments[6].hsps[0].query=="PVPSHELSLENWNQ-VIDTNLTGAFLGSREAIKYFVENDI"
assert record.alignments[6].hsps[0].match=="PVP  ELS  +W+Q  + T+ T  F  S     YF  N I"
assert record.alignments[6].hsps[0].sbjct=="PVPMQELSKVHWSQFFLTTSPTMTFFFSHYLANYFFRNSI"
assert record.alignments[6].hsps[0].query_start==98
assert record.alignments[6].hsps[0].query_end==136
assert record.alignments[6].hsps[0].sbjct_start==220
assert record.alignments[6].hsps[0].sbjct_end==339
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Feb 11, 2000  2:37 PM',)]
assert record.ka_params==[0.315, 0.134, 0.378]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==6165875
assert record.num_sequences==87792
assert record.num_extends==55665
assert record.num_good_extends==148
assert record.num_seqs_better_e==14
assert record.hsps_no_gap==5
assert record.hsps_prelim_gapped==2
assert record.hsps_gapped==7
assert record.query_length==261
assert record.database_length==10666284
assert record.effective_hsp_length==50
assert record.effective_query_length==211
assert record.effective_database_length==6276684
assert record.effective_search_space==1324380324
assert record.effective_search_space_used==1324380324
assert record.frameshift==('50,','0.1')
assert record.threshold==13
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(42,22.0)
assert record.blast_cutoff==(58,27.0)

handle = open('Blast/bt055')
record = parser.parse(handle)
assert record.application=="TBLASTN"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|127420|sp|P19888|MTBA_BACAR MODIFICATION METHYLASE BANI\n(CYTOSINE-SPECIFIC METHYLTRANSFERASE BANI) (M.BANI)"
assert record.query_letters==428
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==0
assert len(record.alignments)==0
assert record.database_name==['data/sts']
assert record.posted_date==[('Feb 11, 2000  2:37 PM',)]
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.ka_params==[0.320, 0.140, 0.403]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==13588598
assert record.num_sequences==87792
assert record.num_extends==162273
assert record.num_good_extends==546
assert record.num_seqs_better_e==0
assert record.hsps_no_gap==0
assert record.hsps_prelim_gapped==0
assert record.hsps_gapped==0
assert record.query_length==428
assert record.database_length==10666284
assert record.effective_hsp_length==48
assert record.effective_query_length==380
assert record.effective_database_length==6452268
assert record.effective_search_space==2451861840
assert record.effective_search_space_used==2451861840
assert record.frameshift==('50,','0.1')
assert record.threshold==13
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.4)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.8)
assert record.blast_cutoff==(61,28.2)

handle = open('Blast/bt056')
record = parser.parse(handle)
assert record.application=="TBLASTX"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|1348853|gb|G26621|G26621 human STS\nSTS_D12006.\x01gi|1396339|gb|G27620|G27620 human STS SHGC-32705."
assert record.query_letters==615
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==19
assert record.descriptions[0].title=="gi|1348853|gb|G26621|G26621 human STS STS_D12006. >gi|1396339|g..."
assert record.descriptions[0].score==398
assert record.descriptions[0].e==1.0000000000000001e-111
assert record.descriptions[1].title=="gi|1348016|gb|G25784|G25784 human STS EST47998."
assert record.descriptions[1].score==302
assert record.descriptions[1].e==9.9999999999999996e-83
assert record.descriptions[2].title=="gi|3403105|gb|G41148|G41148 Z7324 Zebrafish AB Danio rerio STS ..."
assert record.descriptions[2].score==31
assert record.descriptions[2].e==0.600000
assert record.descriptions[3].title=="gi|1234967|emb|Z53521|HSB316XA9 H.sapiens (D22S1166) DNA segmen..."
assert record.descriptions[3].score==30
assert record.descriptions[3].e==0.830000
assert record.descriptions[4].title=="gi|4185670|gb|G42865|G42865 Xq3845 KWOK Homo sapiens STS genomi..."
assert record.descriptions[4].score==30
assert record.descriptions[4].e==1.100000
assert record.descriptions[5].title=="gi|4757419|gb|G49246.1|G49246 stbK116F5_30376 chromosome 22 gen..."
assert record.descriptions[5].score==29
assert record.descriptions[5].e==1.600000
assert record.descriptions[6].title=="gi|6120694|gb|G55375.1|G55375 SHGC-100697 Human Homo sapiens ST..."
assert record.descriptions[6].score==29
assert record.descriptions[6].e==3.000000
assert record.descriptions[7].title=="gi|5225124|gb|G53947.1|G53947 SHGC-85304 Human Homo sapiens STS..."
assert record.descriptions[7].score==29
assert record.descriptions[7].e==3.000000
assert record.descriptions[8].title=="gi|1311530|gb|L77996|HUMSWX945 Human chromosome X STS sWXD945, ..."
assert record.descriptions[8].score==28
assert record.descriptions[8].e==4.100000
assert record.descriptions[9].title=="gi|4631600|dbj|AU046965.1|AU046965 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[9].score==28
assert record.descriptions[9].e==4.100000
assert record.descriptions[10].title=="gi|4631518|dbj|AU046883.1|AU046883 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[10].score==28
assert record.descriptions[10].e==4.100000
assert record.descriptions[11].title=="gi|2226478|gb|G33174|G33174 human STS SHGC-6097 clone pG-2470"
assert record.descriptions[11].score==28
assert record.descriptions[11].e==4.100000
assert record.descriptions[12].title=="gi|2734624|gb|G36957|G36957 SHGC-56642 Human Homo sapiens STS cDNA"
assert record.descriptions[12].score==28
assert record.descriptions[12].e==4.100000
assert record.descriptions[13].title=="gi|859804|gb|G06559|G06559 human STS WI-7401."
assert record.descriptions[13].score==27
assert record.descriptions[13].e==7.700000
assert record.descriptions[14].title=="gi|938611|gb|G08061|G08061 human STS CHLC.GGAA7E02.P7438 clone ..."
assert record.descriptions[14].score==27
assert record.descriptions[14].e==7.700000
assert record.descriptions[15].title=="gi|307789|gb|L18105|HUMUT1736 Human STS UT1736."
assert record.descriptions[15].score==27
assert record.descriptions[15].e==7.700000
assert record.descriptions[16].title=="gi|4492122|gb|G45831.1|G45831 Z4588_1 Zebrafish AB Danio rerio ..."
assert record.descriptions[16].score==27
assert record.descriptions[16].e==7.700000
assert record.descriptions[17].title=="gi|6121804|gb|G56635.1|G56635 SHGC-102032 Human Homo sapiens ST..."
assert record.descriptions[17].score==27
assert record.descriptions[17].e==7.700000
assert record.descriptions[18].title=="gi|4493143|gb|G46852.1|G46852 Z14841_1 Zebrafish AB Danio rerio..."
assert record.descriptions[18].score==27
assert record.descriptions[18].e==7.700000
assert len(record.alignments)==19
assert record.alignments[0].title==">gi|1348853|gb|G26621|G26621 human STS STS_D12006. >gi|1396339|gb|G27620|G27620 human STS SHGC-32705."
assert record.alignments[0].length==615
assert record.alignments[1].title==">gi|1348016|gb|G25784|G25784 human STS EST47998."
assert record.alignments[1].length==617
assert record.alignments[2].title==">gi|3403105|gb|G41148|G41148 Z7324 Zebrafish AB Danio rerio STS genomic"
assert record.alignments[2].length==351
assert record.alignments[3].title==">gi|1234967|emb|Z53521|HSB316XA9 H.sapiens (D22S1166) DNA segment containing (CA) repeat; clone AFMb316xa9; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[3].length==345
assert record.alignments[4].title==">gi|4185670|gb|G42865|G42865 Xq3845 KWOK Homo sapiens STS genomic, sequence tagged site [Homo sapiens]"
assert record.alignments[4].length==1200
assert record.alignments[5].title==">gi|4757419|gb|G49246.1|G49246 stbK116F5_30376 chromosome 22 genomic clone Homo sapiens STS genomic clone 116F5, sequence tagged site"
assert record.alignments[5].length==375
assert record.alignments[6].title==">gi|6120694|gb|G55375.1|G55375 SHGC-100697 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[6].length==460
assert record.alignments[7].title==">gi|5225124|gb|G53947.1|G53947 SHGC-85304 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[7].length==444
assert record.alignments[8].title==">gi|1311530|gb|L77996|HUMSWX945 Human chromosome X STS sWXD945, single read."
assert record.alignments[8].length==196
assert record.alignments[9].title==">gi|4631600|dbj|AU046965.1|AU046965 Rattus norvegicus, OTSUKA clone, 108a02, microsatellite sequence, sequence tagged site"
assert record.alignments[9].length==330
assert record.alignments[10].title==">gi|4631518|dbj|AU046883.1|AU046883 Rattus norvegicus, OTSUKA clone, 085f03, microsatellite sequence, sequence tagged site"
assert record.alignments[10].length==351
assert record.alignments[11].title==">gi|2226478|gb|G33174|G33174 human STS SHGC-6097 clone pG-2470"
assert record.alignments[11].length==299
assert record.alignments[12].title==">gi|2734624|gb|G36957|G36957 SHGC-56642 Human Homo sapiens STS cDNA"
assert record.alignments[12].length==466
assert record.alignments[13].title==">gi|859804|gb|G06559|G06559 human STS WI-7401."
assert record.alignments[13].length==3280
assert record.alignments[14].title==">gi|938611|gb|G08061|G08061 human STS CHLC.GGAA7E02.P7438 clone GGAA7E02"
assert record.alignments[14].length==338
assert record.alignments[15].title==">gi|307789|gb|L18105|HUMUT1736 Human STS UT1736."
assert record.alignments[15].length==355
assert record.alignments[16].title==">gi|4492122|gb|G45831.1|G45831 Z4588_1 Zebrafish AB Danio rerio STS genomic clone Z4588 5', sequence tagged site"
assert record.alignments[16].length==398
assert record.alignments[17].title==">gi|6121804|gb|G56635.1|G56635 SHGC-102032 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[17].length==541
assert record.alignments[18].title==">gi|4493143|gb|G46852.1|G46852 Z14841_1 Zebrafish AB Danio rerio STS genomic clone Z14841 5', sequence tagged site"
assert record.alignments[18].length==291
assert record.alignments[0].hsps[0].score==796
assert record.alignments[0].hsps[0].bits==367
assert record.alignments[0].hsps[0].expect==1e-102
assert record.alignments[0].hsps[1].score==759
assert record.alignments[0].hsps[1].bits==350
assert record.alignments[0].hsps[1].expect==3e-97
assert record.alignments[0].hsps[2].score==387
assert record.alignments[0].hsps[2].bits==180
assert record.alignments[0].hsps[2].expect==9e-91
assert record.alignments[0].hsps[3].score==368
assert record.alignments[0].hsps[3].bits==171
assert record.alignments[0].hsps[3].expect==9e-91
assert record.alignments[0].hsps[4].score==864
assert record.alignments[0].hsps[4].bits==398
assert record.alignments[0].hsps[4].expect==1e-111
assert record.alignments[0].hsps[5].score==846
assert record.alignments[0].hsps[5].bits==390
assert record.alignments[0].hsps[5].expect==1e-109
assert record.alignments[0].hsps[6].score==684
assert record.alignments[0].hsps[6].bits==316
assert record.alignments[0].hsps[6].expect==7e-87
assert len(record.alignments[0].hsps)==7
assert record.alignments[1].hsps[0].score==366
assert record.alignments[1].hsps[0].bits==170
assert record.alignments[1].hsps[0].expect==3e-63
assert record.alignments[1].hsps[1].score==188
assert record.alignments[1].hsps[1].bits==89.0
assert record.alignments[1].hsps[1].expect==3e-63
assert record.alignments[1].hsps[2].score==590
assert record.alignments[1].hsps[2].bits==273
assert record.alignments[1].hsps[2].expect==7e-74
assert record.alignments[1].hsps[3].score==593
assert record.alignments[1].hsps[3].bits==274
assert record.alignments[1].hsps[3].expect==8e-76
assert record.alignments[1].hsps[4].score==53
assert record.alignments[1].hsps[4].bits==27.2
assert record.alignments[1].hsps[4].expect==8e-76
assert record.alignments[1].hsps[5].score==653
assert record.alignments[1].hsps[5].bits==302
assert record.alignments[1].hsps[5].expect==1e-82
assert record.alignments[1].hsps[6].score==598
assert record.alignments[1].hsps[6].bits==276
assert record.alignments[1].hsps[6].expect==5e-75
assert record.alignments[1].hsps[7].score==628
assert record.alignments[1].hsps[7].bits==290
assert record.alignments[1].hsps[7].expect==4e-79
assert len(record.alignments[1].hsps)==8
assert record.alignments[2].hsps[0].score==61
assert record.alignments[2].hsps[0].bits==30.8
assert record.alignments[2].hsps[0].expect==0.60
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==60
assert record.alignments[3].hsps[0].bits==30.4
assert record.alignments[3].hsps[0].expect==0.83
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==59
assert record.alignments[4].hsps[0].bits==29.9
assert record.alignments[4].hsps[0].expect==1.1
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==58
assert record.alignments[5].hsps[0].bits==29.5
assert record.alignments[5].hsps[0].expect==1.6
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==56
assert record.alignments[6].hsps[0].bits==28.6
assert record.alignments[6].hsps[0].expect==3.0
assert len(record.alignments[6].hsps)==1
assert record.alignments[7].hsps[0].score==56
assert record.alignments[7].hsps[0].bits==28.6
assert record.alignments[7].hsps[0].expect==3.0
assert len(record.alignments[7].hsps)==1
assert record.alignments[8].hsps[0].score==55
assert record.alignments[8].hsps[0].bits==28.1
assert record.alignments[8].hsps[0].expect==4.1
assert len(record.alignments[8].hsps)==1
assert record.alignments[9].hsps[0].score==55
assert record.alignments[9].hsps[0].bits==28.1
assert record.alignments[9].hsps[0].expect==4.1
assert len(record.alignments[9].hsps)==1
assert record.alignments[10].hsps[0].score==55
assert record.alignments[10].hsps[0].bits==28.1
assert record.alignments[10].hsps[0].expect==4.1
assert len(record.alignments[10].hsps)==1
assert record.alignments[11].hsps[0].score==55
assert record.alignments[11].hsps[0].bits==28.1
assert record.alignments[11].hsps[0].expect==4.1
assert len(record.alignments[11].hsps)==1
assert record.alignments[12].hsps[0].score==55
assert record.alignments[12].hsps[0].bits==28.1
assert record.alignments[12].hsps[0].expect==4.1
assert len(record.alignments[12].hsps)==1
assert record.alignments[13].hsps[0].score==53
assert record.alignments[13].hsps[0].bits==27.2
assert record.alignments[13].hsps[0].expect==7.7
assert len(record.alignments[13].hsps)==1
assert record.alignments[14].hsps[0].score==53
assert record.alignments[14].hsps[0].bits==27.2
assert record.alignments[14].hsps[0].expect==7.7
assert len(record.alignments[14].hsps)==1
assert record.alignments[15].hsps[0].score==53
assert record.alignments[15].hsps[0].bits==27.2
assert record.alignments[15].hsps[0].expect==7.7
assert len(record.alignments[15].hsps)==1
assert record.alignments[16].hsps[0].score==53
assert record.alignments[16].hsps[0].bits==27.2
assert record.alignments[16].hsps[0].expect==7.7
assert len(record.alignments[16].hsps)==1
assert record.alignments[17].hsps[0].score==53
assert record.alignments[17].hsps[0].bits==27.2
assert record.alignments[17].hsps[0].expect==7.7
assert len(record.alignments[17].hsps)==1
assert record.alignments[18].hsps[0].score==53
assert record.alignments[18].hsps[0].bits==27.2
assert record.alignments[18].hsps[0].expect==7.7
assert len(record.alignments[18].hsps)==1
assert record.alignments[0].hsps[0].identities==(192, 200)
assert record.alignments[0].hsps[0].positives==(192, 200)
assert record.alignments[0].hsps[1].identities==(195, 205)
assert record.alignments[0].hsps[1].positives==(195, 205)
assert record.alignments[0].hsps[2].identities==(74, 74)
assert record.alignments[0].hsps[2].positives==(74, 74)
assert record.alignments[0].hsps[3].identities==(114, 114)
assert record.alignments[0].hsps[3].positives==(114, 114)
assert record.alignments[0].hsps[4].identities==(205, 205)
assert record.alignments[0].hsps[4].positives==(205, 205)
assert record.alignments[0].hsps[5].identities==(196, 196)
assert record.alignments[0].hsps[5].positives==(196, 196)
assert record.alignments[0].hsps[6].identities==(146, 146)
assert record.alignments[0].hsps[6].positives==(146, 146)
assert record.alignments[1].hsps[0].identities==(71, 74)
assert record.alignments[1].hsps[0].positives==(71, 74)
assert record.alignments[1].hsps[1].identities==(42, 67)
assert record.alignments[1].hsps[1].positives==(43, 67)
assert record.alignments[1].hsps[2].identities==(121, 133)
assert record.alignments[1].hsps[2].positives==(121, 133)
assert record.alignments[1].hsps[3].identities==(112, 131)
assert record.alignments[1].hsps[3].positives==(112, 131)
assert record.alignments[1].hsps[4].identities==(9, 13)
assert record.alignments[1].hsps[4].positives==(10, 13)
assert record.alignments[1].hsps[5].identities==(128, 157)
assert record.alignments[1].hsps[5].positives==(132, 157)
assert record.alignments[1].hsps[6].identities==(122, 130)
assert record.alignments[1].hsps[6].positives==(122, 130)
assert record.alignments[1].hsps[7].identities==(119, 131)
assert record.alignments[1].hsps[7].positives==(120, 131)
assert record.alignments[2].hsps[0].identities==(11, 27)
assert record.alignments[2].hsps[0].positives==(18, 27)
assert record.alignments[3].hsps[0].identities==(10, 19)
assert record.alignments[3].hsps[0].positives==(13, 19)
assert record.alignments[4].hsps[0].identities==(10, 24)
assert record.alignments[4].hsps[0].positives==(14, 24)
assert record.alignments[5].hsps[0].identities==(15, 34)
assert record.alignments[5].hsps[0].positives==(17, 34)
assert record.alignments[6].hsps[0].identities==(9, 28)
assert record.alignments[6].hsps[0].positives==(16, 28)
assert record.alignments[7].hsps[0].identities==(10, 24)
assert record.alignments[7].hsps[0].positives==(15, 24)
assert record.alignments[8].hsps[0].identities==(9, 33)
assert record.alignments[8].hsps[0].positives==(19, 33)
assert record.alignments[9].hsps[0].identities==(12, 39)
assert record.alignments[9].hsps[0].positives==(19, 39)
assert record.alignments[10].hsps[0].identities==(8, 26)
assert record.alignments[10].hsps[0].positives==(15, 26)
assert record.alignments[11].hsps[0].identities==(13, 38)
assert record.alignments[11].hsps[0].positives==(19, 38)
assert record.alignments[12].hsps[0].identities==(10, 15)
assert record.alignments[12].hsps[0].positives==(11, 15)
assert record.alignments[13].hsps[0].identities==(9, 21)
assert record.alignments[13].hsps[0].positives==(13, 21)
assert record.alignments[14].hsps[0].identities==(12, 24)
assert record.alignments[14].hsps[0].positives==(13, 24)
assert record.alignments[15].hsps[0].identities==(12, 24)
assert record.alignments[15].hsps[0].positives==(13, 24)
assert record.alignments[16].hsps[0].identities==(7, 22)
assert record.alignments[16].hsps[0].positives==(13, 22)
assert record.alignments[17].hsps[0].identities==(9, 20)
assert record.alignments[17].hsps[0].positives==(13, 20)
assert record.alignments[18].hsps[0].identities==(8, 20)
assert record.alignments[18].hsps[0].positives==(13, 20)
assert record.alignments[0].hsps[0].frame==("+2", "+2")
assert record.alignments[0].hsps[1].frame==("+1", "+1")
assert record.alignments[0].hsps[2].frame==("+3", "+3")
assert record.alignments[0].hsps[3].frame==("+3", "+3")
assert record.alignments[0].hsps[4].frame==("-1", "-1")
assert record.alignments[0].hsps[5].frame==("-3", "-3")
assert record.alignments[0].hsps[6].frame==("-2", "-2")
assert record.alignments[1].hsps[0].frame==("+3", "+3")
assert record.alignments[1].hsps[1].frame==("+3", "+3")
assert record.alignments[1].hsps[2].frame==("+1", "+1")
assert record.alignments[1].hsps[3].frame==("+2", "+2")
assert record.alignments[1].hsps[4].frame==("+2", "+2")
assert record.alignments[1].hsps[5].frame==("-3", "-2")
assert record.alignments[1].hsps[6].frame==("-2", "-1")
assert record.alignments[1].hsps[7].frame==("-1", "-3")
assert record.alignments[2].hsps[0].frame==("+1", "+2")
assert record.alignments[3].hsps[0].frame==("-3", "+1")
assert record.alignments[4].hsps[0].frame==("-3", "-2")
assert record.alignments[5].hsps[0].frame==("-3", "-2")
assert record.alignments[6].hsps[0].frame==("+3", "+1")
assert record.alignments[7].hsps[0].frame==("+1", "-2")
assert record.alignments[8].hsps[0].frame==("+3", "-2")
assert record.alignments[9].hsps[0].frame==("+2", "-3")
assert record.alignments[10].hsps[0].frame==("+1", "+3")
assert record.alignments[11].hsps[0].frame==("+1", "+2")
assert record.alignments[12].hsps[0].frame==("-3", "-3")
assert record.alignments[13].hsps[0].frame==("-1", "-1")
assert record.alignments[14].hsps[0].frame==("-3", "+3")
assert record.alignments[15].hsps[0].frame==("-3", "-1")
assert record.alignments[16].hsps[0].frame==("+1", "+1")
assert record.alignments[17].hsps[0].frame==("-1", "+2")
assert record.alignments[18].hsps[0].frame==("-3", "+1")
assert record.alignments[0].hsps[0].query=="IRMPLHS*DSSFCPL*QEKWECMXXXXXXXXRPKRCLQPHPLNWP*LGLNALMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPLFWWTNNLGTPNPGRKNIQHYEXALX*S*NGFPKLVTHGPGXVXLXLXGLSSFQEFXSTVANPWGX*XXXXXXXXXFXXGXRXXXLXXXXGGCXXVXXVXXXWXXXF"
assert record.alignments[0].hsps[0].match=="IRMPLHS*DSSFCPL*QEKWECM        RPKRCLQPHPLNWP*LGLNALMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPLFWWTNNLGTPNPGRKNIQHYE AL *S*NGFPKLVTHGPG V L L GLSSFQEF STVANPWG *         F  G R   L    GGC  V  V   W   F"
assert record.alignments[0].hsps[0].sbjct=="IRMPLHS*DSSFCPL*QEKWECMQSSQKKQKRPKRCLQPHPLNWP*LGLNALMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPLFWWTNNLGTPNPGRKNIQHYEXALX*S*NGFPKLVTHGPGXVXLXLXGLSSFQEFXSTVANPWGX*XXXXXXXXXFXXGXRXXXLXXXXGGCXXVXXVXXXWXXXF"
assert record.alignments[0].hsps[0].query_start==2
assert record.alignments[0].hsps[0].query_end==601
assert record.alignments[0].hsps[0].sbjct_start==2
assert record.alignments[0].hsps[0].sbjct_end==601
assert record.alignments[0].hsps[1].query=="DQNAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKCTNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL*XSTXVILKWISKAGYTWPWXXXIXFXXAFFXXXXXXXXXXXXXXLXGXVX*XXXXSXGPXXXXVXXPXGWVXXGXFXFXXXXXXFXXXLG"
assert record.alignments[0].hsps[1].match=="DQNAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKCTNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL* ST VILKWISKAGYTWPW   I F  AFF              L G V *    S GP    V  P GWV  G F F      F   LG"
assert record.alignments[0].hsps[1].sbjct=="DQNAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKCTNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL*XSTXVILKWISKAGYTWPWXXXIXFXXAFFXSGVXVNGGKSXGXLXGXVX*XXXXSXGPXXXXVXXPXGWVXXGXFXFXXXXXXFXXXLG"
assert record.alignments[0].hsps[1].query_start==1
assert record.alignments[0].hsps[1].query_end==615
assert record.alignments[0].hsps[1].sbjct_start==1
assert record.alignments[0].hsps[1].sbjct_end==615
assert record.alignments[0].hsps[2].query=="SECPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*MH*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[0].hsps[2].match=="SECPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*MH*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[0].hsps[2].sbjct=="SECPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*MH*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[0].hsps[2].query_start==3
assert record.alignments[0].hsps[2].query_end==224
assert record.alignments[0].hsps[2].sbjct_start==3
assert record.alignments[0].hsps[2].sbjct_end==224
assert record.alignments[0].hsps[3].query=="YSHCSGGQIT*VLQILAGKIYSIMXQHSXNPKMDFQSWLHMALXXSY*XXXGFLXFRSXGQRWQIXGVXXWXGXLXXXXFXGAXGGXGXXSXXVGAXGXXXFXXXGXXXXXXXW"
assert record.alignments[0].hsps[3].match=="YSHCSGGQIT*VLQILAGKIYSIM QHS NPKMDFQSWLHMAL  SY*   GFL FRS GQRWQI GV  W G L    F GA GG G  S  VGA G   F   G       W"
assert record.alignments[0].hsps[3].sbjct=="YSHCSGGQIT*VLQILAGKIYSIMXQHSXNPKMDFQSWLHMALXXSY*XXXGFLXFRSXGQRWQIXGVXXWXGXLXXXXFXGAXGGXGXXSXXVGAXGXXXFXXXGXXXXXXXW"
assert record.alignments[0].hsps[3].query_start==273
assert record.alignments[0].hsps[3].query_end==614
assert record.alignments[0].hsps[3].sbjct_start==273
assert record.alignments[0].hsps[3].sbjct_end==614
assert record.alignments[0].hsps[4].query=="PKXXXKXXXXXXKTKXTXXHPPXWXXNXXXLWPXXXXXXLXNXXX*XPXGFATVDXNS*XEESPXKXNXTXPGPCVTSFGNPF*DYXSAXS*CCIFFLPGFGVPRLFVHQNNGCKGEN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCISAFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGILI"
assert record.alignments[0].hsps[4].match=="PK   K      KTK T  HPP W  N   LWP      L N   * P GFATVD NS* EESP K N T PGPCVTSFGNPF*DY SA S*CCIFFLPGFGVPRLFVHQNNGCKGEN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCISAFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGILI"
assert record.alignments[0].hsps[4].sbjct=="PKXXXKXXXXXXKTKXTXXHPPXWXXNXXXLWPXXXXXXLXNXXX*XPXGFATVDXNS*XEESPXKXNXTXPGPCVTSFGNPF*DYXSAXS*CCIFFLPGFGVPRLFVHQNNGCKGEN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCISAFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGILI"
assert record.alignments[0].hsps[4].query_start==615
assert record.alignments[0].hsps[4].query_end==1
assert record.alignments[0].hsps[4].sbjct_start==615
assert record.alignments[0].hsps[4].sbjct_end==1
assert record.alignments[0].hsps[5].query=="PXXXNXXNPXAPTXXXX*PXPPXAPXKXXXXXXPXXLXTPXICHR*PKLLKXRKPX*XQXDXXRAMCNQLWKSILGLXECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD"
assert record.alignments[0].hsps[5].match=="P   N  NP APT    *P PP AP K      P  L TP ICHR*PKLLK RKP * Q D  RAMCNQLWKSILGL EC FIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD"
assert record.alignments[0].hsps[5].sbjct=="PXXXNXXNPXAPTXXXX*PXPPXAPXKXXXXXXPXXLXTPXICHR*PKLLKXRKPX*XQXDXXRAMCNQLWKSILGLXECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD"
assert record.alignments[0].hsps[5].query_start==589
assert record.alignments[0].hsps[5].query_end==2
assert record.alignments[0].hsps[5].sbjct_start==589
assert record.alignments[0].hsps[5].sbjct_end==2
assert record.alignments[0].hsps[6].query=="EXKKAXXXSIXXXQGHV*PALEIHFRIT*VLXHNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALVHLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGAF*"
assert record.alignments[0].hsps[6].match=="E KKA   SI   QGHV*PALEIHFRIT*VL HNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALVHLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGAF*"
assert record.alignments[0].hsps[6].sbjct=="EXKKAXXXSIXXXQGHV*PALEIHFRIT*VLXHNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALVHLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGAF*"
assert record.alignments[0].hsps[6].query_start==440
assert record.alignments[0].hsps[6].query_end==3
assert record.alignments[0].hsps[6].sbjct_start==440
assert record.alignments[0].hsps[6].sbjct_end==3
assert record.alignments[1].hsps[0].query=="SECPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*MH*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[1].hsps[0].match=="S CPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*  *CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[1].hsps[0].sbjct=="SXCPSTHETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHIL*IGPDLD*XX*CKIP*QKRIGYFKPALFFM*PV"
assert record.alignments[1].hsps[0].query_start==3
assert record.alignments[1].hsps[0].query_end==224
assert record.alignments[1].hsps[0].sbjct_start==3
assert record.alignments[1].hsps[0].sbjct_end==224
assert record.alignments[1].hsps[1].query=="YSHCSGGQIT*VLQILAGKIYSIMXQHSXNPKMDFQSWLHMALXXSY*XXXGFLXFRSXGQRWQIXG"
assert record.alignments[1].hsps[1].match=="YSHCSGGQIT*VLQILAGKIYSIM QHS   K  FQSWLHM     +     F   R  GQR Q  G"
assert record.alignments[1].hsps[1].sbjct=="YSHCSGGQIT*VLQILAGKIYSIMKQHSVILKWIFQSWLHMXCKVLFKFKRPFSFTRGLGQRXQTPG"
assert record.alignments[1].hsps[1].query_start==273
assert record.alignments[1].hsps[1].query_end==473
assert record.alignments[1].hsps[1].sbjct_start==273
assert record.alignments[1].hsps[1].sbjct_end==473
assert record.alignments[1].hsps[2].query=="DQNAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKCTNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL*XSTXVILKWISKAGYTWP"
assert record.alignments[1].hsps[2].match=="DQ APPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTK  NAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL* ST       SKAGYT P"
assert record.alignments[1].hsps[2].sbjct=="DQXAPPLMRLFILSTLTGKVGMYAELSKETKKAKTVPSATSSELALTWTKXXNAKSLDKSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVLVDK*PRYSKSWQEKYTAL*NSTQ*S*NGFSKAGYTXP"
assert record.alignments[1].hsps[2].query_start==1
assert record.alignments[1].hsps[2].query_end==399
assert record.alignments[1].hsps[2].sbjct_start==1
assert record.alignments[1].hsps[2].sbjct_end==399
assert record.alignments[1].hsps[3].query=="IRMPLHS*DSSFCPL*QEKWECMXXXXXXXXRPKRCLQPHPLNWP*LGLNALMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPLFWWTNNLGTPNPGRKNIQHYEXALX*S*NGFPKLVTH"
assert record.alignments[1].hsps[3].match=="IR PLHS*DSSFCPL*QEKWECM        RPKRCLQPHPLNWP*LGL  LMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQF PLQPLFWWTNNLGTPNPGRKNIQHYE AL      FPKLVTH"
assert record.alignments[1].hsps[3].sbjct=="IRXPLHS*DSSFCPL*QEKWECMQSSQKKQKRPKRCLQPHPLNWP*LGLXPLMQNPLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFFPLQPLFWWTNNLGTPNPGRKNIQHYETALSNPKMDFPKLVTH"
assert record.alignments[1].hsps[3].query_start==2
assert record.alignments[1].hsps[3].query_end==394
assert record.alignments[1].hsps[3].sbjct_start==2
assert record.alignments[1].hsps[3].sbjct_end==394
assert record.alignments[1].hsps[4].query=="FQEFXSTVANPWG"
assert record.alignments[1].hsps[4].match=="+Q F ST ANPWG"
assert record.alignments[1].hsps[4].sbjct=="YQGFRSTXANPWG"
assert record.alignments[1].hsps[4].query_start==437
assert record.alignments[1].hsps[4].query_end==475
assert record.alignments[1].hsps[4].sbjct_start==437
assert record.alignments[1].hsps[4].sbjct_end==475
assert record.alignments[1].hsps[5].query=="PXICHR*PKLLKXRKPX*XQXDXXRAMCNQLWKSILGLXECXFIMLYIFPARIWST*VICPPEQWL*RRKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*CI*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHSD"
assert record.alignments[1].hsps[5].match=="P +C R*PK L   K         + MCNQLWK    + EC FIMLYIFPARIWST*VICPPEQWL*R+KLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*  *SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEG SD"
assert record.alignments[1].hsps[5].sbjct=="PGVCXR*PKPLVNEKGLLNLNKTLQXMCNQLWKIHFRITECCFIMLYIFPARIWST*VICPPEQWL*RKKLSSGQKLLRRCGKTGYIKNNAGLK*PMRFCQGILH*XX*SKSGPIQRMWLKAPFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGXSD"
assert record.alignments[1].hsps[5].query_start==472
assert record.alignments[1].hsps[5].query_end==2
assert record.alignments[1].hsps[5].sbjct_start==472
assert record.alignments[1].hsps[5].sbjct_end==2
assert record.alignments[1].hsps[6].query=="GHV*PALEIHFRIT*VLXHNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALVHLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGA"
assert record.alignments[1].hsps[6].match=="GHV*PALE  F   *VL HNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFAL  LVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGA"
assert record.alignments[1].hsps[6].sbjct=="GHV*PALENPF*DY*VLFHNAVYFSCQDLEYLGYLSTRTMAVKEKIEQWSEAAEKMR*NRLHKKQCWFEITYALLSRDFALXXLVQVRANSEDVAEGTVLAFFVSFESSAYIPTFPVKVDKMKSLMSGGA"
assert record.alignments[1].hsps[6].query_start==398
assert record.alignments[1].hsps[6].query_end==9
assert record.alignments[1].hsps[6].sbjct_start==398
assert record.alignments[1].hsps[6].sbjct_end==9
assert record.alignments[1].hsps[7].query=="CVTSFGNPF*DYXSAXS*CCIFFLPGFGVPRLFVHQNNGCKGEN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCISAFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGILI"
assert record.alignments[1].hsps[7].match=="CVTSFG       SA S*CCIFFLPGFGVPRLFVHQNNGCKG+N*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCI  FSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRG LI"
assert record.alignments[1].hsps[7].sbjct=="CVTSFGKSILGLLSAVS*CCIFFLPGFGVPRLFVHQNNGCKGKN*AVVRSC*EDAVKQVT*KTMLV*NNLCAFVKGFCIXGFSPSQGQFRGCG*RHRFGLFCFF*ELCIHSHFSCQSGQNEESHEWRGXLI"
assert record.alignments[1].hsps[7].query_start==393
assert record.alignments[1].hsps[7].query_end==1
assert record.alignments[1].hsps[7].sbjct_start==393
assert record.alignments[1].hsps[7].sbjct_end==1
assert record.alignments[2].hsps[0].query=="LCNLFYRIFSAASDHCSIFSFTAIVLV"
assert record.alignments[2].hsps[0].match=="L NL +R+F++  DHCS+  F   V++"
assert record.alignments[2].hsps[0].sbjct=="LTNLIFRLFTSK*DHCSLKQFNNKVVI"
assert record.alignments[2].hsps[0].query_start==211
assert record.alignments[2].hsps[0].query_end==291
assert record.alignments[2].hsps[0].sbjct_start==2
assert record.alignments[2].hsps[0].sbjct_end==82
assert record.alignments[3].hsps[0].query=="IWST*VICPPEQWL*RRKL"
assert record.alignments[3].hsps[0].match=="+W T  +CPP  W *RR+L"
assert record.alignments[3].hsps[0].sbjct=="LWGTAPLCPPVPWA*RRQL"
assert record.alignments[3].hsps[0].query_start==316
assert record.alignments[3].hsps[0].query_end==260
assert record.alignments[3].hsps[0].sbjct_start==226
assert record.alignments[3].hsps[0].sbjct_end==282
assert record.alignments[4].hsps[0].query=="PFWPFLFLLRALHTFPLFLSKWTK"
assert record.alignments[4].hsps[0].match=="P WP + LL  +H+FP  +  W K"
assert record.alignments[4].hsps[0].sbjct=="PEWPSINLLPGMHSFPRVIPCWEK"
assert record.alignments[4].hsps[0].query_start==106
assert record.alignments[4].hsps[0].query_end==35
assert record.alignments[4].hsps[0].sbjct_start==716
assert record.alignments[4].hsps[0].sbjct_end==645
assert record.alignments[5].hsps[0].query=="PFWPFLFLLRALHTFPLFLSKWTK*RVS*VEGHS"
assert record.alignments[5].hsps[0].match=="P W   FLLR LHT P    K+    V+  EG S"
assert record.alignments[5].hsps[0].sbjct=="PGWTLRFLLRGLHTPPKAAEKFPSQGVAEPEGFS"
assert record.alignments[5].hsps[0].query_start==106
assert record.alignments[5].hsps[0].query_end==5
assert record.alignments[5].hsps[0].sbjct_start==176
assert record.alignments[5].hsps[0].sbjct_end==75
assert record.alignments[6].hsps[0].query=="STHETLHFVHFDRKSGNVCRALKRNKKG"
assert record.alignments[6].hsps[0].match=="+TH    F+H  +K  N+C  +  N++G"
assert record.alignments[6].hsps[0].sbjct=="TTHSMTVFLHVKKKLSNICIYMPENREG"
assert record.alignments[6].hsps[0].query_start==15
assert record.alignments[6].hsps[0].query_end==98
assert record.alignments[6].hsps[0].sbjct_start==82
assert record.alignments[6].hsps[0].sbjct_end==165
assert record.alignments[7].hsps[0].query=="HCFLCNLFYRIFSAASDHCSIFSF"
assert record.alignments[7].hsps[0].match=="HC L ++F RIF    +  S+F+F"
assert record.alignments[7].hsps[0].sbjct=="HCHLSSVFCRIFFTLEESLSLFAF"
assert record.alignments[7].hsps[0].query_start==202
assert record.alignments[7].hsps[0].query_end==273
assert record.alignments[7].hsps[0].sbjct_start==443
assert record.alignments[7].hsps[0].sbjct_end==372
assert record.alignments[8].hsps[0].query=="ETLHFVHFDRKSGNVCRALKRNKKGQNGAFSHI"
assert record.alignments[8].hsps[0].match=="E +H+V F +  G +  + ++ ++ QN   SH+"
assert record.alignments[8].hsps[0].sbjct=="EKIHYVLF*KNCGKIMTS*RQRQENQNNLLSHV"
assert record.alignments[8].hsps[0].query_start==24
assert record.alignments[8].hsps[0].query_end==122
assert record.alignments[8].hsps[0].sbjct_start==141
assert record.alignments[8].hsps[0].sbjct_end==43
assert record.alignments[9].hsps[0].query=="PLTKAHRLFQTSIVFYVTCFTASSQQLLTTAQFSPLQPL"
assert record.alignments[9].hsps[0].match=="PL    RLF + +  ++ C T     L  +  F P++PL"
assert record.alignments[9].hsps[0].sbjct=="PLKTQGRLFYSKVSLFLKCRTTLCLFLNVSEGFXPIEPL"
assert record.alignments[9].hsps[0].query_start==167
assert record.alignments[9].hsps[0].query_end==283
assert record.alignments[9].hsps[0].sbjct_start==211
assert record.alignments[9].hsps[0].sbjct_end==95
assert record.alignments[10].hsps[0].query=="CFLCNLFYRIFSAASDHCSIFSFTAI"
assert record.alignments[10].hsps[0].match=="CF+C L  +I+  A +   +  FT++"
assert record.alignments[10].hsps[0].sbjct=="CFICKLVLKIYLRAEERTQLIEFTSL"
assert record.alignments[10].hsps[0].query_start==205
assert record.alignments[10].hsps[0].query_end==282
assert record.alignments[10].hsps[0].sbjct_start==249
assert record.alignments[10].hsps[0].sbjct_end==326
assert record.alignments[11].hsps[0].query=="KSA*VISNQHCFLCNLFYRIFSAASDHCSIFSFTAIVL"
assert record.alignments[11].hsps[0].match=="+SA +I   HC    +   +FS ASD  S    T ++L"
assert record.alignments[11].hsps[0].sbjct=="QSAWIIGMSHCAWAVIHCLLFSTASDEKSAVDLTGVLL"
assert record.alignments[11].hsps[0].query_start==175
assert record.alignments[11].hsps[0].query_end==288
assert record.alignments[11].hsps[0].sbjct_start==119
assert record.alignments[11].hsps[0].sbjct_end==232
assert record.alignments[12].hsps[0].query=="LKAPFWPFLFLLRAL"
assert record.alignments[12].hsps[0].match=="L  PFW FLF L+AL"
assert record.alignments[12].hsps[0].sbjct=="LSVPFWKFLFYLQAL"
assert record.alignments[12].hsps[0].query_start==115
assert record.alignments[12].hsps[0].query_end==71
assert record.alignments[12].hsps[0].sbjct_start==242
assert record.alignments[12].hsps[0].sbjct_end==198
assert record.alignments[13].hsps[0].query=="CCIFFLPGFGVPRLFVHQNNG"
assert record.alignments[13].hsps[0].match=="C IF++P F    LF+H+  G"
assert record.alignments[13].hsps[0].sbjct=="CFIFYVPDFPWSNLFLHRGRG"
assert record.alignments[13].hsps[0].query_start==339
assert record.alignments[13].hsps[0].query_end==277
assert record.alignments[13].hsps[0].sbjct_start==1396
assert record.alignments[13].hsps[0].sbjct_end==1334
assert record.alignments[14].hsps[0].query=="WLKAPFWPFLFLLRALHTFPLFLS"
assert record.alignments[14].hsps[0].match=="WL  PF PFL  L +L   P  LS"
assert record.alignments[14].hsps[0].sbjct=="WLLFPFLPFLPFLPSLPFLPFLLS"
assert record.alignments[14].hsps[0].query_start==118
assert record.alignments[14].hsps[0].query_end==47
assert record.alignments[14].hsps[0].sbjct_start==153
assert record.alignments[14].hsps[0].sbjct_end==224
assert record.alignments[15].hsps[0].query=="WLKAPFWPFLFLLRALHTFPLFLS"
assert record.alignments[15].hsps[0].match=="WL  PF PFL  L +L   P  LS"
assert record.alignments[15].hsps[0].sbjct=="WLLFPFLPFLPFLPSLPFLPFLLS"
assert record.alignments[15].hsps[0].query_start==118
assert record.alignments[15].hsps[0].query_end==47
assert record.alignments[15].hsps[0].sbjct_start==235
assert record.alignments[15].hsps[0].sbjct_end==164
assert record.alignments[16].hsps[0].query=="CFLCNLFYRIFSAASDHCSIFS"
assert record.alignments[16].hsps[0].match=="CF   + +R+F+    HC+ F+"
assert record.alignments[16].hsps[0].sbjct=="CFALTVIWRVFAGCRPHCATFT"
assert record.alignments[16].hsps[0].query_start==205
assert record.alignments[16].hsps[0].query_end==270
assert record.alignments[16].hsps[0].sbjct_start==13
assert record.alignments[16].hsps[0].sbjct_end==78
assert record.alignments[17].hsps[0].query=="VKGFCISAFSPSQGQFRGCG"
assert record.alignments[17].hsps[0].match=="V GFC+++FSP +    G G"
assert record.alignments[17].hsps[0].sbjct=="VDGFCVTSFSPKKDTHPGSG"
assert record.alignments[17].hsps[0].query_start==174
assert record.alignments[17].hsps[0].query_end==115
assert record.alignments[17].hsps[0].sbjct_start==125
assert record.alignments[17].hsps[0].sbjct_end==184
assert record.alignments[18].hsps[0].query=="NQLWKSILGLXECXFIMLYI"
assert record.alignments[18].hsps[0].match=="NQ+WK  L +  C F+ +Y+"
assert record.alignments[18].hsps[0].sbjct=="NQMWKHSLEVCMCVFVYIYV"
assert record.alignments[18].hsps[0].query_start==388
assert record.alignments[18].hsps[0].query_end==329
assert record.alignments[18].hsps[0].sbjct_start==37
assert record.alignments[18].hsps[0].sbjct_end==96
assert record.database_name==['data/sts']
assert record.posted_date==[('Feb 11, 2000  2:37 PM',)]
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.ka_params==[0.318, 0.135, 0.401]
assert record.matrix=='BLOSUM62'
assert record.num_hits==40473548
assert record.num_sequences==87792
assert record.num_extends==487631
assert record.num_good_extends==13175
assert record.num_seqs_better_e==38
assert record.query_length==205
assert record.database_length==10666284
assert record.effective_hsp_length==46
assert record.effective_query_length==158
assert record.effective_database_length==6627852
assert record.effective_search_space==1047200616
assert record.effective_search_space_used==1047200616
assert record.frameshift==('50,', '0.1')
assert record.threshold==13
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(0,0.0)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(52,26.7)

handle = open('Blast/bt057')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|585505|sp|Q08386|MOPB_RHOCA MOLYBDENUM-PTERIN BINDING\nPROTEIN MOPB"
assert record.query_letters==270
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==13
assert record.descriptions[0].title=="gi|585505|sp|Q08386|MOPB_RHOCA MOLYBDENUM-PTERIN BINDING PROTEI..."
assert record.descriptions[0].score==467
assert record.descriptions[0].e==9.9999999999999999e-132
assert record.descriptions[1].title=="gi|585504|sp|Q08385|MOPA_RHOCA MOLYBDENUM-PTERIN BINDING PROTEI..."
assert record.descriptions[1].score==207
assert record.descriptions[1].e==2.0000000000000001e-53
assert record.descriptions[2].title=="gi|585492|sp|P37733|MODA_AZOVI MOLYBDENUM TRANSPORT PROTEIN MODA"
assert record.descriptions[2].score==145
assert record.descriptions[2].e==9.0000000000000002e-35
assert record.descriptions[3].title=="gi|1709070|sp|P46930|MODE_ECOLI MOLYBDENUM TRANSPORT PROTEIN MODE"
assert record.descriptions[3].score==87
assert record.descriptions[3].e==4.9999999999999999e-17
assert record.descriptions[4].title=="gi|1709071|sp|P45324|MODE_HAEIN MOLYBDENUM TRANSPORT PROTEIN MO..."
assert record.descriptions[4].score==54
assert record.descriptions[4].e==1.9999999999999999e-07
assert record.descriptions[5].title=="gi|585502|sp|P04952|MOP1_CLOPA MOLYBDENUM-PTERIN BINDING PROTEIN I"
assert record.descriptions[5].score==53
assert record.descriptions[5].e==5.9999999999999997e-07
assert record.descriptions[6].title=="gi|127241|sp|P08854|MOP2_CLOPA MOLYBDENUM-PTERIN BINDING PROTEI..."
assert record.descriptions[6].score==52
assert record.descriptions[6].e==0.000001
assert record.descriptions[7].title=="gi|585503|sp|P38366|MOP3_CLOPA MOLYBDENUM-PTERIN BINDING PROTEI..."
assert record.descriptions[7].score==51
assert record.descriptions[7].e==0.000003
assert record.descriptions[8].title=="gi|1170996|sp|P45183|MOP_HAEIN PROBABLE MOLYBDENUM-PTERIN BINDI..."
assert record.descriptions[8].score==46
assert record.descriptions[8].e==0.000050
assert record.descriptions[9].title=="gi|1709069|sp|P09833|MODC_ECOLI MOLYBDENUM TRANSPORT ATP-BINDIN..."
assert record.descriptions[9].score==38
assert record.descriptions[9].e==0.021000
assert record.descriptions[10].title=="gi|585500|sp|P37732|MODD_AZOVI MOLYBDENUM TRANSPORT ATP-BINDING..."
assert record.descriptions[10].score==33
assert record.descriptions[10].e==0.530000
assert record.descriptions[11].title=="gi|2507168|sp|P08838|PT1_BACSU PHOSPHOENOLPYRUVATE-PROTEIN PHOS..."
assert record.descriptions[11].score==30
assert record.descriptions[11].e==4.600000
assert record.descriptions[12].title=="gi|729786|sp|Q05355|HYDL_STRHA PUTATIVE POLYKETIDE HYDROXYLASE"
assert record.descriptions[12].score==29
assert record.descriptions[12].e==7.900000
assert len(record.alignments)==13
assert record.alignments[0].title==">gi|585505|sp|Q08386|MOPB_RHOCA MOLYBDENUM-PTERIN BINDING PROTEIN MOPB"
assert record.alignments[0].length==270
assert record.alignments[1].title==">gi|585504|sp|Q08385|MOPA_RHOCA MOLYBDENUM-PTERIN BINDING PROTEIN MOPA"
assert record.alignments[1].length==265
assert record.alignments[2].title==">gi|585492|sp|P37733|MODA_AZOVI MOLYBDENUM TRANSPORT PROTEIN MODA"
assert record.alignments[2].length==270
assert record.alignments[3].title==">gi|1709070|sp|P46930|MODE_ECOLI MOLYBDENUM TRANSPORT PROTEIN MODE"
assert record.alignments[3].length==262
assert record.alignments[4].title==">gi|1709071|sp|P45324|MODE_HAEIN MOLYBDENUM TRANSPORT PROTEIN MODE HOMOLOG"
assert record.alignments[4].length==255
assert record.alignments[5].title==">gi|585502|sp|P04952|MOP1_CLOPA MOLYBDENUM-PTERIN BINDING PROTEIN I"
assert record.alignments[5].length==68
assert record.alignments[6].title==">gi|127241|sp|P08854|MOP2_CLOPA MOLYBDENUM-PTERIN BINDING PROTEIN II"
assert record.alignments[6].length==68
assert record.alignments[7].title==">gi|585503|sp|P38366|MOP3_CLOPA MOLYBDENUM-PTERIN BINDING PROTEIN III"
assert record.alignments[7].length==68
assert record.alignments[8].title==">gi|1170996|sp|P45183|MOP_HAEIN PROBABLE MOLYBDENUM-PTERIN BINDING PROTEIN"
assert record.alignments[8].length==69
assert record.alignments[9].title==">gi|1709069|sp|P09833|MODC_ECOLI MOLYBDENUM TRANSPORT ATP-BINDING PROTEIN MODC"
assert record.alignments[9].length==352
assert record.alignments[10].title==">gi|585500|sp|P37732|MODD_AZOVI MOLYBDENUM TRANSPORT ATP-BINDING PROTEIN MODD"
assert record.alignments[10].length==380
assert record.alignments[11].title==">gi|2507168|sp|P08838|PT1_BACSU PHOSPHOENOLPYRUVATE-PROTEIN PHOSPHOTRANSFERASE (PHOSPHOTRANSFERASE SYSTEM, ENZYME I)"
assert record.alignments[11].length==570
assert record.alignments[12].title==">gi|729786|sp|Q05355|HYDL_STRHA PUTATIVE POLYKETIDE HYDROXYLASE"
assert record.alignments[12].length==555
assert record.alignments[0].hsps[0].score==1189
assert record.alignments[0].hsps[0].bits==467
assert record.alignments[0].hsps[0].expect==1e-131
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==521
assert record.alignments[1].hsps[0].bits==207
assert record.alignments[1].hsps[0].expect==2e-53
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==362
assert record.alignments[2].hsps[0].bits==145
assert record.alignments[2].hsps[0].expect==9e-35
assert record.alignments[2].hsps[1].score==98
assert record.alignments[2].hsps[1].bits==42.6
assert record.alignments[2].hsps[1].expect==8e-04
assert len(record.alignments[2].hsps)==2
assert record.alignments[3].hsps[0].score==211
assert record.alignments[3].hsps[0].bits==86.6
assert record.alignments[3].hsps[0].expect==5e-17
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==128
assert record.alignments[4].hsps[0].bits==54.3
assert record.alignments[4].hsps[0].expect==2e-07
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==125
assert record.alignments[5].hsps[0].bits==53.1
assert record.alignments[5].hsps[0].expect==6e-07
assert record.alignments[5].hsps[1].score==84
assert record.alignments[5].hsps[1].bits==37.1
assert record.alignments[5].hsps[1].expect==0.036
assert len(record.alignments[5].hsps)==2
assert record.alignments[6].hsps[0].score==123
assert record.alignments[6].hsps[0].bits==52.3
assert record.alignments[6].hsps[0].expect==1e-06
assert record.alignments[6].hsps[1].score==86
assert record.alignments[6].hsps[1].bits==37.9
assert record.alignments[6].hsps[1].expect==0.021
assert len(record.alignments[6].hsps)==2
assert record.alignments[7].hsps[0].score==119
assert record.alignments[7].hsps[0].bits==50.8
assert record.alignments[7].hsps[0].expect==3e-06
assert record.alignments[7].hsps[1].score==83
assert record.alignments[7].hsps[1].bits==36.7
assert record.alignments[7].hsps[1].expect==0.047
assert len(record.alignments[7].hsps)==2
assert record.alignments[8].hsps[0].score==108
assert record.alignments[8].hsps[0].bits==46.5
assert record.alignments[8].hsps[0].expect==5e-05
assert len(record.alignments[8].hsps)==1
assert record.alignments[9].hsps[0].score==86
assert record.alignments[9].hsps[0].bits==37.9
assert record.alignments[9].hsps[0].expect==0.021
assert len(record.alignments[9].hsps)==1
assert record.alignments[10].hsps[0].score==74
assert record.alignments[10].hsps[0].bits==33.2
assert record.alignments[10].hsps[0].expect==0.53
assert len(record.alignments[10].hsps)==1
assert record.alignments[11].hsps[0].score==66
assert record.alignments[11].hsps[0].bits==30.1
assert record.alignments[11].hsps[0].expect==4.6
assert len(record.alignments[11].hsps)==1
assert record.alignments[12].hsps[0].score==64
assert record.alignments[12].hsps[0].bits==29.3
assert record.alignments[12].hsps[0].expect==7.9
assert len(record.alignments[12].hsps)==1
assert record.alignments[0].hsps[0].identities==(247, 270)
assert record.alignments[0].hsps[0].positives==(247, 270)
assert record.alignments[1].hsps[0].identities==(123, 259)
assert record.alignments[1].hsps[0].positives==(155, 259)
assert record.alignments[1].hsps[0].gaps==(13, 259)
assert record.alignments[2].hsps[0].identities==(93, 253)
assert record.alignments[2].hsps[0].positives==(132, 253)
assert record.alignments[2].hsps[0].gaps==(8, 253)
assert record.alignments[2].hsps[1].identities==(33, 99)
assert record.alignments[2].hsps[1].positives==(47, 99)
assert record.alignments[2].hsps[1].gaps==(7, 99)
assert record.alignments[3].hsps[0].identities==(76, 247)
assert record.alignments[3].hsps[0].positives==(114, 247)
assert record.alignments[3].hsps[0].gaps==(17, 247)
assert record.alignments[4].hsps[0].identities==(46, 170)
assert record.alignments[4].hsps[0].positives==(76, 170)
assert record.alignments[4].hsps[0].gaps==(3, 170)
assert record.alignments[5].hsps[0].identities==(25, 64)
assert record.alignments[5].hsps[0].positives==(43, 64)
assert record.alignments[5].hsps[1].identities==(19, 63)
assert record.alignments[5].hsps[1].positives==(36, 63)
assert record.alignments[6].hsps[0].identities==(24, 64)
assert record.alignments[6].hsps[0].positives==(43, 64)
assert record.alignments[6].hsps[1].identities==(21, 63)
assert record.alignments[6].hsps[1].positives==(36, 63)
assert record.alignments[7].hsps[0].identities==(24, 64)
assert record.alignments[7].hsps[0].positives==(43, 64)
assert record.alignments[7].hsps[1].identities==(20, 63)
assert record.alignments[7].hsps[1].positives==(37, 63)
assert record.alignments[8].hsps[0].identities==(19, 67)
assert record.alignments[8].hsps[0].positives==(46, 67)
assert record.alignments[9].hsps[0].identities==(23, 62)
assert record.alignments[9].hsps[0].positives==(37, 62)
assert record.alignments[9].hsps[0].gaps==(1, 62)
assert record.alignments[10].hsps[0].identities==(41, 143)
assert record.alignments[10].hsps[0].positives==(62, 143)
assert record.alignments[10].hsps[0].gaps==(12, 143)
assert record.alignments[11].hsps[0].identities==(32, 141)
assert record.alignments[11].hsps[0].positives==(61, 141)
assert record.alignments[11].hsps[0].gaps==(6, 141)
assert record.alignments[12].hsps[0].identities==(21, 62)
assert record.alignments[12].hsps[0].positives==(29, 62)
assert record.alignments[0].hsps[0].query=="MAATKQGGGDDGRCARGVVLERTGARMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAYGAIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILALP"
assert record.alignments[0].hsps[0].match=="MAATKQGGGDDGRCARGVVLERTGARMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNN                       AGEKLIAAYGAIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILALP"
assert record.alignments[0].hsps[0].sbjct=="MAATKQGGGDDGRCARGVVLERTGARMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNLLAAPVVTAAPGGKAGGGAVLTPAGEKLIAAYGAIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILALP"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==270
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==270
assert record.alignments[1].hsps[0].query=="LERTGA-RMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAYGAIEAGVAKLL-------SSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILALP"
assert record.alignments[1].hsps[0].match=="L+R GA R+G +R+ LL AI R G+I+ AAREVGLSYK AWD V  +NN                       AG+ LIA +G +E  + K L       S+ EK+LN      L  L++RTS RN   C V  V    V A+V + L +G  LTAVIT RSA EM LAPG EV AL+K++FV+LA  G P R+S  NR+ G V  R D P+++E+ LDLG  K+ITA IT  SA+ L L PGV  TAL K+SHVILA+P"
assert record.alignments[1].hsps[0].sbjct=="LQRAGAPRVGGDRIRLLEAIARHGTIAGAAREVGLSYKTAWDAVGTLNNLFEQPLVEAAPGGRTGGNARVTEAGQALIAGFGLLEGALTKALGVLEGGVSAPEKALN-----TLWSLTMRTSNRNTLRCTVTRVTLGAVNAEVELALTDGHSLTAVITERSATEMGLAPGVEVFALIKASFVMLAAGGDPGRISACNRLTGIVAARTDGPVNTEIILDLGNCKSITAVITHTSADALGLAPGVPATALFKASHVILAMP"
assert record.alignments[1].hsps[0].query_start==20
assert record.alignments[1].hsps[0].query_end==270
assert record.alignments[1].hsps[0].sbjct_start==12
assert record.alignments[1].hsps[0].sbjct_end==265
assert record.alignments[2].hsps[0].query=="GARMGAERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAYGAIEAGVAKLLSSFEKSLNLDPA-------EVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILAL"
assert record.alignments[2].hsps[0].match=="G  +   R+ LL AI R GSI+ AA+ V LSYKAAWD +  MNN                        G +++A Y A+E      L    + LN            ++  +S++TSARN +A  V  +    V  +VR+RL    ++ AVIT  SA  + LA G EV ALVKS+ V+L       +L+ RN++ G VI+  + P+++EVTL L  G+++T  +T DS + L L PGV   A  KSS VILA+"
assert record.alignments[2].hsps[0].sbjct=="GTALSDTRIRLLEAIEREGSINRAAKVVPLSYKAAWDAIDTMNNLAPEPLVVRVAGGRQGGGTQLTDYGRRIVAMYRALEIEYQSALDRLSERLNEVTGGDIQAFQRLMHSMSMKTSARNQFAGIVTGLRVGGVDYEVRIRLDAENEIAAVITKASAENLELAIGKEVFALVKSSSVMLT-TEPSLKLTARNQLWGEVIDIHEGPVNNEVTLALPSGRSVTCVVTADSCKALGLAPGVAACAFFKSSSVILAV"
assert record.alignments[2].hsps[0].query_start==24
assert record.alignments[2].hsps[0].query_end==269
assert record.alignments[2].hsps[0].sbjct_start==17
assert record.alignments[2].hsps[0].sbjct_end==268
assert record.alignments[2].hsps[1].query=="AIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAG"
assert record.alignments[2].hsps[1].match=="AI   V  L+ S    L  +P       SL+ +ARN    +V  +    V  +V + L  G+ +T V+TA S   + LAPG    A  KS+ V+LA  G"
assert record.alignments[2].hsps[1].sbjct=="AIGKEVFALVKSSSVMLTTEP-------SLKLTARNQLWGEVIDIHEGPVNNEVTLALPSGRSVTCVVTADSCKALGLAPGVAACAFFKSSSVILAVYG"
assert record.alignments[2].hsps[1].query_start==101
assert record.alignments[2].hsps[1].query_end==199
assert record.alignments[2].hsps[1].sbjct_start==179
assert record.alignments[2].hsps[1].sbjct_end==270
assert record.alignments[3].hsps[0].query=="RVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAY---GAIEAGVAKLLSSFEK-SLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQD-LTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVR----NRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILA"
assert record.alignments[3].hsps[0].match=="R++LL  I  +GSIS  A++ G+SYK+AWD +  MN                         G++LI  Y     I+     +LS  +   LN   A + R  SL+TSARN W   + +   DDV   V + L +G+  L   ITA+S A + L  G EVL L+K+ +V     G+ +  +V     N++ G +          EV + L  G+T+ AT+  +  E   L  G   TA   +  VI+A"
assert record.alignments[3].hsps[0].sbjct=="RISLLKHIALSGSISQGAKDAGISYKSAWDAINEMNQLSEHILVERATGGKGGGGAVLTRYGQRLIQLYDLLAQIQQKAFDVLSDDDALPLNSLLAAISR-FSLQTSARNQWFGTITARDHDDVQQHVDVLLADGKTRLKVAITAQSGARLGLDEGKEVLILLKAPWV-----GITQDEAVAQNADNQLPGIISHIERGAEQCEVLMALPDGQTLCATVPVN--EATSLQQGQNVTAYFNADSVIIA"
assert record.alignments[3].hsps[0].query_start==31
assert record.alignments[3].hsps[0].query_end==268
assert record.alignments[3].hsps[0].sbjct_start==21
assert record.alignments[3].hsps[0].sbjct_end==259
assert record.alignments[4].hsps[0].query=="ERVALLAAIGRTGSISAAAREVGLSYKAAWDGVQAMNNXXXXXXXXXXXXXXXXXXXXXXXAGEKLIAAYGAIEAGVAKLLSSF-EKSLNLDP-AEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRL-GEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLA"
assert record.alignments[4].hsps[0].match=="+RV LL  I + GSI+ AA+   +SYK+AWD ++AMN                          E+L+  Y  +E           ++S+ LD         SL++SARN +  +V      D    V + + G    L   IT +S+A ++L    EV+ + K+ +V ++"
assert record.alignments[4].hsps[0].sbjct=="KRVRLLKEIQQCGSINQAAKNAKVSYKSAWDHLEAMNKISPRPLLERNTGGKNGGGTALTTYAERLLQLYDLLERTQEHAFHILQDESVPLDSLLTATARFSLQSSARNQFFGRVAQQRIIDSRCVVDVNVQGLPTPLQVSITTKSSARLKLITEKEVMLMFKAPWVKIS"
assert record.alignments[4].hsps[0].query_start==30
assert record.alignments[4].hsps[0].query_end==196
assert record.alignments[4].hsps[0].sbjct_start==21
assert record.alignments[4].hsps[0].sbjct_end==190
assert record.alignments[5].hsps[0].query=="LSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVIL"
assert record.alignments[5].hsps[0].match=="+S RN+++G+V+      +++EV L++ GG  IT+ I+ DS E L +  G E TA+IKS+ V++"
assert record.alignments[5].hsps[0].sbjct=="ISARNQLKGKVVGLKKGVITAEVVLEIAGGNKITSIISLDSVEELGVKEGAELTAVIKSTDVMI"
assert record.alignments[5].hsps[0].query_start==204
assert record.alignments[5].hsps[0].query_end==267
assert record.alignments[5].hsps[0].sbjct_start==3
assert record.alignments[5].hsps[0].sbjct_end==66
assert record.alignments[5].hsps[1].query=="SARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLL"
assert record.alignments[5].hsps[1].match=="SARN    KV  +    + A+V + +  G  +T++I+  S  E+ +  G+E+ A++KS  V++"
assert record.alignments[5].hsps[1].sbjct=="SARNQLKGKVVGLKKGVITAEVVLEIAGGNKITSIISLDSVEELGVKEGAELTAVIKSTDVMI"
assert record.alignments[5].hsps[1].query_start==133
assert record.alignments[5].hsps[1].query_end==195
assert record.alignments[5].hsps[1].sbjct_start==4
assert record.alignments[5].hsps[1].sbjct_end==66
assert record.alignments[6].hsps[0].query=="LSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVIL"
assert record.alignments[6].hsps[0].match=="+S RN+++G+V+      +++EV L++ GG  IT+ I+ DS E L +  G E TA++KS+ V++"
assert record.alignments[6].hsps[0].sbjct=="ISARNQLKGKVVGLKKGVVTAEVVLEIAGGNKITSIISLDSVEELGVKEGAELTAVVKSTDVMI"
assert record.alignments[6].hsps[0].query_start==204
assert record.alignments[6].hsps[0].query_end==267
assert record.alignments[6].hsps[0].sbjct_start==3
assert record.alignments[6].hsps[0].sbjct_end==66
assert record.alignments[6].hsps[1].query=="SARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLL"
assert record.alignments[6].hsps[1].match=="SARN    KV  +    V A+V + +  G  +T++I+  S  E+ +  G+E+ A+VKS  V++"
assert record.alignments[6].hsps[1].sbjct=="SARNQLKGKVVGLKKGVVTAEVVLEIAGGNKITSIISLDSVEELGVKEGAELTAVVKSTDVMI"
assert record.alignments[6].hsps[1].query_start==133
assert record.alignments[6].hsps[1].query_end==195
assert record.alignments[6].hsps[1].sbjct_start==4
assert record.alignments[6].hsps[1].sbjct_end==66
assert record.alignments[7].hsps[0].query=="LSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVIL"
assert record.alignments[7].hsps[0].match=="+S RN+++G+V+      +++EV L++ GG  +T+ I+ DS E L +  G E TA+IKS+ V++"
assert record.alignments[7].hsps[0].sbjct=="ISARNQLKGKVVAVKKGLVTAEVVLEIAGGDKVTSIISLDSIEDLGVKEGTELTAVIKSTDVMI"
assert record.alignments[7].hsps[0].query_start==204
assert record.alignments[7].hsps[0].query_end==267
assert record.alignments[7].hsps[0].sbjct_start==3
assert record.alignments[7].hsps[0].sbjct_end==66
assert record.alignments[7].hsps[1].query=="SARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLL"
assert record.alignments[7].hsps[1].match=="SARN    KV +V    V A+V + +  G  +T++I+  S  ++ +  G+E+ A++KS  V++"
assert record.alignments[7].hsps[1].sbjct=="SARNQLKGKVVAVKKGLVTAEVVLEIAGGDKVTSIISLDSIEDLGVKEGTELTAVIKSTDVMI"
assert record.alignments[7].hsps[1].query_start==133
assert record.alignments[7].hsps[1].query_end==195
assert record.alignments[7].hsps[1].sbjct_start==4
assert record.alignments[7].hsps[1].sbjct_end==66
assert record.alignments[8].hsps[0].query=="RLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKSSHVILAL"
assert record.alignments[8].hsps[0].match=="++S RN+++G+V+   +  +++ V +D+GGG  +++T++  + + L+L  G E  A+IK++ V++ +"
assert record.alignments[8].hsps[0].sbjct=="KISARNQLKGKVVSIENGSVNAIVHIDIGGGNVLSSTVSLAAVKELNLEVGKEAYAIIKATSVMVGV"
assert record.alignments[8].hsps[0].query_start==203
assert record.alignments[8].hsps[0].query_end==269
assert record.alignments[8].hsps[0].sbjct_start==2
assert record.alignments[8].hsps[0].sbjct_end==68
assert record.alignments[9].hsps[0].query=="PERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTITATITRDSAEMLDLHPGVETTALIKS"
assert record.alignments[9].hsps[0].match=="P++ S+RN +R +V+   D     EV L++ GGKT+ A I+  + + L + PG+   A IKS"
assert record.alignments[9].hsps[0].sbjct=="PQQTSIRNVLRAKVVNSYDDNGQVEVELEV-GGKTLWARISPWARDELAIKPGLWLYAQIKS"
assert record.alignments[9].hsps[0].query_start==201
assert record.alignments[9].hsps[0].query_end==262
assert record.alignments[9].hsps[0].sbjct_start==287
assert record.alignments[9].hsps[0].sbjct_end==347
assert record.alignments[10].hsps[0].query=="EVLRGLSLRTSARNAWACKVWSVAA--DDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGG-GKTITATITRDSAEMLDLHPGVETTALIKS"
assert record.alignments[10].hsps[0].match=="+++  L L T+        + SV A  DD     R+    G    AV+ AR       APG  +   V +  V LA + + E  S+ N +   V E ++A   + V + L   G  + A ITR S + L + PG    A IK+"
assert record.alignments[10].hsps[0].sbjct=="DIMARLDLPTAFHEDAGVVIESVVAEHDDHYHLTRLAFPGG----AVLVARRPE----APGQRLRLRVHARDVSLANSRI-EDSSITNVLPATVREVVEADTPAHVLVRLEAEGTPLIARITRRSCDQLGIAPGRRMWAQIKA"
assert record.alignments[10].hsps[0].query_start==123
assert record.alignments[10].hsps[0].query_end==262
assert record.alignments[10].hsps[0].sbjct_start==242
assert record.alignments[10].hsps[0].sbjct_end==375
assert record.alignments[11].hsps[0].query=="AAYGAIEAGVAKLLSSFEKSLNLDP-AEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQDLTAVITARSAAEMRLAPGSEVLALVKSNFVLLAGAGVPERLSVRNRVRGRVIERIDAPLSSEVTLDLGGGKTI"
assert record.alignments[11].hsps[0].match=="AA G I+ GV  ++      + +DP AE ++    + +A  A   + W+   ++       + G   +L A I      +  L  G E + L ++ F+ +    +P      +  +  V+ER++       TLD+GG K +"
assert record.alignments[11].hsps[0].sbjct=="AATGTIQNGVTVIVDGINGDVIIDPSAETVKEYEEKHNAYLAQKAE-WAKLVNEPTVS---KDGHHVELAANIGTPDDVKGVLENGGEAVGLYRTEFLYMGRDQLPTEDEQFDAYK-TVLERMEGKSVVVRTLDIGGDKEL"
assert record.alignments[11].hsps[0].query_start==97
assert record.alignments[11].hsps[0].query_end==236
assert record.alignments[11].hsps[0].sbjct_start==207
assert record.alignments[11].hsps[0].sbjct_end==342
assert record.alignments[12].hsps[0].query=="AIEAGVAKLLSSFEKSLNLDPAEVLRGLSLRTSARNAWACKVWSVAADDVAAQVRMRLGEGQ"
assert record.alignments[12].hsps[0].match=="A+E G     S+  +S   DPA V   +  R S  +      + VAAD   + VR +LG GQ"
assert record.alignments[12].hsps[0].sbjct=="AVELGGEIRFSTELQSFEQDPAGVTAVIKSRRSGEHTTVRADYLVAADGPRSPVREQLGIGQ"
assert record.alignments[12].hsps[0].query_start==101
assert record.alignments[12].hsps[0].query_end==162
assert record.alignments[12].hsps[0].sbjct_start==136
assert record.alignments[12].hsps[0].sbjct_end==197
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.ka_params==[0.316, 0.131, 0.361]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==12068104
assert record.num_sequences==82258
assert record.num_extends==396723
assert record.num_good_extends==1066
assert record.num_seqs_better_e==13
assert record.hsps_no_gap==10
assert record.hsps_prelim_gapped==3
assert record.hsps_gapped==18
assert record.query_length==270
assert record.database_length==29652561
assert record.effective_hsp_length==56
assert record.effective_query_length==214
assert record.effective_database_length==25046113
assert record.effective_search_space==5359868182
assert record.effective_search_space_used==5359868182
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.6)
assert record.blast_cutoff==(64,29.3)

handle = open('Blast/bt058')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|730725|sp|Q05362|SCHB_STRHA SCHB PROTEIN"
assert record.query_letters==138
assert record.database=="data/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==0
assert len(record.alignments)==0
assert record.database_name==['data/swissprot']
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.posted_date==[('Feb 2, 2000  9:39 AM',)]
assert record.ka_params==[0.319, 0.139, 0.415]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==8952568
assert record.num_sequences==82258
assert record.num_extends==387403
assert record.num_good_extends==727
assert record.num_seqs_better_e==23
assert record.hsps_no_gap==13
assert record.hsps_prelim_gapped==10
assert record.hsps_gapped==23
assert record.query_length==138
assert record.database_length==29652561
assert record.effective_hsp_length==47
assert record.effective_query_length==91
assert record.effective_database_length==25786435
assert record.effective_search_space==2346565585
assert record.effective_search_space_used==2346565585
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.4)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(61,28.2)

handle = open('Blast/bt059')
record = parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.11'
assert record.date=="Jan-20-2000"
assert record.reference==reference
assert record.query=="gi|129628|sp|P07175|PARA_AGRTU PARA PROTEIN"
assert record.query_letters==222
assert record.database=="/taiyang/jchang/parsers/swissprot"
assert record.database_sequences==82258
assert record.database_letters==29652561
assert len(record.descriptions)==107
assert record.descriptions[0].title=="gi|129628|sp|P07175|PARA_AGRTU PARA PROTEIN"
assert record.descriptions[0].score==427
assert record.descriptions[0].e==9.9999999999999998e-121
assert record.descriptions[1].title=="gi|138456|sp|P06665|VIC1_AGRT6 VIRC1 PROTEIN"
assert record.descriptions[1].score==53
assert record.descriptions[1].e==2.9999999999999999e-07
assert record.descriptions[2].title=="gi|138454|sp|P13459|VIC1_AGRRA VIRC1 PROTEIN"
assert record.descriptions[2].score==52
assert record.descriptions[2].e==6.9999999999999997e-07
assert record.descriptions[3].title=="gi|138455|sp|P07165|VIC1_AGRT5 VIRC1 PROTEIN"
assert record.descriptions[3].score==49
assert record.descriptions[3].e==0.000006
assert record.descriptions[4].title=="gi|2497977|sp|P72190|YCAB_PSEFR HYPOTHETICAL 30.2 KD PROTEIN IN..."
assert record.descriptions[4].score==42
assert record.descriptions[4].e==0.001000
assert record.descriptions[5].title=="gi|586852|sp|P37522|SOJ_BACSU SOJ PROTEIN"
assert record.descriptions[5].score==42
assert record.descriptions[5].e==0.001000
assert record.descriptions[6].title=="gi|132348|sp|P05682|REPA_AGRRA POSSIBLE REPLICATION PROTEIN A"
assert record.descriptions[6].score==41
assert record.descriptions[6].e==0.001000
assert record.descriptions[7].title=="gi|120545|sp|P28373|CHLL_SYNY3 PROTOCHLOROPHYLLIDE REDUCTASE IR..."
assert record.descriptions[7].score==41
assert record.descriptions[7].e==0.001000
assert record.descriptions[8].title=="gi|2496226|sp|Q60283|YZ24_METJA HYPOTHETICAL PROTEIN MJECL24"
assert record.descriptions[8].score==40
assert record.descriptions[8].e==0.003000
assert record.descriptions[9].title=="gi|3024144|sp|Q55900|MIND_SYNY3 SEPTUM SITE-DETERMINING PROTEIN..."
assert record.descriptions[9].score==40
assert record.descriptions[9].e==0.003000
assert record.descriptions[10].title=="gi|1709100|sp|P53382|MRP_MYCLE MRP PROTEIN HOMOLOG"
assert record.descriptions[10].score==40
assert record.descriptions[10].e==0.003000
assert record.descriptions[11].title=="gi|462615|sp|P21590|MRP_ECOLI MRP PROTEIN"
assert record.descriptions[11].score==40
assert record.descriptions[11].e==0.004000
assert record.descriptions[12].title=="gi|3024135|sp|P56346|MIND_CHLVU PUTATIVE SEPTUM SITE-DETERMININ..."
assert record.descriptions[12].score==39
assert record.descriptions[12].e==0.005000
assert record.descriptions[13].title=="gi|1709101|sp|P53383|MRP_SYNY3 MRP PROTEIN HOMOLOG"
assert record.descriptions[13].score==39
assert record.descriptions[13].e==0.009000
assert record.descriptions[14].title=="gi|127097|sp|P18197|MIND_ECOLI SEPTUM SITE-DETERMINING PROTEIN ..."
assert record.descriptions[14].score==39
assert record.descriptions[14].e==0.009000
assert record.descriptions[15].title=="gi|3183081|sp|O24999|MRP_HELPY MRP PROTEIN HOMOLOG"
assert record.descriptions[15].score==38
assert record.descriptions[15].e==0.012000
assert record.descriptions[16].title=="gi|1345782|sp|P48110|CHLL_CYAPA PROTOCHLOROPHYLLIDE REDUCTASE I..."
assert record.descriptions[16].score==38
assert record.descriptions[16].e==0.012000
assert record.descriptions[17].title=="gi|400260|sp|Q01464|MIND_BACSU SEPTUM SITE-DETERMINING PROTEIN ..."
assert record.descriptions[17].score==38
assert record.descriptions[17].e==0.012000
assert record.descriptions[18].title=="gi|3913244|sp|O47041|CHLL_PICAB PROTOCHLOROPHYLLIDE REDUCTASE I..."
assert record.descriptions[18].score==38
assert record.descriptions[18].e==0.016000
assert record.descriptions[19].title=="gi|1168937|sp|P41645|CHLL_PINTH PROTOCHLOROPHYLLIDE REDUCTASE I..."
assert record.descriptions[19].score==38
assert record.descriptions[19].e==0.016000
assert record.descriptions[20].title=="gi|120543|sp|P26181|CHLL_PINCO PROTOCHLOROPHYLLIDE REDUCTASE IR..."
assert record.descriptions[20].score==38
assert record.descriptions[20].e==0.016000
assert record.descriptions[21].title=="gi|401555|sp|P31856|YGI1_PSEPU HYPOTHETICAL 28.9 KD PROTEIN IN ..."
assert record.descriptions[21].score==38
assert record.descriptions[21].e==0.016000
assert record.descriptions[22].title=="gi|120544|sp|Q00237|CHLL_PLEBO PROTOCHLOROPHYLLIDE REDUCTASE IR..."
assert record.descriptions[22].score==38
assert record.descriptions[22].e==0.016000
assert record.descriptions[23].title=="gi|6225723|sp|O33225|MRP_MYCTU MRP PROTEIN HOMOLOG"
assert record.descriptions[23].score==38
assert record.descriptions[23].e==0.021000
assert record.descriptions[24].title=="gi|2496603|sp|P55393|Y4CK_RHISN PUTATIVE REPLICATION PROTEIN A"
assert record.descriptions[24].score==38
assert record.descriptions[24].e==0.021000
assert record.descriptions[25].title=="gi|120542|sp|P06267|CHLL_MARPO PROTOCHLOROPHYLLIDE REDUCTASE IR..."
assert record.descriptions[25].score==38
assert record.descriptions[25].e==0.021000
assert record.descriptions[26].title=="gi|1705820|sp|P54207|CHLL_SYNP7 PROTOCHLOROPHYLLIDE REDUCTASE I..."
assert record.descriptions[26].score==38
assert record.descriptions[26].e==0.021000
assert record.descriptions[27].title=="gi|124472|sp|P07673|INC2_ECOLI INCC PROTEIN"
assert record.descriptions[27].score==38
assert record.descriptions[27].e==0.021000
assert record.descriptions[28].title=="gi|3023485|sp|P56291|CHLL_CHLVU PROTOCHLOROPHYLLIDE REDUCTASE I..."
assert record.descriptions[28].score==37
assert record.descriptions[28].e==0.027000
assert record.descriptions[29].title=="gi|6016572|sp|O78436|MIND_GUITH PUTATIVE SEPTUM SITE-DETERMININ..."
assert record.descriptions[29].score==36
assert record.descriptions[29].e==0.047000
assert record.descriptions[30].title=="gi|6225722|sp|O66946|MRP_AQUAE MRP PROTEIN HOMOLOG"
assert record.descriptions[30].score==36
assert record.descriptions[30].e==0.047000
assert record.descriptions[31].title=="gi|1723296|sp|P50863|MRP_BACSU MRP PROTEIN HOMOLOG"
assert record.descriptions[31].score==36
assert record.descriptions[31].e==0.047000
assert record.descriptions[32].title=="gi|128211|sp|P08625|NIF2_METTL NITROGENASE IRON PROTEIN 2 (NITR..."
assert record.descriptions[32].score==36
assert record.descriptions[32].e==0.047000
assert record.descriptions[33].title=="gi|120541|sp|Q00469|CHLL_CHLRE PROTOCHLOROPHYLLIDE REDUCTASE IR..."
assert record.descriptions[33].score==36
assert record.descriptions[33].e==0.047000
assert record.descriptions[34].title=="gi|585623|sp|Q07733|OPPD_LACLA OLIGOPEPTIDE TRANSPORT ATP-BINDI..."
assert record.descriptions[34].score==36
assert record.descriptions[34].e==0.047000
assert record.descriptions[35].title=="gi|1705819|sp|P51187|CHLL_PORPU PROTOCHLOROPHYLLIDE REDUCTASE I..."
assert record.descriptions[35].score==36
assert record.descriptions[35].e==0.061000
assert record.descriptions[36].title=="gi|6226409|sp|O58667|Y949_PYRHO HYPOTHETICAL PROTEIN PH0949"
assert record.descriptions[36].score==36
assert record.descriptions[36].e==0.081000
assert record.descriptions[37].title=="gi|6226573|sp|P31897|COOC_RHORU CARBON MONOXIDE DEHYDROGENASE A..."
assert record.descriptions[37].score==36
assert record.descriptions[37].e==0.081000
assert record.descriptions[38].title=="gi|544021|sp|P36439|CHLL_POLAC PROTOCHLOROPHYLLIDE REDUCTASE IR..."
assert record.descriptions[38].score==36
assert record.descriptions[38].e==0.081000
assert record.descriptions[39].title=="gi|2499207|sp|Q58289|NIFH_METJA NITROGENASE IRON PROTEIN (NITRO..."
assert record.descriptions[39].score==35
assert record.descriptions[39].e==0.110000
assert record.descriptions[40].title=="gi|2495840|sp|Q57633|Y169_METJA HYPOTHETICAL ATP-BINDING PROTEI..."
assert record.descriptions[40].score==35
assert record.descriptions[40].e==0.110000
assert record.descriptions[41].title=="gi|2496129|sp|Q58233|Y823_METJA HYPOTHETICAL ATP-BINDING PROTEI..."
assert record.descriptions[41].score==35
assert record.descriptions[41].e==0.140000
assert record.descriptions[42].title=="gi|1171022|sp|P45135|MRP_HAEIN MRP PROTEIN HOMOLOG"
assert record.descriptions[42].score==34
assert record.descriptions[42].e==0.180000
assert record.descriptions[43].title=="gi|132518|sp|P22670|RFX1_HUMAN MHC CLASS II REGULATORY FACTOR R..."
assert record.descriptions[43].score==34
assert record.descriptions[43].e==0.180000
assert record.descriptions[44].title=="gi|113780|sp|P08144|AMYA_DROME ALPHA-AMYLASE A PRECURSOR (1,4-A..."
assert record.descriptions[44].score==34
assert record.descriptions[44].e==0.180000
assert record.descriptions[45].title=="gi|1709471|sp|P50980|OPPD_LACLC OLIGOPEPTIDE TRANSPORT ATP-BIND..."
assert record.descriptions[45].score==34
assert record.descriptions[45].e==0.180000
assert record.descriptions[46].title=="gi|6225724|sp|Q9ZE27|MRP_RICPR MRP PROTEIN HOMOLOG"
assert record.descriptions[46].score==34
assert record.descriptions[46].e==0.240000
assert record.descriptions[47].title=="gi|2496034|sp|Q57967|Y547_METJA HYPOTHETICAL ATP-BINDING PROTEI..."
assert record.descriptions[47].score==34
assert record.descriptions[47].e==0.240000
assert record.descriptions[48].title=="gi|128215|sp|P22548|NIF4_CLOPA NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[48].score==34
assert record.descriptions[48].e==0.240000
assert record.descriptions[49].title=="gi|1706104|sp|Q04663|CPSC_STRAG CPSC PROTEIN"
assert record.descriptions[49].score==34
assert record.descriptions[49].e==0.310000
assert record.descriptions[50].title=="gi|548362|sp|P26248|NIF1_AZOCH NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[50].score==34
assert record.descriptions[50].e==0.310000
assert record.descriptions[51].title=="gi|128203|sp|P00459|NIF1_AZOVI NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[51].score==34
assert record.descriptions[51].e==0.310000
assert record.descriptions[52].title=="gi|2497979|sp|Q57731|Y283_METJA HYPOTHETICAL PROTEIN MJ0283"
assert record.descriptions[52].score==33
assert record.descriptions[52].e==0.410000
assert record.descriptions[53].title=="gi|2499205|sp|Q44044|NIFH_ALCFA NITROGENASE IRON PROTEIN (NITRO..."
assert record.descriptions[53].score==33
assert record.descriptions[53].e==0.410000
assert record.descriptions[54].title=="gi|1703363|sp|P52145|ARA2_ECOLI ARSENICAL PUMP-DRIVING ATPASE"
assert record.descriptions[54].score==33
assert record.descriptions[54].e==0.410000
assert record.descriptions[55].title=="gi|132821|sp|P27683|RK24_SPIOL 50S RIBOSOMAL PROTEIN L24, CHLOR..."
assert record.descriptions[55].score==33
assert record.descriptions[55].e==0.410000
assert record.descriptions[56].title=="gi|399100|sp|Q02431|BCHX_RHOSH CHLOROPHYLLIDE REDUCTASE 35.5 KD..."
assert record.descriptions[56].score==33
assert record.descriptions[56].e==0.410000
assert record.descriptions[57].title=="gi|128207|sp|P06118|NIF2_AZOCH NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[57].score==33
assert record.descriptions[57].e==0.410000
assert record.descriptions[58].title=="gi|128267|sp|P00458|NIFH_KLEPN NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[58].score==33
assert record.descriptions[58].e==0.410000
assert record.descriptions[59].title=="gi|128208|sp|P15335|NIF2_AZOVI NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[59].score==33
assert record.descriptions[59].e==0.410000
assert record.descriptions[60].title=="gi|1703299|sp|P54215|AMYA_DROMA ALPHA-AMYLASE A PRECURSOR (1,4-..."
assert record.descriptions[60].score==33
assert record.descriptions[60].e==0.530000
assert record.descriptions[61].title=="gi|114867|sp|P26177|BCHX_RHOCA CHLOROPHYLLIDE REDUCTASE 35.5 KD..."
assert record.descriptions[61].score==33
assert record.descriptions[61].e==0.530000
assert record.descriptions[62].title=="gi|128217|sp|P09555|NIF6_CLOPA NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[62].score==33
assert record.descriptions[62].e==0.530000
assert record.descriptions[63].title=="gi|128277|sp|P06661|NIFH_THIFE NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[63].score==33
assert record.descriptions[63].e==0.530000
assert record.descriptions[64].title=="gi|1709270|sp|P54800|NIF2_METBA NITROGENASE IRON PROTEIN 2 (NIT..."
assert record.descriptions[64].score==32
assert record.descriptions[64].e==0.700000
assert record.descriptions[65].title=="gi|6014740|sp|O33946|CTB1_ACILW MUCONATE CYCLOISOMERASE I 1 (CI..."
assert record.descriptions[65].score==32
assert record.descriptions[65].e==0.920000
assert record.descriptions[66].title=="gi|134731|sp|P08866|SOPA_ECOLI SOPA PROTEIN (PROTEIN A)"
assert record.descriptions[66].score==32
assert record.descriptions[66].e==0.920000
assert record.descriptions[67].title=="gi|114220|sp|P08690|ARA1_ECOLI ARSENICAL PUMP-DRIVING ATPASE"
assert record.descriptions[67].score==32
assert record.descriptions[67].e==0.920000
assert record.descriptions[68].title=="gi|3334342|sp|O29633|SR54_ARCFU PROBABLE SIGNAL RECOGNITION 54 ..."
assert record.descriptions[68].score==32
assert record.descriptions[68].e==1.200000
assert record.descriptions[69].title=="gi|2499206|sp|Q59270|NIFH_CLOCB NITROGENASE IRON PROTEIN (NITRO..."
assert record.descriptions[69].score==32
assert record.descriptions[69].e==1.200000
assert record.descriptions[70].title=="gi|2842582|sp|Q58334|Y924_METJA HYPOTHETICAL ATP-BINDING PROTEI..."
assert record.descriptions[70].score==31
assert record.descriptions[70].e==1.600000
assert record.descriptions[71].title=="gi|128264|sp|P00463|NIFH_BRASP NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[71].score==31
assert record.descriptions[71].e==1.600000
assert record.descriptions[72].title=="gi|128268|sp|P06119|NIFH_METVO NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[72].score==31
assert record.descriptions[72].e==1.600000
assert record.descriptions[73].title=="gi|128263|sp|P06117|NIFH_BRAJA NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[73].score==31
assert record.descriptions[73].e==1.600000
assert record.descriptions[74].title=="gi|2833531|sp|Q58098|Y685_METJA HYPOTHETICAL ATP-BINDING PROTEI..."
assert record.descriptions[74].score==31
assert record.descriptions[74].e==2.100000
assert record.descriptions[75].title=="gi|1709997|sp|P53692|RA18_SCHPO DNA REPAIR PROTEIN RAD18"
assert record.descriptions[75].score==31
assert record.descriptions[75].e==2.100000
assert record.descriptions[76].title=="gi|127221|sp|P22900|MOBD_THIFE MOBD PROTEIN"
assert record.descriptions[76].score==31
assert record.descriptions[76].e==2.100000
assert record.descriptions[77].title=="gi|114863|sp|P26237|BCHL_RHOCA PROTOCHLOROPHILLIDE REDUCTASE 33..."
assert record.descriptions[77].score==31
assert record.descriptions[77].e==2.700000
assert record.descriptions[78].title=="gi|128216|sp|P09554|NIF5_CLOPA NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[78].score==31
assert record.descriptions[78].e==2.700000
assert record.descriptions[79].title=="gi|128204|sp|P00456|NIF1_CLOPA NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[79].score==31
assert record.descriptions[79].e==2.700000
assert record.descriptions[80].title=="gi|128209|sp|P09552|NIF2_CLOPA NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[80].score==31
assert record.descriptions[80].e==2.700000
assert record.descriptions[81].title=="gi|128261|sp|P00457|NIFH_ANASP NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[81].score==31
assert record.descriptions[81].e==2.700000
assert record.descriptions[82].title=="gi|266624|sp|Q00240|NIFH_PLEBO NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[82].score==31
assert record.descriptions[82].e==2.700000
assert record.descriptions[83].title=="gi|6225393|sp|O67066|FTSY_AQUAE CELL DIVISION PROTEIN FTSY HOMOLOG"
assert record.descriptions[83].score==30
assert record.descriptions[83].e==3.500000
assert record.descriptions[84].title=="gi|6166125|sp|O14936|CSKP_HUMAN PERIPHERAL PLASMA MEMBRANE PROT..."
assert record.descriptions[84].score==30
assert record.descriptions[84].e==3.500000
assert record.descriptions[85].title=="gi|1709099|sp|P53381|MRP_CLOPE MRP PROTEIN HOMOLOG"
assert record.descriptions[85].score==30
assert record.descriptions[85].e==3.500000
assert record.descriptions[86].title=="gi|1171710|sp|P46034|NIFH_FRASP NITROGENASE IRON PROTEIN (NITRO..."
assert record.descriptions[86].score==30
assert record.descriptions[86].e==3.500000
assert record.descriptions[87].title=="gi|128212|sp|P26252|NIF2_RHISO NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[87].score==30
assert record.descriptions[87].e==3.500000
assert record.descriptions[88].title=="gi|128206|sp|P26251|NIF1_RHISO NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[88].score==30
assert record.descriptions[88].e==3.500000
assert record.descriptions[89].title=="gi|2495787|sp|Q60392|Y084_METJA HYPOTHETICAL ATP-BINDING PROTEI..."
assert record.descriptions[89].score==30
assert record.descriptions[89].e==4.600000
assert record.descriptions[90].title=="gi|1709283|sp|P52336|NIFH_NOSSN NITROGENASE IRON PROTEIN (NITRO..."
assert record.descriptions[90].score==30
assert record.descriptions[90].e==4.600000
assert record.descriptions[91].title=="gi|731773|sp|P40558|YIA3_YEAST HYPOTHETICAL 31.9 KD PROTEIN IN ..."
assert record.descriptions[91].score==30
assert record.descriptions[91].e==4.600000
assert record.descriptions[92].title=="gi|128266|sp|P08925|NIFH_FRAAL NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[92].score==30
assert record.descriptions[92].e==4.600000
assert record.descriptions[93].title=="gi|732146|sp|P40742|YLXH_BACSU HYPOTHETICAL 33.2 KD PROTEIN IN ..."
assert record.descriptions[93].score==30
assert record.descriptions[93].e==4.600000
assert record.descriptions[94].title=="gi|417362|sp|P33178|NIFH_ANASL NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[94].score==30
assert record.descriptions[94].e==4.600000
assert record.descriptions[95].title=="gi|1709281|sp|P26250|NIFH_NOSCO NITROGENASE IRON PROTEIN (NITRO..."
assert record.descriptions[95].score==30
assert record.descriptions[95].e==4.600000
assert record.descriptions[96].title=="gi|6225290|sp|Q9Z900|EFP_CHLPN ELONGATION FACTOR P (EF-P)"
assert record.descriptions[96].score==29
assert record.descriptions[96].e==6.100000
assert record.descriptions[97].title=="gi|3183458|sp|P75796|YLIA_ECOLI HYPOTHETICAL ABC TRANSPORTER AT..."
assert record.descriptions[97].score==29
assert record.descriptions[97].e==6.100000
assert record.descriptions[98].title=="gi|2499682|sp|Q52312|INC1_ECOLI INCC PROTEIN"
assert record.descriptions[98].score==29
assert record.descriptions[98].e==6.100000
assert record.descriptions[99].title=="gi|1345941|sp|P15555|DAC_STRSQ D-ALANYL-D-ALANINE CARBOXYPEPTID..."
assert record.descriptions[99].score==29
assert record.descriptions[99].e==6.100000
assert record.descriptions[100].title=="gi|128275|sp|P08718|NIF1_RHOCA NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[100].score==29
assert record.descriptions[100].e==6.100000
assert record.descriptions[101].title=="gi|267440|sp|P29811|YAC3_PSEAE HYPOTHETICAL 23.9 KD PROTEIN IN ..."
assert record.descriptions[101].score==29
assert record.descriptions[101].e==6.100000
assert record.descriptions[102].title=="gi|732050|sp|P39342|YJGR_ECOLI HYPOTHETICAL 54.3 KD PROTEIN IN ..."
assert record.descriptions[102].score==29
assert record.descriptions[102].e==6.100000
assert record.descriptions[103].title=="gi|2833520|sp|Q58000|Y580_METJA HYPOTHETICAL PROTEIN MJ0580"
assert record.descriptions[103].score==29
assert record.descriptions[103].e==7.900000
assert record.descriptions[104].title=="gi|128272|sp|P00460|NIFH_RHIME NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[104].score==29
assert record.descriptions[104].e==7.900000
assert record.descriptions[105].title=="gi|128271|sp|P00461|NIFH_RHILT NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[105].score==29
assert record.descriptions[105].e==7.900000
assert record.descriptions[106].title=="gi|128276|sp|P22921|NIFH_RHORU NITROGENASE IRON PROTEIN (NITROG..."
assert record.descriptions[106].score==29
assert record.descriptions[106].e==7.900000
assert len(record.alignments)==0
assert record.database_name==['/taiyang/jchang/parsers/swissprot']
assert record.posted_date==[('May 12, 2000 11:25 AM',)]
assert record.num_letters_in_database==[29652561]
assert record.num_sequences_in_database==[82258]
assert record.ka_params==[0.315, 0.130, 0.349]
assert record.ka_params_gap==[0.270, 0.0470, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==11011949
assert record.num_sequences==82258
assert record.num_extends==406641
assert record.num_good_extends==1349
assert record.num_seqs_better_e==107
assert record.hsps_no_gap==89
assert record.hsps_prelim_gapped==18
assert record.hsps_gapped==131
assert record.query_length==222
assert record.database_length==29652561
assert record.effective_hsp_length==57
assert record.effective_query_length==165
assert record.effective_database_length==24963855
assert record.effective_search_space==4119036075
assert record.effective_search_space_used==4119036075
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.6)
assert record.blast_cutoff==(63,29.0)

handle = open('Blast/bt060')
record = pb_parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.12'
assert record.date=="Apr-21-2000"
assert record.reference==reference
assert record.query=="aaaa"
assert record.query_letters==889
assert record.database=="/dbase/swissprot/main/release/sp;/dbase/swissprot/main/update/spu"
assert record.database_sequences==88201
assert record.database_letters==31957340
assert len(record.rounds)==5
assert len(record.rounds[0].new_seqs)==27
assert record.rounds[0].new_seqs[0].title=="100K_RAT Q62671 rattus norvegicus (rat). 100 kda protein (ec..."
assert record.rounds[0].new_seqs[0].score==1516
assert record.rounds[0].new_seqs[0].e==0.000000
assert record.rounds[0].new_seqs[1].title=="HYDP_DROME P51592 drosophila melanogaster (fruit fly). hyper..."
assert record.rounds[0].new_seqs[1].score==513
assert record.rounds[0].new_seqs[1].e==9.9999999999999991e-146
assert record.rounds[0].new_seqs[2].title=="PUB1_SCHPO Q92462 schizosaccharomyces pombe (fission yeast)...."
assert record.rounds[0].new_seqs[2].score==95
assert record.rounds[0].new_seqs[2].e==7.0000000000000003e-19
assert record.rounds[0].new_seqs[3].title=="RSP5_YEAST P39940 saccharomyces cerevisiae (baker's yeast). ..."
assert record.rounds[0].new_seqs[3].score==93
assert record.rounds[0].new_seqs[3].e==2.9999999999999998e-18
assert record.rounds[0].new_seqs[4].title=="NED4_HUMAN P46934 homo sapiens (human). nedd-4 protein (ec 6..."
assert record.rounds[0].new_seqs[4].score==88
assert record.rounds[0].new_seqs[4].e==8.9999999999999996e-17
assert record.rounds[0].new_seqs[5].title=="NED4_MOUSE P46935 mus musculus (mouse). nedd-4 protein (ec 6..."
assert record.rounds[0].new_seqs[5].score==87
assert record.rounds[0].new_seqs[5].e==9.9999999999999998e-17
assert record.rounds[0].new_seqs[6].title=="URB1_RAT P51593 rattus norvegicus (rat). dna binding protein..."
assert record.rounds[0].new_seqs[6].score==85
assert record.rounds[0].new_seqs[6].e==3.9999999999999999e-16
assert record.rounds[0].new_seqs[7].title=="UE3A_HUMAN Q05086 homo sapiens (human). ubiquitin-protein li..."
assert record.rounds[0].new_seqs[7].score==84
assert record.rounds[0].new_seqs[7].e==1.0000000000000001e-15
assert record.rounds[0].new_seqs[8].title=="HUL5_YEAST P53119 saccharomyces cerevisiae (baker's yeast). ..."
assert record.rounds[0].new_seqs[8].score==83
assert record.rounds[0].new_seqs[8].e==2.0000000000000002e-15
assert record.rounds[0].new_seqs[9].title=="UE3A_MOUSE O08759 mus musculus (mouse). ubiquitin-protein li..."
assert record.rounds[0].new_seqs[9].score==82
assert record.rounds[0].new_seqs[9].e==5.9999999999999997e-15
assert record.rounds[0].new_seqs[10].title=="HUL4_YEAST P40985 saccharomyces cerevisiae (baker's yeast). ..."
assert record.rounds[0].new_seqs[10].score==69
assert record.rounds[0].new_seqs[10].e==3e-11
assert record.rounds[0].new_seqs[11].title=="UFD4_YEAST P33202 saccharomyces cerevisiae (baker's yeast). ..."
assert record.rounds[0].new_seqs[11].score==58
assert record.rounds[0].new_seqs[11].e==7.0000000000000005e-08
assert record.rounds[0].new_seqs[12].title=="Y032_HUMAN Q15034 homo sapiens (human). hypothetical protein..."
assert record.rounds[0].new_seqs[12].score==56
assert record.rounds[0].new_seqs[12].e==2.9999999999999999e-07
assert record.rounds[0].new_seqs[13].title=="TR12_HUMAN Q14669 homo sapiens (human). thyroid receptor int..."
assert record.rounds[0].new_seqs[13].score==50
assert record.rounds[0].new_seqs[13].e==0.000020
assert record.rounds[0].new_seqs[14].title=="PAB1_MOUSE P29341 mus musculus (mouse). polyadenylate-bindin..."
assert record.rounds[0].new_seqs[14].score==49
assert record.rounds[0].new_seqs[14].e==0.000030
assert record.rounds[0].new_seqs[15].title=="PAB1_HUMAN P11940 homo sapiens (human). polyadenylate-bindin..."
assert record.rounds[0].new_seqs[15].score==49
assert record.rounds[0].new_seqs[15].e==0.000050
assert record.rounds[0].new_seqs[16].title=="PABP_XENLA P20965 xenopus laevis (african clawed frog). poly..."
assert record.rounds[0].new_seqs[16].score==48
assert record.rounds[0].new_seqs[16].e==0.000080
assert record.rounds[0].new_seqs[17].title=="PABP_DROME P21187 drosophila melanogaster (fruit fly). polya..."
assert record.rounds[0].new_seqs[17].score==41
assert record.rounds[0].new_seqs[17].e==0.008000
assert record.rounds[0].new_seqs[18].title=="PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). po..."
assert record.rounds[0].new_seqs[18].score==37
assert record.rounds[0].new_seqs[18].e==0.150000
assert record.rounds[0].new_seqs[19].title=="PAB2_ARATH P42731 arabidopsis thaliana (mouse-ear cress). po..."
assert record.rounds[0].new_seqs[19].score==34
assert record.rounds[0].new_seqs[19].e==0.990000
assert record.rounds[0].new_seqs[20].title=="MCM3_MOUSE P25206 mus musculus (mouse). dna replication lice..."
assert record.rounds[0].new_seqs[20].score==34
assert record.rounds[0].new_seqs[20].e==1.700000
assert record.rounds[0].new_seqs[21].title=="PABP_SCHPO P31209 schizosaccharomyces pombe (fission yeast)...."
assert record.rounds[0].new_seqs[21].score==32
assert record.rounds[0].new_seqs[21].e==3.800000
assert record.rounds[0].new_seqs[22].title=="PROB_SERMA P17856 serratia marcescens. glutamate 5-kinase (e..."
assert record.rounds[0].new_seqs[22].score==32
assert record.rounds[0].new_seqs[22].e==6.600000
assert record.rounds[0].new_seqs[23].title=="ST20_CANAL Q92212 candida albicans (yeast). serine/threonine..."
assert record.rounds[0].new_seqs[23].score==31
assert record.rounds[0].new_seqs[23].e==8.600000
assert record.rounds[0].new_seqs[24].title=="KEND_HUMAN O95613 homo sapiens (human). kendrin (kiaa0402). ..."
assert record.rounds[0].new_seqs[24].score==31
assert record.rounds[0].new_seqs[24].e==8.600000
assert record.rounds[0].new_seqs[25].title=="FIXK_RHIME P13295 rhizobium meliloti (sinorhizobium meliloti..."
assert record.rounds[0].new_seqs[25].score==31
assert record.rounds[0].new_seqs[25].e==8.600000
assert record.rounds[0].new_seqs[26].title=="CC24_YEAST P11433 saccharomyces cerevisiae (baker's yeast). ..."
assert record.rounds[0].new_seqs[26].score==31
assert record.rounds[0].new_seqs[26].e==8.600000
assert len(record.rounds[0].alignments)==27
assert record.rounds[0].alignments[0].title==">100K_RAT Q62671 rattus norvegicus (rat). 100 kda protein (ec 6.3.2.-). 7/1999"
assert record.rounds[0].alignments[0].length==889
assert record.rounds[0].alignments[1].title==">HYDP_DROME P51592 drosophila melanogaster (fruit fly). hyperplastic discs protein (hyd protein) (ec 6.3.2.-). 12/1998"
assert record.rounds[0].alignments[1].length==2895
assert record.rounds[0].alignments[2].title==">PUB1_SCHPO Q92462 schizosaccharomyces pombe (fission yeast). ubiquitin--protein ligase pub1 (ec 6.3.2.-). 12/1998"
assert record.rounds[0].alignments[2].length==767
assert record.rounds[0].alignments[3].title==">RSP5_YEAST P39940 saccharomyces cerevisiae (baker's yeast). ubiquitin--protein ligase rsp5 (ec 6.3.2.-). 7/1999"
assert record.rounds[0].alignments[3].length==809
assert record.rounds[0].alignments[4].title==">NED4_HUMAN P46934 homo sapiens (human). nedd-4 protein (ec 6.3.2.-) (kiaa0093) (fragment). 7/1999"
assert record.rounds[0].alignments[4].length==927
assert record.rounds[0].alignments[5].title==">NED4_MOUSE P46935 mus musculus (mouse). nedd-4 protein (ec 6.3.2.-) (fragment). 11/1997"
assert record.rounds[0].alignments[5].length==957
assert record.rounds[0].alignments[6].title==">URB1_RAT P51593 rattus norvegicus (rat). dna binding protein ure-b1 (ec 6.3.2.-). 10/1996"
assert record.rounds[0].alignments[6].length==310
assert record.rounds[0].alignments[7].title==">UE3A_HUMAN Q05086 homo sapiens (human). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap) (human papillomavirus e6-associated protein). 5/2000"
assert record.rounds[0].alignments[7].length==875
assert record.rounds[0].alignments[8].title==">HUL5_YEAST P53119 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul5 (ec 6.3.2.-). 5/2000"
assert record.rounds[0].alignments[8].length==910
assert record.rounds[0].alignments[9].title==">UE3A_MOUSE O08759 mus musculus (mouse). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap). 5/2000"
assert record.rounds[0].alignments[9].length==885
assert record.rounds[0].alignments[10].title==">HUL4_YEAST P40985 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul4 (ec 6.3.2.-). 5/2000"
assert record.rounds[0].alignments[10].length==892
assert record.rounds[0].alignments[11].title==">UFD4_YEAST P33202 saccharomyces cerevisiae (baker's yeast). ubiquitin fusion degradation protein 4 (ub fusion protein 4). 11/1997"
assert record.rounds[0].alignments[11].length==1483
assert record.rounds[0].alignments[12].title==">Y032_HUMAN Q15034 homo sapiens (human). hypothetical protein kiaa0032. 5/2000"
assert record.rounds[0].alignments[12].length==1050
assert record.rounds[0].alignments[13].title==">TR12_HUMAN Q14669 homo sapiens (human). thyroid receptor interacting protein 12 (trip12) (kiaa0045). 12/1998"
assert record.rounds[0].alignments[13].length==1992
assert record.rounds[0].alignments[14].title==">PAB1_MOUSE P29341 mus musculus (mouse). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[0].alignments[14].length==636
assert record.rounds[0].alignments[15].title==">PAB1_HUMAN P11940 homo sapiens (human). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[0].alignments[15].length==636
assert record.rounds[0].alignments[16].title==">PABP_XENLA P20965 xenopus laevis (african clawed frog). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[0].alignments[16].length==633
assert record.rounds[0].alignments[17].title==">PABP_DROME P21187 drosophila melanogaster (fruit fly). polyadenylate-binding protein (poly(a) binding protein) (pabp). 2/1995"
assert record.rounds[0].alignments[17].length==632
assert record.rounds[0].alignments[18].title==">PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 5 (poly(a) binding protein 5) (pabp 5). 11/1995"
assert record.rounds[0].alignments[18].length==668
assert record.rounds[0].alignments[19].title==">PAB2_ARATH P42731 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 2 (poly(a) binding protein 2) (pabp 2). 12/1998"
assert record.rounds[0].alignments[19].length==629
assert record.rounds[0].alignments[20].title==">MCM3_MOUSE P25206 mus musculus (mouse). dna replication licensing factor mcm3 (dna polymerase alpha holoenzyme-associated protein p1) (p1-mcm3). 5/2000"
assert record.rounds[0].alignments[20].length==812
assert record.rounds[0].alignments[21].title==">PABP_SCHPO P31209 schizosaccharomyces pombe (fission yeast). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[0].alignments[21].length==653
assert record.rounds[0].alignments[22].title==">PROB_SERMA P17856 serratia marcescens. glutamate 5-kinase (ec 2.7.2.11) (gamma-glutamyl kinase) (gk). 10/1996"
assert record.rounds[0].alignments[22].length==367
assert record.rounds[0].alignments[23].title==">ST20_CANAL Q92212 candida albicans (yeast). serine/threonine-protein kinase ste20 homolog (ec 2.7.1.-). 5/2000"
assert record.rounds[0].alignments[23].length==1230
assert record.rounds[0].alignments[24].title==">KEND_HUMAN O95613 homo sapiens (human). kendrin (kiaa0402). 5/2000"
assert record.rounds[0].alignments[24].length==3321
assert record.rounds[0].alignments[25].title==">FIXK_RHIME P13295 rhizobium meliloti (sinorhizobium meliloti). nitrogen fixation regulation protein fixk. 5/2000"
assert record.rounds[0].alignments[25].length==211
assert record.rounds[0].alignments[26].title==">CC24_YEAST P11433 saccharomyces cerevisiae (baker's yeast). cell division control protein 24 (calcium regulatory protein). 7/1999"
assert record.rounds[0].alignments[26].length==854
assert len(record.rounds[1].new_seqs)==9
assert record.rounds[1].new_seqs[0].title=="PABP_DROME P21187 drosophila melanogaster (fruit fly). polya..."
assert record.rounds[1].new_seqs[0].score==67
assert record.rounds[1].new_seqs[0].e==1e-10
assert record.rounds[1].new_seqs[1].title=="PABP_SCHPO P31209 schizosaccharomyces pombe (fission yeast)...."
assert record.rounds[1].new_seqs[1].score==44
assert record.rounds[1].new_seqs[1].e==0.001
assert record.rounds[1].new_seqs[2].title=="PAB2_ARATH P42731 arabidopsis thaliana (mouse-ear cress). po..."
assert record.rounds[1].new_seqs[2].score==43
assert record.rounds[1].new_seqs[2].e==0.003
assert record.rounds[1].new_seqs[3].title=="PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). po..."
assert record.rounds[1].new_seqs[3].score==43
assert record.rounds[1].new_seqs[3].e==0.003
assert record.rounds[1].new_seqs[4].title=="PABP_YEAST P04147 saccharomyces cerevisiae polyadenylate-bin..."
assert record.rounds[1].new_seqs[4].score==41
assert record.rounds[1].new_seqs[4].e==0.010
assert record.rounds[1].new_seqs[5].title=="RNR_AQUAE O67834 aquifex aeolicus. ribonuclease r (ec 3.1.-...."
assert record.rounds[1].new_seqs[5].score==33
assert record.rounds[1].new_seqs[5].e==2.2
assert record.rounds[1].new_seqs[6].title=="NIRA_EMENI P28348 emericella nidulans (aspergillus nidulans)..."
assert record.rounds[1].new_seqs[6].score==33
assert record.rounds[1].new_seqs[6].e==2.9
assert record.rounds[1].new_seqs[7].title=="YK44_YEAST P36023 saccharomyces cerevisiae (baker's yeast). ..."
assert record.rounds[1].new_seqs[7].score==31
assert record.rounds[1].new_seqs[7].e==8.4
assert record.rounds[1].new_seqs[8].title=="SYQ_YEAST P13188 saccharomyces cerevisiae (baker's yeast). g..."
assert record.rounds[1].new_seqs[8].score==31
assert record.rounds[1].new_seqs[8].e==8.4
assert len(record.rounds[1].alignments)==26
assert record.rounds[1].alignments[0].title==">100K_RAT Q62671 rattus norvegicus (rat). 100 kda protein (ec 6.3.2.-). 7/1999"
assert record.rounds[1].alignments[0].length==889
assert record.rounds[1].alignments[1].title==">HYDP_DROME P51592 drosophila melanogaster (fruit fly). hyperplastic discs protein (hyd protein) (ec 6.3.2.-). 12/1998"
assert record.rounds[1].alignments[1].length==2895
assert record.rounds[1].alignments[2].title==">NED4_HUMAN P46934 homo sapiens (human). nedd-4 protein (ec 6.3.2.-) (kiaa0093) (fragment). 7/1999"
assert record.rounds[1].alignments[2].length==927
assert record.rounds[1].alignments[3].title==">PUB1_SCHPO Q92462 schizosaccharomyces pombe (fission yeast). ubiquitin--protein ligase pub1 (ec 6.3.2.-). 12/1998"
assert record.rounds[1].alignments[3].length==767
assert record.rounds[1].alignments[4].title==">NED4_MOUSE P46935 mus musculus (mouse). nedd-4 protein (ec 6.3.2.-) (fragment). 11/1997"
assert record.rounds[1].alignments[4].length==957
assert record.rounds[1].alignments[5].title==">RSP5_YEAST P39940 saccharomyces cerevisiae (baker's yeast). ubiquitin--protein ligase rsp5 (ec 6.3.2.-). 7/1999"
assert record.rounds[1].alignments[5].length==809
assert record.rounds[1].alignments[6].title==">UE3A_HUMAN Q05086 homo sapiens (human). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap) (human papillomavirus e6-associated protein). 5/2000"
assert record.rounds[1].alignments[6].length==875
assert record.rounds[1].alignments[7].title==">UE3A_MOUSE O08759 mus musculus (mouse). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap). 5/2000"
assert record.rounds[1].alignments[7].length==885
assert record.rounds[1].alignments[8].title==">HUL4_YEAST P40985 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul4 (ec 6.3.2.-). 5/2000"
assert record.rounds[1].alignments[8].length==892
assert record.rounds[1].alignments[9].title==">URB1_RAT P51593 rattus norvegicus (rat). dna binding protein ure-b1 (ec 6.3.2.-). 10/1996"
assert record.rounds[1].alignments[9].length==310
assert record.rounds[1].alignments[10].title==">TR12_HUMAN Q14669 homo sapiens (human). thyroid receptor interacting protein 12 (trip12) (kiaa0045). 12/1998"
assert record.rounds[1].alignments[10].length==1992
assert record.rounds[1].alignments[11].title==">Y032_HUMAN Q15034 homo sapiens (human). hypothetical protein kiaa0032. 5/2000"
assert record.rounds[1].alignments[11].length==1050
assert record.rounds[1].alignments[12].title==">HUL5_YEAST P53119 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul5 (ec 6.3.2.-). 5/2000"
assert record.rounds[1].alignments[12].length==910
assert record.rounds[1].alignments[13].title==">UFD4_YEAST P33202 saccharomyces cerevisiae (baker's yeast). ubiquitin fusion degradation protein 4 (ub fusion protein 4). 11/1997"
assert record.rounds[1].alignments[13].length==1483
assert record.rounds[1].alignments[14].title==">PAB1_HUMAN P11940 homo sapiens (human). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[1].alignments[14].length==636
assert record.rounds[1].alignments[15].title==">PABP_XENLA P20965 xenopus laevis (african clawed frog). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[1].alignments[15].length==633
assert record.rounds[1].alignments[16].title==">PAB1_MOUSE P29341 mus musculus (mouse). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[1].alignments[16].length==636
assert record.rounds[1].alignments[17].title==">PABP_DROME P21187 drosophila melanogaster (fruit fly). polyadenylate-binding protein (poly(a) binding protein) (pabp). 2/1995"
assert record.rounds[1].alignments[17].length==632
assert record.rounds[1].alignments[18].title==">PABP_SCHPO P31209 schizosaccharomyces pombe (fission yeast). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[1].alignments[18].length==653
assert record.rounds[1].alignments[19].title==">PAB2_ARATH P42731 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 2 (poly(a) binding protein 2) (pabp 2). 12/1998"
assert record.rounds[1].alignments[19].length==629
assert record.rounds[1].alignments[20].title==">PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 5 (poly(a) binding protein 5) (pabp 5). 11/1995"
assert record.rounds[1].alignments[20].length==668
assert record.rounds[1].alignments[21].title==">PABP_YEAST P04147 saccharomyces cerevisiae polyadenylate-binding protein, cytoplasmic and nuclear (pabp) (ars consensus binding protein acbp-67) (polyadenylate tail-binding protein). 2/1996"
assert record.rounds[1].alignments[21].length==576
assert record.rounds[1].alignments[22].title==">RNR_AQUAE O67834 aquifex aeolicus. ribonuclease r (ec 3.1.-.-) (rnase r) (vacb protein homolog). 5/2000"
assert record.rounds[1].alignments[22].length==705
assert record.rounds[1].alignments[23].title==">NIRA_EMENI P28348 emericella nidulans (aspergillus nidulans). nitrogen assimilation transcription factor nira. 4/1993"
assert record.rounds[1].alignments[23].length==892
assert record.rounds[1].alignments[24].title==">YK44_YEAST P36023 saccharomyces cerevisiae (baker's yeast). putative 101.8 kda transcriptional regulatory protein in las1-ccp1 intergenic region. 2/1995"
assert record.rounds[1].alignments[24].length==863
assert record.rounds[1].alignments[25].title==">SYQ_YEAST P13188 saccharomyces cerevisiae (baker's yeast). glutaminyl-trna synthetase (ec 6.1.1.18) (glutamine--trna ligase) (glnrs). 11/1997"
assert record.rounds[1].alignments[25].length==809
assert len(record.rounds[2].new_seqs)==6
assert record.rounds[2].new_seqs[0].title=="PAB2_ARATH P42731 arabidopsis thaliana (mouse-ear cress). po..."
assert record.rounds[2].new_seqs[0].score==48
assert record.rounds[2].new_seqs[0].e==1e-04
assert record.rounds[2].new_seqs[1].title=="PABP_SCHPO P31209 schizosaccharomyces pombe (fission yeast)...."
assert record.rounds[2].new_seqs[1].score==45
assert record.rounds[2].new_seqs[1].e==7e-04
assert record.rounds[2].new_seqs[2].title=="PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). po..."
assert record.rounds[2].new_seqs[2].score==44
assert record.rounds[2].new_seqs[2].e==0.001
assert record.rounds[2].new_seqs[3].title=="PABP_YEAST P04147 saccharomyces cerevisiae polyadenylate-bin..."
assert record.rounds[2].new_seqs[3].score==42
assert record.rounds[2].new_seqs[3].e==0.006
assert record.rounds[2].new_seqs[4].title=="NIRA_EMENI P28348 emericella nidulans (aspergillus nidulans)..."
assert record.rounds[2].new_seqs[4].score==33
assert record.rounds[2].new_seqs[4].e==2.9
assert record.rounds[2].new_seqs[5].title=="RNR_AQUAE O67834 aquifex aeolicus. ribonuclease r (ec 3.1.-...."
assert record.rounds[2].new_seqs[5].score==32
assert record.rounds[2].new_seqs[5].e==3.7
assert len(record.rounds[2].alignments)==24
assert record.rounds[2].alignments[0].title==">100K_RAT Q62671 rattus norvegicus (rat). 100 kda protein (ec 6.3.2.-). 7/1999"
assert record.rounds[2].alignments[0].length==889
assert record.rounds[2].alignments[1].title==">HYDP_DROME P51592 drosophila melanogaster (fruit fly). hyperplastic discs protein (hyd protein) (ec 6.3.2.-). 12/1998"
assert record.rounds[2].alignments[1].length==2895
assert record.rounds[2].alignments[2].title==">PUB1_SCHPO Q92462 schizosaccharomyces pombe (fission yeast). ubiquitin--protein ligase pub1 (ec 6.3.2.-). 12/1998"
assert record.rounds[2].alignments[2].length==767
assert record.rounds[2].alignments[3].title==">NED4_HUMAN P46934 homo sapiens (human). nedd-4 protein (ec 6.3.2.-) (kiaa0093) (fragment). 7/1999"
assert record.rounds[2].alignments[3].length==927
assert record.rounds[2].alignments[4].title==">NED4_MOUSE P46935 mus musculus (mouse). nedd-4 protein (ec 6.3.2.-) (fragment). 11/1997"
assert record.rounds[2].alignments[4].length==957
assert record.rounds[2].alignments[5].title==">RSP5_YEAST P39940 saccharomyces cerevisiae (baker's yeast). ubiquitin--protein ligase rsp5 (ec 6.3.2.-). 7/1999"
assert record.rounds[2].alignments[5].length==809
assert record.rounds[2].alignments[6].title==">UE3A_HUMAN Q05086 homo sapiens (human). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap) (human papillomavirus e6-associated protein). 5/2000"
assert record.rounds[2].alignments[6].length==875
assert record.rounds[2].alignments[7].title==">UE3A_MOUSE O08759 mus musculus (mouse). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap). 5/2000"
assert record.rounds[2].alignments[7].length==885
assert record.rounds[2].alignments[8].title==">HUL4_YEAST P40985 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul4 (ec 6.3.2.-). 5/2000"
assert record.rounds[2].alignments[8].length==892
assert record.rounds[2].alignments[9].title==">TR12_HUMAN Q14669 homo sapiens (human). thyroid receptor interacting protein 12 (trip12) (kiaa0045). 12/1998"
assert record.rounds[2].alignments[9].length==1992
assert record.rounds[2].alignments[10].title==">URB1_RAT P51593 rattus norvegicus (rat). dna binding protein ure-b1 (ec 6.3.2.-). 10/1996"
assert record.rounds[2].alignments[10].length==310
assert record.rounds[2].alignments[11].title==">Y032_HUMAN Q15034 homo sapiens (human). hypothetical protein kiaa0032. 5/2000"
assert record.rounds[2].alignments[11].length==1050
assert record.rounds[2].alignments[12].title==">UFD4_YEAST P33202 saccharomyces cerevisiae (baker's yeast). ubiquitin fusion degradation protein 4 (ub fusion protein 4). 11/1997"
assert record.rounds[2].alignments[12].length==1483
assert record.rounds[2].alignments[13].title==">HUL5_YEAST P53119 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul5 (ec 6.3.2.-). 5/2000"
assert record.rounds[2].alignments[13].length==910
assert record.rounds[2].alignments[14].title==">PABP_DROME P21187 drosophila melanogaster (fruit fly). polyadenylate-binding protein (poly(a) binding protein) (pabp). 2/1995"
assert record.rounds[2].alignments[14].length==632
assert record.rounds[2].alignments[15].title==">PAB1_HUMAN P11940 homo sapiens (human). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[2].alignments[15].length==636
assert record.rounds[2].alignments[16].title==">PABP_XENLA P20965 xenopus laevis (african clawed frog). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[2].alignments[16].length==633
assert record.rounds[2].alignments[17].title==">PAB1_MOUSE P29341 mus musculus (mouse). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[2].alignments[17].length==636
assert record.rounds[2].alignments[18].title==">PAB2_ARATH P42731 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 2 (poly(a) binding protein 2) (pabp 2). 12/1998"
assert record.rounds[2].alignments[18].length==629
assert record.rounds[2].alignments[19].title==">PABP_SCHPO P31209 schizosaccharomyces pombe (fission yeast). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[2].alignments[19].length==653
assert record.rounds[2].alignments[20].title==">PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 5 (poly(a) binding protein 5) (pabp 5). 11/1995"
assert record.rounds[2].alignments[20].length==668
assert record.rounds[2].alignments[21].title==">PABP_YEAST P04147 saccharomyces cerevisiae polyadenylate-binding protein, cytoplasmic and nuclear (pabp) (ars consensus binding protein acbp-67) (polyadenylate tail-binding protein). 2/1996"
assert record.rounds[2].alignments[21].length==576
assert record.rounds[2].alignments[22].title==">NIRA_EMENI P28348 emericella nidulans (aspergillus nidulans). nitrogen assimilation transcription factor nira. 4/1993"
assert record.rounds[2].alignments[22].length==892
assert record.rounds[2].alignments[23].title==">RNR_AQUAE O67834 aquifex aeolicus. ribonuclease r (ec 3.1.-.-) (rnase r) (vacb protein homolog). 5/2000"
assert record.rounds[2].alignments[23].length==705
assert len(record.rounds[3].new_seqs)==4
assert record.rounds[3].new_seqs[0].title=="PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). po..."
assert record.rounds[3].new_seqs[0].score==51
assert record.rounds[3].new_seqs[0].e== 9e-06
assert record.rounds[3].new_seqs[1].title=="PABP_YEAST P04147 saccharomyces cerevisiae polyadenylate-bin..."
assert record.rounds[3].new_seqs[1].score==45
assert record.rounds[3].new_seqs[1].e== 5e-04
assert record.rounds[3].new_seqs[2].title=="RNR_AQUAE O67834 aquifex aeolicus. ribonuclease r (ec 3.1.-...."
assert record.rounds[3].new_seqs[2].score==33
assert record.rounds[3].new_seqs[2].e== 2.9
assert record.rounds[3].new_seqs[3].title=="NIRA_EMENI P28348 emericella nidulans (aspergillus nidulans)..."
assert record.rounds[3].new_seqs[3].score==33
assert record.rounds[3].new_seqs[3].e==2.9
assert len(record.rounds[3].alignments)==24
assert record.rounds[3].alignments[0].title==">100K_RAT Q62671 rattus norvegicus (rat). 100 kda protein (ec 6.3.2.-). 7/1999"
assert record.rounds[3].alignments[0].length==889
assert record.rounds[3].alignments[1].title==">HYDP_DROME P51592 drosophila melanogaster (fruit fly). hyperplastic discs protein (hyd protein) (ec 6.3.2.-). 12/1998"
assert record.rounds[3].alignments[1].length==2895
assert record.rounds[3].alignments[2].title==">PUB1_SCHPO Q92462 schizosaccharomyces pombe (fission yeast). ubiquitin--protein ligase pub1 (ec 6.3.2.-). 12/1998"
assert record.rounds[3].alignments[2].length==767
assert record.rounds[3].alignments[3].title==">NED4_HUMAN P46934 homo sapiens (human). nedd-4 protein (ec 6.3.2.-) (kiaa0093) (fragment). 7/1999"
assert record.rounds[3].alignments[3].length==927
assert record.rounds[3].alignments[4].title==">NED4_MOUSE P46935 mus musculus (mouse). nedd-4 protein (ec 6.3.2.-) (fragment). 11/1997"
assert record.rounds[3].alignments[4].length==957
assert record.rounds[3].alignments[5].title==">RSP5_YEAST P39940 saccharomyces cerevisiae (baker's yeast). ubiquitin--protein ligase rsp5 (ec 6.3.2.-). 7/1999"
assert record.rounds[3].alignments[5].length==809
assert record.rounds[3].alignments[6].title==">UE3A_HUMAN Q05086 homo sapiens (human). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap) (human papillomavirus e6-associated protein). 5/2000"
assert record.rounds[3].alignments[6].length==875
assert record.rounds[3].alignments[7].title==">UE3A_MOUSE O08759 mus musculus (mouse). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap). 5/2000"
assert record.rounds[3].alignments[7].length==885
assert record.rounds[3].alignments[8].title==">HUL4_YEAST P40985 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul4 (ec 6.3.2.-). 5/2000"
assert record.rounds[3].alignments[8].length==892
assert record.rounds[3].alignments[9].title==">TR12_HUMAN Q14669 homo sapiens (human). thyroid receptor interacting protein 12 (trip12) (kiaa0045). 12/1998"
assert record.rounds[3].alignments[9].length==1992
assert record.rounds[3].alignments[10].title==">URB1_RAT P51593 rattus norvegicus (rat). dna binding protein ure-b1 (ec 6.3.2.-). 10/1996"
assert record.rounds[3].alignments[10].length==310
assert record.rounds[3].alignments[11].title==">Y032_HUMAN Q15034 homo sapiens (human). hypothetical protein kiaa0032. 5/2000"
assert record.rounds[3].alignments[11].length==1050
assert record.rounds[3].alignments[12].title==">UFD4_YEAST P33202 saccharomyces cerevisiae (baker's yeast). ubiquitin fusion degradation protein 4 (ub fusion protein 4). 11/1997"
assert record.rounds[3].alignments[12].length==1483
assert record.rounds[3].alignments[13].title==">HUL5_YEAST P53119 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul5 (ec 6.3.2.-). 5/2000"
assert record.rounds[3].alignments[13].length==910
assert record.rounds[3].alignments[14].title==">PAB1_HUMAN P11940 homo sapiens (human). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[3].alignments[14].length==636
assert record.rounds[3].alignments[15].title==">PABP_DROME P21187 drosophila melanogaster (fruit fly). polyadenylate-binding protein (poly(a) binding protein) (pabp). 2/1995"
assert record.rounds[3].alignments[15].length==632
assert record.rounds[3].alignments[16].title==">PABP_XENLA P20965 xenopus laevis (african clawed frog). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[3].alignments[16].length==633
assert record.rounds[3].alignments[17].title==">PAB1_MOUSE P29341 mus musculus (mouse). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[3].alignments[17].length==636
assert record.rounds[3].alignments[18].title==">PABP_SCHPO P31209 schizosaccharomyces pombe (fission yeast). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[3].alignments[18].length==653
assert record.rounds[3].alignments[19].title==">PAB2_ARATH P42731 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 2 (poly(a) binding protein 2) (pabp 2). 12/1998"
assert record.rounds[3].alignments[19].length==629
assert record.rounds[3].alignments[20].title==">PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 5 (poly(a) binding protein 5) (pabp 5). 11/1995"
assert record.rounds[3].alignments[20].length==668
assert record.rounds[3].alignments[21].title==">PABP_YEAST P04147 saccharomyces cerevisiae polyadenylate-binding protein, cytoplasmic and nuclear (pabp) (ars consensus binding protein acbp-67) (polyadenylate tail-binding protein). 2/1996"
assert record.rounds[3].alignments[21].length==576
assert record.rounds[3].alignments[22].title==">RNR_AQUAE O67834 aquifex aeolicus. ribonuclease r (ec 3.1.-.-) (rnase r) (vacb protein homolog). 5/2000"
assert record.rounds[3].alignments[22].length==705
assert record.rounds[3].alignments[23].title==">NIRA_EMENI P28348 emericella nidulans (aspergillus nidulans). nitrogen assimilation transcription factor nira. 4/1993"
assert record.rounds[3].alignments[23].length==892
assert len(record.rounds[4].new_seqs)==2
assert record.rounds[4].new_seqs[0].title=="RNR_AQUAE O67834 aquifex aeolicus. ribonuclease r (ec 3.1.-...."
assert record.rounds[4].new_seqs[0].score==33
assert record.rounds[4].new_seqs[0].e==2.2
assert record.rounds[4].new_seqs[1].title=="NIRA_EMENI P28348 emericella nidulans (aspergillus nidulans)..."
assert record.rounds[4].new_seqs[1].score==33
assert record.rounds[4].new_seqs[1].e==2.9
assert len(record.rounds[4].alignments)==24
assert record.rounds[4].alignments[0].title==">100K_RAT Q62671 rattus norvegicus (rat). 100 kda protein (ec 6.3.2.-). 7/1999"
assert record.rounds[4].alignments[0].length==889
assert record.rounds[4].alignments[1].title==">HYDP_DROME P51592 drosophila melanogaster (fruit fly). hyperplastic discs protein (hyd protein) (ec 6.3.2.-). 12/1998"
assert record.rounds[4].alignments[1].length==2895
assert record.rounds[4].alignments[2].title==">PUB1_SCHPO Q92462 schizosaccharomyces pombe (fission yeast). ubiquitin--protein ligase pub1 (ec 6.3.2.-). 12/1998"
assert record.rounds[4].alignments[2].length==767
assert record.rounds[4].alignments[3].title==">NED4_HUMAN P46934 homo sapiens (human). nedd-4 protein (ec 6.3.2.-) (kiaa0093) (fragment). 7/1999"
assert record.rounds[4].alignments[3].length==927
assert record.rounds[4].alignments[4].title==">NED4_MOUSE P46935 mus musculus (mouse). nedd-4 protein (ec 6.3.2.-) (fragment). 11/1997"
assert record.rounds[4].alignments[4].length==957
assert record.rounds[4].alignments[5].title==">RSP5_YEAST P39940 saccharomyces cerevisiae (baker's yeast). ubiquitin--protein ligase rsp5 (ec 6.3.2.-). 7/1999"
assert record.rounds[4].alignments[5].length==809
assert record.rounds[4].alignments[6].title==">UE3A_HUMAN Q05086 homo sapiens (human). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap) (human papillomavirus e6-associated protein). 5/2000"
assert record.rounds[4].alignments[6].length==875
assert record.rounds[4].alignments[7].title==">UE3A_MOUSE O08759 mus musculus (mouse). ubiquitin-protein ligase e3a (ec 6.3.2.-) (oncogenic protein- associated protein e6-ap). 5/2000"
assert record.rounds[4].alignments[7].length==885
assert record.rounds[4].alignments[8].title==">HUL4_YEAST P40985 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul4 (ec 6.3.2.-). 5/2000"
assert record.rounds[4].alignments[8].length==892
assert record.rounds[4].alignments[9].title==">TR12_HUMAN Q14669 homo sapiens (human). thyroid receptor interacting protein 12 (trip12) (kiaa0045). 12/1998"
assert record.rounds[4].alignments[9].length==1992
assert record.rounds[4].alignments[10].title==">URB1_RAT P51593 rattus norvegicus (rat). dna binding protein ure-b1 (ec 6.3.2.-). 10/1996"
assert record.rounds[4].alignments[10].length==310
assert record.rounds[4].alignments[11].title==">Y032_HUMAN Q15034 homo sapiens (human). hypothetical protein kiaa0032. 5/2000"
assert record.rounds[4].alignments[11].length==1050
assert record.rounds[4].alignments[12].title==">UFD4_YEAST P33202 saccharomyces cerevisiae (baker's yeast). ubiquitin fusion degradation protein 4 (ub fusion protein 4). 11/1997"
assert record.rounds[4].alignments[12].length==1483
assert record.rounds[4].alignments[13].title==">HUL5_YEAST P53119 saccharomyces cerevisiae (baker's yeast). probable ubiquitin--protein ligase hul5 (ec 6.3.2.-). 5/2000"
assert record.rounds[4].alignments[13].length==910
assert record.rounds[4].alignments[14].title==">PABP_DROME P21187 drosophila melanogaster (fruit fly). polyadenylate-binding protein (poly(a) binding protein) (pabp). 2/1995"
assert record.rounds[4].alignments[14].length==632
assert record.rounds[4].alignments[15].title==">PAB1_HUMAN P11940 homo sapiens (human). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[4].alignments[15].length==636
assert record.rounds[4].alignments[16].title==">PAB1_MOUSE P29341 mus musculus (mouse). polyadenylate-binding protein 1 (poly(a) binding protein 1) (pabp 1). 7/1998"
assert record.rounds[4].alignments[16].length==636
assert record.rounds[4].alignments[17].title==">PABP_XENLA P20965 xenopus laevis (african clawed frog). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[4].alignments[17].length==633
assert record.rounds[4].alignments[18].title==">PABP_SCHPO P31209 schizosaccharomyces pombe (fission yeast). polyadenylate-binding protein (poly(a) binding protein) (pabp). 7/1998"
assert record.rounds[4].alignments[18].length==653
assert record.rounds[4].alignments[19].title==">PAB2_ARATH P42731 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 2 (poly(a) binding protein 2) (pabp 2). 12/1998"
assert record.rounds[4].alignments[19].length==629
assert record.rounds[4].alignments[20].title==">PAB5_ARATH Q05196 arabidopsis thaliana (mouse-ear cress). polyadenylate-binding protein 5 (poly(a) binding protein 5) (pabp 5). 11/1995"
assert record.rounds[4].alignments[20].length==668
assert record.rounds[4].alignments[21].title==">PABP_YEAST P04147 saccharomyces cerevisiae polyadenylate-binding protein, cytoplasmic and nuclear (pabp) (ars consensus binding protein acbp-67) (polyadenylate tail-binding protein). 2/1996"
assert record.rounds[4].alignments[21].length==576
assert record.rounds[4].alignments[22].title==">RNR_AQUAE O67834 aquifex aeolicus. ribonuclease r (ec 3.1.-.-) (rnase r) (vacb protein homolog). 5/2000"
assert record.rounds[4].alignments[22].length==705
assert record.rounds[4].alignments[23].title==">NIRA_EMENI P28348 emericella nidulans (aspergillus nidulans). nitrogen assimilation transcription factor nira. 4/1993"
assert record.rounds[4].alignments[23].length==892
assert record.rounds[0].alignments[0].hsps[0].score==3882
assert record.rounds[0].alignments[0].hsps[0].bits==1516
assert record.rounds[0].alignments[0].hsps[0].expect==0.0
assert len(record.rounds[0].alignments[0].hsps)==1
assert record.rounds[0].alignments[1].hsps[0].score==1308
assert record.rounds[0].alignments[1].hsps[0].bits==513
assert record.rounds[0].alignments[1].hsps[0].expect==1e-145
assert record.rounds[0].alignments[1].hsps[1].score==272
assert record.rounds[0].alignments[1].hsps[1].bits==110
assert record.rounds[0].alignments[1].hsps[1].expect==1e-23
assert len(record.rounds[0].alignments[1].hsps)==2
assert record.rounds[0].alignments[2].hsps[0].score==232
assert record.rounds[0].alignments[2].hsps[0].bits==94.8
assert record.rounds[0].alignments[2].hsps[0].expect==7e-19
assert len(record.rounds[0].alignments[2].hsps)==1
assert record.rounds[0].alignments[3].hsps[0].score==227
assert record.rounds[0].alignments[3].hsps[0].bits==92.8
assert record.rounds[0].alignments[3].hsps[0].expect==3e-18
assert len(record.rounds[0].alignments[3].hsps)==1
assert record.rounds[0].alignments[4].hsps[0].score==214
assert record.rounds[0].alignments[4].hsps[0].bits==87.8
assert record.rounds[0].alignments[4].hsps[0].expect==9e-17
assert len(record.rounds[0].alignments[4].hsps)==1
assert record.rounds[0].alignments[5].hsps[0].score==212
assert record.rounds[0].alignments[5].hsps[0].bits==87.0
assert record.rounds[0].alignments[5].hsps[0].expect==1e-16
assert len(record.rounds[0].alignments[5].hsps)==1
assert record.rounds[0].alignments[6].hsps[0].score==208
assert record.rounds[0].alignments[6].hsps[0].bits==85.4
assert record.rounds[0].alignments[6].hsps[0].expect==4e-16
assert len(record.rounds[0].alignments[6].hsps)==1
assert record.rounds[0].alignments[7].hsps[0].score==204
assert record.rounds[0].alignments[7].hsps[0].bits==83.9
assert record.rounds[0].alignments[7].hsps[0].expect==1e-15
assert len(record.rounds[0].alignments[7].hsps)==1
assert record.rounds[0].alignments[8].hsps[0].score==202
assert record.rounds[0].alignments[8].hsps[0].bits==83.1
assert record.rounds[0].alignments[8].hsps[0].expect==2e-15
assert len(record.rounds[0].alignments[8].hsps)==1
assert record.rounds[0].alignments[9].hsps[0].score==198
assert record.rounds[0].alignments[9].hsps[0].bits==81.5
assert record.rounds[0].alignments[9].hsps[0].expect==6e-15
assert len(record.rounds[0].alignments[9].hsps)==1
assert record.rounds[0].alignments[10].hsps[0].score==167
assert record.rounds[0].alignments[10].hsps[0].bits==69.5
assert record.rounds[0].alignments[10].hsps[0].expect==3e-11
assert len(record.rounds[0].alignments[10].hsps)==1
assert record.rounds[0].alignments[11].hsps[0].score==138
assert record.rounds[0].alignments[11].hsps[0].bits==58.2
assert record.rounds[0].alignments[11].hsps[0].expect==7e-08
assert len(record.rounds[0].alignments[11].hsps)==1
assert record.rounds[0].alignments[12].hsps[0].score==133
assert record.rounds[0].alignments[12].hsps[0].bits==56.2
assert record.rounds[0].alignments[12].hsps[0].expect==3e-07
assert len(record.rounds[0].alignments[12].hsps)==1
assert record.rounds[0].alignments[13].hsps[0].score==118
assert record.rounds[0].alignments[13].hsps[0].bits==50.4
assert record.rounds[0].alignments[13].hsps[0].expect==2e-05
assert len(record.rounds[0].alignments[13].hsps)==1
assert record.rounds[0].alignments[14].hsps[0].score==115
assert record.rounds[0].alignments[14].hsps[0].bits==49.2
assert record.rounds[0].alignments[14].hsps[0].expect==3e-05
assert len(record.rounds[0].alignments[14].hsps)==1
assert record.rounds[0].alignments[15].hsps[0].score==114
assert record.rounds[0].alignments[15].hsps[0].bits==48.8
assert record.rounds[0].alignments[15].hsps[0].expect==5e-05
assert len(record.rounds[0].alignments[15].hsps)==1
assert record.rounds[0].alignments[16].hsps[0].score==112
assert record.rounds[0].alignments[16].hsps[0].bits==48.0
assert record.rounds[0].alignments[16].hsps[0].expect==8e-05
assert len(record.rounds[0].alignments[16].hsps)==1
assert record.rounds[0].alignments[17].hsps[0].score==95
assert record.rounds[0].alignments[17].hsps[0].bits==41.4
assert record.rounds[0].alignments[17].hsps[0].expect==0.008
assert len(record.rounds[0].alignments[17].hsps)==1
assert record.rounds[0].alignments[18].hsps[0].score==84
assert record.rounds[0].alignments[18].hsps[0].bits==37.1
assert record.rounds[0].alignments[18].hsps[0].expect==0.15
assert len(record.rounds[0].alignments[18].hsps)==1
assert record.rounds[0].alignments[19].hsps[0].score==77
assert record.rounds[0].alignments[19].hsps[0].bits==34.4
assert record.rounds[0].alignments[19].hsps[0].expect==0.99
assert len(record.rounds[0].alignments[19].hsps)==1
assert record.rounds[0].alignments[20].hsps[0].score==75
assert record.rounds[0].alignments[20].hsps[0].bits==33.6
assert record.rounds[0].alignments[20].hsps[0].expect==1.7
assert len(record.rounds[0].alignments[20].hsps)==1
assert record.rounds[0].alignments[21].hsps[0].score==72
assert record.rounds[0].alignments[21].hsps[0].bits==32.5
assert record.rounds[0].alignments[21].hsps[0].expect==3.8
assert len(record.rounds[0].alignments[21].hsps)==1
assert record.rounds[0].alignments[22].hsps[0].score==70
assert record.rounds[0].alignments[22].hsps[0].bits==31.7
assert record.rounds[0].alignments[22].hsps[0].expect==6.6
assert len(record.rounds[0].alignments[22].hsps)==1
assert record.rounds[0].alignments[23].hsps[0].score==69
assert record.rounds[0].alignments[23].hsps[0].bits==31.3
assert record.rounds[0].alignments[23].hsps[0].expect==8.6
assert len(record.rounds[0].alignments[23].hsps)==1
assert record.rounds[0].alignments[24].hsps[0].score==69
assert record.rounds[0].alignments[24].hsps[0].bits==31.3
assert record.rounds[0].alignments[24].hsps[0].expect==8.6
assert len(record.rounds[0].alignments[24].hsps)==1
assert record.rounds[0].alignments[25].hsps[0].score==69
assert record.rounds[0].alignments[25].hsps[0].bits==31.3
assert record.rounds[0].alignments[25].hsps[0].expect==8.6
assert len(record.rounds[0].alignments[25].hsps)==1
assert record.rounds[0].alignments[26].hsps[0].score==69
assert record.rounds[0].alignments[26].hsps[0].bits==31.3
assert record.rounds[0].alignments[26].hsps[0].expect==8.6
assert record.rounds[1].alignments[0].hsps[0].score==3163
assert record.rounds[1].alignments[0].hsps[0].bits==1236
assert record.rounds[1].alignments[0].hsps[0].expect==0.0
assert len(record.rounds[1].alignments[0].hsps)==1
assert record.rounds[1].alignments[1].hsps[0].score==1810
assert record.rounds[1].alignments[1].hsps[0].bits==709
assert record.rounds[1].alignments[1].hsps[0].expect==0.0
assert record.rounds[1].alignments[1].hsps[1].score==649
assert record.rounds[1].alignments[1].hsps[1].bits==257
assert record.rounds[1].alignments[1].hsps[1].expect==8e-68
assert len(record.rounds[1].alignments[1].hsps)==2
assert record.rounds[1].alignments[2].hsps[0].score==873
assert record.rounds[1].alignments[2].hsps[0].bits==344
assert record.rounds[1].alignments[2].hsps[0].expect==4e-94
assert len(record.rounds[1].alignments[2].hsps)==1
assert record.rounds[1].alignments[3].hsps[0].score==870
assert record.rounds[1].alignments[3].hsps[0].bits==343
assert record.rounds[1].alignments[3].hsps[0].expect==1e-93
assert len(record.rounds[1].alignments[3].hsps)==1
assert record.rounds[1].alignments[4].hsps[0].score==870
assert record.rounds[1].alignments[4].hsps[0].bits==343
assert record.rounds[1].alignments[4].hsps[0].expect==1e-93
assert len(record.rounds[1].alignments[4].hsps)==1
assert record.rounds[1].alignments[5].hsps[0].score==825
assert record.rounds[1].alignments[5].hsps[0].bits==325
assert record.rounds[1].alignments[5].hsps[0].expect==2e-88
assert len(record.rounds[1].alignments[5].hsps)==1
assert record.rounds[1].alignments[6].hsps[0].score==795
assert record.rounds[1].alignments[6].hsps[0].bits==314
assert record.rounds[1].alignments[6].hsps[0].expect==6e-85
assert len(record.rounds[1].alignments[6].hsps)==1
assert record.rounds[1].alignments[7].hsps[0].score==765
assert record.rounds[1].alignments[7].hsps[0].bits==302
assert record.rounds[1].alignments[7].hsps[0].expect==2e-81
assert len(record.rounds[1].alignments[7].hsps)==1
assert record.rounds[1].alignments[8].hsps[0].score==711
assert record.rounds[1].alignments[8].hsps[0].bits==281
assert record.rounds[1].alignments[8].hsps[0].expect==4e-75
assert len(record.rounds[1].alignments[8].hsps)==1
assert record.rounds[1].alignments[9].hsps[0].score==709
assert record.rounds[1].alignments[9].hsps[0].bits==280
assert record.rounds[1].alignments[9].hsps[0].expect==8e-75
assert len(record.rounds[1].alignments[9].hsps)==1
assert record.rounds[1].alignments[10].hsps[0].score==687
assert record.rounds[1].alignments[10].hsps[0].bits==272
assert record.rounds[1].alignments[10].hsps[0].expect==3e-72
assert len(record.rounds[1].alignments[10].hsps)==1
assert record.rounds[1].alignments[11].hsps[0].score==684
assert record.rounds[1].alignments[11].hsps[0].bits==270
assert record.rounds[1].alignments[11].hsps[0].expect==6e-72
assert len(record.rounds[1].alignments[11].hsps)==1
assert record.rounds[1].alignments[12].hsps[0].score==633
assert record.rounds[1].alignments[12].hsps[0].bits==251
assert record.rounds[1].alignments[12].hsps[0].expect==6e-66
assert len(record.rounds[1].alignments[12].hsps)==1
assert record.rounds[1].alignments[13].hsps[0].score==566
assert record.rounds[1].alignments[13].hsps[0].bits==224
assert record.rounds[1].alignments[13].hsps[0].expect==4e-58
assert len(record.rounds[1].alignments[13].hsps)==1
assert record.rounds[1].alignments[14].hsps[0].score==216
assert record.rounds[1].alignments[14].hsps[0].bits==88.6
assert record.rounds[1].alignments[14].hsps[0].expect==5e-17
assert len(record.rounds[1].alignments[14].hsps)==1
assert record.rounds[1].alignments[15].hsps[0].score==215
assert record.rounds[1].alignments[15].hsps[0].bits==88.2
assert record.rounds[1].alignments[15].hsps[0].expect==6e-17
assert len(record.rounds[1].alignments[15].hsps)==1
assert record.rounds[1].alignments[16].hsps[0].score==214
assert record.rounds[1].alignments[16].hsps[0].bits==87.8
assert record.rounds[1].alignments[16].hsps[0].expect==8e-17
assert len(record.rounds[1].alignments[16].hsps)==1
assert record.rounds[1].alignments[17].hsps[0].score==161
assert record.rounds[1].alignments[17].hsps[0].bits==67.2
assert record.rounds[1].alignments[17].hsps[0].expect==1e-10
assert len(record.rounds[1].alignments[17].hsps)==1
assert record.rounds[1].alignments[18].hsps[0].score==102
assert record.rounds[1].alignments[18].hsps[0].bits==44.2
assert record.rounds[1].alignments[18].hsps[0].expect==0.001
assert len(record.rounds[1].alignments[18].hsps)==1
assert record.rounds[1].alignments[19].hsps[0].score==99
assert record.rounds[1].alignments[19].hsps[0].bits==43.0
assert record.rounds[1].alignments[19].hsps[0].expect==0.003
assert len(record.rounds[1].alignments[19].hsps)==1
assert record.rounds[1].alignments[20].hsps[0].score==98
assert record.rounds[1].alignments[20].hsps[0].bits==42.6
assert record.rounds[1].alignments[20].hsps[0].expect==0.003
assert len(record.rounds[1].alignments[20].hsps)==1
assert record.rounds[1].alignments[21].hsps[0].score==94
assert record.rounds[1].alignments[21].hsps[0].bits==41.1
assert record.rounds[1].alignments[21].hsps[0].expect==0.010
assert len(record.rounds[1].alignments[21].hsps)==1
assert record.rounds[1].alignments[22].hsps[0].score==74
assert record.rounds[1].alignments[22].hsps[0].bits==33.3
assert record.rounds[1].alignments[22].hsps[0].expect==2.2
assert len(record.rounds[1].alignments[22].hsps)==1
assert record.rounds[1].alignments[23].hsps[0].score==73
assert record.rounds[1].alignments[23].hsps[0].bits==32.9
assert record.rounds[1].alignments[23].hsps[0].expect==2.9
assert len(record.rounds[1].alignments[23].hsps)==1
assert record.rounds[1].alignments[24].hsps[0].score==69
assert record.rounds[1].alignments[24].hsps[0].bits==31.3
assert record.rounds[1].alignments[24].hsps[0].expect==8.4
assert len(record.rounds[1].alignments[24].hsps)==1
assert record.rounds[1].alignments[25].hsps[0].score==69
assert record.rounds[1].alignments[25].hsps[0].bits==31.3
assert record.rounds[1].alignments[25].hsps[0].expect==8.4
assert record.rounds[2].alignments[0].hsps[0].score==3143
assert record.rounds[2].alignments[0].hsps[0].bits==1228
assert record.rounds[2].alignments[0].hsps[0].expect==0.0
assert len(record.rounds[2].alignments[0].hsps)==1
assert record.rounds[2].alignments[1].hsps[0].score==1792
assert record.rounds[2].alignments[1].hsps[0].bits==702
assert record.rounds[2].alignments[1].hsps[0].expect==0.0
assert record.rounds[2].alignments[1].hsps[1].score==643
assert record.rounds[2].alignments[1].hsps[1].bits==254
assert record.rounds[2].alignments[1].hsps[1].expect==4e-67
assert len(record.rounds[2].alignments[1].hsps)==2
assert record.rounds[2].alignments[2].hsps[0].score==864
assert record.rounds[2].alignments[2].hsps[0].bits==340
assert record.rounds[2].alignments[2].hsps[0].expect==5e-93
assert len(record.rounds[2].alignments[2].hsps)==1
assert record.rounds[2].alignments[3].hsps[0].score==860
assert record.rounds[2].alignments[3].hsps[0].bits==339
assert record.rounds[2].alignments[3].hsps[0].expect==1e-92
assert len(record.rounds[2].alignments[3].hsps)==1
assert record.rounds[2].alignments[4].hsps[0].score==859
assert record.rounds[2].alignments[4].hsps[0].bits==339
assert record.rounds[2].alignments[4].hsps[0].expect==2e-92
assert len(record.rounds[2].alignments[4].hsps)==1
assert record.rounds[2].alignments[5].hsps[0].score==836
assert record.rounds[2].alignments[5].hsps[0].bits==330
assert record.rounds[2].alignments[5].hsps[0].expect==1e-89
assert len(record.rounds[2].alignments[5].hsps)==1
assert record.rounds[2].alignments[6].hsps[0].score==801
assert record.rounds[2].alignments[6].hsps[0].bits==316
assert record.rounds[2].alignments[6].hsps[0].expect==1e-85
assert len(record.rounds[2].alignments[6].hsps)==1
assert record.rounds[2].alignments[7].hsps[0].score==775
assert record.rounds[2].alignments[7].hsps[0].bits==306
assert record.rounds[2].alignments[7].hsps[0].expect==1e-82
assert len(record.rounds[2].alignments[7].hsps)==1
assert record.rounds[2].alignments[8].hsps[0].score==733
assert record.rounds[2].alignments[8].hsps[0].bits==289
assert record.rounds[2].alignments[8].hsps[0].expect==1e-77
assert len(record.rounds[2].alignments[8].hsps)==1
assert record.rounds[2].alignments[9].hsps[0].score==719
assert record.rounds[2].alignments[9].hsps[0].bits==284
assert record.rounds[2].alignments[9].hsps[0].expect==5e-76
assert len(record.rounds[2].alignments[9].hsps)==1
assert record.rounds[2].alignments[10].hsps[0].score==707
assert record.rounds[2].alignments[10].hsps[0].bits==279
assert record.rounds[2].alignments[10].hsps[0].expect==1e-74
assert len(record.rounds[2].alignments[10].hsps)==1
assert record.rounds[2].alignments[11].hsps[0].score==694
assert record.rounds[2].alignments[11].hsps[0].bits==274
assert record.rounds[2].alignments[11].hsps[0].expect==4e-73
assert len(record.rounds[2].alignments[11].hsps)==1
assert record.rounds[2].alignments[12].hsps[0].score==655
assert record.rounds[2].alignments[12].hsps[0].bits==259
assert record.rounds[2].alignments[12].hsps[0].expect==2e-68
assert len(record.rounds[2].alignments[12].hsps)==1
assert record.rounds[2].alignments[13].hsps[0].score==649
assert record.rounds[2].alignments[13].hsps[0].bits==257
assert record.rounds[2].alignments[13].hsps[0].expect==8e-68
assert len(record.rounds[2].alignments[13].hsps)==1
assert record.rounds[2].alignments[14].hsps[0].score==215
assert record.rounds[2].alignments[14].hsps[0].bits==88.2
assert record.rounds[2].alignments[14].hsps[0].expect==6e-17
assert len(record.rounds[2].alignments[14].hsps)==1
assert record.rounds[2].alignments[15].hsps[0].score==213
assert record.rounds[2].alignments[15].hsps[0].bits==87.4
assert record.rounds[2].alignments[15].hsps[0].expect==1e-16
assert len(record.rounds[2].alignments[15].hsps)==1
assert record.rounds[2].alignments[16].hsps[0].score==210
assert record.rounds[2].alignments[16].hsps[0].bits==86.2
assert record.rounds[2].alignments[16].hsps[0].expect==2e-16
assert len(record.rounds[2].alignments[16].hsps)==1
assert record.rounds[2].alignments[17].hsps[0].score==209
assert record.rounds[2].alignments[17].hsps[0].bits==85.9
assert record.rounds[2].alignments[17].hsps[0].expect==3e-16
assert len(record.rounds[2].alignments[17].hsps)==1
assert record.rounds[2].alignments[18].hsps[0].score==111
assert record.rounds[2].alignments[18].hsps[0].bits==47.7
assert record.rounds[2].alignments[18].hsps[0].expect==1e-04
assert len(record.rounds[2].alignments[18].hsps)==1
assert record.rounds[2].alignments[19].hsps[0].score==104
assert record.rounds[2].alignments[19].hsps[0].bits==45.0
assert record.rounds[2].alignments[19].hsps[0].expect==7e-04
assert len(record.rounds[2].alignments[19].hsps)==1
assert record.rounds[2].alignments[20].hsps[0].score==101
assert record.rounds[2].alignments[20].hsps[0].bits==43.8
assert record.rounds[2].alignments[20].hsps[0].expect==0.001
assert len(record.rounds[2].alignments[20].hsps)==1
assert record.rounds[2].alignments[21].hsps[0].score==96
assert record.rounds[2].alignments[21].hsps[0].bits==41.8
assert record.rounds[2].alignments[21].hsps[0].expect==0.006
assert len(record.rounds[2].alignments[21].hsps)==1
assert record.rounds[2].alignments[22].hsps[0].score==73
assert record.rounds[2].alignments[22].hsps[0].bits==32.9
assert record.rounds[2].alignments[22].hsps[0].expect==2.9
assert len(record.rounds[2].alignments[22].hsps)==1
assert record.rounds[2].alignments[23].hsps[0].score==72
assert record.rounds[2].alignments[23].hsps[0].bits==32.5
assert record.rounds[2].alignments[23].hsps[0].expect==3.7
assert record.rounds[3].alignments[0].hsps[0].score==3136
assert record.rounds[3].alignments[0].hsps[0].bits==1225
assert record.rounds[3].alignments[0].hsps[0].expect==0.0
assert len(record.rounds[3].alignments[0].hsps)==1
assert record.rounds[3].alignments[1].hsps[0].score==1776
assert record.rounds[3].alignments[1].hsps[0].bits==696
assert record.rounds[3].alignments[1].hsps[0].expect==0.0
assert record.rounds[3].alignments[1].hsps[1].score==643
assert record.rounds[3].alignments[1].hsps[1].bits==254
assert record.rounds[3].alignments[1].hsps[1].expect==4e-67
assert len(record.rounds[3].alignments[1].hsps)==2
assert record.rounds[3].alignments[2].hsps[0].score==866
assert record.rounds[3].alignments[2].hsps[0].bits==341
assert record.rounds[3].alignments[2].hsps[0].expect==3e-93
assert len(record.rounds[3].alignments[2].hsps)==1
assert record.rounds[3].alignments[3].hsps[0].score==862
assert record.rounds[3].alignments[3].hsps[0].bits==340
assert record.rounds[3].alignments[3].hsps[0].expect==9e-93
assert len(record.rounds[3].alignments[3].hsps)==1
assert record.rounds[3].alignments[4].hsps[0].score==861
assert record.rounds[3].alignments[4].hsps[0].bits==339
assert record.rounds[3].alignments[4].hsps[0].expect==1e-92
assert len(record.rounds[3].alignments[4].hsps)==1
assert record.rounds[3].alignments[5].hsps[0].score==839
assert record.rounds[3].alignments[5].hsps[0].bits==331
assert record.rounds[3].alignments[5].hsps[0].expect==4e-90
assert len(record.rounds[3].alignments[5].hsps)==1
assert record.rounds[3].alignments[6].hsps[0].score==802
assert record.rounds[3].alignments[6].hsps[0].bits==316
assert record.rounds[3].alignments[6].hsps[0].expect==9e-86
assert len(record.rounds[3].alignments[6].hsps)==1
assert record.rounds[3].alignments[7].hsps[0].score==777
assert record.rounds[3].alignments[7].hsps[0].bits==307
assert record.rounds[3].alignments[7].hsps[0].expect==8e-83
assert len(record.rounds[3].alignments[7].hsps)==1
assert record.rounds[3].alignments[8].hsps[0].score==735
assert record.rounds[3].alignments[8].hsps[0].bits==290
assert record.rounds[3].alignments[8].hsps[0].expect==7e-78
assert len(record.rounds[3].alignments[8].hsps)==1
assert record.rounds[3].alignments[9].hsps[0].score==720
assert record.rounds[3].alignments[9].hsps[0].bits==284
assert record.rounds[3].alignments[9].hsps[0].expect==4e-76
assert len(record.rounds[3].alignments[9].hsps)==1
assert record.rounds[3].alignments[10].hsps[0].score==708
assert record.rounds[3].alignments[10].hsps[0].bits==280
assert record.rounds[3].alignments[10].hsps[0].expect==1e-74
assert len(record.rounds[3].alignments[10].hsps)==1
assert record.rounds[3].alignments[11].hsps[0].score==696
assert record.rounds[3].alignments[11].hsps[0].bits==275
assert record.rounds[3].alignments[11].hsps[0].expect==3e-73
assert len(record.rounds[3].alignments[11].hsps)==1
assert record.rounds[3].alignments[12].hsps[0].score==656
assert record.rounds[3].alignments[12].hsps[0].bits==259
assert record.rounds[3].alignments[12].hsps[0].expect==1e-68
assert len(record.rounds[3].alignments[12].hsps)==1
assert record.rounds[3].alignments[13].hsps[0].score==652
assert record.rounds[3].alignments[13].hsps[0].bits==258
assert record.rounds[3].alignments[13].hsps[0].expect==4e-68
assert len(record.rounds[3].alignments[13].hsps)==1
assert record.rounds[3].alignments[14].hsps[0].score==211
assert record.rounds[3].alignments[14].hsps[0].bits==86.6
assert record.rounds[3].alignments[14].hsps[0].expect==2e-16
assert len(record.rounds[3].alignments[14].hsps)==1
assert record.rounds[3].alignments[15].hsps[0].score==209
assert record.rounds[3].alignments[15].hsps[0].bits==85.9
assert record.rounds[3].alignments[15].hsps[0].expect==3e-16
assert len(record.rounds[3].alignments[15].hsps)==1
assert record.rounds[3].alignments[16].hsps[0].score==208
assert record.rounds[3].alignments[16].hsps[0].bits==85.5
assert record.rounds[3].alignments[16].hsps[0].expect==4e-16
assert len(record.rounds[3].alignments[16].hsps)==1
assert record.rounds[3].alignments[17].hsps[0].score==208
assert record.rounds[3].alignments[17].hsps[0].bits==85.5
assert record.rounds[3].alignments[17].hsps[0].expect==4e-16
assert len(record.rounds[3].alignments[17].hsps)==1
assert record.rounds[3].alignments[18].hsps[0].score==189
assert record.rounds[3].alignments[18].hsps[0].bits==78.1
assert record.rounds[3].alignments[18].hsps[0].expect==7e-14
assert len(record.rounds[3].alignments[18].hsps)==1
assert record.rounds[3].alignments[19].hsps[0].score==158
assert record.rounds[3].alignments[19].hsps[0].bits==66.0
assert record.rounds[3].alignments[19].hsps[0].expect==3e-10
assert len(record.rounds[3].alignments[19].hsps)==1
assert record.rounds[3].alignments[20].hsps[0].score==120
assert record.rounds[3].alignments[20].hsps[0].bits==51.2
assert record.rounds[3].alignments[20].hsps[0].expect==9e-06
assert len(record.rounds[3].alignments[20].hsps)==1
assert record.rounds[3].alignments[21].hsps[0].score==105
assert record.rounds[3].alignments[21].hsps[0].bits==45.3
assert record.rounds[3].alignments[21].hsps[0].expect==5e-04
assert len(record.rounds[3].alignments[21].hsps)==1
assert record.rounds[3].alignments[22].hsps[0].score==73
assert record.rounds[3].alignments[22].hsps[0].bits==32.9
assert record.rounds[3].alignments[22].hsps[0].expect==2.9
assert len(record.rounds[3].alignments[22].hsps)==1
assert record.rounds[3].alignments[23].hsps[0].score==73
assert record.rounds[3].alignments[23].hsps[0].bits==32.9
assert record.rounds[3].alignments[23].hsps[0].expect==2.9
assert record.rounds[4].alignments[0].hsps[0].score==3106
assert record.rounds[4].alignments[0].hsps[0].bits==1214
assert record.rounds[4].alignments[0].hsps[0].expect==0.0
assert len(record.rounds[4].alignments[0].hsps)==1
assert record.rounds[4].alignments[1].hsps[0].score==1758
assert record.rounds[4].alignments[1].hsps[0].bits==689
assert record.rounds[4].alignments[1].hsps[0].expect==0.0
assert record.rounds[4].alignments[1].hsps[1].score==644
assert record.rounds[4].alignments[1].hsps[1].bits==255
assert record.rounds[4].alignments[1].hsps[1].expect==3e-67
assert len(record.rounds[4].alignments[1].hsps)==2
assert record.rounds[4].alignments[2].hsps[0].score==867
assert record.rounds[4].alignments[2].hsps[0].bits==342
assert record.rounds[4].alignments[2].hsps[0].expect==2e-93
assert len(record.rounds[4].alignments[2].hsps)==1
assert record.rounds[4].alignments[3].hsps[0].score==865
assert record.rounds[4].alignments[3].hsps[0].bits==341
assert record.rounds[4].alignments[3].hsps[0].expect==4e-93
assert len(record.rounds[4].alignments[3].hsps)==1
assert record.rounds[4].alignments[4].hsps[0].score==864
assert record.rounds[4].alignments[4].hsps[0].bits==340
assert record.rounds[4].alignments[4].hsps[0].expect==5e-93
assert len(record.rounds[4].alignments[4].hsps)==1
assert record.rounds[4].alignments[5].hsps[0].score==840
assert record.rounds[4].alignments[5].hsps[0].bits==331
assert record.rounds[4].alignments[5].hsps[0].expect==3e-90
assert len(record.rounds[4].alignments[5].hsps)==1
assert record.rounds[4].alignments[6].hsps[0].score==802
assert record.rounds[4].alignments[6].hsps[0].bits==316
assert record.rounds[4].alignments[6].hsps[0].expect==9e-86
assert len(record.rounds[4].alignments[6].hsps)==1
assert record.rounds[4].alignments[7].hsps[0].score==777
assert record.rounds[4].alignments[7].hsps[0].bits==307
assert record.rounds[4].alignments[7].hsps[0].expect==8e-83
assert len(record.rounds[4].alignments[7].hsps)==1
assert record.rounds[4].alignments[8].hsps[0].score==735
assert record.rounds[4].alignments[8].hsps[0].bits==290
assert record.rounds[4].alignments[8].hsps[0].expect==7e-78
assert len(record.rounds[4].alignments[8].hsps)==1
assert record.rounds[4].alignments[9].hsps[0].score==720
assert record.rounds[4].alignments[9].hsps[0].bits==284
assert record.rounds[4].alignments[9].hsps[0].expect==4e-76
assert len(record.rounds[4].alignments[9].hsps)==1
assert record.rounds[4].alignments[10].hsps[0].score==709
assert record.rounds[4].alignments[10].hsps[0].bits==280
assert record.rounds[4].alignments[10].hsps[0].expect==8e-75
assert len(record.rounds[4].alignments[10].hsps)==1
assert record.rounds[4].alignments[11].hsps[0].score==697
assert record.rounds[4].alignments[11].hsps[0].bits==275
assert record.rounds[4].alignments[11].hsps[0].expect==2e-73
assert len(record.rounds[4].alignments[11].hsps)==1
assert record.rounds[4].alignments[12].hsps[0].score==657
assert record.rounds[4].alignments[12].hsps[0].bits==260
assert record.rounds[4].alignments[12].hsps[0].expect==9e-69
assert len(record.rounds[4].alignments[12].hsps)==1
assert record.rounds[4].alignments[13].hsps[0].score==651
assert record.rounds[4].alignments[13].hsps[0].bits==258
assert record.rounds[4].alignments[13].hsps[0].expect==5e-68
assert len(record.rounds[4].alignments[13].hsps)==1
assert record.rounds[4].alignments[14].hsps[0].score==196
assert record.rounds[4].alignments[14].hsps[0].bits==80.8
assert record.rounds[4].alignments[14].hsps[0].expect==1e-14
assert len(record.rounds[4].alignments[14].hsps)==1
assert record.rounds[4].alignments[15].hsps[0].score==194
assert record.rounds[4].alignments[15].hsps[0].bits==80.0
assert record.rounds[4].alignments[15].hsps[0].expect==2e-14
assert len(record.rounds[4].alignments[15].hsps)==1
assert record.rounds[4].alignments[16].hsps[0].score==192
assert record.rounds[4].alignments[16].hsps[0].bits==79.2
assert record.rounds[4].alignments[16].hsps[0].expect==3e-14
assert len(record.rounds[4].alignments[16].hsps)==1
assert record.rounds[4].alignments[17].hsps[0].score==190
assert record.rounds[4].alignments[17].hsps[0].bits==78.5
assert record.rounds[4].alignments[17].hsps[0].expect==5e-14
assert len(record.rounds[4].alignments[17].hsps)==1
assert record.rounds[4].alignments[18].hsps[0].score==168
assert record.rounds[4].alignments[18].hsps[0].bits==69.9
assert record.rounds[4].alignments[18].hsps[0].expect==2e-11
assert len(record.rounds[4].alignments[18].hsps)==1
assert record.rounds[4].alignments[19].hsps[0].score==155
assert record.rounds[4].alignments[19].hsps[0].bits==64.8
assert record.rounds[4].alignments[19].hsps[0].expect==7e-10
assert len(record.rounds[4].alignments[19].hsps)==1
assert record.rounds[4].alignments[20].hsps[0].score==154
assert record.rounds[4].alignments[20].hsps[0].bits==64.4
assert record.rounds[4].alignments[20].hsps[0].expect==9e-10
assert len(record.rounds[4].alignments[20].hsps)==1
assert record.rounds[4].alignments[21].hsps[0].score==137
assert record.rounds[4].alignments[21].hsps[0].bits==57.8
assert record.rounds[4].alignments[21].hsps[0].expect==9e-08
assert len(record.rounds[4].alignments[21].hsps)==1
assert record.rounds[4].alignments[22].hsps[0].score==74
assert record.rounds[4].alignments[22].hsps[0].bits==33.3
assert record.rounds[4].alignments[22].hsps[0].expect==2.2
assert len(record.rounds[4].alignments[22].hsps)==1
assert record.rounds[4].alignments[23].hsps[0].score==73
assert record.rounds[4].alignments[23].hsps[0].bits==32.9
assert record.rounds[4].alignments[23].hsps[0].expect==2.9
assert len(record.rounds[4].alignments[23].hsps)==1
assert record.rounds[0].alignments[0].hsps[0].identities==(765, 889)
assert record.rounds[0].alignments[0].hsps[0].positives==(765, 889)
assert record.rounds[0].alignments[1].hsps[0].identities==(281, 634)
assert record.rounds[0].alignments[1].hsps[0].positives==(375, 634)
assert record.rounds[0].alignments[1].hsps[0].gaps==(32, 634)
assert record.rounds[0].alignments[1].hsps[1].identities==(69, 209)
assert record.rounds[0].alignments[1].hsps[1].positives==(107, 209)
assert record.rounds[0].alignments[1].hsps[1].gaps==(26, 209)
assert record.rounds[0].alignments[2].hsps[0].identities==(69, 267)
assert record.rounds[0].alignments[2].hsps[0].positives==(116, 267)
assert record.rounds[0].alignments[2].hsps[0].gaps==(22, 267)
assert record.rounds[0].alignments[3].hsps[0].identities==(71, 267)
assert record.rounds[0].alignments[3].hsps[0].positives==(122, 267)
assert record.rounds[0].alignments[3].hsps[0].gaps==(22, 267)
assert record.rounds[0].alignments[4].hsps[0].identities==(72, 268)
assert record.rounds[0].alignments[4].hsps[0].positives==(122, 268)
assert record.rounds[0].alignments[4].hsps[0].gaps==(23, 268)
assert record.rounds[0].alignments[5].hsps[0].identities==(73, 270)
assert record.rounds[0].alignments[5].hsps[0].positives==(120, 270)
assert record.rounds[0].alignments[5].hsps[0].gaps==(27, 270)
assert record.rounds[0].alignments[6].hsps[0].identities==(65, 250)
assert record.rounds[0].alignments[6].hsps[0].positives==(110, 250)
assert record.rounds[0].alignments[6].hsps[0].gaps==(19, 250)
assert record.rounds[0].alignments[7].hsps[0].identities==(68, 249)
assert record.rounds[0].alignments[7].hsps[0].positives==(117, 249)
assert record.rounds[0].alignments[7].hsps[0].gaps==(12, 249)
assert record.rounds[0].alignments[8].hsps[0].identities==(67, 252)
assert record.rounds[0].alignments[8].hsps[0].positives==(125, 252)
assert record.rounds[0].alignments[8].hsps[0].gaps==(21, 252)
assert record.rounds[0].alignments[9].hsps[0].identities==(67, 243)
assert record.rounds[0].alignments[9].hsps[0].positives==(114, 243)
assert record.rounds[0].alignments[9].hsps[0].gaps==(12, 243)
assert record.rounds[0].alignments[10].hsps[0].identities==(58, 264)
assert record.rounds[0].alignments[10].hsps[0].positives==(114, 264)
assert record.rounds[0].alignments[10].hsps[0].gaps==(16, 264)
assert record.rounds[0].alignments[11].hsps[0].identities==(51, 210)
assert record.rounds[0].alignments[11].hsps[0].positives==(92, 210)
assert record.rounds[0].alignments[11].hsps[0].gaps==(17, 210)
assert record.rounds[0].alignments[12].hsps[0].identities==(60, 255)
assert record.rounds[0].alignments[12].hsps[0].positives==(97, 255)
assert record.rounds[0].alignments[12].hsps[0].gaps==(21, 255)
assert record.rounds[0].alignments[13].hsps[0].identities==(63, 289)
assert record.rounds[0].alignments[13].hsps[0].positives==(116, 289)
assert record.rounds[0].alignments[13].hsps[0].gaps==(34, 289)
assert record.rounds[0].alignments[14].hsps[0].identities==(25, 70)
assert record.rounds[0].alignments[14].hsps[0].positives==(36, 70)
assert record.rounds[0].alignments[15].hsps[0].identities==(25, 70)
assert record.rounds[0].alignments[15].hsps[0].positives==(35, 70)
assert record.rounds[0].alignments[16].hsps[0].identities==(25, 70)
assert record.rounds[0].alignments[16].hsps[0].positives==(34, 70)
assert record.rounds[0].alignments[17].hsps[0].identities==(21, 59)
assert record.rounds[0].alignments[17].hsps[0].positives==(29, 59)
assert record.rounds[0].alignments[18].hsps[0].identities==(17, 57)
assert record.rounds[0].alignments[18].hsps[0].positives==(29, 57)
assert record.rounds[0].alignments[19].hsps[0].identities==(18, 56)
assert record.rounds[0].alignments[19].hsps[0].positives==(28, 56)
assert record.rounds[0].alignments[20].hsps[0].identities==(23, 87)
assert record.rounds[0].alignments[20].hsps[0].positives==(38, 87)
assert record.rounds[0].alignments[20].hsps[0].gaps==(1, 87)
assert record.rounds[0].alignments[21].hsps[0].identities==(18, 61)
assert record.rounds[0].alignments[21].hsps[0].positives==(26, 61)
assert record.rounds[0].alignments[22].hsps[0].identities==(29, 99)
assert record.rounds[0].alignments[22].hsps[0].positives==(40, 99)
assert record.rounds[0].alignments[22].hsps[0].gaps==(14, 99)
assert record.rounds[0].alignments[23].hsps[0].identities==(26, 100)
assert record.rounds[0].alignments[23].hsps[0].positives==(46, 100)
assert record.rounds[0].alignments[23].hsps[0].gaps==(3, 100)
assert record.rounds[0].alignments[24].hsps[0].identities==(23, 90)
assert record.rounds[0].alignments[24].hsps[0].positives==(44, 90)
assert record.rounds[0].alignments[24].hsps[0].gaps==(10, 90)
assert record.rounds[0].alignments[25].hsps[0].identities==(22, 74)
assert record.rounds[0].alignments[25].hsps[0].positives==(39, 74)
assert record.rounds[0].alignments[25].hsps[0].gaps==(10, 74)
assert record.rounds[0].alignments[26].hsps[0].identities==(24, 78)
assert record.rounds[0].alignments[26].hsps[0].positives==(35, 78)
assert record.rounds[0].alignments[26].hsps[0].gaps==(8, 78)
assert record.rounds[1].alignments[0].hsps[0].identities==(765, 889)
assert record.rounds[1].alignments[0].hsps[0].positives==(765, 889)
assert record.rounds[1].alignments[1].hsps[0].identities==(281, 635)
assert record.rounds[1].alignments[1].hsps[0].positives==(374, 635)
assert record.rounds[1].alignments[1].hsps[0].gaps==(34, 635)
assert record.rounds[1].alignments[1].hsps[1].identities==(69, 209)
assert record.rounds[1].alignments[1].hsps[1].positives==(107, 209)
assert record.rounds[1].alignments[1].hsps[1].gaps==(26, 209)
assert record.rounds[1].alignments[2].hsps[0].identities==(72, 268)
assert record.rounds[1].alignments[2].hsps[0].positives==(119, 268)
assert record.rounds[1].alignments[2].hsps[0].gaps==(23, 268)
assert record.rounds[1].alignments[3].hsps[0].identities==(69, 265)
assert record.rounds[1].alignments[3].hsps[0].positives==(115, 265)
assert record.rounds[1].alignments[3].hsps[0].gaps==(18, 265)
assert record.rounds[1].alignments[4].hsps[0].identities==(72, 268)
assert record.rounds[1].alignments[4].hsps[0].positives==(117, 268)
assert record.rounds[1].alignments[4].hsps[0].gaps==(23, 268)
assert record.rounds[1].alignments[5].hsps[0].identities==(69, 265)
assert record.rounds[1].alignments[5].hsps[0].positives==(116, 265)
assert record.rounds[1].alignments[5].hsps[0].gaps==(18, 265)
assert record.rounds[1].alignments[6].hsps[0].identities==(67, 249)
assert record.rounds[1].alignments[6].hsps[0].positives==(116, 249)
assert record.rounds[1].alignments[6].hsps[0].gaps==(12, 249)
assert record.rounds[1].alignments[7].hsps[0].identities==(67, 249)
assert record.rounds[1].alignments[7].hsps[0].positives==(115, 249)
assert record.rounds[1].alignments[7].hsps[0].gaps==(13, 249)
assert record.rounds[1].alignments[8].hsps[0].identities==(58, 271)
assert record.rounds[1].alignments[8].hsps[0].positives==(116, 271)
assert record.rounds[1].alignments[8].hsps[0].gaps==(16, 271)
assert record.rounds[1].alignments[9].hsps[0].identities==(65, 250)
assert record.rounds[1].alignments[9].hsps[0].positives==(109, 250)
assert record.rounds[1].alignments[9].hsps[0].gaps==(19, 250)
assert record.rounds[1].alignments[10].hsps[0].identities==(58, 286)
assert record.rounds[1].alignments[10].hsps[0].positives==(109, 286)
assert record.rounds[1].alignments[10].hsps[0].gaps==(28, 286)
assert record.rounds[1].alignments[11].hsps[0].identities==(58, 256)
assert record.rounds[1].alignments[11].hsps[0].positives==(96, 256)
assert record.rounds[1].alignments[11].hsps[0].gaps==(15, 256)
assert record.rounds[1].alignments[12].hsps[0].identities==(68, 259)
assert record.rounds[1].alignments[12].hsps[0].positives==(126, 259)
assert record.rounds[1].alignments[12].hsps[0].gaps==(21, 259)
assert record.rounds[1].alignments[13].hsps[0].identities==(60, 277)
assert record.rounds[1].alignments[13].hsps[0].positives==(109, 277)
assert record.rounds[1].alignments[13].hsps[0].gaps==(29, 277)
assert record.rounds[1].alignments[14].hsps[0].identities==(25, 70)
assert record.rounds[1].alignments[14].hsps[0].positives==(35, 70)
assert record.rounds[1].alignments[15].hsps[0].identities==(25, 70)
assert record.rounds[1].alignments[15].hsps[0].positives==(34, 70)
assert record.rounds[1].alignments[16].hsps[0].identities==(25, 70)
assert record.rounds[1].alignments[16].hsps[0].positives==(36, 70)
assert record.rounds[1].alignments[17].hsps[0].identities==(22, 72)
assert record.rounds[1].alignments[17].hsps[0].positives==(30, 72)
assert record.rounds[1].alignments[18].hsps[0].identities==(18, 62)
assert record.rounds[1].alignments[18].hsps[0].positives==(26, 62)
assert record.rounds[1].alignments[19].hsps[0].identities==(18, 63)
assert record.rounds[1].alignments[19].hsps[0].positives==(28, 63)
assert record.rounds[1].alignments[20].hsps[0].identities==(17, 60)
assert record.rounds[1].alignments[20].hsps[0].positives==(29, 60)
assert record.rounds[1].alignments[21].hsps[0].identities==(20, 49)
assert record.rounds[1].alignments[21].hsps[0].positives==(23, 49)
assert record.rounds[1].alignments[21].hsps[0].gaps==(6, 49)
assert record.rounds[1].alignments[22].hsps[0].identities==(17, 92)
assert record.rounds[1].alignments[22].hsps[0].positives==(25, 92)
assert record.rounds[1].alignments[22].hsps[0].gaps==(13, 92)
assert record.rounds[1].alignments[23].hsps[0].identities==(22, 95)
assert record.rounds[1].alignments[23].hsps[0].positives==(35, 95)
assert record.rounds[1].alignments[24].hsps[0].identities==(26, 238)
assert record.rounds[1].alignments[24].hsps[0].positives==(59, 238)
assert record.rounds[1].alignments[24].hsps[0].gaps==(45, 238)
assert record.rounds[1].alignments[25].hsps[0].identities==(20, 117)
assert record.rounds[1].alignments[25].hsps[0].positives==(44, 117)
assert record.rounds[1].alignments[25].hsps[0].gaps==(18, 117)
assert record.rounds[2].alignments[0].hsps[0].identities==(765, 889)
assert record.rounds[2].alignments[0].hsps[0].positives==(765, 889)
assert record.rounds[2].alignments[1].hsps[0].identities==(281, 635)
assert record.rounds[2].alignments[1].hsps[0].positives==(374, 635)
assert record.rounds[2].alignments[1].hsps[0].gaps==(34, 635)
assert record.rounds[2].alignments[1].hsps[1].identities==(69, 209)
assert record.rounds[2].alignments[1].hsps[1].positives==(107, 209)
assert record.rounds[2].alignments[1].hsps[1].gaps==(26, 209)
assert record.rounds[2].alignments[2].hsps[0].identities==(69, 265)
assert record.rounds[2].alignments[2].hsps[0].positives==(115, 265)
assert record.rounds[2].alignments[2].hsps[0].gaps==(18, 265)
assert record.rounds[2].alignments[3].hsps[0].identities==(72, 268)
assert record.rounds[2].alignments[3].hsps[0].positives==(119, 268)
assert record.rounds[2].alignments[3].hsps[0].gaps==(23, 268)
assert record.rounds[2].alignments[4].hsps[0].identities==(72, 268)
assert record.rounds[2].alignments[4].hsps[0].positives==(117, 268)
assert record.rounds[2].alignments[4].hsps[0].gaps==(23, 268)
assert record.rounds[2].alignments[5].hsps[0].identities==(69, 265)
assert record.rounds[2].alignments[5].hsps[0].positives==(116, 265)
assert record.rounds[2].alignments[5].hsps[0].gaps==(18, 265)
assert record.rounds[2].alignments[6].hsps[0].identities==(67, 249)
assert record.rounds[2].alignments[6].hsps[0].positives==(116, 249)
assert record.rounds[2].alignments[6].hsps[0].gaps==(12, 249)
assert record.rounds[2].alignments[7].hsps[0].identities==(67, 249)
assert record.rounds[2].alignments[7].hsps[0].positives==(114, 249)
assert record.rounds[2].alignments[7].hsps[0].gaps==(13, 249)
assert record.rounds[2].alignments[8].hsps[0].identities==(58, 271)
assert record.rounds[2].alignments[8].hsps[0].positives==(116, 271)
assert record.rounds[2].alignments[8].hsps[0].gaps==(16, 271)
assert record.rounds[2].alignments[9].hsps[0].identities==(58, 286)
assert record.rounds[2].alignments[9].hsps[0].positives==(109, 286)
assert record.rounds[2].alignments[9].hsps[0].gaps==(28, 286)
assert record.rounds[2].alignments[10].hsps[0].identities==(65, 250)
assert record.rounds[2].alignments[10].hsps[0].positives==(109, 250)
assert record.rounds[2].alignments[10].hsps[0].gaps==(19, 250)
assert record.rounds[2].alignments[11].hsps[0].identities==(58, 256)
assert record.rounds[2].alignments[11].hsps[0].positives==(96, 256)
assert record.rounds[2].alignments[11].hsps[0].gaps==(15, 256)
assert record.rounds[2].alignments[12].hsps[0].identities==(60, 277)
assert record.rounds[2].alignments[12].hsps[0].positives==(109, 277)
assert record.rounds[2].alignments[12].hsps[0].gaps==(29, 277)
assert record.rounds[2].alignments[13].hsps[0].identities==(68, 259)
assert record.rounds[2].alignments[13].hsps[0].positives==(125, 259)
assert record.rounds[2].alignments[13].hsps[0].gaps==(21, 259)
assert record.rounds[2].alignments[14].hsps[0].identities==(22, 72)
assert record.rounds[2].alignments[14].hsps[0].positives==(30, 72)
assert record.rounds[2].alignments[15].hsps[0].identities==(25, 71)
assert record.rounds[2].alignments[15].hsps[0].positives==(35, 71)
assert record.rounds[2].alignments[16].hsps[0].identities==(25, 71)
assert record.rounds[2].alignments[16].hsps[0].positives==(34, 71)
assert record.rounds[2].alignments[17].hsps[0].identities==(25, 71)
assert record.rounds[2].alignments[17].hsps[0].positives==(36, 71)
assert record.rounds[2].alignments[18].hsps[0].identities==(18, 63)
assert record.rounds[2].alignments[18].hsps[0].positives==(28, 63)
assert record.rounds[2].alignments[19].hsps[0].identities==(22, 80)
assert record.rounds[2].alignments[19].hsps[0].positives==(30, 80)
assert record.rounds[2].alignments[19].hsps[0].gaps==(6, 80)
assert record.rounds[2].alignments[20].hsps[0].identities==(17, 62)
assert record.rounds[2].alignments[20].hsps[0].positives==(29, 62)
assert record.rounds[2].alignments[21].hsps[0].identities==(20, 49)
assert record.rounds[2].alignments[21].hsps[0].positives==(23, 49)
assert record.rounds[2].alignments[21].hsps[0].gaps==(6, 49)
assert record.rounds[2].alignments[22].hsps[0].identities==(22, 95)
assert record.rounds[2].alignments[22].hsps[0].positives==(35, 95)
assert record.rounds[2].alignments[23].hsps[0].identities==(17, 92)
assert record.rounds[2].alignments[23].hsps[0].positives==(25, 92)
assert record.rounds[2].alignments[23].hsps[0].gaps==(13, 92)
assert record.rounds[3].alignments[0].hsps[0].identities==(765, 889)
assert record.rounds[3].alignments[0].hsps[0].positives==(765, 889)
assert record.rounds[3].alignments[1].hsps[0].identities==(281, 635)
assert record.rounds[3].alignments[1].hsps[0].positives==(374, 635)
assert record.rounds[3].alignments[1].hsps[0].gaps==(34, 635)
assert record.rounds[3].alignments[1].hsps[1].identities==(69, 209)
assert record.rounds[3].alignments[1].hsps[1].positives==(107, 209)
assert record.rounds[3].alignments[1].hsps[1].gaps==(26, 209)
assert record.rounds[3].alignments[2].hsps[0].identities==(69, 265)
assert record.rounds[3].alignments[2].hsps[0].positives==(115, 265)
assert record.rounds[3].alignments[2].hsps[0].gaps==(18, 265)
assert record.rounds[3].alignments[3].hsps[0].identities==(71, 268)
assert record.rounds[3].alignments[3].hsps[0].positives==(119, 268)
assert record.rounds[3].alignments[3].hsps[0].gaps==(23, 268)
assert record.rounds[3].alignments[4].hsps[0].identities==(71, 268)
assert record.rounds[3].alignments[4].hsps[0].positives==(117, 268)
assert record.rounds[3].alignments[4].hsps[0].gaps==(23, 268)
assert record.rounds[3].alignments[5].hsps[0].identities==(69, 265)
assert record.rounds[3].alignments[5].hsps[0].positives==(116, 265)
assert record.rounds[3].alignments[5].hsps[0].gaps==(18, 265)
assert record.rounds[3].alignments[6].hsps[0].identities==(67, 249)
assert record.rounds[3].alignments[6].hsps[0].positives==(116, 249)
assert record.rounds[3].alignments[6].hsps[0].gaps==(12, 249)
assert record.rounds[3].alignments[7].hsps[0].identities==(67, 249)
assert record.rounds[3].alignments[7].hsps[0].positives==(114, 249)
assert record.rounds[3].alignments[7].hsps[0].gaps==(13, 249)
assert record.rounds[3].alignments[8].hsps[0].identities==(58, 271)
assert record.rounds[3].alignments[8].hsps[0].positives==(116, 271)
assert record.rounds[3].alignments[8].hsps[0].gaps==(16, 271)
assert record.rounds[3].alignments[9].hsps[0].identities==(58, 286)
assert record.rounds[3].alignments[9].hsps[0].positives==(109, 286)
assert record.rounds[3].alignments[9].hsps[0].gaps==(28, 286)
assert record.rounds[3].alignments[10].hsps[0].identities==(65, 250)
assert record.rounds[3].alignments[10].hsps[0].positives==(109, 250)
assert record.rounds[3].alignments[10].hsps[0].gaps==(19, 250)
assert record.rounds[3].alignments[11].hsps[0].identities==(58, 256)
assert record.rounds[3].alignments[11].hsps[0].positives==(96, 256)
assert record.rounds[3].alignments[11].hsps[0].gaps==(15, 256)
assert record.rounds[3].alignments[12].hsps[0].identities==(60, 277)
assert record.rounds[3].alignments[12].hsps[0].positives==(109, 277)
assert record.rounds[3].alignments[12].hsps[0].gaps==(29, 277)
assert record.rounds[3].alignments[13].hsps[0].identities==(68, 259)
assert record.rounds[3].alignments[13].hsps[0].positives==(125, 259)
assert record.rounds[3].alignments[13].hsps[0].gaps==(21, 259)
assert record.rounds[3].alignments[14].hsps[0].identities==(25, 71)
assert record.rounds[3].alignments[14].hsps[0].positives==(35, 71)
assert record.rounds[3].alignments[15].hsps[0].identities==(22, 72)
assert record.rounds[3].alignments[15].hsps[0].positives==(30, 72)
assert record.rounds[3].alignments[16].hsps[0].identities==(25, 71)
assert record.rounds[3].alignments[16].hsps[0].positives==(34, 71)
assert record.rounds[3].alignments[17].hsps[0].identities==(25, 71)
assert record.rounds[3].alignments[17].hsps[0].positives==(36, 71)
assert record.rounds[3].alignments[18].hsps[0].identities==(22, 80)
assert record.rounds[3].alignments[18].hsps[0].positives==(30, 80)
assert record.rounds[3].alignments[18].hsps[0].gaps==(6, 80)
assert record.rounds[3].alignments[19].hsps[0].identities==(18, 63)
assert record.rounds[3].alignments[19].hsps[0].positives==(28, 63)
assert record.rounds[3].alignments[20].hsps[0].identities==(17, 62)
assert record.rounds[3].alignments[20].hsps[0].positives==(29, 62)
assert record.rounds[3].alignments[21].hsps[0].identities==(20, 49)
assert record.rounds[3].alignments[21].hsps[0].positives==(23, 49)
assert record.rounds[3].alignments[21].hsps[0].gaps==(6, 49)
assert record.rounds[3].alignments[22].hsps[0].identities==(17, 92)
assert record.rounds[3].alignments[22].hsps[0].positives==(25, 92)
assert record.rounds[3].alignments[22].hsps[0].gaps==(13, 92)
assert record.rounds[3].alignments[23].hsps[0].identities==(22, 95)
assert record.rounds[3].alignments[23].hsps[0].positives==(35, 95)
assert record.rounds[4].alignments[0].hsps[0].identities==(765, 889)
assert record.rounds[4].alignments[0].hsps[0].positives==(765, 889)
assert record.rounds[4].alignments[1].hsps[0].identities==(281, 635)
assert record.rounds[4].alignments[1].hsps[0].positives==(374, 635)
assert record.rounds[4].alignments[1].hsps[0].gaps==(34, 635)
assert record.rounds[4].alignments[1].hsps[1].identities==(69, 209)
assert record.rounds[4].alignments[1].hsps[1].positives==(107, 209)
assert record.rounds[4].alignments[1].hsps[1].gaps==(26, 209)
assert record.rounds[4].alignments[2].hsps[0].identities==(69, 265)
assert record.rounds[4].alignments[2].hsps[0].positives==(115, 265)
assert record.rounds[4].alignments[2].hsps[0].gaps==(18, 265)
assert record.rounds[4].alignments[3].hsps[0].identities==(71, 268)
assert record.rounds[4].alignments[3].hsps[0].positives==(119, 268)
assert record.rounds[4].alignments[3].hsps[0].gaps==(23, 268)
assert record.rounds[4].alignments[4].hsps[0].identities==(71, 268)
assert record.rounds[4].alignments[4].hsps[0].positives==(117, 268)
assert record.rounds[4].alignments[4].hsps[0].gaps==(23, 268)
assert record.rounds[4].alignments[5].hsps[0].identities==(69, 265)
assert record.rounds[4].alignments[5].hsps[0].positives==(116, 265)
assert record.rounds[4].alignments[5].hsps[0].gaps==(18, 265)
assert record.rounds[4].alignments[6].hsps[0].identities==(67, 249)
assert record.rounds[4].alignments[6].hsps[0].positives==(116, 249)
assert record.rounds[4].alignments[6].hsps[0].gaps==(12, 249)
assert record.rounds[4].alignments[7].hsps[0].identities==(67, 249)
assert record.rounds[4].alignments[7].hsps[0].positives==(114, 249)
assert record.rounds[4].alignments[7].hsps[0].gaps==(13, 249)
assert record.rounds[4].alignments[8].hsps[0].identities==(58, 271)
assert record.rounds[4].alignments[8].hsps[0].positives==(116, 271)
assert record.rounds[4].alignments[8].hsps[0].gaps==(16, 271)
assert record.rounds[4].alignments[9].hsps[0].identities==(58, 286)
assert record.rounds[4].alignments[9].hsps[0].positives==(109, 286)
assert record.rounds[4].alignments[9].hsps[0].gaps==(28, 286)
assert record.rounds[4].alignments[10].hsps[0].identities==(64, 250)
assert record.rounds[4].alignments[10].hsps[0].positives==(107, 250)
assert record.rounds[4].alignments[10].hsps[0].gaps==(19, 250)
assert record.rounds[4].alignments[11].hsps[0].identities==(58, 256)
assert record.rounds[4].alignments[11].hsps[0].positives==(96, 256)
assert record.rounds[4].alignments[11].hsps[0].gaps==(15, 256)
assert record.rounds[4].alignments[12].hsps[0].identities==(60, 277)
assert record.rounds[4].alignments[12].hsps[0].positives==(109, 277)
assert record.rounds[4].alignments[12].hsps[0].gaps==(29, 277)
assert record.rounds[4].alignments[13].hsps[0].identities==(66, 258)
assert record.rounds[4].alignments[13].hsps[0].positives==(120, 258)
assert record.rounds[4].alignments[13].hsps[0].gaps==(19, 258)
assert record.rounds[4].alignments[14].hsps[0].identities==(22, 73)
assert record.rounds[4].alignments[14].hsps[0].positives==(30, 73)
assert record.rounds[4].alignments[15].hsps[0].identities==(25, 73)
assert record.rounds[4].alignments[15].hsps[0].positives==(35, 73)
assert record.rounds[4].alignments[16].hsps[0].identities==(25, 73)
assert record.rounds[4].alignments[16].hsps[0].positives==(36, 73)
assert record.rounds[4].alignments[17].hsps[0].identities==(25, 71)
assert record.rounds[4].alignments[17].hsps[0].positives==(34, 71)
assert record.rounds[4].alignments[18].hsps[0].identities==(22, 80)
assert record.rounds[4].alignments[18].hsps[0].positives==(30, 80)
assert record.rounds[4].alignments[18].hsps[0].gaps==(6, 80)
assert record.rounds[4].alignments[19].hsps[0].identities==(18, 63)
assert record.rounds[4].alignments[19].hsps[0].positives==(28, 63)
assert record.rounds[4].alignments[20].hsps[0].identities==(17, 62)
assert record.rounds[4].alignments[20].hsps[0].positives==(29, 62)
assert record.rounds[4].alignments[21].hsps[0].identities==(20, 49)
assert record.rounds[4].alignments[21].hsps[0].positives==(23, 49)
assert record.rounds[4].alignments[21].hsps[0].gaps==(6, 49)
assert record.rounds[4].alignments[22].hsps[0].identities==(17, 92)
assert record.rounds[4].alignments[22].hsps[0].positives==(25, 92)
assert record.rounds[4].alignments[22].hsps[0].gaps==(13, 92)
assert record.rounds[4].alignments[23].hsps[0].identities==(22, 95)
assert record.rounds[4].alignments[23].hsps[0].positives==(35, 95)
assert record.rounds[0].alignments[0].hsps[0].query=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[0].alignments[0].hsps[0].match=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQ            GI            TSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQ                                      QLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGM                   RARVEEAMELIVAHGRENGA            EKVQENRKRHGSSRSVV            NAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP                   NFGFV"
assert record.rounds[0].alignments[0].hsps[0].sbjct=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQLERKRTRELLELGIDNEDSEHENDDDTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQRLRNRGERDREREREREMRRSSGLRAGSRRDRDRDFRRQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMLLELSPAQLLLLLASEDSLRARVEEAMELIVAHGRENGADSILDLGLLDSSEKVQENRKRHGSSRSVVDMDLDDTDDGDDNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPLYSSKQILKQKLLLAIKTKNFGFV"
assert record.rounds[0].alignments[0].hsps[0].query_start==1
assert record.rounds[0].alignments[0].hsps[0].query_end==889
assert record.rounds[0].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[0].alignments[0].hsps[0].sbjct_end==889
assert record.rounds[0].alignments[1].hsps[0].query=="QCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSL-EVDRDRDLLIQQTMRQLNNHFGR--RCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGN---PSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVA-EQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[0].alignments[1].hsps[0].match=="Q +F G+   + +L  +  L+L+LFGR FM+DVG E GS+L EL GF VKE +FRR MEKLRN Q RDL L +++R+R+ LI QT ++LN  FG   R    P+  +RVKVTFKDEPGEGSGVARSFYT+IA+A L++ K+PNL+ +Q      H+  +                                        L++D RP+ P +  +   P    D L  H Q +GERLYP++ ++    A KITGM                   R +V EA+E+I    +   +               Q ++ +   S  VV            N PLFY PGKRGFY PR G  +  R+N FRNIGR++GLCLLQNEL P+ L RHV+K +L RK+ +HD AFFDP +YES RQ+I  +Q+ + +   + M+L F +DL KEEG G  ELIP G ++ V+ + +Y       ++   + + E+ L A++ G+ DVLP NS+ +LTAED RLL+NG G++NV  LIS+T+FNDES E  +KL +FK+WFWSIVE+M++ ERQ LVYFWT SP+LPASEEGFQP+PS+TIRP DD HLPTANTCISRLY+P                   NFGFV"
assert record.rounds[0].alignments[1].hsps[0].sbjct=="QSNFSGIAPCNFLLSAK--LTLDLFGRPFMDDVGMEHGSVLPELRGFPVKEMRFRRHMEKLRNGQQRDLVLCKLERNRESLIVQTFKELNTQFGNQSRRIQPPITFNRVKVTFKDEPGEGSGVARSFYTSIAEALLASAKIPNLESVQVGTN--HSKYVVPFSSILRTTVSGSSRDQSTLQRRGSNSKILWRSARERKALNLDARPYTPPNSSDNATPESLNDHLSVHLQQIGERLYPKIHSINQTHAPKITGMLLEIPTPQLLSVISSDETLRQKVNEAIEIITFKQKSETSA--------------QSSQPKKSPSVVVVDPVDDD------NEPLFYSPGKRGFYIPRQGFASFERINAFRNIGRLIGLCLLQNELLPLFLQRHVLKYILRRKIKFHDLAFFDPALYESFRQIIQNAQTKEGEETINRMELCFVIDLMKEEGCGNRELIPGGRDV-VSHRVIYSSTSDAIQNIDXIKSQEKALEALKDGVFDVLPDNSMINLTAEDLRLLLNGVGDINVSTLISYTTFNDESSEGPDKLXKFKKWFWSIVEKMNIMERQHLVYFWTGSPALPASEEGFQPLPSVTIRPADDSHLPTANTCISRLYIPLYSSKSILRSKMLMAIKSKNFGFV"
assert record.rounds[0].alignments[1].hsps[0].query_start==263
assert record.rounds[0].alignments[1].hsps[0].query_end==889
assert record.rounds[0].alignments[1].hsps[0].sbjct_start==2287
assert record.rounds[0].alignments[1].hsps[0].sbjct_end==2895
assert record.rounds[0].alignments[1].hsps[1].query=="SARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDE------SLPA--ETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMD---RNCLEVLPTKMSYAANLK"
assert record.rounds[0].alignments[1].hsps[1].match=="++R DF  Y LSLMRSH  EH D LPVLD+ +L+H+AYV  A +Y+++  N     D                I              + A L + + +      S+P+  +  + H FF RS+S   LGC  P  FE+PL  A+PLAD+PHLLQPN+++++LF      + +++  SG      T D    +  +  PT++ ++ +LK"
assert record.rounds[0].alignments[1].hsps[1].sbjct=="NSRRDFFTYCLSLMRSHTSEHRDALPVLDITALRHIAYVLDAFVYYMR--NDSGFYDKQDSISGR--------INNLSPMTESYDTDDELANLEEFNADVQMSASSMPSGSQGTRRHAFFARSESTLTLGCSAPEGFELPLDMAMPLADKPHLLQPNSKRQELFANLPLLVTTNANNSG-----ATNDGDGASIFDYTPTRLGFSNSLK"
assert record.rounds[0].alignments[1].hsps[1].query_start==3
assert record.rounds[0].alignments[1].hsps[1].query_end==200
assert record.rounds[0].alignments[1].hsps[1].sbjct_start==1913
assert record.rounds[0].alignments[1].hsps[1].sbjct_end==2106
assert record.rounds[0].alignments[2].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQ-------PMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[2].hsps[0].match=="P  G N E  LN F+ IGR++GL +               K++L +KV   D    D   Y SL  ++    +   D  FS  D  F   +        ++L PNG NI VT +N  EYV      R+    E+  +A  +G  +++P+  +      +  LL+ G  E++++     T +   S  +     Q  +WFW +++  S  ++  L+ F T +  +P +  GF+       P      +  +   LP A+TC +RL +P"
assert record.rounds[0].alignments[2].hsps[0].sbjct=="PHSGINPE-HLNYFKFIGRVIGLAIFHRRFVDAFFVVSFYKMILQKKVTLQDMESMDAEYYRSLVWILDNDITGVLDLTFSVEDNCFGEVVT-------IDLKPNGRNIEVTEENKREYVDLVTVWRIQKRIEEQFNAFHEGFSELIPQELINVFDERELELLIGGISEIDMEDWKKHTDYRSYSEND-----QIIKWFWELMDEWSNEKKSRLLQFTTGTSRIPVN--GFKDLQGSDGPRKFTIEKAGEPNKLPKAHTCFNRLDLP"
assert record.rounds[0].alignments[2].hsps[0].query_start==606
assert record.rounds[0].alignments[2].hsps[0].query_end==865
assert record.rounds[0].alignments[2].hsps[0].sbjct_start==491
assert record.rounds[0].alignments[2].hsps[0].sbjct_end==742
assert record.rounds[0].alignments[3].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPS------ITIRPPDD-QHLPTANTCISRLYVP"
assert record.rounds[0].alignments[3].hsps[0].match=="P  G N E  LN F+ IGR++GL +             + K++L +KV   D    D  +Y SL  ++  S     D  FSA D  F   +        V+L P+G NI VT  N  EYV  Y + R++   ++   A   G  +++P++ +      +  LL+ G  E++++     T +     + +++++Q   WFW  V      +R  L+ F T +  +P +  GF+ +         TI    + Q LP ++TC +R+ +P"
assert record.rounds[0].alignments[3].hsps[0].sbjct=="PNSGINPE-HLNYFKFIGRVVGLGVFHRRFLDAFFVGALYKMMLRKKVVLQDMEGVDAEVYNSLNWMLENSIDGVLDLTFSADDERFGEVVT-------VDLKPDGRNIEVTDGNKKEYVELYTQWRIVDRVQEQFKAFMDGFNELIPEDLVTVFDERELELLIGGIAEIDIEDWKKHTDY--RGYQESDEVIQ---WFWKCVSEWDNEQRARLLQFTTGTSRIPVN--GFKDLQGSDGPRRFTIEKAGEVQQLPKSHTCFNRVDLP"
assert record.rounds[0].alignments[3].hsps[0].query_start==606
assert record.rounds[0].alignments[3].hsps[0].query_end==865
assert record.rounds[0].alignments[3].hsps[0].sbjct_start==533
assert record.rounds[0].alignments[3].hsps[0].sbjct_end==784
assert record.rounds[0].alignments[4].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[4].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++     +D     + +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N + +    +WFW  V  M   +R  L+ F T +  +P    A   G     S T+      + LP A+TC +RL +P"
assert record.rounds[0].alignments[4].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLHKPITLHDMESVDSEYYNSLRWIL----ENDP----TELDLRFIID---EELFGQTHQHELKNGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSANHQVI----QWFWKAVLMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPEKLPRAHTCFNRLDLP"
assert record.rounds[0].alignments[4].hsps[0].query_start==606
assert record.rounds[0].alignments[4].hsps[0].query_end==865
assert record.rounds[0].alignments[4].hsps[0].sbjct_start==649
assert record.rounds[0].alignments[4].hsps[0].sbjct_end==901
assert record.rounds[0].alignments[5].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIR---PPDDQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[5].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N + +     WFW  V  M   +R  L+ F T +  +P    A   G     S T+     PD   LP A+TC +RL +P"
assert record.rounds[0].alignments[5].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLQKLITLHDMESVDSEYYSSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKTGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSMNHQVI----HWFWKAVWMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPD--KLPRAHTCFNRLDLP"
assert record.rounds[0].alignments[5].hsps[0].query_start==606
assert record.rounds[0].alignments[5].hsps[0].query_end==865
assert record.rounds[0].alignments[5].hsps[0].sbjct_start==679
assert record.rounds[0].alignments[5].hsps[0].sbjct_end==931
assert record.rounds[0].alignments[6].hsps[0].query=="GRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRPPD--DQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[6].hsps[0].match=="GR +   ++   L      R   K +LG+ V + D    D   Y+ L  L+      + D      DL F+ ++ +E G  +V +L PNG NI VT +N  EYV    + RM     + L A  +G  +++PK  +   T ++  LL  G   +++  L S T ++     + +      +WFW  +      +R   + F T +  +P    A+ EG   +    I   D     LP+A+TC ++L +P"
assert record.rounds[0].alignments[6].hsps[0].sbjct=="GRYVAKAVMTTALLECYFTRSFYKHILGKSVRYTDMESEDYHFYQGLVYLL------ENDVSTLGYDLTFSTEV-QEFGVCEVRDLKPNGANILVTEENKKEYVHLVCQMRMTGAIRKQLAAFLEGFYEIIPKRLISIFTEQELELLYTGLPTIDIDDLKSNTEYHKYQSNSIQ-----IQWFWRALRSFDQADRAKFLQFVTGTSKVPLQGFAALEGMNGIQKFQIHRDDRSTDRLPSAHTCFNQLDLP"
assert record.rounds[0].alignments[6].hsps[0].query_start==623
assert record.rounds[0].alignments[6].hsps[0].query_end==865
assert record.rounds[0].alignments[6].hsps[0].sbjct_start==45
assert record.rounds[0].alignments[6].hsps[0].sbjct_end==282
assert record.rounds[0].alignments[7].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDL-TAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[7].hsps[0].match=="F  IG +LGL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    + + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+ L   E+  LL+ G   ++ Q L   T +  + G   + +L   R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[0].alignments[7].hsps[0].sbjct=="FTLIGIVLGLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGNVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVNLYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRDSVL--IREFWEIVHSFTDEQKRLFLQFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[0].alignments[7].hsps[0].query_start==619
assert record.rounds[0].alignments[7].hsps[0].query_end==865
assert record.rounds[0].alignments[7].hsps[0].sbjct_start==612
assert record.rounds[0].alignments[7].hsps[0].sbjct_end==850
assert record.rounds[0].alignments[8].hsps[0].query=="IGRILGLCLLQNELCPITLNRHVIKVLL----GRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGE-VNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPM-PSITIRPPDDQ--HLPTANTCISRLYVP"
assert record.rounds[0].alignments[8].hsps[0].match=="+G+++G CL ++ L  ++     +K LL    G   ++ D   +D V+Y +L +L+  + ++D      ++DL F +D   E     V+LIPNG    VT  NV  YV K  ++++     +P+ A   GL  ++  + +E   + + ++L++G  + +++  L S T +     E+     Q    FW ++      E+ + + F TS P  P   +GF+ + P   IR    +   LPTA+TC++ L +P"
assert record.rounds[0].alignments[8].hsps[0].sbjct=="LGKVVGKCLYEHVLIDVSFADFFLKKLLNYSNGFLSSFSDLGSYDSVLYNNLIKLL--NMTTDE---IKSLDLTFEIDE-PESSAKVVDLIPNGSKTYVTKDNVLLYVTKVTDYKLNKRCFKPVSAFHGGLSVIIAPHWMEMFNSIELQMLISGERDNIDLDDLKSNTEYGGYKEED-----QTIVDFWEVLNEFKFEEKLNFLKFVTSVPQAPL--QGFKALDPKFGIRNAGTEKYRLPTASTCVNLLKLP"
assert record.rounds[0].alignments[8].hsps[0].query_start==622
assert record.rounds[0].alignments[8].hsps[0].query_end==865
assert record.rounds[0].alignments[8].hsps[0].sbjct_start==647
assert record.rounds[0].alignments[8].hsps[0].sbjct_end==885
assert record.rounds[0].alignments[9].hsps[0].query=="ILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDL-TAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[9].hsps[0].match=="ILGL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    S + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+ L   E+  LL+ G   ++ Q L   T +  + G   E ++   R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[0].alignments[9].hsps[0].sbjct=="ILGLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGSVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVISYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRESVV--IREFWEIVHSFTDEQKRLFLLFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[0].alignments[9].hsps[0].query_start==625
assert record.rounds[0].alignments[9].hsps[0].query_end==865
assert record.rounds[0].alignments[9].hsps[0].sbjct_start==628
assert record.rounds[0].alignments[9].hsps[0].sbjct_end==860
assert record.rounds[0].alignments[10].hsps[0].query=="GKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAV------DLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPK-NSLEDLTAEDFRLLVNGCGE---VNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRL"
assert record.rounds[0].alignments[10].hsps[0].match=="GKN++  L  +   G ++GL +  + +  +   + + K L    +++ D++   P    +L +++  ++ +  D      +  +        D    +    VEL  NG N+P+T  N +E+V K+ E  +    E   +    G   V  + NS++   +E+   LV G  E    + + L S T +     +++  +     WFW I+E      ++ L+ F T+S  +PA+     P     +   D   LP A+TC + +"
assert record.rounds[0].alignments[10].hsps[0].sbjct=="GKNSQLEL--YYLFGVVMGLAIFNSTILDLQFPKALYKKLCSEPLSFEDYSELFPETSRNLIKMLNYTEDNFEDVFSLTFETTYRNNNWILNDSKSSKEYVTVELCENGRNVPITQSNKHEFVMKWVEFYLEKSIEPQYNKFVSGFKRVFAECNSIKLFNSEELERLVCGDEEQTKFDFKSLRSVTKYVGGFSDDSRAVC----WFWEIIESWDYPLQKKLLQFVTASDRIPATGISTIPFKISLLGSHDSDDLPLAHTCFNEI"
assert record.rounds[0].alignments[10].hsps[0].query_start==609
assert record.rounds[0].alignments[10].hsps[0].query_end==862
assert record.rounds[0].alignments[10].hsps[0].sbjct_start==607
assert record.rounds[0].alignments[10].hsps[0].sbjct_end==864
assert record.rounds[0].alignments[11].hsps[0].query=="DPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP-ASEEGFQPMPSITIRPPD-----DQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[11].hsps[0].match=="DP++ +SL+ ++    + D +    ++ L F V      G   +ELIP G N  +   NV EY+    +  +    E+ L A  +G   V     +  L  ++  + + G  E +  M   +T+ N E G   +  +     F SI+      ER+  + F T SP LP    +   P  ++ ++  +     D++LP+  TC + L +P"
assert record.rounds[0].alignments[11].hsps[0].sbjct=="DPLLAKSLKYIVA---NKDDNMTLESLSLTFTVP-----GNDDIELIPGGCNKSLNSSNVEEYIHGVIDQILGKGIEKQLKAFIEGFSKVFSYERMLILFPDEL-VDIFGRVEEDWSMATLYTNLNAEHGYTMDSSIIHD--FISIISAFGKHERRLFLQFLTGSPKLPIGGFKSLNPKFTVVLKHAEDGLTADEYLPSVMTCANYLKLP"
assert record.rounds[0].alignments[11].hsps[0].query_start==662
assert record.rounds[0].alignments[11].hsps[0].query_end==865
assert record.rounds[0].alignments[11].hsps[0].sbjct_start==1259
assert record.rounds[0].alignments[11].hsps[0].sbjct_end==1457
assert record.rounds[0].alignments[12].hsps[0].query=="NCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGG--GQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITI----RPPDDQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[12].hsps[0].match=="N F  IG   GL +  + +  +     + K LL  K    D     P    SL++L L     D +  F    L F +  C+E  G   Q +LIP G N+ V   N  E+V  Y  +   +   +   A   G L V     LE     + R ++ G    N + L     +  +       +    + FW       + +++  + F T S  +P        M S+ I        +++LP A+TC + L +P"
assert record.rounds[0].alignments[12].hsps[0].sbjct=="NWFHLIGITCGLAIYNSTVVDLHFPLALYKKLLNVKPGLEDLKELSPTEGRSLQEL-LDYPGEDVEETFC---LNFTI--CRESYGVIEQKKLIPGGDNVTVCKDNRQEFVDAYVNYVFQISVHEWYTAFSSGFLKVCGGKVLELFQPSELRAMMVGNSNYNWEELEETAIYKGDYSATHPTV----KLFWETFHEFPLEKKKKFLLFLTGSDRIP-----IYGMASLQIVIQSTASGEEYLPVAHTCYNLLDLP"
assert record.rounds[0].alignments[12].hsps[0].query_start==617
assert record.rounds[0].alignments[12].hsps[0].query_end==865
assert record.rounds[0].alignments[12].hsps[0].sbjct_start==786
assert record.rounds[0].alignments[12].hsps[0].sbjct_end==1025
assert record.rounds[0].alignments[13].hsps[0].query=="TPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNW--HDFAFFDPVMYES---LRQLILASQSSDADAVFSAMDLAFAVDL-----CKEEGGG---------QVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPM-PSITI--------RPPDDQHLPTANTCISRLYVP"
assert record.rounds[0].alignments[13].hsps[0].match=="T +P    + ++  FR +G+++   ++   L  + L     K +L ++ +   HD    DPV+  S   L  ++   +  + D   +   L +A++      C  E  G          +EL   G +IPVT  N+ EY+R      +     +   + R G   V P + L+    E+   L+ G                 + G   +   +  ++ + I+      +++  + F T SP LP    GF+ + P +TI          PDD  LP+  TC++ L +P"
assert record.rounds[0].alignments[13].hsps[0].sbjct=="TAKPAHIAKVKMK-FRFLGKLMAKAIMDFRLVDLPLGLPFYKWMLRQETSLTSHDLFDIDPVVARSVYHLEDIVRQKKRLEQDKSQTKESLQYALETLTMNGCSVEDLGLDFTLPGFPNIELKKGGKDIPVTIHNLEEYLRLVIFWALNEGVSRQFDSFRDGFESVFPLSHLQYFYPEELDQLLCGSKADTWDAKTLMECCRPDHGYTHDS--RAVKFLFEILSSFDNEQQRLFLQFVTGSPRLPVG--GFRSLNPPLTIVRKTFESTENPDD-FLPSVMTCVNYLKLP"
assert record.rounds[0].alignments[13].hsps[0].query_start==605
assert record.rounds[0].alignments[13].hsps[0].query_end==865
assert record.rounds[0].alignments[13].hsps[0].sbjct_start==1684
assert record.rounds[0].alignments[13].hsps[0].sbjct_end==1966
assert record.rounds[0].alignments[14].hsps[0].query=="PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[0].alignments[14].hsps[0].match=="P    +Q LGERL+P +QAM P+ A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[0].alignments[14].hsps[0].sbjct=="PPQEQKQMLGERLFPLIQAMHPSLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[0].alignments[14].hsps[0].query_start==480
assert record.rounds[0].alignments[14].hsps[0].query_end==549
assert record.rounds[0].alignments[14].hsps[0].sbjct_start==554
assert record.rounds[0].alignments[14].hsps[0].sbjct_end==623
assert record.rounds[0].alignments[15].hsps[0].query=="PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[0].alignments[15].hsps[0].match=="P    +Q LGERL+P +QAM P  A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[0].alignments[15].hsps[0].sbjct=="PPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[0].alignments[15].hsps[0].query_start==480
assert record.rounds[0].alignments[15].hsps[0].query_end==549
assert record.rounds[0].alignments[15].hsps[0].sbjct_start==554
assert record.rounds[0].alignments[15].hsps[0].sbjct_end==623
assert record.rounds[0].alignments[16].hsps[0].query=="PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[0].alignments[16].hsps[0].match=="P    +Q LGERL+P +QAM P  A KITGM                   R +V+EA+ ++ AH  +  A"
assert record.rounds[0].alignments[16].hsps[0].sbjct=="PPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRLKVDEAVAVLQAHQAKEAA"
assert record.rounds[0].alignments[16].hsps[0].query_start==480
assert record.rounds[0].alignments[16].hsps[0].query_end==549
assert record.rounds[0].alignments[16].hsps[0].sbjct_start==552
assert record.rounds[0].alignments[16].hsps[0].sbjct_end==621
assert record.rounds[0].alignments[17].hsps[0].query=="RQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAH"
assert record.rounds[0].alignments[17].hsps[0].match=="+Q LGERLYP ++ M    A KITGM                   +A+VEEA+ ++  H"
assert record.rounds[0].alignments[17].hsps[0].sbjct=="KQILGERLYPMIEHMHANLAGKITGMLLEIENSELLHMIEDQEALKAKVEEAVAVLQVH"
assert record.rounds[0].alignments[17].hsps[0].query_start==485
assert record.rounds[0].alignments[17].hsps[0].query_end==543
assert record.rounds[0].alignments[17].hsps[0].sbjct_start==567
assert record.rounds[0].alignments[17].hsps[0].sbjct_end==625
assert record.rounds[0].alignments[18].hsps[0].query=="HRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[0].alignments[18].hsps[0].match=="H + LG+ LYP V+  +PA  +K+TGM                   +A+V EA++++"
assert record.rounds[0].alignments[18].hsps[0].sbjct=="HPRMLGDHLYPLVEQQEPANPAKVTGMLLEMDQAEILHLLESPEALKAKVSEALDVL"
assert record.rounds[0].alignments[18].hsps[0].query_start==484
assert record.rounds[0].alignments[18].hsps[0].query_end==540
assert record.rounds[0].alignments[18].hsps[0].sbjct_start==590
assert record.rounds[0].alignments[18].hsps[0].sbjct_end==646
assert record.rounds[0].alignments[19].hsps[0].query=="RQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[0].alignments[19].hsps[0].match=="R  LGE LYP V+ ++   A+K+TGM                   +A+V EAM+++"
assert record.rounds[0].alignments[19].hsps[0].sbjct=="RTMLGEVLYPLVEQVEAESAAKVTGMLLEMDQTEVLHLLESPEALKAKVAEAMDVL"
assert record.rounds[0].alignments[19].hsps[0].query_start==485
assert record.rounds[0].alignments[19].hsps[0].query_end==540
assert record.rounds[0].alignments[19].hsps[0].sbjct_start==556
assert record.rounds[0].alignments[19].hsps[0].sbjct_end==611
assert record.rounds[0].alignments[20].hsps[0].query=="LADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSK"
assert record.rounds[0].alignments[20].hsps[0].match=="L  Q  +    AR   +  R  + L   +    K  +  T+D    E     + YA   K V+  + ++KKA ED+S L +E + S+"
assert record.rounds[0].alignments[20].hsps[0].sbjct=="LRSQDSMSSDTARTSPVTARTLETLIRLATAHAKARMSKTVDLQDAEEAVELVQYAY-FKKVLEKEKKRKKASEDESDLEDEEEKSQ"
assert record.rounds[0].alignments[20].hsps[0].query_start==141
assert record.rounds[0].alignments[20].hsps[0].query_end==227
assert record.rounds[0].alignments[20].hsps[0].sbjct_start==597
assert record.rounds[0].alignments[20].hsps[0].sbjct_end==682
assert record.rounds[0].alignments[21].hsps[0].query=="PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[0].alignments[21].hsps[0].match=="P  + +Q LGE LYP+V   +   + KITGM                     RV EA+ ++"
assert record.rounds[0].alignments[21].hsps[0].sbjct=="PEESRKQVLGELLYPKVFVREEKLSGKITGMLLEMPNSELLELLEDDSALNERVNEAIGVL"
assert record.rounds[0].alignments[21].hsps[0].query_start==480
assert record.rounds[0].alignments[21].hsps[0].query_end==540
assert record.rounds[0].alignments[21].hsps[0].sbjct_start==581
assert record.rounds[0].alignments[21].hsps[0].sbjct_end==641
assert record.rounds[0].alignments[22].hsps[0].query=="RWRLSLELFGRVFMEDVGAEP----GSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATT------PMAVHR"
assert record.rounds[0].alignments[22].hsps[0].match=="RW       G + ++D   E     GS L   G  EVK    R E+ ++RN   RDL+  V R         MR +  H  +  +        P+AVHR"
assert record.rounds[0].alignments[22].hsps[0].sbjct=="RWIFGAPPAGEITVDDGAVEAMMARGSSLLPKGIREVKGDFSRGEVIRIRNLTGRDLAHGVSRYN----SDAMRMIAGHHSQEISEILGYEYGPVAVHR"
assert record.rounds[0].alignments[22].hsps[0].query_start==279
assert record.rounds[0].alignments[22].hsps[0].query_end==367
assert record.rounds[0].alignments[22].hsps[0].sbjct_start==267
assert record.rounds[0].alignments[22].hsps[0].sbjct_end==361
assert record.rounds[0].alignments[23].hsps[0].query=="FEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTA--IAQAFLSNEKLPNLDCI"
assert record.rounds[0].alignments[23].hsps[0].match=="+E++++K++   +KLR +++R+L  E+ R R+   +Q  +Q         A+    +         +P  GSG  R    A  IAQ     +K  NL  I"
assert record.rounds[0].alignments[23].hsps[0].sbjct=="YEIQQTKYQEAQQKLREKKARELE-EIQRLREKNERQNRQQETGQNNADTASGGSNIAPPVPVPNKKPPSGSGGGRDAKQAALIAQKKREEKKRKNLQII"
assert record.rounds[0].alignments[23].hsps[0].query_start==309
assert record.rounds[0].alignments[23].hsps[0].query_end==406
assert record.rounds[0].alignments[23].hsps[0].sbjct_start==839
assert record.rounds[0].alignments[23].hsps[0].sbjct_end==937
assert record.rounds[0].alignments[24].hsps[0].query=="LVEVTMDRNCLEVLPTKMSYAANLK-NVMNMQNRQKKAGEDQSMLAEEA---------DSSKPGPSAHDVAAQLKSSLLAEIGLTESEGP"
assert record.rounds[0].alignments[24].hsps[0].match=="L E+  +   L+V+ T+ S    LK  + N+Q  QK+  ++++   E+          + S  GP  H+V+     SL +E+  +++ GP"
assert record.rounds[0].alignments[24].hsps[0].sbjct=="LRELEEENTSLKVIYTRSSEIEELKATIENLQENQKRLQKEKAEEIEQLHEVIEKLQHELSLMGPVVHEVSDSQAGSLQSELLCSQAGGP"
assert record.rounds[0].alignments[24].hsps[0].query_start==176
assert record.rounds[0].alignments[24].hsps[0].query_end==255
assert record.rounds[0].alignments[24].hsps[0].sbjct_start==1683
assert record.rounds[0].alignments[24].hsps[0].sbjct_end==1772
assert record.rounds[0].alignments[25].hsps[0].query=="AFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVA---EQPLHAMRKG-LLDVLPKNSLEDL"
assert record.rounds[0].alignments[25].hsps[0].match=="AF VDLC+ +GGG+       + +P++ Q++ +Y+    E    VV    E+ L A+R    +D++   +L  L"
assert record.rounds[0].alignments[25].hsps[0].sbjct=="AFLVDLCERQGGGR------QLRLPMSRQDIADYLGLTIETVSRVVTKLKERSLIALRDARTIDIMKPEALRSL"
assert record.rounds[0].alignments[25].hsps[0].query_start==691
assert record.rounds[0].alignments[25].hsps[0].query_end==760
assert record.rounds[0].alignments[25].hsps[0].sbjct_start==142
assert record.rounds[0].alignments[25].hsps[0].sbjct_end==209
assert record.rounds[0].alignments[26].hsps[0].query=="KNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAF--FDPVMYESLRQLILASQSSDADAVF"
assert record.rounds[0].alignments[26].hsps[0].match=="K  EA L+  +NI R +      NE    T N  V+K L GR VNW  +    F  ++Y     +   + SS+ +  F"
assert record.rounds[0].alignments[26].hsps[0].sbjct=="KELEAALDISKNIARSI------NENQRRTENHQVVKKLYGRVVNWKGYRISKFGELLYFDKVFISTTNSSSEPEREF"
assert record.rounds[0].alignments[26].hsps[0].query_start==610
assert record.rounds[0].alignments[26].hsps[0].query_end==685
assert record.rounds[0].alignments[26].hsps[0].sbjct_start==435
assert record.rounds[0].alignments[26].hsps[0].sbjct_end==506
assert record.rounds[1].alignments[0].hsps[0].query=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[1].alignments[0].hsps[0].match=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQ            GI            TSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQ                                      QLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGM                   RARVEEAMELIVAHGRENGA            EKVQENRKRHGSSRSVV            NAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP                   NFGFV"
assert record.rounds[1].alignments[0].hsps[0].sbjct=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQLERKRTRELLELGIDNEDSEHENDDDTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQRLRNRGERDREREREREMRRSSGLRAGSRRDRDRDFRRQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMLLELSPAQLLLLLASEDSLRARVEEAMELIVAHGRENGADSILDLGLLDSSEKVQENRKRHGSSRSVVDMDLDDTDDGDDNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPLYSSKQILKQKLLLAIKTKNFGFV"
assert record.rounds[1].alignments[0].hsps[0].query_start==1
assert record.rounds[1].alignments[0].hsps[0].query_end==889
assert record.rounds[1].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[1].alignments[0].hsps[0].sbjct_end==889
assert record.rounds[1].alignments[1].hsps[0].query=="QCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSL-EVDRDRDLLIQQTMRQLNNHFGR--RCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGN---PSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRM-LVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKL-LQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[1].alignments[1].hsps[0].match=="Q +F G+   + +L  +  L+L+LFGR FM+DVG E GS+L EL GF VKE +FRR MEKLRN Q RDL L +++R+R+ LI QT ++LN  FG   R    P+  +RVKVTFKDEPGEGSGVARSFYT+IA+A L++ K+PNL+ +Q      H+  +                                        L++D RP+ P +  +   P    D L  H Q +GERLYP++ ++    A KITGM                   R +V EA+E+I    +   +               Q ++ +   S  VV            N PLFY PGKRGFY PR G  +  R+N FRNIGR++GLCLLQNEL P+ L RHV+K +L RK+ +HD AFFDP +YES RQ+I  +Q+ + +   + M+L F +DL KEEG G  ELIP G ++ V+ + +Y       ++   +   E+ L A++ G+ DVLP NS+ +LTAED RLL+NG G++NV  LIS+T+FNDES E  +KL  +FK+WFWSIVE+M++ ERQ LVYFWT SP+LPASEEGFQP+PS+TIRP DD HLPTANTCISRLY+P                   NFGFV"
assert record.rounds[1].alignments[1].hsps[0].sbjct=="QSNFSGIAPCNFLLSAK--LTLDLFGRPFMDDVGMEHGSVLPELRGFPVKEMRFRRHMEKLRNGQQRDLVLCKLERNRESLIVQTFKELNTQFGNQSRRIQPPITFNRVKVTFKDEPGEGSGVARSFYTSIAEALLASAKIPNLESVQVGTN--HSKYVVPFSSILRTTVSGSSRDQSTLQRRGSNSKILWRSARERKALNLDARPYTPPNSSDNATPESLNDHLSVHLQQIGERLYPKIHSINQTHAPKITGMLLEIPTPQLLSVISSDETLRQKVNEAIEIITFKQKSETSA--------------QSSQPKKSPSVVVVDPVDDD------NEPLFYSPGKRGFYIPRQGFASFERINAFRNIGRLIGLCLLQNELLPLFLQRHVLKYILRRKIKFHDLAFFDPALYESFRQIIQNAQTKEGEETINRMELCFVIDLMKEEGCGNRELIPGGRDV-VSHRVIYSSTSDAIQNIDXIKSQEKALEALKDGVFDVLPDNSMINLTAEDLRLLLNGVGDINVSTLISYTTFNDESSEGPDKLX-KFKKWFWSIVEKMNIMERQHLVYFWTGSPALPASEEGFQPLPSVTIRPADDSHLPTANTCISRLYIPLYSSKSILRSKMLMAIKSKNFGFV"
assert record.rounds[1].alignments[1].hsps[0].query_start==263
assert record.rounds[1].alignments[1].hsps[0].query_end==889
assert record.rounds[1].alignments[1].hsps[0].sbjct_start==2287
assert record.rounds[1].alignments[1].hsps[0].sbjct_end==2895
assert record.rounds[1].alignments[1].hsps[1].query=="SARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDE------SLPA--ETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMD---RNCLEVLPTKMSYAANLK"
assert record.rounds[1].alignments[1].hsps[1].match=="++R DF  Y LSLMRSH  EH D LPVLD+ +L+H+AYV  A +Y+++  N     D                I              + A L + + +      S+P+  +  + H FF RS+S   LGC  P  FE+PL  A+PLAD+PHLLQPN+++++LF      + +++  SG      T D    +  +  PT++ ++ +LK"
assert record.rounds[1].alignments[1].hsps[1].sbjct=="NSRRDFFTYCLSLMRSHTSEHRDALPVLDITALRHIAYVLDAFVYYMR--NDSGFYDKQDSISGR--------INNLSPMTESYDTDDELANLEEFNADVQMSASSMPSGSQGTRRHAFFARSESTLTLGCSAPEGFELPLDMAMPLADKPHLLQPNSKRQELFANLPLLVTTNANNSG-----ATNDGDGASIFDYTPTRLGFSNSLK"
assert record.rounds[1].alignments[1].hsps[1].query_start==3
assert record.rounds[1].alignments[1].hsps[1].query_end==200
assert record.rounds[1].alignments[1].hsps[1].sbjct_start==1913
assert record.rounds[1].alignments[1].hsps[1].sbjct_end==2106
assert record.rounds[1].alignments[2].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[2].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N     Q  +WFW  V  M   +R  L+ F T +  +P    A   G     S T+      + LP A+TC +RL +P"
assert record.rounds[1].alignments[2].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLHKPITLHDMESVDSEYYNSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKNGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSANH----QVIQWFWKAVLMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPEKLPRAHTCFNRLDLP"
assert record.rounds[1].alignments[2].hsps[0].query_start==606
assert record.rounds[1].alignments[2].hsps[0].query_end==865
assert record.rounds[1].alignments[2].hsps[0].sbjct_start==649
assert record.rounds[1].alignments[2].hsps[0].sbjct_end==901
assert record.rounds[1].alignments[3].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIR-PPDDQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[3].hsps[0].match=="P  G N E  LN F+ IGR++GL +               K++L +KV   D    D   Y SL  ++    +   D  FS  D  F   +        ++L PNG NI VT +N  EYV      R+    E+  +A  +G  +++P+  +      +  LL+ G  E++++     T +   S  +     Q  +WFW +++  S  ++  L+ F T +  +P +     +G       TI    +   LP A+TC +RL +P"
assert record.rounds[1].alignments[3].hsps[0].sbjct=="PHSGINPE-HLNYFKFIGRVIGLAIFHRRFVDAFFVVSFYKMILQKKVTLQDMESMDAEYYRSLVWILDNDITGVLDLTFSVEDNCFGEVVT-------IDLKPNGRNIEVTEENKREYVDLVTVWRIQKRIEEQFNAFHEGFSELIPQELINVFDERELELLIGGISEIDMEDWKKHTDYRSYSEND-----QIIKWFWELMDEWSNEKKSRLLQFTTGTSRIPVNGFKDLQGSDGPRKFTIEKAGEPNKLPKAHTCFNRLDLP"
assert record.rounds[1].alignments[3].hsps[0].query_start==606
assert record.rounds[1].alignments[3].hsps[0].query_end==865
assert record.rounds[1].alignments[3].hsps[0].sbjct_start==491
assert record.rounds[1].alignments[3].hsps[0].sbjct_end==742
assert record.rounds[1].alignments[4].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[4].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N     Q   WFW  V  M   +R  L+ F T +  +P    A   G     S T+        LP A+TC +RL +P"
assert record.rounds[1].alignments[4].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLQKLITLHDMESVDSEYYSSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKTGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSMNH----QVIHWFWKAVWMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPDKLPRAHTCFNRLDLP"
assert record.rounds[1].alignments[4].hsps[0].query_start==606
assert record.rounds[1].alignments[4].hsps[0].query_end==865
assert record.rounds[1].alignments[4].hsps[0].sbjct_start==679
assert record.rounds[1].alignments[4].hsps[0].sbjct_end==931
assert record.rounds[1].alignments[5].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRPPDD-QHLPTANTCISRLYVP"
assert record.rounds[1].alignments[5].hsps[0].match=="P  G N E  LN F+ IGR++GL +             + K++L +KV   D    D  +Y SL  ++  S     D  FSA D  F   +        V+L P+G NI VT  N  EYV  Y + R++   ++   A   G  +++P++ +      +  LL+ G  E++++     T +      +     +  +WFW  V      +R  L+ F T +  +P +     +G       TI    + Q LP ++TC +R+ +P"
assert record.rounds[1].alignments[5].hsps[0].sbjct=="PNSGINPE-HLNYFKFIGRVVGLGVFHRRFLDAFFVGALYKMMLRKKVVLQDMEGVDAEVYNSLNWMLENSIDGVLDLTFSADDERFGEVVT-------VDLKPDGRNIEVTDGNKKEYVELYTQWRIVDRVQEQFKAFMDGFNELIPEDLVTVFDERELELLIGGIAEIDIEDWKKHTDYRGYQESD-----EVIQWFWKCVSEWDNEQRARLLQFTTGTSRIPVNGFKDLQGSDGPRRFTIEKAGEVQQLPKSHTCFNRVDLP"
assert record.rounds[1].alignments[5].hsps[0].query_start==606
assert record.rounds[1].alignments[5].hsps[0].query_end==865
assert record.rounds[1].alignments[5].hsps[0].sbjct_start==533
assert record.rounds[1].alignments[5].hsps[0].sbjct_end==784
assert record.rounds[1].alignments[6].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLED-LTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[6].hsps[0].match=="F  IG +LGL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    + + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+     E+  LL+ G   ++ Q L   T +  + G   + +L   R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[1].alignments[6].hsps[0].sbjct=="FTLIGIVLGLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGNVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVNLYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRDSVL--IREFWEIVHSFTDEQKRLFLQFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[1].alignments[6].hsps[0].query_start==619
assert record.rounds[1].alignments[6].hsps[0].query_end==865
assert record.rounds[1].alignments[6].hsps[0].sbjct_start==612
assert record.rounds[1].alignments[6].hsps[0].sbjct_end==850
assert record.rounds[1].alignments[7].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLED-LTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[7].hsps[0].match=="F  IG + GL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    S + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+     E+  LL+ G   ++ Q L   T +  + G   E ++   R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[1].alignments[7].hsps[0].sbjct=="FTLIGIL-GLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGSVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVISYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRESVV--IREFWEIVHSFTDEQKRLFLLFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[1].alignments[7].hsps[0].query_start==619
assert record.rounds[1].alignments[7].hsps[0].query_end==865
assert record.rounds[1].alignments[7].hsps[0].sbjct_start==623
assert record.rounds[1].alignments[7].hsps[0].sbjct_end==860
assert record.rounds[1].alignments[8].hsps[0].query=="YTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAV------DLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPK-NSLEDLTAEDFRLLVNGCGE---VNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYV"
assert record.rounds[1].alignments[8].hsps[0].match=="+    GKN++  L  +   G ++GL +  + +  +   + + K L    +++ D++   P    +L +++  ++ +  D      +  +        D    +    VEL  NG N+P+T  N +E+V K+ E  +    E   +    G   V  + NS++   +E+   LV G  E    + + L S T +     +++  +     WFW I+E      ++ L+ F T+S  +PA+     P     +   D   LP A+TC + + +"
assert record.rounds[1].alignments[8].hsps[0].sbjct=="FDKSKGKNSQLEL--YYLFGVVMGLAIFNSTILDLQFPKALYKKLCSEPLSFEDYSELFPETSRNLIKMLNYTEDNFEDVFSLTFETTYRNNNWILNDSKSSKEYVTVELCENGRNVPITQSNKHEFVMKWVEFYLEKSIEPQYNKFVSGFKRVFAECNSIKLFNSEELERLVCGDEEQTKFDFKSLRSVTKYVGGFSDDSRAVC----WFWEIIESWDYPLQKKLLQFVTASDRIPATGISTIPFKISLLGSHDSDDLPLAHTCFNEICL"
assert record.rounds[1].alignments[8].hsps[0].query_start==604
assert record.rounds[1].alignments[8].hsps[0].query_end==864
assert record.rounds[1].alignments[8].hsps[0].sbjct_start==602
assert record.rounds[1].alignments[8].hsps[0].sbjct_end==866
assert record.rounds[1].alignments[9].hsps[0].query=="GRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRPPD--DQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[9].hsps[0].match=="GR +   ++   L      R   K +LG+ V + D    D   Y+ L  L+        D      DL F+ ++ +E G  +V +L PNG NI VT +N  EYV    + RM     + L A  +G  +++PK  +   T ++  LL  G   +++  L S T ++     + +      +WFW  +      +R   + F T +  +P    A+ EG   +    I   D     LP+A+TC ++L +P"
assert record.rounds[1].alignments[9].hsps[0].sbjct=="GRYVAKAVMTTALLECYFTRSFYKHILGKSVRYTDMESEDYHFYQGLVYLLEN------DVSTLGYDLTFSTEV-QEFGVCEVRDLKPNGANILVTEENKKEYVHLVCQMRMTGAIRKQLAAFLEGFYEIIPKRLISIFTEQELELLYTGLPTIDIDDLKSNTEYHKYQSNSIQ-----IQWFWRALRSFDQADRAKFLQFVTGTSKVPLQGFAALEGMNGIQKFQIHRDDRSTDRLPSAHTCFNQLDLP"
assert record.rounds[1].alignments[9].hsps[0].query_start==623
assert record.rounds[1].alignments[9].hsps[0].query_end==865
assert record.rounds[1].alignments[9].hsps[0].sbjct_start==45
assert record.rounds[1].alignments[9].hsps[0].sbjct_end==282
assert record.rounds[1].alignments[10].hsps[0].query=="TPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNW--HDFAFFDPVMYES---LRQLILASQSSDADAVFSAMDLAFAVDL-----CKEE---------GGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRP------PDDQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[10].hsps[0].match=="T +P    + ++  FR +G+++   ++   L  + L     K +L ++ +   HD    DPV+  S   L  ++   +  + D   +   L +A++      C  E         G   +EL   G +IPVT  N+ EY+R      +     +   + R G   V P + L+    E+   L+ G                 + G   +   +  ++ + I+      +++  + F T SP LP         P   +R         D  LP+  TC++ L +P"
assert record.rounds[1].alignments[10].hsps[0].sbjct=="TAKPAHIAKVKMK-FRFLGKLMAKAIMDFRLVDLPLGLPFYKWMLRQETSLTSHDLFDIDPVVARSVYHLEDIVRQKKRLEQDKSQTKESLQYALETLTMNGCSVEDLGLDFTLPGFPNIELKKGGKDIPVTIHNLEEYLRLVIFWALNEGVSRQFDSFRDGFESVFPLSHLQYFYPEELDQLLCGSKADTWDAKTLMECCRPDHGYTHDS--RAVKFLFEILSSFDNEQQRLFLQFVTGSPRLPVGGFRSLNPPLTIVRKTFESTENPDDFLPSVMTCVNYLKLP"
assert record.rounds[1].alignments[10].hsps[0].query_start==605
assert record.rounds[1].alignments[10].hsps[0].query_end==865
assert record.rounds[1].alignments[10].hsps[0].sbjct_start==1684
assert record.rounds[1].alignments[10].hsps[0].sbjct_end==1966
assert record.rounds[1].alignments[11].hsps[0].query=="EARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGG--QVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPS-ITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[11].hsps[0].match=="    N F  IG   GL +  + +  +     + K LL  K    D     P    SL++L+      D +  F    L F +  C+E  G   Q +LIP G N+ V   N  E+V  Y  +   +   +   A   G L V     LE     + R ++ G    N + L     +  +       +    + FW       + +++  + F T S  +P    G   +   I      +++LP A+TC + L +P"
assert record.rounds[1].alignments[11].hsps[0].sbjct=="FVEHNWFHLIGITCGLAIYNSTVVDLHFPLALYKKLLNVKPGLEDLKELSPTEGRSLQELLDY-PGEDVEETFC---LNFTI--CRESYGVIEQKKLIPGGDNVTVCKDNRQEFVDAYVNYVFQISVHEWYTAFSSGFLKVCGGKVLELFQPSELRAMMVGNSNYNWEELEETAIYKGDYSATHPTV----KLFWETFHEFPLEKKKKFLLFLTGSDRIPI--YGMASLQIVIQSTASGEEYLPVAHTCYNLLDLP"
assert record.rounds[1].alignments[11].hsps[0].query_start==613
assert record.rounds[1].alignments[11].hsps[0].query_end==865
assert record.rounds[1].alignments[11].hsps[0].sbjct_start==782
assert record.rounds[1].alignments[11].hsps[0].sbjct_end==1025
assert record.rounds[1].alignments[12].hsps[0].query=="RLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLL----GRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGE-VNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPM-PSITIRPPDDQ--HLPTANTCISRLYVP"
assert record.rounds[1].alignments[12].hsps[0].match=="+L     +G+++G CL ++ L  ++     +K LL    G   ++ D   +D V+Y +L +L+    ++D      ++DL F +D   E     V+LIPNG    VT  NV  YV K  ++++     +P+ A   GL  ++  + +E   + + ++L++G  + +++  L S T +     E+     Q    FW ++      E+ + + F TS P  P   +GF+ + P   IR    +   LPTA+TC++ L +P"
assert record.rounds[1].alignments[12].hsps[0].sbjct=="KLKYIWFLGKVVGKCLYEHVLIDVSFADFFLKKLLNYSNGFLSSFSDLGSYDSVLYNNLIKLLN--MTTDE---IKSLDLTFEIDE-PESSAKVVDLIPNGSKTYVTKDNVLLYVTKVTDYKLNKRCFKPVSAFHGGLSVIIAPHWMEMFNSIELQMLISGERDNIDLDDLKSNTEYGGYKEED-----QTIVDFWEVLNEFKFEEKLNFLKFVTSVPQAP--LQGFKALDPKFGIRNAGTEKYRLPTASTCVNLLKLP"
assert record.rounds[1].alignments[12].hsps[0].query_start==615
assert record.rounds[1].alignments[12].hsps[0].query_end==865
assert record.rounds[1].alignments[12].hsps[0].sbjct_start==640
assert record.rounds[1].alignments[12].hsps[0].sbjct_end==885
assert record.rounds[1].alignments[13].hsps[0].query=="RPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLG------------RKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP-ASEEGFQPMPSITIRPPD-----DQHLPTANTCISRLYVP"
assert record.rounds[1].alignments[13].hsps[0].match==" P  N E  +  F  +G  +   LL N +     ++   ++L               +         DP++ +SL+ ++      D +    ++ L F V      G   +ELIP G N  +   NV EY+    +  +    E+ L A  +G   V     +  L  ++  + + G  E +  M   +T+ N E G   +  +     F SI+      ER+  + F T SP LP    +   P  ++ ++  +     D++LP+  TC + L +P"
assert record.rounds[1].alignments[13].hsps[0].sbjct=="NPFSNNEKVIELFGYLGTFVARSLLDNRILDFRFSKVFFELLHRMSTPNVTTVPSDVETCLLMIELVDPLLAKSLKYIVAN---KDDNMTLESLSLTFTVP-----GNDDIELIPGGCNKSLNSSNVEEYIHGVIDQILGKGIEKQLKAFIEGFSKVFSYERMLILFPDEL-VDIFGRVEEDWSMATLYTNLNAEHGYTMDSSI--IHDFISIISAFGKHERRLFLQFLTGSPKLPIGGFKSLNPKFTVVLKHAEDGLTADEYLPSVMTCANYLKLP"
assert record.rounds[1].alignments[13].hsps[0].query_start==607
assert record.rounds[1].alignments[13].hsps[0].query_end==865
assert record.rounds[1].alignments[13].hsps[0].sbjct_start==1192
assert record.rounds[1].alignments[13].hsps[0].sbjct_end==1457
assert record.rounds[1].alignments[14].hsps[0].query=="PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[1].alignments[14].hsps[0].match=="P    +Q LGERL+P +QAM P  A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[1].alignments[14].hsps[0].sbjct=="PPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[1].alignments[14].hsps[0].query_start==480
assert record.rounds[1].alignments[14].hsps[0].query_end==549
assert record.rounds[1].alignments[14].hsps[0].sbjct_start==554
assert record.rounds[1].alignments[14].hsps[0].sbjct_end==623
assert record.rounds[1].alignments[15].hsps[0].query=="PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[1].alignments[15].hsps[0].match=="P    +Q LGERL+P +QAM P  A KITGM                   R +V+EA+ ++ AH  +  A"
assert record.rounds[1].alignments[15].hsps[0].sbjct=="PPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRLKVDEAVAVLQAHQAKEAA"
assert record.rounds[1].alignments[15].hsps[0].query_start==480
assert record.rounds[1].alignments[15].hsps[0].query_end==549
assert record.rounds[1].alignments[15].hsps[0].sbjct_start==552
assert record.rounds[1].alignments[15].hsps[0].sbjct_end==621
assert record.rounds[1].alignments[16].hsps[0].query=="PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[1].alignments[16].hsps[0].match=="P    +Q LGERL+P +QAM P+ A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[1].alignments[16].hsps[0].sbjct=="PPQEQKQMLGERLFPLIQAMHPSLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[1].alignments[16].hsps[0].query_start==480
assert record.rounds[1].alignments[16].hsps[0].query_end==549
assert record.rounds[1].alignments[16].hsps[0].sbjct_start==554
assert record.rounds[1].alignments[16].hsps[0].sbjct_end==623
assert record.rounds[1].alignments[17].hsps[0].query=="PDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[1].alignments[17].hsps[0].match=="       +Q LGERLYP ++ M    A KITGM                   +A+VEEA+ ++  H     A"
assert record.rounds[1].alignments[17].hsps[0].sbjct=="NAKPQEQKQILGERLYPMIEHMHANLAGKITGMLLEIENSELLHMIEDQEALKAKVEEAVAVLQVHRVTEPA"
assert record.rounds[1].alignments[17].hsps[0].query_start==478
assert record.rounds[1].alignments[17].hsps[0].query_end==549
assert record.rounds[1].alignments[17].hsps[0].sbjct_start==560
assert record.rounds[1].alignments[17].hsps[0].sbjct_end==631
assert record.rounds[1].alignments[18].hsps[0].query=="PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIV"
assert record.rounds[1].alignments[18].hsps[0].match=="P  + +Q LGE LYP+V   +   + KITGM                     RV EA+ ++ "
assert record.rounds[1].alignments[18].hsps[0].sbjct=="PEESRKQVLGELLYPKVFVREEKLSGKITGMLLEMPNSELLELLEDDSALNERVNEAIGVLQ"
assert record.rounds[1].alignments[18].hsps[0].query_start==480
assert record.rounds[1].alignments[18].hsps[0].query_end==541
assert record.rounds[1].alignments[18].hsps[0].sbjct_start==581
assert record.rounds[1].alignments[18].hsps[0].sbjct_end==642
assert record.rounds[1].alignments[19].hsps[0].query=="PDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[1].alignments[19].hsps[0].match=="       R  LGE LYP V+ ++   A+K+TGM                   +A+V EAM+++"
assert record.rounds[1].alignments[19].hsps[0].sbjct=="NATPEQQRTMLGEVLYPLVEQVEAESAAKVTGMLLEMDQTEVLHLLESPEALKAKVAEAMDVL"
assert record.rounds[1].alignments[19].hsps[0].query_start==478
assert record.rounds[1].alignments[19].hsps[0].query_end==540
assert record.rounds[1].alignments[19].hsps[0].sbjct_start==549
assert record.rounds[1].alignments[19].hsps[0].sbjct_end==611
assert record.rounds[1].alignments[20].hsps[0].query=="LPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[1].alignments[20].hsps[0].match=="   H + LG+ LYP V+  +PA  +K+TGM                   +A+V EA++++"
assert record.rounds[1].alignments[20].hsps[0].sbjct=="PDKHPRMLGDHLYPLVEQQEPANPAKVTGMLLEMDQAEILHLLESPEALKAKVSEALDVL"
assert record.rounds[1].alignments[20].hsps[0].query_start==481
assert record.rounds[1].alignments[20].hsps[0].query_end==540
assert record.rounds[1].alignments[20].hsps[0].sbjct_start==587
assert record.rounds[1].alignments[20].hsps[0].sbjct_end==646
assert record.rounds[1].alignments[21].hsps[0].query=="PASEGNPSDDPDP----LPAHRQALGERLYPRVQA--MQPAFASKITGM"
assert record.rounds[1].alignments[21].hsps[0].match=="P   G P +  D         RQALGE+LY +V A       A KITGM"
assert record.rounds[1].alignments[21].hsps[0].sbjct=="PPQGGFPRNANDNNQFYQQKQRQALGEQLYKKVSAKTSNEEAAGKITGM"
assert record.rounds[1].alignments[21].hsps[0].query_start==468
assert record.rounds[1].alignments[21].hsps[0].query_end==510
assert record.rounds[1].alignments[21].hsps[0].sbjct_start==484
assert record.rounds[1].alignments[21].hsps[0].sbjct_end==532
assert record.rounds[1].alignments[22].hsps[0].query=="SAMDLAFAVDLCKEEGGGQ---VELIPNGVNIPVTPQNVYEYVRKY-AEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCG"
assert record.rounds[1].alignments[22].hsps[0].match=="    L F +D  K         +EL P G            YV      + +    E    A ++G     P  +L  L  E     +    "
assert record.rounds[1].alignments[22].hsps[0].sbjct=="LREQLCFTIDPEKAGDFDDAVAIELTPEGY--------YKLYVHIADVSYYVREGTETDKEAYKRGFTYYFPDRALHML-PEKLSAKLCSLR"
assert record.rounds[1].alignments[22].hsps[0].query_start==686
assert record.rounds[1].alignments[22].hsps[0].query_end==773
assert record.rounds[1].alignments[22].hsps[0].sbjct_start==243
assert record.rounds[1].alignments[22].hsps[0].sbjct_end==325
assert record.rounds[1].alignments[23].hsps[0].query=="LAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKE"
assert record.rounds[1].alignments[23].hsps[0].match=="+A +  S+  GP +      ++ S  A++   E+  P  T  RP  +  G V S      +W    E   R F         S  + L  F+  E"
assert record.rounds[1].alignments[23].hsps[0].sbjct=="VASQYSSTPSGPVSVSAMRAVQRSFSAQLAHNEARQPEPTYLRPVSTSYGPVPSTQSAQEQWYSPTEAQFRAFTAAHSMPTTSAQSPLTTFDTPE"
assert record.rounds[1].alignments[23].hsps[0].query_start==219
assert record.rounds[1].alignments[23].hsps[0].query_end==313
assert record.rounds[1].alignments[23].hsps[0].sbjct_start==702
assert record.rounds[1].alignments[23].hsps[0].sbjct_end==796
assert record.rounds[1].alignments[24].hsps[0].query=="KRGFYTPRP-GKNTEARLNCFRN----IGRILGLCLLQNELC----PITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVN-IPVTPQNVY---EYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKN------------SLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMS"
assert record.rounds[1].alignments[24].hsps[0].match=="KR F+  +   KN +  +  F      I  +  L LL  E+        + +   K       ++ D       + E        S + + D++           +  +        I  G     +   +     E++   +         +     RK     +                +      + +L  NG   ++   +                +      F   V++  "
assert record.rounds[1].alignments[24].hsps[0].sbjct=="KRNFHILKSLSKNEKETMRHFDFEMEWIIDLFSLSLLHGEMIFYEYDCNITK-FYK-------SFQDLWDMVIHISEKCYNYFFNSDALEVDSLTKFYTNRIVEIVANKVLVIVPAFILRGDRFKTIQYADKKKMIEFLYGVSSVYFNEFGFEYYRCFRKMFTAKIAYKILNRSCEKDAWRIILKFLLNELKLEDNGDSYIDYNDMRLN------------DICPIILEFQETVQKYD"
assert record.rounds[1].alignments[24].hsps[0].query_start==600
assert record.rounds[1].alignments[24].hsps[0].query_end==812
assert record.rounds[1].alignments[24].hsps[0].sbjct_start==523
assert record.rounds[1].alignments[24].hsps[0].sbjct_end==740
assert record.rounds[1].alignments[25].hsps[0].query=="QLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPN-GVNIPVTPQNVYEYVRKYAE----------HRMLVVAEQPLHAMRK-------GLLDVLPKNSLEDLTAEDFRLLV"
assert record.rounds[1].alignments[25].hsps[0].match=="+LI+    +        +D AF       E   ++ +  N GV I +T   V  YV +Y +          ++++      +  +++           ++ +  L+ L  +D R L+"
assert record.rounds[1].alignments[25].hsps[0].sbjct=="ELIVNGIINGDLKTSLQVDAAFKYVKANGEASTKMGMNENSGVGIEITEDQVRNYVMQYIQENKERILTERYKLVPGIFADVKNLKELKWADPRSFKPIIDQEVLKLLGPKDERDLI"
assert record.rounds[1].alignments[25].hsps[0].query_start==671
assert record.rounds[1].alignments[25].hsps[0].query_end==769
assert record.rounds[1].alignments[25].hsps[0].sbjct_start==71
assert record.rounds[1].alignments[25].hsps[0].sbjct_end==187
assert record.rounds[2].alignments[0].hsps[0].query=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[2].alignments[0].hsps[0].match=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQ            GI            TSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQ                                      QLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGM                   RARVEEAMELIVAHGRENGA            EKVQENRKRHGSSRSVV            NAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP                   NFGFV"
assert record.rounds[2].alignments[0].hsps[0].sbjct=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQLERKRTRELLELGIDNEDSEHENDDDTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQRLRNRGERDREREREREMRRSSGLRAGSRRDRDRDFRRQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMLLELSPAQLLLLLASEDSLRARVEEAMELIVAHGRENGADSILDLGLLDSSEKVQENRKRHGSSRSVVDMDLDDTDDGDDNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPLYSSKQILKQKLLLAIKTKNFGFV"
assert record.rounds[2].alignments[0].hsps[0].query_start==1
assert record.rounds[2].alignments[0].hsps[0].query_end==889
assert record.rounds[2].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[2].alignments[0].hsps[0].sbjct_end==889
assert record.rounds[2].alignments[1].hsps[0].query=="QCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSL-EVDRDRDLLIQQTMRQLNNHFGR--RCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGN---PSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRM-LVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLL-QFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[2].alignments[1].hsps[0].match=="Q +F G+   + +L  +  L+L+LFGR FM+DVG E GS+L EL GF VKE +FRR MEKLRN Q RDL L +++R+R+ LI QT ++LN  FG   R    P+  +RVKVTFKDEPGEGSGVARSFYT+IA+A L++ K+PNL+ +Q      H+  +                                        L++D RP+ P +  +   P    D L  H Q +GERLYP++ ++    A KITGM                   R +V EA+E+I    +   +               Q ++ +   S  VV            N PLFY PGKRGFY PR G  +  R+N FRNIGR++GLCLLQNEL P+ L RHV+K +L RK+ +HD AFFDP +YES RQ+I  +Q+ + +   + M+L F +DL KEEG G  ELIP G ++ V+ + +Y       ++   +   E+ L A++ G+ DVLP NS+ +LTAED RLL+NG G++NV  LIS+T+FNDES E  +K L +FK+WFWSIVE+M++ ERQ LVYFWT SP+LPASEEGFQP+PS+TIRP DD HLPTANTCISRLY+P                   NFGFV"
assert record.rounds[2].alignments[1].hsps[0].sbjct=="QSNFSGIAPCNFLLSAK--LTLDLFGRPFMDDVGMEHGSVLPELRGFPVKEMRFRRHMEKLRNGQQRDLVLCKLERNRESLIVQTFKELNTQFGNQSRRIQPPITFNRVKVTFKDEPGEGSGVARSFYTSIAEALLASAKIPNLESVQVGTN--HSKYVVPFSSILRTTVSGSSRDQSTLQRRGSNSKILWRSARERKALNLDARPYTPPNSSDNATPESLNDHLSVHLQQIGERLYPKIHSINQTHAPKITGMLLEIPTPQLLSVISSDETLRQKVNEAIEIITFKQKSETSA--------------QSSQPKKSPSVVVVDPVDDD------NEPLFYSPGKRGFYIPRQGFASFERINAFRNIGRLIGLCLLQNELLPLFLQRHVLKYILRRKIKFHDLAFFDPALYESFRQIIQNAQTKEGEETINRMELCFVIDLMKEEGCGNRELIPGGRDV-VSHRVIYSSTSDAIQNIDXIKSQEKALEALKDGVFDVLPDNSMINLTAEDLRLLLNGVGDINVSTLISYTTFNDESSEGPDK-LXKFKKWFWSIVEKMNIMERQHLVYFWTGSPALPASEEGFQPLPSVTIRPADDSHLPTANTCISRLYIPLYSSKSILRSKMLMAIKSKNFGFV"
assert record.rounds[2].alignments[1].hsps[0].query_start==263
assert record.rounds[2].alignments[1].hsps[0].query_end==889
assert record.rounds[2].alignments[1].hsps[0].sbjct_start==2287
assert record.rounds[2].alignments[1].hsps[0].sbjct_end==2895
assert record.rounds[2].alignments[1].hsps[1].query=="SARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDE------SLPA--ETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMD---RNCLEVLPTKMSYAANLK"
assert record.rounds[2].alignments[1].hsps[1].match=="++R DF  Y LSLMRSH  EH D LPVLD+ +L+H+AYV  A +Y+++  N     D                I              + A L + + +      S+P+  +  + H FF RS+S   LGC  P  FE+PL  A+PLAD+PHLLQPN+++++LF      + +++  SG      T D    +  +  PT++ ++ +LK"
assert record.rounds[2].alignments[1].hsps[1].sbjct=="NSRRDFFTYCLSLMRSHTSEHRDALPVLDITALRHIAYVLDAFVYYMR--NDSGFYDKQDSISGR--------INNLSPMTESYDTDDELANLEEFNADVQMSASSMPSGSQGTRRHAFFARSESTLTLGCSAPEGFELPLDMAMPLADKPHLLQPNSKRQELFANLPLLVTTNANNSG-----ATNDGDGASIFDYTPTRLGFSNSLK"
assert record.rounds[2].alignments[1].hsps[1].query_start==3
assert record.rounds[2].alignments[1].hsps[1].query_end==200
assert record.rounds[2].alignments[1].hsps[1].sbjct_start==1913
assert record.rounds[2].alignments[1].hsps[1].sbjct_end==2106
assert record.rounds[2].alignments[2].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRPP-DDQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[2].hsps[0].match=="P  G N E  LN F+ IGR++GL +               K++L +KV   D    D   Y SL  ++    +   D  FS  D  F   +        ++L PNG NI VT +N  EYV      R+    E+  +A  +G  +++P+  +      +  LL+ G  E++++     T +   S  +     Q  +WFW +++  S  ++  L+ F T +  +P +     +G       TI    +   LP A+TC +RL +P"
assert record.rounds[2].alignments[2].hsps[0].sbjct=="PHSGINPE-HLNYFKFIGRVIGLAIFHRRFVDAFFVVSFYKMILQKKVTLQDMESMDAEYYRSLVWILDNDITGVLDLTFSVEDNCFGEVVT-------IDLKPNGRNIEVTEENKREYVDLVTVWRIQKRIEEQFNAFHEGFSELIPQELINVFDERELELLIGGISEIDMEDWKKHTDYRSYSEND-----QIIKWFWELMDEWSNEKKSRLLQFTTGTSRIPVNGFKDLQGSDGPRKFTIEKAGEPNKLPKAHTCFNRLDLP"
assert record.rounds[2].alignments[2].hsps[0].query_start==606
assert record.rounds[2].alignments[2].hsps[0].query_end==865
assert record.rounds[2].alignments[2].hsps[0].sbjct_start==491
assert record.rounds[2].alignments[2].hsps[0].sbjct_end==742
assert record.rounds[2].alignments[3].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[3].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N     Q  +WFW  V  M   +R  L+ F T +  +P    A   G     S T+      + LP A+TC +RL +P"
assert record.rounds[2].alignments[3].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLHKPITLHDMESVDSEYYNSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKNGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSANH----QVIQWFWKAVLMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPEKLPRAHTCFNRLDLP"
assert record.rounds[2].alignments[3].hsps[0].query_start==606
assert record.rounds[2].alignments[3].hsps[0].query_end==865
assert record.rounds[2].alignments[3].hsps[0].sbjct_start==649
assert record.rounds[2].alignments[3].hsps[0].sbjct_end==901
assert record.rounds[2].alignments[4].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[4].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N     Q   WFW  V  M   +R  L+ F T +  +P    A   G     S T+        LP A+TC +RL +P"
assert record.rounds[2].alignments[4].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLQKLITLHDMESVDSEYYSSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKTGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSMNH----QVIHWFWKAVWMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPDKLPRAHTCFNRLDLP"
assert record.rounds[2].alignments[4].hsps[0].query_start==606
assert record.rounds[2].alignments[4].hsps[0].query_end==865
assert record.rounds[2].alignments[4].hsps[0].sbjct_start==679
assert record.rounds[2].alignments[4].hsps[0].sbjct_end==931
assert record.rounds[2].alignments[5].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRPPDD-QHLPTANTCISRLYVP"
assert record.rounds[2].alignments[5].hsps[0].match=="P  G N E  LN F+ IGR++GL +             + K++L +KV   D    D  +Y SL  ++  S     D  FSA D  F   +        V+L P+G NI VT  N  EYV  Y + R++   ++   A   G  +++P++ +      +  LL+ G  E++++     T +      +     +  +WFW  V      +R  L+ F T +  +P +     +G       TI    + Q LP ++TC +R+ +P"
assert record.rounds[2].alignments[5].hsps[0].sbjct=="PNSGINPE-HLNYFKFIGRVVGLGVFHRRFLDAFFVGALYKMMLRKKVVLQDMEGVDAEVYNSLNWMLENSIDGVLDLTFSADDERFGEVVT-------VDLKPDGRNIEVTDGNKKEYVELYTQWRIVDRVQEQFKAFMDGFNELIPEDLVTVFDERELELLIGGIAEIDIEDWKKHTDYRGYQESD-----EVIQWFWKCVSEWDNEQRARLLQFTTGTSRIPVNGFKDLQGSDGPRRFTIEKAGEVQQLPKSHTCFNRVDLP"
assert record.rounds[2].alignments[5].hsps[0].query_start==606
assert record.rounds[2].alignments[5].hsps[0].query_end==865
assert record.rounds[2].alignments[5].hsps[0].sbjct_start==533
assert record.rounds[2].alignments[5].hsps[0].sbjct_end==784
assert record.rounds[2].alignments[6].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLED-LTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[6].hsps[0].match=="F  IG +LGL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    + + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+     E+  LL+ G   ++ Q L   T +  + G   + +L   R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[2].alignments[6].hsps[0].sbjct=="FTLIGIVLGLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGNVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVNLYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRDSVL--IREFWEIVHSFTDEQKRLFLQFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[2].alignments[6].hsps[0].query_start==619
assert record.rounds[2].alignments[6].hsps[0].query_end==865
assert record.rounds[2].alignments[6].hsps[0].sbjct_start==612
assert record.rounds[2].alignments[6].hsps[0].sbjct_end==850
assert record.rounds[2].alignments[7].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLED-LTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[7].hsps[0].match=="F  IG + GL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    S + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+     E+  LL+ G   ++ Q L   T +  + G   E +    R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[2].alignments[7].hsps[0].sbjct=="FTLIGIL-GLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGSVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVISYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRESV--VIREFWEIVHSFTDEQKRLFLLFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[2].alignments[7].hsps[0].query_start==619
assert record.rounds[2].alignments[7].hsps[0].query_end==865
assert record.rounds[2].alignments[7].hsps[0].sbjct_start==623
assert record.rounds[2].alignments[7].hsps[0].sbjct_end==860
assert record.rounds[2].alignments[8].hsps[0].query=="YTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAV------DLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPK-NSLEDLTAEDFRLLVNGCGE---VNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYV"
assert record.rounds[2].alignments[8].hsps[0].match=="+    GKN++  L  +   G ++GL +  + +  +   + + K L    +++ D++   P    +L +++  ++ +  D      +  +        D    +    VEL  NG N+P+T  N +E+V K+ E  +    E   +    G   V  + NS++   +E+   LV G  E    + + L S T +     +++  +     WFW I+E      ++ L+ F T+S  +PA+     P     +   D   LP A+TC + + +"
assert record.rounds[2].alignments[8].hsps[0].sbjct=="FDKSKGKNSQLEL--YYLFGVVMGLAIFNSTILDLQFPKALYKKLCSEPLSFEDYSELFPETSRNLIKMLNYTEDNFEDVFSLTFETTYRNNNWILNDSKSSKEYVTVELCENGRNVPITQSNKHEFVMKWVEFYLEKSIEPQYNKFVSGFKRVFAECNSIKLFNSEELERLVCGDEEQTKFDFKSLRSVTKYVGGFSDDSRAVC----WFWEIIESWDYPLQKKLLQFVTASDRIPATGISTIPFKISLLGSHDSDDLPLAHTCFNEICL"
assert record.rounds[2].alignments[8].hsps[0].query_start==604
assert record.rounds[2].alignments[8].hsps[0].query_end==864
assert record.rounds[2].alignments[8].hsps[0].sbjct_start==602
assert record.rounds[2].alignments[8].hsps[0].sbjct_end==866
assert record.rounds[2].alignments[9].hsps[0].query=="TPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNW--HDFAFFDPVMYES---LRQLILASQSSDADAVFSAMDLAFAVDL-----CKEE---------GGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRP------PDDQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[9].hsps[0].match=="T +P    + ++  FR +G+++   ++   L  + L     K +L ++ +   HD    DPV+  S   L  ++   +  + D   +   L +A++      C  E         G   +EL   G +IPVT  N+ EY+R      +     +   + R G   V P + L+    E+   L+ G                 + G   +   +  ++ + I+      +++  + F T SP LP         P   +R         D  LP+  TC++ L +P"
assert record.rounds[2].alignments[9].hsps[0].sbjct=="TAKPAHIAKVKMK-FRFLGKLMAKAIMDFRLVDLPLGLPFYKWMLRQETSLTSHDLFDIDPVVARSVYHLEDIVRQKKRLEQDKSQTKESLQYALETLTMNGCSVEDLGLDFTLPGFPNIELKKGGKDIPVTIHNLEEYLRLVIFWALNEGVSRQFDSFRDGFESVFPLSHLQYFYPEELDQLLCGSKADTWDAKTLMECCRPDHGYTHDS--RAVKFLFEILSSFDNEQQRLFLQFVTGSPRLPVGGFRSLNPPLTIVRKTFESTENPDDFLPSVMTCVNYLKLP"
assert record.rounds[2].alignments[9].hsps[0].query_start==605
assert record.rounds[2].alignments[9].hsps[0].query_end==865
assert record.rounds[2].alignments[9].hsps[0].sbjct_start==1684
assert record.rounds[2].alignments[9].hsps[0].sbjct_end==1966
assert record.rounds[2].alignments[10].hsps[0].query=="GRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRPPD--DQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[10].hsps[0].match=="GR +   ++   L      R   K +LG+ V + D    D   Y+ L  L+        D      DL F+ ++ +E G  +V +L PNG NI VT +N  EYV    + RM     + L A  +G  +++PK  +   T ++  LL  G   +++  L S T ++     + +      +WFW  +      +R   + F T +  +P    A+ EG   +    I   D     LP+A+TC ++L +P"
assert record.rounds[2].alignments[10].hsps[0].sbjct=="GRYVAKAVMTTALLECYFTRSFYKHILGKSVRYTDMESEDYHFYQGLVYLLEN------DVSTLGYDLTFSTEV-QEFGVCEVRDLKPNGANILVTEENKKEYVHLVCQMRMTGAIRKQLAAFLEGFYEIIPKRLISIFTEQELELLYTGLPTIDIDDLKSNTEYHKYQSNSIQ-----IQWFWRALRSFDQADRAKFLQFVTGTSKVPLQGFAALEGMNGIQKFQIHRDDRSTDRLPSAHTCFNQLDLP"
assert record.rounds[2].alignments[10].hsps[0].query_start==623
assert record.rounds[2].alignments[10].hsps[0].query_end==865
assert record.rounds[2].alignments[10].hsps[0].sbjct_start==45
assert record.rounds[2].alignments[10].hsps[0].sbjct_end==282
assert record.rounds[2].alignments[11].hsps[0].query=="EARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGG--QVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPS-ITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[11].hsps[0].match=="    N F  IG   GL +  + +  +     + K LL  K    D     P    SL++L+      D +  F    L F +  C+E  G   Q +LIP G N+ V   N  E+V  Y  +   +   +   A   G L V     LE     + R ++ G    N + L     +  +       +    + FW       + +++  + F T S  +P    G   +   I      +++LP A+TC + L +P"
assert record.rounds[2].alignments[11].hsps[0].sbjct=="FVEHNWFHLIGITCGLAIYNSTVVDLHFPLALYKKLLNVKPGLEDLKELSPTEGRSLQELLDY-PGEDVEETFC---LNFTI--CRESYGVIEQKKLIPGGDNVTVCKDNRQEFVDAYVNYVFQISVHEWYTAFSSGFLKVCGGKVLELFQPSELRAMMVGNSNYNWEELEETAIYKGDYSATHPTV----KLFWETFHEFPLEKKKKFLLFLTGSDRIPI--YGMASLQIVIQSTASGEEYLPVAHTCYNLLDLP"
assert record.rounds[2].alignments[11].hsps[0].query_start==613
assert record.rounds[2].alignments[11].hsps[0].query_end==865
assert record.rounds[2].alignments[11].hsps[0].sbjct_start==782
assert record.rounds[2].alignments[11].hsps[0].sbjct_end==1025
assert record.rounds[2].alignments[12].hsps[0].query=="RPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLG------------RKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE-EGFQPMPSITIRPPD-----DQHLPTANTCISRLYVP"
assert record.rounds[2].alignments[12].hsps[0].match==" P  N E  +  F  +G  +   LL N +     ++   ++L               +         DP++ +SL+ ++      D +    ++ L F V      G   +ELIP G N  +   NV EY+    +  +    E+ L A  +G   V     +  L  ++  + + G  E +  M   +T+ N E G   +  +     F SI+      ER+  + F T SP LP    +   P  ++ ++  +     D++LP+  TC + L +P"
assert record.rounds[2].alignments[12].hsps[0].sbjct=="NPFSNNEKVIELFGYLGTFVARSLLDNRILDFRFSKVFFELLHRMSTPNVTTVPSDVETCLLMIELVDPLLAKSLKYIVAN---KDDNMTLESLSLTFTVP-----GNDDIELIPGGCNKSLNSSNVEEYIHGVIDQILGKGIEKQLKAFIEGFSKVFSYERMLILFPDEL-VDIFGRVEEDWSMATLYTNLNAEHGYTMDSSI--IHDFISIISAFGKHERRLFLQFLTGSPKLPIGGFKSLNPKFTVVLKHAEDGLTADEYLPSVMTCANYLKLP"
assert record.rounds[2].alignments[12].hsps[0].query_start==607
assert record.rounds[2].alignments[12].hsps[0].query_end==865
assert record.rounds[2].alignments[12].hsps[0].sbjct_start==1192
assert record.rounds[2].alignments[12].hsps[0].sbjct_end==1457
assert record.rounds[2].alignments[13].hsps[0].query=="RLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLL----GRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCG-EVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPM-PSITIRPPDDQ--HLPTANTCISRLYVP"
assert record.rounds[2].alignments[13].hsps[0].match=="+L     +G+++G CL ++ L  ++     +K LL    G   ++ D   +D V+Y +L +L+    ++D      ++DL F +D   E     V+LIPNG    VT  NV  YV K  ++++     +P+ A   GL  ++  + +E   + + ++L++G    +++  L S T +     E+     Q    FW ++      E+ + + F TS P  P   +GF+ + P   IR    +   LPTA+TC++ L +P"
assert record.rounds[2].alignments[13].hsps[0].sbjct=="KLKYIWFLGKVVGKCLYEHVLIDVSFADFFLKKLLNYSNGFLSSFSDLGSYDSVLYNNLIKLLN--MTTDE---IKSLDLTFEIDE-PESSAKVVDLIPNGSKTYVTKDNVLLYVTKVTDYKLNKRCFKPVSAFHGGLSVIIAPHWMEMFNSIELQMLISGERDNIDLDDLKSNTEYGGYKEED-----QTIVDFWEVLNEFKFEEKLNFLKFVTSVPQAP--LQGFKALDPKFGIRNAGTEKYRLPTASTCVNLLKLP"
assert record.rounds[2].alignments[13].hsps[0].query_start==615
assert record.rounds[2].alignments[13].hsps[0].query_end==865
assert record.rounds[2].alignments[13].hsps[0].sbjct_start==640
assert record.rounds[2].alignments[13].hsps[0].sbjct_end==885
assert record.rounds[2].alignments[14].hsps[0].query=="PDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[2].alignments[14].hsps[0].match=="       +Q LGERLYP ++ M    A KITGM                   +A+VEEA+ ++  H     A"
assert record.rounds[2].alignments[14].hsps[0].sbjct=="NAKPQEQKQILGERLYPMIEHMHANLAGKITGMLLEIENSELLHMIEDQEALKAKVEEAVAVLQVHRVTEPA"
assert record.rounds[2].alignments[14].hsps[0].query_start==478
assert record.rounds[2].alignments[14].hsps[0].query_end==549
assert record.rounds[2].alignments[14].hsps[0].sbjct_start==560
assert record.rounds[2].alignments[14].hsps[0].sbjct_end==631
assert record.rounds[2].alignments[15].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[2].alignments[15].hsps[0].match==" P    +Q LGERL+P +QAM P  A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[2].alignments[15].hsps[0].sbjct=="APPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[2].alignments[15].hsps[0].query_start==479
assert record.rounds[2].alignments[15].hsps[0].query_end==549
assert record.rounds[2].alignments[15].hsps[0].sbjct_start==553
assert record.rounds[2].alignments[15].hsps[0].sbjct_end==623
assert record.rounds[2].alignments[16].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[2].alignments[16].hsps[0].match==" P    +Q LGERL+P +QAM P  A KITGM                   R +V+EA+ ++ AH  +  A"
assert record.rounds[2].alignments[16].hsps[0].sbjct=="APPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRLKVDEAVAVLQAHQAKEAA"
assert record.rounds[2].alignments[16].hsps[0].query_start==479
assert record.rounds[2].alignments[16].hsps[0].query_end==549
assert record.rounds[2].alignments[16].hsps[0].sbjct_start==551
assert record.rounds[2].alignments[16].hsps[0].sbjct_end==621
assert record.rounds[2].alignments[17].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[2].alignments[17].hsps[0].match==" P    +Q LGERL+P +QAM P+ A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[2].alignments[17].hsps[0].sbjct=="APPQEQKQMLGERLFPLIQAMHPSLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[2].alignments[17].hsps[0].query_start==479
assert record.rounds[2].alignments[17].hsps[0].query_end==549
assert record.rounds[2].alignments[17].hsps[0].sbjct_start==553
assert record.rounds[2].alignments[17].hsps[0].sbjct_end==623
assert record.rounds[2].alignments[18].hsps[0].query=="PDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[2].alignments[18].hsps[0].match=="       R  LGE LYP V+ ++   A+K+TGM                   +A+V EAM+++"
assert record.rounds[2].alignments[18].hsps[0].sbjct=="NATPEQQRTMLGEVLYPLVEQVEAESAAKVTGMLLEMDQTEVLHLLESPEALKAKVAEAMDVL"
assert record.rounds[2].alignments[18].hsps[0].query_start==478
assert record.rounds[2].alignments[18].hsps[0].query_end==540
assert record.rounds[2].alignments[18].hsps[0].sbjct_start==549
assert record.rounds[2].alignments[18].hsps[0].sbjct_end==611
assert record.rounds[2].alignments[19].hsps[0].query=="PASEGNPSDDPD------PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIV"
assert record.rounds[2].alignments[19].hsps[0].match=="PA  G P           P  + +Q LGE LYP+V   +   + KITGM                     RV EA+ ++ "
assert record.rounds[2].alignments[19].hsps[0].sbjct=="PAVPGMPERFTAADLAAVPEESRKQVLGELLYPKVFVREEKLSGKITGMLLEMPNSELLELLEDDSALNERVNEAIGVLQ"
assert record.rounds[2].alignments[19].hsps[0].query_start==468
assert record.rounds[2].alignments[19].hsps[0].query_end==541
assert record.rounds[2].alignments[19].hsps[0].sbjct_start==563
assert record.rounds[2].alignments[19].hsps[0].sbjct_end==642
assert record.rounds[2].alignments[20].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[2].alignments[20].hsps[0].match=="     H + LG+ LYP V+  +PA  +K+TGM                   +A+V EA++++"
assert record.rounds[2].alignments[20].hsps[0].sbjct=="ASPDKHPRMLGDHLYPLVEQQEPANPAKVTGMLLEMDQAEILHLLESPEALKAKVSEALDVL"
assert record.rounds[2].alignments[20].hsps[0].query_start==479
assert record.rounds[2].alignments[20].hsps[0].query_end==540
assert record.rounds[2].alignments[20].hsps[0].sbjct_start==585
assert record.rounds[2].alignments[20].hsps[0].sbjct_end==646
assert record.rounds[2].alignments[21].hsps[0].query=="PASEGNPSDDPDP----LPAHRQALGERLYPRVQA--MQPAFASKITGM"
assert record.rounds[2].alignments[21].hsps[0].match=="P   G P +  D         RQALGE+LY +V A       A KITGM"
assert record.rounds[2].alignments[21].hsps[0].sbjct=="PPQGGFPRNANDNNQFYQQKQRQALGEQLYKKVSAKTSNEEAAGKITGM"
assert record.rounds[2].alignments[21].hsps[0].query_start==468
assert record.rounds[2].alignments[21].hsps[0].query_end==510
assert record.rounds[2].alignments[21].hsps[0].sbjct_start==484
assert record.rounds[2].alignments[21].hsps[0].sbjct_end==532
assert record.rounds[2].alignments[22].hsps[0].query=="LAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKE"
assert record.rounds[2].alignments[22].hsps[0].match=="+A +  S+  GP +      ++ S  A++   E+  P  T  RP  +  G V S      +W    E   R F         S  + L  F+  E"
assert record.rounds[2].alignments[22].hsps[0].sbjct=="VASQYSSTPSGPVSVSAMRAVQRSFSAQLAHNEARQPEPTYLRPVSTSYGPVPSTQSAQEQWYSPTEAQFRAFTAAHSMPTTSAQSPLTTFDTPE"
assert record.rounds[2].alignments[22].hsps[0].query_start==219
assert record.rounds[2].alignments[22].hsps[0].query_end==313
assert record.rounds[2].alignments[22].hsps[0].sbjct_start==702
assert record.rounds[2].alignments[22].hsps[0].sbjct_end==796
assert record.rounds[2].alignments[23].hsps[0].query=="SAMDLAFAVDLCKEEGGGQ---VELIPNGVNIPVTPQNVYEYVRKY-AEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCG"
assert record.rounds[2].alignments[23].hsps[0].match=="    L F +D  K         +EL P G            YV      + +    E    A ++G     P  +L  L  E     +    "
assert record.rounds[2].alignments[23].hsps[0].sbjct=="LREQLCFTIDPEKAGDFDDAVAIELTPEGY--------YKLYVHIADVSYYVREGTETDKEAYKRGFTYYFPDRALHML-PEKLSAKLCSLR"
assert record.rounds[2].alignments[23].hsps[0].query_start==686
assert record.rounds[2].alignments[23].hsps[0].query_end==773
assert record.rounds[2].alignments[23].hsps[0].sbjct_start==243
assert record.rounds[2].alignments[23].hsps[0].sbjct_end==325
assert record.rounds[3].alignments[0].hsps[0].query=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[3].alignments[0].hsps[0].match=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQ            GI            TSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQ                                      QLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGM                   RARVEEAMELIVAHGRENGA            EKVQENRKRHGSSRSVV            NAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP                   NFGFV"
assert record.rounds[3].alignments[0].hsps[0].sbjct=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQLERKRTRELLELGIDNEDSEHENDDDTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQRLRNRGERDREREREREMRRSSGLRAGSRRDRDRDFRRQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMLLELSPAQLLLLLASEDSLRARVEEAMELIVAHGRENGADSILDLGLLDSSEKVQENRKRHGSSRSVVDMDLDDTDDGDDNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPLYSSKQILKQKLLLAIKTKNFGFV"
assert record.rounds[3].alignments[0].hsps[0].query_start==1
assert record.rounds[3].alignments[0].hsps[0].query_end==889
assert record.rounds[3].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[3].alignments[0].hsps[0].sbjct_end==889
assert record.rounds[3].alignments[1].hsps[0].query=="QCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSL-EVDRDRDLLIQQTMRQLNNHFGR--RCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGN---PSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRM-LVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLL-QFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[3].alignments[1].hsps[0].match=="Q +F G+   + +L  +  L+L+LFGR FM+DVG E GS+L EL GF VKE +FRR MEKLRN Q RDL L +++R+R+ LI QT ++LN  FG   R    P+  +RVKVTFKDEPGEGSGVARSFYT+IA+A L++ K+PNL+ +Q      H+  +                                        L++D RP+ P +  +   P    D L  H Q +GERLYP++ ++    A KITGM                   R +V EA+E+I    +   +               Q ++ +   S  VV            N PLFY PGKRGFY PR G  +  R+N FRNIGR++GLCLLQNEL P+ L RHV+K +L RK+ +HD AFFDP +YES RQ+I  +Q+ + +   + M+L F +DL KEEG G  ELIP G ++ V+ + +Y       ++   +   E+ L A++ G+ DVLP NS+ +LTAED RLL+NG G++NV  LIS+T+FNDES E  +K L +FK+WFWSIVE+M++ ERQ LVYFWT SP+LPASEEGFQP+PS+TIRP DD HLPTANTCISRLY+P                   NFGFV"
assert record.rounds[3].alignments[1].hsps[0].sbjct=="QSNFSGIAPCNFLLSAK--LTLDLFGRPFMDDVGMEHGSVLPELRGFPVKEMRFRRHMEKLRNGQQRDLVLCKLERNRESLIVQTFKELNTQFGNQSRRIQPPITFNRVKVTFKDEPGEGSGVARSFYTSIAEALLASAKIPNLESVQVGTN--HSKYVVPFSSILRTTVSGSSRDQSTLQRRGSNSKILWRSARERKALNLDARPYTPPNSSDNATPESLNDHLSVHLQQIGERLYPKIHSINQTHAPKITGMLLEIPTPQLLSVISSDETLRQKVNEAIEIITFKQKSETSA--------------QSSQPKKSPSVVVVDPVDDD------NEPLFYSPGKRGFYIPRQGFASFERINAFRNIGRLIGLCLLQNELLPLFLQRHVLKYILRRKIKFHDLAFFDPALYESFRQIIQNAQTKEGEETINRMELCFVIDLMKEEGCGNRELIPGGRDV-VSHRVIYSSTSDAIQNIDXIKSQEKALEALKDGVFDVLPDNSMINLTAEDLRLLLNGVGDINVSTLISYTTFNDESSEGPDK-LXKFKKWFWSIVEKMNIMERQHLVYFWTGSPALPASEEGFQPLPSVTIRPADDSHLPTANTCISRLYIPLYSSKSILRSKMLMAIKSKNFGFV"
assert record.rounds[3].alignments[1].hsps[0].query_start==263
assert record.rounds[3].alignments[1].hsps[0].query_end==889
assert record.rounds[3].alignments[1].hsps[0].sbjct_start==2287
assert record.rounds[3].alignments[1].hsps[0].sbjct_end==2895
assert record.rounds[3].alignments[1].hsps[1].query=="SARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDE------SLPA--ETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMD---RNCLEVLPTKMSYAANLK"
assert record.rounds[3].alignments[1].hsps[1].match=="++R DF  Y LSLMRSH  EH D LPVLD+ +L+H+AYV  A +Y+++  N     D                I              + A L + + +      S+P+  +  + H FF RS+S   LGC  P  FE+PL  A+PLAD+PHLLQPN+++++LF      + +++  SG      T D    +  +  PT++ ++ +LK"
assert record.rounds[3].alignments[1].hsps[1].sbjct=="NSRRDFFTYCLSLMRSHTSEHRDALPVLDITALRHIAYVLDAFVYYMR--NDSGFYDKQDSISGR--------INNLSPMTESYDTDDELANLEEFNADVQMSASSMPSGSQGTRRHAFFARSESTLTLGCSAPEGFELPLDMAMPLADKPHLLQPNSKRQELFANLPLLVTTNANNSG-----ATNDGDGASIFDYTPTRLGFSNSLK"
assert record.rounds[3].alignments[1].hsps[1].query_start==3
assert record.rounds[3].alignments[1].hsps[1].query_end==200
assert record.rounds[3].alignments[1].hsps[1].sbjct_start==1913
assert record.rounds[3].alignments[1].hsps[1].sbjct_end==2106
assert record.rounds[3].alignments[2].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRPP-DDQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[2].hsps[0].match=="P  G N E  LN F+ IGR++GL +               K++L +KV   D    D   Y SL  ++    +   D  FS  D  F   +        ++L PNG NI VT +N  EYV      R+    E+  +A  +G  +++P+  +      +  LL+ G  E++++     T +   S  +     Q  +WFW +++  S  ++  L+ F T +  +P +     +G       TI    +   LP A+TC +RL +P"
assert record.rounds[3].alignments[2].hsps[0].sbjct=="PHSGINPE-HLNYFKFIGRVIGLAIFHRRFVDAFFVVSFYKMILQKKVTLQDMESMDAEYYRSLVWILDNDITGVLDLTFSVEDNCFGEVVT-------IDLKPNGRNIEVTEENKREYVDLVTVWRIQKRIEEQFNAFHEGFSELIPQELINVFDERELELLIGGISEIDMEDWKKHTDYRSYSEND-----QIIKWFWELMDEWSNEKKSRLLQFTTGTSRIPVNGFKDLQGSDGPRKFTIEKAGEPNKLPKAHTCFNRLDLP"
assert record.rounds[3].alignments[2].hsps[0].query_start==606
assert record.rounds[3].alignments[2].hsps[0].query_end==865
assert record.rounds[3].alignments[2].hsps[0].sbjct_start==491
assert record.rounds[3].alignments[2].hsps[0].sbjct_end==742
assert record.rounds[3].alignments[3].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[3].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N     Q  +WFW  V  M   +R  L+ F T +  +P +      G     S T+      + LP A+TC +RL +P"
assert record.rounds[3].alignments[3].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLHKPITLHDMESVDSEYYNSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKNGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSANH----QVIQWFWKAVLMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPEKLPRAHTCFNRLDLP"
assert record.rounds[3].alignments[3].hsps[0].query_start==606
assert record.rounds[3].alignments[3].hsps[0].query_end==865
assert record.rounds[3].alignments[3].hsps[0].sbjct_start==649
assert record.rounds[3].alignments[3].hsps[0].sbjct_end==901
assert record.rounds[3].alignments[4].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[4].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N     Q   WFW  V  M   +R  L+ F T +  +P +      G     S T+        LP A+TC +RL +P"
assert record.rounds[3].alignments[4].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLQKLITLHDMESVDSEYYSSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKTGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSMNH----QVIHWFWKAVWMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPDKLPRAHTCFNRLDLP"
assert record.rounds[3].alignments[4].hsps[0].query_start==606
assert record.rounds[3].alignments[4].hsps[0].query_end==865
assert record.rounds[3].alignments[4].hsps[0].sbjct_start==679
assert record.rounds[3].alignments[4].hsps[0].sbjct_end==931
assert record.rounds[3].alignments[5].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRPPDD-QHLPTANTCISRLYVP"
assert record.rounds[3].alignments[5].hsps[0].match=="P  G N E  LN F+ IGR++GL +             + K++L +KV   D    D  +Y SL  ++  S     D  FSA D  F   +        V+L P+G NI VT  N  EYV  Y + R++   ++   A   G  +++P++ +      +  LL+ G  E++++     T +      +     +  +WFW  V      +R  L+ F T +  +P +     +G       TI    + Q LP ++TC +R+ +P"
assert record.rounds[3].alignments[5].hsps[0].sbjct=="PNSGINPE-HLNYFKFIGRVVGLGVFHRRFLDAFFVGALYKMMLRKKVVLQDMEGVDAEVYNSLNWMLENSIDGVLDLTFSADDERFGEVVT-------VDLKPDGRNIEVTDGNKKEYVELYTQWRIVDRVQEQFKAFMDGFNELIPEDLVTVFDERELELLIGGIAEIDIEDWKKHTDYRGYQESD-----EVIQWFWKCVSEWDNEQRARLLQFTTGTSRIPVNGFKDLQGSDGPRRFTIEKAGEVQQLPKSHTCFNRVDLP"
assert record.rounds[3].alignments[5].hsps[0].query_start==606
assert record.rounds[3].alignments[5].hsps[0].query_end==865
assert record.rounds[3].alignments[5].hsps[0].sbjct_start==533
assert record.rounds[3].alignments[5].hsps[0].sbjct_end==784
assert record.rounds[3].alignments[6].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLED-LTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[6].hsps[0].match=="F  IG +LGL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    + + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+     E+  LL+ G   ++ Q L   T +  + G   + +L   R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[3].alignments[6].hsps[0].sbjct=="FTLIGIVLGLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGNVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVNLYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRDSVL--IREFWEIVHSFTDEQKRLFLQFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[3].alignments[6].hsps[0].query_start==619
assert record.rounds[3].alignments[6].hsps[0].query_end==865
assert record.rounds[3].alignments[6].hsps[0].sbjct_start==612
assert record.rounds[3].alignments[6].hsps[0].sbjct_end==850
assert record.rounds[3].alignments[7].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLED-LTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[7].hsps[0].match=="F  IG + GL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    S + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+     E+  LL+ G   ++ Q L   T +  + G   E +    R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[3].alignments[7].hsps[0].sbjct=="FTLIGIL-GLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGSVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVISYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRESV--VIREFWEIVHSFTDEQKRLFLLFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[3].alignments[7].hsps[0].query_start==619
assert record.rounds[3].alignments[7].hsps[0].query_end==865
assert record.rounds[3].alignments[7].hsps[0].sbjct_start==623
assert record.rounds[3].alignments[7].hsps[0].sbjct_end==860
assert record.rounds[3].alignments[8].hsps[0].query=="YTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAV------DLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPK-NSLEDLTAEDFRLLVNGCGE---VNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYV"
assert record.rounds[3].alignments[8].hsps[0].match=="+    GKN++  L  +   G ++GL +  + +  +   + + K L    +++ D++   P    +L +++  ++ +  D      +  +        D    +    VEL  NG N+P+T  N +E+V K+ E  +    E   +    G   V  + NS++   +E+   LV G  E    + + L S T +     +++  +     WFW I+E      ++ L+ F T+S  +PA+     P     +   D   LP A+TC + + +"
assert record.rounds[3].alignments[8].hsps[0].sbjct=="FDKSKGKNSQLEL--YYLFGVVMGLAIFNSTILDLQFPKALYKKLCSEPLSFEDYSELFPETSRNLIKMLNYTEDNFEDVFSLTFETTYRNNNWILNDSKSSKEYVTVELCENGRNVPITQSNKHEFVMKWVEFYLEKSIEPQYNKFVSGFKRVFAECNSIKLFNSEELERLVCGDEEQTKFDFKSLRSVTKYVGGFSDDSRAVC----WFWEIIESWDYPLQKKLLQFVTASDRIPATGISTIPFKISLLGSHDSDDLPLAHTCFNEICL"
assert record.rounds[3].alignments[8].hsps[0].query_start==604
assert record.rounds[3].alignments[8].hsps[0].query_end==864
assert record.rounds[3].alignments[8].hsps[0].sbjct_start==602
assert record.rounds[3].alignments[8].hsps[0].sbjct_end==866
assert record.rounds[3].alignments[9].hsps[0].query=="TPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNW--HDFAFFDPVMYES---LRQLILASQSSDADAVFSAMDLAFAVDL-----CKEE---------GGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRP------PDDQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[9].hsps[0].match=="T +P    + ++  FR +G+++   ++   L  + L     K +L ++ +   HD    DPV+  S   L  ++   +  + D   +   L +A++      C  E         G   +EL   G +IPVT  N+ EY+R      +     +   + R G   V P + L+    E+   L+ G                 + G   +   +  ++ + I+      +++  + F T SP LP         P   +R         D  LP+  TC++ L +P"
assert record.rounds[3].alignments[9].hsps[0].sbjct=="TAKPAHIAKVKMK-FRFLGKLMAKAIMDFRLVDLPLGLPFYKWMLRQETSLTSHDLFDIDPVVARSVYHLEDIVRQKKRLEQDKSQTKESLQYALETLTMNGCSVEDLGLDFTLPGFPNIELKKGGKDIPVTIHNLEEYLRLVIFWALNEGVSRQFDSFRDGFESVFPLSHLQYFYPEELDQLLCGSKADTWDAKTLMECCRPDHGYTHDS--RAVKFLFEILSSFDNEQQRLFLQFVTGSPRLPVGGFRSLNPPLTIVRKTFESTENPDDFLPSVMTCVNYLKLP"
assert record.rounds[3].alignments[9].hsps[0].query_start==605
assert record.rounds[3].alignments[9].hsps[0].query_end==865
assert record.rounds[3].alignments[9].hsps[0].sbjct_start==1684
assert record.rounds[3].alignments[9].hsps[0].sbjct_end==1966
assert record.rounds[3].alignments[10].hsps[0].query=="GRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLP----ASEEGFQPMPSITIRPPD--DQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[10].hsps[0].match=="GR +   ++   L      R   K +LG+ V + D    D   Y+ L  L+        D      DL F+ ++ +E G  +V +L PNG NI VT +N  EYV    + RM     + L A  +G  +++PK  +   T ++  LL  G   +++  L S T ++     + +      +WFW  +      +R   + F T +  +P    A+ EG   +    I   D     LP+A+TC ++L +P"
assert record.rounds[3].alignments[10].hsps[0].sbjct=="GRYVAKAVMTTALLECYFTRSFYKHILGKSVRYTDMESEDYHFYQGLVYLLEN------DVSTLGYDLTFSTEV-QEFGVCEVRDLKPNGANILVTEENKKEYVHLVCQMRMTGAIRKQLAAFLEGFYEIIPKRLISIFTEQELELLYTGLPTIDIDDLKSNTEYHKYQSNSIQ-----IQWFWRALRSFDQADRAKFLQFVTGTSKVPLQGFAALEGMNGIQKFQIHRDDRSTDRLPSAHTCFNQLDLP"
assert record.rounds[3].alignments[10].hsps[0].query_start==623
assert record.rounds[3].alignments[10].hsps[0].query_end==865
assert record.rounds[3].alignments[10].hsps[0].sbjct_start==45
assert record.rounds[3].alignments[10].hsps[0].sbjct_end==282
assert record.rounds[3].alignments[11].hsps[0].query=="EARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGG--QVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPS-ITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[11].hsps[0].match=="    N F  IG   GL +  + +  +     + K LL  K    D     P    SL++L+      D +  F    L F +  C+E  G   Q +LIP G N+ V   N  E+V  Y  +   +   +   A   G L V     LE     + R ++ G    N + L     +  +       +    + FW       + +++  + F T S  +P    G   +   I      +++LP A+TC + L +P"
assert record.rounds[3].alignments[11].hsps[0].sbjct=="FVEHNWFHLIGITCGLAIYNSTVVDLHFPLALYKKLLNVKPGLEDLKELSPTEGRSLQELLDY-PGEDVEETFC---LNFTI--CRESYGVIEQKKLIPGGDNVTVCKDNRQEFVDAYVNYVFQISVHEWYTAFSSGFLKVCGGKVLELFQPSELRAMMVGNSNYNWEELEETAIYKGDYSATHPTV----KLFWETFHEFPLEKKKKFLLFLTGSDRIPI--YGMASLQIVIQSTASGEEYLPVAHTCYNLLDLP"
assert record.rounds[3].alignments[11].hsps[0].query_start==613
assert record.rounds[3].alignments[11].hsps[0].query_end==865
assert record.rounds[3].alignments[11].hsps[0].sbjct_start==782
assert record.rounds[3].alignments[11].hsps[0].sbjct_end==1025
assert record.rounds[3].alignments[12].hsps[0].query=="RPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLG------------RKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE-EGFQPMPSITIRPPD-----DQHLPTANTCISRLYVP"
assert record.rounds[3].alignments[12].hsps[0].match==" P  N E  +  F  +G  +   LL N +     ++   ++L               +         DP++ +SL+ ++      D +    ++ L F V      G   +ELIP G N  +   NV EY+    +  +    E+ L A  +G   V     +  L  ++  + + G  E +  M   +T+ N E G   +  +     F SI+      ER+  + F T SP LP    +   P  ++ ++  +     D++LP+  TC + L +P"
assert record.rounds[3].alignments[12].hsps[0].sbjct=="NPFSNNEKVIELFGYLGTFVARSLLDNRILDFRFSKVFFELLHRMSTPNVTTVPSDVETCLLMIELVDPLLAKSLKYIVAN---KDDNMTLESLSLTFTVP-----GNDDIELIPGGCNKSLNSSNVEEYIHGVIDQILGKGIEKQLKAFIEGFSKVFSYERMLILFPDEL-VDIFGRVEEDWSMATLYTNLNAEHGYTMDSSI--IHDFISIISAFGKHERRLFLQFLTGSPKLPIGGFKSLNPKFTVVLKHAEDGLTADEYLPSVMTCANYLKLP"
assert record.rounds[3].alignments[12].hsps[0].query_start==607
assert record.rounds[3].alignments[12].hsps[0].query_end==865
assert record.rounds[3].alignments[12].hsps[0].sbjct_start==1192
assert record.rounds[3].alignments[12].hsps[0].sbjct_end==1457
assert record.rounds[3].alignments[13].hsps[0].query=="RLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLL----GRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCG-EVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPM-PSITIRPPDDQ--HLPTANTCISRLYVP"
assert record.rounds[3].alignments[13].hsps[0].match=="+L     +G+++G CL ++ L  ++     +K LL    G   ++ D   +D V+Y +L +L+    ++D      ++DL F +D   E     V+LIPNG    VT  NV  YV K  ++++     +P+ A   GL  ++  + +E   + + ++L++G    +++  L S T +     E+     Q    FW ++      E+ + + F TS P  P   +GF+ + P   IR    +   LPTA+TC++ L +P"
assert record.rounds[3].alignments[13].hsps[0].sbjct=="KLKYIWFLGKVVGKCLYEHVLIDVSFADFFLKKLLNYSNGFLSSFSDLGSYDSVLYNNLIKLLN--MTTDE---IKSLDLTFEIDE-PESSAKVVDLIPNGSKTYVTKDNVLLYVTKVTDYKLNKRCFKPVSAFHGGLSVIIAPHWMEMFNSIELQMLISGERDNIDLDDLKSNTEYGGYKEED-----QTIVDFWEVLNEFKFEEKLNFLKFVTSVPQAP--LQGFKALDPKFGIRNAGTEKYRLPTASTCVNLLKLP"
assert record.rounds[3].alignments[13].hsps[0].query_start==615
assert record.rounds[3].alignments[13].hsps[0].query_end==865
assert record.rounds[3].alignments[13].hsps[0].sbjct_start==640
assert record.rounds[3].alignments[13].hsps[0].sbjct_end==885
assert record.rounds[3].alignments[14].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[3].alignments[14].hsps[0].match==" P    +Q LGERL+P +QAM P  A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[3].alignments[14].hsps[0].sbjct=="APPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[3].alignments[14].hsps[0].query_start==479
assert record.rounds[3].alignments[14].hsps[0].query_end==549
assert record.rounds[3].alignments[14].hsps[0].sbjct_start==553
assert record.rounds[3].alignments[14].hsps[0].sbjct_end==623
assert record.rounds[3].alignments[15].hsps[0].query=="PDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[3].alignments[15].hsps[0].match=="       +Q LGERLYP ++ M    A KITGM                   +A+VEEA+ ++  H     A"
assert record.rounds[3].alignments[15].hsps[0].sbjct=="NAKPQEQKQILGERLYPMIEHMHANLAGKITGMLLEIENSELLHMIEDQEALKAKVEEAVAVLQVHRVTEPA"
assert record.rounds[3].alignments[15].hsps[0].query_start==478
assert record.rounds[3].alignments[15].hsps[0].query_end==549
assert record.rounds[3].alignments[15].hsps[0].sbjct_start==560
assert record.rounds[3].alignments[15].hsps[0].sbjct_end==631
assert record.rounds[3].alignments[16].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[3].alignments[16].hsps[0].match==" P    +Q LGERL+P +QAM P  A KITGM                   R +V+EA+ ++ AH  +  A"
assert record.rounds[3].alignments[16].hsps[0].sbjct=="APPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRLKVDEAVAVLQAHQAKEAA"
assert record.rounds[3].alignments[16].hsps[0].query_start==479
assert record.rounds[3].alignments[16].hsps[0].query_end==549
assert record.rounds[3].alignments[16].hsps[0].sbjct_start==551
assert record.rounds[3].alignments[16].hsps[0].sbjct_end==621
assert record.rounds[3].alignments[17].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[3].alignments[17].hsps[0].match==" P    +Q LGERL+P +QAM P+ A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[3].alignments[17].hsps[0].sbjct=="APPQEQKQMLGERLFPLIQAMHPSLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[3].alignments[17].hsps[0].query_start==479
assert record.rounds[3].alignments[17].hsps[0].query_end==549
assert record.rounds[3].alignments[17].hsps[0].sbjct_start==553
assert record.rounds[3].alignments[17].hsps[0].sbjct_end==623
assert record.rounds[3].alignments[18].hsps[0].query=="PASEGNPSDDPD------PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIV"
assert record.rounds[3].alignments[18].hsps[0].match=="PA  G P           P  + +Q LGE LYP+V   +   + KITGM                     RV EA+ ++ "
assert record.rounds[3].alignments[18].hsps[0].sbjct=="PAVPGMPERFTAADLAAVPEESRKQVLGELLYPKVFVREEKLSGKITGMLLEMPNSELLELLEDDSALNERVNEAIGVLQ"
assert record.rounds[3].alignments[18].hsps[0].query_start==468
assert record.rounds[3].alignments[18].hsps[0].query_end==541
assert record.rounds[3].alignments[18].hsps[0].sbjct_start==563
assert record.rounds[3].alignments[18].hsps[0].sbjct_end==642
assert record.rounds[3].alignments[19].hsps[0].query=="PDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[3].alignments[19].hsps[0].match=="       R  LGE LYP V+ ++   A+K+TGM                   +A+V EAM+++"
assert record.rounds[3].alignments[19].hsps[0].sbjct=="NATPEQQRTMLGEVLYPLVEQVEAESAAKVTGMLLEMDQTEVLHLLESPEALKAKVAEAMDVL"
assert record.rounds[3].alignments[19].hsps[0].query_start==478
assert record.rounds[3].alignments[19].hsps[0].query_end==540
assert record.rounds[3].alignments[19].hsps[0].sbjct_start==549
assert record.rounds[3].alignments[19].hsps[0].sbjct_end==611
assert record.rounds[3].alignments[20].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[3].alignments[20].hsps[0].match=="     H + LG+ LYP V+  +PA  +K+TGM                   +A+V EA++++"
assert record.rounds[3].alignments[20].hsps[0].sbjct=="ASPDKHPRMLGDHLYPLVEQQEPANPAKVTGMLLEMDQAEILHLLESPEALKAKVSEALDVL"
assert record.rounds[3].alignments[20].hsps[0].query_start==479
assert record.rounds[3].alignments[20].hsps[0].query_end==540
assert record.rounds[3].alignments[20].hsps[0].sbjct_start==585
assert record.rounds[3].alignments[20].hsps[0].sbjct_end==646
assert record.rounds[3].alignments[21].hsps[0].query=="PASEGNPSDDPDP----LPAHRQALGERLYPRVQAM--QPAFASKITGM"
assert record.rounds[3].alignments[21].hsps[0].match=="P   G P +  D         RQALGE+LY +V A       A KITGM"
assert record.rounds[3].alignments[21].hsps[0].sbjct=="PPQGGFPRNANDNNQFYQQKQRQALGEQLYKKVSAKTSNEEAAGKITGM"
assert record.rounds[3].alignments[21].hsps[0].query_start==468
assert record.rounds[3].alignments[21].hsps[0].query_end==510
assert record.rounds[3].alignments[21].hsps[0].sbjct_start==484
assert record.rounds[3].alignments[21].hsps[0].sbjct_end==532
assert record.rounds[3].alignments[22].hsps[0].query=="SAMDLAFAVDLCKEEGGGQ---VELIPNGVNIPVTPQNVYEYVRKY-AEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCG"
assert record.rounds[3].alignments[22].hsps[0].match=="    L F +D  K         +EL P G            YV      + +    E    A ++G     P  +L  L  E     +    "
assert record.rounds[3].alignments[22].hsps[0].sbjct=="LREQLCFTIDPEKAGDFDDAVAIELTPEGY--------YKLYVHIADVSYYVREGTETDKEAYKRGFTYYFPDRALHML-PEKLSAKLCSLR"
assert record.rounds[3].alignments[22].hsps[0].query_start==686
assert record.rounds[3].alignments[22].hsps[0].query_end==773
assert record.rounds[3].alignments[22].hsps[0].sbjct_start==243
assert record.rounds[3].alignments[22].hsps[0].sbjct_end==325
assert record.rounds[3].alignments[23].hsps[0].query=="LAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKE"
assert record.rounds[3].alignments[23].hsps[0].match=="+A +  S+  GP +      ++ S  A++   E+  P  T  RP  +  G V S      +W    E   R F         S  + L  F+  E"
assert record.rounds[3].alignments[23].hsps[0].sbjct=="VASQYSSTPSGPVSVSAMRAVQRSFSAQLAHNEARQPEPTYLRPVSTSYGPVPSTQSAQEQWYSPTEAQFRAFTAAHSMPTTSAQSPLTTFDTPE"
assert record.rounds[3].alignments[23].hsps[0].query_start==219
assert record.rounds[3].alignments[23].hsps[0].query_end==313
assert record.rounds[3].alignments[23].hsps[0].sbjct_start==702
assert record.rounds[3].alignments[23].hsps[0].sbjct_end==796
assert record.rounds[4].alignments[0].hsps[0].query=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[4].alignments[0].hsps[0].match=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQ            GI            TSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQ                                      QLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGM                   RARVEEAMELIVAHGRENGA            EKVQENRKRHGSSRSVV            NAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP                   NFGFV"
assert record.rounds[4].alignments[0].hsps[0].sbjct=="MMSARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQLERKRTRELLELGIDNEDSEHENDDDTSQSATLNDKDDESLPAETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMDRNCLEVLPTKMSYAANLKNVMNMQNRQKKAGEDQSMLAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSLEVDRDRDLLIQQTMRQLNNHFGRRCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQRLRNRGERDREREREREMRRSSGLRAGSRRDRDRDFRRQLSIDTRPFRPASEGNPSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMLLELSPAQLLLLLASEDSLRARVEEAMELIVAHGRENGADSILDLGLLDSSEKVQENRKRHGSSRSVVDMDLDDTDDGDDNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPLYSSKQILKQKLLLAIKTKNFGFV"
assert record.rounds[4].alignments[0].hsps[0].query_start==1
assert record.rounds[4].alignments[0].hsps[0].query_end==889
assert record.rounds[4].alignments[0].hsps[0].sbjct_start==1
assert record.rounds[4].alignments[0].hsps[0].sbjct_end==889
assert record.rounds[4].alignments[1].hsps[0].query=="QCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKESKFRREMEKLRNQQSRDLSL-EVDRDRDLLIQQTMRQLNNHFGR--RCATTPMAVHRVKVTFKDEPGEGSGVARSFYTAIAQAFLSNEKLPNLDCIQNANKGTHTSLMQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQLSIDTRPFRPASEGN---PSDDPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGAXXXXXXXXXXXXEKVQENRKRHGSSRSVVXXXXXXXXXXXXNAPLFYQPGKRGFYTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRM-LVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLL-QFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVPXXXXXXXXXXXXXXXXXXXNFGFV"
assert record.rounds[4].alignments[1].hsps[0].match=="Q +F G+   + +L  +  L+L+LFGR FM+DVG E GS+L EL GF VKE +FRR MEKLRN Q RDL L +++R+R+ LI QT ++LN  FG   R    P+  +RVKVTFKDEPGEGSGVARSFYT+IA+A L++ K+PNL+ +Q      H+  +                                        L++D RP+ P +  +   P    D L  H Q +GERLYP++ ++    A KITGM                   R +V EA+E+I    +   +               Q ++ +   S  VV            N PLFY PGKRGFY PR G  +  R+N FRNIGR++GLCLLQNEL P+ L RHV+K +L RK+ +HD AFFDP +YES RQ+I  +Q+ + +   + M+L F +DL KEEG G  ELIP G ++ V+ + +Y       ++   +   E+ L A++ G+ DVLP NS+ +LTAED RLL+NG G++NV  LIS+T+FNDES E  +K L +FK+WFWSIVE+M++ ERQ LVYFWT SP+LPASEEGFQP+PS+TIRP DD HLPTANTCISRLY+P                   NFGFV"
assert record.rounds[4].alignments[1].hsps[0].sbjct=="QSNFSGIAPCNFLLSAK--LTLDLFGRPFMDDVGMEHGSVLPELRGFPVKEMRFRRHMEKLRNGQQRDLVLCKLERNRESLIVQTFKELNTQFGNQSRRIQPPITFNRVKVTFKDEPGEGSGVARSFYTSIAEALLASAKIPNLESVQVGTN--HSKYVVPFSSILRTTVSGSSRDQSTLQRRGSNSKILWRSARERKALNLDARPYTPPNSSDNATPESLNDHLSVHLQQIGERLYPKIHSINQTHAPKITGMLLEIPTPQLLSVISSDETLRQKVNEAIEIITFKQKSETSA--------------QSSQPKKSPSVVVVDPVDDD------NEPLFYSPGKRGFYIPRQGFASFERINAFRNIGRLIGLCLLQNELLPLFLQRHVLKYILRRKIKFHDLAFFDPALYESFRQIIQNAQTKEGEETINRMELCFVIDLMKEEGCGNRELIPGGRDV-VSHRVIYSSTSDAIQNIDXIKSQEKALEALKDGVFDVLPDNSMINLTAEDLRLLLNGVGDINVSTLISYTTFNDESSEGPDK-LXKFKKWFWSIVEKMNIMERQHLVYFWTGSPALPASEEGFQPLPSVTIRPADDSHLPTANTCISRLYIPLYSSKSILRSKMLMAIKSKNFGFV"
assert record.rounds[4].alignments[1].hsps[0].query_start==263
assert record.rounds[4].alignments[1].hsps[0].query_end==889
assert record.rounds[4].alignments[1].hsps[0].sbjct_start==2287
assert record.rounds[4].alignments[1].hsps[0].sbjct_end==2895
assert record.rounds[4].alignments[1].hsps[1].query=="SARGDFLNYALSLMRSHNDEHSDVLPVLDVCSLKHVAYVFQALIYWIKAMNQQTTLDTPQXXXXXXXXXXXXGIXXXXXXXXXXXXTSQSATLNDKDDE------SLPA--ETGQNHPFFRRSDSMTFLGCIPPNPFEVPLAEAIPLADQPHLLQPNARKEDLFGRPSQGLYSSSAGSGKCLVEVTMD---RNCLEVLPTKMSYAANLK"
assert record.rounds[4].alignments[1].hsps[1].match=="++R DF  Y LSLMRSH  EH D LPVLD+ +L+H+AYV  A +Y+++  N     D                I              + A L + + +      S+P+  +  + H FF RS+S   LGC  P  FE+PL  A+PLAD+PHLLQPN+++++LF      + +++  SG      T D    +  +  PT++ ++ +LK"
assert record.rounds[4].alignments[1].hsps[1].sbjct=="NSRRDFFTYCLSLMRSHTSEHRDALPVLDITALRHIAYVLDAFVYYMR--NDSGFYDKQDSISGR--------INNLSPMTESYDTDDELANLEEFNADVQMSASSMPSGSQGTRRHAFFARSESTLTLGCSAPEGFELPLDMAMPLADKPHLLQPNSKRQELFANLPLLVTTNANNSG-----ATNDGDGASIFDYTPTRLGFSNSLK"
assert record.rounds[4].alignments[1].hsps[1].query_start==3
assert record.rounds[4].alignments[1].hsps[1].query_end==200
assert record.rounds[4].alignments[1].hsps[1].sbjct_start==1913
assert record.rounds[4].alignments[1].hsps[1].sbjct_end==2106
assert record.rounds[4].alignments[2].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRPP-DDQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[2].hsps[0].match=="P  G N E  LN F+ IGR++GL +               K++L +KV   D    D   Y SL  ++    +   D  FS  D  F   +        ++L PNG NI VT +N  EYV      R+    E+  +A  +G  +++P+  +      +  LL+ G  E++++     T +   S  +     Q  +WFW +++  S  ++  L+ F T +  +P +     +G       TI    +   LP A+TC +RL +P"
assert record.rounds[4].alignments[2].hsps[0].sbjct=="PHSGINPE-HLNYFKFIGRVIGLAIFHRRFVDAFFVVSFYKMILQKKVTLQDMESMDAEYYRSLVWILDNDITGVLDLTFSVEDNCFGEVVT-------IDLKPNGRNIEVTEENKREYVDLVTVWRIQKRIEEQFNAFHEGFSELIPQELINVFDERELELLIGGISEIDMEDWKKHTDYRSYSEND-----QIIKWFWELMDEWSNEKKSRLLQFTTGTSRIPVNGFKDLQGSDGPRKFTIEKAGEPNKLPKAHTCFNRLDLP"
assert record.rounds[4].alignments[2].hsps[0].query_start==606
assert record.rounds[4].alignments[2].hsps[0].query_end==865
assert record.rounds[4].alignments[2].hsps[0].sbjct_start==491
assert record.rounds[4].alignments[2].hsps[0].sbjct_end==742
assert record.rounds[4].alignments[3].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[3].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N     Q  +WFW  V  M   +R  L+ F T +  +P +      G     S T+      + LP A+TC +RL +P"
assert record.rounds[4].alignments[3].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLHKPITLHDMESVDSEYYNSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKNGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSANH----QVIQWFWKAVLMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPEKLPRAHTCFNRLDLP"
assert record.rounds[4].alignments[3].hsps[0].query_start==606
assert record.rounds[4].alignments[3].hsps[0].query_end==865
assert record.rounds[4].alignments[3].hsps[0].sbjct_start==649
assert record.rounds[4].alignments[3].hsps[0].sbjct_end==901
assert record.rounds[4].alignments[4].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV---ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRP-PDDQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[4].hsps[0].match=="P  G   E  L+ F+ IGR+ G+ +   +L      R   K++L + +  HD    D   Y SLR ++    +         +DL F +D   EE  GQ    EL   G  I VT +N  EY+    + R +   ++ + A ++G  +++P++ ++     +  LL+ G G+V+V      T + +    N     Q   WFW  V  M   +R  L+ F T +  +P +      G     S T+        LP A+TC +RL +P"
assert record.rounds[4].alignments[4].hsps[0].sbjct=="PNSGLCNEDHLSYFKFIGRVAGMAVYHGKLLDGFFIRPFYKMMLQKLITLHDMESVDSEYYSSLRWILENDPTE--------LDLRFIID---EELFGQTHQHELKTGGSEIVVTNKNKKEYIYLVIQWRFVNRIQKQMAAFKEGFFELIPQDLIKIFDENELELLMCGLGDVDVNDWREHTKYKNGYSMNH----QVIHWFWKAVWMMDSEKRIRLLQFVTGTSRVPMNGFAELYGSNGPQSFTVEQWGTPDKLPRAHTCFNRLDLP"
assert record.rounds[4].alignments[4].hsps[0].query_start==606
assert record.rounds[4].alignments[4].hsps[0].query_end==865
assert record.rounds[4].alignments[4].hsps[0].sbjct_start==679
assert record.rounds[4].alignments[4].hsps[0].sbjct_end==931
assert record.rounds[4].alignments[5].hsps[0].query=="PRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRPPDD-QHLPTANTCISRLYVP"
assert record.rounds[4].alignments[5].hsps[0].match=="P  G N E  LN F+ IGR++GL +             + K++L +KV   D    D  +Y SL  ++  S     D  FSA D  F   +        V+L P+G NI VT  N  EYV  Y + R++   ++   A   G  +++P++ +      +  LL+ G  E++++     T +      +     +  +WFW  V      +R  L+ F T +  +P +     +G       TI    + Q LP ++TC +R+ +P"
assert record.rounds[4].alignments[5].hsps[0].sbjct=="PNSGINPE-HLNYFKFIGRVVGLGVFHRRFLDAFFVGALYKMMLRKKVVLQDMEGVDAEVYNSLNWMLENSIDGVLDLTFSADDERFGEVVT-------VDLKPDGRNIEVTDGNKKEYVELYTQWRIVDRVQEQFKAFMDGFNELIPEDLVTVFDERELELLIGGIAEIDIEDWKKHTDYRGYQESD-----EVIQWFWKCVSEWDNEQRARLLQFTTGTSRIPVNGFKDLQGSDGPRRFTIEKAGEVQQLPKSHTCFNRVDLP"
assert record.rounds[4].alignments[5].hsps[0].query_start==606
assert record.rounds[4].alignments[5].hsps[0].query_end==865
assert record.rounds[4].alignments[5].hsps[0].sbjct_start==533
assert record.rounds[4].alignments[5].hsps[0].sbjct_end==784
assert record.rounds[4].alignments[6].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLED-LTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[6].hsps[0].match=="F  IG +LGL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    + + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+     E+  LL+ G   ++ Q L   T +  + G   + +L   R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[4].alignments[6].hsps[0].sbjct=="FTLIGIVLGLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGNVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVNLYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRDSVL--IREFWEIVHSFTDEQKRLFLQFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[4].alignments[6].hsps[0].query_start==619
assert record.rounds[4].alignments[6].hsps[0].query_end==865
assert record.rounds[4].alignments[6].hsps[0].sbjct_start==612
assert record.rounds[4].alignments[6].hsps[0].sbjct_end==850
assert record.rounds[4].alignments[7].hsps[0].query=="FRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLED-LTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[7].hsps[0].match=="F  IG + GL +  N +  +     V + L+G+K  + D     PV+Y+SL+ L+    S + D     M + F +      G   + +L  NG  IP+T +N  E+V  Y+++ +    E+   A R+G   V  ++ L+     E+  LL+ G   ++ Q L   T +  + G   E +    R FW IV   +  +++  + F T +   P    G   M  I    PD + LPT++TC + L +P"
assert record.rounds[4].alignments[7].hsps[0].sbjct=="FTLIGIL-GLAIYNNCILDVHFPMVVYRKLMGKKGTFRDLGDSHPVLYQSLKDLLEYEGSVEDD-----MMITFQISQTDLFGNPMMYDLKENGDKIPITNENRKEFVISYSDYILNKSVEKQFKAFRRGFHMVTNESPLKYLFRPEEIELLICGSRNLDFQALEETTEY--DGGYTRESV--VIREFWEIVHSFTDEQKRLFLLFTTGTDRAPVGGLGKLKM-IIAKNGPDTERLPTSHTCFNVLLLP"
assert record.rounds[4].alignments[7].hsps[0].query_start==619
assert record.rounds[4].alignments[7].hsps[0].query_end==865
assert record.rounds[4].alignments[7].hsps[0].sbjct_start==623
assert record.rounds[4].alignments[7].hsps[0].sbjct_end==860
assert record.rounds[4].alignments[8].hsps[0].query=="YTPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAV------DLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPK-NSLEDLTAEDFRLLVNGCGE---VNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQHLPTANTCISRLYV"
assert record.rounds[4].alignments[8].hsps[0].match=="+    GKN++  L  +   G ++GL +  + +  +   + + K L    +++ D++   P    +L +++  ++ +  D      +  +        D    +    VEL  NG N+P+T  N +E+V K+ E  +    E   +    G   V  + NS++   +E+   LV G  E    + + L S T +     +++  +     WFW I+E      ++ L+ F T+S  +PA+     P     +   D   LP A+TC + + +"
assert record.rounds[4].alignments[8].hsps[0].sbjct=="FDKSKGKNSQLEL--YYLFGVVMGLAIFNSTILDLQFPKALYKKLCSEPLSFEDYSELFPETSRNLIKMLNYTEDNFEDVFSLTFETTYRNNNWILNDSKSSKEYVTVELCENGRNVPITQSNKHEFVMKWVEFYLEKSIEPQYNKFVSGFKRVFAECNSIKLFNSEELERLVCGDEEQTKFDFKSLRSVTKYVGGFSDDSRAVC----WFWEIIESWDYPLQKKLLQFVTASDRIPATGISTIPFKISLLGSHDSDDLPLAHTCFNEICL"
assert record.rounds[4].alignments[8].hsps[0].query_start==604
assert record.rounds[4].alignments[8].hsps[0].query_end==864
assert record.rounds[4].alignments[8].hsps[0].sbjct_start==602
assert record.rounds[4].alignments[8].hsps[0].sbjct_end==866
assert record.rounds[4].alignments[9].hsps[0].query=="TPRPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNW--HDFAFFDPVMYES---LRQLILASQSSDADAVFSAMDLAFAVDL-----CKEE---------GGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRP------PDDQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[9].hsps[0].match=="T +P    + ++  FR +G+++   ++   L  + L     K +L ++ +   HD    DPV+  S   L  ++   +  + D   +   L +A++      C  E         G   +EL   G +IPVT  N+ EY+R      +     +   + R G   V P + L+    E+   L+ G                 + G   +   +  ++ + I+      +++  + F T SP LP         P   +R         D  LP+  TC++ L +P"
assert record.rounds[4].alignments[9].hsps[0].sbjct=="TAKPAHIAKVKMK-FRFLGKLMAKAIMDFRLVDLPLGLPFYKWMLRQETSLTSHDLFDIDPVVARSVYHLEDIVRQKKRLEQDKSQTKESLQYALETLTMNGCSVEDLGLDFTLPGFPNIELKKGGKDIPVTIHNLEEYLRLVIFWALNEGVSRQFDSFRDGFESVFPLSHLQYFYPEELDQLLCGSKADTWDAKTLMECCRPDHGYTHDS--RAVKFLFEILSSFDNEQQRLFLQFVTGSPRLPVGGFRSLNPPLTIVRKTFESTENPDDFLPSVMTCVNYLKLP"
assert record.rounds[4].alignments[9].hsps[0].query_start==605
assert record.rounds[4].alignments[9].hsps[0].query_end==865
assert record.rounds[4].alignments[9].hsps[0].sbjct_start==1684
assert record.rounds[4].alignments[9].hsps[0].sbjct_end==1966
assert record.rounds[4].alignments[10].hsps[0].query=="GRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQV-ELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE----EGFQPMPSITIRPPD--DQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[10].hsps[0].match=="GR +   ++   L      R   K +LG+ V + D    D   Y+ L  L+        D      DL F+ ++ +E G  +V +L PNG NI VT +N  EYV    + RM     + L A  +G  +++PK  +   T ++  LL  G   +++  L S T ++     + +      +WFW  +      +R   + F T +  +P       EG   +    I   D     LP+A+TC ++L +P"
assert record.rounds[4].alignments[10].hsps[0].sbjct=="GRYVAKAVMTTALLECYFTRSFYKHILGKSVRYTDMESEDYHFYQGLVYLLEN------DVSTLGYDLTFSTEV-QEFGVCEVRDLKPNGANILVTEENKKEYVHLVCQMRMTGAIRKQLAAFLEGFYEIIPKRLISIFTEQELELLYTGLPTIDIDDLKSNTEYHKYQSNSIQ-----IQWFWRALRSFDQADRAKFLQFVTGTSKVPLQGFAALEGMNGIQKFQIHRDDRSTDRLPSAHTCFNQLDLP"
assert record.rounds[4].alignments[10].hsps[0].query_start==623
assert record.rounds[4].alignments[10].hsps[0].query_end==865
assert record.rounds[4].alignments[10].hsps[0].sbjct_start==45
assert record.rounds[4].alignments[10].hsps[0].sbjct_end==282
assert record.rounds[4].alignments[11].hsps[0].query=="EARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLGRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGG--QVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPS-ITIRPPDDQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[11].hsps[0].match=="    N F  IG   GL +  + +  +     + K LL  K    D     P    SL++L+      D +  F    L F +  C+E  G   Q +LIP G N+ V   N  E+V  Y  +   +   +   A   G L V     LE     + R ++ G    N + L     +  +       +    + FW       + +++  + F T S  +P    G   +   I      +++LP A+TC + L +P"
assert record.rounds[4].alignments[11].hsps[0].sbjct=="FVEHNWFHLIGITCGLAIYNSTVVDLHFPLALYKKLLNVKPGLEDLKELSPTEGRSLQELLDY-PGEDVEETFC---LNFTI--CRESYGVIEQKKLIPGGDNVTVCKDNRQEFVDAYVNYVFQISVHEWYTAFSSGFLKVCGGKVLELFQPSELRAMMVGNSNYNWEELEETAIYKGDYSATHPTV----KLFWETFHEFPLEKKKKFLLFLTGSDRIPI--YGMASLQIVIQSTASGEEYLPVAHTCYNLLDLP"
assert record.rounds[4].alignments[11].hsps[0].query_start==613
assert record.rounds[4].alignments[11].hsps[0].query_end==865
assert record.rounds[4].alignments[11].hsps[0].sbjct_start==782
assert record.rounds[4].alignments[11].hsps[0].sbjct_end==1025
assert record.rounds[4].alignments[12].hsps[0].query=="RPGKNTEARLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLLG------------RKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCGEVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASE-EGFQPMPSITIRPPD-----DQHLPTANTCISRLYVP"
assert record.rounds[4].alignments[12].hsps[0].match==" P  N E  +  F  +G  +   LL N +     ++   ++L               +         DP++ +SL+ ++      D +    ++ L F V      G   +ELIP G N  +   NV EY+    +  +    E+ L A  +G   V     +  L  ++  + + G  E +  M   +T+ N E G   +  +     F SI+      ER+  + F T SP LP    +   P  ++ ++  +     D++LP+  TC + L +P"
assert record.rounds[4].alignments[12].hsps[0].sbjct=="NPFSNNEKVIELFGYLGTFVARSLLDNRILDFRFSKVFFELLHRMSTPNVTTVPSDVETCLLMIELVDPLLAKSLKYIVAN---KDDNMTLESLSLTFTVP-----GNDDIELIPGGCNKSLNSSNVEEYIHGVIDQILGKGIEKQLKAFIEGFSKVFSYERMLILFPDEL-VDIFGRVEEDWSMATLYTNLNAEHGYTMDSSI--IHDFISIISAFGKHERRLFLQFLTGSPKLPIGGFKSLNPKFTVVLKHAEDGLTADEYLPSVMTCANYLKLP"
assert record.rounds[4].alignments[12].hsps[0].query_start==607
assert record.rounds[4].alignments[12].hsps[0].query_end==865
assert record.rounds[4].alignments[12].hsps[0].sbjct_start==1192
assert record.rounds[4].alignments[12].hsps[0].sbjct_end==1457
assert record.rounds[4].alignments[13].hsps[0].query=="RLNCFRNIGRILGLCLLQNELCPITLNRHVIKVLL----GRKVNWHDFAFFDPVMYESLRQLILASQSSDADAVFSAMDLAFAVDLCKEEGGGQVELIPNGVNIPVTPQNVYEYVRKYAEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCG-EVNVQMLISFTSFNDESGENAEKLLQFKRWFWSIVERMSMTERQDLVYFWTSSPSLPASEEGFQPMPSITIRPPDDQ--HLPTANTCISRLYVP"
assert record.rounds[4].alignments[13].hsps[0].match=="+L     +G+++G CL ++ L  ++     +K LL    G   ++ D   +D V+Y +L +L+    ++D      ++DL F +D   E     V+LIPNG    VT  NV  YV K  ++++     +P+ A   GL  ++  + +E   + + ++L++G    +++  L S T +     E+     Q    FW ++      E+ + + F TS P  P         P   IR    +   LPTA+TC++ L +P"
assert record.rounds[4].alignments[13].hsps[0].sbjct=="KLKYIWFLGKVVGKCLYEHVLIDVSFADFFLKKLLNYSNGFLSSFSDLGSYDSVLYNNLIKLLN--MTTDE---IKSLDLTFEIDE-PESSAKVVDLIPNGSKTYVTKDNVLLYVTKVTDYKLNKRCFKPVSAFHGGLSVIIAPHWMEMFNSIELQMLISGERDNIDLDDLKSNTEYGGYKEED-----QTIVDFWEVLNEFKFEEKLNFLKFVTSVPQAPLQGFKALD-PKFGIRNAGTEKYRLPTASTCVNLLKLP"
assert record.rounds[4].alignments[13].hsps[0].query_start==615
assert record.rounds[4].alignments[13].hsps[0].query_end==865
assert record.rounds[4].alignments[13].hsps[0].sbjct_start==640
assert record.rounds[4].alignments[13].hsps[0].sbjct_end==885
assert record.rounds[4].alignments[14].hsps[0].query=="DPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[4].alignments[14].hsps[0].match=="        +Q LGERLYP ++ M    A KITGM                   +A+VEEA+ ++  H     A"
assert record.rounds[4].alignments[14].hsps[0].sbjct=="ANAKPQEQKQILGERLYPMIEHMHANLAGKITGMLLEIENSELLHMIEDQEALKAKVEEAVAVLQVHRVTEPA"
assert record.rounds[4].alignments[14].hsps[0].query_start==477
assert record.rounds[4].alignments[14].hsps[0].query_end==549
assert record.rounds[4].alignments[14].hsps[0].sbjct_start==559
assert record.rounds[4].alignments[14].hsps[0].sbjct_end==631
assert record.rounds[4].alignments[15].hsps[0].query=="DPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[4].alignments[15].hsps[0].match=="   P    +Q LGERL+P +QAM P  A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[4].alignments[15].hsps[0].sbjct=="ASAPPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[4].alignments[15].hsps[0].query_start==477
assert record.rounds[4].alignments[15].hsps[0].query_end==549
assert record.rounds[4].alignments[15].hsps[0].sbjct_start==551
assert record.rounds[4].alignments[15].hsps[0].sbjct_end==623
assert record.rounds[4].alignments[16].hsps[0].query=="DPDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[4].alignments[16].hsps[0].match=="   P    +Q LGERL+P +QAM P+ A KITGM                   R++V+EA+ ++ AH  +  A"
assert record.rounds[4].alignments[16].hsps[0].sbjct=="ASAPPQEQKQMLGERLFPLIQAMHPSLAGKITGMLLEIDNSELLHMLESPESLRSKVDEAVAVLQAHQAKEAA"
assert record.rounds[4].alignments[16].hsps[0].query_start==477
assert record.rounds[4].alignments[16].hsps[0].query_end==549
assert record.rounds[4].alignments[16].hsps[0].sbjct_start==551
assert record.rounds[4].alignments[16].hsps[0].sbjct_end==623
assert record.rounds[4].alignments[17].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIVAHGRENGA"
assert record.rounds[4].alignments[17].hsps[0].match==" P    +Q LGERL+P +QAM P  A KITGM                   R +V+EA+ ++ AH  +  A"
assert record.rounds[4].alignments[17].hsps[0].sbjct=="APPQEQKQMLGERLFPLIQAMHPTLAGKITGMLLEIDNSELLHMLESPESLRLKVDEAVAVLQAHQAKEAA"
assert record.rounds[4].alignments[17].hsps[0].query_start==479
assert record.rounds[4].alignments[17].hsps[0].query_end==549
assert record.rounds[4].alignments[17].hsps[0].sbjct_start==551
assert record.rounds[4].alignments[17].hsps[0].sbjct_end==621
assert record.rounds[4].alignments[18].hsps[0].query=="PASEGNPSDDPD------PLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELIV"
assert record.rounds[4].alignments[18].hsps[0].match=="PA  G P           P  + +Q LGE LYP+V   +   + KITGM                     RV EA+ ++ "
assert record.rounds[4].alignments[18].hsps[0].sbjct=="PAVPGMPERFTAADLAAVPEESRKQVLGELLYPKVFVREEKLSGKITGMLLEMPNSELLELLEDDSALNERVNEAIGVLQ"
assert record.rounds[4].alignments[18].hsps[0].query_start==468
assert record.rounds[4].alignments[18].hsps[0].query_end==541
assert record.rounds[4].alignments[18].hsps[0].sbjct_start==563
assert record.rounds[4].alignments[18].hsps[0].sbjct_end==642
assert record.rounds[4].alignments[19].hsps[0].query=="PDPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[4].alignments[19].hsps[0].match=="       R  LGE LYP V+ ++   A+K+TGM                   +A+V EAM+++"
assert record.rounds[4].alignments[19].hsps[0].sbjct=="NATPEQQRTMLGEVLYPLVEQVEAESAAKVTGMLLEMDQTEVLHLLESPEALKAKVAEAMDVL"
assert record.rounds[4].alignments[19].hsps[0].query_start==478
assert record.rounds[4].alignments[19].hsps[0].query_end==540
assert record.rounds[4].alignments[19].hsps[0].sbjct_start==549
assert record.rounds[4].alignments[19].hsps[0].sbjct_end==611
assert record.rounds[4].alignments[20].hsps[0].query=="DPLPAHRQALGERLYPRVQAMQPAFASKITGMXXXXXXXXXXXXXXXXXXXRARVEEAMELI"
assert record.rounds[4].alignments[20].hsps[0].match=="     H + LG+ LYP V+  +PA  +K+TGM                   +A+V EA++++"
assert record.rounds[4].alignments[20].hsps[0].sbjct=="ASPDKHPRMLGDHLYPLVEQQEPANPAKVTGMLLEMDQAEILHLLESPEALKAKVSEALDVL"
assert record.rounds[4].alignments[20].hsps[0].query_start==479
assert record.rounds[4].alignments[20].hsps[0].query_end==540
assert record.rounds[4].alignments[20].hsps[0].sbjct_start==585
assert record.rounds[4].alignments[20].hsps[0].sbjct_end==646
assert record.rounds[4].alignments[21].hsps[0].query=="PASEGNPSDDPDP----LPAHRQALGERLYPRVQAM--QPAFASKITGM"
assert record.rounds[4].alignments[21].hsps[0].match=="P   G P +  D         RQALGE+LY +V A       A KITGM"
assert record.rounds[4].alignments[21].hsps[0].sbjct=="PPQGGFPRNANDNNQFYQQKQRQALGEQLYKKVSAKTSNEEAAGKITGM"
assert record.rounds[4].alignments[21].hsps[0].query_start==468
assert record.rounds[4].alignments[21].hsps[0].query_end==510
assert record.rounds[4].alignments[21].hsps[0].sbjct_start==484
assert record.rounds[4].alignments[21].hsps[0].sbjct_end==532
assert record.rounds[4].alignments[22].hsps[0].query=="SAMDLAFAVDLCKEEGGGQ---VELIPNGVNIPVTPQNVYEYVRKY-AEHRMLVVAEQPLHAMRKGLLDVLPKNSLEDLTAEDFRLLVNGCG"
assert record.rounds[4].alignments[22].hsps[0].match=="    L F +D  K         +EL P G            YV      + +    E    A ++G     P  +L  L  E     +    "
assert record.rounds[4].alignments[22].hsps[0].sbjct=="LREQLCFTIDPEKAGDFDDAVAIELTPEGY--------YKLYVHIADVSYYVREGTETDKEAYKRGFTYYFPDRALHML-PEKLSAKLCSLR"
assert record.rounds[4].alignments[22].hsps[0].query_start==686
assert record.rounds[4].alignments[22].hsps[0].query_end==773
assert record.rounds[4].alignments[22].hsps[0].sbjct_start==243
assert record.rounds[4].alignments[22].hsps[0].sbjct_end==325
assert record.rounds[4].alignments[23].hsps[0].query=="LAEEADSSKPGPSAHDVAAQLKSSLLAEIGLTESEGPPLTSFRPQCSFMGMVISHDMLLGRWRLSLELFGRVFMEDVGAEPGSILTELGGFEVKE"
assert record.rounds[4].alignments[23].hsps[0].match=="+A +  S+  GP +      ++ S  A++   E+  P  T  RP  +  G V S      +W    E   R F         S  + L  F+  E"
assert record.rounds[4].alignments[23].hsps[0].sbjct=="VASQYSSTPSGPVSVSAMRAVQRSFSAQLAHNEARQPEPTYLRPVSTSYGPVPSTQSAQEQWYSPTEAQFRAFTAAHSMPTTSAQSPLTTFDTPE"
assert record.rounds[4].alignments[23].hsps[0].query_start==219
assert record.rounds[4].alignments[23].hsps[0].query_end==313
assert record.rounds[4].alignments[23].hsps[0].sbjct_start==702
assert record.rounds[4].alignments[23].hsps[0].sbjct_end==796
assert record.database_name==['/dbase/swissprot/main/release/sp', '/dbase/swissprot/main/update/spu']
assert record.posted_date==[('Jun 21, 2000 12:39 PM',), ('Nov 3, 1999  8:09 PM',)]
assert record.num_letters_in_database==[31411157, 546183]
assert record.num_sequences_in_database==[86593, 1608]
assert record.ka_params==[0.318, 0.131, 0.369]
assert record.ka_params_gap==[0.270, 0.0458, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==221707272
assert record.num_sequences==88201
assert record.num_extends==8189338
assert record.num_good_extends==22743
assert record.num_seqs_better_e==125
assert record.hsps_no_gap==105
assert record.hsps_prelim_gapped==20
assert record.hsps_gapped==150
assert record.query_length==889
assert record.database_length==31957340
assert record.effective_hsp_length==55
assert record.effective_query_length==834
assert record.effective_database_length==27106285
assert record.effective_search_space==22606641690
assert record.effective_search_space_used==22606641690
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(69,31.3)

handle = open('Blast/bt062')
record = parser.parse(handle)
assert record.application=="BLASTN"
assert record.version=='2.0.14'
assert record.date=="Jun-29-2000"
assert record.reference==reference
assert record.query=="gi|1348916|gb|G26684|G26684 human STS\nSTS_D11570.\x01gi|1375195|gb|G26945|G26945 human STS SHGC-32699."
assert record.query_letters==285
assert record.database=="data/sts"
assert record.database_sequences==87792
assert record.database_letters==31998854
assert len(record.descriptions)==23
assert record.descriptions[0].title=="gi|1348916|gb|G26684|G26684 human STS STS_D11570. >gi|1375195|g..."
assert record.descriptions[0].score==517
assert record.descriptions[0].e==1e-146
assert record.descriptions[1].title=="gi|4516686|dbj|AU026763.1|AU026763 Rattus norvegicus, OTSUKA cl..."
assert record.descriptions[1].score==32
assert record.descriptions[1].e==1.600000
assert record.descriptions[2].title=="gi|6120827|gb|G55508.1|G55508 SHGC-100856 Human Homo sapiens ST..."
assert record.descriptions[2].score==32
assert record.descriptions[2].e==1.600000
assert record.descriptions[3].title=="gi|720683|gb|G03725|G03725 human STS WI-344."
assert record.descriptions[3].score==32
assert record.descriptions[3].e==1.600000
assert record.descriptions[4].title=="gi|5690111|gb|G54226.1|G54226 B124N23/SP6 Human Chromosome 12 H..."
assert record.descriptions[4].score==30
assert record.descriptions[4].e==6.500000
assert record.descriptions[5].title=="gi|4493307|gb|G47007.1|G47007 Z15259_1 Zebrafish AB Danio rerio..."
assert record.descriptions[5].score==30
assert record.descriptions[5].e==6.500000
assert record.descriptions[6].title=="gi|4491799|gb|G45508.1|G45508 Z24506_1 Zebrafish AB Danio rerio..."
assert record.descriptions[6].score==30
assert record.descriptions[6].e==6.500000
assert record.descriptions[7].title=="gi|6121596|gb|G56277.1|G56277 SHGC-101791 Human Homo sapiens ST..."
assert record.descriptions[7].score==30
assert record.descriptions[7].e==6.500000
assert record.descriptions[8].title=="gi|5222417|gb|G51240.1|G51240 SHGC-80720 Human Homo sapiens STS..."
assert record.descriptions[8].score==30
assert record.descriptions[8].e==6.500000
assert record.descriptions[9].title=="gi|5221977|gb|G50800.1|G50800 SHGC-83850 Human Homo sapiens STS..."
assert record.descriptions[9].score==30
assert record.descriptions[9].e==6.500000
assert record.descriptions[10].title=="gi|5224501|gb|G53324.1|G53324 SHGC-82315 Human Homo sapiens STS..."
assert record.descriptions[10].score==30
assert record.descriptions[10].e==6.500000
assert record.descriptions[11].title=="gi|4529247|gb|G48587.1|G48587 SHGC-82546 Human Homo sapiens STS..."
assert record.descriptions[11].score==30
assert record.descriptions[11].e==6.500000
assert record.descriptions[12].title=="gi|3359917|gb|G40708|G40708 Z8947 Zebrafish AB Danio rerio STS ..."
assert record.descriptions[12].score==30
assert record.descriptions[12].e==6.500000
assert record.descriptions[13].title=="gi|3359244|gb|G40035|G40035 Z13538 Zebrafish AB Danio rerio STS..."
assert record.descriptions[13].score==30
assert record.descriptions[13].e==6.500000
assert record.descriptions[14].title=="gi|1347715|gb|G25483|G25483 human STS EST334642."
assert record.descriptions[14].score==30
assert record.descriptions[14].e==6.500000
assert record.descriptions[15].title=="gi|1244262|gb|G19475|G19475 human STS SHGC-18755."
assert record.descriptions[15].score==30
assert record.descriptions[15].e==6.500000
assert record.descriptions[16].title=="gi|1232611|emb|Z51311|HS302WC9 H.sapiens (D5S2069) DNA segment ..."
assert record.descriptions[16].score==30
assert record.descriptions[16].e==6.500000
assert record.descriptions[17].title=="gi|1223022|gb|G18565|G18565 BMS485 cow Bos taurus STS genomic, ..."
assert record.descriptions[17].score==30
assert record.descriptions[17].e==6.500000
assert record.descriptions[18].title=="gi|1161779|gb|G15890|G15890 human STS CHLC.UTR_01448_M84721.P56..."
assert record.descriptions[18].score==30
assert record.descriptions[18].e==6.500000
assert record.descriptions[19].title=="gi|858803|gb|G05558|G05558 human STS WI-7105."
assert record.descriptions[19].score==30
assert record.descriptions[19].e==6.500000
assert record.descriptions[20].title=="gi|1342455|gb|G22129|G22129 human STS WI-14200."
assert record.descriptions[20].score==30
assert record.descriptions[20].e==6.500000
assert record.descriptions[21].title=="gi|1347001|gb|G24769|G24769 human STS EST129834."
assert record.descriptions[21].score==30
assert record.descriptions[21].e==6.500000
assert record.descriptions[22].title=="gi|605469|gb|L31223|HUMUT821B Human STS UT821, 3' primer bind."
assert record.descriptions[22].score==30
assert record.descriptions[22].e==6.500000
assert len(record.alignments)==23
assert record.alignments[0].title==">gi|1348916|gb|G26684|G26684 human STS STS_D11570. gb|G26945|G26945 human STS SHGC-32699."
assert record.alignments[0].length==285
assert record.alignments[1].title==">gi|4516686|dbj|AU026763.1|AU026763 Rattus norvegicus, OTSUKA clone, OT33.16/752f07, microsatellite sequence, sequence tagged site"
assert record.alignments[1].length==307
assert record.alignments[2].title==">gi|6120827|gb|G55508.1|G55508 SHGC-100856 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[2].length==711
assert record.alignments[3].title==">gi|720683|gb|G03725|G03725 human STS WI-344."
assert record.alignments[3].length==246
assert record.alignments[4].title==">gi|5690111|gb|G54226.1|G54226 B124N23/SP6 Human Chromosome 12 Homo sapiens STS genomic clone RPCI-11-B124N23 SP6, sequence tagged site"
assert record.alignments[4].length==550
assert record.alignments[5].title==">gi|4493307|gb|G47007.1|G47007 Z15259_1 Zebrafish AB Danio rerio STS genomic clone Z15259 5', sequence tagged site"
assert record.alignments[5].length==442
assert record.alignments[6].title==">gi|4491799|gb|G45508.1|G45508 Z24506_1 Zebrafish AB Danio rerio STS genomic clone Z24506 5', sequence tagged site"
assert record.alignments[6].length==272
assert record.alignments[7].title==">gi|6121596|gb|G56277.1|G56277 SHGC-101791 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[7].length==641
assert record.alignments[8].title==">gi|5222417|gb|G51240.1|G51240 SHGC-80720 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[8].length==712
assert record.alignments[9].title==">gi|5221977|gb|G50800.1|G50800 SHGC-83850 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[9].length==422
assert record.alignments[10].title==">gi|5224501|gb|G53324.1|G53324 SHGC-82315 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[10].length==428
assert record.alignments[11].title==">gi|4529247|gb|G48587.1|G48587 SHGC-82546 Human Homo sapiens STS genomic, sequence tagged site"
assert record.alignments[11].length==694
assert record.alignments[12].title==">gi|3359917|gb|G40708|G40708 Z8947 Zebrafish AB Danio rerio STS genomic"
assert record.alignments[12].length==549
assert record.alignments[13].title==">gi|3359244|gb|G40035|G40035 Z13538 Zebrafish AB Danio rerio STS genomic"
assert record.alignments[13].length==536
assert record.alignments[14].title==">gi|1347715|gb|G25483|G25483 human STS EST334642."
assert record.alignments[14].length==407
assert record.alignments[15].title==">gi|1244262|gb|G19475|G19475 human STS SHGC-18755."
assert record.alignments[15].length==400
assert record.alignments[16].title==">gi|1232611|emb|Z51311|HS302WC9 H.sapiens (D5S2069) DNA segment containing (CA) repeat; clone AFM302wc9; single read, sequence tagged site [Homo sapiens]"
assert record.alignments[16].length==374
assert record.alignments[17].title==">gi|1223022|gb|G18565|G18565 BMS485 cow Bos taurus STS genomic, sequence tagged site [Bos taurus]"
assert record.alignments[17].length==181
assert record.alignments[18].title==">gi|1161779|gb|G15890|G15890 human STS CHLC.UTR_01448_M84721.P56085 clone UTR_01448_M84721."
assert record.alignments[18].length==729
assert record.alignments[19].title==">gi|858803|gb|G05558|G05558 human STS WI-7105."
assert record.alignments[19].length==735
assert record.alignments[20].title==">gi|1342455|gb|G22129|G22129 human STS WI-14200."
assert record.alignments[20].length==373
assert record.alignments[21].title==">gi|1347001|gb|G24769|G24769 human STS EST129834."
assert record.alignments[21].length==306
assert record.alignments[22].title==">gi|605469|gb|L31223|HUMUT821B Human STS UT821, 3' primer bind."
assert record.alignments[22].length==127
assert record.alignments[0].hsps[0].score==261
assert record.alignments[0].hsps[0].bits==517
assert record.alignments[0].hsps[0].expect==1e-146
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==16
assert record.alignments[1].hsps[0].bits==32.2
assert record.alignments[1].hsps[0].expect==1.6
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==16
assert record.alignments[2].hsps[0].bits==32.2
assert record.alignments[2].hsps[0].expect==1.6
assert len(record.alignments[2].hsps)==1
assert record.alignments[3].hsps[0].score==16
assert record.alignments[3].hsps[0].bits==32.2
assert record.alignments[3].hsps[0].expect==1.6
assert len(record.alignments[3].hsps)==1
assert record.alignments[4].hsps[0].score==15
assert record.alignments[4].hsps[0].bits==30.2
assert record.alignments[4].hsps[0].expect==6.5
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==15
assert record.alignments[5].hsps[0].bits==30.2
assert record.alignments[5].hsps[0].expect==6.5
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==15
assert record.alignments[6].hsps[0].bits==30.2
assert record.alignments[6].hsps[0].expect==6.5
assert len(record.alignments[6].hsps)==1
assert record.alignments[7].hsps[0].score==15
assert record.alignments[7].hsps[0].bits==30.2
assert record.alignments[7].hsps[0].expect==6.5
assert len(record.alignments[7].hsps)==1
assert record.alignments[8].hsps[0].score==15
assert record.alignments[8].hsps[0].bits==30.2
assert record.alignments[8].hsps[0].expect==6.5
assert len(record.alignments[8].hsps)==1
assert record.alignments[9].hsps[0].score==15
assert record.alignments[9].hsps[0].bits==30.2
assert record.alignments[9].hsps[0].expect==6.5
assert len(record.alignments[9].hsps)==1
assert record.alignments[10].hsps[0].score==15
assert record.alignments[10].hsps[0].bits==30.2
assert record.alignments[10].hsps[0].expect==6.5
assert len(record.alignments[10].hsps)==1
assert record.alignments[11].hsps[0].score==15
assert record.alignments[11].hsps[0].bits==30.2
assert record.alignments[11].hsps[0].expect==6.5
assert len(record.alignments[11].hsps)==1
assert record.alignments[12].hsps[0].score==15
assert record.alignments[12].hsps[0].bits==30.2
assert record.alignments[12].hsps[0].expect==6.5
assert len(record.alignments[12].hsps)==1
assert record.alignments[13].hsps[0].score==15
assert record.alignments[13].hsps[0].bits==30.2
assert record.alignments[13].hsps[0].expect==6.5
assert len(record.alignments[13].hsps)==1
assert record.alignments[14].hsps[0].score==15
assert record.alignments[14].hsps[0].bits==30.2
assert record.alignments[14].hsps[0].expect==6.5
assert len(record.alignments[14].hsps)==1
assert record.alignments[15].hsps[0].score==15
assert record.alignments[15].hsps[0].bits==30.2
assert record.alignments[15].hsps[0].expect==6.5
assert len(record.alignments[15].hsps)==1
assert record.alignments[16].hsps[0].score==15
assert record.alignments[16].hsps[0].bits==30.2
assert record.alignments[16].hsps[0].expect==6.5
assert len(record.alignments[16].hsps)==1
assert record.alignments[17].hsps[0].score==15
assert record.alignments[17].hsps[0].bits==30.2
assert record.alignments[17].hsps[0].expect==6.5
assert len(record.alignments[17].hsps)==1
assert record.alignments[18].hsps[0].score==15
assert record.alignments[18].hsps[0].bits==30.2
assert record.alignments[18].hsps[0].expect==6.5
assert len(record.alignments[18].hsps)==1
assert record.alignments[19].hsps[0].score==15
assert record.alignments[19].hsps[0].bits==30.2
assert record.alignments[19].hsps[0].expect==6.5
assert len(record.alignments[19].hsps)==1
assert record.alignments[20].hsps[0].score==15
assert record.alignments[20].hsps[0].bits==30.2
assert record.alignments[20].hsps[0].expect==6.5
assert len(record.alignments[20].hsps)==1
assert record.alignments[21].hsps[0].score==15
assert record.alignments[21].hsps[0].bits==30.2
assert record.alignments[21].hsps[0].expect==6.5
assert len(record.alignments[21].hsps)==1
assert record.alignments[22].hsps[0].score==15
assert record.alignments[22].hsps[0].bits==30.2
assert record.alignments[22].hsps[0].expect==6.5
assert len(record.alignments[22].hsps)==1
assert record.alignments[0].hsps[0].identities==(285, 285)
assert record.alignments[1].hsps[0].identities==(16, 16)
assert record.alignments[2].hsps[0].identities==(18, 19)
assert record.alignments[3].hsps[0].identities==(16, 16)
assert record.alignments[4].hsps[0].identities==(15, 15)
assert record.alignments[5].hsps[0].identities==(15, 15)
assert record.alignments[6].hsps[0].identities==(15, 15)
assert record.alignments[7].hsps[0].identities==(15, 15)
assert record.alignments[8].hsps[0].identities==(17, 18)
assert record.alignments[9].hsps[0].identities==(15, 15)
assert record.alignments[10].hsps[0].identities==(15, 15)
assert record.alignments[11].hsps[0].identities==(17, 18)
assert record.alignments[12].hsps[0].identities==(15, 15)
assert record.alignments[13].hsps[0].identities==(15, 15)
assert record.alignments[14].hsps[0].identities==(15, 15)
assert record.alignments[15].hsps[0].identities==(15, 15)
assert record.alignments[16].hsps[0].identities==(15, 15)
assert record.alignments[17].hsps[0].identities==(15, 15)
assert record.alignments[18].hsps[0].identities==(18, 19)
assert record.alignments[19].hsps[0].identities==(15, 15)
assert record.alignments[20].hsps[0].identities==(15, 15)
assert record.alignments[21].hsps[0].identities==(15, 15)
assert record.alignments[22].hsps[0].identities==(15, 15)
assert record.alignments[0].hsps[0].strand==("Plus", "Plus")
assert record.alignments[1].hsps[0].strand==("Plus", "Plus")
assert record.alignments[2].hsps[0].strand==("Plus", "Plus")
assert record.alignments[3].hsps[0].strand==("Plus", "Minus")
assert record.alignments[4].hsps[0].strand==("Plus", "Plus")
assert record.alignments[5].hsps[0].strand==("Plus", "Minus")
assert record.alignments[6].hsps[0].strand==("Plus", "Minus")
assert record.alignments[7].hsps[0].strand==("Plus", "Minus")
assert record.alignments[8].hsps[0].strand==("Plus", "Plus")
assert record.alignments[9].hsps[0].strand==("Plus", "Minus")
assert record.alignments[10].hsps[0].strand==("Plus", "Minus")
assert record.alignments[11].hsps[0].strand==("Plus", "Plus")
assert record.alignments[12].hsps[0].strand==("Plus", "Minus")
assert record.alignments[13].hsps[0].strand==("Plus", "Minus")
assert record.alignments[14].hsps[0].strand==("Plus", "Plus")
assert record.alignments[15].hsps[0].strand==("Plus", "Minus")
assert record.alignments[16].hsps[0].strand==("Plus", "Plus")
assert record.alignments[17].hsps[0].strand==("Plus", "Plus")
assert record.alignments[18].hsps[0].strand==("Plus", "Plus")
assert record.alignments[19].hsps[0].strand==("Plus", "Minus")
assert record.alignments[20].hsps[0].strand==("Plus", "Minus")
assert record.alignments[21].hsps[0].strand==("Plus", "Plus")
assert record.alignments[22].hsps[0].strand==("Plus", "Minus")
assert record.alignments[0].hsps[0].query=="gatccctacccttnccgttggtctctntcgctgactcgaggcacctaacatccattcacacccaacacaggccagcgacttctggggctcagccacagacatggtttgtnactnttgagcttctgttcctagagaatcctagaggcttgattggcccaggctgctgtntgtnctggaggcaaagaatccctacctcctaggggtgaaaggaaatnaaaatggaaagttcttgtagcgcaaggcctgacatgggtagctgctcaataaatgctagtntgttatttc"
assert record.alignments[0].hsps[0].match=="|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
assert record.alignments[0].hsps[0].sbjct=="gatccctacccttnccgttggtctctntcgctgactcgaggcacctaacatccattcacacccaacacaggccagcgacttctggggctcagccacagacatggtttgtnactnttgagcttctgttcctagagaatcctagaggcttgattggcccaggctgctgtntgtnctggaggcaaagaatccctacctcctaggggtgaaaggaaatnaaaatggaaagttcttgtagcgcaaggcctgacatgggtagctgctcaataaatgctagtntgttatttc"
assert record.alignments[0].hsps[0].query_start==1
assert record.alignments[0].hsps[0].query_end==285
assert record.alignments[0].hsps[0].sbjct_start==1
assert record.alignments[0].hsps[0].sbjct_end==285
assert record.alignments[1].hsps[0].query=="ggaaagttcttgtagc"
assert record.alignments[1].hsps[0].match=="||||||||||||||||"
assert record.alignments[1].hsps[0].sbjct=="ggaaagttcttgtagc"
assert record.alignments[1].hsps[0].query_start==221
assert record.alignments[1].hsps[0].query_end==236
assert record.alignments[1].hsps[0].sbjct_start==32
assert record.alignments[1].hsps[0].sbjct_end==47
assert record.alignments[2].hsps[0].query=="gaaatnaaaatggaaagtt"
assert record.alignments[2].hsps[0].match=="||||| |||||||||||||"
assert record.alignments[2].hsps[0].sbjct=="gaaataaaaatggaaagtt"
assert record.alignments[2].hsps[0].query_start==210
assert record.alignments[2].hsps[0].query_end==228
assert record.alignments[2].hsps[0].sbjct_start==588
assert record.alignments[2].hsps[0].sbjct_end==606
assert record.alignments[3].hsps[0].query=="ctcaataaatgctagt"
assert record.alignments[3].hsps[0].match=="||||||||||||||||"
assert record.alignments[3].hsps[0].sbjct=="ctcaataaatgctagt"
assert record.alignments[3].hsps[0].query_start==260
assert record.alignments[3].hsps[0].query_end==275
assert record.alignments[3].hsps[0].sbjct_start==178
assert record.alignments[3].hsps[0].sbjct_end==163
assert record.alignments[4].hsps[0].query=="aaaatggaaagttct"
assert record.alignments[4].hsps[0].match=="|||||||||||||||"
assert record.alignments[4].hsps[0].sbjct=="aaaatggaaagttct"
assert record.alignments[4].hsps[0].query_start==216
assert record.alignments[4].hsps[0].query_end==230
assert record.alignments[4].hsps[0].sbjct_start==330
assert record.alignments[4].hsps[0].sbjct_end==344
assert record.alignments[5].hsps[0].query=="ttctgttcctagaga"
assert record.alignments[5].hsps[0].match=="|||||||||||||||"
assert record.alignments[5].hsps[0].sbjct=="ttctgttcctagaga"
assert record.alignments[5].hsps[0].query_start==121
assert record.alignments[5].hsps[0].query_end==135
assert record.alignments[5].hsps[0].sbjct_start==384
assert record.alignments[5].hsps[0].sbjct_end==370
assert record.alignments[6].hsps[0].query=="ggaaagttcttgtag"
assert record.alignments[6].hsps[0].match=="|||||||||||||||"
assert record.alignments[6].hsps[0].sbjct=="ggaaagttcttgtag"
assert record.alignments[6].hsps[0].query_start==221
assert record.alignments[6].hsps[0].query_end==235
assert record.alignments[6].hsps[0].sbjct_start==138
assert record.alignments[6].hsps[0].sbjct_end==124
assert record.alignments[7].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[7].hsps[0].match=="|||||||||||||||"
assert record.alignments[7].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[7].hsps[0].query_start==258
assert record.alignments[7].hsps[0].query_end==272
assert record.alignments[7].hsps[0].sbjct_start==216
assert record.alignments[7].hsps[0].sbjct_end==202
assert record.alignments[8].hsps[0].query=="taaatgctagtntgttat"
assert record.alignments[8].hsps[0].match=="||||||||||| ||||||"
assert record.alignments[8].hsps[0].sbjct=="taaatgctagtttgttat"
assert record.alignments[8].hsps[0].query_start==265
assert record.alignments[8].hsps[0].query_end==282
assert record.alignments[8].hsps[0].sbjct_start==293
assert record.alignments[8].hsps[0].sbjct_end==310
assert record.alignments[9].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[9].hsps[0].match=="|||||||||||||||"
assert record.alignments[9].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[9].hsps[0].query_start==258
assert record.alignments[9].hsps[0].query_end==272
assert record.alignments[9].hsps[0].sbjct_start==32
assert record.alignments[9].hsps[0].sbjct_end==18
assert record.alignments[10].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[10].hsps[0].match=="|||||||||||||||"
assert record.alignments[10].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[10].hsps[0].query_start==258
assert record.alignments[10].hsps[0].query_end==272
assert record.alignments[10].hsps[0].sbjct_start==47
assert record.alignments[10].hsps[0].sbjct_end==33
assert record.alignments[11].hsps[0].query=="taaatgctagtntgttat"
assert record.alignments[11].hsps[0].match=="||||||||||| ||||||"
assert record.alignments[11].hsps[0].sbjct=="taaatgctagtttgttat"
assert record.alignments[11].hsps[0].query_start==265
assert record.alignments[11].hsps[0].query_end==282
assert record.alignments[11].hsps[0].sbjct_start==292
assert record.alignments[11].hsps[0].sbjct_end==309
assert record.alignments[12].hsps[0].query=="aacatccattcacac"
assert record.alignments[12].hsps[0].match=="|||||||||||||||"
assert record.alignments[12].hsps[0].sbjct=="aacatccattcacac"
assert record.alignments[12].hsps[0].query_start==47
assert record.alignments[12].hsps[0].query_end==61
assert record.alignments[12].hsps[0].sbjct_start==479
assert record.alignments[12].hsps[0].sbjct_end==465
assert record.alignments[13].hsps[0].query=="ttctgttcctagaga"
assert record.alignments[13].hsps[0].match=="|||||||||||||||"
assert record.alignments[13].hsps[0].sbjct=="ttctgttcctagaga"
assert record.alignments[13].hsps[0].query_start==121
assert record.alignments[13].hsps[0].query_end==135
assert record.alignments[13].hsps[0].sbjct_start==433
assert record.alignments[13].hsps[0].sbjct_end==419
assert record.alignments[14].hsps[0].query=="ctaacatccattcac"
assert record.alignments[14].hsps[0].match=="|||||||||||||||"
assert record.alignments[14].hsps[0].sbjct=="ctaacatccattcac"
assert record.alignments[14].hsps[0].query_start==45
assert record.alignments[14].hsps[0].query_end==59
assert record.alignments[14].hsps[0].sbjct_start==389
assert record.alignments[14].hsps[0].sbjct_end==403
assert record.alignments[15].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[15].hsps[0].match=="|||||||||||||||"
assert record.alignments[15].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[15].hsps[0].query_start==258
assert record.alignments[15].hsps[0].query_end==272
assert record.alignments[15].hsps[0].sbjct_start==324
assert record.alignments[15].hsps[0].sbjct_end==310
assert record.alignments[16].hsps[0].query=="acagacatggtttgt"
assert record.alignments[16].hsps[0].match=="|||||||||||||||"
assert record.alignments[16].hsps[0].sbjct=="acagacatggtttgt"
assert record.alignments[16].hsps[0].query_start==95
assert record.alignments[16].hsps[0].query_end==109
assert record.alignments[16].hsps[0].sbjct_start==246
assert record.alignments[16].hsps[0].sbjct_end==260
assert record.alignments[17].hsps[0].query=="ctcaataaatgctag"
assert record.alignments[17].hsps[0].match=="|||||||||||||||"
assert record.alignments[17].hsps[0].sbjct=="ctcaataaatgctag"
assert record.alignments[17].hsps[0].query_start==260
assert record.alignments[17].hsps[0].query_end==274
assert record.alignments[17].hsps[0].sbjct_start==145
assert record.alignments[17].hsps[0].sbjct_end==159
assert record.alignments[18].hsps[0].query=="gtagctgctcaataaatgc"
assert record.alignments[18].hsps[0].match=="|||| ||||||||||||||"
assert record.alignments[18].hsps[0].sbjct=="gtaggtgctcaataaatgc"
assert record.alignments[18].hsps[0].query_start==253
assert record.alignments[18].hsps[0].query_end==271
assert record.alignments[18].hsps[0].sbjct_start==698
assert record.alignments[18].hsps[0].sbjct_end==716
assert record.alignments[19].hsps[0].query=="gaaagttcttgtagc"
assert record.alignments[19].hsps[0].match=="|||||||||||||||"
assert record.alignments[19].hsps[0].sbjct=="gaaagttcttgtagc"
assert record.alignments[19].hsps[0].query_start==222
assert record.alignments[19].hsps[0].query_end==236
assert record.alignments[19].hsps[0].sbjct_start==543
assert record.alignments[19].hsps[0].sbjct_end==529
assert record.alignments[20].hsps[0].query=="tgctcaataaatgct"
assert record.alignments[20].hsps[0].match=="|||||||||||||||"
assert record.alignments[20].hsps[0].sbjct=="tgctcaataaatgct"
assert record.alignments[20].hsps[0].query_start==258
assert record.alignments[20].hsps[0].query_end==272
assert record.alignments[20].hsps[0].sbjct_start==33
assert record.alignments[20].hsps[0].sbjct_end==19
assert record.alignments[21].hsps[0].query=="tggaaagttcttgta"
assert record.alignments[21].hsps[0].match=="|||||||||||||||"
assert record.alignments[21].hsps[0].sbjct=="tggaaagttcttgta"
assert record.alignments[21].hsps[0].query_start==220
assert record.alignments[21].hsps[0].query_end==234
assert record.alignments[21].hsps[0].sbjct_start==144
assert record.alignments[21].hsps[0].sbjct_end==158
assert record.alignments[22].hsps[0].query=="acagacatggtttgt"
assert record.alignments[22].hsps[0].match=="|||||||||||||||"
assert record.alignments[22].hsps[0].sbjct=="acagacatggtttgt"
assert record.alignments[22].hsps[0].query_start==95
assert record.alignments[22].hsps[0].query_end==109
assert record.alignments[22].hsps[0].sbjct_start==106
assert record.alignments[22].hsps[0].sbjct_end==92
assert record.database_name==['data/sts']
assert record.num_letters_in_database==[31998854]
assert record.num_sequences_in_database==[87792]
assert record.posted_date==[('Feb 11, 2000  2:37 PM',)]
assert record.ka_params==[1.37, 0.711, 1.31]
assert record.ka_params_gap==[1.37, 0.711, 1.31]
assert record.matrix=='blastn matrix:1 -3'
assert record.gap_penalties==[5,2]
assert record.num_hits==3835
assert record.num_sequences==87792
assert record.num_extends==3835
assert record.num_good_extends==301
assert record.num_seqs_better_e==24
assert record.query_length==285
assert record.database_length==31998854
assert record.effective_hsp_length==17
assert record.effective_query_length==268
assert record.effective_database_length==30506390
assert record.effective_search_space==8175712520
assert record.effective_search_space_used==8175712520
assert record.threshold==0
assert record.window_size==0
assert record.dropoff_1st_pass==(6,11.9)
assert record.gap_x_dropoff==(10,19.8)
assert record.gap_trigger==(12,24.3)
assert record.blast_cutoff==(15,30.2)

handle = open('Blast/bt063')
record = pb_parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.0.14'
assert record.date=="Jun-29-2000"
assert record.reference==reference
assert record.query=="aseq"
assert record.query_letters==82
assert record.database=="SWISS_PROT"
assert record.database_sequences==87404
assert record.database_letters==31753758
assert len(record.rounds)==2
assert len(record.rounds[0].new_seqs)==6
assert record.rounds[0].new_seqs[0].title=="YMZ9_YEAST Q04439 saccharomyces cerevisiae (baker's yeast). ..."
assert record.rounds[0].new_seqs[0].score==27
assert record.rounds[0].new_seqs[0].e==6.000000
assert record.rounds[0].new_seqs[1].title=="DPOL_WHVW6 P11292 woodchuck hepatitis virus w64 (isolate pws..."
assert record.rounds[0].new_seqs[1].score==27
assert record.rounds[0].new_seqs[1].e==7.900000
assert record.rounds[0].new_seqs[2].title=="DPOL_WHV8I P17396 woodchuck hepatitis virus 8 (infectious cl..."
assert record.rounds[0].new_seqs[2].score==27
assert record.rounds[0].new_seqs[2].e==7.900000
assert record.rounds[0].new_seqs[3].title=="DPOL_WHV8 P06275 woodchuck hepatitis virus 8 (whv 8). dna po..."
assert record.rounds[0].new_seqs[3].score==27
assert record.rounds[0].new_seqs[3].e==7.900000
assert record.rounds[0].new_seqs[4].title=="DPOL_WHV7 P12898 woodchuck hepatitis virus 7 (whv 7). dna po..."
assert record.rounds[0].new_seqs[4].score==27
assert record.rounds[0].new_seqs[4].e==7.900000
assert record.rounds[0].new_seqs[5].title=="DPOL_WHV59 P12899 woodchuck hepatitis virus 59 (whv 59). dna..."
assert record.rounds[0].new_seqs[5].score==27
assert record.rounds[0].new_seqs[5].e==7.900000
assert len(record.rounds[0].alignments)==6
assert record.rounds[0].alignments[0].title==">YMZ9_YEAST Q04439 saccharomyces cerevisiae (baker's yeast). hypothetical myosin-like protein in ilv2-ade17 intergenic region. 5/2000"
assert record.rounds[0].alignments[0].length==1219
assert record.rounds[0].alignments[1].title==">DPOL_WHVW6 P11292 woodchuck hepatitis virus w64 (isolate pws23). dna polymerase (ec 2.7.7.7) (fragment). 12/1998"
assert record.rounds[0].alignments[1].length==556
assert record.rounds[0].alignments[2].title==">DPOL_WHV8I P17396 woodchuck hepatitis virus 8 (infectious clone) (whv 8). dna polymerase (ec 2.7.7.7). 12/1998"
assert record.rounds[0].alignments[2].length==884
assert record.rounds[0].alignments[3].title==">DPOL_WHV8 P06275 woodchuck hepatitis virus 8 (whv 8). dna polymerase (ec 2.7.7.7). 12/1998"
assert record.rounds[0].alignments[3].length==883
assert record.rounds[0].alignments[4].title==">DPOL_WHV7 P12898 woodchuck hepatitis virus 7 (whv 7). dna polymerase (ec 2.7.7.7). 12/1998"
assert record.rounds[0].alignments[4].length==884
assert record.rounds[0].alignments[5].title==">DPOL_WHV59 P12899 woodchuck hepatitis virus 59 (whv 59). dna polymerase (ec 2.7.7.7). 12/1998"
assert record.rounds[0].alignments[5].length==884
assert len(record.rounds[1].new_seqs)==7
assert record.rounds[1].new_seqs[0].title=="YMZ9_YEAST Q04439 saccharomyces cerevisiae (baker's yeast). ..."
assert record.rounds[1].new_seqs[0].score==27
assert record.rounds[1].new_seqs[0].e==5.9
assert record.rounds[1].new_seqs[1].title=="DPOL_WHVW6 P11292 woodchuck hepatitis virus w64 (isolate pws..."
assert record.rounds[1].new_seqs[1].score==27
assert record.rounds[1].new_seqs[1].e==5.9
assert record.rounds[1].new_seqs[2].title=="DPOL_WHV8I P17396 woodchuck hepatitis virus 8 (infectious cl..."
assert record.rounds[1].new_seqs[2].score==27
assert record.rounds[1].new_seqs[2].e==5.9
assert record.rounds[1].new_seqs[3].title=="DPOL_WHV8 P06275 woodchuck hepatitis virus 8 (whv 8). dna po..."
assert record.rounds[1].new_seqs[3].score==27
assert record.rounds[1].new_seqs[3].e==5.9
assert record.rounds[1].new_seqs[4].title=="DPOL_WHV7 P12898 woodchuck hepatitis virus 7 (whv 7). dna po..."
assert record.rounds[1].new_seqs[4].score==27
assert record.rounds[1].new_seqs[4].e==5.9
assert record.rounds[1].new_seqs[5].title=="DPOL_WHV59 P12899 woodchuck hepatitis virus 59 (whv 59). dna..."
assert record.rounds[1].new_seqs[5].score==27
assert record.rounds[1].new_seqs[5].e==5.9
assert record.rounds[1].new_seqs[6].title=="DPOL_HPBGS P03161 ground squirrel hepatitis virus (gshv). dn..."
assert record.rounds[1].new_seqs[6].score==27
assert record.rounds[1].new_seqs[6].e==7.8
assert len(record.rounds[1].alignments)==7
assert record.rounds[1].alignments[0].title==">YMZ9_YEAST Q04439 saccharomyces cerevisiae (baker's yeast). hypothetical myosin-like protein in ilv2-ade17 intergenic region. 5/2000"
assert record.rounds[1].alignments[0].length==1219
assert record.rounds[1].alignments[1].title==">DPOL_WHVW6 P11292 woodchuck hepatitis virus w64 (isolate pws23). dna polymerase (ec 2.7.7.7) (fragment). 12/1998"
assert record.rounds[1].alignments[1].length==556
assert record.rounds[1].alignments[2].title==">DPOL_WHV8I P17396 woodchuck hepatitis virus 8 (infectious clone) (whv 8). dna polymerase (ec 2.7.7.7). 12/1998"
assert record.rounds[1].alignments[2].length==884
assert record.rounds[1].alignments[3].title==">DPOL_WHV8 P06275 woodchuck hepatitis virus 8 (whv 8). dna polymerase (ec 2.7.7.7). 12/1998"
assert record.rounds[1].alignments[3].length==883
assert record.rounds[1].alignments[4].title==">DPOL_WHV7 P12898 woodchuck hepatitis virus 7 (whv 7). dna polymerase (ec 2.7.7.7). 12/1998"
assert record.rounds[1].alignments[4].length==884
assert record.rounds[1].alignments[5].title==">DPOL_WHV59 P12899 woodchuck hepatitis virus 59 (whv 59). dna polymerase (ec 2.7.7.7). 12/1998"
assert record.rounds[1].alignments[5].length==884
assert record.rounds[1].alignments[6].title==">DPOL_HPBGS P03161 ground squirrel hepatitis virus (gshv). dna polymerase (ec 2.7.7.7) (a protein). 12/1998"
assert record.rounds[1].alignments[6].length==881
assert record.rounds[0].alignments[0].hsps[0].score==59
assert record.rounds[0].alignments[0].hsps[0].bits==27.4
assert record.rounds[0].alignments[0].hsps[0].expect==6.0
assert len(record.rounds[0].alignments[0].hsps)==1
assert record.rounds[0].alignments[1].hsps[0].score==58
assert record.rounds[0].alignments[1].hsps[0].bits==27.0
assert record.rounds[0].alignments[1].hsps[0].expect==7.9
assert len(record.rounds[0].alignments[1].hsps)==1
assert record.rounds[0].alignments[2].hsps[0].score==58
assert record.rounds[0].alignments[2].hsps[0].bits==27.0
assert record.rounds[0].alignments[2].hsps[0].expect==7.9
assert len(record.rounds[0].alignments[2].hsps)==1
assert record.rounds[0].alignments[3].hsps[0].score==58
assert record.rounds[0].alignments[3].hsps[0].bits==27.0
assert record.rounds[0].alignments[3].hsps[0].expect==7.9
assert len(record.rounds[0].alignments[3].hsps)==1
assert record.rounds[0].alignments[4].hsps[0].score==58
assert record.rounds[0].alignments[4].hsps[0].bits==27.0
assert record.rounds[0].alignments[4].hsps[0].expect==7.9
assert len(record.rounds[0].alignments[4].hsps)==1
assert record.rounds[0].alignments[5].hsps[0].score==58
assert record.rounds[0].alignments[5].hsps[0].bits==27.0
assert record.rounds[0].alignments[5].hsps[0].expect==7.9
assert record.rounds[1].alignments[0].hsps[0].score==59
assert record.rounds[1].alignments[0].hsps[0].bits==27.4
assert record.rounds[1].alignments[0].hsps[0].expect==5.9
assert len(record.rounds[1].alignments[0].hsps)==1
assert record.rounds[1].alignments[1].hsps[0].score==59
assert record.rounds[1].alignments[1].hsps[0].bits==27.4
assert record.rounds[1].alignments[1].hsps[0].expect==5.9
assert len(record.rounds[1].alignments[1].hsps)==1
assert record.rounds[1].alignments[2].hsps[0].score==59
assert record.rounds[1].alignments[2].hsps[0].bits==27.4
assert record.rounds[1].alignments[2].hsps[0].expect==5.9
assert len(record.rounds[1].alignments[2].hsps)==1
assert record.rounds[1].alignments[3].hsps[0].score==59
assert record.rounds[1].alignments[3].hsps[0].bits==27.4
assert record.rounds[1].alignments[3].hsps[0].expect==5.9
assert len(record.rounds[1].alignments[3].hsps)==1
assert record.rounds[1].alignments[4].hsps[0].score==59
assert record.rounds[1].alignments[4].hsps[0].bits==27.4
assert record.rounds[1].alignments[4].hsps[0].expect==5.9
assert len(record.rounds[1].alignments[4].hsps)==1
assert record.rounds[1].alignments[5].hsps[0].score==59
assert record.rounds[1].alignments[5].hsps[0].bits==27.4
assert record.rounds[1].alignments[5].hsps[0].expect==5.9
assert len(record.rounds[1].alignments[5].hsps)==1
assert record.rounds[1].alignments[6].hsps[0].score==58
assert record.rounds[1].alignments[6].hsps[0].bits==27.0
assert record.rounds[1].alignments[6].hsps[0].expect==7.8
assert len(record.rounds[1].alignments[6].hsps)==1
assert record.rounds[0].alignments[0].hsps[0].identities==(10, 22)
assert record.rounds[0].alignments[0].hsps[0].positives==(16, 22)
assert record.rounds[0].alignments[1].hsps[0].identities==(17, 43)
assert record.rounds[0].alignments[1].hsps[0].positives==(20, 43)
assert record.rounds[0].alignments[1].hsps[0].gaps==(10, 43)
assert record.rounds[0].alignments[2].hsps[0].identities==(17, 43)
assert record.rounds[0].alignments[2].hsps[0].positives==(20, 43)
assert record.rounds[0].alignments[2].hsps[0].gaps==(10, 43)
assert record.rounds[0].alignments[3].hsps[0].identities==(17, 43)
assert record.rounds[0].alignments[3].hsps[0].positives==(20, 43)
assert record.rounds[0].alignments[3].hsps[0].gaps==(10, 43)
assert record.rounds[0].alignments[4].hsps[0].identities==(17, 43)
assert record.rounds[0].alignments[4].hsps[0].positives==(20, 43)
assert record.rounds[0].alignments[4].hsps[0].gaps==(10, 43)
assert record.rounds[0].alignments[5].hsps[0].identities==(17, 43)
assert record.rounds[0].alignments[5].hsps[0].positives==(20, 43)
assert record.rounds[0].alignments[5].hsps[0].gaps==(10, 43)
assert record.rounds[1].alignments[0].hsps[0].identities==(10, 22)
assert record.rounds[1].alignments[0].hsps[0].positives==(16, 22)
assert record.rounds[1].alignments[1].hsps[0].identities==(17, 43)
assert record.rounds[1].alignments[1].hsps[0].positives==(20, 43)
assert record.rounds[1].alignments[1].hsps[0].gaps==(10, 43)
assert record.rounds[1].alignments[2].hsps[0].identities==(17, 43)
assert record.rounds[1].alignments[2].hsps[0].positives==(20, 43)
assert record.rounds[1].alignments[2].hsps[0].gaps==(10, 43)
assert record.rounds[1].alignments[3].hsps[0].identities==(17, 43)
assert record.rounds[1].alignments[3].hsps[0].positives==(20, 43)
assert record.rounds[1].alignments[3].hsps[0].gaps==(10, 43)
assert record.rounds[1].alignments[4].hsps[0].identities==(17, 43)
assert record.rounds[1].alignments[4].hsps[0].positives==(20, 43)
assert record.rounds[1].alignments[4].hsps[0].gaps==(10, 43)
assert record.rounds[1].alignments[5].hsps[0].identities==(17, 43)
assert record.rounds[1].alignments[5].hsps[0].positives==(20, 43)
assert record.rounds[1].alignments[5].hsps[0].gaps==(10, 43)
assert record.rounds[1].alignments[6].hsps[0].identities==(17, 43)
assert record.rounds[1].alignments[6].hsps[0].positives==(20, 43)
assert record.rounds[1].alignments[6].hsps[0].gaps==(10, 43)
assert record.rounds[0].alignments[0].hsps[0].query=="NLIALQHIPLSPAGVIAKRPAP"
assert record.rounds[0].alignments[0].hsps[0].match=="++ A QH+P +PA   +K+PAP"
assert record.rounds[0].alignments[0].hsps[0].sbjct=="SIAAAQHVPTAPASRHSKKPAP"
assert record.rounds[0].alignments[0].hsps[0].query_start==52
assert record.rounds[0].alignments[0].hsps[0].query_end==73
assert record.rounds[0].alignments[0].hsps[0].sbjct_start==992
assert record.rounds[0].alignments[0].hsps[0].sbjct_end==1013
assert record.rounds[0].alignments[1].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[0].alignments[1].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[0].alignments[1].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[0].alignments[1].hsps[0].query_start==34
assert record.rounds[0].alignments[1].hsps[0].query_end==66
assert record.rounds[0].alignments[1].hsps[0].sbjct_start==115
assert record.rounds[0].alignments[1].hsps[0].sbjct_end==157
assert record.rounds[0].alignments[2].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[0].alignments[2].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[0].alignments[2].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[0].alignments[2].hsps[0].query_start==34
assert record.rounds[0].alignments[2].hsps[0].query_end==66
assert record.rounds[0].alignments[2].hsps[0].sbjct_start==443
assert record.rounds[0].alignments[2].hsps[0].sbjct_end==485
assert record.rounds[0].alignments[3].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[0].alignments[3].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[0].alignments[3].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[0].alignments[3].hsps[0].query_start==34
assert record.rounds[0].alignments[3].hsps[0].query_end==66
assert record.rounds[0].alignments[3].hsps[0].sbjct_start==442
assert record.rounds[0].alignments[3].hsps[0].sbjct_end==484
assert record.rounds[0].alignments[4].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[0].alignments[4].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[0].alignments[4].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[0].alignments[4].hsps[0].query_start==34
assert record.rounds[0].alignments[4].hsps[0].query_end==66
assert record.rounds[0].alignments[4].hsps[0].sbjct_start==443
assert record.rounds[0].alignments[4].hsps[0].sbjct_end==485
assert record.rounds[0].alignments[5].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[0].alignments[5].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[0].alignments[5].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[0].alignments[5].hsps[0].query_start==34
assert record.rounds[0].alignments[5].hsps[0].query_end==66
assert record.rounds[0].alignments[5].hsps[0].sbjct_start==443
assert record.rounds[0].alignments[5].hsps[0].sbjct_end==485
assert record.rounds[1].alignments[0].hsps[0].query=="NLIALQHIPLSPAGVIAKRPAP"
assert record.rounds[1].alignments[0].hsps[0].match=="++ A QH+P +PA   +K+PAP"
assert record.rounds[1].alignments[0].hsps[0].sbjct=="SIAAAQHVPTAPASRHSKKPAP"
assert record.rounds[1].alignments[0].hsps[0].query_start==52
assert record.rounds[1].alignments[0].hsps[0].query_end==73
assert record.rounds[1].alignments[0].hsps[0].sbjct_start==992
assert record.rounds[1].alignments[0].hsps[0].sbjct_end==1013
assert record.rounds[1].alignments[1].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[1].alignments[1].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[1].alignments[1].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[1].alignments[1].hsps[0].query_start==34
assert record.rounds[1].alignments[1].hsps[0].query_end==66
assert record.rounds[1].alignments[1].hsps[0].sbjct_start==115
assert record.rounds[1].alignments[1].hsps[0].sbjct_end==157
assert record.rounds[1].alignments[2].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[1].alignments[2].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[1].alignments[2].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[1].alignments[2].hsps[0].query_start==34
assert record.rounds[1].alignments[2].hsps[0].query_end==66
assert record.rounds[1].alignments[2].hsps[0].sbjct_start==443
assert record.rounds[1].alignments[2].hsps[0].sbjct_end==485
assert record.rounds[1].alignments[3].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[1].alignments[3].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[1].alignments[3].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[1].alignments[3].hsps[0].query_start==34
assert record.rounds[1].alignments[3].hsps[0].query_end==66
assert record.rounds[1].alignments[3].hsps[0].sbjct_start==442
assert record.rounds[1].alignments[3].hsps[0].sbjct_end==484
assert record.rounds[1].alignments[4].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[1].alignments[4].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[1].alignments[4].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[1].alignments[4].hsps[0].query_start==34
assert record.rounds[1].alignments[4].hsps[0].query_end==66
assert record.rounds[1].alignments[4].hsps[0].sbjct_start==443
assert record.rounds[1].alignments[4].hsps[0].sbjct_end==485
assert record.rounds[1].alignments[5].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[1].alignments[5].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[1].alignments[5].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPISPAAV"
assert record.rounds[1].alignments[5].hsps[0].query_start==34
assert record.rounds[1].alignments[5].hsps[0].query_end==66
assert record.rounds[1].alignments[5].hsps[0].sbjct_start==443
assert record.rounds[1].alignments[5].hsps[0].sbjct_end==485
assert record.rounds[1].alignments[6].hsps[0].query=="IHWPSFY--NVVTGKTLALPNL--------IALQHIPLSPAGV"
assert record.rounds[1].alignments[6].hsps[0].match=="+HWP F   N+ T   L   NL         A  HIP+SPA V"
assert record.rounds[1].alignments[6].hsps[0].sbjct=="VHWPKFAVPNLQTLANLLSTNLQWLSLDVSAAFYHIPVSPAAV"
assert record.rounds[1].alignments[6].hsps[0].query_start==34
assert record.rounds[1].alignments[6].hsps[0].query_end==66
assert record.rounds[1].alignments[6].hsps[0].sbjct_start==440
assert record.rounds[1].alignments[6].hsps[0].sbjct_end==482
assert record.database_name==['/dbase/swissprot/main/release/sp', '/dbase/swissprot/main/update/spu']
assert record.posted_date==[('Jun 21, 2000 12:39 PM',), ('Jul 30, 2000 11:00 PM',)]
assert record.num_letters_in_database==[31411157, 342601]
assert record.num_sequences_in_database==[86593, 811]
assert record.ka_params==[0.318, 0.132, 0.422]
assert record.ka_params_gap==[0.270, 0.0464, 0.230]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==7377465
assert record.num_sequences==87404
assert record.num_extends==229017
assert record.num_good_extends==569
assert record.num_seqs_better_e==13
assert record.hsps_no_gap==2
assert record.hsps_prelim_gapped==11
assert record.hsps_gapped==13
assert record.query_length==82
assert record.database_length==31753758
assert record.effective_hsp_length==44
assert record.effective_query_length==38
assert record.effective_database_length==27907982
assert record.effective_search_space==1060503316
assert record.effective_search_space_used==1060503316
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.8)
assert record.gap_x_dropoff_final==(64,24.9)
assert record.gap_trigger==(41,21.7)
assert record.blast_cutoff==(58,27.0)

handle = open('Blast/bt067')
record = parser.parse(handle)
assert record.application=="BLASTX"
assert record.version=='2.2.1'
assert record.date=="Jul-12-2001"
assert record.reference==reference
assert record.query=="gi|8333950|gb|BE038934.1|BE038934 AB07D09 AB Arabidopsis\nthaliana cDNA 5' similar to cold-regulated protein cor47, mRNA\nsequence"
assert record.query_letters==738
assert record.database=="cold_tolerance_protein.fasta"
assert record.database_sequences==19
assert record.database_letters==7312
assert len(record.descriptions)==2
assert record.descriptions[0].title=="gi|729204|sp|Q06396|CR13_ORYSA 13 KD COLD-INDUCED PROTEIN"
assert record.descriptions[0].score==21
assert record.descriptions[0].e==0.780000
assert record.descriptions[1].title=="gi|729205|sp|Q06397|CR18_ORYSA 18 KD COLD-INDUCED PROTEIN"
assert record.descriptions[1].score==18
assert record.descriptions[1].e==5.100000
assert len(record.alignments)==2
assert record.alignments[0].title==">gi|729204|sp|Q06396|CR13_ORYSA 13 KD COLD-INDUCED PROTEIN"
assert record.alignments[0].length==120
assert record.alignments[1].title==">gi|729205|sp|Q06397|CR18_ORYSA 18 KD COLD-INDUCED PROTEIN"
assert record.alignments[1].length==156
assert record.alignments[0].hsps[0].score==42
assert record.alignments[0].hsps[0].bits==20.8
assert record.alignments[0].hsps[0].expect==0.78
assert len(record.alignments[0].hsps)==1
assert record.alignments[1].hsps[0].score==35
assert record.alignments[1].hsps[0].bits==18.1
assert record.alignments[1].hsps[0].expect==5.1
assert len(record.alignments[1].hsps)==1
assert record.alignments[0].hsps[0].identities==(6, 7)
assert record.alignments[0].hsps[0].positives==(6, 7)
assert record.alignments[1].hsps[0].identities==(5, 7)
assert record.alignments[1].hsps[0].positives==(7, 7)
assert record.alignments[0].hsps[0].frame==("-3", )
assert record.alignments[1].hsps[0].frame==("+3", )
assert record.alignments[0].hsps[0].query=="CCRHYSF"
assert record.alignments[0].hsps[0].match=="CC HYSF"
assert record.alignments[0].hsps[0].sbjct=="CCSHYSF"
assert record.alignments[0].hsps[0].query_start==181
assert record.alignments[0].hsps[0].query_end==161
assert record.alignments[0].hsps[0].sbjct_start==93
assert record.alignments[0].hsps[0].sbjct_end==99
assert record.alignments[1].hsps[0].query=="HCLGVLS"
assert record.alignments[1].hsps[0].match=="HCLG++S"
assert record.alignments[1].hsps[0].sbjct=="HCLGIVS"
assert record.alignments[1].hsps[0].query_start==33
assert record.alignments[1].hsps[0].query_end==53
assert record.alignments[1].hsps[0].sbjct_start==54
assert record.alignments[1].hsps[0].sbjct_end==60
assert record.database_name==['cold_tolerance_protein.fasta']
assert record.num_letters_in_database==[7312]
assert record.num_sequences_in_database==[19]
assert record.posted_date==[('Sep 27, 2001  4:17 PM',)]
assert record.ka_params==[0.318, 0.135, 0.401]
assert record.ka_params_gap==[0.267, 0.0410, 0.140]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==3730
assert record.num_sequences==19
assert record.num_extends==61
assert record.num_good_extends==2
assert record.num_seqs_better_e==4
assert record.hsps_no_gap==2
assert record.hsps_prelim_gapped==0
assert record.hsps_gapped==2
assert record.database_length==7312
assert record.effective_hsp_length==47
assert record.effective_database_length==6419
assert record.effective_search_space_used==1270962
assert record.frameshift==('50,','0.1')
assert record.threshold==12
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.gap_x_dropoff==(38,14.6)
assert record.gap_x_dropoff_final==(64,24.7)
assert record.gap_trigger==(33,18.0)

handle = open('Blast/bt071')
record = pb_parser.parse(handle)
assert record.application=="BLASTP"
assert record.version=='2.2.8'
assert record.date=="Jan-05-2004"
assert record.reference==reference
assert record.query=="gi|16130963|ref|NP_417539.1| RNA polymerase sigma(70)\n[Escherichia coli K12]"
assert record.query_letters==613
assert record.database=="NC_000913"
assert record.database_sequences==4311
assert record.database_letters==1350842
assert len(record.rounds)==2
assert len(record.rounds[0].new_seqs)==4
assert record.rounds[0].new_seqs[0].title=="gi|16130963|ref|NP_417539.1| RNA polymerase sigma(70) [Escherich..."
assert record.rounds[0].new_seqs[0].score==1032
assert record.rounds[0].new_seqs[0].e==0.000000
assert record.rounds[0].new_seqs[1].title=="gi|16130648|ref|NP_417221.1| sigma S factor [Escherichia coli K12]"
assert record.rounds[0].new_seqs[1].score==189
assert record.rounds[0].new_seqs[1].e==3e-49
assert record.rounds[0].new_seqs[2].title=="gi|16131333|ref|NP_417918.1| sigma(32) factor [Escherichia coli ..."
assert record.rounds[0].new_seqs[2].score==84
assert record.rounds[0].new_seqs[2].e==2.0000000000000001e-17
assert record.rounds[0].new_seqs[3].title=="gi|33347603|ref|NP_416432.2| sigma factor 28 [Escherichia coli K12]"
assert record.rounds[0].new_seqs[3].score==63
assert record.rounds[0].new_seqs[3].e==5.0000000000000002e-11
assert len(record.rounds[0].alignments)==0
assert len(record.rounds[1].new_seqs)==0
assert len(record.rounds[1].alignments)==0
assert record.database_name==['NC_000913']
assert record.num_letters_in_database==[1350842]
assert record.num_sequences_in_database==[4311]
assert record.posted_date==[('Apr 8, 2004 12:56 AM',)]
assert record.ka_params==[0.313, 0.141, 0.379]
assert record.ka_params_gap==[0.267, 0.0431, 0.140]
assert record.matrix=='BLOSUM62'
assert record.gap_penalties==[11,1]
assert record.num_hits==2949500
assert record.num_sequences==4311
assert record.num_extends==123966
assert record.num_good_extends==504
assert record.num_seqs_better_e==6
assert record.hsps_no_gap==9
assert record.hsps_prelim_gapped==5
assert record.hsps_gapped==23
assert record.query_length==613
assert record.database_length==1350842
assert record.effective_hsp_length==90
assert record.effective_query_length==523
assert record.effective_database_length==962852
assert record.effective_search_space==503571596
assert record.effective_search_space_used==503571596
assert record.threshold==11
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.2)
assert record.gap_x_dropoff==(38,14.6)
assert record.gap_x_dropoff_final==(64,24.7)
assert record.gap_trigger==(41,21.3)
assert record.blast_cutoff==(64,29.2)

handle = open('Blast/bt102')
record = parser.parse(handle)
assert record.application=="TBLASTN"
assert record.version=='2.2.16'
assert record.date=="Mar-25-2007"
assert record.reference==reference
assert record.query=="gi|1710680|sp|P00579.2|RPOD_ECOLI RecName: Full=RNA polymerase\nsigma factor rpoD; AltName: Full=Sigma-70"
assert record.query_letters==613
assert record.database=="ecoli.nt"
assert record.database_sequences==400
assert record.database_letters==4662239
assert len(record.descriptions)==8
assert record.descriptions[0].title=="gb|AE000388.1|AE000388 Escherichia coli K-12 MG1655 section ..."
assert record.descriptions[0].score==933
assert record.descriptions[0].e==0.0
assert record.descriptions[0].num_alignments==2
assert record.descriptions[1].title=="gb|AE000358.1|AE000358 Escherichia coli K-12 MG1655 section ..."
assert record.descriptions[1].score==240
assert record.descriptions[1].e==1e-063
assert record.descriptions[1].num_alignments==1
assert record.descriptions[2].title=="gb|AE000422.1|AE000422 Escherichia coli K-12 MG1655 section ..."
assert record.descriptions[2].score==71
assert record.descriptions[2].e==2e-016
assert record.descriptions[2].num_alignments==2
assert record.descriptions[3].title=="gb|AE000285.1|AE000285 Escherichia coli K-12 MG1655 section ..."
assert record.descriptions[3].score==42
assert record.descriptions[3].e==1e-006
assert record.descriptions[3].num_alignments==2
assert record.descriptions[4].title=="gb|AE000397.1|AE000397 Escherichia coli K-12 MG1655 section ..."
assert record.descriptions[4].score==32
assert record.descriptions[4].e==0.50
assert record.descriptions[4].num_alignments==1
assert record.descriptions[5].title=="gb|AE000436.1|AE000436 Escherichia coli K-12 MG1655 section ..."
assert record.descriptions[5].score==29
assert record.descriptions[5].e==3.4 
assert record.descriptions[5].num_alignments==1
assert record.descriptions[6].title=="gb|AE000375.1|AE000375 Escherichia coli K-12 MG1655 section ..."
assert record.descriptions[6].score==29
assert record.descriptions[6].e==3.4
assert record.descriptions[6].num_alignments==1
assert record.descriptions[7].title=="gb|AE000405.1|AE000405 Escherichia coli K-12 MG1655 section ..."
assert record.descriptions[7].score==28
assert record.descriptions[7].e==6.3
assert record.descriptions[7].num_alignments==1
assert len(record.alignments)==8
assert record.alignments[0].title==">gb|AE000388.1|AE000388 Escherichia coli K-12 MG1655 section 278 of 400 of the complete genome"
assert record.alignments[0].length==10334
assert record.alignments[1].title==">gb|AE000358.1|AE000358 Escherichia coli K-12 MG1655 section 248 of 400 of the complete genome"
assert record.alignments[1].length==11457
assert record.alignments[2].title==">gb|AE000422.1|AE000422 Escherichia coli K-12 MG1655 section 312 of 400 of the complete genome"
assert record.alignments[2].length==10864
assert record.alignments[3].title==">gb|AE000285.1|AE000285 Escherichia coli K-12 MG1655 section 175 of 400 of the complete genome"
assert record.alignments[3].length==10165
assert record.alignments[4].title==">gb|AE000397.1|AE000397 Escherichia coli K-12 MG1655 section 287 of 400 of the complete genome"
assert record.alignments[4].length==14820
assert record.alignments[5].title==">gb|AE000436.1|AE000436 Escherichia coli K-12 MG1655 section 326 of 400 of the complete genome"
assert record.alignments[5].length==11150
assert record.alignments[6].title==">gb|AE000375.1|AE000375 Escherichia coli K-12 MG1655 section 265 of 400 of the complete genome"
assert record.alignments[6].length==10362
assert record.alignments[7].title==">gb|AE000405.1|AE000405 Escherichia coli K-12 MG1655 section 295 of 400 of the complete genome"
assert record.alignments[7].length==11095
assert record.alignments[0].hsps[0].score==2038
assert record.alignments[0].hsps[0].bits==933
assert record.alignments[0].hsps[0].expect==0.0
assert record.alignments[0].hsps[1].score==880
assert record.alignments[0].hsps[1].bits==404
assert record.alignments[0].hsps[1].expect==0.0
assert len(record.alignments[0].hsps)==2
assert record.alignments[1].hsps[0].score==519
assert record.alignments[1].hsps[0].bits==240
assert record.alignments[1].hsps[0].expect==1e-063
assert len(record.alignments[1].hsps)==1
assert record.alignments[2].hsps[0].score==149
assert record.alignments[2].hsps[0].bits==71.0
assert record.alignments[2].hsps[0].expect==2e-016
assert record.alignments[2].hsps[1].score==78
assert record.alignments[2].hsps[1].bits==38.6
assert record.alignments[2].hsps[1].expect==2e-016
assert len(record.alignments[2].hsps)==2
assert record.alignments[3].hsps[0].score==85
assert record.alignments[3].hsps[0].bits==41.7
assert record.alignments[3].hsps[0].expect==1e-006
assert record.alignments[3].hsps[1].score==51
assert record.alignments[3].hsps[1].bits==26.2
assert record.alignments[3].hsps[1].expect==1e-006
assert len(record.alignments[3].hsps)==2
assert record.alignments[4].hsps[0].score==63
assert record.alignments[4].hsps[0].bits==31.7
assert record.alignments[4].hsps[0].expect==0.50
assert len(record.alignments[4].hsps)==1
assert record.alignments[5].hsps[0].score==57
assert record.alignments[5].hsps[0].bits==29.0
assert record.alignments[5].hsps[0].expect==3.4
assert len(record.alignments[5].hsps)==1
assert record.alignments[6].hsps[0].score==57
assert record.alignments[6].hsps[0].bits==29.0
assert record.alignments[6].hsps[0].expect==3.4
assert len(record.alignments[6].hsps)==1
assert record.alignments[7].hsps[0].score==55
assert record.alignments[7].hsps[0].bits==28.0
assert record.alignments[7].hsps[0].expect==6.3
assert len(record.alignments[7].hsps)==1
assert record.alignments[0].hsps[0].identities==(404, 404)
assert record.alignments[0].hsps[0].positives==(404, 404)
assert record.alignments[0].hsps[1].identities==(177, 187)
assert record.alignments[0].hsps[1].positives==(177, 187)
assert record.alignments[1].hsps[0].identities==(96, 225)
assert record.alignments[1].hsps[0].positives==(157, 225)
assert record.alignments[2].hsps[0].identities==(31, 78)
assert record.alignments[2].hsps[0].positives==(46, 78)
assert record.alignments[2].hsps[1].identities==(14, 35)
assert record.alignments[2].hsps[1].positives==(25, 35)
assert record.alignments[3].hsps[0].identities==(17, 37)
assert record.alignments[3].hsps[0].positives==(26, 37)
assert record.alignments[3].hsps[1].identities==(10, 19)
assert record.alignments[3].hsps[1].positives==(14, 19)
assert record.alignments[4].hsps[0].identities==(12, 37)
assert record.alignments[4].hsps[0].positives==(22, 37)
assert record.alignments[5].hsps[0].identities==(12, 26)
assert record.alignments[5].hsps[0].positives==(17, 26)
assert record.alignments[6].hsps[0].identities==(11, 41)
assert record.alignments[6].hsps[0].positives==(20, 41)
assert record.alignments[7].hsps[0].identities==(12, 50)
assert record.alignments[7].hsps[0].positives==(24, 50)
assert record.alignments[0].hsps[0].frame==("+1", )
assert record.alignments[0].hsps[1].frame==("+1", )
assert record.alignments[1].hsps[0].frame==("-2", )
assert record.alignments[2].hsps[0].frame==("-3", )
assert record.alignments[2].hsps[1].frame==("-3", )
assert record.alignments[3].hsps[0].frame==("-1", )
assert record.alignments[3].hsps[1].frame==("-1", )
assert record.alignments[4].hsps[0].frame==("-3", )
assert record.alignments[5].hsps[0].frame==("+1", )
assert record.alignments[6].hsps[0].frame==("+3", )
assert record.alignments[7].hsps[0].frame==("+2", )
assert record.alignments[0].hsps[0].query=="NSIDPELAREKFAELRAQYVVTRDTIKAKGRSHATAQEEILKLSEVFKQFRLVPKQFDYLVNSMRVMMDRVRTQERLIMKLCVEQCKMPKKNFITLFTGNETSDTWFNAAIAMNKPWSEKLHDVSEEVHRALQKLQQIEEETGLTIEQVKDINRRMSIGEAKARRAKKEMVEANLRLVISIAKKYTNRGLQFLDLIQEGNIGLMKAVDKFEYRRGYKFSTYATWWIRQAITRSIADQARTIRIPVHMIETINKLNRISRQMLQEMGREPTPEELAERMLMPEDKIRKVLKIAKEPISMETPIGDDEDSHLGDFIEDTTLELPLDSATTESLRAATHDVLAGLTAREAKVLRMRFGIDMNTDYTLEEVGKQFDVTRERIRQIEAKALRKLRHPSRSEVLRSFLDD"
assert record.alignments[0].hsps[0].match=="NSIDPELAREKFAELRAQYVVTRDTIKAKGRSHATAQEEILKLSEVFKQFRLVPKQFDYLVNSMRVMMDRVRTQERLIMKLCVEQCKMPKKNFITLFTGNETSDTWFNAAIAMNKPWSEKLHDVSEEVHRALQKLQQIEEETGLTIEQVKDINRRMSIGEAKARRAKKEMVEANLRLVISIAKKYTNRGLQFLDLIQEGNIGLMKAVDKFEYRRGYKFSTYATWWIRQAITRSIADQARTIRIPVHMIETINKLNRISRQMLQEMGREPTPEELAERMLMPEDKIRKVLKIAKEPISMETPIGDDEDSHLGDFIEDTTLELPLDSATTESLRAATHDVLAGLTAREAKVLRMRFGIDMNTDYTLEEVGKQFDVTRERIRQIEAKALRKLRHPSRSEVLRSFLDD"
assert record.alignments[0].hsps[0].sbjct=="NSIDPELAREKFAELRAQYVVTRDTIKAKGRSHATAQEEILKLSEVFKQFRLVPKQFDYLVNSMRVMMDRVRTQERLIMKLCVEQCKMPKKNFITLFTGNETSDTWFNAAIAMNKPWSEKLHDVSEEVHRALQKLQQIEEETGLTIEQVKDINRRMSIGEAKARRAKKEMVEANLRLVISIAKKYTNRGLQFLDLIQEGNIGLMKAVDKFEYRRGYKFSTYATWWIRQAITRSIADQARTIRIPVHMIETINKLNRISRQMLQEMGREPTPEELAERMLMPEDKIRKVLKIAKEPISMETPIGDDEDSHLGDFIEDTTLELPLDSATTESLRAATHDVLAGLTAREAKVLRMRFGIDMNTDYTLEEVGKQFDVTRERIRQIEAKALRKLRHPSRSEVLRSFLDD"
assert record.alignments[0].hsps[0].query_start==210
assert record.alignments[0].hsps[0].query_end==613
assert record.alignments[0].hsps[0].sbjct_start==7345
assert record.alignments[0].hsps[0].sbjct_end==8556
assert record.alignments[0].hsps[1].query=="MEQNPQSQLKLLVTRGKEQGYLTYAEVNDHLPEDIVDSDQIEDIIQMINDMGIQVMEEAPDADDLMLAENTXXXXXXXXXXQVLSSVESEIGRTTDPVRMYMREMGTVELLTREGEIDIAKRIEDGINQVQCSVAEYPEAITYLLEQYDRVEAEEARLSDLITGFVDPNAEEDLAPTATHVGSELSQ"
assert record.alignments[0].hsps[1].match=="MEQNPQSQLKLLVTRGKEQGYLTYAEVNDHLPEDIVDSDQIEDIIQMINDMGIQVMEEAPDADDLMLAENT          QVLSSVESEIGRTTDPVRMYMREMGTVELLTREGEIDIAKRIEDGINQVQCSVAEYPEAITYLLEQYDRVEAEEARLSDLITGFVDPNAEEDLAPTATHVGSELSQ"
assert record.alignments[0].hsps[1].sbjct=="MEQNPQSQLKLLVTRGKEQGYLTYAEVNDHLPEDIVDSDQIEDIIQMINDMGIQVMEEAPDADDLMLAENTADEDAAEAAAQVLSSVESEIGRTTDPVRMYMREMGTVELLTREGEIDIAKRIEDGINQVQCSVAEYPEAITYLLEQYDRVEAEEARLSDLITGFVDPNAEEDLAPTATHVGSELSQ"
assert record.alignments[0].hsps[1].query_start==1
assert record.alignments[0].hsps[1].query_end==187
assert record.alignments[0].hsps[1].sbjct_start==6718
assert record.alignments[0].hsps[1].sbjct_end==7278
assert record.alignments[1].hsps[0].query=="AKKEMVEANLRLVISIAKKYTNRGLQFLDLIQEGNIGLMKAVDKFEYRRGYKFSTYATWWIRQAITRSIADQARTIRIPVHMIETINKLNRISRQMLQEMGREPTPEELAERMLMPEDKIRKVLKIAKEPISMETPIGDDEDSHLGDFIEDTTLELPLDSATTESLRAATHDVLAGLTAREAKVLRMRFGIDMNTDYTLEEVGKQFDVTRERIRQIEAKALRKLR"
assert record.alignments[1].hsps[0].match=="+++ M+E+NLRLV+ IA++Y NRGL  LDLI+EGN+GL++AV+KF+  RG++FSTYATWWIRQ I R+I +Q RTIR+P+H+++ +N   R +R++  ++  EP+ EE+AE++  P D + ++L++ +   S++TP+G D +  L D + D     P D+   + ++ +    L  L A++ +VL  RFG+      TLE+VG++  +TRER+RQI+ + LR+LR"
assert record.alignments[1].hsps[0].sbjct=="SRRRMIESNLRLVVKIARRYGNRGLALLDLIEEGNLGLIRAVEKFDPERGFRFSTYATWWIRQTIERAIMNQTRTIRLPIHIVKELNVYLRTARELSHKLDHEPSAEEIAEQLDKPVDDVSRMLRLNERITSVDTPLGGDSEKALLDILADEKENGPEDTTQDDDMKQSIVKWLFELNAKQREVLARRFGLLGYEAATLEDVGREIGLTRERVRQIQVEGLRRLR"
assert record.alignments[1].hsps[0].query_start==375
assert record.alignments[1].hsps[0].query_end==599
assert record.alignments[1].hsps[0].sbjct_start==2258
assert record.alignments[1].hsps[0].sbjct_end==1584
assert record.alignments[2].hsps[0].query=="AKKEMVEANLRLVISIAKKYTNRGLQFLDLIQEGNIGLMKAVDKFEYRRGYKFSTYATWWIRQAITRSIADQARTIRI"
assert record.alignments[2].hsps[0].match=="A K ++ ++LR V+ IA+ Y   GL   DLIQEGNIGLMKAV +F    G +  ++A  WI+  I   +    R +++"
assert record.alignments[2].hsps[0].sbjct=="AAKTLILSHLRFVVHIARNYAGYGLPQADLIQEGNIGLMKAVRRFNPEVGVRLVSFAVHWIKAEIHEYVLRNWRIVKV"
assert record.alignments[2].hsps[0].query_start==375
assert record.alignments[2].hsps[0].query_end==452
assert record.alignments[2].hsps[0].sbjct_start==2759
assert record.alignments[2].hsps[0].sbjct_end==2526
assert record.alignments[2].hsps[1].query=="IDMNTDYTLEEVGKQFDVTRERIRQIEAKALRKLR"
assert record.alignments[2].hsps[1].match=="+D +   TL+E+  ++ V+ ER+RQ+E  A++KLR"
assert record.alignments[2].hsps[1].sbjct=="LDEDNKSTLQELADRYGVSAERVRQLEKNAMKKLR"
assert record.alignments[2].hsps[1].query_start==565
assert record.alignments[2].hsps[1].query_end==599
assert record.alignments[2].hsps[1].sbjct_start==2171
assert record.alignments[2].hsps[1].sbjct_end==2067
assert record.alignments[3].hsps[0].query=="DLIQEGNIGLMKAVDKFEYRRGYKFSTYATWWIRQAI"
assert record.alignments[3].hsps[0].match=="DL+Q G IGL+ AV++++  +G  F+TYA   IR A+"
assert record.alignments[3].hsps[0].sbjct=="DLLQAGGIGLLNAVERYDALQGTAFTTYAVQRIRGAM"
assert record.alignments[3].hsps[0].query_start==403
assert record.alignments[3].hsps[0].query_end==439
assert record.alignments[3].hsps[0].sbjct_start==1264
assert record.alignments[3].hsps[0].sbjct_end==1154
assert record.alignments[3].hsps[1].query=="QMLQEMGREPTPEELAERM"
assert record.alignments[3].hsps[1].match=="Q+ QE+GR  T  E+AER+"
assert record.alignments[3].hsps[1].sbjct=="QLEQELGRNATETEVAERL"
assert record.alignments[3].hsps[1].query_start==469
assert record.alignments[3].hsps[1].query_end==487
assert record.alignments[3].hsps[1].sbjct_start==1075
assert record.alignments[3].hsps[1].sbjct_end==1019
assert record.alignments[4].hsps[0].query=="ISMETPIGDDEDSHLGDFIEDTTLELPLDSATTESLR"
assert record.alignments[4].hsps[0].match=="I++E    +DE  +LGD++ED    +  D  TT++ +"
assert record.alignments[4].hsps[0].sbjct=="ITLEAARYEDESLNLGDYVEDQIESVTFDRITTQTAK"
assert record.alignments[4].hsps[0].query_start==505
assert record.alignments[4].hsps[0].query_end==541
assert record.alignments[4].hsps[0].sbjct_start==13927
assert record.alignments[4].hsps[0].sbjct_end==13817
assert record.alignments[5].hsps[0].query=="ATTESLRAATHDVLAGLTAREAKVLR"
assert record.alignments[5].hsps[0].match=="ATT +LR + H++  GL  R+ KV R"
assert record.alignments[5].hsps[0].sbjct=="ATTRNLRWSWHEITWGLCLRQTKVCR"
assert record.alignments[5].hsps[0].query_start==535
assert record.alignments[5].hsps[0].query_end==560
assert record.alignments[5].hsps[0].sbjct_start==1762
assert record.alignments[5].hsps[0].sbjct_end==1839
assert record.alignments[6].hsps[0].query=="KFSTYATWWIRQAITRSIADQARTIRIPVHMIETINKLNRI"
assert record.alignments[6].hsps[0].match=="K S+  TWW   A+ R +  +    R+P  +  ++   +RI"
assert record.alignments[6].hsps[0].sbjct=="KLSSITTWWHLAALPRRVRRKPYPPRLPAELTNSMRPKSRI"
assert record.alignments[6].hsps[0].query_start==426
assert record.alignments[6].hsps[0].query_end==466
assert record.alignments[6].hsps[0].sbjct_start==5193
assert record.alignments[6].hsps[0].sbjct_end==5315
assert record.alignments[7].hsps[0].query=="FDYLVNSMRVMMDRVRTQERLIMKLCVEQCKMPKKNFITLFTGNETSDTW"
assert record.alignments[7].hsps[0].match=="F+  VNS  + +  V +Q+++  K   +  K   KN+     G + +D +"
assert record.alignments[7].hsps[0].sbjct=="FEQRVNSDVLTVSTVNSQDQVTQKPLRDSVKQALKNYFAQLNGQDVNDLY"
assert record.alignments[7].hsps[0].query_start==266
assert record.alignments[7].hsps[0].query_end==315
assert record.alignments[7].hsps[0].sbjct_start==1070
assert record.alignments[7].hsps[0].sbjct_end==1219
assert record.database_name==['ecoli.nt']
assert record.num_letters_in_database==[4662239]
assert record.num_sequences_in_database==[400]
assert record.posted_date==[('May 31, 2007 11:31 AM',)]
assert record.ka_params==[0.317, 0.132, 0.363]
assert record.matrix=='BLOSUM62'
assert record.num_sequences==400
assert record.num_hits==2109753
assert record.num_extends==23326
assert record.num_good_extends==120
assert record.num_seqs_better_e==8
assert record.query_length==613
assert record.database_length==1554079
assert record.effective_query_length==570
assert record.effective_database_length==1536879
assert record.effective_search_space==876021030
assert record.effective_search_space_used==876021030
assert record.threshold==13
assert record.window_size==40
assert record.dropoff_1st_pass==(16,7.3)
assert record.blast_cutoff==(52,26.7)
