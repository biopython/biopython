#!/usr/bin/env python
"""Test the Martel-based genbank_format, now found in Bio.expressions.genbank
"""
from Bio.expressions.genbank import *

import martel_support
import Martel

test_list = martel_support.Storage()
add_test = test_list.add_test
add_test_lines = test_list.add_test_lines


header_s = """\
GBPRI8.SEQ           Genetic Sequence Data Bank
                          15 October 2000

                NCBI-GenBank Flat File Release 120.0

                         Primate Sequences (Part 8)

   28652 loci,    36017768 bases, from    28652 reported sequences


"""

add_test("header 1", header, header_s)

add_test_lines("locus", locus_line, """\
LOCUS       HUAF001549 202004 bp    DNA             PRI       28-JUL-1998
LOCUS       HUB384D8   139934 bp    DNA             PRI       03-FEB-2000
LOCUS       HUM12LIPO     239 bp    mRNA            PRI       27-APR-1993
LOCUS       HUM17QYACA    225 bp ss-RNA             PRI       12-APR-1996
LOCUS       HUMALPPC2      53 bp ds-mRNA            PRI       18-DEC-1998
LOCUS       HUMMGAT       648 bp    DNA   circular  PRI       26-MAR-1997
LOCUS       S40706        895 bp                    PRI       08-MAY-1993
LOCUS       HUMMTTSF1    1040 bp    mRNA  circular  PRI       26-JUL-1995
LOCUS       18_Char_LOCUS_Name 99999999999 bp ss-snRNA circular HTC DD-MMM-YYYY
LOCUS       AB000383                  5423 bp    DNA   circular VRL 05-FEB-1999
LOCUS       NM_007355               2459 bp    mRNA    linear   PRI 02-NOV-2000
LOCUS       AL123456  4411529 bp          circular  BCT       07-JUL-1998
""")

add_test("definition 1", definition_block, """\
DEFINITION  Human Chromosome 16 BAC clone CIT987SK-A-270G1, complete sequence.
""")

add_test("definition 2", definition_block, """\
DEFINITION  Homo sapiens chromosome 22q13 BAC clone CIT987SK-384D8 complete
            sequence.
""")

add_test("definition 3", definition_block, """\
DEFINITION  Human cytoplasmic acid phosphatase B (ACP1) gene, exons 2-3S;
            cytoplasmic acid phosphatase A (ACP1) gene; cytoplasmic acid
            phosphatase C (ACP1) gene.
""")

add_test("accession 1", accession_block, """\
ACCESSION   L38283 AF001549
""")
add_test("accession 2", accession_block, """\
ACCESSION   M24309 K01883 M11832
""")
add_test("accession 3", accession_block, """\
ACCESSION   M24317 K01883 M11831 M11832 M11833 M11834 M11835 M11836 M11837
            M11838 M11839
""")
add_test("accession 4", accession_block, """\
ACCESSION   U01317 J00179 J00093 J00094 J00096 J00158 J00159 J00160 J00161
            J00162 J00163 J00164 J00165 J00166 J00167 J00168 J00169 J00170
            J00171 J00172 J00173 J00174 J00175 J00177 J00178 K01239 K01890
            K02544 M18047 M19067 M24868 M24886 X00423 X00424 X00672
""")

add_test_lines("version", version_line, """\
VERSION     Z96953.1  GI:2791932
VERSION     AC000062.1  GI:1669374
VERSION     D14718.1  GI:285928
VERSION     AL123456
""")

# The documentation says:
#    Each line in theEach line in the keywords field ends in a
#    semicolon; keywords field ends in a semicolon;
# This shows it isn't true
add_test("keywords 1", keywords_block, """\
KEYWORDS    HTG.
""")
add_test("keywords 2", keywords_block, """\
KEYWORDS    .
""")
add_test("keywords 3", keywords_block, """\
KEYWORDS    immunoglobulin heavy chain subgroup VH-III.
""")
add_test("keywords 4", keywords_block, """\
KEYWORDS    cell surface glycoprotein; class II gene; integral membrane
            protein; major histocompatibility complex.
""")

add_test_lines("segment", segment_line, """\
SEGMENT     1 of 39
SEGMENT     987 of 989
""")

add_test("source 1", source_block, """\
SOURCE      Homo sapiens DNA.
""")
add_test("source 2", source_block, """\
SOURCE      Homo sapiens (tissue library: Subclones in pOT2 from P1 clone H11)
            DNA.
""")
add_test("source 3", source_block, """\
SOURCE      Homo sapiens (clone: clone SEL226a) (clone library: 17q YAC
            (303G8)) fetus brain/ thymus/ spleen/ testis/ total fetus cDNA to
            other RNA.
""")

add_test("organism 1", organism_block, """\
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Homo.
""")
add_test("organism 2", organism_block, """\
  ORGANISM  Pan paniscus
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Pan.
""")
add_test("organism 3", organism_block, """\
  ORGANISM  Mitochondrion Varecia variegata rubra
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Primates; Strepsirhini; Lemuridae; Varecia.
""")
add_test("organism 4", organism_block, """\
  ORGANISM  Macaca mulatta
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Primates; Catarrhini; Cercopithecidae;
            Cercopithecinae; Macaca.
""")

add_test("reference_line 1", reference_line, """\
REFERENCE   3  (bases 1 to 2722)
""")
add_test("reference_line 2", reference_line, """\
REFERENCE   2  (bases 1 to 60; 1561 to 3000; 3241 to 3300; 6721 to 6780; 7201
            to 7265)
""")
add_test("reference_line 3", reference_line, """\
REFERENCE   14 (sites)
""")
add_test("reference_line 4", reference_line, """\
REFERENCE   2
""")

add_test("authors 1", authors_block, """\
  AUTHORS   Adams,M.D.
""")
add_test("authors 2", authors_block, """\
  AUTHORS   Loftus,B.J., Kim,U.J., Sneddon,V.P., Kalush,F., Brandon,R.,
            Fuhrmann,J., Mason,T., Crosby,M.L., Barnstead,M., Cronin,L.,
            Deslattes Mays,A., Cao,Y., Xu,R.X., Kang,H.L., Mitchell,S.,
            Eichler,E.E., Harris,P.C., Venter,J.C. and Adams,M.D.
""")
add_test("authors 3", authors_block, """\
  AUTHORS   Mahraoui,L., Takeda,J., Mesonero,J., Chantret,I., Dussaulx,E.,
            Bell,G.I. and Brot-Laroche,E.
""")

add_test("title 1", title_block, """\
  TITLE     Regulation of expression of the human fructose transporter (GLUT5)
            by cyclic AMP
""")
add_test("title 2", title_block, """\
  TITLE     Direct Submission
""")

add_test("journal 1", journal_block, """\
  JOURNAL   Genomics 60 (3), 295-308 (1999)
""")
add_test("journal 2", journal_block, """\
  JOURNAL   Unpublished
""")
add_test("journal 3", journal_block, """\
  JOURNAL   Submitted (26-JUN-1996) The Institute for Genomic Research, 9712
            Medical Center Dr, Rockville, Maryland 20850, USA
""")
add_test("journal 4", journal_block, """\
  JOURNAL   Biochem. J. 301 (Pt 1), 169-175 (1994)
""")
add_test("journal 5", journal_block, """\
  JOURNAL   Submitted (24-JAN-1994) Jun Takeda, Howard Hughes Medical
            Institute, The University of Chicago, 5841 S. Maryland Ave.,
            Chicago, IL 60637, USA
""")

add_test_lines("medline", medline_line, """\
  MEDLINE   92355116
  MEDLINE   93329075
  MEDLINE   86313652
""")
add_test_lines("medline", pubmed_line, """\
   PUBMED   10493829
   PUBMED   9384599
   PUBMED   8486766
""")

add_test("remark 1", remark_block, """\
  REMARK    sequence updated by GenBank staff 05/21/99 to correspond to
            published sequence
""")
add_test("remark 2", remark_block, """\
  REMARK    Sequence update by submitter
""")
add_test("remark 3", remark_block, """\
  REMARK    Erratum:[J Biol Chem 1988 Aug 5;263(22):11016]
""")
add_test("remark 4", remark_block, """\
  REMARK    This entry shows the sequence encoding a part of the heavy chain of
            the human IgG1 antibody LuNm01 (NOT the LuNm03 antibody as stated
            in Fig. 2 of Reference 2). The light chain variable domain of
            LuNm01 is encoded by GenBank Accession Number M97803.
""")

add_test("reference 1", reference, """\
REFERENCE   1  (bases 1 to 139934)
  AUTHORS   Loftus,B.J., Kim,U.J., Sneddon,V.P., Kalush,F., Brandon,R.,
            Fuhrmann,J., Mason,T., Crosby,M.L., Barnstead,M., Cronin,L.,
            Deslattes Mays,A., Cao,Y., Xu,R.X., Kang,H.L., Mitchell,S.,
            Eichler,E.E., Harris,P.C., Venter,J.C. and Adams,M.D.
  TITLE     Genome duplications and other features in 12 Mb of DNA sequence
            from human chromosome 16p and 16q
  JOURNAL   Genomics 60 (3), 295-308 (1999)
  MEDLINE   99425270
   PUBMED   10493829
""")
add_test("reference 2", reference, """\
REFERENCE   1  (bases 1 to 422)
  AUTHORS   Tikka,L., Elomaa,O., Pihlajaniemi,T. and Tryggvason,K.
  TITLE     Human alpha 1 (XIII) collagen gene. Multiple forms of the gene
            transcripts are generated through complex alternative splicing of
            several short exons
  JOURNAL   J. Biol. Chem. 266 (26), 17713-17719 (1991)
  MEDLINE   91373404
""")
add_test("reference 3", reference, """\
REFERENCE   2  (bases 1 to 270)
  AUTHORS   Pimtanothai,N. and Hurley,C.K.
  TITLE     Direct Submission
  JOURNAL   Submitted (27-AUG-1999) Microbiology & Immunology, Georgetown
            University, 3970 Reservoir Rd. N.W., Washington, D.C. 20007, USA
""")
add_test("reference 3", reference, """\
REFERENCE   1  (bases 1 to 159)
  AUTHORS   Weber,J.L.
  JOURNAL   Unpublished (1993)
""")
add_test("reference 4", reference, """\
REFERENCE   1  (bases 1 to 1381)
  AUTHORS   Rosen,G.D., Birkenmeier,T.M. and Dean,D.C.
  TITLE     Characterization of the alpha 4 integrin gene promoter
  JOURNAL   Proc. Natl. Acad. Sci. U.S.A. 88 (10), 4094-4098 (1991)
  MEDLINE   91239513
  REMARK    sequence updated by GenBank staff 05/21/99 to correspond to
            published sequence
""")

add_test("comment 1", comment_block, """\
COMMENT     On May 21, 1999 this sequence version replaced gi:177901.
""")
add_test("comment 2", comment_block, """\
COMMENT     [1]  revises [2].
            Printed copy of sequence for [2],[1] kindly provided by W.L.Miller,
            08-DEC-1987.
""")
add_test("comment 3", comment_block, """\
COMMENT     A draft entry and printed copy of the sequence in [1] were kindly
            provided by A.Yoshida, 30-MAY-1986.
            The other human class I ADH1 alpha subunit sequence is found under
            accession M11307.
""")

add_test("features_line", features_line, """\
FEATURES             Location/Qualifiers
""")

add_test("feature_key_location 1", feature_key_line, """\
     source          1..1400
""")
add_test("feature_key_location 2", feature_key_line, """\
     mRNA            5..>93
""")
add_test("feature_key_location 3", feature_key_line, """\
     gene            join(M24308.1:23..375,M24309.1:1..132,M24310.1:1..169,
                     M24311.1:1..118,M24312.1:1..250,M24313.1:1..291,
                     M24314.1:1..166,M24315.1:1..169,1..631)
""")
add_test("feature_key_location 4", feature_key_line, """\
     CDS             join(L38283.1:602..619,L38284.1:432..533,
                     L38285.1:265..403,L38286.1:99..186,L38286.1:284..503,
                     L38287.1:23..283,L38288.1:139..274,L38289.1:122..260,
                     155..179)
""")
add_test("feature_key_location 5", feature_key_line, """\
     intron          order(M81117.1:182..196,1..15)
""")

add_test("qualifier_block 1", qualifier, """\
                     /gene="ADH5"
""")
add_test("qualifier_block 2", qualifier, """\
                     /product="alcohol dehydrogenase"
""")
add_test("qualifier_block 3", qualifier, """\
                     /translation="MANEVIKCKAAVAWEAGKPLSIEEIEVAPPKAHEVRIKIIATAV
                     CHTDAYTLSGADPEGCFPVILGHEGAGIVESVGEGVTKLKAGDTVIPLYIPQCGECKF
                     CLNPKTNLCQKIRVTQGKGLMPDGTSRFTCKGKTILHYMGTSTFSEYTVVADISVAKI
                     DPLAPLDKVCLLGCGISTGYGAAVNTAKLEPGSVCAVFGLGGVGLAVIMGCKVAGASR
                     IIGVDINKDKFARAKEFGATECINPQDFSKPIQEVLIEMTDGGVDYSFECIGNVKVMR
                     AALEACHKGWGVSVVVGVAASGEEIATRPFQLVTGRTWKGTAFGGWKSVESVPKLVSE
                     YMSKKIKVDEFVTHNLSFDEINKAFELMHSGKSIRTVVKI"
""")
add_test("qualifier_block 4", qualifier, """\
                     /note="this region is complementary to part of the first
                     monomer unit of a consensus Alu repeat; putative"
""")
add_test("qualifier_block 5", qualifier, """\
                     /citation=[2]
""")
add_test("qualifier_block 6", qualifier, """\
                     /partial
""")
add_test("qualifier_block 7", qualifier, """\
                     /clone="lambda [537,488,659,130,475, and 476]"
""")

add_test("feature_key_block 1", feature, """\
     source          1..1624
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
""")
add_test("feature_key_block 2", feature, """\
     repeat_region   complement(21..183)
                     /partial
                     /note="this region is complementary to part of the first
                     monomer unit of a consensus Alu repeat; putative"
                     /citation=[2]
                     /rpt_family="Alu"
""")
add_test("feature_key_block 3", feature, """\
     GC_signal       complement(1175..1180)
                     /note="putative"
                     /citation=[1]
                     /citation=[2]
""")
add_test("feature_key_block 4", feature, """\
     exon            118..358
                     /gene="ADP-ribosylation factor 3"
                     /number=2
""")
add_test("features 1", features_line + Martel.Rep1(feature), """\
FEATURES             Location/Qualifiers
     source          1..144
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
                     /germline
                     /map="Chromosome 4"
     exon            1..144
                     /note="5' altered transcript of AF-4 gene. The sequence
                     joined to AF-4 cDNA sequence at 436nt."
""")
add_test("features 2", features_line + Martel.Rep1(feature), """\
FEATURES             Location/Qualifiers
     source          1..8588
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
                     /cell_line="HuH-7"
                     /clone="lambda [537,488,659,130,475, and 476]"
                     /tissue_type="hepatoma"
     CDS             130..8481
                     /codon_start=1
                     /evidence=not_experimental
                     /product="alpha-fetoprotein enhancer binding protein"
                     /protein_id="BAA01095.1"
                     /db_xref="GI:219430"
                     /translation="MRLGGGQLVSEELMNLGESFIQTNDPSLKLFQCAVCNKFTTDNL
                     DMLGLHMNVERSLSEDEWKAVMGDSYQCKLCRYNTQLKANFQLHCKTDKHVQKYQLVA
                     HIKEGGKANEWRLKCVAIGNPVHLKCNACDYYTNSLEKLRLHTVNSRHEASLKLYKHL
                     QQHESGVEGESCYYHCVLCNYSTKAKLNLIQHVRSMKHQRSESLRKLQRLQKGLPEED
                     EDLGQIFTIRRCPSTDPEEAIEDVEGPSETAADPEELAKDQEGGASSSQAEKELTDSP
                     ATSKRISFPGSSESPLSSKRPKTAEEIKPEQMYQCPYCKYSNADVNRLRVHAMTQHSV
                     QPMLRCPLCQDMLNNKIHLQLHLTHLHSVAPDCVEKLIMTVTTPEMVMPSSMFLPAAV
                     PPLSSSSTVTSSSCSTSGVQPSMPTDDYSEESDTDLSQKSDGPASPVEGPKDPSCPKD
                     SGLTSVGTDTFRL"
     variation       2316
                     /note="c in lambda 537; t in lambda 488"
                     /replace="t"
     variation       4338
                     /note="t in lambda 659; c in lambda 488"
                     /replace="c"
     variation       4378
                     /note="deletion in lambda 488"
                     /replace=""
     polyA_signal    8497..8502
                     /evidence=not_experimental
     polyA_site      8521
                     /note="'poly(A) tail is attached to the cytosine residue
                     at 8521 inlambda 476'"
                     /evidence=not_experimental
     polyA_signal    8529..8534
                     /evidence=not_experimental
     polyA_site      8588
                     /note="'poly(A) tail is attatched to the cytosine residue
                     at 8588 inlambda 475'"
                     /evidence=not_experimental
""")

add_test("features 3", features_line + Martel.Rep1(feature), """\
FEATURES             Location/Qualifiers
     source          1..159
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
                     /sex="Male"
                     /tissue_type="Blood"
     primer_bind     1..24
                     /standard_name="PCR primer"
                     /evidence=experimental
     repeat_region   29..78
                     /rpt_family="(dC-dA)n.(dG-dT)n"
                     /rpt_type=tandem
                     /evidence=experimental
     primer_bind     complement(87..105)
                     /standard_name="PCR primer"
                     /evidence=experimental
""")

add_test_lines("base_count", base_count_line, """\
BASE COUNT       53 a     82 c     96 g     39 t
BASE COUNT       83 a     60 c     45 g     83 t      1 others
BASE COUNT       93 a     70 c     48 g     74 t      6 others
BASE COUNT        0 a      1 c      1 g      1 t
""")

# No spaces / spaces / "." (even though doc says it isn't needed) / text
add_test_lines("origin", origin_line, """\
ORIGIN
ORIGIN      
ORIGIN      .
ORIGIN      1 bp upstream of AvaII site; chromosome 6p23-p21.1.
""")

# the 's' is in HUM1PLYSAT ; the 'r' is in HUM3PFEPO
#  (present to make sure you don't think [atcgn] is sufficient :)
add_test_lines("sequence", Martel.Rep1(sequence_line), """\
        1 aactcggtac taggaaaact cctattttaa aatccagccc tgagtgggaa gatttgggaa
   201901 ggcagctctg cgacctcctc cctgcagcgg cccaggtggg aactcagcca gggcagcggc
   201961 gggggtcaca gtccccgcct gggacttcct atctgtcgaa gctt
      241 ttttnagntt caccatcata tatntccctc tgtaattnaa cctcttgcca caccactg
      121 tatgcctctg accatgggnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn
      301 gagttcgcga ccagcctggc caacatgnta aaascccgtg ctctacta
      301 aaggagaggg aggccagagg acgggtcctt ttgggagttt tggrggctgg taacagctgc
""")

record_s1 = """\
LOCUS       HUM20N        209 bp    DNA             PRI       12-JUN-1993
DEFINITION  Human (clone: pHyTM1/20(N)) DNA sequence.
ACCESSION   L15245
VERSION     L15245.1  GI:291592
KEYWORDS    .
SOURCE      Homo sapiens DNA.
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Homo.
REFERENCE   1  (bases 1 to 209)
  AUTHORS   Farr,C.J., Stevanovic,M., Thomson,E.J., Goodfellow,P.N. and
            Cooke,H.J.
  TITLE     Telomere-associated chromosome fragmentation: applications in
            genome manipulation and analysis
  JOURNAL   Nature Genet. 2, 275-282 (1992)
  MEDLINE   93265159
REFERENCE   2  (bases 1 to 209)
  AUTHORS   Bayne,R.A.L., Taggart,M., Farr,C. and Cooke,H.J.
  TITLE     Sequence of end-clones recovered from TC8 hybrid deletion panel:
            (R) or (N) denotes sequence from Eco RI or Nde I end of insert
            where the Nde I end is proximal to the plasmid insertion site
  JOURNAL   Unpublished (1993)
FEATURES             Location/Qualifiers
     source          1..209
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
                     /cell_line="HTM18TC8 (a CHO/Human X hybrid line)"
BASE COUNT       83 a     23 c     45 g     58 t
ORIGIN      
        1 agtatcatag ccacacctaa cttcaaggag taggaaagta cgatcctact gtgaccttct
       61 aaggaggaga accagaaatg ttgatgaaca gcaattatta atgagaatgt ataattgaga
      121 aaaaaggtac aaaaagagta ggtataaaga gagtattttt caaagatgtg tgtgtacttg
      181 gatttatatt ataaatgtat ggacaatat
//
"""
                     
record_s2 = """\
LOCUS       HUM21DC4Z    3540 bp    DNA             PRI       10-MAY-1995
DEFINITION  Homo sapiens (subclone 10_b2 from P1 H21) DNA sequence.
ACCESSION   L42095
VERSION     L42095.1  GI:804735
KEYWORDS    Interleukin growth hormone cluster on chromosome 5 (5q31).
SOURCE      Homo sapiens (tissue library: Subclones in Homo sapiens from P1
            clone H21) DNA.
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Homo.
REFERENCE   1  (bases 1 to 3540)
  AUTHORS   Martin,C.H., Bondoc,M.M., Chiang,A., Cloutier,T., Davis,C.A.,
            Ericsson,C.L., Jaklevic,M.A., Kim,R.J., Lee,M.T., Li,M.,
            Mayeda,C.A., Steiert-El Kheir,A. and Palazzolo,M.J.
  TITLE     Sequencing of the interleukin growth hormone cluster on chromosome
            5 (5q31) of homo sapiens
  JOURNAL   Unpublished (1995)
COMMENT     Sequence submitted by:
            Human Genome Center and
            Drosophila Genome Center
            Lawrence Berkeley Laboratory
            Berkeley, CA 94720
            e-mail: seq@genome.lbl.gov
            This subclone overlaps H21 1_b3 and H21 10_d2.
            The P1, from which this subclone is derived, is adjacent to P1
            (1857) and  (5005).
FEATURES             Location/Qualifiers
     source          1..3540
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
                     /tissue_lib="Subclones in Homo sapiens from P1 clone H21"
BASE COUNT      967 a    694 c    608 g   1271 t
ORIGIN      
        1 atgttactac tgtaattgtt tgtgtgccac aaaccatcca catataagag gtgaacttaa
       61 tccattaacg tgtgtgtcct gactgcttta ctgacctgcc attcccgtct ctctccctct
      121 ccttggaacc tgattgcctg agacacaata atatggaaat taggccaatt agtaacccta
      181 caacagcccc taagtgttta agcgaaagaa gagtcaaaca tctcgtttta aatcaaaaac
      241 tagaaatgat taagcttagt tgagaaaagc atgtcaaaat ccaaaacagg ttgaaagtta
//
"""
record_s3 = """\
LOCUS       HUM28HLA1     270 bp    DNA             PRI       29-SEP-1999
DEFINITION  Homo sapiens MHC class I antigen (HLA-A) gene, HLA-A*0213 allele,
            exon 2.
ACCESSION   AF181101
VERSION     AF181101.1  GI:5932391
KEYWORDS    .
SEGMENT     1 of 2
SOURCE      human.
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Homo.
REFERENCE   1  (bases 1 to 270)
  AUTHORS   Steiner,N.K., Hurley,C.K. and Koester,R.P.
  TITLE     Novel HLA-A allele
  JOURNAL   Unpublished
REFERENCE   2  (bases 1 to 270)
  AUTHORS   Steiner,N.K., Hurley,C.K. and Koester,R.P.
  TITLE     Direct Submission
  JOURNAL   Submitted (25-AUG-1999) Microbiology and Immunology, Georgetown
            University Medical Center, 3970 Reservoir Road,NW, Washington, D.C.
            20007, USA
FEATURES             Location/Qualifiers
     source          1..270
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
                     /clone="GN00286"
     exon            1..270
                     /gene="HLA-A"
                     /number=2
BASE COUNT       53 a     82 c     96 g     39 t
ORIGIN
        1 gctctcactc catgaggtat ttcttcacat ccgtgtcccg gcccggccgc ggggagcccc
       61 gcttcatcgc agtgggctac gtggacgaca cgcagttcgt gcggttcgac agcgacgccg
      121 cgagccagag gatggagccg cgggcgccgt ggatagagca ggagggtccg gagtattggg
      181 acggggagac acggaaagtg aaggcccact cacagactca ccgagtggac ctggggaccc
      241 tgcgcggcta ctacaaccag agcgaggccg
//
"""
add_test("record 1", record, record_s1)
add_test("record 2", record, record_s2)
add_test("record 3", record, record_s3)

add_test("format", format, record_s1+record_s2+record_s3)
add_test("ncbi_format", ncbi_format, header_s + record_s1+record_s2+record_s3)

def run_tests(args):
    test_list.test()

def dump():
    test_list.dump()

if __name__ == "__main__":
    test_list.test()
