from testformats.swissprot38 import *


import support

test_list = support.Storage()
add_test = test_list.add_test
add_test_lines = test_list.add_test_lines

add_test_lines("ID line", ID, """\
ID   100K_RAT       STANDARD;      PRT;   889 AA.
ID   CYC_BOVIN      STANDARD;      PRT;   104 AA.
ID   GIA2_GIALA     STANDARD;      PRT;   296 AA.
""")

add_test_lines("AC line", AC, """\
AC   Q62671;
AC   P00321; P05348;
""")

add_test("AC (block)", AC_block, """\
AC   Q62671; Q05349; Q05351; Q05352; Q05353; Q05354; Q05355; Q05356;
AC   Q92671; Q95349; Q95351; Q95352; Q95353; Q95354; Q95355; Q95356;
AC   Q98763;
""")

add_test("date 1", DT_created + DT_seq_update + DT_ann_update, """\
DT   01-OCT-1996 (Rel. 34, Created)
DT   01-OCT-1996 (Rel. 34, Last sequence update)
DT   01-NOV-1997 (Rel. 35, Last annotation update)
""")

add_test("date 2", DT_created + DT_seq_update + DT_ann_update, """\
DT   01-AUG-1988 (Rel. 08, Created)
DT   01-JAN-1990 (Rel. 13, Last sequence update)
DT   15-APR-1999 (Rel. 38, Last annotation update)
""")

add_test_lines("DE (single line)", DE, """\
DE   100 KD PROTEIN (EC 6.3.2.-).
DE   10 KD PROTEIN PRECURSOR (CLONE PSAS10).
""")
add_test("DE (muliline) 1", DE_block, """\
DE   14-3-3 PROTEIN BETA/ALPHA (PROTEIN KINASE C INHIBITOR PROTEIN-1)
DE   (KCIP-1).
""")
add_test("DE (muliline) 2", DE_block, """\
DE   ANNEXIN V (LIPOCORTIN V) (ENDONEXIN II) (CALPHOBINDIN I) (CBP-I)
DE   (PLACENTAL ANTICOAGULANT PROTEIN I) (PAP-I) (PP4) (THROMBOPLASTIN
DE   INHIBITOR) (VASCULAR ANTICOAGULANT-ALPHA) (VAC-ALPHA) (ANCHORIN CII).
""")

add_test_lines("GN (single line)", GN, """\
GN   HAG3.
GN   REX-1.
GN   HNS OR DRDX OR OSMZ OR BGLY.
GN   GVPA AND (GVPB OR GVPA2).
""")

# from CALM_HUMAN
add_test("GN (block)", GN_block, """\
GN   (CALM1 OR CAM1 OR CALM OR CAM) AND (CALM2 OR CAM2 OR CAMB) AND
GN   (CALM3 OR CAM3 OR CAMC).
""")

add_test_lines("OS (single line)", OS, """\
OS   Helianthus annuus (Common sunflower).
OS   Escherichia coli.
OS   Homo sapiens (Human).
OS   Acer spicatum (Moose maple) (Mountain maple).
OS   Rous sarcoma virus (strain Schmidt-Ruppin).
""")

add_test("OS (block) 1", OS_block, """\
OS   Oncorhynchus nerka (Sockeye salmon), and
OS   Oncorhynchus masou (Cherry salmon) (Masu salmon).
""")

add_test("OS (block) 2", OS_block, """\
OS   Mus musculus (Mouse), Rattus norvegicus (Rat), and
OS   Bos taurus (Bovine).
""")

add_test_lines("OG (single line)", OG, """\
OG   Chloroplast.
OG   Cyanelle.
OG   Mitochondrion.
OG   Plasmid name.
OG   Plasmid IncI1 ColIb.
""")

add_test("OG (block)", OG_block, """\
OG   Plasmid pDGO100, Plasmid IncQ pIE723, Plasmid pBP201, and
OG   Plasmid IncM pBWH1.
""")

add_test("OG (block)", OG_block, """\
OG   Plasmid R6-5, Plasmid IncFII NR1, and
OG   Plasmid IncFII R1-19 (R1 drd-19).
""")
add_test_lines("OC (single line)", OC, """\
OC   Eukaryota; Alveolata; Apicomplexa; Haemosporida; Plasmodium.
OC   Eukaryota; Entamoebidae; Entamoeba.
""")

add_test("OC (block) 1", OC_block, """\
OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
OC   euphyllophytes; Spermatophyta; Magnoliophyta; eudicotyledons;
OC   core eudicots; Asteridae; euasterids I; Solanales; Solanaceae;
OC   Solanum.
""")

add_test("OC (block) 2", OC_block, """\
OC   Eukaryota; Metazoa; Chordata; Vertebrata; Mammalia; Eutheria;
OC   Primates; Catarrhini; Hominidae; Homo.
""")

add_test_lines("RN", RN, """\
RN   [1]
RN   [2]
RN   [3]
RN   [23]
RN   [876543]
""")

add_test_lines("RP", RP, """\
RP   SEQUENCE FROM N.A.
RP   SEQUENCE FROM N.A., AND SEQUENCE OF 12-35.
RP   SEQUENCE OF 34-56; 67-73 AND 123-345, AND DISULFIDE BONDS.
RP   REVISIONS TO 67-89.
RP   STRUCTURE BY NMR.
RP   X-RAY CRYSTALLOGRAPHY (1.8 ANGSTROMS).
RP   CHARACTERIZATION.
RP   MUTAGENESIS OF TYR-56.
RP   REVIEW.
RP   VARIANT ALA-58.
RP   VARIANTS XLI LEU-341; ARG-372 AND TYR-446.
""")

add_test_lines("RC (single line)", RC, """\
RC   STRAIN=SPRAGUE-DAWLEY; TISSUE=LIVER;
RC   STRAIN=HOLSTEIN; TISSUE=MAMMARY GLAND, AND LYMPH NODE;
RC   SPECIES=RAT; STRAIN=WISTAR;
RC   PLASMID=INCFII R100;
""")

add_test("RC (block)", RC_block, """\
RC   STRAIN=MVZ CATALOG 172969, 172970, 174109, 174110, 174229, AND 174230;
RC   TISSUE=LIVER;
""")

add_test_lines("RX (single line)", RX, """\
RX   MEDLINE; 91002678.
RX   MEDLINE; 93144687.
""")

add_test("RA (block) 1", RA_block, """\
RA   SMITH H. JR., VON BRAUN M.T. III;
""")

add_test("RA (block) 2", RA_block, """\
RA   YANOFSKY C., PLATT T., CRAWFORD I.P., NICHOLS B.P., CHRISTIE G.E.,
RA   HOROWITZ H., VAN CLEEMPUT M., WU A.M.;
""")
add_test("RT (single line)", RT, """\
RT   "Organization of the sunflower 11S storage protein gene family.";
""")

add_test("RT (block) 1", RT_block, """\
RT   "New insulin-like proteins with atypical disulfide bond pattern
RT   characterized in Caenorhabditis elegans by comparative sequence
RT   analysis and homology modeling.";
""")

add_test("RT (block) 2", RT_block, """\
RT   "Stored mRNA in cotyledons of Vigna unguiculata seeds: nucleotide
RT   sequence of cloned cDNA for a stored mRNA and induction of its
RT   synthesis by precocious germination.";
""")
add_test_lines("RL (single line)", RL, """\
RL   J. Mol. Biol. 168:321-331(1983).
RL   Nucleic Acids Res. 27:0-0(1999).
RL   Thesis (1972), University of Geneva, Switzerland.
""")

add_test("RL (block) 1", RL_block, """\
RL   (In) Boyer P.D. (eds.);
RL   The enzymes (3rd ed.), pp.11:397-547, Academic Press, New York (1975).
""")

add_test("RL (block) 2", RL_block, """\
RL   (In) Rich D.H., Gross E. (eds.);
RL   Proceedings of the 7th american peptide symposium, pp.69-72,
RL   Pierce Chemical Co., Rockford Il. (1981).
""")

add_test("RL (block) 3", RL_block, """\
RL   (In) Magnusson S., Ottesen M., Foltmann B., Dano K.,
RL   Neurath H. (eds.);
RL   Regulatory proteolytic enzymes and their inhibitors, pp.163-172,
RL   Pergamon Press, New York (1978).
""")

add_test("RL (block) 4", RL_block, """\
RL   (In) Plant Gene Register PGR98-023.
RL   (In) Worm Breeder's Gazette 15(3):34(1998).
""")

add_test("reference 1", reference, """\
RN   [1]
RP   SEQUENCE FROM N.A.
RC   STRAIN=WISTAR; TISSUE=TESTIS;
RX   MEDLINE; 92253337.
RA   MUELLER D., REHBEIN M., BAUMEISTER H., RICHTER D.;
RT   "Molecular characterization of a novel rat protein structurally
RT   related to poly(A) binding proteins and the 70K protein of the U1
RT   small nuclear ribonucleoprotein particle (snRNP).";
RL   Nucleic Acids Res. 20:1471-1475(1992).
""")

add_test("reference 2", reference, """\
RN   [2]
RP   ERRATUM.
RA   MUELLER D., REHBEIN M., BAUMEISTER H., RICHTER D.;
RL   Nucleic Acids Res. 20:2624-2624(1992).
""")
s1 = """\
CC   -!- FUNCTION: E3 UBIQUITIN-PROTEIN LIGASE WHICH ACCEPTS UBIQUITIN FROM
CC       AN E2 UBIQUITIN-CONJUGATING ENZYME IN THE FORM OF A THIOESTER AND
CC       THEN DIRECTLY TRANSFERS THE UBIQUITIN TO TARGETED SUBSTRATES (BY
CC       SIMILARITY). THIS PROTEIN MAY BE INVOLVED IN MATURATION AND/OR
CC       POST-TRANSCRIPTIONAL REGULATION OF MRNA.
"""
s2 = """\
CC   -!- TISSUE SPECIFICITY: HIGHEST LEVELS FOUND IN TESTIS. ALSO PRESENT
CC       IN LIVER, KIDNEY, LUNG AND BRAIN.
"""
s3 = """\
CC   -!- DEVELOPMENTAL STAGE: IN EARLY POST-NATAL LIFE, EXPRESSION IN
CC       THE TESTIS INCREASES TO REACH A MAXIMUM AROUND DAY 28.
"""
add_test("single comment 1", single_comment, s1)
add_test("single comment 2", single_comment, s2)
add_test("single comment 3", single_comment, s3)

copyright = """\
CC   --------------------------------------------------------------------------
CC   This SWISS-PROT entry is copyright. It is produced through a collaboration
CC   between  the Swiss Institute of Bioinformatics  and the  EMBL outstation -
CC   the European Bioinformatics Institute.  There are no  restrictions on  its
CC   use  by  non-profit  institutions as long  as its content  is  in  no  way
CC   modified and this statement is not removed.  Usage  by  and for commercial
CC   entities requires a license agreement (See http://www.isb-sib.ch/announce/
CC   or send an email to license@isb-sib.ch).
CC   --------------------------------------------------------------------------
"""                     


add_test("set of comments 1", comment, s1)
add_test("set of comments 2", comment, s1+s2)
add_test("set of comments 3", comment, s1+s3)
add_test("set of comments 4", comment, s1+s2+s3)
add_test("set of comments 5", comment, s1+copyright)
add_test("set of comments 6", comment, s1+s2+s3+copyright)


s = """\
DR   AARHUS/GHENT-2DPAGE; 8006; IEF.
DR   DICTYDB; DD01047; MYOA.
DR   ECO2DBASE; G052.0; 6TH EDITION.
DR   ECOGENE; EG10054; ARAC.
DR   FLYBASE; FBgn0000055; Adh.
DR   GCRDB; GCR_0087; -.
DR   HIV; K02013; NEF$BRU.
DR   HSC-2DPAGE; P47985; HUMAN.
DR   HSSP; P00438; 1DOB.
DR   MAIZEDB; 25342; -.
DR   MAIZE-2DPAGE; P80607; COLEOPTILE.
DR   MENDEL; 2596; AMAhy;psbA;1.
DR   MGD; MGI:87920; ADFP.
DR   MGD; MGI:95401; EPB4.1.
DR   MIM; 249900; -.
DR   PDB; 3ADK; 16-APR-88.
DR   PIR; A02768; R5EC7.
DR   REBASE; RB00005; EcoRI.
DR   SGD; L0000008; AAR2.
DR   STYGENE; SG10312; PROV.
DR   SUBTILIST; BG10774; OPPD.
DR   SWISS-2DPAGE; P10599; HUMAN.
DR   TIGR; MJ0125; -.
DR   TRANSFAC; T00141; -.
DR   WORMPEP; ZK637.7; CE00437.
DR   YEPD; 4270; -.
DR   ZFIN; ZDB-GENE-980526-290; hoxa1.
DR   EMBL; Y00312; CAA68412.1; -.
DR   EMBL; L29151; AAA99430.1; ALT_INIT.
DR   EMBL; L20562; AAA26884.1; ALT_TERM.
DR   EMBL; X56420; CAA39814.1; ALT_FRAME.
DR   EMBL; M28482; AAA26378.1; ALT_SEQ.
DR   EMBL; M63397; AAA51662.1; -.
DR   EMBL; M63395; AAA51662.1; JOINED.
DR   EMBL; M63396; AAA51662.1; JOINED.
DR   EMBL; J04126; -; NOT_ANNOTATED_CDS.
DR   PROSITE; PS00107; PROTEIN_KINASE_ATP; 1.
DR   PROSITE; PS00028; ZINC_FINGER_C2H2; 6.
DR   PROSITE; PS00237; G_PROTEIN_RECEPTOR; FALSE_NEG.
DR   PROSITE; PS01128; SHIKIMATE_KINASE; PARTIAL.
DR   PROSITE; PS00383; TYR_PHOSPHATASE_1; UNKNOWN_1.
DR   PFAM; PF00017; SH2; 1.
DR   PFAM; PF00008; EGF; 8.
DR   PFAM; PF00595; PDZ; PARTIAL.
"""

add_test_lines("DR", DR, s)
add_test("DR (block)", DR_block, s)

add_test_lines("KW (single line)", KW, """\
KW   Oxidoreductase; Acetylation.
KW   Acetylation; Oxidoreductase.
KW   Ubiquitin conjugation; Ligase.
KW   Signal.
KW   Seed storage protein; Multigene family; Signal.
""")

add_test("KW (block) 1", KW_block, """\
KW   Brain; Neurone; Phosphorylation; Acetylation; Multigene family;
KW   3D-structure.
""")

add_test("KW (block) 2", KW_block, """\
KW   Steroidogenesis; Oxidoreductase; NAD; Isomerase; Mitochondrion;
KW   Multigene family; Multifunctional enzyme; Transmembrane;
KW   Endoplasmic reticulum.
""")

add_test("KW (block) 3", KW_block, """\
KW   Hydrolase; Ligase; Oxidoreductase; NADP; Multifunctional enzyme;
KW   One-carbon metabolism; ATP-binding; Purine biosynthesis;
KW   Amino-acid biosynthesis; Methionine biosynthesis;
KW   Histidine biosynthesis.
""")

add_test_lines("FT range / single line", FT_range, """\
FT   DOMAIN       77     88       ASP/GLU-RICH (ACIDIC).
FT   DOMAIN      127    150       PRO-RICH.
FT   DOMAIN      420    439       ARG/GLU-RICH (MIXED CHARGE).
FT   BINDING     858    858       UBIQUITIN (BY SIMILARITY).
FT   DOMAIN       43     57       PRO/THR-RICH.
FT   SIGNAL       <1      8       BY SIMILARITY.
FT   NON_TER       1      1
FT   DISULFID     56     67
FT   CARBOHYD    114    114       POTENTIAL.
FT   CONFLICT    102    102       D -> S (IN REF. 2).
FT   CONFLICT    105    105       MISSING (IN REF. 3).
FT   CHAIN         ?     75       10 KD PROTEIN.
FT   SIGNAL        1    ?24       POTENTIAL.
FT   PROPEP      ?25    ?31       POTENTIAL.
FT   SIGNAL        1      ?
FT   INIT_MET      0      0
""")

add_test("FT w/ continuation 1", FT, """\
FT   MOD_RES       9      9       AMIDATION (G-10 PROVIDE AMIDE GROUP)
FT                                (BY SIMILARITY).
""")

add_test("FT w/ continuation 2", FT, """\
FT   DOMAIN      131    296       13.5 X 12 AA TANDEM REPEATS OF E-E-T-Q-K-
FT                                T-V-E-P-E-Q-T.
""")

add_test("FT w/ continuation 3", FT, """\
FT   VARIANT      33     33       F -> Y (IN A*0205, A*0206, A*0208, A*0210
FT                                AND A*0221).
FT                                /FTId=VAR_004334.
""")


add_test("feature (block)", feature_block, """\
FT   DOMAIN       77     88       ASP/GLU-RICH (ACIDIC).
FT   DOMAIN      127    150       PRO-RICH.
FT   DOMAIN      420    439       ARG/GLU-RICH (MIXED CHARGE).
FT   DOMAIN      131    296       13.5 X 12 AA TANDEM REPEATS OF E-E-T-Q-K-
FT                                T-V-E-P-E-Q-T.
FT   BINDING     858    858       UBIQUITIN (BY SIMILARITY).
FT   DOMAIN       43     57       PRO/THR-RICH.
FT   SIGNAL       <1      8       BY SIMILARITY.
FT   NON_TER       1      1
FT   DISULFID     56     67
FT   CARBOHYD    114    114       POTENTIAL.
FT   VARIANT      33     33       F -> Y (IN A*0205, A*0206, A*0208, A*0210
FT                                AND A*0221).
FT                                /FTId=VAR_004334.
FT   CONFLICT    102    102       D -> S (IN REF. 2).
FT   CONFLICT    105    105       MISSING (IN REF. 3).
FT   CHAIN         ?     75       10 KD PROTEIN.
FT   SIGNAL        1    ?24       POTENTIAL.
FT   PROPEP      ?25    ?31       POTENTIAL.
FT   MOD_RES       9      9       AMIDATION (G-10 PROVIDE AMIDE GROUP)
FT                                (BY SIMILARITY).
FT   SIGNAL        1      ?
FT   INIT_MET      0      0
""")

add_test_lines("SQ header", SQ, """\
SQ   SEQUENCE   889 AA;  100368 MW;  DD7E6C7A CRC32;
SQ   SEQUENCE   111 AA;  12416 MW;  103BBA8B CRC32;
SQ   SEQUENCE   29 AA;  2900 MW;  BA38C516 CRC32;
SQ   SEQUENCE   1707 AA;  194328 MW;  31FDA77C CRC32;
""")

add_test_lines("SQ_data", SQ_data, """\
     ISFTSFNDES GENAEKLLQF KRWFWSIVER MSMTERQDLV YFWTSSPSLP ASEEGFQPMP
     SITIRPPDDQ HLPTANTCIS RLYVPLYSSK QILKQKLLLA IKTKNFGFV
     SITIRPPDDQ HLP
     A
""")

add_test("sequence 1", sequence, """\
SQ   SEQUENCE   889 AA;  100368 MW;  DD7E6C7A CRC32;
     MMSARGDFLN YALSLMRSHN DEHSDVLPVL DVCSLKHVAY VFQALIYWIK AMNQQTTLDT
     PQLERKRTRE LLELGIDNED SEHENDDDTS QSATLNDKDD ESLPAETGQN HPFFRRSDSM
     S
""")
add_test("sequence 2", sequence, """\
SQ   SEQUENCE   4 AA;  408 MW;  34BC4AD8 CRC32;
     GFAD
""")

add_test("end", end, """\
//
""")


record1 = """ID   100K_RAT       STANDARD;      PRT;   889 AA.
AC   Q62671;
DT   01-NOV-1997 (Rel. 35, Created)
DT   01-NOV-1997 (Rel. 35, Last sequence update)
DT   15-JUL-1999 (Rel. 38, Last annotation update)
DE   100 KD PROTEIN (EC 6.3.2.-).
OS   Rattus norvegicus (Rat).
OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Mammalia;
OC   Eutheria; Rodentia; Sciurognathi; Muridae; Murinae; Rattus.
RN   [1]
RP   SEQUENCE FROM N.A.
RC   STRAIN=WISTAR; TISSUE=TESTIS;
RX   MEDLINE; 92253337.
RA   MUELLER D., REHBEIN M., BAUMEISTER H., RICHTER D.;
RT   "Molecular characterization of a novel rat protein structurally
RT   related to poly(A) binding proteins and the 70K protein of the U1
RT   small nuclear ribonucleoprotein particle (snRNP).";
RL   Nucleic Acids Res. 20:1471-1475(1992).
RN   [2]
RP   ERRATUM.
RA   MUELLER D., REHBEIN M., BAUMEISTER H., RICHTER D.;
RL   Nucleic Acids Res. 20:2624-2624(1992).
CC   -!- FUNCTION: E3 UBIQUITIN-PROTEIN LIGASE WHICH ACCEPTS UBIQUITIN FROM
CC       AN E2 UBIQUITIN-CONJUGATING ENZYME IN THE FORM OF A THIOESTER AND
CC       THEN DIRECTLY TRANSFERS THE UBIQUITIN TO TARGETED SUBSTRATES (BY
CC       SIMILARITY). THIS PROTEIN MAY BE INVOLVED IN MATURATION AND/OR
CC       POST-TRANSCRIPTIONAL REGULATION OF MRNA.
CC   -!- TISSUE SPECIFICITY: HIGHEST LEVELS FOUND IN TESTIS. ALSO PRESENT
CC       IN LIVER, KIDNEY, LUNG AND BRAIN.
CC   -!- DEVELOPMENTAL STAGE: IN EARLY POST-NATAL LIFE, EXPRESSION IN
CC       THE TESTIS INCREASES TO REACH A MAXIMUM AROUND DAY 28.
CC   -!- MISCELLANEOUS: A CYSTEINE RESIDUE IS REQUIRED FOR
CC       UBIQUITIN-THIOLESTER FORMATION.
CC   -!- SIMILARITY: CONTAINS AN HECT-TYPE E3 UBIQUITIN-PROTEIN LIGASE
CC       DOMAIN.
CC   -!- SIMILARITY: A CENTRAL REGION (AA 485-514) IS SIMILAR TO THE
CC       C-TERMINAL DOMAINS OF MAMMALIAN AND YEAST POLY (A) RNA BINDING
CC       PROTEINS (PABP).
CC   -!- SIMILARITY: THE C-TERMINAL HALF SHOWS HIGH SIMILARITY TO
CC       DROSOPHILA HYPERPLASMIC DISC PROTEIN AND SOME, TO HUMAN E6-AP.
CC   -!- SIMILARITY: CONTAINS MIXED-CHARGE DOMAINS SIMILAR TO RNA-BINDING
CC       PROTEINS.
CC   --------------------------------------------------------------------------
CC   This SWISS-PROT entry is copyright. It is produced through a collaboration
CC   between  the Swiss Institute of Bioinformatics  and the  EMBL outstation -
CC   the European Bioinformatics Institute.  There are no  restrictions on  its
CC   use  by  non-profit  institutions as long  as its content  is  in  no  way
CC   modified and this statement is not removed.  Usage  by  and for commercial
CC   entities requires a license agreement (See http://www.isb-sib.ch/announce/
CC   or send an email to license@isb-sib.ch).
CC   --------------------------------------------------------------------------
DR   EMBL; X64411; CAA45756.1; -.
DR   PFAM; PF00632; HECT; 1.
DR   PFAM; PF00658; PABP; 1.
DR   PROSITE; PS00107; PROTEIN_KINASE_ATP; 1.
DR   AARHUS/GHENT-2DPAGE; 8006; IEF.
DR   DICTYDB; DD01047; MYOA.
KW   Ubiquitin conjugation; G-protein coupled receptor; Transmembrane;
KW   Glycoprotein; Ligase.
FT   DOMAIN       77     88       ASP/GLU-RICH (ACIDIC).
FT   DOMAIN      127    150       PRO-RICH.
FT   DOMAIN      420    439       ARG/GLU-RICH (MIXED CHARGE).
FT   DOMAIN      448    457       ARG/ASP-RICH (MIXED CHARGE).
FT   DOMAIN      485    514       PABP-LIKE.
FT   DOMAIN      579    590       ASP/GLU-RICH (ACIDIC).
FT   DOMAIN      786    889       HECT DOMAIN.
FT   DOMAIN      827    847       PRO-RICH.
FT   BINDING     858    858       UBIQUITIN (BY SIMILARITY).
SQ   SEQUENCE   889 AA;  100368 MW;  DD7E6C7A CRC32;
     MMSARGDFLN YALSLMRSHN DEHSDVLPVL DVCSLKHVAY VFQALIYWIK AMNQQTTLDT
     PQLERKRTRE LLELGIDNED SEHENDDDTS QSATLNDKDD ESLPAETGQN HPFFRRSDSM
     TFLGCIPPNP FEVPLAEAIP LADQPHLLQP NARKEDLFGR PSQGLYSSSA GSGKCLVEVT
     MDRNCLEVLP TKMSYAANLK NVMNMQNRQK KAGEDQSMLA EEADSSKPGP SAHDVAAQLK
     SSLLAEIGLT ESEGPPLTSF RPQCSFMGMV ISHDMLLGRW RLSLELFGRV FMEDVGAEPG
     SILTELGGFE VKESKFRREM EKLRNQQSRD LSLEVDRDRD LLIQQTMRQL NNHFGRRCAT
     TPMAVHRVKV TFKDEPGEGS GVARSFYTAI AQAFLSNEKL PNLDCIQNAN KGTHTSLMQR
     LRNRGERDRE REREREMRRS SGLRAGSRRD RDRDFRRQLS IDTRPFRPAS EGNPSDDPDP
     LPAHRQALGE RLYPRVQAMQ PAFASKITGM LLELSPAQLL LLLASEDSLR ARVEEAMELI
     VAHGRENGAD SILDLGLLDS SEKVQENRKR HGSSRSVVDM DLDDTDDGDD NAPLFYQPGK
     RGFYTPRPGK NTEARLNCFR NIGRILGLCL LQNELCPITL NRHVIKVLLG RKVNWHDFAF
     FDPVMYESLR QLILASQSSD ADAVFSAMDL AFAVDLCKEE GGGQVELIPN GVNIPVTPQN
     VYEYVRKYAE HRMLVVAEQP LHAMRKGLLD VLPKNSLEDL TAEDFRLLVN GCGEVNVQML
     ISFTSFNDES GENAEKLLQF KRWFWSIVER MSMTERQDLV YFWTSSPSLP ASEEGFQPMP
     SITIRPPDDQ HLPTANTCIS RLYVPLYSSK QILKQKLLLA IKTKNFGFV
//
"""

record2 = """\
ID   12KD_FRAAN     STANDARD;      PRT;   111 AA.
AC   Q05349;
DT   01-OCT-1996 (Rel. 34, Created)
DT   01-OCT-1996 (Rel. 34, Last sequence update)
DT   01-NOV-1997 (Rel. 35, Last annotation update)
DE   AUXIN-REPRESSED 12.5 KD PROTEIN.
OS   Fragaria ananassa (Strawberry).
OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
OC   euphyllophytes; Spermatophyta; Magnoliophyta; eudicotyledons;
OC   core eudicots; Rosidae; eurosids I; Rosales; Rosaceae; Fragaria.
RN   [1]
RP   SEQUENCE FROM N.A.
RC   STRAIN=CV. OZARK BEAUTY; TISSUE=FLOWER;
RX   MEDLINE; 91329668.
RA   REDDY A.S.N., POOVAIAH B.W.;
RT   "Molecular cloning and sequencing of a cDNA for an auxin-repressed
RT   mRNA: correlation between fruit growth and repression of the
RT   auxin-regulated gene.";
RL   Plant Mol. Biol. 14:127-136(1990).
CC   -!- FUNCTION: E3 UBIQUITIN-PROTEIN LIGASE WHICH ACCEPTS UBIQUITIN FROM
CC       AN E2 UBIQUITIN-CONJUGATING ENZYME IN THE FORM OF A THIOESTER AND
CC       THEN DIRECTLY TRANSFERS THE UBIQUITIN TO TARGETED SUBSTRATES (BY
CC       SIMILARITY). THIS PROTEIN MAY BE INVOLVED IN MATURATION AND/OR
CC       POST-TRANSCRIPTIONAL REGULATION OF MRNA.
CC   -!- TISSUE SPECIFICITY: HIGHEST LEVELS FOUND IN TESTIS. ALSO PRESENT
CC       IN LIVER, KIDNEY, LUNG AND BRAIN.
CC   -!- DEVELOPMENTAL STAGE: IN EARLY POST-NATAL LIFE, EXPRESSION IN
CC       THE TESTIS INCREASES TO REACH A MAXIMUM AROUND DAY 28.
CC   -!- MISCELLANEOUS: A CYSTEINE RESIDUE IS REQUIRED FOR
CC       UBIQUITIN-THIOLESTER FORMATION.
CC   -!- SIMILARITY: CONTAINS AN HECT-TYPE E3 UBIQUITIN-PROTEIN LIGASE
CC       DOMAIN.
CC   -!- SIMILARITY: A CENTRAL REGION (AA 485-514) IS SIMILAR TO THE
CC       C-TERMINAL DOMAINS OF MAMMALIAN AND YEAST POLY (A) RNA BINDING
CC       PROTEINS (PABP).
CC   -!- SIMILARITY: THE C-TERMINAL HALF SHOWS HIGH SIMILARITY TO
CC       DROSOPHILA HYPERPLASMIC DISC PROTEIN AND SOME, TO HUMAN E6-AP.
CC   -!- SIMILARITY: CONTAINS MIXED-CHARGE DOMAINS SIMILAR TO RNA-BINDING
CC       PROTEINS.
CC   --------------------------------------------------------------------------
CC   This SWISS-PROT entry is copyright. It is produced through a collaboration
CC   between  the Swiss Institute of Bioinformatics  and the  EMBL outstation -
CC   the European Bioinformatics Institute.  There are no  restrictions on  its
CC   use  by  non-profit  institutions as long  as its content  is  in  no  way
CC   modified and this statement is not removed.  Usage  by  and for commercial
CC   entities requires a license agreement (See http://www.isb-sib.ch/announce/
CC   or send an email to license@isb-sib.ch).
CC   --------------------------------------------------------------------------
DR   EMBL; X52429; CAA36676.1; -.
DR   EMBL; X64411; CAA45756.1; -.
DR   PFAM; PF00632; HECT; 1.
DR   PFAM; PF00658; PABP; 1.
DR   PROSITE; PS00107; PROTEIN_KINASE_ATP; 1.
DR   AARHUS/GHENT-2DPAGE; 8006; IEF.
DR   DICTYDB; DD01047; MYOA.
KW   Ubiquitin conjugation; G-protein coupled receptor; Transmembrane;
KW   Glycoprotein; Ligase.
FT   DOMAIN       43     57       PRO/THR-RICH.
SQ   SEQUENCE   111 AA;  12416 MW;  103BBA8B CRC32;
     MVLLDKLWDD IVAGPQPERG LGMLRKVPQP LNLKDEGESS KITMPTTPTT PVTPTTPISA
     RKDNVWRSVF HPGSNLSSKT MGNQVFDSPQ PNSPTVYDWM YSGETRSKHH R
//
"""

add_test("record 1", record, record1)
add_test("record 2", record, record2)
add_test("format", format, record1 + record2)


def test():
    test_list.test()

def dump():
    test_list.dump()
       
if __name__ == "__main__":
    test()

