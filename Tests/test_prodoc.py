# Copyright 1999 by Cayte Lindner.  All rights reserved.
# Copyright 2009 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest
from Bio.ExPASy import Prodoc




class TestProdocRead(unittest.TestCase):

    def test_read_pdoc00100(self):
        "Reading Prodoc record PDOC00100"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00100.txt')
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00100")
        self.assertEqual(len(record.prosite_refs), 4)
        self.assertEqual(record.prosite_refs[0], ("PS00107", "PROTEIN_KINASE_ATP"))
        self.assertEqual(record.prosite_refs[1], ("PS00108", "PROTEIN_KINASE_ST"))
        self.assertEqual(record.prosite_refs[2], ("PS00109", "PROTEIN_KINASE_TYR"))
        self.assertEqual(record.prosite_refs[3], ("PS50011", "PROTEIN_KINASE_DOM"))
        self.assertEqual(record.text, """\
******************************************
* Protein kinases signatures and profile *
******************************************

Eukaryotic  protein kinases [1 to 5]  are  enzymes  that   belong  to  a  very
extensive family of  proteins which share a conserved catalytic core common to
both serine/threonine and  tyrosine protein kinases.  There  are  a  number of
conserved regions in the catalytic domain of protein kinases. We have selected
two of these regions to build signature patterns.  The  first region, which is
located in the N-terminal extremity of the catalytic domain, is a glycine-rich
stretch of residues in the vicinity  of a lysine residue, which has been shown
to be involved in ATP binding.   The second  region,  which is  located in the
central part of the  catalytic  domain,  contains  a  conserved  aspartic acid
residue  which is important for the catalytic activity  of  the enzyme [6]; we
have derived  two signature patterns for that region: one specific for serine/
threonine kinases  and  the  other  for  tyrosine kinases. We also developed a
profile which is based on the alignment in [1] and covers the entire catalytic
domain.

-Consensus pattern: [LIV]-G-{P}-G-{P}-[FYWMGSTNH]-[SGA]-{PW}-[LIVCAT]-{PD}-x-
                    [GSTACLIVMFY]-x(5,18)-[LIVMFYWCSTAR]-[AIVP]-[LIVMFAGCKR]-K
                    [K binds ATP]
-Sequences known to belong to this class detected by the pattern: the majority
 of known  protein  kinases  but it fails to find a number of them, especially
 viral kinases  which  are  quite  divergent in this region and are completely
 missed by this pattern.
-Other sequence(s) detected in Swiss-Prot: 42.

-Consensus pattern: [LIVMFYC]-x-[HY]-x-D-[LIVMFY]-K-x(2)-N-[LIVMFYCT](3)
                    [D is an active site residue]
-Sequences known to belong to this class detected by the pattern: Most serine/
 threonine  specific protein  kinases  with  10 exceptions (half of them viral
 kinases) and  also  Epstein-Barr  virus BGLF4 and Drosophila ninaC which have
 respectively Ser and Arg instead of the conserved Lys and which are therefore
 detected by the tyrosine kinase specific pattern described below.
-Other sequence(s) detected in Swiss-Prot: 1.

-Consensus pattern: [LIVMFYC]-{A}-[HY]-x-D-[LIVMFY]-[RSTAC]-{D}-{PF}-N-
                    [LIVMFYC](3)
                    [D is an active site residue]
-Sequences known to belong to this class detected by the pattern: ALL tyrosine
 specific protein  kinases  with  the  exception of human ERBB3 and mouse blk.
 This pattern    will    also    detect    most    bacterial    aminoglycoside
 phosphotransferases [8,9]  and  herpesviruses ganciclovir kinases [10]; which
 are proteins structurally and evolutionary related to protein kinases.
-Other sequence(s) detected in Swiss-Prot: 17.

-Sequences known to belong to this class detected by the profile: ALL,  except
 for three  viral  kinases.  This  profile  also  detects  receptor  guanylate
 cyclases (see   <PDOC00430>)  and  2-5A-dependent  ribonucleases.    Sequence
 similarities between  these  two  families  and the eukaryotic protein kinase
 family have been noticed before. It also detects Arabidopsis thaliana kinase-
 like protein TMKL1 which seems to have lost its catalytic activity.
-Other sequence(s) detected in Swiss-Prot: 4.

-Note: If a protein  analyzed  includes the two protein kinase signatures, the
 probability of it being a protein kinase is close to 100%
-Note: Eukaryotic-type protein  kinases  have  also  been found in prokaryotes
 such as Myxococcus xanthus [11] and Yersinia pseudotuberculosis.
-Note: The  patterns  shown  above has been updated since their publication in
 [7].

-Expert(s) to contact by email:
           Hunter T.; hunter@salk-sc2.sdsc.edu
           Quinn A.M.; quinn@biomed.med.yale.edu

-Last update: April 2006 / Pattern revised.

""")

        self.assertEqual(len(record.references), 11)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Hanks S.K., Hunter T.")
        self.assertEqual(record.references[ 0].citation, """\
"Protein kinases 6. The eukaryotic protein kinase superfamily: kinase
(catalytic) domain structure and classification."
FASEB J. 9:576-596(1995).
PubMed=7768349""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Hunter T.")
        self.assertEqual(record.references[ 1].citation, """\
"Protein kinase classification."
Methods Enzymol. 200:3-37(1991).
PubMed=1835513""")
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Hanks S.K., Quinn A.M.")
        self.assertEqual(record.references[ 2].citation, """\
"Protein kinase catalytic domain sequence database: identification of
conserved features of primary structure and classification of family
members."
Methods Enzymol. 200:38-62(1991).
PubMed=1956325""")
        self.assertEqual(record.references[ 3].number, "4")
        self.assertEqual(record.references[ 3].authors, "Hanks S.K.")
        self.assertEqual(record.references[ 3].citation, 'Curr. Opin. Struct. Biol. 1:369-383(1991).')
        self.assertEqual(record.references[ 4].number, "5")
        self.assertEqual(record.references[ 4].authors, "Hanks S.K., Quinn A.M., Hunter T.")
        self.assertEqual(record.references[ 4].citation, """\
"The protein kinase family: conserved features and deduced phylogeny
of the catalytic domains."
Science 241:42-52(1988).
PubMed=3291115""")
        self.assertEqual(record.references[ 5].number, "6")
        self.assertEqual(record.references[ 5].authors, "Knighton D.R., Zheng J.H., Ten Eyck L.F., Ashford V.A., Xuong N.-H., Taylor S.S., Sowadski J.M.")
        self.assertEqual(record.references[ 5].citation, """\
"Crystal structure of the catalytic subunit of cyclic adenosine
monophosphate-dependent protein kinase."
Science 253:407-414(1991).
PubMed=1862342""")
        self.assertEqual(record.references[ 6].number, "7")
        self.assertEqual(record.references[ 6].authors, "Bairoch A., Claverie J.-M.")
        self.assertEqual(record.references[ 6].citation, """\
"Sequence patterns in protein kinases."
Nature 331:22-22(1988).
PubMed=3340146; DOI=10.1038/331022a0""")
        self.assertEqual(record.references[ 7].number, "8")
        self.assertEqual(record.references[ 7].authors, "Benner S.")
        self.assertEqual(record.references[ 7].citation, 'Nature 329:21-21(1987).')
        self.assertEqual(record.references[ 8].number, "9")
        self.assertEqual(record.references[ 8].authors, "Kirby R.")
        self.assertEqual(record.references[ 8].citation, """\
"Evolutionary origin of aminoglycoside phosphotransferase resistance
genes."
J. Mol. Evol. 30:489-492(1990).
PubMed=2165531""")
        self.assertEqual(record.references[ 9].number, "10")
        self.assertEqual(record.references[ 9].authors, "Littler E., Stuart A.D., Chee M.S.")
        self.assertEqual(record.references[ 9].citation, 'Nature 358:160-162(1992).')
        self.assertEqual(record.references[10].number, "11")
        self.assertEqual(record.references[10].authors, "Munoz-Dorado J., Inouye S., Inouye M.")
        self.assertEqual(record.references[10].citation, 'Cell 67:995-1006(1991).')

    def test_read_pdoc00113(self):
        "Reading Prodoc record PDOC00113"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00113.txt')
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00113")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS00123", "ALKALINE_PHOSPHATASE"))
        self.assertEqual(record.text, """\
************************************
* Alkaline phosphatase active site *
************************************

Alkaline phosphatase (EC 3.1.3.1) (ALP) [1] is a zinc and magnesium-containing
metalloenzyme  which hydrolyzes phosphate esters, optimally at high pH.  It is
found in nearly  all living organisms,  with the exception of some plants.  In
Escherichia coli, ALP (gene phoA) is found in the periplasmic space.  In yeast
it (gene  PHO8)  is  found  in  lysosome-like vacuoles and in mammals, it is a
glycoprotein attached to the membrane by a GPI-anchor.

In mammals, four different isozymes are currently known [2]. Three of them are
tissue-specific:  the  placental,  placental-like (germ cell)   and intestinal
isozymes.  The fourth form is  tissue non-specific and was previously known as
the liver/bone/kidney isozyme.

Streptomyces' species  involved  in  the  synthesis  of  streptomycin (SM), an
antibiotic, express  a  phosphatase (EC 3.1.3.39) (gene strK) which is  highly
related to ALP.   It specifically cleaves  both  streptomycin-6-phosphate and,
more slowly, streptomycin-3"-phosphate.

A serine is involved   in the catalytic activity of ALP. The region around the
active site serine is relatively well conserved and can be used as a signature
pattern.

-Consensus pattern: [IV]-x-D-S-[GAS]-[GASC]-[GAST]-[GA]-T
                    [S is the active site residue]
-Sequences known to belong to this class detected by the pattern: ALL.
-Other sequence(s) detected in Swiss-Prot: 3.
-Last update: June 1994 / Text revised.

""")

        self.assertEqual(len(record.references), 3)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Trowsdale J., Martin D., Bicknell D., Campbell I.")
        self.assertEqual(record.references[ 0].citation, """\
"Alkaline phosphatases."
Biochem. Soc. Trans. 18:178-180(1990).
PubMed=2379681""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Manes T., Glade K., Ziomek C.A., Millan J.L.")
        self.assertEqual(record.references[ 1].citation, """\
"Genomic structure and comparison of mouse tissue-specific alkaline
phosphatase genes."
Genomics 8:541-554(1990).
PubMed=2286375""")
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Mansouri K., Piepersberg W.")
        self.assertEqual(record.references[ 2].citation, """\
"Genetics of streptomycin production in Streptomyces griseus:
nucleotide sequence of five genes, strFGHIK, including a phosphatase
gene."
Mol. Gen. Genet. 228:459-469(1991).
PubMed=1654502""")

    def test_read_pdoc00144(self):
        "Reading Prodoc record PDOC00144"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00144.txt')
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00144")
        self.assertEqual(len(record.prosite_refs), 2)
        self.assertEqual(record.prosite_refs[0], ("PS00159", "ALDOLASE_KDPG_KHG_1"))
        self.assertEqual(record.prosite_refs[1], ("PS00160", "ALDOLASE_KDPG_KHG_2"))
        self.assertEqual(record.text, """\
*************************************************
* KDPG and KHG aldolases active site signatures *
*************************************************

4-hydroxy-2-oxoglutarate aldolase (EC 4.1.3.16)  (KHG-aldolase)  catalyzes the
interconversion of  4-hydroxy-2-oxoglutarate  into  pyruvate  and  glyoxylate.
Phospho-2-dehydro-3-deoxygluconate  aldolase   (EC 4.1.2.14)   (KDPG-aldolase)
catalyzes the interconversion of  6-phospho-2-dehydro-3-deoxy-D-gluconate into
pyruvate and glyceraldehyde 3-phosphate.

These two enzymes are structurally and functionally related [1]. They are both
homotrimeric proteins of approximately 220 amino-acid residues. They are class
I aldolases whose catalytic mechanism involves  the formation of a Schiff-base
intermediate  between  the  substrate  and the epsilon-amino group of a lysine
residue. In both enzymes, an arginine is required for catalytic activity.

We developed  two signature patterns for these enzymes. The first one contains
the active  site  arginine  and the second, the lysine involved in the Schiff-
base formation.

-Consensus pattern: G-[LIVM]-x(3)-E-[LIV]-T-[LF]-R
                    [R is the active site residue]
-Sequences known to belong to this class detected by the pattern: ALL,  except
 for Bacillus  subtilis  KDPG-aldolase  which  has  Thr  instead of Arg in the
 active site.
-Other sequence(s) detected in Swiss-Prot: NONE.

-Consensus pattern: G-x(3)-[LIVMF]-K-[LF]-F-P-[SA]-x(3)-G
                    [K is involved in Schiff-base formation]
-Sequences known to belong to this class detected by the pattern: ALL.
-Other sequence(s) detected in Swiss-Prot: NONE.

-Last update: November 1997 / Patterns and text revised.

""")

        self.assertEqual(len(record.references), 1)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Vlahos C.J., Dekker E.E.")
        self.assertEqual(record.references[ 0].citation, """\
"The complete amino acid sequence and identification of the
active-site arginine peptide of Escherichia coli
2-keto-4-hydroxyglutarate aldolase."
J. Biol. Chem. 263:11683-11691(1988).
PubMed=3136164""")

    def test_read_pdoc00149(self):
        "Reading Prodoc record PDOC00149"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00149.txt')
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00149")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS00165", "DEHYDRATASE_SER_THR"))
        self.assertEqual(record.text, """\
*********************************************************************
* Serine/threonine dehydratases pyridoxal-phosphate attachment site *
*********************************************************************

Serine and threonine  dehydratases [1,2]  are  functionally  and  structurally
related pyridoxal-phosphate dependent enzymes:

 - L-serine dehydratase (EC 4.3.1.17) and D-serine  dehydratase  (EC 4.3.1.18)
   catalyze the dehydratation of L-serine (respectively D-serine) into ammonia
   and pyruvate.
 - Threonine dehydratase  (EC 4.3.1.19) (TDH) catalyzes  the  dehydratation of
   threonine into  alpha-ketobutarate  and  ammonia.  In Escherichia coli  and
   other microorganisms,  two  classes  of  TDH  are  known  to  exist. One is
   involved in  the  biosynthesis of isoleucine, the other in hydroxamino acid
   catabolism.

Threonine synthase  (EC 4.2.3.1) is  also  a  pyridoxal-phosphate  enzyme,  it
catalyzes the  transformation of  homoserine-phosphate into threonine.  It has
been shown [3] that  threonine  synthase  is  distantly related to the serine/
threonine dehydratases.

In all these enzymes, the pyridoxal-phosphate group is  attached  to a  lysine
residue.  The sequence around  this residue is sufficiently conserved to allow
the derivation  of  a  pattern  specific  to serine/threonine dehydratases and
threonine synthases.

-Consensus pattern: [DESH]-x(4,5)-[STVG]-{EVKD}-[AS]-[FYI]-K-[DLIFSA]-[RLVMF]-
                    [GA]-[LIVMGA]
                    [The K is the pyridoxal-P attachment site]
-Sequences known to belong to this class detected by the pattern: ALL.
-Other sequence(s) detected in Swiss-Prot: 17.

-Note: Some   bacterial L-serine dehydratases - such as those from Escherichia
 coli - are iron-sulfur proteins [4] and do not belong to this family.

-Last update: December 2004 / Pattern and text revised.

""")

        self.assertEqual(len(record.references), 4)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Ogawa H., Gomi T., Konishi K., Date T., Nakashima H., Nose K., Matsuda Y., Peraino C., Pitot H.C., Fujioka M.")
        self.assertEqual(record.references[ 0].citation, """\
"Human liver serine dehydratase. cDNA cloning and sequence homology
with hydroxyamino acid dehydratases from other sources."
J. Biol. Chem. 264:15818-15823(1989).
PubMed=2674117""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Datta P., Goss T.J., Omnaas J.R., Patil R.V.")
        self.assertEqual(record.references[ 1].citation, """\
"Covalent structure of biodegradative threonine dehydratase of
Escherichia coli: homology with other dehydratases."
Proc. Natl. Acad. Sci. U.S.A. 84:393-397(1987).
PubMed=3540965""")
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Parsot C.")
        self.assertEqual(record.references[ 2].citation, """\
"Evolution of biosynthetic pathways: a common ancestor for threonine
synthase, threonine dehydratase and D-serine dehydratase."
EMBO J. 5:3013-3019(1986).
PubMed=3098560""")
        self.assertEqual(record.references[ 3].number, "4")
        self.assertEqual(record.references[ 3].authors, "Grabowski R., Hofmeister A.E.M., Buckel W.")
        self.assertEqual(record.references[ 3].citation, """\
"Bacterial L-serine dehydratases: a new family of enzymes containing
iron-sulfur clusters."
Trends Biochem. Sci. 18:297-300(1993).
PubMed=8236444""")

    def test_read_pdoc00340(self):
        "Reading Prodoc record PDOC00340"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00340.txt')
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00340")
        self.assertEqual(len(record.prosite_refs), 3)
        self.assertEqual(record.prosite_refs[0], ("PS00406", "ACTINS_1"))
        self.assertEqual(record.prosite_refs[1], ("PS00432", "ACTINS_2"))
        self.assertEqual(record.prosite_refs[2], ("PS01132", "ACTINS_ACT_LIKE"))
        self.assertEqual(record.text, """\
*********************
* Actins signatures *
*********************

Actins [1 to 4] are highly conserved contractile  proteins that are present in
all eukaryotic cells. In vertebrates there are three groups of actin isoforms:
alpha, beta and gamma.  The alpha actins are found in muscle tissues and are a
major constituent of the contractile apparatus.  The beta and gamma actins co-
exists in most cell  types as  components of the cytoskeleton and as mediators
of internal cell motility.  In plants [5]  there  are  many isoforms which are
probably involved  in  a  variety of  functions such as cytoplasmic streaming,
cell shape determination,  tip growth,  graviperception, cell wall deposition,
etc.

Actin exists either in a monomeric form (G-actin) or in a polymerized form (F-
actin). Each actin monomer  can  bind a molecule of ATP;  when  polymerization
occurs, the ATP is hydrolyzed.

Actin is a protein of from 374 to 379 amino acid  residues.  The  structure of
actin has been highly conserved in the course of evolution.

Recently some  divergent  actin-like  proteins have been identified in several
species. These proteins are:

 - Centractin  (actin-RPV)  from mammals, fungi (yeast ACT5, Neurospora crassa
   ro-4) and  Pneumocystis  carinii  (actin-II).  Centractin  seems  to  be  a
   component of  a  multi-subunit  centrosomal complex involved in microtubule
   based vesicle motility. This subfamily is also known as ARP1.
 - ARP2  subfamily  which  includes  chicken ACTL, yeast ACT2, Drosophila 14D,
   C.elegans actC.
 - ARP3  subfamily  which includes actin 2 from mammals, Drosophila 66B, yeast
   ACT4 and fission yeast act2.
 - ARP4  subfamily  which includes yeast ACT3 and Drosophila 13E.

We developed  three  signature  patterns. The first two are specific to actins
and span  positions  54 to 64 and 357 to 365. The last signature picks up both
actins and  the actin-like proteins and corresponds to positions 106 to 118 in
actins.

-Consensus pattern: [FY]-[LIV]-[GV]-[DE]-E-[ARV]-[QLAH]-x(1,2)-[RKQ](2)-[GD]
-Sequences known to belong to this class detected by the pattern: ALL,  except
 for the actin-like proteins and 10 actins.
-Other sequence(s) detected in Swiss-Prot: NONE.

-Consensus pattern: W-[IVC]-[STAK]-[RK]-x-[DE]-Y-[DNE]-[DE]
-Sequences known to belong to this class detected by the pattern: ALL,  except
 for the actin-like proteins and 9 actins.
-Other sequence(s) detected in Swiss-Prot: NONE.

-Consensus pattern: [LM]-[LIVMA]-T-E-[GAPQ]-x-[LIVMFYWHQPK]-[NS]-[PSTAQ]-x(2)-
                    N-[KR]
-Sequences known to belong to this class detected by the pattern: ALL,  except
 for 5 actins.
-Other sequence(s) detected in Swiss-Prot: NONE.

-Last update: December 2004 / Patterns and text revised.

""")

        self.assertEqual(len(record.references), 5)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Sheterline P., Clayton J., Sparrow J.C.")
        self.assertEqual(record.references[ 0].citation, '(In) Actins, 3rd Edition, Academic Press Ltd, London, (1996).')
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Pollard T.D., Cooper J.A.")
        self.assertEqual(record.references[ 1].citation, 'Annu. Rev. Biochem. 55:987-1036(1986).')
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Pollard T.D.")
        self.assertEqual(record.references[ 2].citation, """\
"Actin."
Curr. Opin. Cell Biol. 2:33-40(1990).
PubMed=2183841""")
        self.assertEqual(record.references[ 3].number, "4")
        self.assertEqual(record.references[ 3].authors, "Rubenstein P.A.")
        self.assertEqual(record.references[ 3].citation, """\
"The functional importance of multiple actin isoforms."
BioEssays 12:309-315(1990).
PubMed=2203335""")
        self.assertEqual(record.references[ 4].number, "5")
        self.assertEqual(record.references[ 4].authors, "Meagher R.B., McLean B.G.")
        self.assertEqual(record.references[ 4].citation, 'Cell Motil. Cytoskeleton 16:164-166(1990).')

    def test_read_pdoc00424(self):
        "Reading Prodoc record PDOC00424"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00424.txt',)
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00424")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS00488", "PAL_HISTIDASE"))
        self.assertEqual(record.text, """\
**********************************************************
* Phenylalanine and histidine ammonia-lyases active site *
**********************************************************

Phenylalanine ammonia-lyase (EC 4.3.1.5) (PAL) is  a  key  enzyme of plant and
fungi  phenylpropanoid  metabolism  which is involved in the biosynthesis of a
wide  variety  of secondary metabolites such  as  flavanoids,   furanocoumarin
phytoalexins and  cell  wall  components.  These compounds have many important
roles in plants during normal growth and in responses to environmental stress.
PAL catalyzes  the  removal  of  an  ammonia  group from phenylalanine to form
trans-cinnamate.

Histidine ammonia-lyase (EC 4.3.1.3) (histidase)  catalyzes  the first step in
histidine degradation, the removal of  an  ammonia  group  from  histidine  to
produce urocanic acid.

The two types of enzymes are functionally and  structurally related [1].  They
are the only enzymes  which are known to have the modified amino acid dehydro-
alanine (DHA) in their active site. A serine residue has been shown [2,3,4] to
be the  precursor  of  this  essential electrophilic moiety. The region around
this active  site  residue  is  well  conserved and can be used as a signature
pattern.

-Consensus pattern: [GS]-[STG]-[LIVM]-[STG]-[SAC]-S-G-[DH]-L-x-[PN]-L-[SA]-
                    x(2,3)-[SAGVTL]
                    [S is the active site residue]
-Sequences known to belong to this class detected by the pattern: ALL.
-Other sequence(s) detected in Swiss-Prot: NONE.
-Last update: April 2006 / Pattern revised.

""")

        self.assertEqual(len(record.references), 4)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Taylor R.G., Lambert M.A., Sexsmith E., Sadler S.J., Ray P.N., Mahuran D.J., McInnes R.R.")
        self.assertEqual(record.references[ 0].citation, """\
"Cloning and expression of rat histidase. Homology to two bacterial
histidases and four phenylalanine ammonia-lyases."
J. Biol. Chem. 265:18192-18199(1990).
PubMed=2120224""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Langer M., Reck G., Reed J., Retey J.")
        self.assertEqual(record.references[ 1].citation, """\
"Identification of serine-143 as the most likely precursor of
dehydroalanine in the active site of histidine ammonia-lyase. A study
of the overexpressed enzyme by site-directed mutagenesis."
Biochemistry 33:6462-6467(1994).
PubMed=8204579""")
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Schuster B., Retey J.")
        self.assertEqual(record.references[ 2].citation, """\
"Serine-202 is the putative precursor of the active site
dehydroalanine of phenylalanine ammonia lyase. Site-directed
mutagenesis studies on the enzyme from parsley (Petroselinum crispum
L.)."
FEBS Lett. 349:252-254(1994).
PubMed=8050576""")
        self.assertEqual(record.references[ 3].number, "4")
        self.assertEqual(record.references[ 3].authors, "Taylor R.G., McInnes R.R.")
        self.assertEqual(record.references[ 3].citation, """\
"Site-directed mutagenesis of conserved serines in rat histidase.
Identification of serine 254 as an essential active site residue."
J. Biol. Chem. 269:27473-27477(1994).
PubMed=7961661""")

    def test_read_pdoc00472(self):
        "Reading Prodoc record PDOC00472"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00472.txt')
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00472")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS00546", "CYSTEINE_SWITCH"))
        self.assertEqual(record.text, """\
*****************************
* Matrixins cysteine switch *
*****************************

Mammalian extracellular matrix metalloproteinases (EC 3.4.24.-), also known as
matrixins [1] (see <PDOC00129>), are zinc-dependent enzymes. They are secreted
by cells  in an inactive form (zymogen) that differs from the mature enzyme by
the presence  of  an  N-terminal propeptide. A highly conserved octapeptide is
found two  residues  downstream  of the C-terminal end of the propeptide. This
region has been shown to be  involved  in  autoinhibition  of matrixins [2,3];
a cysteine  within the octapeptide chelates  the  active  site  zinc ion, thus
inhibiting the  enzyme.  This  region has been called the 'cysteine switch' or
'autoinhibitor region'.

A cysteine switch has been found in the following zinc proteases:

 - MMP-1 (EC 3.4.24.7) (interstitial collagenase).
 - MMP-2 (EC 3.4.24.24) (72 Kd gelatinase).
 - MMP-3 (EC 3.4.24.17) (stromelysin-1).
 - MMP-7 (EC 3.4.24.23) (matrilysin).
 - MMP-8 (EC 3.4.24.34) (neutrophil collagenase).
 - MMP-9 (EC 3.4.24.35) (92 Kd gelatinase).
 - MMP-10 (EC 3.4.24.22) (stromelysin-2).
 - MMP-11 (EC 3.4.24.-) (stromelysin-3).
 - MMP-12 (EC 3.4.24.65) (macrophage metalloelastase).
 - MMP-13 (EC 3.4.24.-) (collagenase 3).
 - MMP-14 (EC 3.4.24.-) (membrane-type matrix metalliproteinase 1).
 - MMP-15 (EC 3.4.24.-) (membrane-type matrix metalliproteinase 2).
 - MMP-16 (EC 3.4.24.-) (membrane-type matrix metalliproteinase 3).
 - Sea urchin hatching enzyme (EC 3.4.24.12) (envelysin) [4].
 - Chlamydomonas reinhardtii gamete lytic enzyme (GLE) [5].

-Consensus pattern: P-R-C-[GN]-x-P-[DR]-[LIVSAPKQ]
                    [C chelates the zinc ion]
-Sequences known to belong to this class detected by the pattern: ALL,  except
 for cat MMP-7 and mouse MMP-11.
-Other sequence(s) detected in Swiss-Prot: NONE.
-Last update: November 1997 / Pattern and text revised.

""")

        self.assertEqual(len(record.references), 5)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Woessner J.F. Jr.")
        self.assertEqual(record.references[ 0].citation, """\
"Matrix metalloproteinases and their inhibitors in connective tissue
remodeling."
FASEB J. 5:2145-2154(1991).
PubMed=1850705""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Sanchez-Lopez R., Nicholson R., Gesnel M.C., Matrisian L.M., Breathnach R.")
        self.assertEqual(record.references[ 1].citation, 'J. Biol. Chem. 263:11892-11899(1988).')
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Park A.J., Matrisian L.M., Kells A.F., Pearson R., Yuan Z.Y., Navre M.")
        self.assertEqual(record.references[ 2].citation, """\
"Mutational analysis of the transin (rat stromelysin) autoinhibitor
region demonstrates a role for residues surrounding the 'cysteine
switch'."
J. Biol. Chem. 266:1584-1590(1991).
PubMed=1988438""")
        self.assertEqual(record.references[ 3].number, "4")
        self.assertEqual(record.references[ 3].authors, "Lepage T., Gache C.")
        self.assertEqual(record.references[ 3].citation, """\
"Early expression of a collagenase-like hatching enzyme gene in the
sea urchin embryo."
EMBO J. 9:3003-3012(1990).
PubMed=2167841""")
        self.assertEqual(record.references[ 4].number, "5")
        self.assertEqual(record.references[ 4].authors, "Kinoshita T., Fukuzawa H., Shimada T., Saito T., Matsuda Y.")
        self.assertEqual(record.references[ 4].citation, """\
"Primary structure and expression of a gamete lytic enzyme in
Chlamydomonas reinhardtii: similarity of functional domains to matrix
metalloproteases."
Proc. Natl. Acad. Sci. U.S.A. 89:4693-4697(1992).
PubMed=1584806""")

    def test_read_pdoc00640(self):
        "Reading Prodoc record PDOC00640"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00640.txt',)
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00640")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS00812", "GLYCOSYL_HYDROL_F8"))
        self.assertEqual(record.text, """\
******************************************
* Glycosyl hydrolases family 8 signature *
******************************************

The microbial degradation  of cellulose and  xylans requires  several types of
enzymes such as endoglucanases (EC 3.2.1.4),  cellobiohydrolases (EC 3.2.1.91)
(exoglucanases), or xylanases (EC 3.2.1.8) [1,2].  Fungi and bacteria produces
a spectrum of cellulolytic  enzymes (cellulases)  and  xylanases which, on the
basis of sequence similarities,  can be classified into families. One of these
families is known as the cellulase family D [3] or as  the glycosyl hydrolases
family 8  [4,E1].  The  enzymes  which  are  currently known to belong to this
family are listed below.

 - Acetobacter xylinum endonuclease cmcAX.
 - Bacillus strain KSM-330 acidic endonuclease K (Endo-K).
 - Cellulomonas josui endoglucanase 2 (celB).
 - Cellulomonas uda endoglucanase.
 - Clostridium cellulolyticum endoglucanases C (celcCC).
 - Clostridium thermocellum endoglucanases A (celA).
 - Erwinia chrysanthemi minor endoglucanase y (celY).
 - Bacillus circulans beta-glucanase (EC 3.2.1.73).
 - Escherichia coli hypothetical protein yhjM.

The most conserved region in  these enzymes is  a stretch of about 20 residues
that contains  two conserved aspartate. The first asparatate is thought [5] to
act as the nucleophile in the catalytic mechanism. We have used this region as
a signature pattern.

-Consensus pattern: A-[ST]-D-[AG]-D-x(2)-[IM]-A-x-[SA]-[LIVM]-[LIVMG]-x-A-
                    x(3)-[FW]
                    [The first D is an active site residue]
-Sequences known to belong to this class detected by the pattern: ALL.
-Other sequence(s) detected in Swiss-Prot: NONE.

-Expert(s) to contact by email:
           Henrissat B.; bernie@afmb.cnrs-mrs.fr

-Last update: November 1997 / Text revised.

""")

        self.assertEqual(len(record.references), 6)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Beguin P.")
        self.assertEqual(record.references[ 0].citation, """\
"Molecular biology of cellulose degradation."
Annu. Rev. Microbiol. 44:219-248(1990).
PubMed=2252383; DOI=10.1146/annurev.mi.44.100190.001251""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Gilkes N.R., Henrissat B., Kilburn D.G., Miller R.C. Jr., Warren R.A.J.")
        self.assertEqual(record.references[ 1].citation, """\
"Domains in microbial beta-1, 4-glycanases: sequence conservation,
function, and enzyme families."
Microbiol. Rev. 55:303-315(1991).
PubMed=1886523""")
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Henrissat B., Claeyssens M., Tomme P., Lemesle L., Mornon J.-P.")
        self.assertEqual(record.references[ 2].citation, """\
"Cellulase families revealed by hydrophobic cluster analysis."
Gene 81:83-95(1989).
PubMed=2806912""")
        self.assertEqual(record.references[ 3].number, "4")
        self.assertEqual(record.references[ 3].authors, "Henrissat B.")
        self.assertEqual(record.references[ 3].citation, """\
"A classification of glycosyl hydrolases based on amino acid sequence
similarities."
Biochem. J. 280:309-316(1991).
PubMed=1747104""")
        self.assertEqual(record.references[ 4].number, "5")
        self.assertEqual(record.references[ 4].authors, "Alzari P.M., Souchon H., Dominguez R.")
        self.assertEqual(record.references[ 4].citation, """\
"The crystal structure of endoglucanase CelA, a family 8 glycosyl
hydrolase from Clostridium thermocellum."
Structure 4:265-275(1996).
PubMed=8805535""")
        self.assertEqual(record.references[ 5].number, "E1")
        self.assertEqual(record.references[ 5].authors, "")
        self.assertEqual(record.references[ 5].citation, 'http://www.expasy.org/cgi-bin/lists?glycosid.txt')

    def test_read_pdoc00787(self):
        "Reading Prodoc record PDOC00787"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00787.txt')
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00787")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS01027", "GLYCOSYL_HYDROL_F39"))
        self.assertEqual(record.text, """\
******************************************************
* Glycosyl hydrolases family 39 putative active site *
******************************************************

It has  been  shown  [1,E1]  that  the  following  glycosyl  hydrolases can be
classified into a single family on the basis of sequence similarities:

 - Mammalian lysosomal alpha-L-iduronidase (EC 3.2.1.76).
 - Caldocellum  saccharolyticum  and  Thermoanaerobacter saccharolyticum beta-
   xylosidase (EC 3.2.1.37) (gene xynB).

The best  conserved  regions  in  these  enzymes is  located in the N-terminal
section. It   contains  a  glutamic  acid  residue  which,  on  the  basis  of
similarities with other  families of glycosyl hydrolases [2], probably acts as
the proton donor in the catalytic mechanism. We use this region as a signature
pattern.

-Consensus pattern: W-x-F-E-x-W-N-E-P-[DN]
                    [The second E may be the active site residue]
-Sequences known to belong to this class detected by the pattern: ALL.
-Other sequence(s) detected in Swiss-Prot: NONE.

-Expert(s) to contact by email:
           Henrissat B.; bernie@afmb.cnrs-mrs.fr

-Last update: May 2004 / Text revised.

""")

        self.assertEqual(len(record.references), 3)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Henrissat B., Bairoch A.")
        self.assertEqual(record.references[ 0].citation, """\
"New families in the classification of glycosyl hydrolases based on
amino acid sequence similarities."
Biochem. J. 293:781-788(1993).
PubMed=8352747""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Henrissat B., Callebaut I., Fabrega S., Lehn P., Mornon J.-P., Davies G.")
        self.assertEqual(record.references[ 1].citation, """\
"Conserved catalytic machinery and the prediction of a common fold for
several families of glycosyl hydrolases."
Proc. Natl. Acad. Sci. U.S.A. 92:7090-7094(1995).
PubMed=7624375""")
        self.assertEqual(record.references[ 2].number, "E1")
        self.assertEqual(record.references[ 2].authors, '')
        self.assertEqual(record.references[ 2].citation, "http://www.expasy.org/cgi-bin/lists?glycosid.txt")

    def test_read_pdoc0933(self):
        "Reading Prodoc record PDOC00933"
        filename = os.path.join( 'Prosite', 'Doc', 'pdoc00933.txt')
        handle = open(filename)
        record = Prodoc.read(handle)
        handle.close()

        self.assertEqual(record.accession, "PDOC00933")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS01213", "GLOBIN_FAM_2"))
        self.assertEqual(record.text, """\
**********************************************
* Protozoan/cyanobacterial globins signature *
**********************************************

Globins are heme-containing  proteins involved in  binding and/or transporting
oxygen [1]. Almost all globins belong to a large family (see <PDOC00793>), the
only exceptions  are  the  following proteins which form a family of their own
[2,3,4]:

 - Monomeric  hemoglobins  from the protozoan Paramecium caudatum, Tetrahymena
   pyriformis and Tetrahymena thermophila.
 - Cyanoglobins  from  the  cyanobacteria Nostoc commune and Synechocystis PCC
   6803.
 - Globins  LI637  and  LI410  from  the chloroplast of the alga Chlamydomonas
   eugametos.
 - Mycobacterium tuberculosis globins glbN and glbO.

These proteins  contain a conserved histidine which could be involved in heme-
binding. As a signature pattern, we use a conserved region that ends with this
residue.

-Consensus pattern: F-[LF]-x(4)-[GE]-G-[PAT]-x(2)-[YW]-x-[GSE]-[KRQAE]-x(1,5)-
                    [LIVM]-x(3)-H
                    [The H may be a heme ligand]
-Sequences known to belong to this class detected by the pattern: ALL.
-Other sequence(s) detected in Swiss-Prot: NONE.
-Last update: April 2006 / Pattern revised.

""")

        self.assertEqual(len(record.references), 4)
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Concise Encyclopedia Biochemistry, Second Edition, Walter de Gruyter, Berlin New-York (1988).")
        self.assertEqual(record.references[ 0].citation, '')
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Takagi T.")
        self.assertEqual(record.references[ 1].citation, 'Curr. Opin. Struct. Biol. 3:413-418(1993).')
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Couture M., Chamberland H., St-Pierre B., Lafontaine J., Guertin M.")
        self.assertEqual(record.references[ 2].citation, """\
"Nuclear genes encoding chloroplast hemoglobins in the unicellular
green alga Chlamydomonas eugametos."
Mol. Gen. Genet. 243:185-197(1994).
PubMed=8177215""")
        self.assertEqual(record.references[ 3].number, "4")
        self.assertEqual(record.references[ 3].authors, "Couture M., Das T.K., Savard P.Y., Ouellet Y., Wittenberg J.B., Wittenberg B.A., Rousseau D.L., Guertin M.")
        self.assertEqual(record.references[ 3].citation, """\
"Structural investigations of the hemoglobin of the cyanobacterium
Synechocystis PCC6803 reveal a unique distal heme pocket."
Eur. J. Biochem. 267:4770-4780(2000).
PubMed=10903511""")

class TestProdocParse(unittest.TestCase):

    def test_parse_pdoc(self):
        "Parsing an excerpt of prosite.doc" 
        filename = os.path.join( 'Prosite', 'Doc', 'prosite.excerpt.doc')
        handle = open(filename)
        records = Prodoc.parse(handle)

        # Testing the first parsed record
        record = records.next()
        self.assertEqual(record.accession, "PDOC00000")
        self.assertEqual(len(record.prosite_refs), 0)
        self.assertEqual(record.text, """\
**********************************
*** PROSITE documentation file ***
**********************************

Release 20.43 of 10-Feb-2009.

PROSITE is developed by the Swiss Institute of Bioinformatics (SIB) under
the responsability of Amos Bairoch and Nicolas Hulo.

This release was prepared by: Nicolas Hulo, Virginie Bulliard, Petra
Langendijk-Genevaux and Christian Sigrist with the help of Edouard
de Castro, Lorenzo Cerutti, Corinne Lachaize and Amos Bairoch.


See: http://www.expasy.org/prosite/
Email: prosite@expasy.org

Acknowledgements:

 - To all those mentioned in this document who have reviewed the entry(ies)
   for which they are listed as experts. With specific thanks to Rein Aasland,
   Mark Boguski, Peer Bork, Josh Cherry, Andre Chollet, Frank Kolakowski,
   David Landsman, Bernard Henrissat, Eugene Koonin, Steve Henikoff, Manuel
   Peitsch and Jonathan Reizer.
 - Jim Apostolopoulos is the author of the PDOC00699 entry.
 - Brigitte Boeckmann is the author of the PDOC00691, PDOC00703, PDOC00829,
   PDOC00796, PDOC00798, PDOC00799, PDOC00906, PDOC00907, PDOC00908,
   PDOC00912, PDOC00913, PDOC00924, PDOC00928, PDOC00929, PDOC00955,
   PDOC00961, PDOC00966, PDOC00988 and PDOC50020 entries.
 - Jean-Louis Boulay is the author of the PDOC01051, PDOC01050, PDOC01052,
   PDOC01053 and PDOC01054 entries.
 - Ryszard Brzezinski is the author of the PDOC60000 entry.
 - Elisabeth Coudert is the author of the PDOC00373 entry.
 - Kirill Degtyarenko is the author of the PDOC60001 entry.
 - Christian Doerig is the author of the PDOC01049 entry.
 - Kay Hofmann is the author of the PDOC50003, PDOC50006, PDOC50007 and
   PDOC50017 entries.
 - Chantal Hulo is the author of the PDOC00987 entry.
 - Karine Michoud is the author of the PDOC01044 and PDOC01042 entries.
 - Yuri Panchin is the author of the PDOC51013 entry.
 - S. Ramakumar is the author of the PDOC51052, PDOC60004, PDOC60010,
   PDOC60011, PDOC60015, PDOC60016, PDOC60018, PDOC60020, PDOC60021,
   PDOC60022, PDOC60023, PDOC60024, PDOC60025, PDOC60026, PDOC60027,
   PDOC60028, PDOC60029 and PDOC60030 entries.
 - Keith Robison is the author of the PDOC00830 and PDOC00861 entries.

   ------------------------------------------------------------------------
   PROSITE is copyright.   It  is  produced  by  the  Swiss  Institute   of
   Bioinformatics (SIB). There are no restrictions on its use by non-profit
   institutions as long as its  content is in no way modified. Usage by and
   for commercial  entities requires a license agreement.   For information
   about  the  licensing  scheme   send  an  email to license@isb-sib.ch or
   see: http://www.expasy.org/prosite/prosite_license.htm.
   ------------------------------------------------------------------------

""")

        # Testing the second parsed record"
        record = records.next()
        self.assertEqual(record.accession, "PDOC00001")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS00001", "ASN_GLYCOSYLATION"))
        self.assertEqual(record.text, """\
************************
* N-glycosylation site *
************************

It has been known for a long time [1] that potential N-glycosylation sites are
specific to the consensus sequence Asn-Xaa-Ser/Thr.  It must be noted that the
presence of the consensus  tripeptide  is  not sufficient  to conclude that an
asparagine residue is glycosylated, due to  the fact that the  folding of  the
protein plays an important  role in the  regulation of N-glycosylation [2]. It
has been shown [3] that  the  presence of proline between Asn and Ser/Thr will
inhibit N-glycosylation; this  has  been confirmed by a recent [4] statistical
analysis of glycosylation sites, which also  shows that about 50% of the sites
that have a proline C-terminal to Ser/Thr are not glycosylated.

It must also  be noted that there  are  a few  reported cases of glycosylation
sites with the pattern Asn-Xaa-Cys; an  experimentally demonstrated occurrence
of such a non-standard site is found in the plasma protein C [5].

-Consensus pattern: N-{P}-[ST]-{P}
                    [N is the glycosylation site]
-Last update: May 1991 / Text revised.

""")
        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Marshall R.D.")
        self.assertEqual(record.references[ 0].citation, """\
"Glycoproteins."
Annu. Rev. Biochem. 41:673-702(1972).
PubMed=4563441; DOI=10.1146/annurev.bi.41.070172.003325""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Pless D.D., Lennarz W.J.")
        self.assertEqual(record.references[ 1].citation, """\
"Enzymatic conversion of proteins to glycoproteins."
Proc. Natl. Acad. Sci. U.S.A. 74:134-138(1977).
PubMed=264667""")
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Bause E.")
        self.assertEqual(record.references[ 2].citation, """\
"Structural requirements of N-glycosylation of proteins. Studies with
proline peptides as conformational probes."
Biochem. J. 209:331-336(1983).
PubMed=6847620""")
        self.assertEqual(record.references[ 3].number, "4")
        self.assertEqual(record.references[ 3].authors, "Gavel Y., von Heijne G.")
        self.assertEqual(record.references[ 3].citation, """\
"Sequence differences between glycosylated and non-glycosylated
Asn-X-Thr/Ser acceptor sites: implications for protein engineering."
Protein Eng. 3:433-442(1990).
PubMed=2349213""")
        self.assertEqual(record.references[ 4].number, "5")
        self.assertEqual(record.references[ 4].authors, "Miletich J.P., Broze G.J. Jr.")
        self.assertEqual(record.references[ 4].citation, """\
"Beta protein C is not glycosylated at asparagine 329. The rate of
translation may influence the frequency of usage at
asparagine-X-cysteine sites."
J. Biol. Chem. 265:11397-11404(1990).
PubMed=1694179""")

        # Testing the third parsed record" 
        record = records.next()
        self.assertEqual(record.accession, "PDOC00004")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS00004", "CAMP_PHOSPHO_SITE"))
        self.assertEqual(record.text, """\
****************************************************************
* cAMP- and cGMP-dependent protein kinase phosphorylation site *
****************************************************************

There has been a  number of studies  relative to the  specificity of cAMP- and
cGMP-dependent protein kinases [1,2,3].  Both types of kinases appear to share
a preference  for  the  phosphorylation  of serine or threonine residues found
close to at least  two consecutive N-terminal  basic residues. It is important
to note that there are quite a number of exceptions to this rule.

-Consensus pattern: [RK](2)-x-[ST]
                    [S or T is the phosphorylation site]
-Last update: June 1988 / First entry.

""")

        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Fremisco J.R., Glass D.B., Krebs E.G.")
        self.assertEqual(record.references[ 0].citation, """\
J. Biol. Chem. 255:4240-4245(1980).""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Glass D.B., Smith S.B.")
        self.assertEqual(record.references[ 1].citation, """\
"Phosphorylation by cyclic GMP-dependent protein kinase of a synthetic
peptide corresponding to the autophosphorylation site in the enzyme."
J. Biol. Chem. 258:14797-14803(1983).
PubMed=6317673""")
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Glass D.B., el-Maghrabi M.R., Pilkis S.J.")
        self.assertEqual(record.references[ 2].citation, """\
"Synthetic peptides corresponding to the site phosphorylated in
6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase as substrates of
cyclic nucleotide-dependent protein kinases."
J. Biol. Chem. 261:2987-2993(1986).
PubMed=3005275""")

        # Testing the fourth parsed record"
        record = records.next()
        self.assertEqual(record.accession, "PDOC60030")
        self.assertEqual(len(record.prosite_refs), 1)
        self.assertEqual(record.prosite_refs[0], ("PS60030", "BACTERIOCIN_IIA"))
        self.assertEqual(record.text, """\
******************************************
* Bacteriocin class IIa family signature *
******************************************

Many Gram-positive  bacteria  produce  ribosomally  synthesized  antimicrobial
peptides, often  termed  bacteriocins. One important and well studied class of
bacteriocins is the class IIa or pediocin-like bacteriocins produced by lactic
acid bacteria.  All  class  IIa  bacteriocins  are produced by food-associated
strains, isolated  from  a  variety of food products of industrial and natural
origins, including  meat  products,  dairy  products and vegetables. Class IIa
bacteriocins are all cationic, display anti-Listeria activity, and kill target
cells by permeabilizing the cell membrane [1-3].

Class IIa  bacteriocins  contain  between  37  and 48 residues. Based on their
primary structures,  the  peptide  chains  of  class  IIa  bacteriocins may be
divided roughly into two regions: a hydrophilic, cationic and highly conserved
N-terminal region,  and  a  less  conserved hydrophobic/amphiphilic C-terminal
region. The  N-terminal  region  contains  the conserved Y-G-N-G-V/L 'pediocin
box' motif  and  two conserved cysteine residues joined by a disulfide bridge.
It forms  a  three-stranded antiparallel beta-sheet supported by the conserved
disulfide bridge  (see <PDB:1OG7>). This cationic N-terminal beta-sheet domain
mediates binding of the class IIa bacteriocin to the target cell membrane. The
C-terminal region forms a hairpin-like domain (see <PDB:1OG7>) that penetrates
into the  hydrophobic  part  of  the  target  cell membrane, thereby mediating
leakage through  the  membrane.  The  two domains are joined by a hinge, which
enables movement of the domains relative to each other [2,3].

Some proteins  known  to belong to the class IIa bacteriocin family are listed
below:

 - Pediococcus acidilactici pediocin PA-1.
 - Leuconostoc mesenteroides mesentericin Y105.
 - Carnobacterium piscicola carnobacteriocin B2.
 - Lactobacillus sake sakacin P.
 - Enterococcus faecium enterocin A.
 - Enterococcus faecium enterocin P.
 - Leuconostoc gelidum leucocin A.
 - Lactobacillus curvatus curvacin A.
 - Listeria innocua listeriocin 743A.

The pattern  we  developed  for  the  class  IIa bacteriocin family covers the
'pediocin box' motif.

-Conserved pattern: Y-G-N-G-[VL]-x-C-x(4)-C
-Sequences known to belong to this class detected by the pattern: ALL.
-Other sequence(s) detected in Swiss-Prot: NONE.

-Expert(s) to contact by email:
           Ramakumar S.; ramak@physics.iisc.ernet.in

-Last update: March 2006 / First entry.

""")

        self.assertEqual(record.references[ 0].number, "1")
        self.assertEqual(record.references[ 0].authors, "Ennahar S., Sonomoto K., Ishizaki A.")
        self.assertEqual(record.references[ 0].citation, """\
"Class IIa bacteriocins from lactic acid bacteria: antibacterial
activity and food preservation."
J. Biosci. Bioeng. 87:705-716(1999).
PubMed=16232543""")
        self.assertEqual(record.references[ 1].number, "2")
        self.assertEqual(record.references[ 1].authors, "Johnsen L., Fimland G., Nissen-Meyer J.")
        self.assertEqual(record.references[ 1].citation, """\
"The C-terminal domain of pediocin-like antimicrobial peptides (class
IIa bacteriocins) is involved in specific recognition of the
C-terminal part of cognate immunity proteins and in determining the
antimicrobial spectrum."
J. Biol. Chem. 280:9243-9250(2005).
PubMed=15611086; DOI=10.1074/jbc.M412712200""")
        self.assertEqual(record.references[ 2].number, "3")
        self.assertEqual(record.references[ 2].authors, "Fimland G., Johnsen L., Dalhus B., Nissen-Meyer J.")
        self.assertEqual(record.references[ 2].citation, """\
"Pediocin-like antimicrobial peptides (class IIa bacteriocins) and
their immunity proteins: biosynthesis, structure, and mode of
action."
J. Pept. Sci. 11:688-696(2005).
PubMed=16059970; DOI=10.1002/psc.699""")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
