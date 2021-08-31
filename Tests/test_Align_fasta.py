# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.fasta module."""
import unittest

from Bio.Align.fasta_m8 import AlignmentIterator

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.emboss."
    ) from None


class TestFasta(unittest.TestCase):

    def test_m8CB(self):
        # Alignment file obtained by running
        # fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg
        # in the fasta36 source distribution
        path = "Fasta/output_m8CB.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.commandline, "fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        alignment = alignments[0]
        return
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 112)
        self.assertEqual(alignment.annotations["similarity"], 112)
        self.assertEqual(alignment.annotations["gaps"], 19)
        self.assertAlmostEqual(alignment.annotations["score"], 591.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||",
        )

    def test_m8CC(self):
        # Alignment file obtained by running
        # fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg
        # in the fasta36 source distribution
        path = "Fasta/output_m8CC.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.commandline, "fasta36 -q -m 8CC seq/mgstm1.aa seq/prot_test.lseg")
            alignments = list(alignments)
        return
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 112)
        self.assertEqual(alignment.annotations["similarity"], 112)
        self.assertEqual(alignment.annotations["gaps"], 19)
        self.assertAlmostEqual(alignment.annotations["score"], 591.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)

fasta_alignments = """\
# fasta36 -q seq/mgstm1.aa seq/prot_test.lseg
FASTA searches a protein or DNA sequence data bank
 version 36.3.8h May, 2020
Please cite:
 W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

Query: seq/mgstm1.aa
  1>>>sp|P10649|GSTM1_MOUSE Glutathione S-transferase Mu 1; GST 1-1; GST class-mu 1; Glutathione S-transferase GT8.7; pmGT10 - 218 aa
Library: seq/prot_test.lseg
     2267 residues in    12 sequences

Statistics: (shuffled [472]) MLE statistics: Lambda= 0.1773;  K=0.007457
 statistics sampled from 4 (4) to 472 sequences
Algorithm: FASTA (3.8 Nov 2011) [optimized]
Parameters: BL50 matrix (15:-5), open/ext: -10/-2
 ktup: 2, E-join: 1 (0.667), E-opt: 0.2 (0.333), width:  16
 Scan time:  0.010

The best scores are:                                      opt bits E(12)
sp|P09488|GSTM1_HUMAN Glutathione S-transferase Mu ( 218) 1242 326.7 2.6e-93
sp|P00502|GSTA1_RAT Glutathione S-transferase alph ( 222)  237 69.6 6.5e-16
sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Ho ( 142)   51 22.0   0.088
sp|P00517|KAPCA_BOVIN cAMP-dependent protein kinas ( 351)   43 20.0    0.86
sp|P14960|RBS_GUITH Ribulose bisphosphate carboxyl ( 139)   36 18.2     1.2
sp|P01593|KV101_HUMAN Ig kappa chain V-I region AG ( 108)   31 16.9     2.1
sp|P99998|CYC_PANTR Cytochrome c OS=Pan troglodyte ( 105)   30 16.6     2.4
sp|P02585|TNNC2_HUMAN Troponin C, skeletal muscle  ( 160)   30 16.6     3.5
sp|P03435|HEMA_I75A3 Hemagglutinin OS=Influenza A  ( 567)   37 18.4     3.5
sp|P60615|NXL1A_BUNMU Alpha-bungarotoxin isoform A (  95)   26 15.6       4
sp|P00193|FER_PEPAS Ferredoxin OS=Peptostreptococc (  54)   22 14.6     4.5
sp|P01834|IGKC_HUMAN Ig kappa chain C region OS=Ho ( 106)   23 14.9     6.5

>>sp|P09488|GSTM1_HUMAN Glutathione S-transferase Mu 1 O  (218 aa)
 initn: 1242 init1: 1242 opt: 1242  Z-score: 1727.0  bits: 326.7 E(12): 2.6e-93
Smith-Waterman score: 1242; 78.0% identity (95.4% similar) in 218 aa overlap (1-218:1-218)

               10        20        30        40        50        60
sp|P10 MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNL
       ::::::::..:::.: ::.:::::::::.::.::::::::.:::::::::::::::::::
sp|P09 MPMILGYWDIRGLAHAIRLLLEYTDSSYEEKKYTMGDAPDYDRSQWLNEKFKLGLDFPNL
               10        20        30        40        50        60

               70        80        90       100       110       120
sp|P10 PYLIDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDF
       ::::::.:::::::::: :.::::.: ::::::.::.::.:::.::..::: :.::::.:
sp|P09 PYLIDGAHKITQSNAILCYIARKHNLCGETEEEKIRVDILENQTMDNHMQLGMICYNPEF
               70        80        90       100       110       120

              130       140       150       160       170       180
sp|P10 EKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPN
       :: ::..:. .:::.:::::::::::::::.:.:.::::.::.:: .:.:::::::::::
sp|P09 EKLKPKYLEELPEKLKLYSEFLGKRPWFAGNKITFVDFLVYDVLDLHRIFEPKCLDAFPN
              130       140       150       160       170       180

              190       200       210        
sp|P10 LRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK
       :.::..:::::.::::::::::..  :.::::: :.::
sp|P09 LKDFISRFEGLEKISAYMKSSRFLPRPVFSKMAVWGNK
              190       200       210        

>>sp|P00502|GSTA1_RAT Glutathione S-transferase alpha-1   (222 aa)
 initn: 204 init1:  73 opt: 237  Z-score: 337.5  bits: 69.6 E(12): 6.5e-16
Smith-Waterman score: 237; 27.4% identity (57.0% similar) in 223 aa overlap (4-218:6-218)

                 10        20        30        40        50        
sp|P10   MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKL--GLD
            .: :.:.::  . :: ::  .   .:::         : .:    ::.:   .: 
sp|P00 MSGKPVLHYFNARGRMECIRWLLAAAGVEFDEK---------FIQSPEDLEKLKKDGNLM
               10        20        30                 40        50 

         60         70        80        90       100        110    
sp|P10 FPNLPYL-IDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMD-TRMQLIML
       : ..:.. :::  :..:. ::: :.: :. : :.  .::   :.  . ..: :.: . ..
sp|P00 FDQVPMVEIDG-MKLAQTRAILNYIATKYDLYGKDMKERALIDMYTEGILDLTEMIMQLV
              60         70        80        90       100       110

          120         130       140         150       160       170
sp|P10 CYNPDFEKQKPEFLK--TIPEKMKLYSEFLGK--RPWFAGDKVTYVDFLAYDILDQYRMF
          :: .. :  . :  :  . .  . . : .  . ...:...: ::.   ..:   . :
sp|P00 ICPPDQKEAKTALAKDRTKNRYLPAFEKVLKSHGQDYLVGNRLTRVDIHLLELLLYVEEF
              120       130       140       150       160       170

              180       190       200       210            
sp|P10 EPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK    
       . . : .:: :. : .:. .: ... ... .     :. .:. . . :    
sp|P00 DASLLTSFPLLKAFKSRISSLPNVKKFLQPGSQRKLPVDAKQIEEARKIFKF
              180       190       200       210       220  

>>sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo s  (142 aa)
 initn:  40 init1:  40 opt:  51  Z-score: 83.9  bits: 22.0 E(12): 0.088
Smith-Waterman score: 51; 25.6% identity (69.2% similar) in 39 aa overlap (177-214:36-73)

        150       160       170       180       190       200      
sp|P10 WFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKS-SRYIA
                                     .::. . .. .:. :.. :: .:. .. .:
sp|P69 ADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFD-LSHGSAQVKGHGKKVA
          10        20        30        40         50        60    

         210                                                       
sp|P10 TPIFSKMAHWSNK                                               
         . . .::                                                   
sp|P69 DALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHA
           70        80        90       100       110       120    

>>sp|P00517|KAPCA_BOVIN cAMP-dependent protein kinase ca  (351 aa)
 initn:  43 init1:  43 opt:  43  Z-score: 65.7  bits: 20.0 E(12): 0.86
Smith-Waterman score: 54; 23.6% identity (47.2% similar) in 72 aa overlap (137-206:229-300)

        110       120       130       140       150       160      
sp|P10 TRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQ
                                     .:    :  :.:: . .   . ..   .  
sp|P00 LCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRF
      200       210       220       230       240       250        

        170       180         190       200       210              
sp|P10 YRMFEPKCLDAFPNLR--DFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK      
          :     : . ::   :.  :: .::.    .:. ...::                  
sp|P00 PSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPK
      260       270       280       290       300       310        

sp|P00 FKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF
      320       330       340       350 

>>sp|P14960|RBS_GUITH Ribulose bisphosphate carboxylase   (139 aa)
 initn:  56 init1:  36 opt:  36  Z-score: 63.3  bits: 18.2 E(12):  1.2
Smith-Waterman score: 36; 57.1% identity (85.7% similar) in 7 aa overlap (7-13:47-53)

                                       10        20        30      
sp|P10                         MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMG
                                     ::.. ::                       
sp|P14 EQIVKQIQYAISKNWALNVEWTDDPHPRNAYWDLWGLPLFGIKDPAAVMFEINACRKAKP
         20        30        40        50        60        70      

         40        50        60        70        80        90      
sp|P10 DAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARKHHLDGETEEERIR
                                                                   
sp|P14 ACYVKVNAFDNSRGVESCCLSFIVQRPTSNEPGFQLIRSEVDSRNIRYTIQSYASTRPEG
         80        90       100       110       120       130      

>>sp|P01593|KV101_HUMAN Ig kappa chain V-I region AG OS=  (108 aa)
 initn:  31 init1:  31 opt:  31  Z-score: 58.3  bits: 16.9 E(12):  2.1
Smith-Waterman score: 33; 36.0% identity (54.0% similar) in 50 aa overlap (150-190:16-64)

     120       130       140       150       160          170      
sp|P10 FEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDI---LDQYRMFE---PK
                                     ::.:: .   . ::   :. :..     ::
sp|P01                DIQMTQSPSSLSASVGDRVTITCQASQDINHYLNWYQQGPKKAPK
                              10        20        30        40     

             180        190       200       210                    
sp|P10 CL--DAFPNLRDFL-ARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK            
        :  ::  ::.  . .:: :                                        
sp|P01 ILIYDA-SNLETGVPSRFSGSGFGTDFTFTISGLQPEDIATYYCQQYDTLPRTFGQGTKL
          50         60        70        80        90       100    

>>sp|P99998|CYC_PANTR Cytochrome c OS=Pan troglodytes GN  (105 aa)
 initn:  30 init1:  30 opt:  30  Z-score: 57.2  bits: 16.6 E(12):  2.4
Smith-Waterman score: 36; 26.5% identity (54.4% similar) in 68 aa overlap (129-193:28-88)

      100       110       120       130       140          150     
sp|P10 IVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWF---AGDKVTY
                                     :: :.   :...  :. : .   :..:   
sp|P99    MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKG
                  10        20        30        40        50       

         160       170       180       190       200       210     
sp|P10 VDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHW
       . . . : : .: . .::    .:. . .   : :.::                      
sp|P99 I-IWGEDTLMEY-LENPK--KYIPGTKMI---FVGIKKKEERADLIAYLKKATNE     
         60         70          80           90       100          

          
sp|P10 SNK

>>sp|P02585|TNNC2_HUMAN Troponin C, skeletal muscle OS=H  (160 aa)
 initn:  30 init1:  30 opt:  30  Z-score: 53.9  bits: 16.6 E(12):  3.5
Smith-Waterman score: 41; 25.9% identity (63.0% similar) in 54 aa overlap (44-92:13-62)

            20        30        40        50        60        70   
sp|P10 THPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKIT--
                                     :. .  .:: ..:.    .  ::.  :.  
sp|P02                   MTDQQAEARSYLSEEMIAEFKAAFDM----FDADGGGDISVK
                                 10        20            30        

              80           90       100       110       120        
sp|P10 QSNAILRYLAR---KHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFL
       . ....:.:..   :..::.  ::                                    
sp|P02 ELGTVMRMLGQTPTKEELDAIIEEVDEDGSGTIDFEEFLVMMVRQMKEDAKGKSEEELAE
       40        50        60        70        80        90        

>>sp|P03435|HEMA_I75A3 Hemagglutinin OS=Influenza A viru  (567 aa)
 initn:  37 init1:  37 opt:  37  Z-score: 53.7  bits: 18.4 E(12):  3.5
Smith-Waterman score: 50; 28.2% identity (64.1% similar) in 39 aa overlap (74-111:399-437)

            50        60        70        80         90       100  
sp|P10 SQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARK-HHLDGETEEERIRADIVEN
                                     : ...   .: :... :  : . : . .:.
sp|P03 GFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEK
      370       380       390       400       410       420        

            110       120       130       140       150       160  
sp|P10 QVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYD
        : ::...:                                                   
sp|P03 YVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEDMGNGCFKIYH
      430       440       450       460       470       480        

>>sp|P60615|NXL1A_BUNMU Alpha-bungarotoxin isoform A31 O  (95 aa)
 initn:  26 init1:  26 opt:  26  Z-score: 52.4  bits: 15.6 E(12):    4
Smith-Waterman score: 30; 54.5% identity (63.6% similar) in 11 aa overlap (115-125:86-94)

           90       100       110       120       130       140    
sp|P10 HLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGK
                                     : ::   ::.:                   
sp|P60 SRGKVVELGCAATCPSKKPYEEVTCCSTDKC-NPH-PKQRPG                  
          60        70        80          90                       

          150       160       170       180       190       200    
sp|P10 RPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYI

>>sp|P00193|FER_PEPAS Ferredoxin OS=Peptostreptococcus a  (54 aa)
 initn:  22 init1:  22 opt:  22  Z-score: 51.3  bits: 14.6 E(12):  4.5
Smith-Waterman score: 25; 50.0% identity (100.0% similar) in 4 aa overlap (171-174:15-18)

              150       160       170       180       190       200
sp|P10 FLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKS
                                     .:.:                          
sp|P00                 AYVINDSCIACGACKPECPVNIQQGSIYAIDADSCIDCGSCASV
                               10        20        30        40    

              210        
sp|P10 SRYIATPIFSKMAHWSNK
                         
sp|P00 CPVGAPNPED        
           50            

>>sp|P01834|IGKC_HUMAN Ig kappa chain C region OS=Homo s  (106 aa)
 initn:  23 init1:  23 opt:  23  Z-score: 47.4  bits: 14.9 E(12):  6.5
Smith-Waterman score: 35; 18.3% identity (53.5% similar) in 71 aa overlap (58-124:12-82)

        30        40        50        60        70          80     
sp|P10 YDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKIT--QSNAILRYLARKHH
                                     :.   : .:. ...   .:   :    . .
sp|P01                    TVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWK
                                  10        20        30        40 

          90       100       110         120       130       140   
sp|P10 LDGETEEERIRADIVENQVMDTRMQL--IMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLG
       .:.  .    . ...:..  :. ..:   .   . :.::.:                   
sp|P01 VDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSF
              50        60        70        80        90       100 

           150       160       170       180       190       200   
sp|P10 KRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRY
                                                                   
sp|P01 NRGEC                                                       
                                                                   



218 residues in 1 query   sequences
2267 residues in 12 library sequences
 Tcomplib [36.3.8h May, 2020] (8 proc in memory [0G])
 start: Tue Aug 31 10:28:51 2021 done: Tue Aug 31 10:28:51 2021
 Total Scan time:  0.010 Total Display time:  0.020

Function used was FASTA [36.3.8h May, 2020]
"""
