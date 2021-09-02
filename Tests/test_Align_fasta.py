# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.fasta module."""
import unittest

from Bio.Seq import Seq
from Bio.Align.fasta_m8 import AlignmentIterator

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.emboss."
    ) from None


class TestFasta(unittest.TestCase):

    query = "MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK"
    targets = {
        "sp|P09488|GSTM1_HUMAN": "MPMILGYWDIRGLAHAIRLLLEYTDSSYEEKKYTMGDAPDYDRSQWLNEKFKLGLDFPNLPYLIDGAHKITQSNAILCYIARKHNLCGETEEEKIRVDILENQTMDNHMQLGMICYNPEFEKLKPKYLEELPEKLKLYSEFLGKRPWFAGNKITFVDFLVYDVLDLHRIFEPKCLDAFPNLKDFISRFEGLEKISAYMKSSRFLPRPVFSKMAVWGNK",
        "sp|P00502|GSTA1_RAT": "MSGKPVLHYFNARGRMECIRWLLAAAGVEFDEKFIQSPEDLEKLKKDGNLMFDQVPMVEIDGMKLAQTRAILNYIATKYDLYGKDMKERALIDMYTEGILDLTEMIMQLVICPPDQKEAKTALAKDRTKNRYLPAFEKVLKSHGQDYLVGNRLTRVDIHLLELLLYVEEFDASLLTSFPLLKAFKSRISSLPNVKKFLQPGSQRKLPVDAKQIEEARKIFKF",
        "sp|P69905|HBA_HUMAN": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "sp|P00517|KAPCA_BOVIN": "MGNAAAAKKGSEQESVKEFLAKAKEDFLKKWENPAQNTAHLDQFERIKTLGTGSFGRVMLVKHMETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPKFKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF",
        "sp|P14960|RBS_GUITH": "MRLTQGAFSFLPDLTDEQIVKQIQYAISKNWALNVEWTDDPHPRNAYWDLWGLPLFGIKDPAAVMFEINACRKAKPACYVKVNAFDNSRGVESCCLSFIVQRPTSNEPGFQLIRSEVDSRNIRYTIQSYASTRPEGERY",
        "sp|P01593|KV101_HUMAN": "DIQMTQSPSSLSASVGDRVTITCQASQDINHYLNWYQQGPKKAPKILIYDASNLETGVPsrfsgsgfgtdftftisgLQPEDIATYYCQQYDTLPRTFGQGTKLEIKR",
        "sp|P99998|CYC_PANTR": "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE",
        "sp|P02585|TNNC2_HUMAN": "MTDQQAEARSYLSEEMIAEfkaafdmfdadgggdISVKELGTVMRMLGQTPTKEELDAIIEEVDEDGSGTIDFEEFLVMMVRQMKEDAKGKSEEELAECFRIFDRNADGYIDPEELAEIFRASGEHVTDEEIESLMKDGDKNNDGRIDFDEFLKMMEGVQ",
        "sp|P03435|HEMA_I75A3": "MKTIIALSYIFCLVFAQDLPGNDNNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRILDGINCTLIDALLGDPHCDGFQNEKWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFINEGFNWTGVTQNGGSSACKRGPDSGFFSRLNWLYKSGSTYPVQNVTMPNNDNSDKLYIWGVHHPSTDKEQTNLYVQASGKVTVSTKRSQQTIIPNVGSRPWVRGLSSRISIYWTIVKPGDILVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQKGNIRCNICI",
        "sp|P60615|NXL1A_BUNMU": "MKtllltlvvvtIVCLDLGYTIVCHTTATSPISAVTCPPGENLCYRKMWCDAFCSSRGKVVELGCAATCPSKKPYEEVTCCSTDKCNPHPKQRPG",
        "sp|P00193|FER_PEPAS": "AYVINDSCIACGACKPECPVNIQQGSIYAIDADSCIDCGSCASVCPVGAPNPED",
        "sp|P01834|IGKC_HUMAN": "TVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDstyslsstltlsKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
    }

    def test_m8CB(self):
        # Alignment file obtained by running
        # fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg
        # in the fasta36 source distribution
        path = "Fasta/output_m8CB.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(
                alignments.commandline,
                "fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg",
            )
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # sp|P10649|GSTM1_MOUSE   sp|P09488|GSTM1_HUMAN
        alignment = alignments[0]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 218)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 218)
        self.assertEqual(alignment.shape, (2, 218))
        self.assertEqual(alignment.sequences[0].id, "sp|P09488|GSTM1_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 48)
        self.assertAlmostEqual(alignment.annotations["evalue"], 6.1e-78)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 275.6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 218], [0, 218]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), len(target))
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = Seq(target)
        alignment.sequences[1].seq = Seq(query)
        self.assertEqual(
            alignment[0],
            "MPMILGYWDIRGLAHAIRLLLEYTDSSYEEKKYTMGDAPDYDRSQWLNEKFKLGLDFPNLPYLIDGAHKITQSNAILCYIARKHNLCGETEEEKIRVDILENQTMDNHMQLGMICYNPEFEKLKPKYLEELPEKLKLYSEFLGKRPWFAGNKITFVDFLVYDVLDLHRIFEPKCLDAFPNLKDFISRFEGLEKISAYMKSSRFLPRPVFSKMAVWGNK",
        )
        self.assertEqual(
            alignment[1],
            "MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P00502|GSTA1_RAT
        alignment = alignments[1]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 205)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 205)
        self.assertEqual(alignment.shape, (2, 223))
        self.assertEqual(alignment.sequences[0].id, "sp|P00502|GSTA1_RAT")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 144)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.1e-13)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 60.7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            5,
                            33,
                            33,
                            46,
                            48,
                            58,
                            59,
                            62,
                            62,
                            101,
                            102,
                            125,
                            127,
                            142,
                            144,
                            218,
                        ],
                        [
                            3,
                            31,
                            40,
                            53,
                            53,
                            63,
                            63,
                            66,
                            67,
                            106,
                            106,
                            129,
                            129,
                            144,
                            144,
                            218,
                        ],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 218)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = Seq(target)
        alignment.sequences[1].seq = Seq(query)
        self.assertEqual(
            alignment[0],
            "VLHYFNARGRMECIRWLLAAAGVEFDEK---------FIQSPEDLEKLKKDGNLMFDQVPMVEIDG-MKLAQTRAILNYIATKYDLYGKDMKERALIDMYTEGILDLTEMIMQLVICPPDQKEAKTALAKDRTKNRYLPAFEKVLKSHGQDYLVGNRLTRVDIHLLELLLYVEEFDASLLTSFPLLKAFKSRISSLPNVKKFLQPGSQRKLPVDAKQIEEARK",
        )
        self.assertEqual(
            alignment[1],
            "ILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKL--GLDFPNLPYL-IDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMD-TRMQLIMLCYNPDFEKQKPEFLK--TIPEKMKLYSEFLGK--RPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P69905|HBA_HUMAN

    def test_m8CC(self):
        # Alignment file obtained by running
        # fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg
        # in the fasta36 source distribution
        path = "Fasta/output_m8CC.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(
                alignments.commandline,
                "fasta36 -q -m 8CC seq/mgstm1.aa seq/prot_test.lseg",
            )
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
sp|P10 WFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKS-SRYIA
sp|P10 TPIFSKMAHWSNK    

sp|P69 ADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFD-LSHGSAQVKGHGKKVA
sp|P69 DALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHA


sp|P10 TRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQ
sp|P10 YRMFEPKCLDAFPNLR--DFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK      
                                     .:    :  :.:: . .   . ..   .  
sp|P00 LCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRF
sp|P00 PSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPK
sp|P00 FKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF

sp|P10                         MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMG
sp|P10 DAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARKHHLDGETEEERIR

sp|P14 EQIVKQIQYAISKNWALNVEWTDDPHPRNAYWDLWGLPLFGIKDPAAVMFEINACRKAKP
sp|P14 ACYVKVNAFDNSRGVESCCLSFIVQRPTSNEPGFQLIRSEVDSRNIRYTIQSYASTRPEG


sp|P10 FEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDI---LDQYRMFE---PK
sp|P10 CL--DAFPNLRDFL-ARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK            

sp|P01                DIQMTQSPSSLSASVGDRVTITCQASQDINHYLNWYQQGPKKAPK
sp|P01 ILIYDA-SNLETGVPSRFSGSGFGTDFTFTISGLQPEDIATYYCQQYDTLPRTFGQGTKL


sp|P10 IVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWF---AGDKVTY
sp|P10 VDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHW
sp|P10 SNK

sp|P99    MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKG
sp|P99 I-IWGEDTLMEY-LENPK--KYIPGTKMI---FVGIKKKEERADLIAYLKKATNE     


sp|P10 THPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKIT--
sp|P10 QSNAILRYLAR---KHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFL

sp|P02                   MTDQQAEARSYLSEEMIAEFKAAFDM----FDADGGGDISVK
sp|P02 ELGTVMRMLGQTPTKEELDAIIEEVDEDGSGTIDFEEFLVMMVRQMKEDAKGKSEEELAE


sp|P10 SQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARK-HHLDGETEEERIRADIVEN
sp|P10 QVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYD

sp|P03 GFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEK
sp|P03 YVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEDMGNGCFKIYH


sp|P10 HLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGK
sp|P10 RPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYI
sp|P60 SRGKVVELGCAATCPSKKPYEEVTCCSTDKC-NPH-PKQRPG                  


sp|P10 FLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKS
sp|P10 SRYIATPIFSKMAHWSNK

sp|P00                 AYVINDSCIACGACKPECPVNIQQGSIYAIDADSCIDCGSCASV
sp|P00 CPVGAPNPED        


sp|P10 YDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKIT--QSNAILRYLARKHH
sp|P10 LDGETEEERIRADIVENQVMDTRMQL--IMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLG
sp|P10 KRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRY

sp|P01                    TVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWK
sp|P01 VDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSF
sp|P01 NRGEC                                                       
"""
