# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unit tests for the Bio.AlignIO.FastaIO module."""
import unittest

from Bio.AlignIO import FastaIO


class FastaIOTests(unittest.TestCase):
    """Test FastaIO module."""

    def test_output001(self):
        """Check output001.m10 file."""
        fasta_file = "Fasta/output001.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 4)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 108)
            self.assertEqual(
                alignments[0][0].seq,
                "SGSNT-RRRAISRPVRLTAEED---QE"
                "IRKRAAECGKTVSGFLRAAALGKKVNS"
                "LTDDRVLKEVM-----RLGALQKKLFI"
                "DGKRVGDREYAEVLIAITEYHRALLSR",
            )
            self.assertEqual(alignments[0][0].id, "gi|10955263|ref|NP_052604.1|")
            self.assertEqual(alignments[0][0].annotations["original_length"], 107)
            self.assertEqual(
                alignments[0][1].seq,
                "AGSGAPRRRGSGLASRISEQSEALLQE"
                "AAKHAAEFGRS------EVDTEHLLLA"
                "LADSDVVKTILGQFKIKVDDLKRQIES"
                "EAKR-GDKPF-EGEIGVSPRVKDALSR",
            )
            self.assertEqual(alignments[0][1].id, "gi|152973457|ref|YP_001338508.1|")
            self.assertEqual(alignments[0][1].annotations["original_length"], 931)
            self.assertEqual(len(alignments[1]), 2)
            self.assertEqual(alignments[1].get_alignment_length(), 64)
            self.assertEqual(
                alignments[1][0].seq,
                "AAECGKTVSGFLRAAALGKKVNSLTDD"
                "RVLKEV-MRLGALQKKLFIDGKRVGDR"
                "EYAEVLIAIT",
            )
            self.assertEqual(alignments[1][0].id, "gi|10955263|ref|NP_052604.1|")
            self.assertEqual(alignments[1][0].annotations["original_length"], 107)
            self.assertEqual(
                alignments[1][1].seq,
                "ASRQGCTVGG--KMDSVQDKASDKDKE"
                "RVMKNINIMWNALSKNRLFDG----NK"
                "ELKEFIMTLT",
            )
            self.assertEqual(alignments[1][1].id, "gi|152973588|ref|YP_001338639.1|")
            self.assertEqual(alignments[1][1].annotations["original_length"], 459)
            self.assertEqual(len(alignments[2]), 2)
            self.assertEqual(alignments[2].get_alignment_length(), 38)
            self.assertEqual(
                alignments[2][0].seq, "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE"
            )
            self.assertEqual(alignments[2][0].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[2][0].annotations["original_length"], 126)
            self.assertEqual(
                alignments[2][1].seq, "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ"
            )
            self.assertEqual(alignments[2][1].id, "gi|152973462|ref|YP_001338513.1|")
            self.assertEqual(alignments[2][1].annotations["original_length"], 101)
            self.assertEqual(len(alignments[3]), 2)
            self.assertEqual(alignments[3].get_alignment_length(), 43)
            self.assertEqual(
                alignments[3][0].seq, "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC"
            )
            self.assertEqual(alignments[3][0].id, "gi|10955265|ref|NP_052606.1|")
            self.assertEqual(alignments[3][0].annotations["original_length"], 346)
            self.assertEqual(
                alignments[3][1].seq, "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC"
            )
            self.assertEqual(alignments[3][1].id, "gi|152973545|ref|YP_001338596.1|")
            self.assertEqual(alignments[3][1].annotations["original_length"], 242)

    def test_output002(self):
        """Check output002.m10 file."""
        fasta_file = "Fasta/output002.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 6)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 88)
            self.assertEqual(
                alignments[0][0].seq,
                "SGSNTRRRAISRPVR--LTAEED"
                "QEIRKRAAECG-KTVSGFLRAAA"
                "LGKKVNSLTDDRVLKEVMRLGAL"
                "QKKLFIDGKRVGDREYAEV",
            )
            self.assertEqual(alignments[0][0].id, "gi|10955263|ref|NP_052604.1|")
            self.assertEqual(alignments[0][0].annotations["original_length"], 107)
            self.assertEqual(
                alignments[0][1].seq,
                "SQRSTRRKPENQPTRVILFNKPY"
                "DVLPQFTDEAGRKTLKEFIPVQG"
                "VYAAGRLDRDSEGLLVLTNNGAL"
                "QARLTQPGKRTGKIYYVQV",
            )
            self.assertEqual(alignments[0][1].id, "gi|162139799|ref|NP_309634.2|")
            self.assertEqual(alignments[0][1].annotations["original_length"], 207)
            self.assertEqual(len(alignments[1]), 2)
            self.assertEqual(alignments[1].get_alignment_length(), 53)
            self.assertEqual(
                alignments[1][0].seq,
                "EIRKRAAECGKTVSGFLRAAA-LGKKV----NSLTDDRVLKEVMRLGALQKKL",
            )
            self.assertEqual(alignments[1][0].id, "gi|10955263|ref|NP_052604.1|")
            self.assertEqual(alignments[1][0].annotations["original_length"], 107)
            self.assertEqual(
                alignments[1][1].seq,
                "EIKPRGTSKGEAIAAFMQEAPFIGRTPVFLGDDLTDESGFAVVNRLGGMSVKI",
            )
            self.assertEqual(alignments[1][1].id, "gi|15831859|ref|NP_310632.1|")
            self.assertEqual(alignments[1][1].annotations["original_length"], 266)
            self.assertEqual(len(alignments[2]), 2)
            self.assertEqual(alignments[2].get_alignment_length(), 92)
            self.assertEqual(
                alignments[2][0].seq,
                "SEFFSKIESDLKKKKSKGDVFFD"
                "LIIPNG-----GKKDRYVYTSFN"
                "GEKFSSYTLNKVTKTDEYNDLSE"
                "LSASFFKKNFDKINVNLLSKATS",
            )
            self.assertEqual(alignments[2][0].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[2][0].annotations["original_length"], 126)
            self.assertEqual(
                alignments[2][1].seq,
                "TELNSELAKAMKVDAQRG-AFVS"
                "QVLPNSSAAKAGIKAGDVITSL"
                "NGKPISSFAALRA-QVGTMPVG"
                "SKLTLGLLRDG-KQVNVNLELQ"
                "QSS",
            )
            self.assertEqual(alignments[2][1].id, "gi|15829419|ref|NP_308192.1|")
            self.assertEqual(alignments[2][1].annotations["original_length"], 474)
            self.assertEqual(len(alignments[3]), 2)
            self.assertEqual(alignments[3].get_alignment_length(), 73)
            self.assertEqual(
                alignments[3][0].seq,
                "FFDLIIPNGGKKDRYVYTSFNGE"
                "KFSSYTLNKVTKTDEYNDLSELS"
                "ASFFKKNFDKINVNLLSKATSFA"
                "LKKG",
            )
            self.assertEqual(alignments[3][0].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[3][0].annotations["original_length"], 126)
            self.assertEqual(
                alignments[3][1].seq,
                "LFDLFLKNDAMHDPMVNESYC-E"
                "TFGWVSKENLARMKE---LTYKA"
                "NDVLKKLFDDAGLILVDFKLEFG"
                "LYKG",
            )
            self.assertEqual(alignments[3][1].id, "gi|15832592|ref|NP_311365.1|")
            self.assertEqual(alignments[3][1].annotations["original_length"], 237)
            self.assertEqual(len(alignments[4]), 2)
            self.assertEqual(alignments[4].get_alignment_length(), 63)
            self.assertEqual(
                alignments[4][0].seq,
                "VDIKK-ETIESELHSKLPKSIDKIHEDIKKQLSCSLI--MKKID-VEMEDYSTYCFSALRAIE",
            )
            self.assertEqual(alignments[4][0].id, "gi|10955265|ref|NP_052606.1|")
            self.assertEqual(alignments[4][0].annotations["original_length"], 346)
            self.assertEqual(
                alignments[4][1].seq,
                "IDPKKIEQIARQVHESMPKGIREFGEDVEKKIRQTLQAQLTRLDLVSREEFDVQTQVLLRTRE",
            )
            self.assertEqual(alignments[4][1].id, "gi|38704138|ref|NP_311957.2|")
            self.assertEqual(alignments[4][1].annotations["original_length"], 111)
            self.assertEqual(len(alignments[5]), 2)
            self.assertEqual(alignments[5].get_alignment_length(), 157)
            self.assertEqual(
                alignments[5][0].seq,
                "QYIMTTSNGDRVRAKIYKRGSIQ"
                "FQGKYLQIASLINDFMCSILNMK"
                "EIVEQKNKEFNVDI---KKETI-"
                "ESELHSKLPKSIDKIHEDIKKQL"
                "SCSLIMKKIDV-EMEDYSTYCFS"
                "ALRA-IEGFIYQILNDVCNPSSS"
                "KNLGEYFTENKPKYIIREI",
            )
            self.assertEqual(alignments[5][0].id, "gi|10955265|ref|NP_052606.1|")
            self.assertEqual(alignments[5][0].annotations["original_length"], 346)
            self.assertEqual(
                alignments[5][1].seq,
                "EFIRLLSDHDQFEKDQISELTVA"
                "ANALKLEVAK--NNY-----NMK"
                "YSFDTQTERRMIELIREQKDLIP"
                "EKYLHQSGIKKL-KLHED---EF"
                "SSLLVDAERQVLEGSSFVLCCGE"
                "KINSTISELLSKKITDLTHPTES"
                "FTLSEYFSYDVYEEIFKKV",
            )
            self.assertEqual(alignments[5][1].id, "gi|15833861|ref|NP_312634.1|")
            self.assertEqual(alignments[5][1].annotations["original_length"], 330)

    def test_output003(self):
        """Check output003.m10 file."""
        fasta_file = "Fasta/output003.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 3)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 55)
            self.assertEqual(
                alignments[0][0].seq,
                "ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGGGLRLSASTKTVDQLVRIAA",
            )
            self.assertEqual(alignments[0][0].id, "gi|152973837|ref|YP_001338874.1|")
            self.assertEqual(alignments[0][0].annotations["original_length"], 183)
            self.assertEqual(
                alignments[0][1].seq,
                "VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGKKVNSLTDDRVLKEVMRLGA",
            )
            self.assertEqual(alignments[0][1].id, "gi|10955263|ref|NP_052604.1|")
            self.assertEqual(alignments[0][1].annotations["original_length"], 107)
            self.assertEqual(len(alignments[1]), 2)
            self.assertEqual(alignments[1].get_alignment_length(), 22)
            self.assertEqual(alignments[1][0].seq, "DDAEHLFRTLSSR-LDALQDGN")
            self.assertEqual(alignments[1][0].id, "gi|152973840|ref|YP_001338877.1|")
            self.assertEqual(alignments[1][0].annotations["original_length"], 63)
            self.assertEqual(alignments[1][1].seq, "DDRANLFEFLSEEGITITEDNN")
            self.assertEqual(alignments[1][1].id, "gi|10955265|ref|NP_052606.1|")
            self.assertEqual(alignments[1][1].annotations["original_length"], 346)
            self.assertEqual(len(alignments[2]), 2)
            self.assertEqual(alignments[2].get_alignment_length(), 63)
            self.assertEqual(
                alignments[2][0].seq,
                "VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIRHLKQEIEHLNTLLLSKDSHIDSLKQAM",
            )
            self.assertEqual(alignments[2][0].id, "gi|152973841|ref|YP_001338878.1|")
            self.assertEqual(alignments[2][0].annotations["original_length"], 133)
            self.assertEqual(
                alignments[2][1].seq,
                "VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATSF-ALKKGI",
            )
            self.assertEqual(alignments[2][1].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[2][1].annotations["original_length"], 126)

    def test_output004(self):
        """Check output004.m10 file."""
        fasta_file = "Fasta/output004.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 1)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 102)
            self.assertEqual(
                alignments[0][0].seq,
                "AAAAAAGATAAAAAATATCAAAT"
                "AGAAGCAATAAAAAATAAAGATA"
                "AAACTTTATTTATTGTCTATGCT"
                "ACTGATATTTATAGCCCGAGCGA"
                "ATTTTTCTCA",
            )
            self.assertEqual(alignments[0][0].id, "ref|NC_002127.1|:c1351-971")
            self.assertEqual(alignments[0][0].annotations["original_length"], 381)
            self.assertEqual(
                alignments[0][1].seq,
                "AGAGAAAATAAAACAAGTAATAA"
                "AATATTAATGGAAAAAATAAATT"
                "CTTGTTTATTTAGACCTGATTCT"
                "AATCACTTTTCTTGCCCGGAGTC"
                "ATTTTTGACA",
            )
            self.assertEqual(alignments[0][1].id, "ref|NC_002695.1|:1970775-1971404")
            self.assertEqual(alignments[0][1].annotations["original_length"], 630)

    def test_output005(self):
        """Check output005.m10 file."""
        fasta_file = "Fasta/output005.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 1)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 110)
            self.assertEqual(
                alignments[0][0].seq,
                "IKNKDKTLFIVYAT-DIYSPSEF"
                "FSKIESDLKKKKSKGDV--FFDL"
                "IIPNGGKKD--RYVYTSFNGEKF"
                "SSYTLNKVTKTDEYNDL--SELS"
                "ASFFKKNFDKINVNLLSK",
            )
            self.assertEqual(alignments[0][0].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[0][0].annotations["original_length"], 126)
            self.assertEqual(
                alignments[0][1].seq,
                "IKDELPVAFCSWASLDLECEVKY"
                "INDVTSLYAKDWMSGERKWFIDW"
                "IAPFGHNMELYKYMRKKYPYELF"
                "RAIRLDESSKTGKIAEFHGGGID"
                "KKLASKIFRQYHHELMSE",
            )
            self.assertEqual(alignments[0][1].id, "gi|10955282|ref|NP_052623.1|")
            self.assertEqual(alignments[0][1].annotations["original_length"], 163)

    def test_output006(self):
        """Check output006.m10 file."""
        fasta_file = "Fasta/output006.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 1)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 131)
            self.assertEqual(
                alignments[0][0].seq,
                "GCAACGCTTCAAGAACTGGAATT"
                "AGGAACCGTGACAACGATTAATG"
                "AGGAGATTTATGAAGAGGGTTCT"
                "TCGATTTTAGGCCAATCGGAAGG"
                "AATTATGTAGCAAGTCCATCAGA"
                "AAATGGAAGAAGTCAT",
            )
            self.assertEqual(alignments[0][0].id, "query")
            self.assertEqual(alignments[0][0].annotations["original_length"], 131)
            self.assertEqual(
                alignments[0][1].seq,
                "GCAACGCTTCAAGAACTGGAATT"
                "AGGAACCGTGACAACGATTAATG"
                "AGGAGATTTATGAAGAGGGTTCT"
                "TCGATTTTAGGCCAATCGGAAGG"
                "AATTATGTAGCAAGTCCATCAGA"
                "AAATGGAAGTAGTCAT",
            )
            self.assertEqual(alignments[0][1].id, "gi|116660610|gb|EG558221.1|EG558221")
            self.assertEqual(alignments[0][1].annotations["original_length"], 573)

    def test_output007(self):
        """Check output007.m10 file."""
        fasta_file = "Fasta/output007.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 9)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 108)
            self.assertEqual(
                alignments[0][0].seq,
                "SGSNT-RRRAISRPVRLTAEED-"
                "--QEIRKRAAECGKTVSGFLRAA"
                "ALGKKVNSLTDDRVLKEVM----"
                "-RLGALQKKLFIDGKRVGDREYA"
                "EVLIAITEYHRALLSR",
            )
            self.assertEqual(alignments[0][0].id, "gi|10955263|ref|NP_052604.1|")
            self.assertEqual(alignments[0][0].annotations["original_length"], 107)
            self.assertEqual(
                alignments[0][1].seq,
                "AGSGAPRRRGSGLASRISEQSEA"
                "LLQEAAKHAAEFGRS------EV"
                "DTEHLLLALADSDVVKTILGQFK"
                "IKVDDLKRQIESEAKR-GDKPF-"
                "EGEIGVSPRVKDALSR",
            )
            self.assertEqual(alignments[0][1].id, "gi|152973457|ref|YP_001338508.1|")
            self.assertEqual(alignments[0][1].annotations["original_length"], 931)
            self.assertEqual(len(alignments[1]), 2)
            self.assertEqual(alignments[1].get_alignment_length(), 64)
            self.assertEqual(
                alignments[1][0].seq,
                "AAECGKTVSGFLRAAALGKKVNS"
                "LTDDRVLKEV-MRLGALQKKLFI"
                "DGKRVGDREYAEVLIAIT",
            )
            self.assertEqual(alignments[1][0].id, "gi|10955263|ref|NP_052604.1|")
            self.assertEqual(alignments[1][0].annotations["original_length"], 107)
            self.assertEqual(
                alignments[1][1].seq,
                "ASRQGCTVGG--KMDSVQDKASD"
                "KDKERVMKNINIMWNALSKNRLF"
                "DG----NKELKEFIMTLT",
            )
            self.assertEqual(alignments[1][1].id, "gi|152973588|ref|YP_001338639.1|")
            self.assertEqual(alignments[1][1].annotations["original_length"], 459)
            self.assertEqual(len(alignments[2]), 2)
            self.assertEqual(alignments[2].get_alignment_length(), 45)
            self.assertEqual(
                alignments[2][0].seq, "EIRKRAAECGKTVSGFLRAAA-----LGKKVNSLTDDRVLKEVMR"
            )
            self.assertEqual(alignments[2][0].id, "gi|10955263|ref|NP_052604.1|")
            self.assertEqual(alignments[2][0].annotations["original_length"], 107)
            self.assertEqual(
                alignments[2][1].seq, "ELVKLIADMGISVRALLRKNVEPYEELGLEEDKFTDDQLIDFMLQ"
            )
            self.assertEqual(alignments[2][1].id, "gi|152973480|ref|YP_001338531.1|")
            self.assertEqual(alignments[2][1].annotations["original_length"], 141)
            self.assertEqual(len(alignments[3]), 2)
            self.assertEqual(alignments[3].get_alignment_length(), 38)
            self.assertEqual(
                alignments[3][0].seq, "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE"
            )
            self.assertEqual(alignments[3][0].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[3][0].annotations["original_length"], 126)
            self.assertEqual(
                alignments[3][1].seq, "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ"
            )
            self.assertEqual(alignments[3][1].id, "gi|152973462|ref|YP_001338513.1|")
            self.assertEqual(alignments[3][1].annotations["original_length"], 101)
            self.assertEqual(len(alignments[4]), 2)
            self.assertEqual(alignments[4].get_alignment_length(), 11)
            self.assertEqual(alignments[4][0].seq, "FFDLIIPNGGK")
            self.assertEqual(alignments[4][0].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[4][0].annotations["original_length"], 126)
            self.assertEqual(alignments[4][1].seq, "FFDLVIENPGK")
            self.assertEqual(alignments[4][1].id, "gi|152973509|ref|YP_001338560.1|")
            self.assertEqual(alignments[4][1].annotations["original_length"], 448)
            self.assertEqual(len(alignments[5]), 2)
            self.assertEqual(alignments[5].get_alignment_length(), 40)
            self.assertEqual(
                alignments[5][0].seq, "DKTLFIVYATDIYSPSE-FFSKIESDLKKKKSKGD-VFFD"
            )
            self.assertEqual(alignments[5][0].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[5][0].annotations["original_length"], 126)
            self.assertEqual(
                alignments[5][1].seq, "ESVVFILMAGFAMSVCYLFFSVLEKVINARKSKDESIYHD"
            )
            self.assertEqual(alignments[5][1].id, "gi|152973581|ref|YP_001338632.1|")
            self.assertEqual(alignments[5][1].annotations["original_length"], 84)
            self.assertEqual(len(alignments[6]), 2)
            self.assertEqual(alignments[6].get_alignment_length(), 30)
            self.assertEqual(alignments[6][0].seq, "ASFFKKNFDKINVNLLSKATSFALKKGIPI")
            self.assertEqual(alignments[6][0].id, "gi|10955264|ref|NP_052605.1|")
            self.assertEqual(alignments[6][0].annotations["original_length"], 126)
            self.assertEqual(alignments[6][1].seq, "ASFSKEEQDKVAVDKVAADVAWQERMNKPV")
            self.assertEqual(alignments[6][1].id, "gi|152973536|ref|YP_001338587.1|")
            self.assertEqual(alignments[6][1].annotations["original_length"], 84)
            self.assertEqual(len(alignments[7]), 2)
            self.assertEqual(alignments[7].get_alignment_length(), 43)
            self.assertEqual(
                alignments[7][0].seq, "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC"
            )
            self.assertEqual(alignments[7][0].id, "gi|10955265|ref|NP_052606.1|")
            self.assertEqual(alignments[7][0].annotations["original_length"], 346)
            self.assertEqual(
                alignments[7][1].seq, "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC"
            )
            self.assertEqual(alignments[7][1].id, "gi|152973545|ref|YP_001338596.1|")
            self.assertEqual(alignments[7][1].annotations["original_length"], 242)
            self.assertEqual(len(alignments[8]), 2)
            self.assertEqual(alignments[8].get_alignment_length(), 64)
            self.assertEqual(
                alignments[8][0].seq,
                "ISGTYKGIDFLIKLMPSGGNTTI"
                "GRASGQNNTYFDEIALIIKENCL"
                "Y--SDTKNFEYTIPKFSD",
            )
            self.assertEqual(alignments[8][0].id, "gi|10955265|ref|NP_052606.1|")
            self.assertEqual(alignments[8][0].annotations["original_length"], 346)
            self.assertEqual(
                alignments[8][1].seq,
                "IDGVITAFD-LRTGMNISKDKVV"
                "AQIQGMDPVW---ISAAVPESIA"
                "YLLKDTSQFEISVPAYPD",
            )
            self.assertEqual(alignments[8][1].id, "gi|152973505|ref|YP_001338556.1|")
            self.assertEqual(alignments[8][1].annotations["original_length"], 430)

    def test_output008(self):
        """Check output008.m10 file."""
        fasta_file = "Fasta/output008.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 12)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 65)
            self.assertEqual(
                alignments[0][0].seq,
                "LQHRHPHQQQQQQQQQQQQQQQQ"
                "QQQQQQQQQQQH---HHHHHHHL"
                "LQDAYMQQYQHATQQQQML",
            )
            self.assertEqual(alignments[0][0].id, "sp|Q9NSY1|BMP2K_HUMAN")
            self.assertEqual(alignments[0][0].annotations["original_length"], 1161)
            self.assertEqual(
                alignments[0][1].seq,
                "IPHQLPHALRHRPAQEAAHASQL"
                "HPAQPGCGQPLHGLWRLHHHPVY"
                "LYAWILRLRGHGMQSGGLL",
            )
            self.assertEqual(alignments[0][1].id, "gi|283855822|gb|GQ290312.1|")
            self.assertEqual(alignments[0][1].annotations["original_length"], 983)
            self.assertEqual(len(alignments[1]), 2)
            self.assertEqual(alignments[1].get_alignment_length(), 201)
            self.assertEqual(
                alignments[1][0].seq,
                "GPEIL---LGQ-GPPQQPPQQHR"
                "VLQQLQQGDWRLQQLH-------"
                "LQHRHPHQQQQQQQQQQQQQQQQ"
                "QQQQQQQQQQQH-----HHHHHH"
                "-HLLQDAYMQQYQHATQQQQMLQ"
                "QQF-LMHSVYQPQPSASQYPTMM"
                "PQYQQAFFQQQMLAQHQPSQQQA"
                "SPEYLTSPQEFSPALVSYTSSLP"
                "A-QVGTIMDSSYSANRS",
            )
            self.assertEqual(alignments[1][0].id, "sp|Q9NSY1|BMP2K_HUMAN")
            self.assertEqual(alignments[1][0].annotations["original_length"], 1161)
            self.assertEqual(
                alignments[1][1].seq,
                "GPELLRALLQQNGCGTQPLRVPT"
                "VLPG*AMAVLHAGRLHVPAHRAW"
                "LPHQLPHALRHGPAQEAAHASQL"
                "HPAQPGRG*PLHGLRWLHHHPLH"
                "/PLCMDTLSLGPQDAIWRASLPH"
                "WAVKLPCGLWWSWPLSGTWWCVS"
                "P*ATSA------LGRTMP*WASL"
                "SPGSWHWPALHPPSLVGPGTSLK"
                "ACSVHAGSTTTHSSQKS",
            )
            self.assertEqual(alignments[1][1].id, "gi|57163782|ref|NM_001009242.1|")
            self.assertEqual(alignments[1][1].annotations["original_length"], 1047)
            self.assertEqual(len(alignments[2]), 2)
            self.assertEqual(alignments[2].get_alignment_length(), 348)
            self.assertEqual(
                alignments[2][0].seq,
                "MNGTEGPNFYVPFSNATGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVLGGFT"
                "STLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLAGWSRYIPEGLQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIIIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTIPAFFAKSAA"
                "IYNPVIYIMMNKQFRNCMLTTIC"
                "CGKNPLGDDEASATVSKTETSQV"
                "APA",
            )
            self.assertEqual(alignments[2][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[2][0].annotations["original_length"], 348)
            self.assertEqual(
                alignments[2][1].seq,
                "MNGTEGPNFYVPFSNKTGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVFGGFT"
                "TTLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLVGWSRYIPEGMQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIVIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTLPAFFAKSSS"
                "IYNPVIYIMMNKQFRNCMLTTLC"
                "CGKNPLGDDEASTTGSKTETSQV"
                "APA",
            )
            self.assertEqual(alignments[2][1].id, "gi|57163782|ref|NM_001009242.1|")
            self.assertEqual(alignments[2][1].annotations["original_length"], 1047)
            self.assertEqual(len(alignments[3]), 2)
            self.assertEqual(alignments[3].get_alignment_length(), 348)
            self.assertEqual(
                alignments[3][0].seq,
                "MNGTEGPNFYVPFSNATGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVLGGFT"
                "STLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLAGWSRYIPEGLQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIIIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTIPAFFAKSAA"
                "IYNPVIYIMMNKQFRNCMLTTIC"
                "CGKNPLGDDEASATVSKTETSQV"
                "APA",
            )
            self.assertEqual(alignments[3][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[3][0].annotations["original_length"], 348)
            self.assertEqual(
                alignments[3][1].seq,
                "MNGTEGPNFYVPFSNKTGVVRSP"
                "FEAPQYYLAEPWQFSMLAAYMFL"
                "LIMLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVFGGFT"
                "TTLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLVGWSRYIPEGMQ"
                "CSCGIDYYTPHEETNNESFVIYM"
                "FVVHFIIPLIVIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWLPYAGVAFYIF"
                "THQGSDFGPIFMTIPAFFAKTSA"
                "VYNPVIYIMMNKQFRNCMVTTLC"
                "CGKNPLGDDEASTTVSKTETSQV"
                "APA",
            )
            self.assertEqual(alignments[3][1].id, "gi|18148870|dbj|AB062417.1|")
            self.assertEqual(alignments[3][1].annotations["original_length"], 1047)
            self.assertEqual(len(alignments[4]), 2)
            self.assertEqual(alignments[4].get_alignment_length(), 326)
            self.assertEqual(
                alignments[4][0].seq,
                "VPFSNATGVVRSPFEYPQYYLAE"
                "PWQFSMLAAYMFLLIVLGFPINF"
                "LTLYVTVQHKKLRTPLNYILLNL"
                "AVADLFMVLGGFTSTLYTSLHGY"
                "FVFGPTGCNLEGFFATLGGEIAL"
                "WSLVVLAIERYVVVCKPMSNFRF"
                "GENHAIMGVAFTWVMALACAAPP"
                "LAGWSRYIPEGLQCSCGIDYYTL"
                "KPEVNNESFVIYMFVVHFTIPMI"
                "IIFFCYGQLVFTVKEAAAQQQES"
                "ATTQKAEKEVTRMVIIMVIAFLI"
                "CWVPYASVAFYIFTHQGSNFGPI"
                "FMTIPAFFAKSAAIYNPVIYIMM"
                "NKQFRNCMLTTICCGKNPLGDDE"
                "ASAT",
            )
            self.assertEqual(alignments[4][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[4][0].annotations["original_length"], 348)
            self.assertEqual(
                alignments[4][1].seq,
                "VPFSNKTGVVRSPFEYPQYYLAE"
                "PWQFSMLAAYMFLLIVLGFPINF"
                "LTLYVTVQHKKLRTPLNYILLNL"
                "AVANLFMVFGGFTTTLYTSMHGY"
                "FVFGATGCNLEGFFATLGGEIAL"
                "WSLVVLAIERYVVVCKPMSNFRF"
                "GENHAIMGLAFTWVMALACAAPP"
                "LAGWSRYIPEGMQCSCGIDYYTL"
                "KPEVNNESFVIYMFVVHFTIPMI"
                "VIFFCYGQLVFTVKEAAAQQQES"
                "ATTQKAEKEVTRMVIIMVVAFLI"
                "CWLPYASVAFYIFTHQGSNFGPV"
                "FMTIPAFFAKSSSIYNPVIYIMM"
                "NKQFRNCMLTTLCCGKNPLGDDE"
                "ASTT",
            )
            self.assertEqual(alignments[4][1].id, "gi|283855822|gb|GQ290312.1|")
            self.assertEqual(alignments[4][1].annotations["original_length"], 983)
            self.assertEqual(len(alignments[5]), 2)
            self.assertEqual(alignments[5].get_alignment_length(), 354)
            self.assertEqual(
                alignments[5][0].seq,
                "MNGTEGPNFYVPFSNATGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVLGGFT"
                "STLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLAGWSRYIPEGLQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIIIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTIPAFFAKSAA"
                "IYNPVIYIMMNKQFRNCMLTTIC"
                "CGKNPLGDDEAS-ATVSKTE---"
                "--TSQVAPA",
            )
            self.assertEqual(alignments[5][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[5][0].annotations["original_length"], 348)
            self.assertEqual(
                alignments[5][1].seq,
                "MNGTEGPNFYIPMSNKTGVVRSP"
                "FEYPQYYLAEPWQYSILCAYMFL"
                "LILLGFPINFMTLYVTIQHKKLR"
                "TPLNYILLNLAFANHFMVLCGFT"
                "VTMYSSMNGYFILGATGCYVEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFSENHAVMGVAFTW"
                "IMALSCAVPPLLGWSRYIPEGMQ"
                "CSCGVDYYTLKPEVNNESFVIYM"
                "FVVHFTIPLIIIFFCYGRLVCTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVVFFLICWVPYASVAFFIF"
                "SNQGSEFGPIFMTVPAFFAKSSS"
                "IYNPVIYIMLNKQFRNCMITTLC"
                "CGKNPFGEDDASSAATSKTEASS"
                "VSSSQVSPA",
            )
            self.assertEqual(alignments[5][1].id, "gi|2734705|gb|U59921.1|BBU59921")
            self.assertEqual(alignments[5][1].annotations["original_length"], 1574)
            self.assertEqual(len(alignments[6]), 2)
            self.assertEqual(alignments[6].get_alignment_length(), 347)
            self.assertEqual(
                alignments[6][0].seq,
                "MNGTEGPNFYVPFSNATGVVRSP"
                "FEYPQYYLAEPWQFSMLAAYMFL"
                "LIVLGFPINFLTLYVTVQHKKLR"
                "TPLNYILLNLAVADLFMVLGGFT"
                "STLYTSLHGYFVFGPTGCNLEGF"
                "FATLGGEIALWSLVVLAIERYVV"
                "VCKPMSNFRFGENHAIMGVAFTW"
                "VMALACAAPPLAGWSRYIPEGLQ"
                "CSCGIDYYTLKPEVNNESFVIYM"
                "FVVHFTIPMIIIFFCYGQLVFTV"
                "KEAAAQQQESATTQKAEKEVTRM"
                "VIIMVIAFLICWVPYASVAFYIF"
                "THQGSNFGPIFMTIPAFFAKSAA"
                "IYNPVIYIMMNKQFRNCMLTTIC"
                "CGKNPLGD-DEASATVSKTETSQ"
                "VA",
            )
            self.assertEqual(alignments[6][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[6][0].annotations["original_length"], 348)
            self.assertEqual(
                alignments[6][1].seq,
                "MNGTEGPNFYIPMSNATGVVRSP"
                "FEYPQYYLAEPWAFSALSAYMFF"
                "LIIAGFPINFLTLYVTIEHKKLR"
                "TPLNYILLNLAVADLFMVFGGFT"
                "TTMYTSMHGYFVFGPTGCNIEGF"
                "FATLGGEIALWCLVVLAIERWMV"
                "VCKPVTNFRFGESHAIMGVMVTW"
                "TMALACALPPLFGWSRYIPEGLQ"
                "CSCGIDYYTRAPGINNESFVIYM"
                "FTCHFSIPLAVISFCYGRLVCTV"
                "KEAAAQQQESETTQRAEREVTRM"
                "VVIMVISFLVCWVPYASVAWYIF"
                "THQGSTFGPIFMTIPSFFAKSSA"
                "LYNPMIYICMNKQFRHCMITTLC"
                "CGKNPFEEEDGASATSSKTEASS"
                "VS",
            )
            self.assertEqual(alignments[6][1].id, "gi|12583664|dbj|AB043817.1|")
            self.assertEqual(alignments[6][1].annotations["original_length"], 1344)
            self.assertEqual(len(alignments[7]), 2)
            self.assertEqual(alignments[7].get_alignment_length(), 111)
            self.assertEqual(
                alignments[7][0].seq,
                "VPFSNATGVVRSPFEYPQYYLAE"
                "PWQFSMLAAYMFLLIVLGFPINF"
                "LTLYVTVQHKKLRTPLNYILLNL"
                "AVADLFMVLGGFTSTLYTSLHGY"
                "FVFGPTGCNLEGFFATLGG",
            )
            self.assertEqual(alignments[7][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[7][0].annotations["original_length"], 348)
            self.assertEqual(
                alignments[7][1].seq,
                "VPFSNKTGVVRSPFEHPQYYLAE"
                "PWQFSMLAAYMFLLIVLGFPINF"
                "LTLYVTVQHKKLRTPLNYILLNL"
                "AVADLFMVFGGFTTTLYTSLHGY"
                "FVFGPTGCNLEGFFATLGG",
            )
            self.assertEqual(alignments[7][1].id, "gi|283855845|gb|GQ290303.1|")
            self.assertEqual(alignments[7][1].annotations["original_length"], 4301)
            self.assertEqual(len(alignments[8]), 2)
            self.assertEqual(alignments[8].get_alignment_length(), 172)
            self.assertEqual(
                alignments[8][0].seq,
                "RYIPEGLQCSCGIDYYTLKPEVN"
                "NESFVIYMFVVHFTIPMIIIFFC"
                "YGQLVFTVKE-------------"
                "-----------------------"
                "AAAQQQESATTQKAEKEVTRMVI"
                "IMVIAFLICWVPYASVAFYIFTH"
                "QGSNFGPIFMTIPAFFAKSAAIY"
                "NPVIYIMMNKQ",
            )
            self.assertEqual(alignments[8][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[8][0].annotations["original_length"], 348)
            self.assertEqual(
                alignments[8][1].seq,
                "RYIPEGMQCSCGIDYYTLKPEVN"
                "NESFVIYMFVVHFTIPMIVIFFC"
                "YGQLVFTVKEVRSCVGHWGHAH*"
                "VNGAQLHSQSCHSLDT*PCVPA\\"
                "AAAQQQESATTQKAEKEVTRMVI"
                "IMVIAFLICWLPYAGVAFYIFTH"
                "QGSNFGPIFMTLPAFFAKSSSIY"
                "NPVIYIMMNKQ",
            )
            self.assertEqual(alignments[8][1].id, "gi|283855845|gb|GQ290303.1|")
            self.assertEqual(alignments[8][1].annotations["original_length"], 4301)
            self.assertEqual(len(alignments[9]), 2)
            self.assertEqual(alignments[9].get_alignment_length(), 73)
            self.assertEqual(
                alignments[9][0].seq,
                "LGGEIALWSLVVLAIERYVVVCK"
                "PMSNFRFGENHAIMGVAFTWVMA"
                "LACAAPPLAGWSR--YIPEGLQC"
                "SCGI",
            )
            self.assertEqual(alignments[9][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[9][0].annotations["original_length"], 348)
            self.assertEqual(
                alignments[9][1].seq,
                "LAGEIALWSLVVLAIERYVVVCK"
                "PMSNFRFGENHAIMGLALTWVMA"
                "LACAAPPLVGWSR*WH*TEG-KC"
                "L*GL",
            )
            self.assertEqual(alignments[9][1].id, "gi|283855845|gb|GQ290303.1|")
            self.assertEqual(alignments[9][1].annotations["original_length"], 4301)
            self.assertEqual(len(alignments[10]), 2)
            self.assertEqual(alignments[10].get_alignment_length(), 30)
            self.assertEqual(alignments[10][0].seq, "IMMNKQFRNCMLTTICCGKNPLGDDEASAT")
            self.assertEqual(alignments[10][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[10][0].annotations["original_length"], 348)
            self.assertEqual(alignments[10][1].seq, "MLLAFQFRNCMLTTLCCGKNPLGDDEASTT")
            self.assertEqual(alignments[10][1].id, "gi|283855845|gb|GQ290303.1|")
            self.assertEqual(alignments[10][1].annotations["original_length"], 4301)
            self.assertEqual(len(alignments[11]), 2)
            self.assertEqual(alignments[11].get_alignment_length(), 31)
            self.assertEqual(alignments[11][0].seq, "AQQQESATTQKAEKEVTRMVIIMVIAFLICW")
            self.assertEqual(alignments[11][0].id, "sp|P08100|OPSD_HUMAN")
            self.assertEqual(alignments[11][0].annotations["original_length"], 348)
            self.assertEqual(alignments[11][1].seq, "SQQIRNATTMMMTMRVTSFSAFWVVADSCCW")
            self.assertEqual(alignments[11][1].id, "gi|283855822|gb|GQ290312.1|")
            self.assertEqual(alignments[11][1].annotations["original_length"], 983)

    def test_output009(self):
        """Check output009.m10 file."""
        fasta_file = "Fasta/output009.m10"
        with open(fasta_file) as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 7)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 22)
            self.assertEqual(alignments[0][0].seq, "TGATGTTCTGTTTCTAAAACAG")
            self.assertEqual(alignments[0][0].id, "gi|255708421:1-99")
            self.assertEqual(alignments[0][0].annotations["original_length"], 99)
            self.assertEqual(alignments[0][1].seq, "TGATTTTTTTTGTCTAAAACAG")
            self.assertEqual(alignments[0][1].id, "gi|23308614|ref|NM_152952.1|")
            self.assertEqual(alignments[0][1].annotations["original_length"], 5188)
            self.assertEqual(len(alignments[1]), 2)
            self.assertEqual(alignments[1].get_alignment_length(), 14)
            self.assertEqual(alignments[1][0].seq, "AGAAGGAAAAAAAA")
            self.assertEqual(alignments[1][0].id, "gi|156718121:2361-2376")
            self.assertEqual(alignments[1][0].annotations["original_length"], 16)
            self.assertEqual(alignments[1][1].seq, "AGAACTAAAAAAAA")
            self.assertEqual(alignments[1][1].id, "gi|47271416|ref|NM_131257.2|")
            self.assertEqual(alignments[1][1].annotations["original_length"], 597)
            self.assertEqual(len(alignments[2]), 2)
            self.assertEqual(alignments[2].get_alignment_length(), 14)
            self.assertEqual(alignments[2][0].seq, "AGAAGGAAAAAAAA")
            self.assertEqual(alignments[2][0].id, "gi|156718121:2361-2376")
            self.assertEqual(alignments[2][0].annotations["original_length"], 16)
            self.assertEqual(alignments[2][1].seq, "AGAAGGTATAAAAA")
            self.assertEqual(alignments[2][1].id, "gi|332859474|ref|XM_001156938.2|")
            self.assertEqual(alignments[2][1].annotations["original_length"], 762)
            self.assertEqual(len(alignments[3]), 2)
            self.assertEqual(alignments[3].get_alignment_length(), 14)
            self.assertEqual(alignments[3][0].seq, "TTTTTTTCCTTCTT")
            self.assertEqual(alignments[3][0].id, "gi|156718121:2361-2376")
            self.assertEqual(alignments[3][0].annotations["original_length"], 16)
            self.assertEqual(alignments[3][1].seq, "TTTTTTTACATCTT")
            self.assertEqual(alignments[3][1].id, "gi|332211534|ref|XM_003254825.1|")
            self.assertEqual(alignments[3][1].annotations["original_length"], 805)
            self.assertEqual(len(alignments[4]), 2)
            self.assertEqual(alignments[4].get_alignment_length(), 14)
            self.assertEqual(alignments[4][0].seq, "AAGAAGGAAAAAAA")
            self.assertEqual(alignments[4][0].id, "gi|156718121:2361-2376")
            self.assertEqual(alignments[4][0].annotations["original_length"], 16)
            self.assertEqual(alignments[4][1].seq, "AATAAGTAAAAAAA")
            self.assertEqual(alignments[4][1].id, "gi|23308614|ref|NM_152952.1|")
            self.assertEqual(alignments[4][1].annotations["original_length"], 5188)
            self.assertEqual(len(alignments[5]), 2)
            self.assertEqual(alignments[5].get_alignment_length(), 14)
            self.assertEqual(alignments[5][0].seq, "TTTTTTTCCTTCTT")
            self.assertEqual(alignments[5][0].id, "gi|156718121:2361-2376")
            self.assertEqual(alignments[5][0].annotations["original_length"], 16)
            self.assertEqual(alignments[5][1].seq, "TTTTTTTACATCTT")
            self.assertEqual(alignments[5][1].id, "gi|297689475|ref|XM_002822130.1|")
            self.assertEqual(alignments[5][1].annotations["original_length"], 1158)
            self.assertEqual(len(alignments[6]), 2)
            self.assertEqual(alignments[6].get_alignment_length(), 14)
            self.assertEqual(alignments[6][0].seq, "AAGAAGGAAAAAAA")
            self.assertEqual(alignments[6][0].id, "gi|156718121:2361-2376")
            self.assertEqual(alignments[6][0].annotations["original_length"], 16)
            self.assertEqual(alignments[6][1].seq, "AAGAAGGTAAAAGA")
            self.assertEqual(alignments[6][1].id, "gi|297689475|ref|XM_002822130.1|")
            self.assertEqual(alignments[6][1].annotations["original_length"], 1158)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
