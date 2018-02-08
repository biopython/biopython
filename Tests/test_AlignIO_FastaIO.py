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
        fasta_file = 'Fasta/output001.m10'
        with open(fasta_file, "r") as handle:
            alignments = list(FastaIO.FastaM10Iterator(handle))
            self.assertEqual(len(alignments), 4)
            self.assertEqual(len(alignments[0]), 2)
            self.assertEqual(alignments[0].get_alignment_length(), 108)
            self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSR gi|10955263|ref|NP_052604.1| 107")
            self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKIKVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSR gi|152973457|ref|YP_001338508.1| 931")
            self.assertEqual(len(alignments[1]), 2)
            self.assertEqual(alignments[1].get_alignment_length(), 64)
            self.assertEqual(str(alignments[1][0].seq) + " " + alignments[1][0].id + " " + str(alignments[1][0].annotations["original_length"]), "AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAIT gi|10955263|ref|NP_052604.1| 107")
            self.assertEqual(str(alignments[1][1].seq) + " " + alignments[1][1].id + " " + str(alignments[1][1].annotations["original_length"]), "ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLT gi|152973588|ref|YP_001338639.1| 459")
            self.assertEqual(len(alignments[2]), 2)
            self.assertEqual(alignments[2].get_alignment_length(), 38)
            self.assertEqual(str(alignments[2][0].seq) + " " + alignments[2][0].id + " " + str(alignments[2][0].annotations["original_length"]), "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE gi|10955264|ref|NP_052605.1| 126")
            self.assertEqual(str(alignments[2][1].seq) + " " + alignments[2][1].id + " " + str(alignments[2][1].annotations["original_length"]), "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ gi|152973462|ref|YP_001338513.1| 101")
            self.assertEqual(len(alignments[3]), 2)
            self.assertEqual(alignments[3].get_alignment_length(), 43)
            self.assertEqual(str(alignments[3][0].seq) + " " + alignments[3][0].id + " " + str(alignments[3][0].annotations["original_length"]), "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC gi|10955265|ref|NP_052606.1| 346")
            self.assertEqual(str(alignments[3][1].seq) + " " + alignments[3][1].id + " " + str(alignments[3][1].annotations["original_length"]), "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC gi|152973545|ref|YP_001338596.1| 242")

    def test_output002(self):
            """Check output002.m10 file."""
            fasta_file = 'Fasta/output002.m10'
            with open(fasta_file, "r") as handle:
                alignments = list(FastaIO.FastaM10Iterator(handle))
                self.assertEqual(len(alignments), 6)
                self.assertEqual(len(alignments[0]), 2)
                self.assertEqual(alignments[0].get_alignment_length(), 88)
                self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "SGSNTRRRAISRPVR--LTAEEDQEIRKRAAECG-KTVSGFLRAAALGKKVNSLTDDRVLKEVMRLGALQKKLFIDGKRVGDREYAEV gi|10955263|ref|NP_052604.1| 107")
                self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "SQRSTRRKPENQPTRVILFNKPYDVLPQFTDEAGRKTLKEFIPVQGVYAAGRLDRDSEGLLVLTNNGALQARLTQPGKRTGKIYYVQV gi|162139799|ref|NP_309634.2| 207")
                self.assertEqual(len(alignments[1]), 2)
                self.assertEqual(alignments[1].get_alignment_length(), 53)
                self.assertEqual(str(alignments[1][0].seq) + " " + alignments[1][0].id + " " + str(alignments[1][0].annotations["original_length"]), "EIRKRAAECGKTVSGFLRAAA-LGKKV----NSLTDDRVLKEVMRLGALQKKL gi|10955263|ref|NP_052604.1| 107")
                self.assertEqual(str(alignments[1][1].seq) + " " + alignments[1][1].id + " " + str(alignments[1][1].annotations["original_length"]), "EIKPRGTSKGEAIAAFMQEAPFIGRTPVFLGDDLTDESGFAVVNRLGGMSVKI gi|15831859|ref|NP_310632.1| 266")
                self.assertEqual(len(alignments[2]), 2)
                self.assertEqual(alignments[2].get_alignment_length(), 92)
                self.assertEqual(str(alignments[2][0].seq) + " " + alignments[2][0].id + " " + str(alignments[2][0].annotations["original_length"]), "SEFFSKIESDLKKKKSKGDVFFDLIIPNG-----GKKDRYVYTSFNGEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATS gi|10955264|ref|NP_052605.1| 126")
                self.assertEqual(str(alignments[2][1].seq) + " " + alignments[2][1].id + " " + str(alignments[2][1].annotations["original_length"]), "TELNSELAKAMKVDAQRG-AFVSQVLPNSSAAKAGIKAGDVITSLNGKPISSFAALRA-QVGTMPVGSKLTLGLLRDG-KQVNVNLELQQSS gi|15829419|ref|NP_308192.1| 474")
                self.assertEqual(len(alignments[3]), 2)
                self.assertEqual(alignments[3].get_alignment_length(), 73)
                self.assertEqual(str(alignments[3][0].seq) + " " + alignments[3][0].id + " " + str(alignments[3][0].annotations["original_length"]), "FFDLIIPNGGKKDRYVYTSFNGEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATSFALKKG gi|10955264|ref|NP_052605.1| 126")
                self.assertEqual(str(alignments[3][1].seq) + " " + alignments[3][1].id + " " + str(alignments[3][1].annotations["original_length"]), "LFDLFLKNDAMHDPMVNESYC-ETFGWVSKENLARMKE---LTYKANDVLKKLFDDAGLILVDFKLEFGLYKG gi|15832592|ref|NP_311365.1| 237")
                self.assertEqual(len(alignments[4]), 2)
                self.assertEqual(alignments[4].get_alignment_length(), 63)
                self.assertEqual(str(alignments[4][0].seq) + " " + alignments[4][0].id + " " + str(alignments[4][0].annotations["original_length"]), "VDIKK-ETIESELHSKLPKSIDKIHEDIKKQLSCSLI--MKKID-VEMEDYSTYCFSALRAIE gi|10955265|ref|NP_052606.1| 346")
                self.assertEqual(str(alignments[4][1].seq) + " " + alignments[4][1].id + " " + str(alignments[4][1].annotations["original_length"]), "IDPKKIEQIARQVHESMPKGIREFGEDVEKKIRQTLQAQLTRLDLVSREEFDVQTQVLLRTRE gi|38704138|ref|NP_311957.2| 111")
                self.assertEqual(len(alignments[5]), 2)
                self.assertEqual(alignments[5].get_alignment_length(), 157)
                self.assertEqual(str(alignments[5][0].seq) + " " + alignments[5][0].id + " " + str(alignments[5][0].annotations["original_length"]), "QYIMTTSNGDRVRAKIYKRGSIQFQGKYLQIASLINDFMCSILNMKEIVEQKNKEFNVDI---KKETI-ESELHSKLPKSIDKIHEDIKKQLSCSLIMKKIDV-EMEDYSTYCFSALRA-IEGFIYQILNDVCNPSSSKNLGEYFTENKPKYIIREI gi|10955265|ref|NP_052606.1| 346")
                self.assertEqual(str(alignments[5][1].seq) + " " + alignments[5][1].id + " " + str(alignments[5][1].annotations["original_length"]), "EFIRLLSDHDQFEKDQISELTVAANALKLEVAK--NNY-----NMKYSFDTQTERRMIELIREQKDLIPEKYLHQSGIKKL-KLHED---EFSSLLVDAERQVLEGSSFVLCCGEKINSTISELLSKKITDLTHPTESFTLSEYFSYDVYEEIFKKV gi|15833861|ref|NP_312634.1| 330")

    def test_output003(self):
            """Check output003.m10 file."""
            fasta_file = 'Fasta/output003.m10'
            with open(fasta_file, "r") as handle:
                alignments = list(FastaIO.FastaM10Iterator(handle))
                self.assertEqual(len(alignments), 3)
                self.assertEqual(len(alignments[0]), 2)
                self.assertEqual(alignments[0].get_alignment_length(), 55)
                self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGGGLRLSASTKTVDQLVRIAA gi|152973837|ref|YP_001338874.1| 183")
                self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGKKVNSLTDDRVLKEVMRLGA gi|10955263|ref|NP_052604.1| 107")
                self.assertEqual(len(alignments[1]), 2)
                self.assertEqual(alignments[1].get_alignment_length(), 22)
                self.assertEqual(str(alignments[1][0].seq) + " " + alignments[1][0].id + " " + str(alignments[1][0].annotations["original_length"]), "DDAEHLFRTLSSR-LDALQDGN gi|152973840|ref|YP_001338877.1| 63")
                self.assertEqual(str(alignments[1][1].seq) + " " + alignments[1][1].id + " " + str(alignments[1][1].annotations["original_length"]), "DDRANLFEFLSEEGITITEDNN gi|10955265|ref|NP_052606.1| 346")
                self.assertEqual(len(alignments[2]), 2)
                self.assertEqual(alignments[2].get_alignment_length(), 63)
                self.assertEqual(str(alignments[2][0].seq) + " " + alignments[2][0].id + " " + str(alignments[2][0].annotations["original_length"]), "VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIRHLKQEIEHLNTLLLSKDSHIDSLKQAM gi|152973841|ref|YP_001338878.1| 133")
                self.assertEqual(str(alignments[2][1].seq) + " " + alignments[2][1].id + " " + str(alignments[2][1].annotations["original_length"]), "VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATSF-ALKKGI gi|10955264|ref|NP_052605.1| 126")

    def test_output004(self):
            """Check output004.m10 file."""
            fasta_file = 'Fasta/output004.m10'
            with open(fasta_file, "r") as handle:
                alignments = list(FastaIO.FastaM10Iterator(handle))
                self.assertEqual(len(alignments), 1)
                self.assertEqual(len(alignments[0]), 2)
                self.assertEqual(alignments[0].get_alignment_length(), 102)
                self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "AAAAAAGATAAAAAATATCAAATAGAAGCAATAAAAAATAAAGATAAAACTTTATTTATTGTCTATGCTACTGATATTTATAGCCCGAGCGAATTTTTCTCA ref|NC_002127.1|:c1351-971 381")
                self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "AGAGAAAATAAAACAAGTAATAAAATATTAATGGAAAAAATAAATTCTTGTTTATTTAGACCTGATTCTAATCACTTTTCTTGCCCGGAGTCATTTTTGACA ref|NC_002695.1|:1970775-1971404 630")

    def test_output005(self):
            """Check output005.m10 file."""
            fasta_file = 'Fasta/output005.m10'
            with open(fasta_file, "r") as handle:
                alignments = list(FastaIO.FastaM10Iterator(handle))
                self.assertEqual(len(alignments), 1)
                self.assertEqual(len(alignments[0]), 2)
                self.assertEqual(alignments[0].get_alignment_length(), 110)
                self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "IKNKDKTLFIVYAT-DIYSPSEFFSKIESDLKKKKSKGDV--FFDLIIPNGGKKD--RYVYTSFNGEKFSSYTLNKVTKTDEYNDL--SELSASFFKKNFDKINVNLLSK gi|10955264|ref|NP_052605.1| 126")
                self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "IKDELPVAFCSWASLDLECEVKYINDVTSLYAKDWMSGERKWFIDWIAPFGHNMELYKYMRKKYPYELFRAIRLDESSKTGKIAEFHGGGIDKKLASKIFRQYHHELMSE gi|10955282|ref|NP_052623.1| 163")

    def test_output006(self):
            """Check output006.m10 file."""
            fasta_file = 'Fasta/output006.m10'
            with open(fasta_file, "r") as handle:
                alignments = list(FastaIO.FastaM10Iterator(handle))
                self.assertEqual(len(alignments), 1)
                self.assertEqual(len(alignments[0]), 2)
                self.assertEqual(alignments[0].get_alignment_length(), 131)
                self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACAACGATTAATGAGGAGATTTATGAAGAGGGTTCTTCGATTTTAGGCCAATCGGAAGGAATTATGTAGCAAGTCCATCAGAAAATGGAAGAAGTCAT query 131")
                self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACAACGATTAATGAGGAGATTTATGAAGAGGGTTCTTCGATTTTAGGCCAATCGGAAGGAATTATGTAGCAAGTCCATCAGAAAATGGAAGTAGTCAT gi|116660610|gb|EG558221.1|EG558221 573")

    def test_output007(self):
            """Check output007.m10 file."""
            fasta_file = 'Fasta/output007.m10'
            with open(fasta_file, "r") as handle:
                alignments = list(FastaIO.FastaM10Iterator(handle))
                self.assertEqual(len(alignments), 9)
                self.assertEqual(len(alignments[0]), 2)
                self.assertEqual(alignments[0].get_alignment_length(), 108)
                self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSR gi|10955263|ref|NP_052604.1| 107")
                self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKIKVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSR gi|152973457|ref|YP_001338508.1| 931")
                self.assertEqual(len(alignments[1]), 2)
                self.assertEqual(alignments[1].get_alignment_length(), 64)
                self.assertEqual(str(alignments[1][0].seq) + " " + alignments[1][0].id + " " + str(alignments[1][0].annotations["original_length"]), "AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAIT gi|10955263|ref|NP_052604.1| 107")
                self.assertEqual(str(alignments[1][1].seq) + " " + alignments[1][1].id + " " + str(alignments[1][1].annotations["original_length"]), "ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLT gi|152973588|ref|YP_001338639.1| 459")
                self.assertEqual(len(alignments[2]), 2)
                self.assertEqual(alignments[2].get_alignment_length(), 45)
                self.assertEqual(str(alignments[2][0].seq) + " " + alignments[2][0].id + " " + str(alignments[2][0].annotations["original_length"]), "EIRKRAAECGKTVSGFLRAAA-----LGKKVNSLTDDRVLKEVMR gi|10955263|ref|NP_052604.1| 107")
                self.assertEqual(str(alignments[2][1].seq) + " " + alignments[2][1].id + " " + str(alignments[2][1].annotations["original_length"]), "ELVKLIADMGISVRALLRKNVEPYEELGLEEDKFTDDQLIDFMLQ gi|152973480|ref|YP_001338531.1| 141")
                self.assertEqual(len(alignments[3]), 2)
                self.assertEqual(alignments[3].get_alignment_length(), 38)
                self.assertEqual(str(alignments[3][0].seq) + " " + alignments[3][0].id + " " + str(alignments[3][0].annotations["original_length"]), "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE gi|10955264|ref|NP_052605.1| 126")
                self.assertEqual(str(alignments[3][1].seq) + " " + alignments[3][1].id + " " + str(alignments[3][1].annotations["original_length"]), "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ gi|152973462|ref|YP_001338513.1| 101")
                self.assertEqual(len(alignments[4]), 2)
                self.assertEqual(alignments[4].get_alignment_length(), 11)
                self.assertEqual(str(alignments[4][0].seq) + " " + alignments[4][0].id + " " + str(alignments[4][0].annotations["original_length"]), "FFDLIIPNGGK gi|10955264|ref|NP_052605.1| 126")
                self.assertEqual(str(alignments[4][1].seq) + " " + alignments[4][1].id + " " + str(alignments[4][1].annotations["original_length"]), "FFDLVIENPGK gi|152973509|ref|YP_001338560.1| 448")
                self.assertEqual(len(alignments[5]), 2)
                self.assertEqual(alignments[5].get_alignment_length(), 40)
                self.assertEqual(str(alignments[5][0].seq) + " " + alignments[5][0].id + " " + str(alignments[5][0].annotations["original_length"]), "DKTLFIVYATDIYSPSE-FFSKIESDLKKKKSKGD-VFFD gi|10955264|ref|NP_052605.1| 126")
                self.assertEqual(str(alignments[5][1].seq) + " " + alignments[5][1].id + " " + str(alignments[5][1].annotations["original_length"]), "ESVVFILMAGFAMSVCYLFFSVLEKVINARKSKDESIYHD gi|152973581|ref|YP_001338632.1| 84")
                self.assertEqual(len(alignments[6]), 2)
                self.assertEqual(alignments[6].get_alignment_length(), 30)
                self.assertEqual(str(alignments[6][0].seq) + " " + alignments[6][0].id + " " + str(alignments[6][0].annotations["original_length"]), "ASFFKKNFDKINVNLLSKATSFALKKGIPI gi|10955264|ref|NP_052605.1| 126")
                self.assertEqual(str(alignments[6][1].seq) + " " + alignments[6][1].id + " " + str(alignments[6][1].annotations["original_length"]), "ASFSKEEQDKVAVDKVAADVAWQERMNKPV gi|152973536|ref|YP_001338587.1| 84")
                self.assertEqual(len(alignments[7]), 2)
                self.assertEqual(alignments[7].get_alignment_length(), 43)
                self.assertEqual(str(alignments[7][0].seq) + " " + alignments[7][0].id + " " + str(alignments[7][0].annotations["original_length"]), "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC gi|10955265|ref|NP_052606.1| 346")
                self.assertEqual(str(alignments[7][1].seq) + " " + alignments[7][1].id + " " + str(alignments[7][1].annotations["original_length"]), "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC gi|152973545|ref|YP_001338596.1| 242")
                self.assertEqual(len(alignments[8]), 2)
                self.assertEqual(alignments[8].get_alignment_length(), 64)
                self.assertEqual(str(alignments[8][0].seq) + " " + alignments[8][0].id + " " + str(alignments[8][0].annotations["original_length"]), "ISGTYKGIDFLIKLMPSGGNTTIGRASGQNNTYFDEIALIIKENCLY--SDTKNFEYTIPKFSD gi|10955265|ref|NP_052606.1| 346")
                self.assertEqual(str(alignments[8][1].seq) + " " + alignments[8][1].id + " " + str(alignments[8][1].annotations["original_length"]), "IDGVITAFD-LRTGMNISKDKVVAQIQGMDPVW---ISAAVPESIAYLLKDTSQFEISVPAYPD gi|152973505|ref|YP_001338556.1| 430")

    def test_output008(self):
            """Check output008.m10 file."""
            fasta_file = 'Fasta/output008.m10'
            with open(fasta_file, "r") as handle:
                alignments = list(FastaIO.FastaM10Iterator(handle))
                self.assertEqual(len(alignments), 12)
                self.assertEqual(len(alignments[0]), 2)
                self.assertEqual(alignments[0].get_alignment_length(), 65)
                self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH---HHHHHHHLLQDAYMQQYQHATQQQQML sp|Q9NSY1|BMP2K_HUMAN 1161")
                self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "IPHQLPHALRHRPAQEAAHASQLHPAQPGCGQPLHGLWRLHHHPVYLYAWILRLRGHGMQSGGLL gi|283855822|gb|GQ290312.1| 983")
                self.assertEqual(len(alignments[1]), 2)
                self.assertEqual(alignments[1].get_alignment_length(), 201)
                self.assertEqual(str(alignments[1][0].seq) + " " + alignments[1][0].id + " " + str(alignments[1][0].annotations["original_length"]), "GPEIL---LGQ-GPPQQPPQQHRVLQQLQQGDWRLQQLH-------LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH-----HHHHHH-HLLQDAYMQQYQHATQQQQMLQQQF-LMHSVYQPQPSASQYPTMMPQYQQAFFQQQMLAQHQPSQQQASPEYLTSPQEFSPALVSYTSSLPA-QVGTIMDSSYSANRS sp|Q9NSY1|BMP2K_HUMAN 1161")
                self.assertEqual(str(alignments[1][1].seq) + " " + alignments[1][1].id + " " + str(alignments[1][1].annotations["original_length"]), "GPELLRALLQQNGCGTQPLRVPTVLPG*AMAVLHAGRLHVPAHRAWLPHQLPHALRHGPAQEAAHASQLHPAQPGRG*PLHGLRWLHHHPLH/PLCMDTLSLGPQDAIWRASLPHWAVKLPCGLWWSWPLSGTWWCVSP*ATSA------LGRTMP*WASLSPGSWHWPALHPPSLVGPGTSLKACSVHAGSTTTHSSQKS gi|57163782|ref|NM_001009242.1| 1047")
                self.assertEqual(len(alignments[2]), 2)
                self.assertEqual(alignments[2].get_alignment_length(), 348)
                self.assertEqual(str(alignments[2][0].seq) + " " + alignments[2][0].id + " " + str(alignments[2][0].annotations["original_length"]), "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[2][1].seq) + " " + alignments[2][1].id + " " + str(alignments[2][1].annotations["original_length"]), "MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTLPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASTTGSKTETSQVAPA gi|57163782|ref|NM_001009242.1| 1047")
                self.assertEqual(len(alignments[3]), 2)
                self.assertEqual(alignments[3].get_alignment_length(), 348)
                self.assertEqual(str(alignments[3][0].seq) + " " + alignments[3][0].id + " " + str(alignments[3][0].annotations["original_length"]), "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[3][1].seq) + " " + alignments[3][1].id + " " + str(alignments[3][1].annotations["original_length"]), "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA gi|18148870|dbj|AB062417.1| 1047")
                self.assertEqual(len(alignments[4]), 2)
                self.assertEqual(alignments[4].get_alignment_length(), 326)
                self.assertEqual(str(alignments[4][0].seq) + " " + alignments[4][0].id + " " + str(alignments[4][0].annotations["original_length"]), "VPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASAT sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[4][1].seq) + " " + alignments[4][1].id + " " + str(alignments[4][1].annotations["original_length"]), "VPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSMHGYFVFGATGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLAFTWVMALACAAPPLAGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVVAFLICWLPYASVAFYIFTHQGSNFGPVFMTIPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASTT gi|283855822|gb|GQ290312.1| 983")
                self.assertEqual(len(alignments[5]), 2)
                self.assertEqual(alignments[5].get_alignment_length(), 354)
                self.assertEqual(str(alignments[5][0].seq) + " " + alignments[5][0].id + " " + str(alignments[5][0].annotations["original_length"]), "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEAS-ATVSKTE-----TSQVAPA sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[5][1].seq) + " " + alignments[5][1].id + " " + str(alignments[5][1].annotations["original_length"]), "MNGTEGPNFYIPMSNKTGVVRSPFEYPQYYLAEPWQYSILCAYMFLLILLGFPINFMTLYVTIQHKKLRTPLNYILLNLAFANHFMVLCGFTVTMYSSMNGYFILGATGCYVEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFSENHAVMGVAFTWIMALSCAVPPLLGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFTIPLIIIFFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVVFFLICWVPYASVAFFIFSNQGSEFGPIFMTVPAFFAKSSSIYNPVIYIMLNKQFRNCMITTLCCGKNPFGEDDASSAATSKTEASSVSSSQVSPA gi|2734705|gb|U59921.1|BBU59921 1574")
                self.assertEqual(len(alignments[6]), 2)
                self.assertEqual(alignments[6].get_alignment_length(), 347)
                self.assertEqual(str(alignments[6][0].seq) + " " + alignments[6][0].id + " " + str(alignments[6][0].annotations["original_length"]), "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGD-DEASATVSKTETSQVA sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[6][1].seq) + " " + alignments[6][1].id + " " + str(alignments[6][1].annotations["original_length"]), "MNGTEGPNFYIPMSNATGVVRSPFEYPQYYLAEPWAFSALSAYMFFLIIAGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVFGPTGCNIEGFFATLGGEIALWCLVVLAIERWMVVCKPVTNFRFGESHAIMGVMVTWTMALACALPPLFGWSRYIPEGLQCSCGIDYYTRAPGINNESFVIYMFTCHFSIPLAVISFCYGRLVCTVKEAAAQQQESETTQRAEREVTRMVVIMVISFLVCWVPYASVAWYIFTHQGSTFGPIFMTIPSFFAKSSALYNPMIYICMNKQFRHCMITTLCCGKNPFEEEDGASATSSKTEASSVS gi|12583664|dbj|AB043817.1| 1344")
                self.assertEqual(len(alignments[7]), 2)
                self.assertEqual(alignments[7].get_alignment_length(), 111)
                self.assertEqual(str(alignments[7][0].seq) + " " + alignments[7][0].id + " " + str(alignments[7][0].annotations["original_length"]), "VPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGG sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[7][1].seq) + " " + alignments[7][1].id + " " + str(alignments[7][1].annotations["original_length"]), "VPFSNKTGVVRSPFEHPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGG gi|283855845|gb|GQ290303.1| 4301")
                self.assertEqual(len(alignments[8]), 2)
                self.assertEqual(alignments[8].get_alignment_length(), 172)
                self.assertEqual(str(alignments[8][0].seq) + " " + alignments[8][0].id + " " + str(alignments[8][0].annotations["original_length"]), "RYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKE------------------------------------AAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQ sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[8][1].seq) + " " + alignments[8][1].id + " " + str(alignments[8][1].annotations["original_length"]), "RYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEVRSCVGHWGHAH*VNGAQLHSQSCHSLDT*PCVPA\AAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSNFGPIFMTLPAFFAKSSSIYNPVIYIMMNKQ gi|283855845|gb|GQ290303.1| 4301")
                self.assertEqual(len(alignments[9]), 2)
                self.assertEqual(alignments[9].get_alignment_length(), 73)
                self.assertEqual(str(alignments[9][0].seq) + " " + alignments[9][0].id + " " + str(alignments[9][0].annotations["original_length"]), "LGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSR--YIPEGLQCSCGI sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[9][1].seq) + " " + alignments[9][1].id + " " + str(alignments[9][1].annotations["original_length"]), "LAGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWVMALACAAPPLVGWSR*WH*TEG-KCL*GL gi|283855845|gb|GQ290303.1| 4301")
                self.assertEqual(len(alignments[10]), 2)
                self.assertEqual(alignments[10].get_alignment_length(), 30)
                self.assertEqual(str(alignments[10][0].seq) + " " + alignments[10][0].id + " " + str(alignments[10][0].annotations["original_length"]), "IMMNKQFRNCMLTTICCGKNPLGDDEASAT sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[10][1].seq) + " " + alignments[10][1].id + " " + str(alignments[10][1].annotations["original_length"]), "MLLAFQFRNCMLTTLCCGKNPLGDDEASTT gi|283855845|gb|GQ290303.1| 4301")
                self.assertEqual(len(alignments[11]), 2)
                self.assertEqual(alignments[11].get_alignment_length(), 31)
                self.assertEqual(str(alignments[11][0].seq) + " " + alignments[11][0].id + " " + str(alignments[11][0].annotations["original_length"]), "AQQQESATTQKAEKEVTRMVIIMVIAFLICW sp|P08100|OPSD_HUMAN 348")
                self.assertEqual(str(alignments[11][1].seq) + " " + alignments[11][1].id + " " + str(alignments[11][1].annotations["original_length"]), "SQQIRNATTMMMTMRVTSFSAFWVVADSCCW gi|283855822|gb|GQ290312.1| 983")

    def test_output009(self):
            """Check output009.m10 file."""
            fasta_file = 'Fasta/output009.m10'
            with open(fasta_file, "r") as handle:
                alignments = list(FastaIO.FastaM10Iterator(handle))
                self.assertEqual(len(alignments), 7)
                self.assertEqual(len(alignments[0]), 2)
                self.assertEqual(alignments[0].get_alignment_length(), 22)
                self.assertEqual(str(alignments[0][0].seq) + " " + alignments[0][0].id + " " + str(alignments[0][0].annotations["original_length"]), "TGATGTTCTGTTTCTAAAACAG gi|255708421:1-99 99")
                self.assertEqual(str(alignments[0][1].seq) + " " + alignments[0][1].id + " " + str(alignments[0][1].annotations["original_length"]), "TGATTTTTTTTGTCTAAAACAG gi|23308614|ref|NM_152952.1| 5188")
                self.assertEqual(len(alignments[1]), 2)
                self.assertEqual(alignments[1].get_alignment_length(), 14)
                self.assertEqual(str(alignments[1][0].seq) + " " + alignments[1][0].id + " " + str(alignments[1][0].annotations["original_length"]), "AGAAGGAAAAAAAA gi|156718121:2361-2376 16")
                self.assertEqual(str(alignments[1][1].seq) + " " + alignments[1][1].id + " " + str(alignments[1][1].annotations["original_length"]), "AGAACTAAAAAAAA gi|47271416|ref|NM_131257.2| 597")
                self.assertEqual(len(alignments[2]), 2)
                self.assertEqual(alignments[2].get_alignment_length(), 14)
                self.assertEqual(str(alignments[2][0].seq) + " " + alignments[2][0].id + " " + str(alignments[2][0].annotations["original_length"]), "AGAAGGAAAAAAAA gi|156718121:2361-2376 16")
                self.assertEqual(str(alignments[2][1].seq) + " " + alignments[2][1].id + " " + str(alignments[2][1].annotations["original_length"]), "AGAAGGTATAAAAA gi|332859474|ref|XM_001156938.2| 762")
                self.assertEqual(len(alignments[3]), 2)
                self.assertEqual(alignments[3].get_alignment_length(), 14)
                self.assertEqual(str(alignments[3][0].seq) + " " + alignments[3][0].id + " " + str(alignments[3][0].annotations["original_length"]), "TTTTTTTCCTTCTT gi|156718121:2361-2376 16")
                self.assertEqual(str(alignments[3][1].seq) + " " + alignments[3][1].id + " " + str(alignments[3][1].annotations["original_length"]), "TTTTTTTACATCTT gi|332211534|ref|XM_003254825.1| 805")
                self.assertEqual(len(alignments[4]), 2)
                self.assertEqual(alignments[4].get_alignment_length(), 14)
                self.assertEqual(str(alignments[4][0].seq) + " " + alignments[4][0].id + " " + str(alignments[4][0].annotations["original_length"]), "AAGAAGGAAAAAAA gi|156718121:2361-2376 16")
                self.assertEqual(str(alignments[4][1].seq) + " " + alignments[4][1].id + " " + str(alignments[4][1].annotations["original_length"]), "AATAAGTAAAAAAA gi|23308614|ref|NM_152952.1| 5188")
                self.assertEqual(len(alignments[5]), 2)
                self.assertEqual(alignments[5].get_alignment_length(), 14)
                self.assertEqual(str(alignments[5][0].seq) + " " + alignments[5][0].id + " " + str(alignments[5][0].annotations["original_length"]), "TTTTTTTCCTTCTT gi|156718121:2361-2376 16")
                self.assertEqual(str(alignments[5][1].seq) + " " + alignments[5][1].id + " " + str(alignments[5][1].annotations["original_length"]), "TTTTTTTACATCTT gi|297689475|ref|XM_002822130.1| 1158")
                self.assertEqual(len(alignments[6]), 2)
                self.assertEqual(alignments[6].get_alignment_length(), 14)
                self.assertEqual(str(alignments[6][0].seq) + " " + alignments[6][0].id + " " + str(alignments[6][0].annotations["original_length"]), "AAGAAGGAAAAAAA gi|156718121:2361-2376 16")
                self.assertEqual(str(alignments[6][1].seq) + " " + alignments[6][1].id + " " + str(alignments[6][1].annotations["original_length"]), "AAGAAGGTAAAAGA gi|297689475|ref|XM_002822130.1| 1158")


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
