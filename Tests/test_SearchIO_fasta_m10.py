# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO FastaIO parsers."""

import os
import unittest

from Bio.SearchIO import parse


# test case files are in the Blast directory
TEST_DIR = 'Fasta'
FMT = 'fasta-m10'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class Fasta34Cases(unittest.TestCase):

    def test_output002(self):
        """Test parsing fasta34 output (output002.m10)"""

        m10_file = get_file('output002.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(3, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|10955263|ref|NP_052604.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('34.26', qresult.version)
        self.assertEqual('NC_002695.faa', qresult.target)
        self.assertEqual(107, qresult.seq_len)
        self.assertEqual('plasmid mobilization [Escherichia coli O157:H7 s 107 aa', qresult.description)
        self.assertEqual(2, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|162139799|ref|NP_309634.2|', hit.id)
        self.assertEqual('23S rRNA pseudouridine synthase E [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(207, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(55, hsp.initn_score)
        self.assertEqual(55, hsp.init1_score)
        self.assertEqual(77, hsp.opt_score)
        self.assertEqual(110.8, hsp.z_score)
        self.assertEqual(26.5, hsp.bitscore)
        self.assertEqual(1.2, hsp.evalue)
        self.assertEqual(77, hsp.sw_score)
        self.assertAlmostEqual(28.4, hsp.ident_pct)
        self.assertAlmostEqual(54.5, hsp.pos_pct)
        self.assertEqual(88, hsp.aln_span)
        self.assertEqual(4, hsp.query_start)
        self.assertEqual(89, hsp.query_end)
        self.assertEqual('SGSNTRRRAISRPVR--LTAEEDQEIRKRAAECG-KTVSGFLRAAALGKKVNSLTDDRVLKEVMRLGALQKKLFIDGKRVGDREYAEV', str(hsp.query.seq))
        self.assertEqual(15, hsp.hit_start)
        self.assertEqual(103, hsp.hit_end)
        self.assertEqual('SQRSTRRKPENQPTRVILFNKPYDVLPQFTDEAGRKTLKEFIPVQGVYAAGRLDRDSEGLLVLTNNGALQARLTQPGKRTGKIYYVQV', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|15831859|ref|NP_310632.1|', hit.id)
        self.assertEqual('trehalose-6-phosphate phosphatase [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(266, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(43, hsp.initn_score)
        self.assertEqual(43, hsp.init1_score)
        self.assertEqual(69, hsp.opt_score)
        self.assertEqual(98.6, hsp.z_score)
        self.assertEqual(24.6, hsp.bitscore)
        self.assertEqual(5.8, hsp.evalue)
        self.assertEqual(69, hsp.sw_score)
        self.assertAlmostEqual(28.3, hsp.ident_pct)
        self.assertAlmostEqual(66.0, hsp.pos_pct)
        self.assertEqual(53, hsp.aln_span)
        self.assertEqual(26, hsp.query_start)
        self.assertEqual(74, hsp.query_end)
        self.assertEqual('EIRKRAAECGKTVSGFLRAAA-LGKKV----NSLTDDRVLKEVMRLGALQKKL', str(hsp.query.seq))
        self.assertEqual(166, hsp.hit_start)
        self.assertEqual(219, hsp.hit_end)
        self.assertEqual('EIKPRGTSKGEAIAAFMQEAPFIGRTPVFLGDDLTDESGFAVVNRLGGMSVKI', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|10955264|ref|NP_052605.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('34.26', qresult.version)
        self.assertEqual('NC_002695.faa', qresult.target)
        self.assertEqual(126, qresult.seq_len)
        self.assertEqual('hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa', qresult.description)
        self.assertEqual(2, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|15829419|ref|NP_308192.1|', hit.id)
        self.assertEqual('serine endoprotease [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(474, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(64, hsp.initn_score)
        self.assertEqual(40, hsp.init1_score)
        self.assertEqual(77, hsp.opt_score)
        self.assertEqual(105.8, hsp.z_score)
        self.assertEqual(27.0, hsp.bitscore)
        self.assertEqual(2.3, hsp.evalue)
        self.assertEqual(77, hsp.sw_score)
        self.assertAlmostEqual(25.0, hsp.ident_pct)
        self.assertAlmostEqual(62.0, hsp.pos_pct)
        self.assertEqual(92, hsp.aln_span)
        self.assertEqual(30, hsp.query_start)
        self.assertEqual(117, hsp.query_end)
        self.assertEqual('SEFFSKIESDLKKKKSKGDVFFDLIIPNG-----GKKDRYVYTSFNGEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATS', str(hsp.query.seq))
        self.assertEqual(295, hsp.hit_start)
        self.assertEqual(384, hsp.hit_end)
        self.assertEqual('TELNSELAKAMKVDAQRG-AFVSQVLPNSSAAKAGIKAGDVITSLNGKPISSFAALRA-QVGTMPVGSKLTLGLLRDG-KQVNVNLELQQSS', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|15832592|ref|NP_311365.1|', hit.id)
        self.assertEqual('phosphoribosylaminoimidazole-succinocarboxamide synthase [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(237, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(73, hsp.initn_score)
        self.assertEqual(45, hsp.init1_score)
        self.assertEqual(74, hsp.opt_score)
        self.assertEqual(105.5, hsp.z_score)
        self.assertEqual(26.0, hsp.bitscore)
        self.assertEqual(2.4, hsp.evalue)
        self.assertEqual(74, hsp.sw_score)
        self.assertAlmostEqual(27.4, hsp.ident_pct)
        self.assertAlmostEqual(58.9, hsp.pos_pct)
        self.assertEqual(73, hsp.aln_span)
        self.assertEqual(50, hsp.query_start)
        self.assertEqual(123, hsp.query_end)
        self.assertEqual('FFDLIIPNGGKKDRYVYTSFNGEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATSFALKKG', str(hsp.query.seq))
        self.assertEqual(116, hsp.hit_start)
        self.assertEqual(185, hsp.hit_end)
        self.assertEqual('LFDLFLKNDAMHDPMVNESYC-ETFGWVSKENLARMKE---LTYKANDVLKKLFDDAGLILVDFKLEFGLYKG', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('gi|10955265|ref|NP_052606.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('34.26', qresult.version)
        self.assertEqual('NC_002695.faa', qresult.target)
        self.assertEqual(346, qresult.seq_len)
        self.assertEqual('hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa', qresult.description)
        self.assertEqual(2, len(qresult))
        # third qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|38704138|ref|NP_311957.2|', hit.id)
        self.assertEqual('hypothetical protein ECs3930 [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(111, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(50, hsp.initn_score)
        self.assertEqual(50, hsp.init1_score)
        self.assertEqual(86, hsp.opt_score)
        self.assertEqual(117.5, hsp.z_score)
        self.assertEqual(28.6, hsp.bitscore)
        self.assertEqual(0.51, hsp.evalue)
        self.assertEqual(86, hsp.sw_score)
        self.assertAlmostEqual(30.2, hsp.ident_pct)
        self.assertAlmostEqual(63.5, hsp.pos_pct)
        self.assertEqual(63, hsp.aln_span)
        self.assertEqual(187, hsp.query_start)
        self.assertEqual(246, hsp.query_end)
        self.assertEqual('VDIKK-ETIESELHSKLPKSIDKIHEDIKKQLSCSLI--MKKID-VEMEDYSTYCFSALRAIE', str(hsp.query.seq))
        self.assertEqual(13, hsp.hit_start)
        self.assertEqual(76, hsp.hit_end)
        self.assertEqual('IDPKKIEQIARQVHESMPKGIREFGEDVEKKIRQTLQAQLTRLDLVSREEFDVQTQVLLRTRE', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)
        # third qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|15833861|ref|NP_312634.1|', hit.id)
        self.assertEqual('hypothetical protein ECs4607 [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(330, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(32, hsp.initn_score)
        self.assertEqual(32, hsp.init1_score)
        self.assertEqual(87, hsp.opt_score)
        self.assertEqual(112.7, hsp.z_score)
        self.assertEqual(29.2, hsp.bitscore)
        self.assertEqual(0.95, hsp.evalue)
        self.assertEqual(87, hsp.sw_score)
        self.assertAlmostEqual(21.0, hsp.ident_pct)
        self.assertAlmostEqual(58.0, hsp.pos_pct)
        self.assertEqual(157, hsp.aln_span)
        self.assertEqual(130, hsp.query_start)
        self.assertEqual(281, hsp.query_end)
        self.assertEqual('QYIMTTSNGDRVRAKIYKRGSIQFQGKYLQIASLINDFMCSILNMKEIVEQKNKEFNVDI---KKETI-ESELHSKLPKSIDKIHEDIKKQLSCSLIMKKIDV-EMEDYSTYCFSALRA-IEGFIYQILNDVCNPSSSKNLGEYFTENKPKYIIREI', str(hsp.query.seq))
        self.assertEqual(9, hsp.hit_start)
        self.assertEqual(155, hsp.hit_end)
        self.assertEqual('EFIRLLSDHDQFEKDQISELTVAANALKLEVAK--NNY-----NMKYSFDTQTERRMIELIREQKDLIPEKYLHQSGIKKL-KLHED---EFSSLLVDAERQVLEGSSFVLCCGEKINSTISELLSKKITDLTHPTESFTLSEYFSYDVYEEIFKKV', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

    def test_output003(self):
        """Test parsing fasta34 output (output003.m10)"""

        m10_file = get_file('output003.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(5, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|152973837|ref|YP_001338874.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('34.26', qresult.version)
        self.assertEqual('NC_002127.faa', qresult.target)
        self.assertEqual(183, qresult.seq_len)
        self.assertEqual('hypothetical protein KPN_pKPN7p10262 [Klebsiella pneumoniae subsp. pneumonia 183 aa', qresult.description)
        self.assertEqual(1, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|10955263|ref|NP_052604.1|', hit.id)
        self.assertEqual('plasmid mobilization [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(107, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(43, hsp.initn_score)
        self.assertEqual(43, hsp.init1_score)
        self.assertEqual(45, hsp.opt_score)
        self.assertEqual(64.1, hsp.z_score)
        self.assertEqual(17.7, hsp.bitscore)
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(59, hsp.sw_score)
        self.assertAlmostEqual(25.5, hsp.ident_pct)
        self.assertAlmostEqual(67.3, hsp.pos_pct)
        self.assertEqual(55, hsp.aln_span)
        self.assertEqual(86, hsp.query_start)
        self.assertEqual(141, hsp.query_end)
        self.assertEqual('ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGGGLRLSASTKTVDQLVRIAA', str(hsp.query.seq))
        self.assertEqual(17, hsp.hit_start)
        self.assertEqual(69, hsp.hit_end)
        self.assertEqual('VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGKKVNSLTDDRVLKEVMRLGA', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|152973838|ref|YP_001338875.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('34.26', qresult.version)
        self.assertEqual('NC_002127.faa', qresult.target)
        self.assertEqual(76, qresult.seq_len)
        self.assertEqual('hypothetical protein KPN_pKPN7p10263 [Klebsiella pneumoniae subsp. pneumonia 76 aa', qresult.description)
        self.assertEqual(0, len(qresult))

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('gi|152973839|ref|YP_001338876.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('34.26', qresult.version)
        self.assertEqual('NC_002127.faa', qresult.target)
        self.assertEqual(112, qresult.seq_len)
        self.assertEqual('hypothetical protein KPN_pKPN7p10264 [Klebsiella pneumoniae subsp. pneumonia 112 aa', qresult.description)
        self.assertEqual(0, len(qresult))

        # test fourth qresult
        qresult = qresults[3]
        self.assertEqual('gi|152973840|ref|YP_001338877.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('34.26', qresult.version)
        self.assertEqual('NC_002127.faa', qresult.target)
        self.assertEqual(63, qresult.seq_len)
        self.assertEqual('RNA one modulator-like protein [Klebsiella pneumoniae subsp. pneumoniae  63 aa', qresult.description)
        self.assertEqual(1, len(qresult))
        # fourth qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|10955265|ref|NP_052606.1|', hit.id)
        self.assertEqual('hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(346, hit.seq_len)
        self.assertEqual(1, len(hit))
        # fourth qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(35, hsp.initn_score)
        self.assertEqual(35, hsp.init1_score)
        self.assertEqual(38, hsp.opt_score)
        self.assertEqual(71.3, hsp.z_score)
        self.assertEqual(19.2, hsp.bitscore)
        self.assertEqual(0.11, hsp.evalue)
        self.assertEqual(38, hsp.sw_score)
        self.assertAlmostEqual(36.4, hsp.ident_pct)
        self.assertAlmostEqual(63.6, hsp.pos_pct)
        self.assertEqual(22, hsp.aln_span)
        self.assertEqual(42, hsp.query_start)
        self.assertEqual(63, hsp.query_end)
        self.assertEqual('DDAEHLFRTLSSR-LDALQDGN', str(hsp.query.seq))
        self.assertEqual(101, hsp.hit_start)
        self.assertEqual(123, hsp.hit_end)
        self.assertEqual('DDRANLFEFLSEEGITITEDNN', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

        # test fifth qresult
        qresult = qresults[4]
        self.assertEqual('gi|152973841|ref|YP_001338878.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('34.26', qresult.version)
        self.assertEqual('NC_002127.faa', qresult.target)
        self.assertEqual(133, qresult.seq_len)
        self.assertEqual('Excl1 protein [Klebsiella pneumoniae subsp. pneumoniae  133 aa', qresult.description)
        self.assertEqual(1, len(qresult))
        # fifth qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|10955264|ref|NP_052605.1|', hit.id)
        self.assertEqual('hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(126, hit.seq_len)
        self.assertEqual(1, len(hit))
        # fifth qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(37, hsp.initn_score)
        self.assertEqual(37, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(80.0, hsp.z_score)
        self.assertEqual(20.4, hsp.bitscore)
        self.assertEqual(0.036, hsp.evalue)
        self.assertEqual(57, hsp.sw_score)
        self.assertAlmostEqual(25.4, hsp.ident_pct)
        self.assertAlmostEqual(65.1, hsp.pos_pct)
        self.assertEqual(63, hsp.aln_span)
        self.assertEqual(48, hsp.query_start)
        self.assertEqual(109, hsp.query_end)
        self.assertEqual('VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIRHLKQEIEHLNTLLLSKDSHIDSLKQAM', str(hsp.query.seq))
        self.assertEqual(65, hsp.hit_start)
        self.assertEqual(124, hsp.hit_end)
        self.assertEqual('VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDKINVNLLSKATSF-ALKKGI', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)


class Fasta35Cases(unittest.TestCase):

    def test_output001(self):
        """Test parsing fasta35 output (output001.m10)"""

        m10_file = get_file('output001.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(3, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|10955263|ref|NP_052604.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('35.03', qresult.version)
        self.assertEqual('NC_009649.faa', qresult.target)
        self.assertEqual(107, qresult.seq_len)
        self.assertEqual('plasmid mobilization [Escherichia coli O157:H7 s 107 aa', qresult.description)
        self.assertEqual(2, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|152973457|ref|YP_001338508.1|', hit.id)
        self.assertEqual('ATPase with chaperone activity, ATP-binding subunit [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(931, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(65, hsp.initn_score)
        self.assertEqual(43, hsp.init1_score)
        self.assertEqual(71, hsp.opt_score)
        self.assertEqual(92.7, hsp.z_score)
        self.assertEqual(25.3, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(71, hsp.sw_score)
        self.assertAlmostEqual(25.0, hsp.ident_pct)
        self.assertAlmostEqual(57.4, hsp.pos_pct)
        self.assertEqual(108, hsp.aln_span)
        self.assertEqual(4, hsp.query_start)
        self.assertEqual(103, hsp.query_end)
        self.assertEqual('SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSR', str(hsp.query.seq))
        self.assertEqual(95, hsp.hit_start)
        self.assertEqual(195, hsp.hit_end)
        self.assertEqual('AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKIKVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSR', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|152973588|ref|YP_001338639.1|', hit.id)
        self.assertEqual('F pilus assembly protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(459, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(33, hsp.initn_score)
        self.assertEqual(33, hsp.init1_score)
        self.assertEqual(63, hsp.opt_score)
        self.assertEqual(87.7, hsp.z_score)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.81, hsp.evalue)
        self.assertEqual(63, hsp.sw_score)
        self.assertAlmostEqual(26.6, hsp.ident_pct)
        self.assertAlmostEqual(65.6, hsp.pos_pct)
        self.assertEqual(64, hsp.aln_span)
        self.assertEqual(31, hsp.query_start)
        self.assertEqual(94, hsp.query_end)
        self.assertEqual('AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAIT', str(hsp.query.seq))
        self.assertEqual(190, hsp.hit_start)
        self.assertEqual(248, hsp.hit_end)
        self.assertEqual('ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLT', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|10955264|ref|NP_052605.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('35.03', qresult.version)
        self.assertEqual('NC_009649.faa', qresult.target)
        self.assertEqual(126, qresult.seq_len)
        self.assertEqual('hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa', qresult.description)
        self.assertEqual(1, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|152973462|ref|YP_001338513.1|', hit.id)
        self.assertEqual('hypothetical protein KPN_pKPN3p05904 [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(101, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(50, hsp.initn_score)
        self.assertEqual(50, hsp.init1_score)
        self.assertEqual(58, hsp.opt_score)
        self.assertEqual(91.6, hsp.z_score)
        self.assertEqual(22.2, hsp.bitscore)
        self.assertEqual(0.49, hsp.evalue)
        self.assertEqual(58, hsp.sw_score)
        self.assertAlmostEqual(28.9, hsp.ident_pct)
        self.assertAlmostEqual(63.2, hsp.pos_pct)
        self.assertEqual(38, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(38, hsp.query_end)
        self.assertEqual('MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE', str(hsp.query.seq))
        self.assertEqual(43, hsp.hit_start)
        self.assertEqual(81, hsp.hit_end)
        self.assertEqual('IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('gi|10955265|ref|NP_052606.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('35.03', qresult.version)
        self.assertEqual('NC_009649.faa', qresult.target)
        self.assertEqual(346, qresult.seq_len)
        self.assertEqual('hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa', qresult.description)
        self.assertEqual(1, len(qresult))
        # third qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|152973545|ref|YP_001338596.1|', hit.id)
        self.assertEqual('putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(242, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(52, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(70, hsp.opt_score)
        self.assertEqual(94.0, hsp.z_score)
        self.assertEqual(25.3, hsp.bitscore)
        self.assertEqual(0.36, hsp.evalue)
        self.assertEqual(70, hsp.sw_score)
        self.assertAlmostEqual(27.9, hsp.ident_pct)
        self.assertAlmostEqual(65.1, hsp.pos_pct)
        self.assertEqual(43, hsp.aln_span)
        self.assertEqual(196, hsp.query_start)
        self.assertEqual(238, hsp.query_end)
        self.assertEqual('SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC', str(hsp.query.seq))
        self.assertEqual(51, hsp.hit_start)
        self.assertEqual(94, hsp.hit_end)
        self.assertEqual('SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

    def test_output004(self):
        """Test parsing fasta35 output (output004.m10)"""

        m10_file = get_file('output004.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(3, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('ref|NC_002127.1|:413-736', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('35.04', qresult.version)
        self.assertEqual('NC_002695.ffn', qresult.target)
        self.assertEqual(324, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('ref|NC_002127.1|:c1351-971', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('35.04', qresult.version)
        self.assertEqual('NC_002695.ffn', qresult.target)
        self.assertEqual(381, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(1, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|NC_002695.1|:1970775-1971404', hit.id)
        self.assertEqual('', hit.description)
        self.assertEqual(630, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(54, hsp.initn_score)
        self.assertEqual(54, hsp.init1_score)
        self.assertEqual(91, hsp.opt_score)
        self.assertEqual(139.3, hsp.z_score)
        self.assertEqual(35.2, hsp.bitscore)
        self.assertEqual(0.045, hsp.evalue)
        self.assertAlmostEqual(57.8, hsp.ident_pct)
        self.assertAlmostEqual(57.8, hsp.pos_pct)
        self.assertEqual(102, hsp.aln_span)
        self.assertEqual(3, hsp.query_start)
        self.assertEqual(105, hsp.query_end)
        self.assertEqual('AAAAAAGATAAAAAATATCAAATAGAAGCAATAAAAAATAAAGATAAAACTTTATTTATTGTCTATGCTACTGATATTTATAGCCCGAGCGAATTTTTCTCA', str(hsp.query.seq))
        self.assertEqual(312, hsp.hit_start)
        self.assertEqual(414, hsp.hit_end)
        self.assertEqual('AGAGAAAATAAAACAAGTAATAAAATATTAATGGAAAAAATAAATTCTTGTTTATTTAGACCTGATTCTAATCACTTTTCTTGCCCGGAGTCATTTTTGACA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual({'similarity': ': : :: :::::: :  : : : :  :  :::  :::: : : ::     ::::::::      :: ::: : :  ::: : :::::     ::::::  ::'}, hsp.aln_annotation)

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('ref|NC_002127.1|:c2388-1348', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('35.04', qresult.version)
        self.assertEqual('NC_002695.ffn', qresult.target)
        self.assertEqual(1041, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(0, len(qresult))

    def test_output005(self):
        """Test parsing ssearch35 output (output005.m10)"""

        m10_file = get_file('output005.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(3, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|10955263|ref|NP_052604.1|', qresult.id)
        self.assertEqual('ssearch', qresult.program)
        self.assertEqual('35.03', qresult.version)
        self.assertEqual('NC_002128.faa', qresult.target)
        self.assertEqual(107, qresult.seq_len)
        self.assertEqual('plasmid mobilization [Escherichia coli O157:H7 s 107 aa', qresult.description)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|10955264|ref|NP_052605.1|', qresult.id)
        self.assertEqual('ssearch', qresult.program)
        self.assertEqual('35.03', qresult.version)
        self.assertEqual('NC_002128.faa', qresult.target)
        self.assertEqual(126, qresult.seq_len)
        self.assertEqual('hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa', qresult.description)
        self.assertEqual(1, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|10955282|ref|NP_052623.1|', hit.id)
        self.assertEqual('hemolysin C [Escherichia coli O157:H7 str. Sakai]', hit.description)
        self.assertEqual(163, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(69, hsp.opt_score)
        self.assertEqual(108.8, hsp.z_score)
        self.assertEqual(26.0, hsp.bitscore)
        self.assertEqual(0.025, hsp.evalue)
        self.assertEqual(69, hsp.sw_score)
        self.assertAlmostEqual(20.9, hsp.ident_pct)
        self.assertAlmostEqual(55.5, hsp.pos_pct)
        self.assertEqual(110, hsp.aln_span)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(114, hsp.query_end)
        self.assertEqual('IKNKDKTLFIVYAT-DIYSPSEFFSKIESDLKKKKSKGDV--FFDLIIPNGGKKD--RYVYTSFNGEKFSSYTLNKVTKTDEYNDL--SELSASFFKKNFDKINVNLLSK', str(hsp.query.seq))
        self.assertEqual(38, hsp.hit_start)
        self.assertEqual(148, hsp.hit_end)
        self.assertEqual('IKDELPVAFCSWASLDLECEVKYINDVTSLYAKDWMSGERKWFIDWIAPFGHNMELYKYMRKKYPYELFRAIRLDESSKTGKIAEFHGGGIDKKLASKIFRQYHHELMSE', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({}, hsp.aln_annotation)

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('gi|10955265|ref|NP_052606.1|', qresult.id)
        self.assertEqual('ssearch', qresult.program)
        self.assertEqual('35.03', qresult.version)
        self.assertEqual('NC_002128.faa', qresult.target)
        self.assertEqual(346, qresult.seq_len)
        self.assertEqual('hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa', qresult.description)
        self.assertEqual(0, len(qresult))

    def test_output006(self):
        """Test parsing fasta35 output (output006.m10)"""

        m10_file = get_file('output006.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(1, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('query', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('35.04', qresult.version)
        self.assertEqual('orchid_cds.txt', qresult.target)
        self.assertEqual(131, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(1, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|116660610|gb|EG558221.1|EG558221', hit.id)
        self.assertEqual('CR03001A07 Root CR03 cDNA library Catharanthus roseus cDNA clone CR03001A07 5\', mRNA sequence', hit.description)
        self.assertEqual(573, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(646, hsp.initn_score)
        self.assertEqual(646, hsp.init1_score)
        self.assertEqual(646, hsp.opt_score)
        self.assertEqual(712.3, hsp.z_score)
        self.assertEqual(139.6, hsp.bitscore)
        self.assertEqual(7.2e-38, hsp.evalue)
        self.assertAlmostEqual(99.2, hsp.ident_pct)
        self.assertAlmostEqual(99.2, hsp.pos_pct)
        self.assertEqual(131, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(131, hsp.query_end)
        self.assertEqual('GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACAACGATTAATGAGGAGATTTATGAAGAGGGTTCTTCGATTTTAGGCCAATCGGAAGGAATTATGTAGCAAGTCCATCAGAAAATGGAAGAAGTCAT', str(hsp.query.seq))
        self.assertEqual(359, hsp.hit_start)
        self.assertEqual(490, hsp.hit_end)
        self.assertEqual('GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACAACGATTAATGAGGAGATTTATGAAGAGGGTTCTTCGATTTTAGGCCAATCGGAAGGAATTATGTAGCAAGTCCATCAGAAAATGGAAGTAGTCAT', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual({'similarity': ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ::::::'}, hsp.aln_annotation)


class Fasta36Cases(unittest.TestCase):

    def test_output007(self):
        """Test parsing fasta36 output (output007.m10)"""

        m10_file = get_file('output007.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(3, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|10955263|ref|NP_052604.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.4', qresult.version)
        self.assertEqual('NC_009649.faa', qresult.target)
        self.assertEqual(107, qresult.seq_len)
        self.assertEqual('plasmid mobilization [Escherichia coli O157:H7 s', qresult.description)
        self.assertEqual(3, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|152973457|ref|YP_001338508.1|', hit.id)
        self.assertEqual('ATPase with chaperone activity, ATP-binding subunit [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(931, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(97, hsp.initn_score)
        self.assertEqual(43, hsp.init1_score)
        self.assertEqual(71, hsp.opt_score)
        self.assertEqual(109.6, hsp.z_score)
        self.assertEqual(28.5, hsp.bitscore)
        self.assertEqual(0.048, hsp.evalue)
        self.assertEqual(71, hsp.sw_score)
        self.assertAlmostEqual(25.0, hsp.ident_pct)
        self.assertAlmostEqual(57.4, hsp.pos_pct)
        self.assertEqual(108, hsp.aln_span)
        self.assertEqual(4, hsp.query_start)
        self.assertEqual(103, hsp.query_end)
        self.assertEqual('SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSR', str(hsp.query.seq))
        self.assertEqual(95, hsp.hit_start)
        self.assertEqual(195, hsp.hit_end)
        self.assertEqual('AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKIKVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSR', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': '.::..-:::. .   :.. . .---::  :.::: :..------ .   . . .:.:. :.: ..-----..  :....  ..::-::. .-:  :...   .  :::'}, hsp.aln_annotation)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|152973588|ref|YP_001338639.1|', hit.id)
        self.assertEqual('F pilus assembly protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(459, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(66, hsp.initn_score)
        self.assertEqual(33, hsp.init1_score)
        self.assertEqual(63, hsp.opt_score)
        self.assertEqual(101.4, hsp.z_score)
        self.assertEqual(25.9, hsp.bitscore)
        self.assertEqual(0.14, hsp.evalue)
        self.assertEqual(63, hsp.sw_score)
        self.assertAlmostEqual(26.6, hsp.ident_pct)
        self.assertAlmostEqual(65.6, hsp.pos_pct)
        self.assertEqual(64, hsp.aln_span)
        self.assertEqual(31, hsp.query_start)
        self.assertEqual(94, hsp.query_end)
        self.assertEqual('AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAIT', str(hsp.query.seq))
        self.assertEqual(190, hsp.hit_start)
        self.assertEqual(248, hsp.hit_end)
        self.assertEqual('ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLT', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': ':.. : ::.:--.  ..  :...   .::.:..-.  .::.:. ..::----..:  : ....:'}, hsp.aln_annotation)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|152973480|ref|YP_001338531.1|', hit.id)
        self.assertEqual('Arsenate reductase (Arsenical pump modifier) [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(141, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(45, hsp.initn_score)
        self.assertEqual(37, hsp.init1_score)
        self.assertEqual(51, hsp.opt_score)
        self.assertEqual(89.6, hsp.z_score)
        self.assertEqual(22.0, hsp.bitscore)
        self.assertEqual(0.63, hsp.evalue)
        self.assertEqual(51, hsp.sw_score)
        self.assertAlmostEqual(26.7, hsp.ident_pct)
        self.assertAlmostEqual(62.2, hsp.pos_pct)
        self.assertEqual(45, hsp.aln_span)
        self.assertEqual(26, hsp.query_start)
        self.assertEqual(66, hsp.query_end)
        self.assertEqual('EIRKRAAECGKTVSGFLRAAA-----LGKKVNSLTDDRVLKEVMR', str(hsp.query.seq))
        self.assertEqual(42, hsp.hit_start)
        self.assertEqual(87, hsp.hit_end)
        self.assertEqual('ELVKLIADMGISVRALLRKNVEPYEELGLEEDKFTDDQLIDFMLQ', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': ':. :  :. : .: ..::  .-----:: . ...:::...  ...'}, hsp.aln_annotation)

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|10955264|ref|NP_052605.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.4', qresult.version)
        self.assertEqual('NC_009649.faa', qresult.target)
        self.assertEqual(126, qresult.seq_len)
        self.assertEqual('hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s', qresult.description)
        self.assertEqual(4, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|152973462|ref|YP_001338513.1|', hit.id)
        self.assertEqual('hypothetical protein KPN_pKPN3p05904 [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(101, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(78, hsp.initn_score)
        self.assertEqual(50, hsp.init1_score)
        self.assertEqual(58, hsp.opt_score)
        self.assertEqual(100.8, hsp.z_score)
        self.assertEqual(23.9, hsp.bitscore)
        self.assertEqual(0.15, hsp.evalue)
        self.assertEqual(58, hsp.sw_score)
        self.assertAlmostEqual(28.9, hsp.ident_pct)
        self.assertAlmostEqual(63.2, hsp.pos_pct)
        self.assertEqual(38, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(38, hsp.query_end)
        self.assertEqual('MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE', str(hsp.query.seq))
        self.assertEqual(43, hsp.hit_start)
        self.assertEqual(81, hsp.hit_end)
        self.assertEqual('IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': '.:::   ..  .::..:::..      :  .:..: ..'}, hsp.aln_annotation)
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|152973509|ref|YP_001338560.1|', hit.id)
        self.assertEqual('probable sensor kinase (silver resistance) [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(448, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(73, hsp.initn_score)
        self.assertEqual(56, hsp.init1_score)
        self.assertEqual(56, hsp.opt_score)
        self.assertEqual(89.9, hsp.z_score)
        self.assertEqual(24.0, hsp.bitscore)
        self.assertEqual(0.6, hsp.evalue)
        self.assertEqual(56, hsp.sw_score)
        self.assertAlmostEqual(72.7, hsp.ident_pct)
        self.assertAlmostEqual(81.8, hsp.pos_pct)
        self.assertEqual(11, hsp.aln_span)
        self.assertEqual(50, hsp.query_start)
        self.assertEqual(61, hsp.query_end)
        self.assertEqual('FFDLIIPNGGK', str(hsp.query.seq))
        self.assertEqual(407, hsp.hit_start)
        self.assertEqual(418, hsp.hit_end)
        self.assertEqual('FFDLVIENPGK', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': '::::.: : ::'}, hsp.aln_annotation)
        # second qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|152973581|ref|YP_001338632.1|', hit.id)
        self.assertEqual('inner membrane protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(84, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(61, hsp.initn_score)
        self.assertEqual(46, hsp.init1_score)
        self.assertEqual(48, hsp.opt_score)
        self.assertEqual(88.5, hsp.z_score)
        self.assertEqual(21.3, hsp.bitscore)
        self.assertEqual(0.72, hsp.evalue)
        self.assertEqual(48, hsp.sw_score)
        self.assertAlmostEqual(30.0, hsp.ident_pct)
        self.assertAlmostEqual(67.5, hsp.pos_pct)
        self.assertEqual(40, hsp.aln_span)
        self.assertEqual(15, hsp.query_start)
        self.assertEqual(53, hsp.query_end)
        self.assertEqual('DKTLFIVYATDIYSPSE-FFSKIESDLKKKKSKGD-VFFD', str(hsp.query.seq))
        self.assertEqual(44, hsp.hit_start)
        self.assertEqual(84, hsp.hit_end)
        self.assertEqual('ESVVFILMAGFAMSVCYLFFSVLEKVINARKSKDESIYHD', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': '....::..:   .:   -::: .:. .. .::: .-.. :'}, hsp.aln_annotation)
        # second qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('gi|152973536|ref|YP_001338587.1|', hit.id)
        self.assertEqual('putative inner membrane protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(84, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(63, hsp.initn_score)
        self.assertEqual(42, hsp.init1_score)
        self.assertEqual(48, hsp.opt_score)
        self.assertEqual(88.5, hsp.z_score)
        self.assertEqual(21.3, hsp.bitscore)
        self.assertEqual(0.72, hsp.evalue)
        self.assertEqual(48, hsp.sw_score)
        self.assertAlmostEqual(26.7, hsp.ident_pct)
        self.assertAlmostEqual(66.7, hsp.pos_pct)
        self.assertEqual(30, hsp.aln_span)
        self.assertEqual(96, hsp.query_start)
        self.assertEqual(126, hsp.query_end)
        self.assertEqual('ASFFKKNFDKINVNLLSKATSFALKKGIPI', str(hsp.query.seq))
        self.assertEqual(6, hsp.hit_start)
        self.assertEqual(36, hsp.hit_end)
        self.assertEqual('ASFSKEEQDKVAVDKVAADVAWQERMNKPV', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': '::: :.. ::. :. ..  ...  . . :.'}, hsp.aln_annotation)

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('gi|10955265|ref|NP_052606.1|', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.4', qresult.version)
        self.assertEqual('NC_009649.faa', qresult.target)
        self.assertEqual(346, qresult.seq_len)
        self.assertEqual('hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s', qresult.description)
        self.assertEqual(2, len(qresult))
        # third qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|152973545|ref|YP_001338596.1|', hit.id)
        self.assertEqual('putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(242, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(72, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(70, hsp.opt_score)
        self.assertEqual(110.9, hsp.z_score)
        self.assertEqual(28.4, hsp.bitscore)
        self.assertEqual(0.041, hsp.evalue)
        self.assertEqual(70, hsp.sw_score)
        self.assertAlmostEqual(27.9, hsp.ident_pct)
        self.assertAlmostEqual(65.1, hsp.pos_pct)
        self.assertEqual(43, hsp.aln_span)
        self.assertEqual(196, hsp.query_start)
        self.assertEqual(238, hsp.query_end)
        self.assertEqual('SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEMEDYSTYC', str(hsp.query.seq))
        self.assertEqual(51, hsp.hit_start)
        self.assertEqual(94, hsp.hit_end)
        self.assertEqual('SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLISSHGEYC', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': ':...: . . :  ::.: : .:: -. . . .:. . ... ::'}, hsp.aln_annotation)
        # third qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|152973505|ref|YP_001338556.1|', hit.id)
        self.assertEqual('putative membrane fusion protein SilB [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]', hit.description)
        self.assertEqual(430, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(95, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(90.1, hsp.z_score)
        self.assertEqual(25.4, hsp.bitscore)
        self.assertEqual(0.59, hsp.evalue)
        self.assertEqual(57, hsp.sw_score)
        self.assertAlmostEqual(23.4, hsp.ident_pct)
        self.assertAlmostEqual(60.9, hsp.pos_pct)
        self.assertEqual(64, hsp.aln_span)
        self.assertEqual(39, hsp.query_start)
        self.assertEqual(101, hsp.query_end)
        self.assertEqual('ISGTYKGIDFLIKLMPSGGNTTIGRASGQNNTYFDEIALIIKENCLY--SDTKNFEYTIPKFSD', str(hsp.query.seq))
        self.assertEqual(221, hsp.hit_start)
        self.assertEqual(281, hsp.hit_end)
        self.assertEqual('IDGVITAFD-LRTGMNISKDKVVAQIQGMDPVW---ISAAVPESIAYLLKDTSQFEISVPAYPD', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': ':.:.  ..:-:   :  . . .... .:.. ..---:.  . :.  :--.::..:: ..: . :'}, hsp.aln_annotation)

    def test_output008(self):
        """Test parsing tfastx36 output (output008.m10)"""

        m10_file = get_file('output008.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(4, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('sp|Q9BS26|ERP44_HUMAN', qresult.id)
        self.assertEqual('tfastx', qresult.program)
        self.assertEqual('36.3.4', qresult.version)
        self.assertEqual('rhodopsin_nucs.fasta', qresult.target)
        self.assertEqual(406, qresult.seq_len)
        self.assertEqual('Endoplasmic reticulum resident protein 44 OS=Homo sapiens GN=ERP44', qresult.description)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('sp|Q9NSY1|BMP2K_HUMAN', qresult.id)
        self.assertEqual('tfastx', qresult.program)
        self.assertEqual('36.3.4', qresult.version)
        self.assertEqual('rhodopsin_nucs.fasta', qresult.target)
        self.assertEqual(1161, qresult.seq_len)
        self.assertEqual('BMP-2-inducible protein kinase OS=Homo sapiens GN=BMP2K', qresult.description)
        self.assertEqual(2, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|283855822|gb|GQ290312.1|', hit.id)
        self.assertEqual('Myotis ricketti voucher GQX10 rhodopsin (RHO) mRNA, partial cds', hit.description)
        self.assertEqual(983, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(106, hsp.initn_score)
        self.assertEqual(53, hsp.init1_score)
        self.assertEqual(70, hsp.opt_score)
        self.assertEqual(88.0, hsp.z_score)
        self.assertEqual(28.0, hsp.bitscore)
        self.assertEqual(0.026, hsp.evalue)
        self.assertEqual(70, hsp.sw_score)
        self.assertAlmostEqual(23.1, hsp.ident_pct)
        self.assertAlmostEqual(46.2, hsp.pos_pct)
        self.assertEqual(65, hsp.aln_span)
        self.assertEqual(452, hsp.query_start)
        self.assertEqual(514, hsp.query_end)
        self.assertEqual('LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH---HHHHHHHLLQDAYMQQYQHATQQQQML', str(hsp.query.seq))
        self.assertEqual(122, hsp.hit_start)
        self.assertEqual(317, hsp.hit_end)
        self.assertEqual('IPHQLPHALRHRPAQEAAHASQLHPAQPGCGQPLHGLWRLHHHPVYLYAWILRLRGHGMQSGGLL', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': '. :. ::  ...  :.  . .: .  :    :  :---. :::   :    ..   :. :.  .:'}, hsp.aln_annotation)
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|57163782|ref|NM_001009242.1|', hit.id)
        self.assertEqual('Felis catus rhodopsin (RHO), mRNA', hit.description)
        self.assertEqual(1047, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(105, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(68, hsp.opt_score)
        self.assertEqual(85.8, hsp.z_score)
        self.assertEqual(27.7, hsp.bitscore)
        self.assertEqual(0.034, hsp.evalue)
        self.assertEqual(72, hsp.sw_score)
        self.assertAlmostEqual(23.5, hsp.ident_pct)
        self.assertAlmostEqual(45.0, hsp.pos_pct)
        self.assertEqual(201, hsp.aln_span)
        self.assertEqual(417, hsp.query_start)
        self.assertEqual(599, hsp.query_end)
        self.assertEqual('GPEIL---LGQ-GPPQQPPQQHRVLQQLQQGDWRLQQLH-------LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH-----HHHHHH-HLLQDAYMQQYQHATQQQQMLQQQF-LMHSVYQPQPSASQYPTMMPQYQQAFFQQQMLAQHQPSQQQASPEYLTSPQEFSPALVSYTSSLPA-QVGTIMDSSYSANRS', str(hsp.query.seq))
        self.assertEqual(14, hsp.hit_start)
        self.assertEqual(595, hsp.hit_end)
        self.assertEqual('GPELLRALLQQNGCGTQPLRVPTVLPG*AMAVLHAGRLHVPAHRAWLPHQLPHALRHGPAQEAAHASQLHPAQPGRG*PLHGLRWLHHHPLH/PLCMDTLSLGPQDAIWRASLPHWAVKLPCGLWWSWPLSGTWWCVSP*ATSA------LGRTMP*WASLSPGSWHWPALHPPSLVGPGTSLKACSVHAGSTTTHSSQKS', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': ':::.:---: :-:   :: .   ::    ..  .  .::-------: :. ::  ..   :.  . .: .  :  .    :-----:::  :- : .:.     : :  . .. .   -:  ...   : .. .  . :   .:------:.. .:   . ::     :    :.::.  .:: :-.: .   ...:...:'}, hsp.aln_annotation)

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('sp|P06213|INSR_HUMAN', qresult.id)
        self.assertEqual('tfastx', qresult.program)
        self.assertEqual('36.3.4', qresult.version)
        self.assertEqual('rhodopsin_nucs.fasta', qresult.target)
        self.assertEqual(1382, qresult.seq_len)
        self.assertEqual('Insulin receptor OS=Homo sapiens GN=INSR', qresult.description)
        self.assertEqual(0, len(qresult))

        # test fourth qresult
        qresult = qresults[3]
        self.assertEqual('sp|P08100|OPSD_HUMAN', qresult.id)
        self.assertEqual('tfastx', qresult.program)
        self.assertEqual('36.3.4', qresult.version)
        self.assertEqual('rhodopsin_nucs.fasta', qresult.target)
        self.assertEqual(348, qresult.seq_len)
        self.assertEqual('Rhodopsin OS=Homo sapiens GN=RHO', qresult.description)
        self.assertEqual(6, len(qresult))
        # fourth qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|57163782|ref|NM_001009242.1|', hit.id)
        self.assertEqual('Felis catus rhodopsin (RHO), mRNA', hit.description)
        self.assertEqual(1047, hit.seq_len)
        self.assertEqual(1, len(hit))
        # fourth qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(2298, hsp.initn_score)
        self.assertEqual(2298, hsp.init1_score)
        self.assertEqual(2298, hsp.opt_score)
        self.assertEqual(3150.5, hsp.z_score)
        self.assertEqual(593.0, hsp.bitscore)
        self.assertEqual(6.7e-173, hsp.evalue)
        self.assertEqual(2298, hsp.sw_score)
        self.assertAlmostEqual(96.6, hsp.ident_pct)
        self.assertAlmostEqual(99.4, hsp.pos_pct)
        self.assertEqual(348, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(348, hsp.query_end)
        self.assertEqual('MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA', str(hsp.query.seq))
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(1044, hsp.hit_end)
        self.assertEqual('MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTLPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASTTGSKTETSQVAPA', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': '::::::::::::::: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.::::.:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.:::::::::.::::::::::::::::::::::::::::::::::.:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.:::::::..:::::::::::::::::::::.:::::::::::::.: :::::::::::'}, hsp.aln_annotation)
        # fourth qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|18148870|dbj|AB062417.1|', hit.id)
        self.assertEqual('Synthetic construct Bos taurus gene for rhodopsin, complete cds', hit.description)
        self.assertEqual(1047, hit.seq_len)
        self.assertEqual(1, len(hit))
        # fourth qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(2237, hsp.initn_score)
        self.assertEqual(2237, hsp.init1_score)
        self.assertEqual(2237, hsp.opt_score)
        self.assertEqual(3067.2, hsp.z_score)
        self.assertEqual(577.6, hsp.bitscore)
        self.assertEqual(2.9e-168, hsp.evalue)
        self.assertEqual(2237, hsp.sw_score)
        self.assertAlmostEqual(93.4, hsp.ident_pct)
        self.assertAlmostEqual(98.6, hsp.pos_pct)
        self.assertEqual(348, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(348, hsp.query_end)
        self.assertEqual('MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA', str(hsp.query.seq))
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(1044, hsp.hit_end)
        self.assertEqual('MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # fourth qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|283855822|gb|GQ290312.1|', hit.id)
        self.assertEqual('Myotis ricketti voucher GQX10 rhodopsin (RHO) mRNA, partial cds', hit.description)
        self.assertEqual(983, hit.seq_len)
        self.assertEqual(2, len(hit))
        self.assertEqual({'similarity': '::::::::::::::: ::::::::: ::::::::::::::::::::::.::::::::::::::::::::::::::::::::::::::.::::.:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.:::::::::.:::::::::: . :.:::::::::::::: ::.:.:::::::::::::::::::::::::::::::::::::::::::::::.:::.:::::::::::.::::::::::::::..:.:::::::::::::::::.::.:::::::::::::.:::::::::::::'}, hsp.aln_annotation)
        # fourth qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(2138, hsp.initn_score)
        self.assertEqual(2138, hsp.init1_score)
        self.assertEqual(2143, hsp.opt_score)
        self.assertEqual(2939.0, hsp.z_score)
        self.assertEqual(553.8, hsp.bitscore)
        self.assertEqual(4.1e-161, hsp.evalue)
        self.assertEqual(2143, hsp.sw_score)
        self.assertAlmostEqual(95.1, hsp.ident_pct)
        self.assertAlmostEqual(99.4, hsp.pos_pct)
        self.assertEqual(326, hsp.aln_span)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(336, hsp.query_end)
        self.assertEqual('VPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASAT', str(hsp.query.seq))
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(978, hsp.hit_end)
        self.assertEqual('VPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSMHGYFVFGATGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLAFTWVMALACAAPPLAGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVVAFLICWLPYASVAFYIFTHQGSNFGPVFMTIPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASTT', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual({'similarity': '::::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.::::.::::.:::::.::::::: :::::::::::::::::::::::::::::::::::::::::::::::::.:::::::::::::::::::::::::.::::::::::::::::::::::::::::::::::.::::::::::::::::::::::::::::::::::::::::.::::::.:::::::::::::::::::.:::::::::::..:::::::::::::::::::::.:::::::::::::.:'}, hsp.aln_annotation)
        # fourth qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        self.assertEqual(74, hsp.initn_score)
        self.assertEqual(58, hsp.init1_score)
        self.assertEqual(59, hsp.opt_score)
        self.assertEqual(91.3, hsp.z_score)
        self.assertEqual(26.9, hsp.bitscore)
        self.assertEqual(0.017, hsp.evalue)
        self.assertEqual(59, hsp.sw_score)
        self.assertAlmostEqual(35.5, hsp.ident_pct)
        self.assertAlmostEqual(61.3, hsp.pos_pct)
        self.assertEqual(31, hsp.aln_span)
        self.assertEqual(234, hsp.query_start)
        self.assertEqual(265, hsp.query_end)
        self.assertEqual('AQQQESATTQKAEKEVTRMVIIMVIAFLICW', str(hsp.query.seq))
        self.assertEqual(674, hsp.hit_start)
        self.assertEqual(767, hsp.hit_end)
        self.assertEqual('SQQIRNATTMMMTMRVTSFSAFWVVADSCCW', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # fourth qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('gi|2734705|gb|U59921.1|BBU59921', hit.id)
        self.assertEqual('Bufo bufo rhodopsin mRNA, complete cds', hit.description)
        self.assertEqual(1574, hit.seq_len)
        self.assertEqual(1, len(hit))
        # fourth qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(2080, hsp.initn_score)
        self.assertEqual(2031, hsp.init1_score)
        self.assertEqual(2057, hsp.opt_score)
        self.assertEqual(2819.6, hsp.z_score)
        self.assertEqual(532.4, hsp.bitscore)
        self.assertEqual(1.8e-154, hsp.evalue)
        self.assertEqual(2057, hsp.sw_score)
        self.assertAlmostEqual(83.3, hsp.ident_pct)
        self.assertAlmostEqual(95.2, hsp.pos_pct)
        self.assertEqual(354, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(348, hsp.query_end)
        self.assertEqual('MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEAS-ATVSKTE-----TSQVAPA', str(hsp.query.seq))
        self.assertEqual(41, hsp.hit_start)
        self.assertEqual(1103, hsp.hit_end)
        self.assertEqual('MNGTEGPNFYIPMSNKTGVVRSPFEYPQYYLAEPWQYSILCAYMFLLILLGFPINFMTLYVTIQHKKLRTPLNYILLNLAFANHFMVLCGFTVTMYSSMNGYFILGATGCYVEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFSENHAVMGVAFTWIMALSCAVPPLLGWSRYIPEGMQCSCGVDYYTLKPEVNNESFVIYMFVVHFTIPLIIIFFCYGRLVCTVKEAAAQQQESATTQKAEKEVTRMVIIMVVFFLICWVPYASVAFFIFSNQGSEFGPIFMTVPAFFAKSSSIYNPVIYIMLNKQFRNCMITTLCCGKNPFGEDDASSAATSKTEASSVSSSQVSPA', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # fourth qresult, fifth hit
        hit = qresult[4]
        self.assertEqual('gi|12583664|dbj|AB043817.1|', hit.id)
        self.assertEqual('Conger myriaster conf gene for fresh water form rod opsin, complete cds', hit.description)
        self.assertEqual(1344, hit.seq_len)
        self.assertEqual(1, len(hit))
        # fourth qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual(1975, hsp.initn_score)
        self.assertEqual(1951, hsp.init1_score)
        self.assertEqual(1993, hsp.opt_score)
        self.assertEqual(2732.8, hsp.z_score)
        self.assertEqual(516.1, hsp.bitscore)
        self.assertEqual(1.2e-149, hsp.evalue)
        self.assertEqual(1993, hsp.sw_score)
        self.assertAlmostEqual(81.6, hsp.ident_pct)
        self.assertAlmostEqual(93.9, hsp.pos_pct)
        self.assertEqual(347, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(346, hsp.query_end)
        self.assertEqual('MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGD-DEASATVSKTETSQVA', str(hsp.query.seq))
        self.assertEqual(22, hsp.hit_start)
        self.assertEqual(1063, hsp.hit_end)
        self.assertEqual('MNGTEGPNFYIPMSNATGVVRSPFEYPQYYLAEPWAFSALSAYMFFLIIAGFPINFLTLYVTIEHKKLRTPLNYILLNLAVADLFMVFGGFTTTMYTSMHGYFVFGPTGCNIEGFFATLGGEIALWCLVVLAIERWMVVCKPVTNFRFGESHAIMGVMVTWTMALACALPPLFGWSRYIPEGLQCSCGIDYYTRAPGINNESFVIYMFTCHFSIPLAVISFCYGRLVCTVKEAAAQQQESETTQRAEREVTRMVVIMVISFLVCWVPYASVAWYIFTHQGSTFGPIFMTIPSFFAKSSALYNPMIYICMNKQFRHCMITTLCCGKNPFEEEDGASATSSKTEASSVS', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # fourth qresult, sixth hit
        hit = qresult[5]
        self.assertEqual('gi|283855845|gb|GQ290303.1|', hit.id)
        self.assertEqual('Cynopterus brachyotis voucher 20020434 rhodopsin (RHO) gene, exons 1 through 5 and partial cds', hit.description)
        self.assertEqual(4301, hit.seq_len)
        self.assertEqual(4, len(hit))
        # fourth qresult, sixth hit, first hsp
        hsp = qresult[5].hsps[0]
        self.assertEqual(2094, hsp.initn_score)
        self.assertEqual(723, hsp.init1_score)
        self.assertEqual(723, hsp.opt_score)
        self.assertEqual(992.9, hsp.z_score)
        self.assertEqual(195.8, hsp.bitscore)
        self.assertEqual(1e-52, hsp.evalue)
        self.assertEqual(723, hsp.sw_score)
        self.assertAlmostEqual(96.4, hsp.ident_pct)
        self.assertAlmostEqual(99.1, hsp.pos_pct)
        self.assertEqual(111, hsp.aln_span)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(121, hsp.query_end)
        self.assertEqual('VPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGG', str(hsp.query.seq))
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(333, hsp.hit_end)
        self.assertEqual('VPFSNKTGVVRSPFEHPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGG', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # fourth qresult, sixth hit, second hsp
        hsp = qresult[5].hsps[1]
        self.assertEqual(1411, hsp.initn_score)
        self.assertEqual(499, hsp.init1_score)
        self.assertEqual(501, hsp.opt_score)
        self.assertEqual(992.9, hsp.z_score)
        self.assertEqual(195.8, hsp.bitscore)
        self.assertEqual(1e-52, hsp.evalue)
        self.assertEqual(783, hsp.sw_score)
        self.assertAlmostEqual(75.4, hsp.ident_pct)
        self.assertAlmostEqual(79.5, hsp.pos_pct)
        self.assertEqual(172, hsp.aln_span)
        self.assertEqual(176, hsp.query_start)
        self.assertEqual(312, hsp.query_end)
        self.assertEqual('RYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKE------------------------------------AAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQ', str(hsp.query.seq))
        self.assertEqual(2854, hsp.hit_start)
        self.assertEqual(3368, hsp.hit_end)
        self.assertEqual(r'RYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEVRSCVGHWGHAH*VNGAQLHSQSCHSLDT*PCVPA\AAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSNFGPIFMTLPAFFAKSSSIYNPVIYIMMNKQ', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # fourth qresult, sixth hit, third hsp
        hsp = qresult[5].hsps[2]
        self.assertEqual(431, hsp.initn_score)
        self.assertEqual(379, hsp.init1_score)
        self.assertEqual(388, hsp.opt_score)
        self.assertEqual(992.9, hsp.z_score)
        self.assertEqual(195.8, hsp.bitscore)
        self.assertEqual(1e-52, hsp.evalue)
        self.assertEqual(388, hsp.sw_score)
        self.assertAlmostEqual(80.8, hsp.ident_pct)
        self.assertAlmostEqual(90.4, hsp.pos_pct)
        self.assertEqual(73, hsp.aln_span)
        self.assertEqual(118, hsp.query_start)
        self.assertEqual(189, hsp.query_end)
        self.assertEqual('LGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSR--YIPEGLQCSCGI', str(hsp.query.seq))
        self.assertEqual(1403, hsp.hit_start)
        self.assertEqual(1619, hsp.hit_end)
        self.assertEqual('LAGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWVMALACAAPPLVGWSR*WH*TEG-KCL*GL', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # fourth qresult, sixth hit, fourth hsp
        hsp = qresult[5].hsps[3]
        self.assertEqual(213, hsp.initn_score)
        self.assertEqual(171, hsp.init1_score)
        self.assertEqual(176, hsp.opt_score)
        self.assertEqual(992.9, hsp.z_score)
        self.assertEqual(195.8, hsp.bitscore)
        self.assertEqual(1e-52, hsp.evalue)
        self.assertEqual(176, hsp.sw_score)
        self.assertAlmostEqual(76.7, hsp.ident_pct)
        self.assertAlmostEqual(93.3, hsp.pos_pct)
        self.assertEqual(30, hsp.aln_span)
        self.assertEqual(306, hsp.query_start)
        self.assertEqual(336, hsp.query_end)
        self.assertEqual('IMMNKQFRNCMLTTICCGKNPLGDDEASAT', str(hsp.query.seq))
        self.assertEqual(4206, hsp.hit_start)
        self.assertEqual(4296, hsp.hit_end)
        self.assertEqual('MLLAFQFRNCMLTTLCCGKNPLGDDEASTT', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)

    def test_output009(self):
        """Test parsing fasta36 output (output009.m10)"""

        m10_file = get_file('output009.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(3, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('mrnalib.fasta', qresult.target)
        self.assertEqual(15, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|255708421:1-99', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('mrnalib.fasta', qresult.target)
        self.assertEqual(99, qresult.seq_len)
        self.assertEqual('Mus musculus myoglobin (Mb), transcript varia', qresult.description)
        self.assertEqual(1, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|23308614|ref|NM_152952.1|', hit.id)
        self.assertEqual('Danio rerio cytoglobin 1 (cygb1), mRNA', hit.description)
        self.assertEqual(5188, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(124, hsp.initn_score)
        self.assertEqual(74, hsp.init1_score)
        self.assertEqual(74, hsp.opt_score)
        self.assertEqual(70.3, hsp.z_score)
        self.assertEqual(23.6, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(81.8, hsp.ident_pct)
        self.assertAlmostEqual(81.8, hsp.pos_pct)
        self.assertEqual(22, hsp.aln_span)
        self.assertEqual(7, hsp.query_start)
        self.assertEqual(29, hsp.query_end)
        self.assertEqual('TGATGTTCTGTTTCTAAAACAG', str(hsp.query.seq))
        self.assertEqual(3483, hsp.hit_start)
        self.assertEqual(3505, hsp.hit_end)
        self.assertEqual('TGATTTTTTTTGTCTAAAACAG', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('gi|156718121:2361-2376', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('mrnalib.fasta', qresult.target)
        self.assertEqual(16, qresult.seq_len)
        self.assertEqual('Bos taurus nucleoporin 43kDa (NU', qresult.description)
        self.assertEqual(5, len(qresult))
        # third qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|47271416|ref|NM_131257.2|', hit.id)
        self.assertEqual('Danio rerio hemoglobin alpha adult-1 (hbaa1), mRNA', hit.description)
        self.assertEqual(597, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(52, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(52, hsp.opt_score)
        self.assertEqual(87.3, hsp.z_score)
        self.assertEqual(21.0, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(85.7, hsp.ident_pct)
        self.assertAlmostEqual(85.7, hsp.pos_pct)
        self.assertEqual(14, hsp.aln_span)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('AGAAGGAAAAAAAA', str(hsp.query.seq))
        self.assertEqual(572, hsp.hit_start)
        self.assertEqual(586, hsp.hit_end)
        self.assertEqual('AGAACTAAAAAAAA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # third qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|332859474|ref|XM_001156938.2|', hit.id)
        self.assertEqual('PREDICTED: Pan troglodytes myoglobin, transcript variant 11 (MB), mRNA', hit.description)
        self.assertEqual(762, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(52, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(52, hsp.opt_score)
        self.assertEqual(85.3, hsp.z_score)
        self.assertEqual(20.9, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(85.7, hsp.ident_pct)
        self.assertAlmostEqual(85.7, hsp.pos_pct)
        self.assertEqual(14, hsp.aln_span)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('AGAAGGAAAAAAAA', str(hsp.query.seq))
        self.assertEqual(84, hsp.hit_start)
        self.assertEqual(98, hsp.hit_end)
        self.assertEqual('AGAAGGTATAAAAA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # third qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|332211534|ref|XM_003254825.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys hemoglobin subunit gamma-2-like (LOC100581638), mRNA', hit.description)
        self.assertEqual(805, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(52, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(52, hsp.opt_score)
        self.assertEqual(84.9, hsp.z_score)
        self.assertEqual(20.9, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(85.7, hsp.ident_pct)
        self.assertAlmostEqual(85.7, hsp.pos_pct)
        self.assertEqual(14, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(15, hsp.query_end)
        self.assertEqual('TTTTTTTCCTTCTT', str(hsp.query.seq))
        self.assertEqual(633, hsp.hit_start)
        self.assertEqual(647, hsp.hit_end)
        self.assertEqual('TTTTTTTACATCTT', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)
        # third qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('gi|23308614|ref|NM_152952.1|', hit.id)
        self.assertEqual('Danio rerio cytoglobin 1 (cygb1), mRNA', hit.description)
        self.assertEqual(5188, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(52, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(52, hsp.opt_score)
        self.assertEqual(70.2, hsp.z_score)
        self.assertEqual(20.9, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(85.7, hsp.ident_pct)
        self.assertAlmostEqual(85.7, hsp.pos_pct)
        self.assertEqual(14, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(15, hsp.query_end)
        self.assertEqual('AAGAAGGAAAAAAA', str(hsp.query.seq))
        self.assertEqual(3547, hsp.hit_start)
        self.assertEqual(3561, hsp.hit_end)
        self.assertEqual('AATAAGTAAAAAAA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # third qresult, fifth hit
        hit = qresult[4]
        self.assertEqual('gi|297689475|ref|XM_002822130.1|', hit.id)
        self.assertEqual('PREDICTED: Pongo abelii hemoglobin subunit gamma-like (LOC100439631), mRNA', hit.description)
        self.assertEqual(1158, hit.seq_len)
        self.assertEqual(2, len(hit))
        # third qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual(52, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(52, hsp.opt_score)
        self.assertEqual(82.0, hsp.z_score)
        self.assertEqual(20.9, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(85.7, hsp.ident_pct)
        self.assertAlmostEqual(85.7, hsp.pos_pct)
        self.assertEqual(14, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(15, hsp.query_end)
        self.assertEqual('TTTTTTTCCTTCTT', str(hsp.query.seq))
        self.assertEqual(983, hsp.hit_start)
        self.assertEqual(997, hsp.hit_end)
        self.assertEqual('TTTTTTTACATCTT', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)
        # third qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        self.assertEqual(52, hsp.initn_score)
        self.assertEqual(52, hsp.init1_score)
        self.assertEqual(52, hsp.opt_score)
        self.assertEqual(82.0, hsp.z_score)
        self.assertEqual(20.9, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(85.7, hsp.ident_pct)
        self.assertAlmostEqual(85.7, hsp.pos_pct)
        self.assertEqual(14, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(15, hsp.query_end)
        self.assertEqual('AAGAAGGAAAAAAA', str(hsp.query.seq))
        self.assertEqual(19, hsp.hit_start)
        self.assertEqual(33, hsp.hit_end)
        self.assertEqual('AAGAAGGTAAAAGA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)

    def test_output010(self):
        """Test parsing fasta36 output (output010.m10)"""

        m10_file = get_file('output010.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(1, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('mrnalib.fasta', qresult.target)
        self.assertEqual(15, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(0, len(qresult))

    def test_output011(self):
        """Test parsing fasta36 output (output011.m10)"""

        m10_file = get_file('output011.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(1, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|255708421:1-99', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('mrnalib.fasta', qresult.target)
        self.assertEqual(99, qresult.seq_len)
        self.assertEqual('Mus musculus myoglobin (Mb), transcript varia', qresult.description)
        self.assertEqual(5, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|284005422|ref|NM_001171502.1|', hit.id)
        self.assertEqual('Oryctolagus cuniculus hemoglobin, zeta (HBZ_2), mRNA', hit.description)
        self.assertEqual(429, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(79, hsp.initn_score)
        self.assertEqual(66, hsp.init1_score)
        self.assertEqual(99, hsp.opt_score)
        self.assertEqual(90.0, hsp.z_score)
        self.assertEqual(23.6, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(73.8, hsp.ident_pct)
        self.assertAlmostEqual(73.8, hsp.pos_pct)
        self.assertEqual(42, hsp.aln_span)
        self.assertEqual(33, hsp.query_start)
        self.assertEqual(75, hsp.query_end)
        self.assertEqual('CAACATCCAGAGGACTGTCATCCTTGTCCCTGTGGGTGAGGG', str(hsp.query.seq))
        self.assertEqual(11, hsp.hit_start)
        self.assertEqual(52, hsp.hit_end)
        self.assertEqual('CAAGAGCGAGAGGACCATCAT-CATGTCCCTCTGGGACAAGG', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|284005386|ref|NM_001171415.1|', hit.id)
        self.assertEqual('Oryctolagus cuniculus zeta globin (HBZ0), mRNA', hit.description)
        self.assertEqual(429, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(73, hsp.initn_score)
        self.assertEqual(66, hsp.init1_score)
        self.assertEqual(99, hsp.opt_score)
        self.assertEqual(90.0, hsp.z_score)
        self.assertEqual(23.6, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(73.8, hsp.ident_pct)
        self.assertAlmostEqual(73.8, hsp.pos_pct)
        self.assertEqual(42, hsp.aln_span)
        self.assertEqual(33, hsp.query_start)
        self.assertEqual(75, hsp.query_end)
        self.assertEqual('CAACATCCAGAGGACTGTCATCCTTGTCCCTGTGGGTGAGGG', str(hsp.query.seq))
        self.assertEqual(11, hsp.hit_start)
        self.assertEqual(52, hsp.hit_end)
        self.assertEqual('CAAGAGCGAGAGGACCATCAT-CATGTCCCTCTGGGACAAGG', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|284005381|ref|NM_001171414.1|', hit.id)
        self.assertEqual('Oryctolagus cuniculus hemoglobin subunit zeta (RA_M008_JSM295ECF), mRNA', hit.description)
        self.assertEqual(429, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(90, hsp.initn_score)
        self.assertEqual(66, hsp.init1_score)
        self.assertEqual(99, hsp.opt_score)
        self.assertEqual(90.0, hsp.z_score)
        self.assertEqual(23.6, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(73.8, hsp.ident_pct)
        self.assertAlmostEqual(73.8, hsp.pos_pct)
        self.assertEqual(42, hsp.aln_span)
        self.assertEqual(33, hsp.query_start)
        self.assertEqual(75, hsp.query_end)
        self.assertEqual('CAACATCCAGAGGACTGTCATCCTTGTCCCTGTGGGTGAGGG', str(hsp.query.seq))
        self.assertEqual(11, hsp.hit_start)
        self.assertEqual(52, hsp.hit_end)
        self.assertEqual('CAAGAGCGAGAGGACCATCAT-CATGTCCCTCTGGGACAAGG', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('gi|23308614|ref|NM_152952.1|', hit.id)
        self.assertEqual('Danio rerio cytoglobin 1 (cygb1), mRNA', hit.description)
        self.assertEqual(5188, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(115, hsp.initn_score)
        self.assertEqual(78, hsp.init1_score)
        self.assertEqual(80, hsp.opt_score)
        self.assertEqual(70.4, hsp.z_score)
        self.assertEqual(23.6, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(80.0, hsp.ident_pct)
        self.assertAlmostEqual(80.0, hsp.pos_pct)
        self.assertEqual(25, hsp.aln_span)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(37, hsp.query_end)
        self.assertEqual('TTAGAAACAGAACATCATCTTCAAC', str(hsp.query.seq))
        self.assertEqual(415, hsp.hit_start)
        self.assertEqual(440, hsp.hit_end)
        self.assertEqual('TGACAAACTCAACACCATCTTCAAC', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, fifth hit
        hit = qresult[4]
        self.assertEqual('gi|291415427|ref|XM_002723908.1|', hit.id)
        self.assertEqual('PREDICTED: Oryctolagus cuniculus zeta globin-like (LOC100357627), mRNA', hit.description)
        self.assertEqual(423, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual(91, hsp.initn_score)
        self.assertEqual(66, hsp.init1_score)
        self.assertEqual(99, hsp.opt_score)
        self.assertEqual(90.0, hsp.z_score)
        self.assertEqual(23.6, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(73.8, hsp.ident_pct)
        self.assertAlmostEqual(73.8, hsp.pos_pct)
        self.assertEqual(42, hsp.aln_span)
        self.assertEqual(33, hsp.query_start)
        self.assertEqual(75, hsp.query_end)
        self.assertEqual('CAACATCCAGAGGACTGTCATCCTTGTCCCTGTGGGTGAGGG', str(hsp.query.seq))
        self.assertEqual(11, hsp.hit_start)
        self.assertEqual(52, hsp.hit_end)
        self.assertEqual('CAAGAGCGAGAGGACCATCAT-CATGTCCCTCTGGGACAAGG', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)

    def test_output012(self):
        """Test parsing fasta36 output (output012.m10)"""

        m10_file = get_file('output012.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(1, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|156718121:2361-2376', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('mrnalib.fasta', qresult.target)
        self.assertEqual(16, qresult.seq_len)
        self.assertEqual('Bos taurus nucleoporin 43kDa (NU', qresult.description)
        self.assertEqual(8, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|77681427|ref|NM_001016495.2|', hit.id)
        self.assertEqual('Xenopus (Silurana) tropicalis hemoglobin, epsilon 1 (hbe1), mRNA', hit.description)
        self.assertEqual(554, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(88.1, hsp.z_score)
        self.assertEqual(21.0, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('AAGAAGGAAAAAAAA', str(hsp.query.seq))
        self.assertEqual(529, hsp.hit_start)
        self.assertEqual(544, hsp.hit_end)
        self.assertEqual('AAGAATAAAAAAAAA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|53749657|ref|NM_182940.2|', hit.id)
        self.assertEqual('Danio rerio hemoglobin alpha embryonic-1 (hbae1), mRNA', hit.description)
        self.assertEqual(564, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(87.9, hsp.z_score)
        self.assertEqual(21.0, hsp.bitscore)
        self.assertEqual(1.4, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('AAGAAGGAAAAAAAA', str(hsp.query.seq))
        self.assertEqual(531, hsp.hit_start)
        self.assertEqual(546, hsp.hit_end)
        self.assertEqual('AAGAAAAAAAAAAAA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|158508631|ref|NM_030206.4|', hit.id)
        self.assertEqual('Mus musculus cytoglobin (Cygb), mRNA', hit.description)
        self.assertEqual(2331, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(62, hsp.initn_score)
        self.assertEqual(62, hsp.init1_score)
        self.assertEqual(62, hsp.opt_score)
        self.assertEqual(76.3, hsp.z_score)
        self.assertEqual(20.9, hsp.bitscore)
        self.assertEqual(1.5, hsp.evalue)
        self.assertAlmostEqual(87.5, hsp.ident_pct)
        self.assertAlmostEqual(87.5, hsp.pos_pct)
        self.assertEqual(16, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('GAAGAAGGAAAAAAAA', str(hsp.query.seq))
        self.assertEqual(1781, hsp.hit_start)
        self.assertEqual(1797, hsp.hit_end)
        self.assertEqual('GAAGAAAAAAAAAAAA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('gi|296217287|ref|XM_002754912.1|', hit.id)
        self.assertEqual('PREDICTED: Callithrix jacchus hemoglobin subunit gamma-like (LOC100389093), mRNA', hit.description)
        self.assertEqual(627, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(86.3, hsp.z_score)
        self.assertEqual(20.8, hsp.bitscore)
        self.assertEqual(1.5, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(15, hsp.query_end)
        self.assertEqual('GAAGAAGGAAAAAAA', str(hsp.query.seq))
        self.assertEqual(131, hsp.hit_start)
        self.assertEqual(146, hsp.hit_end)
        self.assertEqual('GAAGAAGGTAAAACA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, fifth hit
        hit = qresult[4]
        self.assertEqual('gi|332211540|ref|XM_003254828.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys hemoglobin subunit gamma-1-like, transcript variant 2 (LOC100582529), mRNA', hit.description)
        self.assertEqual(640, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(86.0, hsp.z_score)
        self.assertEqual(20.8, hsp.bitscore)
        self.assertEqual(1.5, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('TTTTTTTTCCTTCTT', str(hsp.query.seq))
        self.assertEqual(589, hsp.hit_start)
        self.assertEqual(604, hsp.hit_end)
        self.assertEqual('TTTTTTTTACATCTT', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)
        # first qresult, sixth hit
        hit = qresult[5]
        self.assertEqual('gi|147903656|ref|NM_001086277.1|', hit.id)
        self.assertEqual('Xenopus laevis hemoglobin, alpha 2 (hba2), mRNA', hit.description)
        self.assertEqual(677, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, sixth hit, first hsp
        hsp = qresult[5].hsps[0]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(85.2, hsp.z_score)
        self.assertEqual(20.8, hsp.bitscore)
        self.assertEqual(1.6, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('AAGAAGGAAAAAAAA', str(hsp.query.seq))
        self.assertEqual(642, hsp.hit_start)
        self.assertEqual(657, hsp.hit_end)
        self.assertEqual('AAAAAAGAAAAAAAA', str(hsp.hit.seq))
        self.assertEqual(1, hsp.query_strand)
        # first qresult, seventh hit
        hit = qresult[6]
        self.assertEqual('gi|380013536|ref|XM_003690762.1|', hit.id)
        self.assertEqual('PREDICTED: Apis florea globin-like (LOC100870092), mRNA', hit.description)
        self.assertEqual(707, hit.seq_len)
        self.assertEqual(3, len(hit))
        # first qresult, seventh hit, first hsp
        hsp = qresult[6].hsps[0]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(84.6, hsp.z_score)
        self.assertEqual(20.7, hsp.bitscore)
        self.assertEqual(1.7, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('TTTTTTTTCCTTCTT', str(hsp.query.seq))
        self.assertEqual(6, hsp.hit_start)
        self.assertEqual(21, hsp.hit_end)
        self.assertEqual('TTTTTTTTTTTTCTT', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)
        # first qresult, seventh hit, second hsp
        hsp = qresult[6].hsps[1]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(84.6, hsp.z_score)
        self.assertEqual(20.7, hsp.bitscore)
        self.assertEqual(1.7, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('TTTTTTTTCCTTCTT', str(hsp.query.seq))
        self.assertEqual(45, hsp.hit_start)
        self.assertEqual(60, hsp.hit_end)
        self.assertEqual('TTTTTTTTTTTTCTT', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)
        # first qresult, seventh hit, third hsp
        hsp = qresult[6].hsps[2]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(84.6, hsp.z_score)
        self.assertEqual(20.7, hsp.bitscore)
        self.assertEqual(1.7, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(15, hsp.query_end)
        self.assertEqual('TTTTTTTCCTTCTTC', str(hsp.query.seq))
        self.assertEqual(644, hsp.hit_start)
        self.assertEqual(659, hsp.hit_end)
        self.assertEqual('TTTTTTTCCTCCTCC', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)
        # first qresult, eighth hit
        hit = qresult[7]
        self.assertEqual('gi|332211538|ref|XM_003254827.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys hemoglobin subunit gamma-1-like, transcript variant 1 (LOC100582529), mRNA', hit.description)
        self.assertEqual(713, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, eighth hit, first hsp
        hsp = qresult[7].hsps[0]
        self.assertEqual(57, hsp.initn_score)
        self.assertEqual(57, hsp.init1_score)
        self.assertEqual(57, hsp.opt_score)
        self.assertEqual(84.5, hsp.z_score)
        self.assertEqual(20.7, hsp.bitscore)
        self.assertEqual(1.7, hsp.evalue)
        self.assertAlmostEqual(86.7, hsp.ident_pct)
        self.assertAlmostEqual(86.7, hsp.pos_pct)
        self.assertEqual(15, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(16, hsp.query_end)
        self.assertEqual('TTTTTTTTCCTTCTT', str(hsp.query.seq))
        self.assertEqual(662, hsp.hit_start)
        self.assertEqual(677, hsp.hit_end)
        self.assertEqual('TTTTTTTTACATCTT', str(hsp.hit.seq))
        self.assertEqual(-1, hsp.query_strand)

    def test_output013(self):
        """Test parsing fasta36 output (output013.m10)"""

        m10_file = get_file('output013.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(3, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('protlib.fasta', qresult.target)
        self.assertEqual(16, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('sp|Q9Y2H6|68-133', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('protlib.fasta', qresult.target)
        self.assertEqual(66, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(1, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|291391832|ref|XP_002712264.1|', hit.id)
        self.assertEqual('PREDICTED: titin [Oryctolagus cuniculus]', hit.description)
        self.assertEqual(33406, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(98, hsp.initn_score)
        self.assertEqual(98, hsp.init1_score)
        self.assertEqual(109, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(30.2, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(109, hsp.sw_score)
        self.assertAlmostEqual(26.8, hsp.ident_pct)
        self.assertAlmostEqual(54.9, hsp.pos_pct)
        self.assertEqual(71, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(66, hsp.query_end)
        self.assertEqual('PNGSVPPIY-----VPPGYAPQVIEDNGVRRVVVVPQAPEFHPGSHTVLHRSPHPPLPGFIPVPTMMPPPP', str(hsp.query.seq))
        self.assertEqual(10704, hsp.hit_start)
        self.assertEqual(10775, hsp.hit_end)
        self.assertEqual('PEKKVPPAVPKKPEAPPAKVPEAPKEVVPEKKIAVPKKPEVPPAKVPEVPKKPVIEEKPVIPVPKKVESPP', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('sp|Q9Y2H6|265-345', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('protlib.fasta', qresult.target)
        self.assertEqual(81, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(4, len(qresult))
        # third qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|260806189|ref|XP_002597967.1|', hit.id)
        self.assertEqual('hypothetical protein BRAFLDRAFT_79792 [Branchiostoma floridae]', hit.description)
        self.assertEqual(23830, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(220, hsp.initn_score)
        self.assertEqual(62, hsp.init1_score)
        self.assertEqual(92, hsp.opt_score)
        self.assertEqual(97.4, hsp.z_score)
        self.assertEqual(30.5, hsp.bitscore)
        self.assertEqual(0.32, hsp.evalue)
        self.assertEqual(92, hsp.sw_score)
        self.assertAlmostEqual(31.6, hsp.ident_pct)
        self.assertAlmostEqual(60.8, hsp.pos_pct)
        self.assertEqual(79, hsp.aln_span)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('LSNIVKPVASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVG-EETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(22430, hsp.hit_start)
        self.assertEqual(22499, hsp.hit_end)
        self.assertEqual('VSNI-RPAASDISPHTLTLTWDTP------EDDGGSLITSYVVEMFDVS---DGKWQTLTTTCRRPPYPVKGLNPSATY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # third qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|348553521|ref|XP_003462575.1|', hit.id)
        self.assertEqual('PREDICTED: receptor-type tyrosine-protein phosphatase F isoform 1 [Cavia porcellus]', hit.description)
        self.assertEqual(1899, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(104, hsp.initn_score)
        self.assertEqual(75, hsp.init1_score)
        self.assertEqual(75, hsp.opt_score)
        self.assertEqual(96.6, hsp.z_score)
        self.assertEqual(26.7, hsp.bitscore)
        self.assertEqual(0.36, hsp.evalue)
        self.assertEqual(75, hsp.sw_score)
        self.assertAlmostEqual(32.4, hsp.ident_pct)
        self.assertAlmostEqual(64.9, hsp.pos_pct)
        self.assertEqual(37, hsp.aln_span)
        self.assertEqual(43, hsp.query_start)
        self.assertEqual(80, hsp.query_end)
        self.assertEqual('YEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH', str(hsp.query.seq))
        self.assertEqual(542, hsp.hit_start)
        self.assertEqual(579, hsp.hit_end)
        self.assertEqual('YELVYWAAEEEGQQRKVTFDPTSSYTLEDLKPDTLYH', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # third qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|348553523|ref|XP_003462576.1|', hit.id)
        self.assertEqual('PREDICTED: receptor-type tyrosine-protein phosphatase F isoform 2 [Cavia porcellus]', hit.description)
        self.assertEqual(1908, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(104, hsp.initn_score)
        self.assertEqual(75, hsp.init1_score)
        self.assertEqual(75, hsp.opt_score)
        self.assertEqual(96.6, hsp.z_score)
        self.assertEqual(26.7, hsp.bitscore)
        self.assertEqual(0.36, hsp.evalue)
        self.assertEqual(75, hsp.sw_score)
        self.assertAlmostEqual(32.4, hsp.ident_pct)
        self.assertAlmostEqual(64.9, hsp.pos_pct)
        self.assertEqual(37, hsp.aln_span)
        self.assertEqual(43, hsp.query_start)
        self.assertEqual(80, hsp.query_end)
        self.assertEqual('YEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH', str(hsp.query.seq))
        self.assertEqual(542, hsp.hit_start)
        self.assertEqual(579, hsp.hit_end)
        self.assertEqual('YELVYWAAEEEGQQRKVTFDPTSSYTLEDLKPDTLYH', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # third qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('gi|221124183|ref|XP_002154464.1|', hit.id)
        self.assertEqual('PREDICTED: similar to FAD104 [Hydra magnipapillata]', hit.description)
        self.assertEqual(860, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(85, hsp.initn_score)
        self.assertEqual(66, hsp.init1_score)
        self.assertEqual(70, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(25.3, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(70, hsp.sw_score)
        self.assertAlmostEqual(27.1, hsp.ident_pct)
        self.assertAlmostEqual(58.6, hsp.pos_pct)
        self.assertEqual(70, hsp.aln_span)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(615, hsp.hit_start)
        self.assertEqual(673, hsp.hit_end)
        self.assertEqual('ASSISYHSIKLKWGHQSS-------KKSI-----LNHTLQMQNKSGSFNTVYSGMDTSFTLSKLKELTPY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)

    def test_output014(self):
        """Test parsing fasta36 output (output014.m10)"""

        m10_file = get_file('output014.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(1, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('protlib.fasta', qresult.target)
        self.assertEqual(16, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(0, len(qresult))

    def test_output015(self):
        """Test parsing fasta36 output (output015.m10)"""

        m10_file = get_file('output015.m10')
        qresults = list(parse(m10_file, FMT))
        self.assertEqual(1, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('sp|Q9Y2H6|68-133', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('protlib.fasta', qresult.target)
        self.assertEqual(66, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(2, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|194762369|ref|XP_001963317.1|', hit.id)
        self.assertEqual('GF14002 [Drosophila ananassae]', hit.description)
        self.assertEqual(1761, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(88, hsp.initn_score)
        self.assertEqual(68, hsp.init1_score)
        self.assertEqual(85, hsp.opt_score)
        self.assertEqual(95.3, hsp.z_score)
        self.assertEqual(26.0, hsp.bitscore)
        self.assertEqual(0.42, hsp.evalue)
        self.assertEqual(85, hsp.sw_score)
        self.assertAlmostEqual(31.0, hsp.ident_pct)
        self.assertAlmostEqual(49.3, hsp.pos_pct)
        self.assertEqual(71, hsp.aln_span)
        self.assertEqual(5, hsp.query_start)
        self.assertEqual(66, hsp.query_end)
        self.assertEqual('PPIY----VPPGYA---PQVIEDNGVRRVVVVPQAPEFH---PGSHTVLHRSPHPPLPGFIPVPTMMPPPP', str(hsp.query.seq))
        self.assertEqual(128, hsp.hit_start)
        self.assertEqual(195, hsp.hit_end)
        self.assertEqual('PPLLQQTATPPQGAQIVPPVCALHHPQQQLALMAAMQHHHPLPPPHA-LHHAPLPPPP---PLPLNPGPPP', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|77812697|ref|NP_035782.3|', hit.id)
        self.assertEqual('titin isoform N2-A [Mus musculus]', hit.description)
        self.assertEqual(33467, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(104, hsp.initn_score)
        self.assertEqual(92, hsp.init1_score)
        self.assertEqual(106, hsp.opt_score)
        self.assertEqual(94.9, hsp.z_score)
        self.assertEqual(30.2, hsp.bitscore)
        self.assertEqual(0.45, hsp.evalue)
        self.assertEqual(106, hsp.sw_score)
        self.assertAlmostEqual(29.4, hsp.ident_pct)
        self.assertAlmostEqual(57.4, hsp.pos_pct)
        self.assertEqual(68, hsp.aln_span)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(66, hsp.query_end)
        self.assertEqual('PNGSVPPIY--VPPGYAPQVIEDNGVRRVVVVPQAPEFHPGSHTVLHRSPHPPLPGFIPVPTMMPPPP', str(hsp.query.seq))
        self.assertEqual(10780, hsp.hit_start)
        self.assertEqual(10848, hsp.hit_end)
        self.assertEqual('PEKKVPPKKPEAPPAKVPEVPKEVVTEKKVAVPKKPEVPPAKVPEVPKKPVIEEKPAIPVVEKVASPP', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)

    def test_output016(self):
        """Test parsing fasta36 output (output016.m10)"""

        m10_file = get_file('output016.m10')
        # only check 8, 10 is too many
        qresults = list(parse(m10_file, FMT))[:8]
        self.assertEqual(1, len(qresults))
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('sp|Q9Y2H6|265-345', qresult.id)
        self.assertEqual('fasta', qresult.program)
        self.assertEqual('36.3.5c', qresult.version)
        self.assertEqual('protlib.fasta', qresult.target)
        self.assertEqual(81, qresult.seq_len)
        self.assertEqual('', qresult.description)
        self.assertEqual(17, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|167518632|ref|XP_001743656.1|', hit.id)
        self.assertEqual('hypothetical protein [Monosiga brevicollis MX1]', hit.description)
        self.assertEqual(1145, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(88, hsp.initn_score)
        self.assertEqual(68, hsp.init1_score)
        self.assertEqual(68, hsp.opt_score)
        self.assertEqual(97.5, hsp.z_score)
        self.assertEqual(26.1, hsp.bitscore)
        self.assertEqual(0.32, hsp.evalue)
        self.assertEqual(68, hsp.sw_score)
        self.assertAlmostEqual(43.5, hsp.ident_pct)
        self.assertAlmostEqual(65.2, hsp.pos_pct)
        self.assertEqual(23, hsp.aln_span)
        self.assertEqual(56, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('YKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(424, hsp.hit_start)
        self.assertEqual(447, hsp.hit_end)
        self.assertEqual('FRPVYTGIDTNYKVVDLTPNCDY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|9507013|ref|NP_062122.1|', hit.id)
        self.assertEqual('receptor-type tyrosine-protein phosphatase F precursor [Rattus norvegicus]', hit.description)
        self.assertEqual(1898, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(83, hsp.initn_score)
        self.assertEqual(43, hsp.init1_score)
        self.assertEqual(72, hsp.opt_score)
        self.assertEqual(96.2, hsp.z_score)
        self.assertEqual(26.6, hsp.bitscore)
        self.assertEqual(0.37, hsp.evalue)
        self.assertEqual(72, hsp.sw_score)
        self.assertAlmostEqual(26.8, hsp.ident_pct)
        self.assertAlmostEqual(54.9, hsp.pos_pct)
        self.assertEqual(71, hsp.aln_span)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('VASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(325, hsp.hit_start)
        self.assertEqual(385, hsp.hit_end)
        self.assertEqual('VVTETTATSVTLTWD------SGNTEPVS---FYG--IQYRAAGTDGPFQEVDGVASTRYSIGGLSPFSEY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(98, hsp.initn_score)
        self.assertEqual(70, hsp.init1_score)
        self.assertEqual(70, hsp.opt_score)
        self.assertEqual(96.2, hsp.z_score)
        self.assertEqual(26.6, hsp.bitscore)
        self.assertEqual(0.37, hsp.evalue)
        self.assertEqual(70, hsp.sw_score)
        self.assertAlmostEqual(32.4, hsp.ident_pct)
        self.assertAlmostEqual(62.2, hsp.pos_pct)
        self.assertEqual(37, hsp.aln_span)
        self.assertEqual(43, hsp.query_start)
        self.assertEqual(80, hsp.query_end)
        self.assertEqual('YEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH', str(hsp.query.seq))
        self.assertEqual(542, hsp.hit_start)
        self.assertEqual(579, hsp.hit_end)
        self.assertEqual('YELVYWAAEDEGQQHKVTFDPTSSYTLEDLKPDTLYH', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('gi|115648048|ref|NP_035343.2|', hit.id)
        self.assertEqual('receptor-type tyrosine-protein phosphatase F precursor [Mus musculus]', hit.description)
        self.assertEqual(1898, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(98, hsp.initn_score)
        self.assertEqual(70, hsp.init1_score)
        self.assertEqual(73, hsp.opt_score)
        self.assertEqual(96.2, hsp.z_score)
        self.assertEqual(26.6, hsp.bitscore)
        self.assertEqual(0.37, hsp.evalue)
        self.assertEqual(73, hsp.sw_score)
        self.assertAlmostEqual(25.6, hsp.ident_pct)
        self.assertAlmostEqual(54.9, hsp.pos_pct)
        self.assertEqual(82, hsp.aln_span)
        self.assertEqual(7, hsp.query_start)
        self.assertEqual(80, hsp.query_end)
        self.assertEqual('PVASDIQARTVVLTWSPPSSL-INGETDESS-----VP---ELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH', str(hsp.query.seq))
        self.assertEqual(497, hsp.hit_start)
        self.assertEqual(579, hsp.hit_end)
        self.assertEqual('PPSPTIQVKTQQGVPAQPADFQANAESDTRIQLSWLLPPQERIVKYELVYWAAEDEGQQHKVTFDPTSSYTLEDLKPDTLYH', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        self.assertEqual(76, hsp.initn_score)
        self.assertEqual(43, hsp.init1_score)
        self.assertEqual(72, hsp.opt_score)
        self.assertEqual(96.2, hsp.z_score)
        self.assertEqual(26.6, hsp.bitscore)
        self.assertEqual(0.37, hsp.evalue)
        self.assertEqual(72, hsp.sw_score)
        self.assertAlmostEqual(26.8, hsp.ident_pct)
        self.assertAlmostEqual(54.9, hsp.pos_pct)
        self.assertEqual(71, hsp.aln_span)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('VASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(325, hsp.hit_start)
        self.assertEqual(385, hsp.hit_end)
        self.assertEqual('VVTETTATSVTLTWD------SGNTEPVS---FYG--IQYRAAGTDGPFQEVDGVASTRYSIGGLSPFSEY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('gi|354481005|ref|XP_003502693.1|', hit.id)
        self.assertEqual('PREDICTED: LOW QUALITY PROTEIN: receptor-type tyrosine-protein phosphatase F-like [Cricetulus griseus]', hit.description)
        self.assertEqual(1898, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(98, hsp.initn_score)
        self.assertEqual(70, hsp.init1_score)
        self.assertEqual(70, hsp.opt_score)
        self.assertEqual(96.2, hsp.z_score)
        self.assertEqual(26.6, hsp.bitscore)
        self.assertEqual(0.37, hsp.evalue)
        self.assertEqual(70, hsp.sw_score)
        self.assertAlmostEqual(32.4, hsp.ident_pct)
        self.assertAlmostEqual(62.2, hsp.pos_pct)
        self.assertEqual(37, hsp.aln_span)
        self.assertEqual(43, hsp.query_start)
        self.assertEqual(80, hsp.query_end)
        self.assertEqual('YEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH', str(hsp.query.seq))
        self.assertEqual(542, hsp.hit_start)
        self.assertEqual(579, hsp.hit_end)
        self.assertEqual('YELVYWAAEDEGQQHKVTFDPTSSYTLEDLKPDTVYH', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, fifth hit
        hit = qresult[4]
        self.assertEqual('gi|328789682|ref|XP_003251305.1|', hit.id)
        self.assertEqual('PREDICTED: LOW QUALITY PROTEIN: twitchin [Apis mellifera]', hit.description)
        self.assertEqual(8619, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual(70, hsp.initn_score)
        self.assertEqual(70, hsp.init1_score)
        self.assertEqual(78, hsp.opt_score)
        self.assertEqual(95.2, hsp.z_score)
        self.assertEqual(28.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(78, hsp.sw_score)
        self.assertAlmostEqual(28.6, hsp.ident_pct)
        self.assertAlmostEqual(54.3, hsp.pos_pct)
        self.assertEqual(70, hsp.aln_span)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(4760, hsp.hit_start)
        self.assertEqual(4823, hsp.hit_end)
        self.assertEqual('ASDVHAEGCTLTWKPP------EDDGGQPIDKYVVEKMDEATGRWVPAGETD-GPQTSLQVEGLTPGHKY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, sixth hit
        hit = qresult[5]
        self.assertEqual('gi|260828627|ref|XP_002609264.1|', hit.id)
        self.assertEqual('hypothetical protein BRAFLDRAFT_124749 [Branchiostoma floridae]', hit.description)
        self.assertEqual(4389, hit.seq_len)
        self.assertEqual(7, len(hit))
        # first qresult, sixth hit, first hsp
        hsp = qresult[5].hsps[0]
        self.assertEqual(81, hsp.initn_score)
        self.assertEqual(73, hsp.init1_score)
        self.assertEqual(97, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(27.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(97, hsp.sw_score)
        self.assertAlmostEqual(21.4, hsp.ident_pct)
        self.assertAlmostEqual(67.1, hsp.pos_pct)
        self.assertEqual(70, hsp.aln_span)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(2241, hsp.hit_start)
        self.assertEqual(2302, hsp.hit_end)
        self.assertEqual('ANAVDSQSIRINWQPPTE-PNGN--------VLGYNIFYTTEGESGNNQQTVGPDDTTYVIEGLRPATQY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, sixth hit, second hsp
        hsp = qresult[5].hsps[1]
        self.assertEqual(177, hsp.initn_score)
        self.assertEqual(55, hsp.init1_score)
        self.assertEqual(90, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(27.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(90, hsp.sw_score)
        self.assertAlmostEqual(30.6, hsp.ident_pct)
        self.assertAlmostEqual(56.9, hsp.pos_pct)
        self.assertEqual(72, hsp.aln_span)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('VASDIQA-RTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(2818, hsp.hit_start)
        self.assertEqual(2881, hsp.hit_end)
        self.assertEqual('VTADGQAPDTVVVTWQSPAET-NGD--------LLGYYIYYQVVGSTETSQAETGPDETTYSISGLRPATEY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, sixth hit, third hsp
        hsp = qresult[5].hsps[2]
        self.assertEqual(196, hsp.initn_score)
        self.assertEqual(61, hsp.init1_score)
        self.assertEqual(84, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(27.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(84, hsp.sw_score)
        self.assertAlmostEqual(27.8, hsp.ident_pct)
        self.assertAlmostEqual(56.9, hsp.pos_pct)
        self.assertEqual(72, hsp.aln_span)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('VASDIQA-RTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(3300, hsp.hit_start)
        self.assertEqual(3363, hsp.hit_end)
        self.assertEqual('VTAEGQAPDTITVTWQSPAET-NGD--------LLGYYIYYQVVGSTEDVRAEAGPEETTYSISGLRPATEY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, sixth hit, fourth hsp
        hsp = qresult[5].hsps[3]
        self.assertEqual(79, hsp.initn_score)
        self.assertEqual(49, hsp.init1_score)
        self.assertEqual(83, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(27.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(83, hsp.sw_score)
        self.assertAlmostEqual(27.9, hsp.ident_pct)
        self.assertAlmostEqual(57.4, hsp.pos_pct)
        self.assertEqual(68, hsp.aln_span)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(80, hsp.query_end)
        self.assertEqual('IQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDYH', str(hsp.query.seq))
        self.assertEqual(3686, hsp.hit_start)
        self.assertEqual(3747, hsp.hit_end)
        self.assertEqual('IDSTTIELQWMPPSP------DEQN-GVIKGYKILYKKVGEEGENEEDAGLLDLMYTLSDLEKWTEYN', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, sixth hit, fifth hsp
        hsp = qresult[5].hsps[4]
        self.assertEqual(100, hsp.initn_score)
        self.assertEqual(50, hsp.init1_score)
        self.assertEqual(81, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(27.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(81, hsp.sw_score)
        self.assertAlmostEqual(25.7, hsp.ident_pct)
        self.assertAlmostEqual(57.1, hsp.pos_pct)
        self.assertEqual(70, hsp.aln_span)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(3398, hsp.hit_start)
        self.assertEqual(3459, hsp.hit_end)
        self.assertEqual('ASSLGSEAIEVSWQPPPQS-NGE--------ILGYRLHYQIVGEESASTQEVEGYETFYLLRGLRPVTEY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, sixth hit, sixth hsp
        hsp = qresult[5].hsps[5]
        self.assertEqual(178, hsp.initn_score)
        self.assertEqual(58, hsp.init1_score)
        self.assertEqual(81, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(27.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(81, hsp.sw_score)
        self.assertAlmostEqual(27.1, hsp.ident_pct)
        self.assertAlmostEqual(55.7, hsp.pos_pct)
        self.assertEqual(70, hsp.aln_span)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('ASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(2145, hsp.hit_start)
        self.assertEqual(2206, hsp.hit_end)
        self.assertEqual('ATPVDPRTVRVEWQPPQQ-PNGE--------IQGYNIYYRTTESDEDALQQAGAQDIFLTLTGLSPFTEY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, sixth hit, seventh hsp
        hsp = qresult[5].hsps[6]
        self.assertEqual(102, hsp.initn_score)
        self.assertEqual(48, hsp.init1_score)
        self.assertEqual(79, hsp.opt_score)
        self.assertEqual(95.1, hsp.z_score)
        self.assertEqual(27.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(79, hsp.sw_score)
        self.assertAlmostEqual(29.4, hsp.ident_pct)
        self.assertAlmostEqual(54.4, hsp.pos_pct)
        self.assertEqual(68, hsp.aln_span)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(79, hsp.query_end)
        self.assertEqual('IQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGE-ETNITLNDLKPAMDY', str(hsp.query.seq))
        self.assertEqual(3497, hsp.hit_start)
        self.assertEqual(3555, hsp.hit_end)
        self.assertEqual('VEPTTITVDWQPPLE-INGV--------LLGYKVIYMPENA-AEFSTVELGPAELSTMLLDLEPATTY', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, seventh hit
        hit = qresult[6]
        self.assertEqual('gi|119220552|ref|NP_689957.3|', hit.id)
        self.assertEqual('protein sidekick-1 isoform 1 [Homo sapiens]', hit.description)
        self.assertEqual(2213, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, seventh hit, first hsp
        hsp = qresult[6].hsps[0]
        self.assertEqual(87, hsp.initn_score)
        self.assertEqual(51, hsp.init1_score)
        self.assertEqual(70, hsp.opt_score)
        self.assertEqual(95.0, hsp.z_score)
        self.assertEqual(26.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(70, hsp.sw_score)
        self.assertAlmostEqual(29.9, hsp.ident_pct)
        self.assertAlmostEqual(58.2, hsp.pos_pct)
        self.assertEqual(67, hsp.aln_span)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(73, hsp.query_end)
        self.assertEqual('VASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYV-GEETNITL-NDL', str(hsp.query.seq))
        self.assertEqual(775, hsp.hit_start)
        self.assertEqual(835, hsp.hit_end)
        self.assertEqual('VASGRTNQSIMVQWQPPP-----ETEHNGV--LRGYILRYRLAGLPGEYQQRNITSPEVNYCLVTDL', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, eighth hit
        hit = qresult[7]
        self.assertEqual('gi|332864595|ref|XP_518946.3|', hit.id)
        self.assertEqual('PREDICTED: protein sidekick-1 [Pan troglodytes]', hit.description)
        self.assertEqual(2213, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, eighth hit, first hsp
        hsp = qresult[7].hsps[0]
        self.assertEqual(68, hsp.initn_score)
        self.assertEqual(45, hsp.init1_score)
        self.assertEqual(76, hsp.opt_score)
        self.assertEqual(95.0, hsp.z_score)
        self.assertEqual(26.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(76, hsp.sw_score)
        self.assertAlmostEqual(32.5, hsp.ident_pct)
        self.assertAlmostEqual(61.0, hsp.pos_pct)
        self.assertEqual(77, hsp.aln_span)
        self.assertEqual(4, hsp.query_start)
        self.assertEqual(80, hsp.query_end)
        self.assertEqual('IVKPVASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYVGEE-TNITLNDLKPAMDYH', str(hsp.query.seq))
        self.assertEqual(674, hsp.hit_start)
        self.assertEqual(740, hsp.hit_end)
        self.assertEqual('LASPNSS--HSHAVVLSWVRP---FDGNS-----PILY-YIVELSENNSPWKVHLSNVGPEMTGITVSGLTPARTYQ', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)
        # first qresult, eighth hit, second hsp
        hsp = qresult[7].hsps[1]
        self.assertEqual(87, hsp.initn_score)
        self.assertEqual(51, hsp.init1_score)
        self.assertEqual(70, hsp.opt_score)
        self.assertEqual(95.0, hsp.z_score)
        self.assertEqual(26.6, hsp.bitscore)
        self.assertEqual(0.43, hsp.evalue)
        self.assertEqual(70, hsp.sw_score)
        self.assertAlmostEqual(29.9, hsp.ident_pct)
        self.assertAlmostEqual(58.2, hsp.pos_pct)
        self.assertEqual(67, hsp.aln_span)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(73, hsp.query_end)
        self.assertEqual('VASDIQARTVVLTWSPPSSLINGETDESSVPELYGYEVLISSTGKDGKYKSVYV-GEETNITL-NDL', str(hsp.query.seq))
        self.assertEqual(775, hsp.hit_start)
        self.assertEqual(835, hsp.hit_end)
        self.assertEqual('VASGRTNQSIMVQWQPPP-----ETEHNGV--LRGYILRYRLAGLPGEYQQRNITSPEVNYCLVTDL', str(hsp.hit.seq))
        self.assertEqual(0, hsp.query_strand)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
