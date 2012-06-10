# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HmmerIO parsers."""


import os
import unittest

from Bio.SearchIO import parse, read

# test case files are in the Blast directory
TEST_DIR = 'Hmmer'
FMT = 'hmmer-text'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HmmerscanCases(unittest.TestCase):

    def test_hmm001(self):
        "Test parsing hmmerscan 3.0 (text_hmm001)"

        xml_file = get_file('text_hmm001.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(0, len(qresult))
        
        # test second result
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|4885477|ref|NP_005359.1|', qresult.id)
        self.assertEqual('myoglobin [Homo sapiens]', qresult.desc)
        self.assertEqual(154, qresult.seq_len)
        self.assertEqual(1, len(qresult))

        hit = qresult[0]
        self.assertEqual('Globin', hit.id)
        self.assertEqual('Globin', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(6e-21, hit.evalue)
        self.assertEqual(74.6, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))
        
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(74.0, hsp.bitscore)
        self.assertEqual(0.2, hsp.bias)
        self.assertEqual(6.7e-25, hsp.evalue_cond)
        self.assertEqual(9.2e-21, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(106, hsp.hit_to)
        self.assertEqual('[.', hsp.hit_endtype)
        self.assertEqual(6, hsp.query_from)
        self.assertEqual(111, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(6, hsp.env_from)
        self.assertEqual(112, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('HHHHHHHHHHHHCHHHHHHHHHHHHHHHHHSGGGGGGGCCCTTTT.HHHHHTSCHHHHHHHHHHHHHHHHHHCTTSHHHHHHHHHHHHHHHHTT-.--HHHHCCHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('qkalvkaswekvkanaeeigaeilkrlfkaypdtkklFkkfgdls.aedlksspkfkahakkvlaaldeavknldnddnlkaalkklgarHakrg.vdpanfklfgeal', \
                hsp.hit.seq.tostring())
        self.assertEqual('EWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKsEDEMKASEDLKKHGATVLTALGGILKK---KGHHEAEIKPLAQSHATKHkIPVKYLEFISECI', \
                hsp.query.seq.tostring())
        self.assertEqual('5789*********************************************************************...6899***********************999998', \
                hsp.alignment_annotation['PP'])

        # test third result
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|126362951:116-221', qresult.id)
        self.assertEqual('leukocyte immunoglobulin-like receptor subfamily B member 1 isoform 2 precursor [Homo sapiens]', qresult.desc)
        self.assertEqual(106, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual('Ig_3', hit.id)
        self.assertEqual('Immunoglobulin domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(1.4e-09, hit.evalue)
        self.assertEqual(38.2, hit.bitscore)
        self.assertEqual(0.4, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(37.6, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(3e-13, hsp.evalue_cond)
        self.assertEqual(2.1e-09, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(72, hsp.hit_to)
        self.assertEqual('[.', hsp.hit_endtype)
        self.assertEqual(8, hsp.query_from)
        self.assertEqual(83, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(8, hsp.env_from)
        self.assertEqual(87, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.94, hsp.acc_avg)
        self.assertEqual('kPvisvspsptvtsggnvtLtCsaeggpppptisWy.....ietppelqgsegssssestLtissvtsedsgtYtCva', \
                hsp.hit.seq.tostring())
        self.assertEqual('KPTLSAQPSPVVNSGGNVILQCDSQVA--FDGFSLCkegedEHPQCLNSQPHARGSSRAIFSVGPVSPSRRWWYRCYA', \
                hsp.query.seq.tostring())
        self.assertEqual('8************************99..78888888****************************************9', \
                hsp.alignment_annotation['PP'])

        hit = qresult[-1]
        self.assertEqual('Ig_2', hit.id)
        self.assertEqual('Immunoglobulin domain', hit.desc)
        self.assertEqual(3.5e-05, hit.evalue)
        self.assertEqual(23.7, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.1, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(6.2e-09, hsp.evalue_cond)
        self.assertEqual(4.3e-05, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(79, hsp.hit_to)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(8, hsp.query_from)
        self.assertEqual(103, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(8, hsp.env_from)
        self.assertEqual(103, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.71, hsp.acc_avg)
        self.assertEqual('kpvlvapp.svvtegenvtLtCsapgnptprvqwykdg.vels......qsqnq........lfipnvsaedsgtYtCra....rnseggktstsveltv', \
                hsp.hit.seq.tostring())
        self.assertEqual('KPTLSAQPsPVVNSGGNVILQCDSQVA-FDGFSLCKEGeDEHPqclnsqP---HargssraiFSVGPVSPSRRWWYRCYAydsnSPYEWSLPSDLLELLV', \
                hsp.query.seq.tostring())
        self.assertEqual('799998885779*************85.899***9988655554443320...134455543444669************88443344588888888766', \
                hsp.alignment_annotation['PP'])

        # test fourth result
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', qresult.id)
        self.assertEqual('exportin-5 [Homo sapiens]', qresult.desc)
        self.assertEqual(1204, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual('Xpo1', hit.id)
        self.assertEqual('Exportin 1-like protein', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(7.8e-34, hit.evalue)
        self.assertEqual(116.6, hit.bitscore)
        self.assertEqual(7.8, hit.bias)
        self.assertEqual(2.8, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))
        
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(116.1, hsp.bitscore)
        self.assertEqual(3.4, hsp.bias)
        self.assertEqual(1.6e-37, hsp.evalue_cond)
        self.assertEqual(1.1e-33, hsp.evalue)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(147, hsp.hit_to)
        self.assertEqual('.]', hsp.hit_endtype)
        self.assertEqual(109, hsp.query_from)
        self.assertEqual(270, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(108, hsp.env_from)
        self.assertEqual(270, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual('HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL', \
                hsp.hit.seq.tostring())
        self.assertEqual('NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL', \
                hsp.query.seq.tostring())
        self.assertEqual('89******************************************************99..79*********************************99*****************************************8889*********************8', \
                hsp.alignment_annotation['PP'])

        hsp = hit[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(-1.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.35, hsp.evalue_cond)
        self.assertEqual(2.4e+03, hsp.evalue)
        self.assertEqual(111, hsp.hit_from)
        self.assertEqual(138, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(498, hsp.query_from)
        self.assertEqual(524, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(495, hsp.env_from)
        self.assertEqual(528, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual('HHCTTS-CHHCHCS.HHHHHCHHCCSCC', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('swvswidiglivnspllsllfqlLndpe', \
                hsp.hit.seq.tostring())
        self.assertEqual('SFVQWEAMTLFLES-VITQMFRTLNREE', \
                hsp.query.seq.tostring())
        self.assertEqual('899*********98.8888899998776', \
                hsp.alignment_annotation['PP'])

        hit = qresult[-1]
        self.assertEqual('IBN_N', hit.id)
        self.assertEqual('Importin-beta N-terminal domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(0.0039, hit.evalue)
        self.assertEqual(16.9, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.7, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))
        
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(14.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.8e-06, hsp.evalue_cond)
        self.assertEqual(0.033, hsp.evalue)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(35, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(32, hsp.env_from)
        self.assertEqual(99, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.87, hsp.acc_avg)
        self.assertEqual('HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll', \
                hsp.hit.seq.tostring())
        self.assertEqual('FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI', \
                hsp.query.seq.tostring())
        self.assertEqual('56788886699*********.6555899******************........999999****99999887', \
                hsp.alignment_annotation['PP'])

        hsp = hit[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(-3.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.2, hsp.evalue_cond)
        self.assertEqual(8e+03, hsp.evalue)
        self.assertEqual(56, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(167, hsp.query_from)
        self.assertEqual(185, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(164, hsp.env_from)
        self.assertEqual(186, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.85, hsp.acc_avg)
        self.assertEqual('HCS-HHHHHHHHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('qqlpeeekelIrnnllnll', \
                hsp.hit.seq.tostring())
        self.assertEqual('QTLPPQRRRDIQQTLTQNM', \
                hsp.query.seq.tostring())
        self.assertEqual('6899*******99998865', \
                hsp.alignment_annotation['PP'])

        # test fifth result
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', qresult.id)
        self.assertEqual('POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]', qresult.desc)
        self.assertEqual(352, qresult.seq_len)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('Pou', hit.id)
        self.assertEqual('Pou domain - N-terminal to homeobox domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(7e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(5e-40, hsp.evalue_cond)
        self.assertEqual(1.4e-36, hsp.evalue)
        self.assertEqual(2, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('.]', hsp.hit_endtype)
        self.assertEqual(132, hsp.query_from)
        self.assertEqual(204, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(130, hsp.env_from)
        self.assertEqual(204, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae', \
                hsp.hit.seq.tostring())
        self.assertEqual('KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD', \
                hsp.query.seq.tostring())
        self.assertEqual('67899******************************************************************96', \
                hsp.alignment_annotation['PP'])

        hit = qresult[1]
        self.assertEqual('Homeobox', hit.id)
        self.assertEqual('Homeobox domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(2.1e-18, hit.evalue)
        self.assertEqual(65.5, hit.bitscore)
        self.assertEqual(1.1, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(64.6, hsp.bitscore)
        self.assertEqual(0.7, hsp.bias)
        self.assertEqual(1.5e-21, hsp.evalue_cond)
        self.assertEqual(4.1e-18, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(56, hsp.hit_to)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(223, hsp.query_from)
        self.assertEqual(279, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(223, hsp.env_from)
        self.assertEqual(279, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual('SS--SS--HHHHHHHHHHCCTSSS--HHHHHHHHHH----HHHHHHHHHHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('rrkRttftkeqleeLeelFeknrypsaeereeLAkklgLterqVkvWFqNrRakekk', \
                hsp.hit.seq.tostring())
        self.assertEqual('KRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQKGKR', \
                hsp.query.seq.tostring())
        self.assertEqual('79****************************************************997', \
                hsp.alignment_annotation['PP'])

        hit = qresult[2]
        self.assertEqual('HTH_31', hit.id)
        self.assertEqual('Helix-turn-helix domain', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.012, hit.evalue)
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(12.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.7e-05, hsp.evalue_cond)
        self.assertEqual(0.16, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(34, hsp.hit_to)
        self.assertEqual('[.', hsp.hit_endtype)
        self.assertEqual(140, hsp.query_from)
        self.assertEqual(180, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(140, hsp.env_from)
        self.assertEqual(183, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.96, hsp.acc_avg)
        self.assertEqual('aLGarLralReraGLtqeevAerlg......vSastlsrlE', \
                hsp.hit.seq.tostring())
        self.assertEqual('QFAKLLKQKRITLGYTQADVGLTLGvlfgkvFSQTTICRFE', \
                hsp.query.seq.tostring())
        self.assertEqual('6999***********************************99', \
                hsp.alignment_annotation['PP'])

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(0.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.19, hsp.evalue_cond)
        self.assertEqual(5.2e+02, hsp.evalue)
        self.assertEqual(38, hsp.hit_from)
        self.assertEqual(61, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(244, hsp.query_from)
        self.assertEqual(267, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(242, hsp.env_from)
        self.assertEqual(269, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual('rgrpsaavlaalaralgldpaera', \
                hsp.hit.seq.tostring())
        self.assertEqual('CPKPSLQQITHIANQLGLEKDVVR', \
                hsp.query.seq.tostring())
        self.assertEqual('678**************9988765', \
                hsp.alignment_annotation['PP'])

        hit = qresult[3]
        self.assertEqual('Homeobox_KN', hit.id)
        self.assertEqual('Homeobox KN domain', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.039, hit.evalue)
        self.assertEqual(13.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(1.6, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(12.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(3.5e-05, hsp.evalue_cond)
        self.assertEqual(0.095, hsp.evalue)
        self.assertEqual(6, hsp.hit_from)
        self.assertEqual(38, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(243, hsp.query_from)
        self.assertEqual(275, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(240, hsp.env_from)
        self.assertEqual(276, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.91, hsp.acc_avg)
        self.assertEqual('hnPYPskevkeelakqTglsrkqidnWFiNaRr', \
                hsp.hit.seq.tostring())
        self.assertEqual('KCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQ', \
                hsp.query.seq.tostring())
        self.assertEqual('56779*************************996', \
                hsp.alignment_annotation['PP'])

        hit = qresult[4]
        self.assertEqual('DUF521', hit.id)
        self.assertEqual('Protein of unknown function (DUF521)', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.14, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(9.4e-05, hsp.evalue_cond)
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(272, hsp.hit_from)
        self.assertEqual(333, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(220, hsp.query_from)
        self.assertEqual(279, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(196, hsp.env_from)
        self.assertEqual(293, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.77, hsp.acc_avg)
        self.assertEqual('adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk', \
                hsp.hit.seq.tostring())
        self.assertEqual('QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR', \
                hsp.query.seq.tostring())
        self.assertEqual('345666667778888899************************99..9999999988876554', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(5, counter)

    def test_hmm002(self):
        "Test parsing hmmerscan 3.0 (text_hmm002)"

        xml_file = get_file('text_hmm002.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_hmm003(self):
        "Test parsing hmmerscan 3.0 (text_hmm003)"

        xml_file = get_file('text_hmm003.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|4885477|ref|NP_005359.1|', qresult.id)
        self.assertEqual('myoglobin [Homo sapiens]', qresult.desc)
        self.assertEqual(154, qresult.seq_len)
        self.assertEqual(1, len(qresult))

        hit = qresult[0]
        self.assertEqual('Globin', hit.id)
        self.assertEqual('Globin', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(6e-21, hit.evalue)
        self.assertEqual(74.6, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))
        
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(74.0, hsp.bitscore)
        self.assertEqual(0.2, hsp.bias)
        self.assertEqual(6.7e-25, hsp.evalue_cond)
        self.assertEqual(9.2e-21, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(106, hsp.hit_to)
        self.assertEqual('[.', hsp.hit_endtype)
        self.assertEqual(6, hsp.query_from)
        self.assertEqual(111, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(6, hsp.env_from)
        self.assertEqual(112, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('HHHHHHHHHHHHCHHHHHHHHHHHHHHHHHSGGGGGGGCCCTTTT.HHHHHTSCHHHHHHHHHHHHHHHHHHCTTSHHHHHHHHHHHHHHHHTT-.--HHHHCCHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('qkalvkaswekvkanaeeigaeilkrlfkaypdtkklFkkfgdls.aedlksspkfkahakkvlaaldeavknldnddnlkaalkklgarHakrg.vdpanfklfgeal', \
                hsp.hit.seq.tostring())
        self.assertEqual('EWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKsEDEMKASEDLKKHGATVLTALGGILKK---KGHHEAEIKPLAQSHATKHkIPVKYLEFISECI', \
                hsp.query.seq.tostring())
        self.assertEqual('5789*********************************************************************...6899***********************999998', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_hmm004(self):
        "Test parsing hmmerscan 3.0 (text_hmm004)"

        xml_file = get_file('text_hmm004.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|126362951:116-221', qresult.id)
        self.assertEqual('leukocyte immunoglobulin-like receptor subfamily B member 1 isoform 2 precursor [Homo sapiens]', qresult.desc)
        self.assertEqual(106, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual('Ig_3', hit.id)
        self.assertEqual('Immunoglobulin domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(1.4e-09, hit.evalue)
        self.assertEqual(38.2, hit.bitscore)
        self.assertEqual(0.4, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(37.6, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(3e-13, hsp.evalue_cond)
        self.assertEqual(2.1e-09, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(72, hsp.hit_to)
        self.assertEqual('[.', hsp.hit_endtype)
        self.assertEqual(8, hsp.query_from)
        self.assertEqual(83, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(8, hsp.env_from)
        self.assertEqual(87, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.94, hsp.acc_avg)
        self.assertEqual('kPvisvspsptvtsggnvtLtCsaeggpppptisWy.....ietppelqgsegssssestLtissvtsedsgtYtCva', \
                hsp.hit.seq.tostring())
        self.assertEqual('KPTLSAQPSPVVNSGGNVILQCDSQVA--FDGFSLCkegedEHPQCLNSQPHARGSSRAIFSVGPVSPSRRWWYRCYA', \
                hsp.query.seq.tostring())
        self.assertEqual('8************************99..78888888****************************************9', \
                hsp.alignment_annotation['PP'])

        hit = qresult[-1]
        self.assertEqual('Ig_2', hit.id)
        self.assertEqual('Immunoglobulin domain', hit.desc)
        self.assertEqual(3.5e-05, hit.evalue)
        self.assertEqual(23.7, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.1, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(6.2e-09, hsp.evalue_cond)
        self.assertEqual(4.3e-05, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(79, hsp.hit_to)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(8, hsp.query_from)
        self.assertEqual(103, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(8, hsp.env_from)
        self.assertEqual(103, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.71, hsp.acc_avg)
        self.assertEqual('kpvlvapp.svvtegenvtLtCsapgnptprvqwykdg.vels......qsqnq........lfipnvsaedsgtYtCra....rnseggktstsveltv', \
                hsp.hit.seq.tostring())
        self.assertEqual('KPTLSAQPsPVVNSGGNVILQCDSQVA-FDGFSLCKEGeDEHPqclnsqP---HargssraiFSVGPVSPSRRWWYRCYAydsnSPYEWSLPSDLLELLV', \
                hsp.query.seq.tostring())
        self.assertEqual('799998885779*************85.899***9988655554443320...134455543444669************88443344588888888766', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_hmm005(self):
        "Test parsing hmmerscan 3.0 (text_hmm005)"

        xml_file = get_file('text_hmm005.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first result
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', qresult.id)
        self.assertEqual('exportin-5 [Homo sapiens]', qresult.desc)
        self.assertEqual(1204, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual('Xpo1', hit.id)
        self.assertEqual('Exportin 1-like protein', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(7.8e-34, hit.evalue)
        self.assertEqual(116.6, hit.bitscore)
        self.assertEqual(7.8, hit.bias)
        self.assertEqual(2.8, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))
        
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(116.1, hsp.bitscore)
        self.assertEqual(3.4, hsp.bias)
        self.assertEqual(1.6e-37, hsp.evalue_cond)
        self.assertEqual(1.1e-33, hsp.evalue)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(147, hsp.hit_to)
        self.assertEqual('.]', hsp.hit_endtype)
        self.assertEqual(109, hsp.query_from)
        self.assertEqual(270, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(108, hsp.env_from)
        self.assertEqual(270, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual('HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL', \
                hsp.hit.seq.tostring())
        self.assertEqual('NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL', \
                hsp.query.seq.tostring())
        self.assertEqual('89******************************************************99..79*********************************99*****************************************8889*********************8', \
                hsp.alignment_annotation['PP'])

        hsp = hit[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(-1.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.35, hsp.evalue_cond)
        self.assertEqual(2.4e+03, hsp.evalue)
        self.assertEqual(111, hsp.hit_from)
        self.assertEqual(138, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(498, hsp.query_from)
        self.assertEqual(524, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(495, hsp.env_from)
        self.assertEqual(528, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual('HHCTTS-CHHCHCS.HHHHHCHHCCSCC', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('swvswidiglivnspllsllfqlLndpe', \
                hsp.hit.seq.tostring())
        self.assertEqual('SFVQWEAMTLFLES-VITQMFRTLNREE', \
                hsp.query.seq.tostring())
        self.assertEqual('899*********98.8888899998776', \
                hsp.alignment_annotation['PP'])

        hit = qresult[-1]
        self.assertEqual('IBN_N', hit.id)
        self.assertEqual('Importin-beta N-terminal domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(0.0039, hit.evalue)
        self.assertEqual(16.9, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.7, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))
        
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(14.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.8e-06, hsp.evalue_cond)
        self.assertEqual(0.033, hsp.evalue)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(35, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(32, hsp.env_from)
        self.assertEqual(99, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.87, hsp.acc_avg)
        self.assertEqual('HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll', \
                hsp.hit.seq.tostring())
        self.assertEqual('FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI', \
                hsp.query.seq.tostring())
        self.assertEqual('56788886699*********.6555899******************........999999****99999887', \
                hsp.alignment_annotation['PP'])

        hsp = hit[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(-3.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.2, hsp.evalue_cond)
        self.assertEqual(8e+03, hsp.evalue)
        self.assertEqual(56, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(167, hsp.query_from)
        self.assertEqual(185, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(164, hsp.env_from)
        self.assertEqual(186, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.85, hsp.acc_avg)
        self.assertEqual('HCS-HHHHHHHHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('qqlpeeekelIrnnllnll', \
                hsp.hit.seq.tostring())
        self.assertEqual('QTLPPQRRRDIQQTLTQNM', \
                hsp.query.seq.tostring())
        self.assertEqual('6899*******99998865', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_hmm006(self):
        "Test parsing hmmerscan 3.0 (text_hmm006)"

        xml_file = get_file('text_hmm006.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', qresult.id)
        self.assertEqual('POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]', qresult.desc)
        self.assertEqual(352, qresult.seq_len)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('Pou', hit.id)
        self.assertEqual('Pou domain - N-terminal to homeobox domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(7e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(5e-40, hsp.evalue_cond)
        self.assertEqual(1.4e-36, hsp.evalue)
        self.assertEqual(2, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('.]', hsp.hit_endtype)
        self.assertEqual(132, hsp.query_from)
        self.assertEqual(204, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(130, hsp.env_from)
        self.assertEqual(204, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae', \
                hsp.hit.seq.tostring())
        self.assertEqual('KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD', \
                hsp.query.seq.tostring())
        self.assertEqual('67899******************************************************************96', \
                hsp.alignment_annotation['PP'])

        hit = qresult[1]
        self.assertEqual('Homeobox', hit.id)
        self.assertEqual('Homeobox domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(2.1e-18, hit.evalue)
        self.assertEqual(65.5, hit.bitscore)
        self.assertEqual(1.1, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(64.6, hsp.bitscore)
        self.assertEqual(0.7, hsp.bias)
        self.assertEqual(1.5e-21, hsp.evalue_cond)
        self.assertEqual(4.1e-18, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(56, hsp.hit_to)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(223, hsp.query_from)
        self.assertEqual(279, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(223, hsp.env_from)
        self.assertEqual(279, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual('SS--SS--HHHHHHHHHHCCTSSS--HHHHHHHHHH----HHHHHHHHHHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('rrkRttftkeqleeLeelFeknrypsaeereeLAkklgLterqVkvWFqNrRakekk', \
                hsp.hit.seq.tostring())
        self.assertEqual('KRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQKGKR', \
                hsp.query.seq.tostring())
        self.assertEqual('79****************************************************997', \
                hsp.alignment_annotation['PP'])

        hit = qresult[2]
        self.assertEqual('HTH_31', hit.id)
        self.assertEqual('Helix-turn-helix domain', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.012, hit.evalue)
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(12.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.7e-05, hsp.evalue_cond)
        self.assertEqual(0.16, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(34, hsp.hit_to)
        self.assertEqual('[.', hsp.hit_endtype)
        self.assertEqual(140, hsp.query_from)
        self.assertEqual(180, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(140, hsp.env_from)
        self.assertEqual(183, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.96, hsp.acc_avg)
        self.assertEqual('aLGarLralReraGLtqeevAerlg......vSastlsrlE', \
                hsp.hit.seq.tostring())
        self.assertEqual('QFAKLLKQKRITLGYTQADVGLTLGvlfgkvFSQTTICRFE', \
                hsp.query.seq.tostring())
        self.assertEqual('6999***********************************99', \
                hsp.alignment_annotation['PP'])

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(0.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.19, hsp.evalue_cond)
        self.assertEqual(5.2e+02, hsp.evalue)
        self.assertEqual(38, hsp.hit_from)
        self.assertEqual(61, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(244, hsp.query_from)
        self.assertEqual(267, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(242, hsp.env_from)
        self.assertEqual(269, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual('rgrpsaavlaalaralgldpaera', \
                hsp.hit.seq.tostring())
        self.assertEqual('CPKPSLQQITHIANQLGLEKDVVR', \
                hsp.query.seq.tostring())
        self.assertEqual('678**************9988765', \
                hsp.alignment_annotation['PP'])

        hit = qresult[3]
        self.assertEqual('Homeobox_KN', hit.id)
        self.assertEqual('Homeobox KN domain', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.039, hit.evalue)
        self.assertEqual(13.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(1.6, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(12.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(3.5e-05, hsp.evalue_cond)
        self.assertEqual(0.095, hsp.evalue)
        self.assertEqual(6, hsp.hit_from)
        self.assertEqual(38, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(243, hsp.query_from)
        self.assertEqual(275, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(240, hsp.env_from)
        self.assertEqual(276, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.91, hsp.acc_avg)
        self.assertEqual('hnPYPskevkeelakqTglsrkqidnWFiNaRr', \
                hsp.hit.seq.tostring())
        self.assertEqual('KCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQ', \
                hsp.query.seq.tostring())
        self.assertEqual('56779*************************996', \
                hsp.alignment_annotation['PP'])

        hit = qresult[4]
        self.assertEqual('DUF521', hit.id)
        self.assertEqual('Protein of unknown function (DUF521)', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.14, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(9.4e-05, hsp.evalue_cond)
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(272, hsp.hit_from)
        self.assertEqual(333, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(220, hsp.query_from)
        self.assertEqual(279, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(196, hsp.env_from)
        self.assertEqual(293, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.77, hsp.acc_avg)
        self.assertEqual('adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk', \
                hsp.hit.seq.tostring())
        self.assertEqual('QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR', \
                hsp.query.seq.tostring())
        self.assertEqual('345666667778888899************************99..9999999988876554', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_hmm007(self):
        "Test parsing hmmerscan 3.0 (text_hmm007)"

        xml_file = get_file('text_hmm007.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', qresult.id)
        self.assertEqual('POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]', qresult.desc)
        self.assertEqual(352, qresult.seq_len)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual('Pou', hit.id)
        self.assertEqual('Pou domain - N-terminal to homeobox domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(7e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(5e-40, hsp.evalue_cond)
        self.assertEqual(1.4e-36, hsp.evalue)
        self.assertEqual(2, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('.]', hsp.hit_endtype)
        self.assertEqual(132, hsp.query_from)
        self.assertEqual(204, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(130, hsp.env_from)
        self.assertEqual(204, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)

        hit = qresult[1]
        self.assertEqual('Homeobox', hit.id)
        self.assertEqual('Homeobox domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(2.1e-18, hit.evalue)
        self.assertEqual(65.5, hit.bitscore)
        self.assertEqual(1.1, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(64.6, hsp.bitscore)
        self.assertEqual(0.7, hsp.bias)
        self.assertEqual(1.5e-21, hsp.evalue_cond)
        self.assertEqual(4.1e-18, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(56, hsp.hit_to)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(223, hsp.query_from)
        self.assertEqual(279, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(223, hsp.env_from)
        self.assertEqual(279, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)

        hit = qresult[2]
        self.assertEqual('HTH_31', hit.id)
        self.assertEqual('Helix-turn-helix domain', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.012, hit.evalue)
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(12.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.7e-05, hsp.evalue_cond)
        self.assertEqual(0.16, hsp.evalue)
        self.assertEqual(0, hsp.hit_from)
        self.assertEqual(34, hsp.hit_to)
        self.assertEqual('[.', hsp.hit_endtype)
        self.assertEqual(140, hsp.query_from)
        self.assertEqual(180, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(140, hsp.env_from)
        self.assertEqual(183, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.96, hsp.acc_avg)

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(0.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.19, hsp.evalue_cond)
        self.assertEqual(5.2e+02, hsp.evalue)
        self.assertEqual(38, hsp.hit_from)
        self.assertEqual(61, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(244, hsp.query_from)
        self.assertEqual(267, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(242, hsp.env_from)
        self.assertEqual(269, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)

        hit = qresult[3]
        self.assertEqual('Homeobox_KN', hit.id)
        self.assertEqual('Homeobox KN domain', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.039, hit.evalue)
        self.assertEqual(13.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(1.6, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(12.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(3.5e-05, hsp.evalue_cond)
        self.assertEqual(0.095, hsp.evalue)
        self.assertEqual(6, hsp.hit_from)
        self.assertEqual(38, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(243, hsp.query_from)
        self.assertEqual(275, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(240, hsp.env_from)
        self.assertEqual(276, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.91, hsp.acc_avg)

        hit = qresult[4]
        self.assertEqual('DUF521', hit.id)
        self.assertEqual('Protein of unknown function (DUF521)', hit.desc)
        self.assertFalse(hit.is_in_inclusion)
        self.assertEqual(0.14, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(9.4e-05, hsp.evalue_cond)
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(272, hsp.hit_from)
        self.assertEqual(333, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(220, hsp.query_from)
        self.assertEqual(279, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(196, hsp.env_from)
        self.assertEqual(293, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.77, hsp.acc_avg)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_hmm008(self):
        "Test parsing hmmerscan 3.0 (text_hmm008)"

        xml_file = get_file('text_hmm008.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first result
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmscan', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/Pfam-A.hmm', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', qresult.id)
        self.assertEqual('exportin-5 [Homo sapiens]', qresult.desc)
        self.assertEqual(1204, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual('Xpo1', hit.id)
        self.assertEqual('Exportin 1-like protein', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(7.8e-34, hit.evalue)
        self.assertEqual(116.6, hit.bitscore)
        self.assertEqual(7.8, hit.bias)
        self.assertEqual(2.8, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))
        
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(116.1, hsp.bitscore)
        self.assertEqual(3.4, hsp.bias)
        self.assertEqual(1.6e-37, hsp.evalue_cond)
        self.assertEqual(1.1e-33, hsp.evalue)
        self.assertEqual(1, hsp.hit_from)
        self.assertEqual(147, hsp.hit_to)
        self.assertEqual('.]', hsp.hit_endtype)
        self.assertEqual(109, hsp.query_from)
        self.assertEqual(270, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(108, hsp.env_from)
        self.assertEqual(270, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual('HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL', \
                hsp.hit.seq.tostring())
        self.assertEqual('NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL', \
                hsp.query.seq.tostring())
        self.assertEqual('89******************************************************99..79*********************************99*****************************************8889*********************8', \
                hsp.alignment_annotation['PP'])

        hsp = hit[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(-1.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.35, hsp.evalue_cond)
        self.assertEqual(2.4e+03, hsp.evalue)
        self.assertEqual(111, hsp.hit_from)
        self.assertEqual(138, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(498, hsp.query_from)
        self.assertEqual(524, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(495, hsp.env_from)
        self.assertEqual(528, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual('HHCTTS-CHHCHCS.HHHHHCHHCCSCC', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('swvswidiglivnspllsllfqlLndpe', \
                hsp.hit.seq.tostring())
        self.assertEqual('SFVQWEAMTLFLES-VITQMFRTLNREE', \
                hsp.query.seq.tostring())
        self.assertEqual('899*********98.8888899998776', \
                hsp.alignment_annotation['PP'])

        hit = qresult[-1]
        self.assertEqual('IBN_N', hit.id)
        self.assertEqual('Importin-beta N-terminal domain', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(0.0039, hit.evalue)
        self.assertEqual(16.9, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.7, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))
        
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(14.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.8e-06, hsp.evalue_cond)
        self.assertEqual(0.033, hsp.evalue)
        self.assertEqual(3, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(35, hsp.query_from)
        self.assertEqual(97, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(32, hsp.env_from)
        self.assertEqual(99, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.87, hsp.acc_avg)
        self.assertEqual('HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll', \
                hsp.hit.seq.tostring())
        self.assertEqual('FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI', \
                hsp.query.seq.tostring())
        self.assertEqual('56788886699*********.6555899******************........999999****99999887', \
                hsp.alignment_annotation['PP'])

        hsp = hit[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_in_inclusion)
        self.assertEqual(-3.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.2, hsp.evalue_cond)
        self.assertEqual(8e+03, hsp.evalue)
        self.assertEqual(56, hsp.hit_from)
        self.assertEqual(74, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(167, hsp.query_from)
        self.assertEqual(185, hsp.query_to)
        self.assertEqual('..', hsp.query_endtype)
        self.assertEqual(164, hsp.env_from)
        self.assertEqual(186, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.85, hsp.acc_avg)
        self.assertEqual('HCS-HHHHHHHHHHHHHHH', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('qqlpeeekelIrnnllnll', \
                hsp.hit.seq.tostring())
        self.assertEqual('QTLPPQRRRDIQQTLTQNM', \
                hsp.query.seq.tostring())
        self.assertEqual('6899*******99998865', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

class HmmersearchCases(unittest.TestCase):

    def test_hmm009(self):
        "Test parsing hmmersearch 3.0 (text_hmm009)"

        xml_file = get_file('text_hmm009.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmsearch', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/uniprot_sprot.fasta', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('globins4', qresult.id)
        self.assertEqual(149, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(1, counter)

    def test_hmm010(self):
        "Test parsing hmmersearch 3.0 (text_hmm010)"

        xml_file = get_file('text_hmm010.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmsearch', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/uniprot_sprot.fasta', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('Pkinase', qresult.id)
        self.assertEqual('PF00069.17', qresult.acc)
        self.assertEqual('Protein kinase domain', qresult.desc)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual(7, len(qresult))

        hit = qresult[0]
        self.assertEqual('sp|Q9WUT3|KS6A2_MOUSE', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase alpha-2 OS', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(58, hsp.hit_from)
        self.assertEqual(317, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(58, hsp.env_from)
        self.assertEqual(317, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF', \
                hsp.hit.seq.tostring())
        self.assertEqual('67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997', \
                hsp.alignment_annotation['PP'])

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.5e-77, hsp.evalue_cond)
        self.assertEqual(1.1e-72, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(414, hsp.hit_from)
        self.assertEqual(671, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(414, hsp.env_from)
        self.assertEqual(671, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI', \
                hsp.hit.seq.tostring())
        self.assertEqual('7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7', \
                hsp.alignment_annotation['PP'])

        hit = qresult[-1]
        self.assertEqual('sp|P18654|KS6A3_MOUSE', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase alpha-3 OS', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(5e-144, hit.evalue)
        self.assertEqual(483.2, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(240.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(8.7e-75, hsp.evalue_cond)
        self.assertEqual(6.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(67, hsp.hit_from)
        self.assertEqual(326, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(67, hsp.env_from)
        self.assertEqual(326, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF', \
                hsp.hit.seq.tostring())
        self.assertEqual('67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886', \
                hsp.alignment_annotation['PP'])

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(241.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.1e-75, hsp.evalue_cond)
        self.assertEqual(3.9e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(421, hsp.hit_from)
        self.assertEqual(678, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(421, hsp.env_from)
        self.assertEqual(678, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI', \
                hsp.hit.seq.tostring())
        self.assertEqual('7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        #self.assertRaises(StopIteration, qresults.next, )
        #self.assertEqual(1, counter)

    def test_hmm011(self):
        "Test parsing hmmersearch 3.0 (text_hmm011)"

        xml_file = get_file('text_hmm011.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmsearch', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/uniprot_sprot.fasta', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('Pkinase', qresult.id)
        self.assertEqual('PF00069.17', qresult.acc)
        self.assertEqual('Protein kinase domain', qresult.desc)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual(7, len(qresult))

        hit = qresult[0]
        self.assertEqual('sp|Q9WUT3|KS6A2_MOUSE', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase alpha-2 OS', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(58, hsp.hit_from)
        self.assertEqual(317, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(58, hsp.env_from)
        self.assertEqual(317, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.5e-77, hsp.evalue_cond)
        self.assertEqual(1.1e-72, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(414, hsp.hit_from)
        self.assertEqual(671, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(414, hsp.env_from)
        self.assertEqual(671, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)

        hit = qresult[-1]
        self.assertEqual('sp|P18654|KS6A3_MOUSE', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase alpha-3 OS', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(5e-144, hit.evalue)
        self.assertEqual(483.2, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(240.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(8.7e-75, hsp.evalue_cond)
        self.assertEqual(6.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(67, hsp.hit_from)
        self.assertEqual(326, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(67, hsp.env_from)
        self.assertEqual(326, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(241.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.1e-75, hsp.evalue_cond)
        self.assertEqual(3.9e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(421, hsp.hit_from)
        self.assertEqual(678, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(421, hsp.env_from)
        self.assertEqual(678, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)

        # test if we've properly finished iteration
        #self.assertRaises(StopIteration, qresults.next, )
        #self.assertEqual(1, counter)

    def test_hmm012(self):
        "Test parsing hmmersearch 3.0 (text_hmm012)"

        xml_file = get_file('text_hmm012.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmsearch', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/uniprot_sprot.fasta', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('Pkinase', qresult.id)
        self.assertEqual('PF00069.17', qresult.acc)
        self.assertEqual('Protein kinase domain', qresult.desc)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual(7, len(qresult))

        hit = qresult[0]
        self.assertEqual('sp|Q9WUT3|KS6A2_MOUSE', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(58, hsp.hit_from)
        self.assertEqual(317, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(58, hsp.env_from)
        self.assertEqual(317, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF', \
                hsp.hit.seq.tostring())
        self.assertEqual('67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997', \
                hsp.alignment_annotation['PP'])

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.5e-77, hsp.evalue_cond)
        self.assertEqual(1.1e-72, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(414, hsp.hit_from)
        self.assertEqual(671, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(414, hsp.env_from)
        self.assertEqual(671, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI', \
                hsp.hit.seq.tostring())
        self.assertEqual('7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7', \
                hsp.alignment_annotation['PP'])

        hit = qresult[-1]
        self.assertEqual('sp|P18654|KS6A3_MOUSE', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(5e-144, hit.evalue)
        self.assertEqual(483.2, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(240.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(8.7e-75, hsp.evalue_cond)
        self.assertEqual(6.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(67, hsp.hit_from)
        self.assertEqual(326, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(67, hsp.env_from)
        self.assertEqual(326, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF', \
                hsp.hit.seq.tostring())
        self.assertEqual('67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886', \
                hsp.alignment_annotation['PP'])

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(241.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.1e-75, hsp.evalue_cond)
        self.assertEqual(3.9e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(421, hsp.hit_from)
        self.assertEqual(678, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(421, hsp.env_from)
        self.assertEqual(678, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI', \
                hsp.hit.seq.tostring())
        self.assertEqual('7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        #self.assertRaises(StopIteration, qresults.next, )
        #self.assertEqual(1, counter)

    def test_hmm013(self):
        "Test parsing hmmersearch 3.0 (text_hmm013)"

        xml_file = get_file('text_hmm013.out')
        qresults = parse(xml_file, FMT)
        counter = 0

        # test first qresult
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmsearch', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/uniprot_sprot.fasta', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('globins4', qresult.id)
        self.assertEqual(149, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults.next()
        counter += 1

        self.assertEqual('hmmsearch', qresult.program)
        self.assertEqual('/home/bow/db/hmmer/uniprot_sprot.fasta', qresult.target)
        self.assertEqual('3.0', qresult.version)
        self.assertEqual('Pkinase', qresult.id)
        self.assertEqual('PF00069.17', qresult.acc)
        self.assertEqual('Protein kinase domain', qresult.desc)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual(7, len(qresult))

        hit = qresult[0]
        self.assertEqual('sp|Q9WUT3|KS6A2_MOUSE', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase alpha-2 OS', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(58, hsp.hit_from)
        self.assertEqual(317, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(58, hsp.env_from)
        self.assertEqual(317, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF', \
                hsp.hit.seq.tostring())
        self.assertEqual('67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997', \
                hsp.alignment_annotation['PP'])

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.5e-77, hsp.evalue_cond)
        self.assertEqual(1.1e-72, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(414, hsp.hit_from)
        self.assertEqual(671, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(414, hsp.env_from)
        self.assertEqual(671, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI', \
                hsp.hit.seq.tostring())
        self.assertEqual('7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7', \
                hsp.alignment_annotation['PP'])

        hit = qresult[-1]
        self.assertEqual('sp|P18654|KS6A3_MOUSE', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase alpha-3 OS', hit.desc)
        self.assertTrue(hit.is_in_inclusion)
        self.assertEqual(5e-144, hit.evalue)
        self.assertEqual(483.2, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(240.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(8.7e-75, hsp.evalue_cond)
        self.assertEqual(6.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(67, hsp.hit_from)
        self.assertEqual(326, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(67, hsp.env_from)
        self.assertEqual(326, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF', \
                hsp.hit.seq.tostring())
        self.assertEqual('67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886', \
                hsp.alignment_annotation['PP'])

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_in_inclusion)
        self.assertEqual(241.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.1e-75, hsp.evalue_cond)
        self.assertEqual(3.9e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_from)
        self.assertEqual(259, hsp.query_to)
        self.assertEqual('[]', hsp.query_endtype)
        self.assertEqual(421, hsp.hit_from)
        self.assertEqual(678, hsp.hit_to)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(421, hsp.env_from)
        self.assertEqual(678, hsp.env_to)
        self.assertEqual('..', hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual('EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG', \
                hsp.alignment_annotation['CS'])
        self.assertEqual('yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl', \
                hsp.query.seq.tostring())
        self.assertEqual('YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI', \
                hsp.hit.seq.tostring())
        self.assertEqual('7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7', \
                hsp.alignment_annotation['PP'])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, qresults.next, )
        self.assertEqual(2, counter)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
