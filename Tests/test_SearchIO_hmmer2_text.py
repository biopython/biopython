# Copyright 2012 by Kai Blin.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from os import path
from Bio.SearchIO import parse

class HmmpfamTests(unittest.TestCase):
    def test_hmmpfam_21(self):
        """Test parsing hmmpfam 2.1 file (text_21_hmmpfam_001.out)"""
        results = parse(path.join("Hmmer", "text_21_hmmpfam_001.out"), "hmmer2-text")
        res = results.next()
        self.assertEqual('roa1_drome', res.id)
        self.assertEqual('<unknown description>', res.description)
        self.assertEqual('hmmpfam', res.program)
        self.assertEqual('2.1.1', res.version)
        self.assertEqual('pfam', res.target)
        self.assertEqual(1, len(res))

        hit = res[0]
        self.assertEqual('SEED', hit.id)
        self.assertEqual('<unknown description>', hit.description)
        self.assertAlmostEqual(146.1, hit.bitscore)
        self.assertAlmostEqual(6.3e-40, hit.evalue)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(77, hsp.hit_end)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(32, hsp.query_start)
        self.assertEqual(103, hsp.query_end)
        self.assertEqual('..', hsp.query_endtype)
        self.assertAlmostEqual(71.2, hsp.bitscore)
        self.assertAlmostEqual(2.2e-17, hsp.evalue)
        self.assertEqual('lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG.kelggrklrv',
                         str(hsp.hit.seq))
        self.assertEqual('LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP',
                         str(hsp.query.seq))

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(77, hsp.hit_end)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(123, hsp.query_start)
        self.assertEqual(194, hsp.query_end)
        self.assertEqual('..', hsp.query_endtype)
        self.assertAlmostEqual(75.5, hsp.bitscore)
        self.assertAlmostEqual(1.1e-18, hsp.evalue)
        self.assertEqual('lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv',
                         str(hsp.hit.seq))
        self.assertEqual('LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV',
                         str(hsp.query.seq))

    def test_hmmpfam_22(self):
        """Test parsing hmmpfam 2.2 file (text_22_hmmpfam_001.out)"""
        results = parse(path.join("Hmmer", "text_22_hmmpfam_001.out"), "hmmer2-text")
        res = results.next()
        self.assertEqual('gi|1522636|gb|AAC37060.1|', res.id)
        self.assertEqual('M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]', res.description)
        self.assertEqual('[none]', res.accession)
        self.assertEqual('hmmpfam', res.program)
        self.assertEqual('2.2g', res.version)
        self.assertEqual('Pfam', res.target)
        self.assertEqual(1, len(res))

        hit = res[0]
        self.assertEqual('Methylase_M', hit.id)
        self.assertEqual('Type I restriction modification system, M', hit.description)
        self.assertAlmostEqual(-105.2, hit.bitscore)
        self.assertAlmostEqual(0.0022, hit.evalue)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(279, hsp.hit_end)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(279, hsp.query_start)
        self.assertEqual(481, hsp.query_end)
        self.assertEqual('..', hsp.query_endtype)
        self.assertAlmostEqual(-105.2, hsp.bitscore)
        self.assertAlmostEqual(0.0022, hsp.evalue)
        self.assertEqual('lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerrieieerktdtesepsldyakledqyeqlededlekedfyqkkGvFilPsqlFwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdldfnsnkLgskaqarnetLtelidlfselelgtPmHNG.dfeelgikDlfGDaYEYLLgkFAeneGKsGGeFYTPqeVSkLiaeiLtigqpsegdfsIYDPAcGSGSLllqaskflgehdgkrnaisyYGQEsn',
                         str(hsp.hit.seq))
        self.assertEqual('NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------GIV---------PRDLLRRTYEDY---KKSNVLI-NYYDAY-L----KPLFYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSNNV--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILGYVYEKLINILAEKGQKGLGAYYTPDEITSYIAKNT-IEPIVVE----------------RFKEIIK--NWKINDINF----ST',
                         str(hsp.query.seq))


    def test_hmmpfam_23(self):
        """Test parsing hmmpfam 2.3 file (text_23_hmmpfam_001.out)"""
        results = parse(path.join("Hmmer", "text_23_hmmpfam_001.out"), "hmmer2-text")
        res = results.next()
        self.assertEqual('gi|90819130|dbj|BAE92499.1|', res.id)
        self.assertEqual('glutamate synthase [Porphyra yezoensis]', res.description)
        self.assertEqual('[none]', res.accession)
        self.assertEqual('hmmpfam', res.program)
        self.assertEqual('2.3.2', res.version)
        self.assertEqual('../Shared/Pfam_fs', res.target)
        self.assertEqual(54, len(res))

        hit = res[0]
        self.assertEqual('Glu_synthase', hit.id)
        self.assertEqual('Conserved region in glutamate synthas', hit.description)
        self.assertAlmostEqual(858.6, hit.bitscore)
        self.assertAlmostEqual(3.6e-255, hit.evalue)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(296, hsp.hit_start)
        self.assertEqual(323, hsp.hit_end)
        self.assertEqual('..', hsp.hit_endtype)
        self.assertEqual(649, hsp.query_start)
        self.assertEqual(676, hsp.query_end)
        self.assertEqual('..', hsp.query_endtype)
        self.assertAlmostEqual(1.3, hsp.bitscore)
        self.assertAlmostEqual(3, hsp.evalue)
        self.assertEqual('lPwelgLaevhqtLvengLRdrVsLia',
                         str(hsp.hit.seq))
        self.assertEqual('IPPLLAVGAVHHHLINKGLRQEASILV',
                         str(hsp.query.seq))

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(412, hsp.hit_end)
        self.assertEqual('[]', hsp.hit_endtype)
        self.assertEqual(829, hsp.query_start)
        self.assertEqual(1216, hsp.query_end)
        self.assertEqual('..', hsp.query_endtype)
        self.assertAlmostEqual(857.3, hsp.bitscore)
        self.assertAlmostEqual(9e-255, hsp.evalue)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
