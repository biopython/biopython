# Copyright 2012 by Kai Blin.
# Revisions copyright 2012 by Wibowo Arindrarto
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO hmmer2 text module."""

import unittest
from os import path

from Bio.SearchIO import parse, read


class HmmpfamTests(unittest.TestCase):

    fmt = "hmmer2-text"

    def test_hmmpfam_21(self):
        """Test parsing hmmpfam 2.1 file (text_21_hmmpfam_001.out)."""
        results = parse(path.join("Hmmer", "text_21_hmmpfam_001.out"), self.fmt)
        res = next(results)
        self.assertEqual("roa1_drome", res.id)
        self.assertEqual("<unknown description>", res.description)
        self.assertEqual("hmmpfam", res.program)
        self.assertEqual("2.1.1", res.version)
        self.assertEqual("pfam", res.target)
        self.assertEqual(1, len(res))

        hit = res[0]
        self.assertEqual("SEED", hit.id)
        self.assertEqual("<unknown description>", hit.description)
        self.assertAlmostEqual(146.1, hit.bitscore)
        self.assertAlmostEqual(6.3e-40, hit.evalue)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(77, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(32, hsp.query_start)
        self.assertEqual(103, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertAlmostEqual(71.2, hsp.bitscore)
        self.assertAlmostEqual(2.2e-17, hsp.evalue)
        self.assertEqual(
            "lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG.kelggrklrv",
            hsp.hit.seq,
        )
        self.assertEqual(
            "lf+g+L + +t+e Lk++F+k G iv++ +++D     + t++s+Gf+F+++  ++  + A +    +++++gr+++ ",
            str(hsp.aln_annotation["similarity"]),
        )
        self.assertEqual(
            "LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP",
            hsp.query.seq,
        )

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(77, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(123, hsp.query_start)
        self.assertEqual(194, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertAlmostEqual(75.5, hsp.bitscore)
        self.assertAlmostEqual(1.1e-18, hsp.evalue)
        self.assertEqual(
            "lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv",
            hsp.hit.seq,
        )
        self.assertEqual(
            "lfVg L  d +e+ ++d+F++fG iv+i+iv+D     ketgk +GfaFVeF++++ ++k +     ++l+g+ + v",
            str(hsp.aln_annotation["similarity"]),
        )
        self.assertEqual(
            "LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV",
            hsp.query.seq,
        )

    def test_hmmpfam_22(self):
        """Test parsing hmmpfam 2.2 file (text_22_hmmpfam_001.out)."""
        results = parse(path.join("Hmmer", "text_22_hmmpfam_001.out"), self.fmt)
        res = next(results)
        self.assertEqual("gi|1522636|gb|AAC37060.1|", res.id)
        self.assertEqual(
            "M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]",
            res.description,
        )
        self.assertEqual("[none]", res.accession)
        self.assertEqual("hmmpfam", res.program)
        self.assertEqual("2.2g", res.version)
        self.assertEqual("Pfam", res.target)
        self.assertEqual(1, len(res))

        hit = res[0]
        self.assertEqual("Methylase_M", hit.id)
        self.assertEqual("Type I restriction modification system, M", hit.description)
        self.assertAlmostEqual(-105.2, hit.bitscore)
        self.assertAlmostEqual(0.0022, hit.evalue)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(279, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(279, hsp.query_start)
        self.assertEqual(481, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertAlmostEqual(-105.2, hsp.bitscore)
        self.assertAlmostEqual(0.0022, hsp.evalue)
        self.assertEqual(
            "lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerrieieerktdtesepsldyakledqyeql"
            "ededlekedfyqkkGvFilPsqlFwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdld"
            "fnsnkLgskaqarnetLtelidlfselelgtPmHNG.dfeelgikDlfGDaYEYLLgkFAeneGKsGGeFYTPq"
            "eVSkLiaeiLtigqpsegdfsIYDPAcGSGSLllqaskflgehdgkrnaisyYGQEsn",
            hsp.hit.seq,
        )
        self.assertEqual(
            " ++EL+++  av+   R              L+F K++ dk      +i+         p +   + +++y   "
            "++   ++ ++y ++      + lF++++   e ++  ++++ + +    ++      + +       Glf +++"
            "+  ++ +s+   +ne ++e+i+ +++ +++     G++ +el   D++G +YE L+   Ae   K+ G +YTP "
            "e++  ia+ + i+  ++                  +++ ++    k+n+i +    s+",
            str(hsp.aln_annotation["similarity"]),
        )
        self.assertEqual(
            "NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------GIV---------PRDLLRRTYEDY---"
            "KKSNVLI-NYYDAY-L----KPLFYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSNN"
            "V--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILGYVYEKLINILAEKGQKGLGAYYTPD"
            "EITSYIAKNT-IEPIVVE----------------RFKEIIK--NWKINDINF----ST",
            hsp.query.seq,
        )

    def test_hmmpfam_23(self):
        """Test parsing hmmpfam 2.3 file (text_23_hmmpfam_001.out)."""
        results = parse(path.join("Hmmer", "text_23_hmmpfam_001.out"), self.fmt)
        res = next(results)
        self.assertEqual("gi|90819130|dbj|BAE92499.1|", res.id)
        self.assertEqual("glutamate synthase [Porphyra yezoensis]", res.description)
        self.assertEqual("[none]", res.accession)
        self.assertEqual("hmmpfam", res.program)
        self.assertEqual("2.3.2", res.version)
        self.assertEqual("../Shared/Pfam_fs", res.target)
        self.assertEqual(54, len(res))

        hit = res[0]
        self.assertEqual("Glu_synthase", hit.id)
        self.assertEqual("Conserved region in glutamate synthas", hit.description)
        self.assertAlmostEqual(858.6, hit.bitscore)
        self.assertAlmostEqual(3.6e-255, hit.evalue)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(296, hsp.hit_start)
        self.assertEqual(323, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(649, hsp.query_start)
        self.assertEqual(676, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertAlmostEqual(1.3, hsp.bitscore)
        self.assertAlmostEqual(3, hsp.evalue)
        self.assertEqual("lPwelgLaevhqtLvengLRdrVsLia", hsp.hit.seq)
        self.assertEqual(
            "+P  l++ +vh  L++ gLR + s+ +", str(hsp.aln_annotation["similarity"])
        )
        self.assertEqual("IPPLLAVGAVHHHLINKGLRQEASILV", hsp.query.seq)

        hsp = hit[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(412, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(829, hsp.query_start)
        self.assertEqual(1216, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertAlmostEqual(857.3, hsp.bitscore)
        self.assertAlmostEqual(9e-255, hsp.evalue)

    def test_hmmpfam_23_no_match(self):
        """Test parsing hmmpfam 2.3 file (text_23_hmmpfam_002.out)."""
        results = parse(path.join("Hmmer", "text_23_hmmpfam_002.out"), self.fmt)
        res = next(results)

        self.assertEqual("SEQ0001", res.id)
        self.assertEqual(0, len(res.hits))

        res = next(results)

        self.assertEqual("SEQ0002", res.id)
        self.assertEqual(0, len(res.hits))

    def test_hmmpfam_23_missing_consensus(self):
        """Test parsing hmmpfam 2.3 file (text_23_hmmpfam_003.out)."""
        results = parse(path.join("Hmmer", "text_23_hmmpfam_003.out"), self.fmt)
        res = next(results)

        self.assertEqual("small_input", res.id)
        self.assertEqual("[none]", res.description)
        self.assertEqual("[none]", res.accession)
        self.assertEqual("hmmpfam", res.program)
        self.assertEqual("2.3.2", res.version)
        self.assertEqual(
            "antismash/specific_modules/lantipeptides/ClassIVLanti.hmm", res.target
        )
        self.assertEqual(1, len(res))

        hit = res[0]
        self.assertEqual("ClassIVLanti", hit.id)
        self.assertEqual("Class-IV", hit.description)
        self.assertAlmostEqual(-79.3, hit.bitscore)
        self.assertAlmostEqual(1, hit.evalue)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(66, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(5, hsp.query_start)
        self.assertEqual(20, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertAlmostEqual(-79.3, hsp.bitscore)
        self.assertAlmostEqual(1, hsp.evalue)
        self.assertEqual(len(hsp.query.seq), len(hsp.hit.seq))
        self.assertEqual(len(hsp.query.seq), len(hsp.aln_annotation["similarity"]))
        self.assertEqual(
            "msEEqLKAFiAKvqaDtsLqEqLKaEGADvvaiAKAaGFtitteDLnahiqakeLsdeeLEgvaGg",
            hsp.hit.seq,
        )
        self.assertEqual(
            "        F+                           G  +t   Ln                   ",
            str(hsp.aln_annotation["similarity"]),
        )
        self.assertEqual(
            "-------CFL---------------------------GCLVTNWVLNRS-----------------",
            hsp.query.seq,
        )

    def test_hmmpfam_23_break_in_end_of_seq(self):
        """Test parsing hmmpfam 2.3 file with a line break in the end of seq marker.

        file (text_23_hmmpfam_004.out)
        """
        results = parse(path.join("Hmmer", "text_23_hmmpfam_004.out"), self.fmt)
        res = next(results)
        self.assertEqual("PKSI-KS", res[0].id)
        self.assertEqual("PKSI-FK", res[1].id)

    def test_hmmpfam_24(self):
        """Test parsing hmmpfam 2.4 file (text_24_hmmpfam_001.out)."""
        results = list(parse(path.join("Hmmer", "text_24_hmmpfam_001.out"), self.fmt))
        self.assertEqual(5, len(results))

        # first qresult
        res = results[0]
        self.assertEqual("random_s00", res.id)
        self.assertEqual("[none]", res.accession)
        self.assertEqual("[none]", res.description)
        self.assertEqual("hmmpfam", res.program)
        self.assertEqual("2.4i", res.version)
        self.assertEqual("/home/bow/db/hmmer/Pfam_fs", res.target)
        self.assertEqual(0, len(res))

        # fourth qresult
        res = results[3]
        self.assertEqual("gi|22748937|ref|NP_065801.1|", res.id)
        self.assertEqual("[none]", res.accession)
        self.assertEqual("exportin-5 [Homo sapiens]", res.description)
        self.assertEqual("hmmpfam", res.program)
        self.assertEqual("2.4i", res.version)
        self.assertEqual("/home/bow/db/hmmer/Pfam_fs", res.target)
        self.assertEqual(33, len(res))

        # fourth qresult, first hit
        hit = res[0]
        self.assertEqual("Xpo1", hit.id)
        self.assertEqual("Exportin 1-like protein", hit.description)
        self.assertAlmostEqual(170.1, hit.bitscore)
        self.assertAlmostEqual(5.1e-48, hit.evalue)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        # fourth qresult, first hit, first hsp
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertAlmostEqual(170.1, hsp.bitscore)
        self.assertAlmostEqual(5.1e-148, hsp.evalue)
        self.assertEqual(108, hsp.query_start)
        self.assertEqual(271, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual("ENHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQG--", hsp.query.seq[:40])
        self.assertEqual(
            "+++++  L+++++e++k+ewP++Wp+ + +l  l++++  ",
            str(hsp.aln_annotation["similarity"])[:40],
        )
        self.assertEqual(
            "WVSMSHITA-ENCkLLEILCLLL----NEQELQLGAAECL", hsp.query.seq[-40:]
        )
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(178, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual("pkflrnKLalalaelakqewPsnWpsffpdlvsllsssss", hsp.hit.seq[:40])
        self.assertEqual(
            "W+++++i + ++++ll++l+ lL    +  +l++ A+eCL",
            str(hsp.aln_annotation["similarity"])[-40:],
        )
        self.assertEqual("Wipiglianvnpi.llnllfslLsgpesdpdlreaAveCL", hsp.hit.seq[-40:])

        # fourth qresult, second from last hit
        hit = res[-2]
        self.assertEqual("Rad50_zn_hook", hit.id)
        self.assertEqual("Rad50 zinc hook motif", hit.description)
        self.assertAlmostEqual(2.2, hit.bitscore)
        self.assertAlmostEqual(9.2, hit.evalue)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        # fourth qresult, second from last hit, first hsp
        hsp = hit[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertAlmostEqual(0.8, hsp.bitscore)
        self.assertAlmostEqual(22, hsp.evalue)
        self.assertEqual(20, hsp.query_start)
        self.assertEqual(47, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual("MDPNSTQRYRLEALKFCEEFKE-KCPIC", hsp.query.seq)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(28, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual("galesekaelkkaieeleeeesscCPvC", hsp.hit.seq)

        # fourth qresult, second from last hit, last hsp
        hsp = hit[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertAlmostEqual(1.3, hsp.bitscore)
        self.assertAlmostEqual(16, hsp.evalue)
        self.assertEqual(789, hsp.query_start)
        self.assertEqual(811, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual("EMLAKMAEPFTKALDMLDAEKS", hsp.query.seq)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(22, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual("galesekaelkkaieeleeees", hsp.hit.seq)


class HmmsearchTests(unittest.TestCase):

    fmt = "hmmer2-text"

    def test_hmmsearch_20(self):
        """Test parsing hmmsearch 2.0 file (text_20_hmmsearch_001.out)."""
        res = read(path.join("Hmmer", "text_20_hmmsearch_001.out"), self.fmt)

        # first query
        self.assertEqual("SEED", res.id)
        self.assertEqual("<unknown description>", res.description)
        self.assertEqual("hmmsearch", res.program)
        self.assertEqual("2.0", res.version)
        self.assertEqual("HMM.dbtemp.29591", res.target)
        self.assertEqual(751, len(res))

        # first hit
        hit = res[0]
        self.assertEqual("PAB2_ARATH", hit.id)
        self.assertEqual("P42731 POLYADENYLATE-BINDING PROTEIN 2 (PO", hit.description)
        self.assertAlmostEqual(393.8, hit.bitscore)
        self.assertAlmostEqual(6.1e-114, hit.evalue)
        self.assertEqual(4, hit.domain_obs_num)
        self.assertEqual(4, len(hit))

        # first hit, first hsp
        hsp = hit[0]
        self.assertEqual("SEED", hsp.query_id)
        self.assertEqual("<unknown description>", hsp.query_description)
        self.assertEqual("PAB2_ARATH", hsp.hit_id)
        self.assertEqual(
            "P42731 POLYADENYLATE-BINDING PROTEIN 2 (PO", hsp.hit_description
        )
        self.assertEqual(3, hsp.domain_index)
        self.assertAlmostEqual(109.1, hsp.bitscore)
        self.assertAlmostEqual(3e-28, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(77, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(216, hsp.hit_start)
        self.assertEqual(287, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)

        # first hit, last hsp
        hsp = hit[-1]
        self.assertEqual("SEED", hsp.query_id)
        self.assertEqual("<unknown description>", hsp.query_description)
        self.assertEqual("PAB2_ARATH", hsp.hit_id)
        self.assertEqual(
            "P42731 POLYADENYLATE-BINDING PROTEIN 2 (PO", hsp.hit_description
        )
        self.assertEqual(1, hsp.domain_index)
        self.assertAlmostEqual(92.1, hsp.bitscore)
        self.assertAlmostEqual(3.9e-23, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(77, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(37, hsp.hit_start)
        self.assertEqual(109, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)

        # last hit
        hit = res[-1]
        self.assertEqual("O00369", hit.id)
        self.assertEqual("O00369 L1 ELEMENT L1.20 P40 AND PUTATIVE P", hit.description)
        self.assertAlmostEqual(-23.8, hit.bitscore)
        self.assertAlmostEqual(9.9e02, hit.evalue)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        # last hit, first hsp
        hsp = hit[0]
        self.assertEqual("SEED", hsp.query_id)
        self.assertEqual("<unknown description>", hsp.query_description)
        self.assertEqual("O00369", hsp.hit_id)
        self.assertEqual(
            "O00369 L1 ELEMENT L1.20 P40 AND PUTATIVE P", hsp.hit_description
        )
        self.assertEqual(1, hsp.domain_index)
        self.assertAlmostEqual(-23.8, hsp.bitscore)
        self.assertAlmostEqual(9.9e02, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(77, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(180, hsp.hit_start)
        self.assertEqual(249, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)

    def test_hmmsearch_22(self):
        """Test parsing hmmsearch 2.2 file (text_22_hmmsearch_001.out)."""
        res = read(path.join("Hmmer", "text_22_hmmsearch_001.out"), self.fmt)

        # first query
        self.assertEqual("Peptidase_C1", res.id)
        self.assertEqual("Papain family cysteine protease", res.description)
        self.assertEqual("hmmsearch", res.program)
        self.assertEqual("2.2g", res.version)
        self.assertEqual("cysprot1b.fa", res.target)
        self.assertEqual(4, len(res))

        # first hit
        hit = res[0]
        self.assertEqual("CATL_RAT", hit.id)
        self.assertEqual("<unknown description>", hit.description)
        self.assertAlmostEqual(449.4, hit.bitscore)
        self.assertAlmostEqual(2e-135, hit.evalue)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        # first hit, first hsp
        hsp = hit[0]
        self.assertEqual("Peptidase_C1", hsp.query_id)
        self.assertEqual("Papain family cysteine protease", hsp.query_description)
        self.assertEqual("CATL_RAT", hit.id)
        self.assertEqual("<unknown description>", hit.description)
        self.assertEqual(1, hsp.domain_index)
        self.assertAlmostEqual(449.4, hsp.bitscore)
        self.assertAlmostEqual(2e-135, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(337, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual("lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgr", hsp.query.seq[:40])
        self.assertEqual(
            "IVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi", hsp.query.seq[-40:]
        )
        self.assertEqual(337, len(hsp.query.seq))
        self.assertEqual(
            "+P+++DWRe kg  VtpVK+QG qCGSCWAFSa g lEg+",
            str(hsp.aln_annotation["similarity"])[:40],
        )
        self.assertEqual(
            "+VKNSWG++WG++GY++ia+++n    n+CG+a+ asypi",
            str(hsp.aln_annotation["similarity"])[-40:],
        )
        self.assertEqual(113, hsp.hit_start)
        self.assertEqual(332, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual("IPKTVDWRE-KG-CVTPVKNQG-QCGSCWAFSASGCLEGQ", hsp.hit.seq[:40])
        self.assertEqual("LVKNSWGKEWGMDGYIKIAKDRN----NHCGLATAASYPI", hsp.hit.seq[-40:])
        self.assertEqual(337, len(hsp.hit.seq))

        # last hit
        hit = res[-1]
        self.assertEqual("PAPA_CARPA", hit.id)
        self.assertEqual("<unknown description>", hit.description)
        self.assertAlmostEqual(337.7, hit.bitscore)
        self.assertAlmostEqual(9e-102, hit.evalue)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        # last hit, last hsp
        hsp = hit[-1]
        self.assertEqual("Peptidase_C1", hsp.query_id)
        self.assertEqual("Papain family cysteine protease", hsp.query_description)
        self.assertEqual("PAPA_CARPA", hit.id)
        self.assertEqual("<unknown description>", hit.description)
        self.assertEqual(1, hsp.domain_index)
        self.assertAlmostEqual(337.7, hsp.bitscore)
        self.assertAlmostEqual(9e-102, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(337, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual("lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgr", hsp.query.seq[:40])
        self.assertEqual(
            "IVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi", hsp.query.seq[-40:]
        )
        self.assertEqual(337, len(hsp.query.seq))
        self.assertEqual(
            "+Pe +DWR+ kg aVtpVK+QG +CGSCWAFSav ++Eg+",
            str(hsp.aln_annotation["similarity"])[:40],
        )
        self.assertEqual(
            "++KNSWGt WGEnGY+ri+Rg+++s ++ CG+ ++  yp+",
            str(hsp.aln_annotation["similarity"])[-40:],
        )
        self.assertEqual(133, hsp.hit_start)
        self.assertEqual(343, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual("IPEYVDWRQ-KG-AVTPVKNQG-SCGSCWAFSAVVTIEGI", hsp.hit.seq[:40])
        self.assertEqual("LIKNSWGTGWGENGYIRIKRGTGNS-YGVCGLYTSSFYPV", hsp.hit.seq[-40:])
        self.assertEqual(337, len(hsp.hit.seq))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
