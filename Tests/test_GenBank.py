# Copyright 2001-2004 by Brad Chapman.  All rights reserved.
# Revisions copyright 2007-2016 by Peter Cock. All rights reserved.
# Revisions copyright 2013 by Kai Blin. All rights reserved.
# Revisions copyright 2015-2016 by Peter Cock.
# Revisions copyright 2019 by Sergio Valqui.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for the GenBank module."""


import os
import sys
import unittest
import warnings
from datetime import datetime

from io import StringIO

from Bio import BiopythonWarning
from Bio import BiopythonParserWarning

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UndefinedSequenceError

# GenBank stuff to test:
from Bio import GenBank


class TestBasics(unittest.TestCase):
    def do_comparison(self, good_record, test_record):
        """Compare two records to see if they are the same.

        This compares the two GenBank records line by line.
        """
        good_handle = StringIO(good_record)
        test_handle = StringIO(test_record)
        while True:
            good_line = good_handle.readline()
            test_line = test_handle.readline()
            if not good_line and not test_line:
                break
            self.assertTrue(good_line, f"Extra info in Test: {test_line!r}")
            self.assertTrue(test_line, f"Extra info in Expected: {good_line!r}")
            test_normalized = " ".join(x for x in test_line.split() if x)
            good_normalized = " ".join(x for x in good_line.split() if x)
            self.assertEqual(test_normalized, good_normalized)

    def test_write_format(self):
        """Test writing to the difference formats."""
        # We only test writing on a subset of the examples:
        filenames = [
            "noref.gb",
            "cor6_6.gb",
            "iro.gb",
            "pri1.gb",
            "arab1.gb",
            "extra_keywords.gb",
            "one_of.gb",
            "origin_line.gb",
        ]
        # don't test writing on protein_refseq, since it is horribly nasty
        # don't test writing on the CONTIG refseq, because the wrapping of
        # locations won't work exactly
        # don't test writing on blank_seq because it lacks a sequence type
        # don't test dbsource_wrap because it is a junky RefSeq file
        record_parser = GenBank.RecordParser(debug_level=0)
        for filename in filenames:
            path = os.path.join("GenBank", filename)
            with open(path) as cur_handle, open(path) as compare_handle:
                iterator = GenBank.Iterator(cur_handle, record_parser)
                compare_iterator = GenBank.Iterator(compare_handle)
                while True:
                    cur_rec = next(iterator)
                    compare_record = next(compare_iterator)
                    if cur_rec is None or compare_record is None:
                        break
                    output_record = str(cur_rec) + "\n"
                    self.do_comparison(compare_record, output_record)

    def test_cleaning_features(self):
        """Test the ability to clean up feature values."""
        gb_parser = GenBank.FeatureParser(
            feature_cleaner=GenBank.utils.FeatureValueCleaner()
        )
        path = "GenBank/arab1.gb"
        with open(path) as handle:
            iterator = GenBank.Iterator(handle, gb_parser)
            first_record = next(iterator)
        # test for cleaning of translation
        translation_feature = first_record.features[1]
        test_trans = translation_feature.qualifiers["translation"][0]
        self.assertNotIn(" ", test_trans, "Did not clean spaces out of the translation")
        self.assertNotIn(
            "\012", test_trans, "Did not clean newlines out of the translation"
        )

    def test_ensembl_locus(self):
        """Test the ENSEMBL locus line."""
        line = "LOCUS       HG531_PATCH 1000000 bp DNA HTG 18-JUN-2011\n"
        s = GenBank.Scanner.GenBankScanner()
        c = GenBank._FeatureConsumer(True)
        s._feed_first_line(c, line)
        self.assertEqual(c.data.name, "HG531_PATCH")
        self.assertEqual(c._expected_size, 1000000)
        line = "LOCUS       HG531_PATCH 759984 bp DNA HTG 18-JUN-2011\n"
        s = GenBank.Scanner.GenBankScanner()
        c = GenBank._FeatureConsumer(True)
        s._feed_first_line(c, line)
        self.assertEqual(c.data.name, "HG531_PATCH")
        self.assertEqual(c._expected_size, 759984)
        line = "LOCUS       HG506_HG1000_1_PATCH 814959 bp DNA HTG 18-JUN-2011\n"
        s = GenBank.Scanner.GenBankScanner()
        c = GenBank._FeatureConsumer(True)
        s._feed_first_line(c, line)
        self.assertEqual(c.data.name, "HG506_HG1000_1_PATCH")
        self.assertEqual(c._expected_size, 814959)
        line = "LOCUS       HG506_HG1000_1_PATCH 1219964 bp DNA HTG 18-JUN-2011\n"
        s = GenBank.Scanner.GenBankScanner()
        c = GenBank._FeatureConsumer(True)
        s._feed_first_line(c, line)
        self.assertEqual(c.data.name, "HG506_HG1000_1_PATCH")
        self.assertEqual(c._expected_size, 1219964)


class TestRecordParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec_parser = GenBank.RecordParser(debug_level=0)

    def perform_record_parser_test(
        self,
        record,
        length,
        locus,
        definition,
        accession,
        titles,
        features,
        tls=None,
        tsa=None,
    ):
        self.assertEqual(len(record.sequence), length)
        self.assertEqual(record.locus, locus)
        self.assertEqual(record.definition, definition)
        self.assertEqual(record.accession, accession)
        self.assertEqual(
            tuple(reference.title for reference in record.references), titles
        )
        self.assertEqual(len(record.features), len(features))
        for feature1, feature2 in zip(record.features, features):
            self.assertEqual(feature1.key, feature2[0])
            self.assertEqual(feature1.location, feature2[1])
            self.assertEqual(len(feature1.qualifiers), len(feature2[2]))
            for qualifier, (key, value) in zip(feature1.qualifiers, feature2[2]):
                self.assertEqual(qualifier.key, key)
                self.assertEqual(qualifier.value, value)
        if tls:
            self.assertEqual(tls, record.tls)
        if tsa:
            self.assertEqual(tsa, record.tsa)

    def test_record_parser_01(self):
        path = "GenBank/noref.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 1622
        locus = "NM_006141"
        definition = (
            "Homo sapiens dynein, cytoplasmic, light intermediate polypeptide 2 "
            "(DNCLI2), mRNA"
        )
        accession = ["NM_006141"]
        titles = ()
        features = [
            (
                "source",
                "1..1622",
                (
                    ("/organism=", '"Homo sapiens"'),
                    ("/db_xref=", '"taxon:9606"'),
                    ("/map=", '"16"'),
                ),
            ),
            (
                "gene",
                "1..1622",
                (
                    ("/gene=", '"DNCLI2"'),
                    ("/note=", '"LIC2"'),
                    ("/db_xref=", '"LocusID:1783"'),
                ),
            ),
            (
                "CDS",
                "7..1485",
                (
                    ("/gene=", '"DNCLI2"'),
                    (
                        "/note=",
                        '"similar to R. norvegicus and G. gallus dynein light '
                        "intermediate chain 2, Swiss-Prot Accession Numbers Q62698 and "
                        'Q90828, respectively"',
                    ),
                    ("/codon_start=", "1"),
                    ("/db_xref=", '"LocusID:1783"'),
                    (
                        "/product=",
                        '"dynein, cytoplasmic, light intermediate polypeptide 2"',
                    ),
                    ("/protein_id=", '"NP_006132.1"'),
                    ("/db_xref=", '"GI:5453634"'),
                    (
                        "/translation=",
                        '"MAPVGVEKKLLLGPNGPAVAAAGDLTSEEEEGQSLWSSILSEVSTRARSKLPSGKNILVFG'
                        "EDGSGKTTLMTKLQGAEHGKKGRGLEYLYLSVHDEDRDDHTRCNVWILDGDLYHKGLLKFAV"
                        "SAESLPETLVIFVADMSRPWTVMESLQKWASVLREHIDKMKIPPEKMRELERKFVKDFQDYM"
                        "EPEEGCQGSPQRRGPLTSGSDEENVALPLGDNVLTHNLGIPVLVVCTKCDAVSVLEKEHDYR"
                        "DEHLDFIQSHLRRFCLQYGAALIYTSVKEEKNLDLLYKYIVHKTYGFHFTTPALVVEKDAVF"
                        "IPAGWDNEKKIAILHENFTTVKPEDAYEDFIVKPPVRKLVHDKELAAEDEQVFLMKQQSLLA"
                        "KQPATPTRASESPARGPSGSPRTQGRGGPASVPSSSPGTSVKKPDPNIKNNAASEGVLASFF"
                        'NSLLSKKTGSPGSPGAGGVQSTAKKSGQKTVLSNVQEELDRMTRKPDSMVTNSSTENEA"',
                    ),
                ),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_02(self):
        path = "GenBank/cor6_6.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
            length = 513
            locus = "ATCOR66M"
            definition = "A.thaliana cor6.6 mRNA"
            accession = ["X55053"]
            titles = (
                "Direct Submission",
                "cDNA sequence analysis and expression of two cold-regulated genes of "
                "Arabidopsis thaliana",
            )
            features = [
                (
                    "source",
                    "1..513",
                    (
                        ("/organism=", '"Arabidopsis thaliana"'),
                        ("/strain=", '"Columbia"'),
                        ("/db_xref=", '"taxon:3702"'),
                    ),
                ),
                ("gene", "50..250", (("/gene=", '"cor6.6"'),)),
                (
                    "CDS",
                    "50..250",
                    (
                        ("/gene=", '"cor6.6"'),
                        ("/note=", '"cold regulated"'),
                        ("/codon_start=", "1"),
                        ("/protein_id=", '"CAA38894.1"'),
                        ("/db_xref=", '"GI:16230"'),
                        ("/db_xref=", '"SWISS-PROT:P31169"'),
                        (
                            "/translation=",
                            '"MSETNKNAFQAGQAAGKAEEKSNVLLDKAKDAAAAAGASAQQAGKSISDAAVGGVNF'
                            'VKDKTGLNK"',
                        ),
                    ),
                ),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )
            record = next(records)
            length = 880
            locus = "ATKIN2"
            definition = "A.thaliana kin2 gene"
            accession = ["X62281"]
            titles = (
                "Direct Submission",
                "Structure and expression of kin2, one of two cold- and ABA-induced "
                "genes of Arabidopsis thaliana",
            )
            features = [
                (
                    "source",
                    "1..880",
                    (
                        ("/organism=", '"Arabidopsis thaliana"'),
                        ("/strain=", '"ssp. L. Heynh, Colombia"'),
                        ("/db_xref=", '"taxon:3702"'),
                    ),
                ),
                ("TATA_signal", "9..20", ()),
                ("exon", "44..160", (("/gene=", '"kin2"'), ("/number=", "1"))),
                ("prim_transcript", "44..>579", (("/gene=", '"kin2"'),)),
                ("mRNA", "join(44..160,320..390,504..>579)", (("/gene=", '"kin2"'),)),
                ("gene", "44..579", (("/gene=", '"kin2"'),)),
                (
                    "CDS",
                    "join(104..160,320..390,504..579)",
                    (
                        ("/gene=", '"kin2"'),
                        ("/codon_start=", "1"),
                        ("/protein_id=", '"CAA44171.1"'),
                        ("/db_xref=", '"GI:16354"'),
                        ("/db_xref=", '"SWISS-PROT:P31169"'),
                        (
                            "/translation=",
                            '"MSETNKNAFQAGQAAGKAERRRAMFCWTRPRMLLLQLELPRNRAGKSISDAAVGGVN'
                            'FVKDKTGLNK"',
                        ),
                    ),
                ),
                ("intron", "161..319", (("/gene=", '"kin2"'), ("/number=", "1"))),
                ("exon", "320..390", (("/gene=", '"kin2"'), ("/number=", "2"))),
                ("intron", "391..503", (("/gene=", '"kin2"'), ("/number=", "2"))),
                ("exon", "504..>579", (("/gene=", '"kin2"'), ("/number=", "3"))),
                ("polyA_signal", "620..625", ()),
                ("polyA_signal", "641..646", ()),
                ("polyA_site", "785", ()),
                ("polyA_site", "800", ()),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )
            record = next(records)
            length = 441
            locus = "BNAKINI"
            definition = "Rapeseed Kin1 protein (kin1) mRNA, complete cds"
            accession = ["M81224"]
            titles = ("Nucleotide sequence of a winter B. napus Kin 1 cDNA",)
            features = [
                (
                    "source",
                    "1..441",
                    (
                        ("/organism=", '"Brassica napus"'),
                        ("/cultivar=", '"Jet neuf"'),
                        ("/db_xref=", '"taxon:3708"'),
                        ("/dev_stage=", '"cold induced"'),
                        ("/tissue_type=", '"leaf"'),
                    ),
                ),
                ("gene", "34..300", (("/gene=", '"kin1"'),)),
                (
                    "CDS",
                    "34..231",
                    (
                        ("/gene=", '"kin1"'),
                        ("/codon_start=", "1"),
                        ("/evidence=", "experimental"),
                        ("/protein_id=", '"AAA32993.1"'),
                        ("/db_xref=", '"GI:167146"'),
                        (
                            "/translation=",
                            '"MADNKQSFQAGQASGRAEEKGNVLMDKVKDAATAAGASAQTAGQKITEAAGGAVNLV'
                            'KEKTGMNK"',
                        ),
                    ),
                ),
                (
                    "polyA_signal",
                    "241..247",
                    (("/gene=", '"kin1"'), ("/note=", '"putative"')),
                ),
                (
                    "polyA_signal",
                    "294..300",
                    (("/gene=", '"kin1"'), ("/note=", '"putative"')),
                ),
                ("polyA_site", "441", (("/gene=", '"kin1"'),)),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )
            record = next(records)
            length = 206
            locus = "ARU237582"
            definition = "Armoracia rusticana csp14 gene (partial), exons 2-3"
            accession = ["AJ237582"]
            titles = ("", "Direct Submission")
            features = [
                (
                    "source",
                    "1..206",
                    (
                        ("/organism=", '"Armoracia rusticana"'),
                        ("/db_xref=", '"taxon:3704"'),
                        ("/country=", '"Russia:Bashkortostan"'),
                    ),
                ),
                ("mRNA", "join(<1..48,143..>206)", (("/gene=", '"csp14"'),)),
                ("exon", "1..48", (("/gene=", '"csp14"'), ("/number=", "2"))),
                ("gene", "1..206", (("/gene=", '"csp14"'),)),
                (
                    "CDS",
                    "join(<1..48,143..>206)",
                    (
                        ("/gene=", '"csp14"'),
                        ("/codon_start=", "2"),
                        ("/product=", '"cold shock protein"'),
                        ("/protein_id=", '"CAB39890.1"'),
                        ("/db_xref=", '"GI:4538893"'),
                        ("/translation=", '"DKAKDAAAAAGASAQQAGKNISDAAAGGVNFVKEKTG"'),
                    ),
                ),
                ("intron", "49..142", (("/gene=", '"csp14"'), ("/number=", "2"))),
                ("exon", "143..206", (("/gene=", '"csp14"'), ("/number=", "3"))),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )
            record = next(records)
            length = 282
            locus = "BRRBIF72"
            definition = "Brassica rapa (clone bif72) kin mRNA, complete cds"
            accession = ["L31939"]
            titles = ("Nucleotide sequences of kin gene in chinese cabbage",)
            features = [
                (
                    "source",
                    "1..282",
                    (
                        ("/organism=", '"Brassica rapa"'),
                        ("/db_xref=", '"taxon:3711"'),
                        ("/dev_stage=", '"flower"'),
                    ),
                ),
                ("gene", "24..221", (("/gene=", '"kin"'),)),
                (
                    "CDS",
                    "24..221",
                    (
                        ("/gene=", '"kin"'),
                        ("/codon_start=", "1"),
                        ("/protein_id=", '"AAA91051.1"'),
                        ("/db_xref=", '"GI:1209262"'),
                        (
                            "/translation=",
                            '"MADNKQSFQAGQAAGRAEEKGNVLLMDKVKDAATAAGALQTAGQKITEAAGGAVNLV'
                            'KEKTGMNK"',
                        ),
                    ),
                ),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )
            record = next(records)
            length = 497
            locus = "AF297471"
            definition = "Brassica napus BN28a (BN28a) gene, complete cds"
            accession = ["AF297471"]
            titles = (
                "BN28a, a low temperature-induced gene of Brassica napus",
                "Direct Submission",
            )
            features = [
                (
                    "source",
                    "1..497",
                    (
                        ("/organism=", '"Brassica napus"'),
                        ("/cultivar=", '"Cascade"'),
                        ("/db_xref=", '"taxon:3708"'),
                    ),
                ),
                (
                    "mRNA",
                    "join(<1..54,241..309,423..>497)",
                    (("/gene=", '"BN28a"'), ("/product=", '"BN28a"')),
                ),
                ("gene", "<1..>497", (("/gene=", '"BN28a"'),)),
                (
                    "CDS",
                    "join(1..54,241..309,423..497)",
                    (
                        ("/gene=", '"BN28a"'),
                        (
                            "/note=",
                            '"low temperature-induced; similar to Brassica napus Kin1 '
                            'in Accession Number M81224"',
                        ),
                        ("/codon_start=", "1"),
                        ("/product=", '"BN28a"'),
                        ("/protein_id=", '"AAG13407.1"'),
                        ("/db_xref=", '"GI:10121869"'),
                        (
                            "/translation=",
                            '"MADNKQSFQAGQAAGRAEEKGNVLMDKVKDAATAAGASAQTAGQKITEAAGGAVNLV'
                            'KEKTGMNK"',
                        ),
                    ),
                ),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )

    def test_record_parser_03(self):
        path = "GenBank/iro.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 1326
        locus = "IRO125195"
        definition = "Homo sapiens mRNA full length insert cDNA clone EUROIMAGE 125195"
        accession = ["AL109817"]
        titles = (
            "The European IMAGE consortium for integrated Molecular analysis "
            "of human gene transcripts",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..1326",
                (
                    ("/organism=", '"Homo sapiens"'),
                    ("/db_xref=", '"taxon:9606"'),
                    ("/chromosome=", '"21"'),
                    ("/clone=", '"IMAGE cDNA clone 125195"'),
                    ("/clone_lib=", '"Soares fetal liver spleen 1NFLS"'),
                    (
                        "/note=",
                        '"contains Alu repeat; '
                        "likely to be be derived from unprocessed nuclear RNA "
                        "or genomic DNA; "
                        "encodes putative exons identical to FTCD; "
                        "formimino transferase cyclodeaminase; "
                        "formimino transferase (EC 2.1.2.5) "
                        '/formimino tetrahydro folate cyclodeaminase (EC 4.3.1.4)"',
                    ),
                ),
            ),
            ("gene", "341..756", (("/gene=", '"FTCD"'),)),
            ("exon", "341..384", (("/gene=", '"FTCD"'), ("/number=", "1"))),
            ("intron", "385..617", (("/gene=", '"FTCD"'), ("/number=", "1"))),
            ("exon", "618..756", (("/gene=", '"FTCD"'), ("/number=", "2"))),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_04(self):
        path = "GenBank/pri1.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 741
        locus = "HUGLUT1"
        definition = "Human fructose transporter (GLUT5) gene, promoter and exon 1"
        accession = ["U05344"]
        titles = (
            "Regulation of expression of the human fructose transporter (GLUT5) by "
            "cyclic AMP",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..741",
                (
                    ("/organism=", '"Homo sapiens"'),
                    ("/db_xref=", '"taxon:9606"'),
                    ("/chromosome=", '"1"'),
                    ("/map=", '"1p31"'),
                    ("/clone=", '"lambda hGT5-157"'),
                    ("/tissue_type=", '"liver"'),
                    (
                        "/clone_lib=",
                        '"partial Hae III/Alu I fetal human liver library in lambda '
                        'Ch4A of Maniatis"',
                    ),
                    ("/dev_stage=", '"fetal"'),
                ),
            ),
            ("repeat_region", "1..73", (("/rpt_family=", '"Alu"'),)),
            ("promoter", "1..513", ()),
            ("5'UTR", "514..609", (("/gene=", '"GLUT5"'),)),
            (
                "exon",
                "514..642",
                (
                    ("/gene=", '"GLUT5"'),
                    ("/number=", "1"),
                    ("/product=", '"fructose transporter"'),
                ),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_05(self):
        path = "GenBank/arab1.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 86436
        locus = "AC007323"
        definition = (
            "Genomic sequence for Arabidopsis thaliana BAC T25K16 from chromosome I, "
            "complete sequence"
        )
        accession = ["AC007323"]
        titles = (
            "Genomic sequence for Arabidopsis thaliana BAC T25K16 from chromosome I",
            "Direct Submission",
            "Direct Submission",
            "Direct Submission",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..86436",
                (
                    ("/organism=", '"Arabidopsis thaliana"'),
                    ("/db_xref=", '"taxon:3702"'),
                    ("/chromosome=", '"1"'),
                    ("/clone=", '"T25K16"'),
                ),
            ),
            (
                "CDS",
                "join(3462..3615,3698..3978,4077..4307,4408..4797,4876..5028,5141..5332)",
                (
                    (
                        "/note=",
                        '"containing similarity to NAM-like proteins gi|3695378"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.1"'),
                    ("/protein_id=", '"AAF26460.1"'),
                    ("/db_xref=", '"GI:6715633"'),
                    (
                        "/translation=",
                        '"MEDQVGFGFRPNDEELVGHYLRNKIEGNTSRDVEVAISEVNICSYDPWNLRFQSKYKSRDA'
                        "MWYFFSRRENNKGNRQSRTTVSGKWKLTGESVEVKDQWGFCSEGFRGKIGHKRVLVFLDGRY"
                        "PDKTKSDWVIHEFHYDLLPEHQKLCNVTLFRFSSYFRLSLLSPMFYTDELMCLPPEILQRTY"
                        "VICRLEYKGDDADILSAYAIDPTPAFVPNMTSSAGSVVNQSRQRNSGSYNTYSEYDSANHGQ"
                        "QFNENSNIMQQQPLQGSFNPLLEYDFANHGGQWLSDYIDLQQQVPYLAPYENESEMIWKHVI"
                        "EENFEFLVDERTSMQQHYSDHRPKKPVSGVLPDDSSDTETGSMIFEDTSSSTDSVGSSDEPG"
                        "HTRIDDIPSLNIIEPLHNYKAQEQPKQQSKEKVISSQKSECEWKMAEDSIKIPPSTNTVKQS"
                        'WIVLENAQWNYLKNMIIGVLLFISVISWIILVG"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join("
                "6617..6953,7266..7351,7464..7603,7916..7998,8087..8166,8273..8368"
                "))",
                (
                    ("/note=", '"hypothetical protein"'),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.2"'),
                    ("/protein_id=", '"AAF26477.1"'),
                    ("/db_xref=", '"GI:6715650"'),
                    (
                        "/translation=",
                        '"MAASEHRCVGCGFRVKSLFIQYSPGNIRLMKCGNCKEVADEYIECERMVCFNHFLSLFGPK'
                        "VYRHVLYNAINPATVNIQVKNYFNSTSRCVVGEIHRQTYLKSPELIIDRSLLLRKSDEESSF"
                        "SDSPVLLSIKVLIGVLSANAAFIISFAIATKGLLNEVSRESLLLQVWEFPMSVIFFVDILLL"
                        "TSNSMALKGQTFKMFSMQIVFCCCYFGISQCKFVFKPVMTESTMTRCIAVCLIAHLIRFLVG"
                        'QIFEPTIFLIQIGSLLQYMSYFFRIV"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(11566..12642)",
                (
                    ("/note=", '"putative RAP2.8 protein gi|3695373"'),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.3"'),
                    ("/protein_id=", '"AAF26476.1"'),
                    ("/db_xref=", '"GI:6715649"'),
                    (
                        "/translation=",
                        '"MDLSLAPTTTTSSDQEQDRDQELTSNIGASSSSGPSGNNNNLPMMMIPPPEKEHMFDKVVT'
                        "PSDVGKLNRLVIPKQHAERYFPLDSSNNQNGTLLNFQDRNGKMWRFRYSYWNSSQSYVMTKG"
                        "WSRFVKEKKLDAGDIVSFQRGIGDESERSKLYIDWRHRPDMSLVQAHQFGNFGFNFNFPTTS"
                        "QYSNRFHPLPEYNSVPIHRGLNIGNHQRSYYNTQRQEFVGYGYGNLAGRCYYTGSPLDHRNI"
                        "VGSEPLVIDSVPVVPGRLTPVMLPPLPPPPSTAGKRLRLFGVNMECGNDYNQQEESWLVPRG"
                        'EIGASSSSSSALRLNLSTDHDDDNDDGDDGDDDQFAKKGKSSLSLNFNP"',
                    ),
                ),
            ),
            (
                "CDS",
                "join("
                "23221..24174,24244..24357,24412..24664,24743..25137,25226..25445,"
                "25527..25711,25783..25905,25994..26478,26564..26730,26814..26983,"
                "27074..27235,27320..27415,27505..28133,28314..28507,28592..28782,"
                "28862..30013,30112..30518,30604..30781"
                ")",
                (
                    (
                        "/note=",
                        '"similar to UFD1 protein emb|CAB10321.1; similar to ESTs '
                        'gb|H36434, gb|AI996152.1"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.4"'),
                    ("/protein_id=", '"AAF26461.1"'),
                    ("/db_xref=", '"GI:6715634"'),
                    (
                        "/translation=",
                        '"MVMEDEPREATIKPSYWLDACEDISCDLIDDLVSEFDPSSVAVNESTDENGVINDFFGGID'
                        "HILDSIKNGGGLPNNGVSDTNSQINEVTVTPQVIAKETVKENGLQKNGGKRDEFSKEEGDKD"
                        "RKRARVCSYQSERSNLSGRGHVNNSREGDRFMNRKRTRNWDEAGNNKKKRECNNYRRDGRDR"
                        "EVRGYWERDKVGSNELVYRSGTWEADHERDVKKVSGGNRECDVKAEENKSKPEERKEKVVEE"
                        "QARRYQLDVLEQAKAKNTIAFLETGAGKTLIAILLIKSVHKDLMSQNRKMLSVFLVPKVPLV"
                        "YQVPPNKKHQAEVIRNQTCFQVGHYCGEMGQDFWDSRRWQREFESKQFLKLTSFFLFSSTQV"
                        "LVMTAQILLNILRHSIIRMETIDLLILDECHHAVKKHPYSLVMSEFYHTTPKDKRPAIFGMT"
                        "ASPVNLKGVSSQVDCAIKIRNLETKLDSTVCTIKDRKELEKHVPMPSEIVVEYDKAATMWSL"
                        "HETIKQMIAAVEEAAQASSRKSKWQFMGARDAGAKDELRQVYGVSERTESDGAANLIHKLRA"
                        "INYTLAELGQWCAYKVGQSFLSALQSDERVNFQVDVKFQESYLSEVVSLLQCELLEGAAAEK"
                        "VAAEVGKPENGNAHDEMEEGELPDDPVVSGGEHVDEVIGAAVADGKVTPKVQSLIKLLLKYQ"
                        "HTADFRAIVFVERVVAALVLPKVRIKVFAELPSLSFIRCASMIGHNNSQEMKSSQMQDTISK"
                        "FRDGHVTLLVATSVAEEGLDIRQCNVVMRFDLAKTVLAYIQSRGRARKPGSDYILMVERYIK"
                        "SFKNYILIFVTTGHQISTDMSTCVTCRGNVSHAAFLRNARNSEETLRKEAIERTDLSHLKDT"
                        "SRLISIDAVPGTVYKVEATGAMVSLNSAVGLVHFYCSQLPGDRYAILRPEFSMEKHEKPGGH"
                        "TEYSCRLQLPCNAPFEILEGPVCSSMRLAQQVDIIVSACKKLHEMGAFTDMLLPDKGSGQDA"
                        "EKADQDDEGEPVPGTARHREFYPEGVADVLKGEWVSSGKEVCESSKLFHLYMYNVRCVDFGS"
                        "SKDPFLSEVSEFAILFGNELDAEVLSMSMDLYVARAMITKASLAFKGSLDITENQLSSLKKF"
                        "HVRLMSIVLDVDVEPSTTPWDPAKAYLFVPVTDNTSMEPIKGINWELVEKITKTTAWDNPLQ"
                        "RARPDVYLGTNERTLGGDRREYGFGKLRHNIVFGQKSHPTYGIRGAVASFDVVRASGLLPVR"
                        "DAFEKEVEEDLSKGKLMMADGCMVAEDLIGKIVTAAHSGKRFYVDSICYDMSAETSFPRKEG"
                        "YLGPLEYNTYADYYKQKIYVVQDRLFFYFLHNLRLLRLYKSSSIMLFIRYGVDLNCKQQPLI"
                        "KGRGVSYCKNLLSPRFEQSGESETVLDKTYYVFLPPELCVVHPLSGSLIRGAQRLPSIMRRV"
                        "ESMLLAVQLKNLISYPIPTSKILEALTAASCQETFCYERAELLGDAYLKWVVSRFLFLKYPQ"
                        "KHEGQLTRMRQQMVSNMVLYQFALVKGLQSYIQADRFAPSRWSAPGVPPVFDEDTKDGGSSF"
                        "FDEEQKPVSEENSDVFEDGEMEDGELEGDLSSYRVLSSKTLADVVEALIGVYYVEGGKIAAN"
                        "HLMKWIGIHVEDDPDEVDGTLKNVNVPESVLKSIDFVGLERALKYEFKEKGLLVEAITHASR"
                        "PSSGVSCYQRLEFVGDAVLDHLITRHLFFTYTSLPPGRLTDLRAAAVNNENFARVAVKHKLH"
                        "LYLRHGSSALEKQVNKIKKQSILFSKSFKCLTVWLLFVFQIREFVKEVQTESSKPGFNSFGL"
                        "GDCKAPKVLGDIVESIAGAIFLDSGKDTTAAWKVFQPLLQPMVTPETLPMHPVRELQERCQQ"
                        "QAEGLEYKASRSGNTATVEVFIDGVQVGVAQNPQKKMAQKLAARNALAALKEKEIAESKEKH"
                        "INNGNAGEDQGENENGNKKNGHQPFTRQTLNDICLRKNWPMPSYRCVKEGGPAHAKRFTFGV"
                        'RVNTSDRGWTDECIGEPMPSVKKAKDSAAVLLLELLNKTFS"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join("
                "31084..31126,31223..31304,31341..31515,31635..31700,31790..31897,"
                "31984..32049,32133..32161,32249..32372"
                "))",
                (
                    (
                        "/note=",
                        '"putative inorganic pyrophosphatase gi|3510259; similar to '
                        'ESTs gb|T42316, gb|AI994042.1, gb|AI994013.1, emb|Z29202"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.5"'),
                    ("/protein_id=", '"AAF26475.1"'),
                    ("/db_xref=", '"GI:6715648"'),
                    (
                        "/translation=",
                        '"MSEETKDNQRLQRPAPRLNERILSSLSRRSVAAHPWHDLEIGPGAPQIFNVVVEITKGSKV'
                        "KYELDKKTGLIKVDRILYSSVVYPHNYGFVPRTLCEDNDPIDVLVIMQEPVLPGCFLRARAI"
                        "GLMPMIDQGEKDDKIIAVCVDDPEYKHYTDIKELPPHRLSEIRRFFEDCILFLQCSSLFISI"
                        'DLSTNKKNENKEVAVNDFLPSESAVEAIQYSMDLYAEYILHTLRR"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join("
                "33694..34029,34103..35173,35269..35349,35432..35701,36326..36387,"
                "36512..36623,36725..36763"
                "))",
                (
                    (
                        "/note=",
                        '"putative late elongated hypocotyl emb|CAA07004; similar to '
                        'ESTS gb|AI993521.1, gb|AA650979"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.6"'),
                    ("/protein_id=", '"AAF26474.1"'),
                    ("/db_xref=", '"GI:6715647"'),
                    (
                        "/translation=",
                        '"MDTNTSGEELLAKARKPYTITKQRERWTEDEHERFLEALRLYGRAWQRIEEHIGTKTAVQI'
                        "RSHAQKFFTKFGKAHSFWFTFQLEKEAEVKGIPVCQALDIEIPPPRPKRKPNTPYPRKPGNN"
                        "GTSSSQVSSAKDAKLVSSASSSQLNQAFLDLEKMPFSEKTSTGKENQDENCSGVSTVNKYPL"
                        "PTKVSGDIETSKTSTVDNAVQDVPKKNKDKDGNDGTTVHSMQNYPWHFHADIVNGNIAKCPQ"
                        "NHPSGMVSQDFMFHPMREETHGHANLQATTASATTTASHQAFPACHSQDDYRSFLQISSTFS"
                        "NLIMSTLLQNPAAHAAATFAASVWPYASVGNSGDSSTPMSSSPPSITAIAAATVAAATAWWA"
                        "SHGLLPVCAPAPITCVPFSTVAVPTPAMTEMDTVENTQPFEKQNTALQDQNLASKSPASSSD"
                        "DSDETGVTKLNADSKTNDDKIEEVVVTAAVHDSNTAQKKNLVDRSSCGSNTPSGSDAETDAL"
                        "DKMEKDKEDVKETDENQPDVIELNNRKIKMRDNNSNNNATTDSWKEVSEEGRIAFQALFARE"
                        "RLPQSFSPPQVAENVNRKQSDTSMPLAPNFKSQDSCAADQEGVVMIGVGTCKSLKTRQTGFK"
                        'PYKRCSMEVKESQVGNINNQSDEKVCKRLRLEGEAST"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join("
                "38600..38756,38838..38989,39111..39516,39915..40031,40377..40579"
                "))",
                (
                    (
                        "/note=",
                        '"similar to Medicago truncatula MtN2 gi|3193308; similar to '
                        'EST gb|H77065"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.7"'),
                    ("/protein_id=", '"AAF26473.1"'),
                    ("/db_xref=", '"GI:6715646"'),
                    (
                        "/translation=",
                        '"MAGDMQGVRVVEKYSPVIVMVMSNVAMGSVNALVKKALDVGVNHMVIGAYRMAISALILVP'
                        "FAYVLERASLMQFFFLLGLSYTSATVSCALVSMLPAITFALALIFRTENVKILKTKAGMLKV"
                        "IGTLICISGALFLTFYKGPQISNSHSHSHGGASHNNNDQDKANNWLLGCLYLTIGTVLLSLW"
                        "MLFQGTLSIKYPCKYSSTCLMSIFAAFQCALLSLYKSRDVNDWIIDDRFVITVIIYAGVVGQ"
                        "AMTTVATTWGIKKLGAVFASAFFPLTLISATLFDFLILHTPLYLGSVIGSLVTITGLYMFLW"
                        'GKNKETESSTALSSGMDNEAQYTTPNKDNDSKSPV"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join("
                "45150..45261,45343..45656,45719..45847,46075..46313,47448..47684,"
                "47777..48554,48638..48868"
                "))",
                (
                    (
                        "/note=",
                        '"putative pyruvate dehydrogenase E1 alpha subunit gi|2454182; '
                        "similar to ESTs emb|Z48417, gb|AW039459.1, gb|T15146, "
                        'emb|Z48416, gb|AF066871, gb|T76832, gb|AI996061.1"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.8"'),
                    ("/protein_id=", '"AAF26472.1"'),
                    ("/db_xref=", '"GI:6715645"'),
                    (
                        "/translation=",
                        '"MATAFAPTKLTATVPLHGSHENRLLLPIRLAPPSSFLGSTRSLSLRRLNHSNATRRSPVVS'
                        "VQEVVKEKQSTNNTSLLITKEEGLELYEDMILGRSFEDMCAQMYYRGKMFGFVHLYNGQEAV"
                        "STGFIKLLTKSDSVVSTYRDHVHALSKGVSARAVMSELFGKVTGCCRGQGGSMHMFSKEHNM"
                        "LGGFAFIGEGIPVATGAAFSSKYRREVLKQDCDDVTVAFFGDGTCNNGQFFECLNMAALYKL"
                        "PIIFVVENNLWAIGMSHLRATSDPEIWKKGPAFGMPGVHVDGMDVLKVREVAKEAVTRARRG"
                        "EGPTLVECETYRFRGHSLADPDELRDAAEKAKYAARDPIAALKKYLIENKLAKEAELKSIEK"
                        "KIDELVEEAVEFADASPQPGRSQLLENVFADPKGFGIGPDGRYRSQPLQIKVSSSELSVLDE"
                        "EKEEEVVKGEAEPNKDSVVSKAEPVKKPRPCELYVCNIPRSYDIAQLLDMFQPFGTVISVEV"
                        "VSRNPQTGESRGSGYVTMGSINSAKIAIASLDGTVRARETKKQEVGGREMRVRYSVDMNPGT"
                        "RRNPEVLNSTPKKILMYESQHKVYVGNLPWFTQPDGLRNHFSKFGTIVSTRVLHDRKTGRNR"
                        'VFAFLSFTSGEERDAALSFNGTVNNMKVAESSSEKVSRRVSRKPTVLLLLQRHLLDTNNV"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join(49986..50039,50121..50333,50585..50656))",
                (
                    (
                        "/note=",
                        '"similar to acidic ribosomal protein p1 gi|2252857; '
                        'similar to ESTs gb|T42111, gb|AI099979, gb|AA728491"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.9"'),
                    ("/protein_id=", '"AAF26471.1"'),
                    ("/db_xref=", '"GI:6715644"'),
                    (
                        "/translation=",
                        '"MSTVGELACSYAVMILEDEGIAITADKIATLVKAAGVSIESYWPMLFAKMAEKRNVTDLIM'
                        'NVGAGGGGGAPVAAAAPAAGGGAAAAPAAEEKKKDEPAEESDGDLGFGLFD"',
                    ),
                ),
            ),
            (
                "CDS",
                "join("
                "51941..52048,52136..52432,52640..52885,53186..53326,53405..54196"
                ")",
                (
                    ("/note=", '"hypothetical protein"'),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.10"'),
                    ("/protein_id=", '"AAF26462.1"'),
                    ("/db_xref=", '"GI:6715635"'),
                    (
                        "/translation=",
                        '"MGKKNGSSSWLTAVKRAFRSPTKKDHSNDVEEDEEKKREKRRWFRKPATQESPVKSSGISP'
                        "PAPQEDSLNVNSKPSPETAPSYATTTPPSNAGKPPSAVVPIATSASKTLAPRRIYYARENYA"
                        "AVVIQTSFRGYLARRALRALKGLVKLQALVRGHNVRKQAKMTLRCMQALVRVQSRVLDQRKR"
                        "LSHDGSRKSAFSDSHAVFESRYLQDLSDRQSMSREGSSAAEDWDDRPHTIDAVKVMLQRRRD"
                        "TALRHDKTNLSQAFSQKMWRTVGNQSTEGHHEVELEEERPKWLDRWMATRPWDKRASSRASV"
                        "DQRVSVKTVEIDTSQPYSRTGAGSPSRGQRPSSPSRTSHHYQSRNNFSATPSPAKSRPILIR"
                        "SASPRCQRDPREDRDRAAYSYTSNTPSLRSNYSFTARSGCSISTTMVNNASLLPNYMASTES"
                        "AKARIRSHSAPRQRPSTPERDRAGLVKKRLSYPVPPPAEYEDNNSLRSPSFKSVAGSHFGGM"
                        'LEQQSNYSSCCTESNGVEISPASTSDFRNWLR"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(57094..58680)",
                (
                    (
                        "/note=",
                        '"putative fatty acid elongase 3-ketoacyl-coA synthase 1 '
                        "gi|4091810; similar to ESTs gb|T42377, gb|N96054, gb|T44368, "
                        'gb|AI999379.1, emb|Z26005"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.11"'),
                    ("/protein_id=", '"AAF26470.1"'),
                    ("/db_xref=", '"GI:6715643"'),
                    (
                        "/translation=",
                        '"MERTNSIEMDRERLTAEMAFRDSSSAVIRIRRRLPDLLTSVKLKYVKLGLHNSCNVTTILF'
                        "FLIILPLTGTVLVQLTGLTFDTFSELWSNQAVQLDTATRLTCLVFLSFVLTLYVANRSKPVY"
                        "LVDFSCYKPEDERKISVDSFLTMTEENGSFTDDTVQFQQRISNRAGLGDETYLPRGITSTPP"
                        "KLNMSEARAEAEAVMFGALDSLFEKTGIKPAEVGILIVNCSLFNPTPSLSAMIVNHYKMRED"
                        "IKSYNLGGMGCSAGLISIDLANNLLKANPNSYAVVVSTENITLNWYFGNDRSMLLCNCIFRM"
                        "GGAAILLSNRRQDRKKSKYSLVNVVRTHKGSDDKNYNCVYQKEDERGTIGVSLARELMSVAG"
                        "DALKTNITTLGPMVLPLSEQLMFLISLVKRKMFKLKVKPYIPDFKLAFEHFCIHAGGRAVLD"
                        "EVQKNLDLKDWHMEPSRMTLHRFGNTSSSSLWYEMAYTEAKGRVKAGDRLWQIAFGSGFKCN"
                        'SAVWKALRPVSTEEMTGNAWAGSIDQYPVKVVQ"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join(59508..59665,61670..61826,63133..63513))",
                (
                    ("/note=", '"hypothetical protein"'),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.12"'),
                    ("/protein_id=", '"AAF26469.1"'),
                    ("/db_xref=", '"GI:6715642"'),
                    (
                        "/translation=",
                        '"MEKRSDSESVEILGDWDSPPPEERIVMVSVPTSPESDYARSNQPKEIESRVSDKETASASG'
                        "EVAARRVLPPWMDPSYEWGGGKWKVDGRKNKNKKEKEKEKEEIIPFKEIIEALLGNSGDKVQ"
                        "QDNKVFEVAPSLHVVELRKTGDDTLEFHKVYFRFNLYQPVQLPLILFVVIRFSMLKIIHYHQ"
                        'FTMAHIKEFVCMWDTHLYKEITNLNIWDTLSSTLVLAIWTVNASHE"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join("
                "64100..64177,64272..64358,64453..64509,64603..64719,64812..64919,"
                "65033..65158,65265..65354,65435..65566,65809..65862,65964..66044,"
                "66152..66259,66380..66451,66537..66599,67026..67214"
                "))",
                (
                    (
                        "/note=",
                        '"similar to wpk4 protein kinase dbj|BAA34675; similar to ESTs '
                        'dbj|AB015122, gb|AI997157.1"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.13"'),
                    ("/protein_id=", '"AAF26468.1"'),
                    ("/db_xref=", '"GI:6715641"'),
                    (
                        "/translation=",
                        '"MSGSRRKATPASRTRVGNYEMGRTLGEGSFAKVKYAKNTVTGDQAAIKILDREKVFRHKMV'
                        "EQLKREISTMKLIKHPNVVEIIEVMASKTKIYIVLELVNGGELFDKIAQQGRLKEDEARRYF"
                        "QQLINAVDYCHSRGVYHRDLKPENLILDANGVLKVSDFGLSAFSRQVREDGLLHTACGTPNY"
                        "VAPEVLSDKGYDGAAADVWSCGVILFVLMAGYLPFDEPNLMTLYKRVRICKAEFSCPPWFSQ"
                        "GAKRVIKRILEPNPITRISIAELLEDEWFKKGYKPPSFDQDDEDITIDDVDAAFSNSKECLV"
                        "TEKKEKPVSMNAFELISSSSEFSLENLFEKQAQLVKKETRFTSQRSASEIMSKMEETAKPLG"
                        "FNVRKDNYKIKMKGDKSGRKGQLSVATEVFEVAPSLHVVELRKTGGDTLEFHKVCDSFYKNF"
                        'SSGLKDVVWNTDAAAEEQKQ"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join(69831..69987,70534..70670,70743..71357,71644..71700))",
                (
                    (
                        "/note=",
                        '"similar to ataxia-telangiectasia group D protein pir|A49618"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.14"'),
                    ("/protein_id=", '"AAF26467.1"'),
                    ("/db_xref=", '"GI:6715640"'),
                    (
                        "/translation=",
                        '"MVSDLPLDEDDIALLKSPYCDDGGDEDVNSAPNIFTYDNVPLKKRHYLGTSDTFRSFEPLN'
                        "EHACIVCDIADDGVVPCSGNECPLAVHRKCVELDCEDPATFYCPYCWFKEQATRSTALRTRG"
                        "VAAAKTLVQYGCSELRSGDIVMTRENSQLENGSDNSLPMQLHENLHQLQELVKHLKARNSQL"
                        "DESTDQFIDMEKSCGEAYAVVNDQPKRVLWTVNEEKMLREGVEKFSDTINKNMPWKKILEMG"
                        "KGIFHTTRNSSDLKDKWRNMVRIIILIWLRSRLTSSSSSQRSEIKMERERNAGVMKKMSPTG"
                        'TIQRLEFVGWYL"',
                    ),
                ),
            ),
            (
                "CDS",
                "join("
                "72285..72371,72789..72865,72989..73097,73190..73442,73524..73585"
                ")",
                (
                    (
                        "/note=",
                        '"similar to SYT gi|2252866; similar to ESTs emb|F14390, '
                        'gb|H36066, emb|F14391"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.15"'),
                    ("/protein_id=", '"AAF26463.1"'),
                    ("/db_xref=", '"GI:6715636"'),
                    (
                        "/translation=",
                        '"MQQQQSPQMFPMVPSIPPANNITTEQIQKYLDENKKLIMAIMENQNLGKLAECAQYQALLQ'
                        "KNLMYLAAIADAQPPPPTPGPSPSTAVAAQMATPHSGMQPPSYFMQHPQASPAGIFAPRGPL"
                        "QFGSPLQFQDPQQQQQIHQQAMQGHMGIRPMGMTNNGMQHAMQQPETGLGGNVGLRGGKQDG"
                        'ADGQGKDDGK"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join(73807..73990,74036..74145))",
                (
                    (
                        "/note=",
                        '"similar to stress-induced protein OZI1 precursor pir|S59544; '
                        'similar to EST gb|AI995719.1"',
                    ),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.16"'),
                    ("/protein_id=", '"AAF26466.1"'),
                    ("/db_xref=", '"GI:6715639"'),
                    (
                        "/translation=",
                        '"MASGGKAKYIIGALIGSFGISYIFDKVISDNKIFGGKDDLNGYLLVKISGTTPGTVSNKEW'
                        'WAATDEKFQAWPRTAGPPVVMNPISRQNFIVKTRPE"',
                    ),
                ),
            ),
            (
                "CDS",
                "join(75335..76249,76516..76653,76733..76982,77015..77148)",
                (
                    ("/note=", '"putative reverse transcriptase gb|AAD17395"'),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.17"'),
                    ("/protein_id=", '"AAF26464.1"'),
                    ("/db_xref=", '"GI:6715637"'),
                    (
                        "/translation=",
                        '"MKEDRRLPHKRDAFQFLKTKAAYVIVIVLTYAFGYFSAYHYHQPLQQQLPPSTTAVETTKP'
                        "QVCSIDNFRVTTPCGNLVPPELIRQTVIDRIFNGTSPYIDFPPPHAKKFLRPKRIKGWGSYG"
                        "AVFENLIRRVKPKTIVEVGSFLGASAIHMANLTRRLGLEETQILCVDDFRGWPGFRDRFKDM"
                        "ALVNGDVLLMYQFMQNVVISDFSGSILPVPFSTGSALEKLCEWGVTADLVEIDAGHDFNSAW"
                        "ADINRAVRILRPGGVIFGHDYFTAADNRGVRRAVNLFAEINRLKVKTDGQHWVIDSVKVINK"
                        "GTRFAISKTVAKIKEDANQWFFAQVLENQDLVNEQAVHISVKVLRGFLRDEHGKVLIHARRS"
                        "FASVHSKLDATFLCWQWAMESMKSLRVDKIIFASEDNDLIGAVTRLPSWPSYKFQIHFLLGE"
                        'LIRSSNLGAHLIAKSVTMEDRRQSYVATGFPFWLKHLFEKERSIA"',
                    ),
                ),
            ),
            (
                "CDS",
                "complement(join(82723..82738,82751..83373,83586..84581))",
                (
                    ("/note=", '"putative cytochrome P450 gi|3831440"'),
                    ("/codon_start=", "1"),
                    ("/evidence=", "not_experimental"),
                    ("/product=", '"T25K16.18"'),
                    ("/protein_id=", '"AAF26465.1"'),
                    ("/db_xref=", '"GI:6715638"'),
                    (
                        "/translation=",
                        '"MFSLNMRTEIESLWVFALASKFNIYMQQHFASLLVAIAITWFTITIVFWSTPGGPAWGKYF'
                        "FTRRFISLDYNRKYKNLIPGPRGFPLVGSMSLRSSHVAHQRIASVAEMSNAKRLMAFSLGDT"
                        "KVVVTCHPAVAKEILNSSVFADRPVDETAYGLMFNRAMGFAPNGTYWRTLRRLGSNHLFNPK"
                        "QIKQSEDQRRVIATQMVNAFARNPKSACAVRDLLKTASLCNMMGLVFGREYELESNNNLESE"
                        "CLKGLVEEGYDLLGTLNWTDHLPWLAGLDFQQIRFRCSQLVPKVNLLLSRIIHEQRAATGNF"
                        "LDMLLSLQGSEKLSESDMVAVLWEMIFRGTDTVAVLVEWVLARIVMHPKVQLTVHDELDRVV"
                        "GRSRTVDESDLPSLTYLTAMIKEVLRLHPPGPLLSWARLSITDTSVDGYHVPAGTTAMVNMW"
                        "AIARDPHVWEDPLEFKPERFVAKEGEAEFSVFGSDLRLAPFGSGKRVCPGKNLGLTTVSFWV"
                        'ATLLHEFEWLPSVEANPPDLSEVLRLSCEMACPLIVNVSSRRKIIAWMF"',
                    ),
                ),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_06(self):
        path = "GenBank/protein_refseq.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 182
        locus = "NP_034640"
        definition = "interferon beta, fibroblast [Mus musculus]"
        accession = ["NP_034640"]
        titles = (
            "structure and expression of a cloned cdna for mouse interferon-beta",
        )
        features = [
            (
                "source",
                "1..182",
                (
                    ("/organism=", '"Mus musculus"'),
                    ("/db_xref=", '"taxon:10090"'),
                    ("/chromosome=", '"4"'),
                    ("/map=", '"4 42.6 cM"'),
                ),
            ),
            ("Protein", "1..182", (("/product=", '"interferon beta, fibroblast"'),)),
            ("sig_peptide", "1..21", ()),
            (
                "Region",
                "1..182",
                (
                    ("/region_name=", '"Interferon alpha/beta domain"'),
                    ("/db_xref=", '"CDD:pfam00143"'),
                    ("/note=", '"interferon"'),
                ),
            ),
            ("mat_peptide", "22..182", (("/product=", '"ifn-beta"'),)),
            (
                "Region",
                "56..170",
                (
                    ("/region_name=", '"Interferon alpha, beta and delta."'),
                    ("/db_xref=", '"CDD:IFabd"'),
                    ("/note=", '"IFabd"'),
                ),
            ),
            (
                "CDS",
                "1..182",
                (
                    ("/gene=", '"Ifnb"'),
                    ("/db_xref=", '"LocusID:15977"'),
                    ("/db_xref=", '"MGD:MGI:107657"'),
                    ("/coded_by=", '"NM_010510.1:21..569"'),
                ),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_07(self):
        path = "GenBank/extra_keywords.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 154329
        locus = "DMBR25B3"
        definition = "Drosophila melanogaster BAC clone BACR25B3"
        accession = ["AL138972"]
        titles = (
            "Sequencing the distal X chromosome of Drosophila melanogaster",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..154329",
                (
                    ("/organism=", '"Drosophila melanogaster"'),
                    ("/db_xref=", '"taxon:7227"'),
                    ("/clone=", '"BAC BACR25B3"'),
                ),
            ),
            (
                "gene",
                "complement(22148..27773)",
                (("/gene=", '"EG:BACR25B3.11"'), ("/note=", "")),
            ),
            (
                "CDS",
                "complement(join("
                "22148..22299,22375..22791,22860..23560,23630..24555,24616..24888,"
                "25024..25178,26677..27009,27623..27773"
                "))",
                (
                    ("/gene=", '"EG:BACR25B3.11"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genefinder'', version:''084'',"
                        " score:''105.71''); /prediction=(method:''genscan'',"
                        " version:''1.0''); /match=(desc:''BASEMENT MEMBRANE-SPECIFIC"
                        " HEPARAN SULFATE PROTEOGLYCAN CORE PROTEIN PRECURSOR (HSPG)"
                        " (PERLECAN) (PLC)'', species:''Homo sapiens (Human)'',"
                        " ranges:(query:24292..24549,"
                        " target:SWISS-PROT::P98160:3713..3628, score:''201.00''),"
                        " (query:24016..24291, target:SWISS-PROT::P98160:3815..3724,"
                        " score:''139.00''), (query:23857..24006,"
                        " target:SWISS-PROT::P98160:3866..3817, score:''99.00''),"
                        " (query:24052..24327, target:SWISS-PROT::P98160:4059..3968,"
                        " score:''143.00''), (query:24046..24312,"
                        " target:SWISS-PROT::P98160:4341..4253, score:''116.00''),"
                        " (query:23806..23901, target:SWISS-PROT::P98160:4177..4146,"
                        " score:''76.00''), (query:23203..23382,"
                        " target:SWISS-PROT::P98160:4062..4003, score:''116.00''),"
                        " (query:22523..22777, target:SWISS-PROT::P98160:4288..4204,"
                        " score:''112.00''), (query:22235..22300,"
                        " target:SWISS-PROT::P98160:4358..4337, score:''64.00'')),"
                        " method:''blastx'', version:''1.4.9'');"
                        " /match=(desc:''GM03359.5prime GM Drosophila melanogaster"
                        " ovary BlueScript Drosophila melanogaster cDNA clone GM03359"
                        " 5prime, mRNA sequence'', species:''Drosophila melanogaster"
                        " (fruit fly)'', ranges:(query:25024..25235,"
                        " target:EMBL::AA801707:438..227, score:''1024.00''),"
                        " (query:24851..24898, target:EMBL::AA801707:476..429,"
                        " score:''204.00'')), method:''blastn'', version:''1.4.9'');"
                        " /match=(desc:''LD08615.5prime LD Drosophila melanogaster"
                        " embryo BlueScript Drosophila melanogaster cDNA clone LD08615"
                        " 5prime, mRNA sequence'', species:''Drosophila melanogaster"
                        " (fruit fly)'', ranges:(query:24629..24727,"
                        " target:EMBL::AA264808:99..1, score:''495.00''),"
                        " (query:24417..24566, target:EMBL::AA264808:250..101,"
                        " score:''687.00''), (query:24048..24420,"
                        " target:EMBL::AA264808:618..246, score:''1847.00''),"
                        " (query:23986..24036, target:EMBL::AA264808:678..628,"
                        " score:''237.00'')), method:''blastn'', version:''1.4.9'');"
                        " /match=(desc:''HL02745.5prime HL Drosophila melanogaster head"
                        " BlueScript Drosophila melanogaster cDNA clone HL02745 5prime,"
                        " mRNA sequence'', species:''Drosophila melanogaster (fruit"
                        " fly)'', ranges:(query:23944..24045,"
                        " target:EMBL::AA697546:103..2, score:''510.00''),"
                        " (query:23630..23943, target:EMBL::AA697546:416..103,"
                        " score:''1570.00''), (query:23419..23561,"
                        " target:EMBL::AA697546:558..416, score:''715.00''),"
                        " (query:23306..23417, target:EMBL::AA697546:670..559,"
                        " score:''524.00''), (query:23280..23316,"
                        " target:EMBL::AA697546:695..659, score:''167.00'')),"
                        " method:''blastn'', version:''1.4.9'');"
                        " /match=(desc:''GM08137.5prime GM Drosophila melanogaster"
                        " ovary BlueScript Drosophila melanogaster cDNA clone GM08137"
                        " 5prime, mRNA sequence'', species:''Drosophila melanogaster"
                        " (fruit fly)'', ranges:(query:23235..23278,"
                        " target:EMBL::AA696682:44..1, score:''139.00''),"
                        " (query:22986..23251, target:EMBL::AA696682:294..29,"
                        " score:''1321.00'')), method:''blastn'', version:''1.4.9'')\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72284.1"'),
                    ("/db_xref=", '"GI:6946669"'),
                    (
                        "/translation=",
                        '"MACNCNQSMIYQSNERRDYNCPGAPQYPYNRFKGGVSLKDTPCMVLYICADFKSSKLSSAK'
                        "PIISGPATTRAPAISYVCQPNDFKCVSHPHTCVRANMVCDGIYDCTDHSDEFNCIAGKGSGK"
                        "SESNSGSGSFKRWKKSPEQGRRSLAKAVKNRKLRKRSFAKSRDYSLKLDDQSSNLRAGESTD"
                        "VECYSSDDTYTDVVWERSDGAPLSNNVRQVGNRLVISNVSPSDAGNYVCKCKTDEGDLYTTS"
                        "YKLEVEDQPHELKSSKIVYAKVGANADLQCGADESRQPTYRWSRQYGQLQAGRSLMNEKLSL"
                        "DSVQANDAGTYICTAQYADGETADFPNILVVTGAIPQFRQEPRSYMSFPTLPNSSFKFNFEL"
                        "TFRPENGDGLLLFNGQTRGSGDYIALSLKDRYAEFRFDFGGKPMLVRAEEPLALNEWHTVRV"
                        "SRFKRDGYIQVDEQHPVAFPTLQQIPQLDLIEDLYIGGVPNWELLPADAVSQQVGFVGCISR"
                        "LTLQGRTVELIREAKYKEGITDCRPCAQGPCQNKGVCLESQTEQAYTCICQPGWTGRDCAIE"
                        "GTQCTPGVCGAGRCENTENDMECLCPLNRSGDRCQYNEILNEHSLNFKGNSFAAYGTPKVTK"
                        "VNITLSVRPASLEDSVILYTAESTLPSGDYLALVLRGGHAELLINTAARLDPVVVRSAEPLP"
                        "LNRWTRIEIRRRLGEGILRVGDGPERKAKAPGSDRILSLKTHLYVGGYDRSTVKVNRDVNIT"
                        "KGFDGCISRLYNFQKPVNLLADIKDAANIQSCGETNMIGGDEDSDNEPPVPPPTPDVHENEL"
                        "QPYAMAPCASDPCENGGSCSEQEDVAVCSCPFGFSGKHCQEHLQLGFNASFRGDGYVELNRS"
                        "HFQPALEQSYTSMGIVFTTNKPNGLLFWWGQEAGEEYTGQDFIAAAVVDGYVEYSMRLDGEE"
                        "AVIRNSDIRVDNGERHIVIAKRDENTAILEVDRMLHSGETRPTSKKSMKLPGNVFVGGAPDL"
                        'EVFTGFRYKHNLNGCIVVVEGETVGQINLSSAAVNGVNANVCPA"',
                    ),
                ),
            ),
            ("gene", "complement(29926..33978)", (("/gene=", '"EG:BACR25B3.10"'),)),
            (
                "CDS",
                "complement(join("
                "29926..30108,30270..30519,30617..31076,31197..31591,31659..31836,"
                "32324..32634,32686..33289,33533..33713,33817..33978"
                "))",
                (
                    ("/gene=", '"EG:BACR25B3.10"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genefinder'', version:''084'',"
                        " score:''98.50''); /prediction=(method:''genscan'',"
                        " version:''1.0''); /match=(desc:''BASEMENT MEMBRANE-SPECIFIC"
                        " HEPARAN SULFATE PROTEOGLYCAN CORE PROTEIN PRECURSOR (HSPG)"
                        " (PERLECAN) (PLC)'', species:''Homo sapiens (Human)'',"
                        " ranges:(query:33540..33716,"
                        " target:SWISS-PROT::P98160:2716..2658, score:''113.00''),"
                        " (query:32859..32963, target:SWISS-PROT::P98160:3341..3307,"
                        " score:''63.00''), (query:33150..33215,"
                        " target:SWISS-PROT::P98160:3530..3509, score:''73.00''),"
                        " (query:32973..33089, target:SWISS-PROT::P98160:3588..3550,"
                        " score:''71.00''), (query:32358..32567,"
                        " target:SWISS-PROT::P98160:3650..3581, score:''107.00''),"
                        " (query:31222..31323, target:SWISS-PROT::P98160:2620..2587,"
                        " score:''80.00''), (query:31489..31572,"
                        " target:SWISS-PROT::P98160:3387..3360, score:''72.00''),"
                        " (query:31495..31593, target:SWISS-PROT::P98160:3575..3543,"
                        " score:''60.00'')), method:''blastx'', version:''1.4.9'');"
                        " /match=(desc:''GM02481.5prime GM Drosophila melanogaster"
                        " ovary BlueScript Drosophila melanogaster cDNA clone GM02481"
                        " 5prime, mRNA sequence'', species:''Drosophila melanogaster"
                        " (fruit fly)'', ranges:(query:30008..30036,"
                        " target:EMBL::AA695253:29..1, score:''145.00''),"
                        " (query:29549..30004, target:EMBL::AA695253:487..32,"
                        " score:''2262.00'')), method:''blastn'', version:''1.4.9'')\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72285.1"'),
                    ("/db_xref=", '"GI:6946670"'),
                    (
                        "/translation=",
                        '"MFLATLDTNDPTDIGTEDPVLTQIIVSIQKPEITIVPVGGSMTLSCSGRMRWSNSPVIVNW'
                        "YKENSRLPENVEVQGGNLYLYDLQVSDSGVYICQAVNNETASVFKDTVSITITKKDQLSPAE"
                        "IVNLPSHVTFEEYVNNEIICEVLGNPAPRVTWARVDGHADAQSTRTYDNRLIFDSPRKSDEG"
                        "RYRCQAENDQNRDEKYVIVYVQSNPPQPPPQQDRLYITPEEINGLAGESFQLNCQFTSVASL"
                        "RYDWSHNGRSLSSSPARNVEIRGNTLEVRDASESDSGVYTCVAYDVRTRRNFTESARVNIDR"
                        "REEQPFGVLMRMMILTDSLINHSNKPIIESLEQNILIIQGEDYSITCEASGSPYPSIKWAKV"
                        "HDFMPENVHISGNVLTIYGARFENRGVYSCVAENDHGSDLSSTSIDIEPRERPSVKIVSAPL"
                        "QTFSVGAPASLYCTVEGIPDPTVEWVRVDGQPLSPRHKIQSPGYMVIDDIQLEDSGDYECRA"
                        "KNIVGEATGVATITVQEPTLVQIIPDNRDLRLTEGDELSLTCVGSGVPNPEVEWVNEMALKR"
                        "DLYSPPSNTAILKIYRVTKADAGIYTCHGKNEAGSDEAHVRVEVQERRGDIGGVDDDSDRDP"
                        "INYNPPQQQNPGIHQPGSNQLLATDIGDNVTLTCDMFQPLNTRWERVDGAPLPRNAYTIKNR"
                        "LEIVRVEQQNLGQYRCNGIGRDGNVKTYFVKELVLMPLPRIRFYPNIPLTVEAGQNLDVHCQ"
                        "VENVRPEDVHWSTDNNRPLPSSVRIVGSVLRFVSITQAAAGEYRCSAFNQYGNRSQIARVAV"
                        "KKPADFHQVPQSQLQRHREGENIQLQCTVTDQYGVRAQDNVEFNWFRDDRRPLPNNARTDSQ"
                        'ILVLTNLRPEDAGRYICNSYDVDRGQQLPEVSIDLQVLSE"',
                    ),
                ),
            ),
            ("gene", "complement(36119..56153)", (("/gene=", '"EG:BACR25B3.1"'),)),
            (
                "CDS",
                "complement(join("
                "36119..37213,37281..39517,39656..40042,40345..40434,40519..40612,"
                "40681..40814,41546..41620,41855..42085,42188..42415,42751..42876,"
                "43604..43837,44241..44438,44812..44928,45148..45233,45661..45793,"
                "45976..46125,46518..46688,47222..47315,47683..47831,48411..48878,"
                "49437..49562,49763..49876,49971..50102,50319..50441,50827..50937,"
                "52849..52966,56031..56153))",
                (
                    ("/gene=", '"EG:BACR25B3.1"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genscan'', version:''1.0''); "
                        "/prediction=(method:''genefinder'', version:''084''); "
                        "/match=(desc:''LOW-DENSITY LIPOPROTEIN RECEPTOR-RELATED PROTEIN PRECURSOR (LRP)'', "
                        "species:''Caenorhabditis elegans'', "
                        "ranges:"
                        "(query:50831..50941, target:SWISS-PROT::Q04833:1221..1185, score:''95.00''), "
                        "(query:50840..51025, target:SWISS-PROT::Q04833:2865..2804, score:''102.00''), "
                        "(query:50828..50935, target:SWISS-PROT::Q04833:3788..3753, score:''119.00''), "
                        "(query:50323..50394, target:SWISS-PROT::Q04833:3706..3683, score:''77.00''), "
                        "(query:50326..50433, target:SWISS-PROT::Q04833:1263..1228, score:''120.00''), "
                        "(query:49948..50079, target:SWISS-PROT::Q04833:2917..2874, score:''88.00''), "
                        "(query:49432..49587, target:SWISS-PROT::Q04833:4085..4034, score:''102.00''), "
                        "(query:49429..49560, target:SWISS-PROT::Q04833:3915..3872, score:''97.00''), "
                        "(query:48622..48720, target:SWISS-PROT::Q04833:1302..1270, score:''99.00''), "
                        "(query:47698..47799, target:SWISS-PROT::Q04833:3996..3963, score:''88.00''), "
                        "(query:47686..47775, target:SWISS-PROT::Q04833:3835..3806, score:''59.00''), "
                        "(query:47692..47787, target:SWISS-PROT::Q04833:4041..4010, score:''83.00''), "
                        "(query:47229..47315, target:SWISS-PROT::Q04833:3742..3714, score:''88.00''), "
                        "(query:47220..47312, target:SWISS-PROT::Q04833:3829..3799, score:''67.00''), "
                        "(query:47232..47318, target:SWISS-PROT::Q04833:3866..3838, score:''78.00''), "
                        "(query:46552..46656, target:SWISS-PROT::Q04833:1344..1310, score:''95.00''), "
                        "(query:46543..46650, target:SWISS-PROT::Q04833:3951..3916, score:''98.00''), "
                        "(query:45983..46129, target:SWISS-PROT::Q04833:2870..2822, score:''82.00''), "
                        "(query:45971..46096, target:SWISS-PROT::Q04833:4089..4048, score:''82.00''), "
                        "(query:45678..45764, target:SWISS-PROT::Q04833:3666..3638, score:''80.00''), "
                        "(query:45128..45238, target:SWISS-PROT::Q04833:94..58, score:''100.00''), "
                        "(query:45158..45268, target:SWISS-PROT::Q04833:3990..3954, score:''80.00''), "
                        "(query:44263..44379, target:SWISS-PROT::Q04833:85..47, score:''77.00''), "
                        "(query:44251..44367, target:SWISS-PROT::Q04833:3995..3957, score:''100.00''), "
                        "(query:43605..43688, target:SWISS-PROT::Q04833:2994..2967, score:''84.00''), "
                        "(query:42764..42877, target:SWISS-PROT::Q04833:2951..2914, score:''77.00''), "
                        "(query:42180..42377, target:SWISS-PROT::Q04833:260..195, score:''148.00''), "
                        "(query:42234..42419, target:SWISS-PROT::Q04833:3199..3138, score:''106.00''), "
                        "(query:39807..40013, target:SWISS-PROT::Q04833:2901..2833, score:''167.00''), "
                        "(query:39645..39857, target:SWISS-PROT::Q04833:3138..3068, score:''151.00''), "
                        "(query:39846..40046, target:SWISS-PROT::Q04833:3241..3175, score:''132.00''), "
                        "(query:39654..39866, target:SWISS-PROT::Q04833:3913..3843, score:''201.00'')), "
                        "method:''blastx'', version:''1.4.9''); "
                        "/match=(desc:''LOW-DENSITY LIPOPROTEIN RECEPTOR-RELATED PROTEIN 2 PRECURSOR (MEGALIN) (GLYCOPROTEIN 330)'', "
                        "species:''Homo sapiens (Human)'', "
                        "ranges:"
                        "(query:50834..50935, target:SWISS-PROT::P98164:2733..2700, score:''99.00''), "
                        "(query:50840..50947, target:SWISS-PROT::P98164:3063..3028, score:''94.00''), "
                        "(query:50831..50926, target:SWISS-PROT::P98164:3918..3887, score:''102.00''), "
                        "(query:50326..50433, target:SWISS-PROT::P98164:1222..1187, score:''107.00''), "
                        "(query:50302..50394, target:SWISS-PROT::P98164:3762..3732, score:''91.00''), "
                        "(query:49773..49904, target:SWISS-PROT::P98164:2939..2896, score:''90.00''), "
                        "(query:49438..49578, target:SWISS-PROT::P98164:217..171, score:''116.00''), "
                        "(query:49429..49545, target:SWISS-PROT::P98164:3796..3758, score:''108.00''), "
                        "(query:48622..48720, target:SWISS-PROT::P98164:3544..3512, score:''94.00''), "
                        "(query:48595..48708, target:SWISS-PROT::P98164:3720..3683, score:''86.00''), "
                        "(query:47701..47814, target:SWISS-PROT::P98164:2817..2780, score:''90.00''), "
                        "(query:47692..47799, target:SWISS-PROT::P98164:3674..3639, score:''60.00''), "
                        "(query:47217..47366, target:SWISS-PROT::P98164:3716..3667, score:''96.00''), "
                        "(query:46543..46647, target:SWISS-PROT::P98164:1101..1067, score:''107.00''), "
                        "(query:46552..46656, target:SWISS-PROT::P98164:3873..3839, score:''84.00''), "
                        "(query:45989..46126, target:SWISS-PROT::P98164:3832..3787, score:''98.00''), "
                        "(query:45149..45274, target:SWISS-PROT::P98164:2775..2734, score:''99.00''), "
                        "(query:44780..44893, target:SWISS-PROT::P98164:268..231, score:''76.00''), "
                        "(query:44813..44905, target:SWISS-PROT::P98164:1223..1193, score:''73.00''), "
                        "(query:44251..44361, target:SWISS-PROT::P98164:3630..3594, score:''119.00''), "
                        "(query:43602..43700, target:SWISS-PROT::P98164:179..147, score:''97.00''), "
                        "(query:43674..43781, target:SWISS-PROT::P98164:191..156, score:''90.00''), "
                        "(query:43584..43685, target:SWISS-PROT::P98164:1107..1074, score:''89.00''), "
                        "(query:42758..42865, target:SWISS-PROT::P98164:1264..1229, score:''79.00''), "
                        "(query:42204..42413, target:SWISS-PROT::P98164:2810..2741, score:''136.00''), "
                        "(query:42189..42377, target:SWISS-PROT::P98164:3027..2965, score:''125.00''), "
                        "(query:42186..42293, target:SWISS-PROT::P98164:3110..3075, score:''109.00''), "
                        "(query:42198..42389, target:SWISS-PROT::P98164:3584..3521, score:''137.00''), "
                        "(query:42309..42422, target:SWISS-PROT::P98164:3793..3756, score:''95.00''), "
                        "(query:39654..39791, target:SWISS-PROT::P98164:63..18, score:''132.00''), "
                        "(query:39786..40049, target:SWISS-PROT::P98164:1183..1096, score:''230.00''), "
                        "(query:39657..39890, target:SWISS-PROT::P98164:3109..3032, score:''200.00''), "
                        "(query:39780..39983, target:SWISS-PROT::P98164:3756..3689, score:''194.00''), "
                        "(query:39618..39761, target:SWISS-PROT::P98164:3845..3798, score:''105.00''), "
                        "(query:39651..39779, target:SWISS-PROT::P98164:3964..3922, score:''128.00'')), "
                        "method:''blastx'', version:''1.4.9''); "
                        "/match=(desc:''GM06086.5prime GM Drosophila melanogaster ovary BlueScript Drosophila melanogaster cDNA clone GM06086 5prime, mRNA sequence'', "
                        "species:''Drosophila melanogaster (fruit fly)'', "
                        "ranges:"
                        "(query:50852..51290, target:EMBL::AA802674:672..234, score:''2195.00'')), "
                        "method:''blastn'', version:''1.4.9''); "
                        "/match=(desc:''SD04592.5prime SD Drosophila melanogaster Schneider L2 cell culture pOT2 Drosophila melanogaster cDNA clone SD04592 5prime, mRNA sequence'', "
                        "species:''Drosophila melanogaster (fruit fly)'', "
                        "ranges:(query:37280..37708, target:EMBL::AI532939:429..1, score:''2136.00''), "
                        "(query:37097..37217, target:EMBL::AI532939:545..425, score:''569.00'')), "
                        "method:''blastn'', version:''1.4.9''); "
                        "/match=(desc:''GH03622.5prime GH Drosophila melanogaster head pOT2 Drosophila melanogaster cDNA clone GH03622 5prime, mRNA sequence'', "
                        "species:''Drosophila melanogaster (fruit fly)'', "
                        "ranges:(query:36446..37075, target:EMBL::AI063674:1..630, score:''3150.00'')), "
                        "method:''blastn'', version:''1.4.9''); "
                        "EST embl|AA802674|AA802674 comes from the 5' UTR\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72286.1"'),
                    ("/db_xref=", '"GI:6946671"'),
                    (
                        "/translation=",
                        '"MLLLQLLLQLLLLGKLLLGKTPPTVFGFRLLFAAFRFPLSLHFPHRMHDHFFVRGDTHSCG'
                        "WKNSTTFTIRISAIYRYLNQCQANEFRCNNGDCIDARKRCNNVSDCSEGEDENEECPAACSG"
                        "MEYQCRDGTRCISVSQQCDGHSDCSDGDDEEHCDGIVPKLRYTCPKGKFTCRDLSCISIVHR"
                        "CDGRADCPNDRSDEEGCPCLYDKWQCDDGTCIAKELLCNGNIDCPEDISDERYCEGGYDSEE"
                        "CRFDEFHCGTGECIPMRQVCDNIYDCNDYSDEVNCVEGEEEDRVGIPIGHQPWRPASKHDDW"
                        "LHEMDTSEYQVYQPSNVYEKANSQNPCASNQFRCTTSNVCIPLHLRCDGFYHCNDMSDEKSC"
                        "EQYQRHTTTRRPLTLATPTSRITTQGPGLLERRNTTTATEASRWPWATKTTTIATTTSNPIT"
                        "TVGVANSPPQTCLENIEFACHNRDCISIESVCDGIPDCGRNEDEDDALCKCSGDKYKCQRGG"
                        "GCIPKSQVCDGKPQCHDRSDESACHLHGRLNKTRLGVKCLESQYQCGDGSCISGYKRCNGIH"
                        "DCADASDEYNCIYDYEDTYDTDPNNNPLNECDILEFECDYSQCLPLEKKCDGYADCEDMSDE"
                        "LECQSYTDHCLESEFECDSYCLPRDQLCNGIPNCQDGSDERNCTFCREDAYLCNTGECVADN"
                        "QRCNGIADCADGSDERHCARIYCPPNKLACNGTCVSRRIKCDGIRDCLDGYDEMYCPETNNH"
                        "YPTQNVNVIRPKLGPNPIPKSCRPHEWQCANLECIDSSLQCNEIKDCSDGSDEELSVCFGTA"
                        "TTRLKPSDCSPEQFYCDESCYNRSVRCNGHVDCSDGSDEVGCSLPCPQHQCPSGRCYTESER"
                        "CDRHRHCEDGSDEANCTAILCKDNEFLCFDRQFCINATQQCDGYYDCRDFSDEQNCIGCYAN"
                        "QFRCNNGDCVSGSAPCNGYSECSDHSDELNCGGTQECLPNQFRCNSGQCVSSSVRCNGRTDC"
                        "QDSSDEQNCGHRHTEVSQGLETTGVFTTSTTSTTAMTPLRIICPPTSFKCENGPCISLGLKC"
                        "NGRVDCPYDGSDEADCGQISNDIDPADSNDRRPNQLNLKTYPDSQIIKESREVIFRCRDEGP"
                        "ARAKVKWSRPGGRPLPPGFTDRNGRLEIPNIRVEDAGTYVCEAVGYASYIPGQQVTVNLNVE"
                        "RSWGENKYEEIRSNRIRYGTVPHIDLEFFGLDNDVGSRPESACTEYQATCMNGECIDKSSIC"
                        "DGNPDCSDASDEQSCSLGLKCQPNQFMCSNSKCVDRTWRCDGENDCGDNSDETSCDPEPSGA"
                        "PCRYNEFQCRSGHCIPKSFQCDNVPDCTDGTDEVGCMAPLPIRPPPQSVSLLEYEVLELTCV"
                        "ATGTPTPTIVWRLNWGHVPDKCESKSYGGTGTLRCPDMRPQDSGAYSCEIINTRGTHFVNPD"
                        "TIVTVRPVRTDVCEAGFFNMLARKAEECVQCFCFGVAKACDSANLFTYAIHPPILSHRVVSV"
                        "ELSPLRQIVINEAAPGQDLLTLLHGVQFRATNVHFSGRETPYLALPADYMGNQLKSYGGNLR"
                        "YEVNYRGSGRPVNGPDVIITGNRFTLTYRVRTQPGQNNRVSIPFVPGGWQKPDGRKASREEI"
                        "MMILANVDNILIRLGYLDSTAREVDLINIALDSAGTADKGLGSASLVEKCQCPPGYVGDSCE"
                        "SCASGYVRQPGGPWLGHCVPFIPDSCPSGTYGDPRRGVPCKECPCPLTGSNNFASGCQQSPD"
                        "GDVVCRCNEGYTGRRCEQCAAGYQGNPLAAGGICRRIPDTSCNVDGTYSVHSNGTCQCKDSV"
                        "IGEQCDTCKSKSFHLNSFTYTGCIECFCSGVGLDCDSSTWYRDQVTSTFGRSRVDHGFVLVT"
                        "NYMQPTPDTVPVSMAAEPNALSFIGSADQSGNTLYWSLPAAFLGNKLSSYGGKLTYTLSYSP"
                        "LPNGIMSRNSAPDVVIKSGEDLRLIHYRKSQVVPSVANTYSVEIKESAWQRGDEVVANREHV"
                        "LMALSDITAIYIKATYTTSTKEASLRQVTLDVATPTNLGTPRAVEVEQCRCPEGYLGLSCEQ"
                        "CAPGYARDPEGGIYLGLCRPCECNGHSKYCNSDTGDCEECSDNTEGPSCERCAAGYVGDATR"
                        "GTIYDCQPDEGYPIPSPPAPGNQTLECTAYCQIEGIYDCRGNECLCKRNVIGDQCDQCRPGT"
                        "YGLSAQNQDGCKECYCSGLASQCRSAALYRQLIPVDFILNAPLITDESGAVQDTENLIPDIS"
                        "RNMYTYTHTSYLPKYWSLRGSVLGNQLFSYGGRLSYSLIVESYGNYERGHDIVLIGNGLKLI"
                        "WSRPDGNENQEEYNVRLHEDEQWTRQDRESARPASRSDFMTVLSDLQHILIRATPRVPTQST"
                        "SIGNVILESAVTTRTPGATHASDIELCQCPSGYVGTSCESCAPLHYRDASGSCSLCPCDVSN"
                        'TESCDLVSGGYVECRCKARWKGDRCREIGE"',
                    ),
                ),
            ),
            ("gene", "complement(70720..75241)", (("/gene=", '"EG:BACR25B3.2"'),)),
            (
                "CDS",
                "complement(join("
                "70720..70988,71424..71621,72605..72768,72839..73016,73086..73559,"
                "75217..75241"
                "))",
                (
                    ("/gene=", '"EG:BACR25B3.2"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genefinder'', version:''084'',"
                        " score:''41.82''); /prediction=(method:''genscan'',"
                        " version:''1.0'')\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72287.1"'),
                    ("/db_xref=", '"GI:6946672"'),
                    (
                        "/translation=",
                        '"MANSKVVAHDESLQGINDSEWQLMGDDIDDGLLDDVDETLKPMETKSEEEDLPTGNWFSQS'
                        "VHRVRRSINRLFGSDDNQERGRRQQRERSQRNRDAINRQKELRRRQKEDHNRWKQMRMERQL"
                        "EKQRLVKRTNHVVFNRATDPRKRASDLYDENEASGYHEEDTTLYRTYFVVNEPYDNEYRDRE"
                        "SVQFQNLQKLLDDDLRNFFHSNYEGNDDEEQEIRSTLERVEPTNDNFKIRVQLRIELPTSVN"
                        "DFGSKLQQQLNVYNRIENLSAATDGVFSFTESSDIEEEAIDVTLPQEEVEGSGSDDSSCRGD"
                        "ATFTCPRSGKTICDEMRCDREIQCPDGEDEEYCNYPNVCTEDQFKCDDKCLELKKRCDGSID"
                        "CLDQTDEAGCINAPEPEPEPEPEPEPEPESEPEAEPEPEPEPEPESEPEQEPEPQVPEANGK"
                        'FY"',
                    ),
                ),
            ),
            ("gene", "121867..127124", (("/gene=", '"EG:BACR25B3.3"'),)),
            (
                "CDS",
                "join("
                "121867..122046,122174..122630,123672..123823,124063..124320,"
                "124392..124688,124755..125018,125094..125254,125317..125576,"
                "126793..127124)",
                (
                    ("/gene=", '"EG:BACR25B3.3"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genscan'', version:''1.0'',"
                        " score:''174.91''); /prediction=(method:''genefinder'',"
                        " version:''084''); /match=(desc:''PROBABLE G PROTEIN-COUPLED"
                        " RECEPTOR C13B9.4 IN CHROMOSOME III'',"
                        " species:''Caenorhabditis elegans'',"
                        " ranges:(query:123671..123775,"
                        " target:SWISS-PROT::Q09460:107..141, score:''80.00''),"
                        " (query:123743..123829, target:SWISS-PROT::Q09460:235..263,"
                        " score:''72.00''), (query:124072..124332,"
                        " target:SWISS-PROT::Q09460:265..351, score:''161.00''),"
                        " (query:124392..124691, target:SWISS-PROT::Q09460:349..448,"
                        " score:''206.00''), (query:124755..124958,"
                        " target:SWISS-PROT::Q09460:448..515, score:''123.00''),"
                        " (query:124764..125027, target:SWISS-PROT::Q09460:454..541,"
                        " score:''108.00'')), method:''blastx'', version:''1.4.9'');"
                        " /match=(desc:''CALCITONIN RECEPTOR PRECURSOR (CT-R)'',"
                        " species:''Sus scrofa (Pig)'', ranges:(query:124165..124236,"
                        " target:SWISS-PROT::P25117:191..214, score:''54.00''),"
                        " (query:124392..124580, target:SWISS-PROT::P25117:233..295,"
                        " score:''118.00''), (query:124725..124886,"
                        " target:SWISS-PROT::P25117:318..371, score:''127.00'')),"
                        " method:''blastx'', version:''1.4.9'')\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72288.1"'),
                    ("/db_xref=", '"GI:6946673"'),
                    (
                        "/translation=",
                        '"MGAGNRKSETKTKTEAEIEIEMERDQFSIAANACMSMGPMLISKDKAPCSGGRVRHADSLH'
                        "IYYAVDGKMTLLSNILDCGGCISAQRFTRLLRQSGSSGPSPSAPTAGTFESKSMLEPTSSHS"
                        "LATGRVPLLHDFDASTTESPGTYVLDGVARVAQLALEPTVMDALPDSDTEQVLGNLNSSAPW"
                        "NLTLASAAATNFENCSALFVNYTLPQTEFAIRKCELDGRWGSRPNATEVNPPGWTDYGPCYK"
                        "PEIIRLMQQMGSKDFDAYIDIARRTRTLEIVGLCLSLFALIVSLLIFCTFRSLRNNRTKIHK"
                        "NLFVAMVLQVIIRLTLYLDQFRRGNKEAATNTSLSVIENTPYLCEASYVLLEYARTAMFMWM"
                        "FIEGLYLHNMVTVAVFQGSFPLKFFSRLGWCVPILMTTVWARCTVMYMDTSLGECLWNYNLT"
                        "PYYWILEGPRLAVILLNFCFLVNIIRVLVMKLRQSQASDIEQTRKAVRAAIVLLPLLGITNL"
                        "LHQLAPLKTATNFAVWSYGTHFLTSFQGFFIALIYCFLNGEVRAVLLKSLATQLSVRGHPEW"
                        "APKRASMYSGAYNTAPDTDAVQPAGDPSATGKRISPPNKRLNGRKPSSASIVMIHEPQQRQR"
                        "LMPRLQNKAREKGKDRVEKTDAEAEPDPTISHIHSKEAGSARSRTRGSKWIMGICFRGQMCD"
                        "AGLAKDAANIHDVANAADVDACSGSNNNYHNINNNNGSQNNNSIHCNHRDDDKVKGESQSDF"
                        'KEPSNTNAESLVHLALFTAHTSNTQNNTHRNTIFTPIRRRNCS"',
                    ),
                ),
            ),
            ("gene", "complement(128489..129414)", (("/gene=", '"EG:BACR25B3.4"'),)),
            (
                "CDS",
                "complement(join("
                "128489..128715,128777..129140,129196..129313,129374..129414"
                "))",
                (
                    ("/gene=", '"EG:BACR25B3.4"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genefinder'', version:''084'',"
                        " score:''61.35''); /prediction=(method:''genscan'',"
                        " version:''1.0''); /match=(desc:''VACUOLAR PROTON-ATPASE"
                        " SUBUNIT D'', species:''Oryctolagus cuniculus (Rabbit)'',"
                        " ranges:(query:129190..129324, target:SPTREMBL::O97755:55..11,"
                        " score:''130.00''), (query:128778..129176,"
                        " target:SPTREMBL::O97755:174..42, score:''472.00''),"
                        " (query:128546..128716, target:SPTREMBL::O97755:231..175,"
                        " score:''169.00'')), method:''blastx'', version:''1.4.9'');"
                        " /match=(desc:''VACUOLAR ATP SYNTHASE SUBUNIT D (EC 3.6.1.34)"
                        " (V-ATPASE D SUBUNIT) (V- ATPASE 28 KD ACCESSORY PROTEIN)'',"
                        " species:''Bos taurus (Bovine)'',"
                        " ranges:(query:129190..129324,"
                        " target:SWISS-PROT::P39942:55..11, score:''130.00''),"
                        " (query:128778..129176, target:SWISS-PROT::P39942:174..42,"
                        " score:''471.00''), (query:128546..128716,"
                        " target:SWISS-PROT::P39942:231..175, score:''173.00'')),"
                        " method:''blastx'', version:''1.4.9'');"
                        " /match=(desc:''GH28048.5prime GH Drosophila melanogaster head"
                        " pOT2 Drosophila melanogaster cDNA clone GH28048 5prime, mRNA"
                        " sequence'', species:''Drosophila melanogaster (fruit fly)'',"
                        " ranges:(query:129196..129317, target:EMBL::AI517334:233..112,"
                        " score:''412.00''), (query:128777..129145,"
                        " target:EMBL::AI517334:597..229, score:''1251.00'')),"
                        " method:''blastn'', version:''1.4.9'');"
                        " /match=(desc:''GH07112.5prime GH Drosophila melanogaster head"
                        " pOT2 Drosophila melanogaster cDNA clone GH07112 5prime, mRNA"
                        " sequence'', species:''Drosophila melanogaster (fruit fly)'',"
                        " ranges:(query:129196..129317, target:EMBL::AI108302:223..102,"
                        " score:''412.00''), (query:128777..129145,"
                        " target:EMBL::AI108302:587..219, score:''1251.00''),"
                        " (query:128636..128716, target:EMBL::AI108302:667..587,"
                        " score:''243.00'')), method:''blastn'', version:''1.4.9'')\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72289.1"'),
                    ("/db_xref=", '"GI:6946674"'),
                    (
                        "/translation=",
                        '"MAAKDRLPIFPSRGAQTLMKSRLAGATKGHGLLKKKADALQMRFRLILGKIIETKTLMGQV'
                        "MKEAAFSLAEVKFTTGDINQIVLQNVTKAQIKIRTKKDNVAGVTLPIFEPYTDGVDTYELAG"
                        "LARGGQQLAKLKKNYQSAVRLLVQLASLQTSFVTLDDVIKVTNRRVNAIEHVIIPRINRTIE"
                        "YIISELDELEREEFYRLKKIQDKKREARKASDKLRAEQRLLGQMAEAQEVQNILDEDGDEDL"
                        'LF"',
                    ),
                ),
            ),
            ("gene", "132240..132926", (("/gene=", '"EG:BACR25B3.5"'),)),
            (
                "CDS",
                "132240..132926",
                (
                    ("/gene=", '"EG:BACR25B3.5"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genefinder'', version:''084'',"
                        " score:''48.06''); /prediction=(method:''genscan'',"
                        " version:''1.0'', score:''132.90'');"
                        " /match=(desc:''N-ACETYLTRANSFERASE'', species:''Drosophila"
                        " melanogaster (Fruit fly)'', ranges:(query:132249..132326,"
                        " target:SPTREMBL::Q94521:60..85, score:''64.00''),"
                        " (query:132600..132842, target:SPTREMBL::Q94521:171..251,"
                        " score:''105.00'')), method:''blastx'', version:''1.4.9'');"
                        " EST embl|AI063093|AI063093 comes from the 3' UTR\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72290.1"'),
                    ("/db_xref=", '"GI:6946675"'),
                    (
                        "/translation=",
                        '"MEYKMIAPEHSEQVMEHLRRNFFADEPLNKAAGLCQNGSSCPALEAHCAEAIQHRMSVMAV'
                        "DAKEKDTLKIVGVVLNGILKPGDTAKALSKLDCNDDADFRKIFDLLHRHNLKHNLFEHFDVD"
                        "CMFDVRILSVDSCYRGQGIANELVKRSVAVAKKNGFRLLKADATGIFSQKIFRSHGFEVFSE"
                        'QPYSKYTDENGKVILPVEAPHIKLQQLYKAICADDQDEKKQSL"',
                    ),
                ),
            ),
            ("gene", "complement(133492..134407)", (("/gene=", '"EG:BACR25B3.6"'),)),
            (
                "CDS",
                "complement(join("
                "133492..133595,133663..133748,133867..134135,134198..134407"
                "))",
                (
                    ("/gene=", '"EG:BACR25B3.6"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genscan'', version:''1.0'',"
                        " score:''119.22''); /prediction=(method:''genefinder'',"
                        " version:''084''); /match=(desc:''LD41675.5prime LD Drosophila"
                        " melanogaster embryo pOT2 Drosophila melanogaster cDNA clone"
                        " LD41675 5prime, mRNA sequence'', species:''Drosophila"
                        " melanogaster (fruit fly)'', ranges:(query:134192..134531,"
                        " target:EMBL::AI515958:340..1, score:''1691.00''),"
                        " (query:133879..134139, target:EMBL::AI515958:591..331,"
                        " score:''1305.00'')), method:''blastn'', version:''1.4.9'')\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72291.1"'),
                    ("/db_xref=", '"GI:6946676"'),
                    (
                        "/translation=",
                        '"MNGLPPSKHYNLTHYQQRYNWDCGLSCIIMILSAQQREQLLGNFDAVCGEEGFGSSTWTID'
                        "LCYLLMRYQVRHEYFTQTLGIDPNYAQHTYYSKIIDKDERRVTRKFKDARAHGLRVEQRTVD"
                        "MEVILRHLARHGPVILLTNASLLTCEVCKRNVLEKFGYAGHYVVLCGYDMAAQKLFYHNPEV"
                        'HDGHICRCLIESMDTARRAYGTDEDIIFIYEKKETRE"',
                    ),
                ),
            ),
            ("gene", "135479..136829", (("/gene=", '"EG:BACR25B3.7"'),)),
            (
                "CDS",
                "join(135479..135749,135961..136586,136641..136829)",
                (
                    ("/gene=", '"EG:BACR25B3.7"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genefinder'', version:''084'',"
                        " score:''66.07''); /prediction=(method:''genscan'',"
                        " version:''1.0'', score:''145.64'');"
                        " /match=(desc:''HYPOTHETICAL 40.4 KD TRP-ASP REPEATS"
                        " CONTAINING PROTEIN C14B1.4 IN CHROMOSOME III'',"
                        " species:''Caenorhabditis elegans'',"
                        " ranges:(query:135548..135748,"
                        " target:SWISS-PROT::Q17963:39..105, score:''120.00''),"
                        " (query:135957..136586, target:SWISS-PROT::Q17963:105..314,"
                        " score:''899.00''), (query:136641..136823,"
                        " target:SWISS-PROT::Q17963:315..375, score:''219.00'')),"
                        " method:''blastx'', version:''1.4.9'');"
                        " /match=(desc:''LD30385.5prime LD Drosophila melanogaster"
                        " embryo pOT2 Drosophila melanogaster cDNA clone LD30385"
                        " 5prime, mRNA sequence'', species:''Drosophila melanogaster"
                        " (fruit fly)'', ranges:(query:135288..135749,"
                        " target:EMBL::AA950546:102..563, score:''2301.00''),"
                        " (query:135956..136047, target:EMBL::AA950546:559..650,"
                        " score:''442.00'')), method:''blastn'', version:''1.4.9'');"
                        " /match=(desc:''LD10938.5prime LD Drosophila melanogaster"
                        " embryo BlueScript Drosophila melanogaster cDNA clone LD10938"
                        " 5prime, mRNA sequence'', species:''Drosophila melanogaster"
                        " (fruit fly)'', ranges:(query:136108..136288,"
                        " target:EMBL::AA392005:776..596, score:''212.00'')),"
                        " method:''blastn'', version:''1.4.9'')\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72292.1"'),
                    ("/db_xref=", '"GI:6946677"'),
                    (
                        "/translation=",
                        '"MVPIGAVHGGHPGVVHPPQQPLPTAPSGPNSLQPNSVGQPGATTSSNSSASNKSSLSVKPN'
                        "YTLKFTLAGHTKAVSAVKFSPNGEWLASSSADKLIKIWGAYDGKFEKTISGHKLGISDVAWS"
                        "SDSRLLVSGSDDKTLKVWELSTGKSLKTLKGHSNYVFCCNFNPQSNLIVSGSFDESVRIWDV"
                        "RTGKCLKTLPAHSDPVSAVHFNRDGSLIVSSSYDGLCRIWDTASGQCLKTLIDDDNPPVSFV"
                        "KFSPNGKYILAATLDNTLKLWDYSKGKCLKTYTGHKNEKYCIFANFSVTGGKWIVSGSEDNM"
                        'VYIWNLQSKEVVQKLQGHTDTVLCTACHPTENIIASAALENDKTIKLWKSDT"',
                    ),
                ),
            ),
            ("gene", "145403..147087", (("/gene=", '"EG:BACR25B3.8"'),)),
            (
                "CDS",
                "join(145403..146203,146515..147087)",
                (
                    ("/gene=", '"EG:BACR25B3.8"'),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72293.1"'),
                    ("/db_xref=", '"GI:6946678"'),
                    (
                        "/translation=",
                        '"MNSTTKHLLHCTLLITVIVTFEVFSGGIKIDENSFTLVDPWTEYGQLATVLLYLLRFLTLL'
                        "TLPQVLFNFCGLVFYNAFPEKVVLKGSPLLAPFICIRVVTRGDFPDLVKTNVLRNMNTCLDT"
                        "GLENFLIEVVTDKAVNLSQHRRIREIVVPKEYKTRTGALFKSRALQYCLEDNVNVLNDSDWI"
                        "VHLDEETLLTENSVRGIINFVLDGKHPFGQGLITYANENVVNWLTTLADSFRVSDDMGKLRL"
                        "QFKLFHKPLFSWKGSYVVTQVSAERSVSFDNGIDGSVAEDCFFAMRAFSQGYTFNFIEGEMY"
                        "EKSPFTLLDFLQQRKRWLQGILLVVHSKMIPFKHKLLLGISVYSWVTMPLSTSNIIFAALYP"
                        "IPCPNLVDFVCAFIAAINIYMYVFGVIKSFSLYRFGLFRFLACVLGAVCTIPVNVVIENVAV"
                        'IWGLVGKKHKFYVVQKDVRVLETV"',
                    ),
                ),
            ),
            ("gene", "complement(148860..152785)", (("/gene=", '"EG:BACR25B3.9"'),)),
            (
                "CDS",
                "complement(join("
                "148860..148905,148966..149462,149546..151809,151881..152032,"
                "152106..152785"
                "))",
                (
                    ("/gene=", '"EG:BACR25B3.9"'),
                    (
                        "/note=",
                        "\"/prediction=(method:''genscan'', version:''1.0'');"
                        " /prediction=(method:''genefinder'', version:''084'');"
                        " /match=(desc:''HYPOTHETICAL 135.8 KD PROTEIN'',"
                        " species:''Drosophila melanogaster (Fruit fly)'',"
                        " ranges:(query:152096..152785, target:SPTREMBL::Q9XZ29:230..1,"
                        " score:''1147.00''), (query:151882..152043,"
                        " target:SPTREMBL::Q9XZ29:277..224, score:''250.00''),"
                        " (query:149546..151816, target:SPTREMBL::Q9XZ29:1032..276,"
                        " score:''3735.00''), (query:148953..149465,"
                        " target:SPTREMBL::Q9XZ29:1202..1032, score:''890.00''),"
                        " (query:148863..148907, target:SPTREMBL::Q9XZ29:1212..1198,"
                        " score:''76.00'')), method:''blastx'', version:''1.4.9'');"
                        " /match=(desc:''LD21815.5prime LD Drosophila melanogaster"
                        " embryo pOT2 Drosophila melanogaster cDNA clone LD21815 5prime"
                        " similar to L19117: Drosophila melanogaster (chromosome X"
                        " 3A6-8) kinesin-like protein of 3A (klp3A) mRNA sequence'',"
                        " species:''Drosophila melanogaster (fruit fly)'',"
                        " ranges:(query:152482..152787, target:EMBL::AA816942:460..155,"
                        " score:''1485.00''), (query:152401..152483,"
                        " target:EMBL::AA816942:540..458, score:''397.00'')),"
                        " method:''blastn'', version:''1.4.9'')\"",
                    ),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72294.1"'),
                    ("/db_xref=", '"GI:6946679"'),
                    (
                        "/translation=",
                        '"MSSEDPSCVAVALRVRPLVQSELDRGCRIAVERSADGAPQVTVNRNESYTYNYVFDIDDSQ'
                        "KDLFETCVQAKVKKLLNGYNVTILAYGQTGSGKTYTMGTAFNGVLDDHVGVIPRAVHDIFTA"
                        "IAEMQSEFRFAVTCSFVELYQEQFYDLFSSKTRDKATVDIREVKNRIIMPGLTELVVTSAQQ"
                        "VTDHLIRGSAGRAVAATAMNETSSRSHAIFTLTLVATKLDGKQSVTTSRFNLVDLAGSERCS"
                        "KTLASGDRFKEGVNINKGLLALGNVINALGSGQAAGYIPYRQSKLTRLLQDSLGGNSITLMI"
                        "ACVSPADYNVAETLSTLRYADRALQIKNKPVVNLDPHAAEVNMLKDVIQKLRVELLSGGKMS"
                        "SSLISAVGAAGLGAIPCEESLAGSMANAAEIQRLKEQVRTLQDRNRKLQQELHQSLLDLTEK"
                        "EMRAHIAEQAHDKLRSHVSELKNKLDQREQAQFGNENTNGDNEMRDFSLLVNRVHVELQRTQ"
                        "EELESQGHESRQRLSSRSHTEGGESGGDEVHEMLHSHSEEYTNKQMNFAGELRNINRQLDLK"
                        "QELHERIMRNFSRLDSDDEDVKLRLCNQKIDDLEAERRDLMDQLRNIKSKDISAKLAEERRK"
                        "RLQLLEQEISDLRRKLITQANLLKIRDKEREKIQNLSTEIRTMKESKVKLIRAMRGESEKFR"
                        "QWKMVREKELTQLKSKDRKMQSEIVRQQTLHSKQRQVLKRKCEEALAANKRLKDALERQASA"
                        "QAQRHKYKDNGGSAAGSSNANAKTDSWVDRELEIILSLIDAEHSLEQLMEDRAVINNHYHLL"
                        "QQEKTSDPAEAAEQARILASLEEELEMRNAQISDLQQKVCPTDLDSRIRSLAEGVQSLGESR"
                        "TVSKQLLKTLVQQRRLQASSLNEQRTTLDELRAQLLDAQQQEDAASKRLRLLQSQHEEQMLA"
                        "QQRAYEEKVSVLIRTANQRWAEARSPAEDQQRNQILEELLSSREALQQELDKLRAKNKSKSK"
                        "AVKSEPQDLDDSFQIVDGNETVVLSDVSDDPDWVPSTSKSKRIQSDSRNVISPPEKQDANVT"
                        "SLGNSSIQSLNSTSATEDGKRCKGCKCRTKCTTKRCGCLSGNNACSETCVCKSNCRNPLNLK"
                        "DHASQCGDGDGQKDETEDADKSDDDGDDEPQTSKENAVKFVTPEAPGKVVASPKQTLQEPKA"
                        'AATPLMNSNVVEDINGPKLAKMSGLAFDTPKRKFF"',
                    ),
                ),
            ),
            (
                "gene",
                "join(153490..154269,AL121804.2:41..610,AL121804.2:672..1487)",
                (("/gene=", '"EG:BACR7C10.3"'),),
            ),
            (
                "CDS",
                "join(153490..154269,AL121804.2:41..610,AL121804.2:672..1487)",
                (
                    ("/gene=", '"EG:BACR7C10.3"'),
                    ("/codon_start=", "1"),
                    ("/protein_id=", '"CAB72295.1"'),
                    ("/db_xref=", '"GI:6946680"'),
                    (
                        "/translation=",
                        '"MEEEAPRFNVLEEAFNGNGNGCANVEATQSAILKVLTRVNRFQMRVRKHIEDNYTEFLPNN'
                        "TSPDIFLEESGSLNREIHDMLENLGSEGLDALDEANVKMAGNGRQLREILLGLGVSEHVLRI"
                        "DELFQCVEEAKATKDYLVLLDLVGRLRAFIYGDDSVDGDAQVATPEVRRIFKALECYETIKV"
                        "KYHVQAYMLQQSLQERFDRLVQLQCKSFPTSRCVTLQVSRDQTQLQDIVQALFQEPYNPARL"
                        "AEFLLDNCIEPVIMRPVMADYSEEADGGTYVRLSLSYATKEPSSAHVRPNYKQVLENLRLLL"
                        "HTLAGINCSVSRDQHVFGIIGDHVKDKMLKLLVDECLIPAVPESTEEYQTSTLCEDVAQLEQ"
                        "LLVDSFIINPEQDRALGQFVEKYETYYRNRMYRRVLETAREIIQRDLQDMVLVAPNNHSAEV"
                        "ANDPFLFPRCMISKSAQDFVKLMDRILRQPTDKLGDQEADPIAGVISIMLHTYINEVPKVHR"
                        "KLLESIPQQAVLFHNNCMFFTHWVAQHANKGIESLAALAKTLQATGQQHFRVQVDYQSSILM"
                        "GIMQEFEFESTHTLGSGPLKLVRQCLRQLELLKNVWANVLPETVYNATFCELINTFVAELIR"
                        "RVFTLRHISAQMACELSDLIDVVLQRAPTLFREPNEVVQVLSWLKLQQLKAMLNASLMEITE"
                        'LWGDGVGPLTASYKSDEIKHLIRALFQDTDWRAKAITQIV"',
                    ),
                ),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_08(self):
        path = "GenBank/one_of.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 2509
        locus = "HSTMPO1"
        definition = "Human thymopoietin (TMPO) gene, exon 1"
        accession = ["U18266"]
        titles = (
            "Structure and mapping of the human thymopoietin (TMPO) gene and"
            " relationship of TMPO beta to rat lamin-associated polypeptide 2",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..2509",
                (
                    ("/organism=", '"Homo sapiens"'),
                    ("/db_xref=", '"taxon:9606"'),
                    ("/chromosome=", '"12"'),
                    ("/map=", '"12q22; 64% (% distance from centromere to telomere)"'),
                    ("/clone=", '"P1.516 (DMPC-HFFno.1B-0943F)"'),
                    (
                        "/clone_lib=",
                        '"DuPont Merck Hum Fibroblast P1 Library no.1 Series B'
                        ' (compressed) (Genome Systems Inc)"',
                    ),
                ),
            ),
            ("5'UTR", "one-of(1888,1901)..2200", (("/gene=", '"TMPO"'),)),
            (
                "gene",
                "join("
                "1888..2509,U18267.1:1..270,U18268.1:1..309,U18270.1:1..6905,"
                "U18269.1:1..128,U18271.1:1..3234"
                ")",
                (("/gene=", '"TMPO"'),),
            ),
            (
                "exon",
                "one-of(1888,1901)..2479",
                (("/gene=", '"TMPO"'), ("/number=", "1")),
            ),
            (
                "CDS",
                "join("
                "2201..2479,U18267.1:120..246,U18268.1:130..288,U18270.1:4691..4788,"
                "U18269.1:82..>128"
                ")",
                (
                    ("/gene=", '"TMPO"'),
                    ("/codon_start=", "1"),
                    ("/product=", '"thymopoietin beta"'),
                    ("/protein_id=", '"AAB60434.1"'),
                    ("/db_xref=", '"GI:885684"'),
                    (
                        "/translation=",
                        '"MPEFLEDPSVLTKDKLKSELVANNVTLPAGEQRKDVYVQLYLQHLTARNRPPLPAGTNSKG'
                        "PPDFSSDEEREPTPVLGSGAAAAGRSRAAVGRKATKKTDKPRQEDKDDLDVTELTNEDLLDQ"
                        "LVKYGVNPGPIVGTTRKLYEKKLLKLREQGTESRSSTPLPTISSSAENTRQNGSNDSDRYSD"
                        'NEEDSKIELKLEKREPLKGRAKTPVTLKQRRVEHNQSYSQAGITETEWTSGS"',
                    ),
                ),
            ),
            (
                "CDS",
                "join("
                "2201..2479,U18267.1:120..246,U18268.1:130..288,U18270.1:39..1558"
                ")",
                (
                    ("/gene=", '"TMPO"'),
                    ("/codon_start=", "1"),
                    ("/product=", '"thymopoietin alpha"'),
                    ("/protein_id=", '"AAB60433.1"'),
                    ("/db_xref=", '"GI:885683"'),
                    (
                        "/translation=",
                        '"MPEFLEDPSVLTKDKLKSELVANNVTLPAGEQRKDVYVQLYLQHLTARNRPPLPAGTNSKG'
                        "PPDFSSDEEREPTPVLGSGAAAAGRSRAAVGRKATKKTDKPRQEDKDDLDVTELTNEDLLDQ"
                        "LVKYGVNPGPIVGTTRKLYEKKLLKLREQGTESRSSTPLPTISSSAENTRQNGSNDSDRYSD"
                        "NEEGKKKEHKKVKSTRDIVPFSELGTTPSGGGFFQGISFPEISTRPPLGSTELQAAKKVHTS"
                        "KGDLPREPLVATNLPGRGQLQKLASERNLFISCKSSHDRCLEKSSSSSSQPEHSAMLVSTAA"
                        "SPSLIKETTTGYYKDIVENICGREKSGIQPLCPERSHISDQSPLSSKRKALEESESSQLISP"
                        "PLAQAIRDYVNSLLVQGGVGSLPGTSNSMPPLDVENIQKRIDQSKFQETEFLSPPRKVPRLS"
                        "EKSVEERDSGSFVAFQNIPGSELMSSFAKTVVSHSLTTLGLEVAKQSQHDKIDASELSFPFH"
                        "ESILKVIEEEWQQVDRQLPSLACKYPVSSREATQILSVPKVDDEILGFISEATPLGGIQAAS"
                        "TESCNQQLDLALCRAYEAAASALQIATHTAFVAKAMQADISEAAQILSSDPSRTHQALGILS"
                        "KTYDAASYICEAAFDEVKMAAHTMGNATVGRRYLWLKDCKINLASKNKLASTPFKGGTLFGG"
                        'EVCKVIKKRGNKH"',
                    ),
                ),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_09(self):
        path = "GenBank/NT_019265.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 0
        locus = "NT_019265"
        definition = "Homo sapiens chromosome 1 working draft sequence segment"
        accession = ["NT_019265"]
        titles = ("Direct Submission",)
        features = [
            (
                "source",
                "1..1250660",
                (
                    ("/organism=", '"Homo sapiens"'),
                    ("/db_xref=", '"taxon:9606"'),
                    ("/chromosome=", '"1"'),
                ),
            ),
            (
                "source",
                "1..3290",
                (
                    ("/note=", '"Accession AL391218 sequenced by The Sanger Centre"'),
                    ("/organism=", '"Homo sapiens"'),
                    ("/db_xref=", '"taxon:9606"'),
                    ("/clone=", '"RP11-13G5"'),
                ),
            ),
            (
                "misc_feature",
                "215902..365470",
                (
                    ("/standard_name=", '"RP11-242F24"'),
                    ("/note=", '"FISH-mapped clone"'),
                ),
            ),
            (
                "variation",
                "217508",
                (
                    ("/allele=", '"T"'),
                    ("/allele=", '"C"'),
                    ("/db_xref=", '"dbSNP:811400"'),
                ),
            ),
            (
                "mRNA",
                "join("
                "342430..342515,363171..363300,365741..365814,376398..376499,"
                "390169..390297,391257..391379,392606..392679,398230..398419,"
                "399082..399167,399534..399650,405844..405913,406704..406761,"
                "406868..407010,407962..408091,408508..409092"
                ")",
                (
                    ("/gene=", '"FLJ10737"'),
                    ("/product=", '"hypothetical protein FLJ10737"'),
                    ("/transcript_id=", '"XM_057697.1"'),
                    ("/db_xref=", '"LocusID:55735"'),
                ),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_10(self):
        path = "GenBank/origin_line.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 180
        locus = "NC_002678"
        definition = "Mesorhizobium loti, complete genome (edited)"
        accession = ["NC_002678"]
        titles = (
            "Complete genome structure of the nitrogen-fixing symbiotic bacterium"
            " Mesorhizobium loti",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..180",
                (
                    ("/organism=", '"Mesorhizobium loti"'),
                    ("/strain=", '"MAFF303099"'),
                    ("/db_xref=", '"taxon:381"'),
                ),
            ),
            ("gene", "20..120", (("/gene=", '"fake"'),)),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_11(self):
        path = "GenBank/blank_seq.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 360
        locus = "NP_001832"
        definition = "cannabinoid receptor 2 (macrophage) [Homo sapiens]"
        accession = ["NP_001832"]
        titles = (
            "Molecular characterization of a peripheral receptor for cannabinoids",
            "Expression of central and peripheral cannabinoid receptors in human "
            "immune tissues and leukocyte subpopulations",
            "Molecular cloning, expression and function of the murine CB2 peripheral "
            "cannabinoid receptor",
        )
        features = [
            (
                "source",
                "1..360",
                (
                    ("/organism=", '"Homo sapiens"'),
                    ("/db_xref=", '"taxon:9606"'),
                    ("/chromosome=", '"1"'),
                    ("/map=", '"1p36.11"'),
                ),
            ),
            (
                "Protein",
                "1..360",
                (("/product=", '"cannabinoid receptor 2 (macrophage)"'),),
            ),
            (
                "Region",
                "50..299",
                (
                    ("/region_name=", '"7 transmembrane receptor (rhodopsin family)"'),
                    ("/db_xref=", '"CDD:pfam00001"'),
                    ("/note=", '"7tm_1"'),
                ),
            ),
            (
                "CDS",
                "1..360",
                (
                    ("/pseudo", ""),
                    ("/gene=", '"CNR2"'),
                    ("/db_xref=", '"LocusID:1269"'),
                    ("/db_xref=", '"MIM:605051"'),
                    ("/coded_by=", '"NM_001841.1:127..1209"'),
                ),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_12(self):
        path = "GenBank/dbsource_wrap.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 64
        locus = "SCX3_BUTOC"
        definition = "Neurotoxin III"
        accession = ["P01485"]
        titles = (
            "Neurotoxins from the venoms of two scorpions: Buthus occitanus tunetanus"
            " and Buthus occitanus mardochei",
        )
        features = [
            (
                "source",
                "1..64",
                (
                    ("/organism=", '"Buthus occitanus tunetanus"'),
                    ("/db_xref=", '"taxon:6871"'),
                ),
            ),
            ("Protein", "1..64", (("/product=", '"Neurotoxin III"'),)),
            (
                "Bond",
                "bond(12,63)",
                (("/bond_type=", '"disulfide"'), ("/note=", '"BY SIMILARITY."')),
            ),
            (
                "Bond",
                "bond(16,36)",
                (("/bond_type=", '"disulfide"'), ("/note=", '"BY SIMILARITY."')),
            ),
            (
                "Bond",
                "bond(22,46)",
                (("/bond_type=", '"disulfide"'), ("/note=", '"BY SIMILARITY."')),
            ),
            (
                "Bond",
                "bond(26,48)",
                (("/bond_type=", '"disulfide"'), ("/note=", '"BY SIMILARITY."')),
            ),
            ("Site", "64", (("/site_type=", '"amidation"'),)),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_13(self):
        path = "GenBank/gbvrl1_start.seq"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
            length = 2007
            locus = "AB000048"
            definition = (
                "Feline panleukopenia virus DNA for nonstructural protein 1,"
                " complete cds"
            )
            accession = ["AB000048"]
            titles = (
                "Evolutionary pattern of feline panleukopenia virus differs from that"
                " of canine parvovirus",
                "Direct Submission",
            )
            features = [
                (
                    "source",
                    "1..2007",
                    (
                        ("/organism=", '"Feline panleukopenia virus"'),
                        ("/mol_type=", '"genomic DNA"'),
                        ("/isolate=", '"483"'),
                        ("/db_xref=", '"taxon:10786"'),
                        ("/lab_host=", '"Felis domesticus"'),
                    ),
                ),
                (
                    "CDS",
                    "1..2007",
                    (
                        ("/codon_start=", "1"),
                        ("/product=", '"nonstructural protein 1"'),
                        ("/protein_id=", '"BAA19009.1"'),
                        ("/db_xref=", '"GI:1769754"'),
                        (
                            "/translation=",
                            '"MSGNQYTEEVMEGVNWLKKHAEDEAFSFVFKCDNVQLNGKDVRWNNYTKPIQNEELT'
                            "SLIRGAQTAMDQTEEEEMDWESEVDSLAKKQVQTFDALIKKCLFEVFVSKNIEPNECV"
                            "WFIQHEWGKDQGWHCHVLLHSKNLQQATGKWLRRQMNMYWSRWLVTLCSINLTPTEKI"
                            "KLREIAEDSEWVTILTYRHKQTKKDYVKMVHFGNMIAYYFLTKKKIVHMTKESGYFLS"
                            "TDSGWKFNFMKYQDRHTVSTLYTEQMKPETVETTVTTAQETKRGRIQTKKEVSIKCTL"
                            "RDLVSKRVTSPEDWMMLQPDSYIEMMAQPGGENLLKNTLEICTLTLARTKTAFELILE"
                            "KADNTKLTNFDLANSRTCQIFRMHGWNWIKVCHAIACVLNRQGGKRNTVLFHGPASTG"
                            "KSIIAQAIAQAVGNVGCYNAANVNFPFNDCTNKNLIWVEEAGNFGQQVNQFKAICSGQ"
                            "TIRIDQKGKGSKQIEPTPVIMTTNENITIVRIGCEERPEHTQPIRDRMLNIKLVCKLP"
                            "GDFGLVDKEEWPLICAWLVKHGYQSTMANYTHHWGKVPEWDENWAEPKIQEGINSPGC"
                            "KDLETQAASNPQSQDHVLTPLTPDVVDLALEPWSTPDTPIAETANQQSNQLGVTHKDV"
                            'QASPTWSEIEADLRAIFTSEQLEEDFRDDLD"',
                        ),
                    ),
                ),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )

            record = next(records)
            length = 2007
            locus = "AB000049"
            definition = (
                "Feline panleukopenia virus DNA for nonstructural protein 1,"
                " complete cds"
            )
            accession = ["AB000049"]
            titles = (
                "Evolutionary pattern of feline panleukopenia virus differs that of"
                " canine parvovirus",
                "Direct Submission",
            )
            features = [
                (
                    "source",
                    "1..2007",
                    (
                        ("/organism=", '"Feline panleukopenia virus"'),
                        ("/mol_type=", '"genomic DNA"'),
                        ("/isolate=", '"94-1"'),
                        ("/db_xref=", '"taxon:10786"'),
                        ("/lab_host=", '"Felis domesticus"'),
                    ),
                ),
                (
                    "CDS",
                    "1..2007",
                    (
                        ("/codon_start=", "1"),
                        ("/product=", '"nonstructural protein 1"'),
                        ("/protein_id=", '"BAA19010.1"'),
                        ("/db_xref=", '"GI:1769756"'),
                        (
                            "/translation=",
                            '"MSGNQYTEEVMEGVNWLKKHAEDEAFSFVFKCDNVQLNGKDVRWNNYTKPIQNEELT'
                            "SLIRGAQTAMDQTEEEEMDWESEVDSLAKKQVQTFDALIKKCLFEVFVSKNIEPNECV"
                            "WFIQHEWGKDQGWHCHVLLHSKNLQQATGKWLRRQMNMYWSRWLVTLCSINLTPTEKI"
                            "KLREIAEDSEWVTILTYRHKQTKKDYVKMVHFGNMIAYYFLTKKKIVHMTKESGYFLS"
                            "TDSGWKFNFMKYQDRHTVSTLYTEQMKPETVETTVTTAQETKRGRIQTKKEVSIKCTL"
                            "RDLVSKRVTSPEDWMMLQPDSYIEMMAQPGGENLLKNTLEICTLTLARTKTAFELILE"
                            "KADNTKLTNFDLANSRTCQIFRMHGWNWIKVCHAIACVLNRQGGKRNTVLFHGPASTG"
                            "KSIIAQAIAQAVGNVGCYNAANVNFPFNDCTNKNLIWVEEAGNFGQQVNQFKAICSGQ"
                            "TIRIDQKGKGSKQIEPTPVIMTTNENITIVRIGCEERPEHTQPIRDRMLNIKLVCKLP"
                            "GDFGLVDKEEWPLICAWLVKHGYQSTMANYTHHWGKVPEWDENWAEPKIQEGINSPGC"
                            "KDLETQAASNPQSQDHVLTPLTPDVVDLALEPWSTPDTPIAETANQQSNQLGVTHKDV"
                            'QASPTWSEIEADLRAIFTSEQLEEDFRDDLD"',
                        ),
                    ),
                ),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )

            record = next(records)
            length = 1755
            locus = "AB000050"
            definition = (
                "Feline panleukopenia virus DNA for capsid protein 2, complete cds"
            )
            accession = ["AB000050"]
            titles = (
                "Evolutionary pattern of feline panleukopenia virus differs from that"
                " of canine parvovirus",
                "Direct Submission",
            )
            features = [
                (
                    "source",
                    "1..1755",
                    (
                        ("/organism=", '"Feline panleukopenia virus"'),
                        ("/mol_type=", '"genomic DNA"'),
                        ("/isolate=", '"94-1"'),
                        ("/db_xref=", '"taxon:10786"'),
                        ("/lab_host=", '"Felis domesticus"'),
                    ),
                ),
                (
                    "CDS",
                    "1..1755",
                    (
                        ("/codon_start=", "1"),
                        ("/product=", '"capsid protein 2"'),
                        ("/protein_id=", '"BAA19011.1"'),
                        ("/db_xref=", '"GI:1769758"'),
                        (
                            "/translation=",
                            '"MSDGAVQPDGGQPAVRNERATGSGNGSGGGGGGGSGGVGISTGTFNNQTEFKFLENG'
                            "WVEITANSSRLVHLNMPESENYKRVVVNNMDKTAVKGNMALDDTHVQIVTPWSLVDAN"
                            "AWGVWFNPGDWQLIVNTMSELHLVSFEQEIFNVVLKTVSESATQPPTKVYNNDLTASL"
                            "MVALDSNNTMPFTPAAMRSETLGFYPWKPTIPTPWRYYFQWDRTLIPSHTGTSGTPTN"
                            "VYHGTDPDDVQFYTIENSVPVHLLRTGDEFATGTFFFDCKPCRLTHTWQTNRALGLPP"
                            "FLNSLPQSEGATNFGDIGVQQDKRRGVTQMGNTDYITEATIMRPAEVGYSAPYYSFEA"
                            "STQGPFKTPIAAGRGGAQTDENQAADGDPRYAFGRQHGQKTTTTGETPERFTYIAHQD"
                            "TGRYPEGDWIQNINFNLPVTNDNVLLPTDPIGGKTGINYTNIFNTYGPLTALNNVPPV"
                            "YPNGQIWDKEFDTDLKPRLHVNAPFVCQNNCPGQLFVKVAPNLTNEYDPDASANMSRI"
                            "VTYSDFWWKGKLVFKAKLRASHTWNPIQQMSINVDNQFNYVPNNIGAMKIVYEKSQLA"
                            'PRKLY"',
                        ),
                    ),
                ),
            ]
            self.perform_record_parser_test(
                record, length, locus, definition, accession, titles, features
            )

    def test_record_parser_14(self):
        path = "GenBank/NC_005816.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Premature end of file in sequence data
                record = next(records)
        length = 9609
        locus = "NC_005816"
        definition = (
            "Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete"
            " sequence"
        )
        accession = ["NC_005816"]
        titles = (
            "Genetics of metabolic variations between Yersinia pestis biovars and the"
            " proposal of a new biovar, microtus",
            "Complete genome sequence of Yersinia pestis strain 91001, an isolate"
            " avirulent to humans",
            "Direct Submission",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..9609",
                (
                    ("/organism=", '"Yersinia pestis biovar Microtus str. 91001"'),
                    ("/mol_type=", '"genomic DNA"'),
                    ("/strain=", '"91001"'),
                    ("/db_xref=", '"taxon:229193"'),
                    ("/plasmid=", '"pPCP1"'),
                    ("/biovar=", '"Microtus"'),
                ),
            ),
            ("repeat_region", "1..1954", ()),
            (
                "gene",
                "87..1109",
                (("/locus_tag=", '"YP_pPCP01"'), ("/db_xref=", '"GeneID:2767718"')),
            ),
            (
                "CDS",
                "87..1109",
                (
                    ("/locus_tag=", '"YP_pPCP01"'),
                    (
                        "/note=",
                        '"similar to corresponding CDS from previously sequenced pPCP'
                        " plasmid of Yersinia pestis KIM (AF053945) and CO92"
                        " (AL109969), also many transposase entries for insertion"
                        " sequence IS100 of Yersinia pestis. Contains IS21-like element"
                        ' transposase, HTH domain (Interpro|IPR007101)"',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"putative transposase"'),
                    ("/protein_id=", '"NP_995567.1"'),
                    ("/db_xref=", '"GI:45478712"'),
                    ("/db_xref=", '"GeneID:2767718"'),
                    (
                        "/translation=",
                        '"MVTFETVMEIKILHKQGMSSRAIARELGISRNTVKRYLQAKSEPPKYTPRPAVASLLDEYR'
                        "DYIRQRIADAHPYKIPATVIAREIRDQGYRGGMTILRAFIRSLSVPQEQEPAVRFETEPGRQ"
                        "MQVDWGTMRNGRSPLHVFVAVLGYSRMLYIEFTDNMRYDTLETCHRNAFRFFGGVPREVLYD"
                        "NMKTVVLQRDAYQTGQHRFHPSLWQFGKEMGFSPRLCRPFRAQTKGKVERMVQYTRNSFYIP"
                        "LMTRLRPMGITVDVETANRHGLRWLHDVANQRKHETIQARPCDRWLEEQQSMLALPPEKKEY"
                        'DVHLDENLVNFDKHPLHHPLSIYDSFCRGVA"',
                    ),
                ),
            ),
            (
                "misc_feature",
                "87..959",
                (
                    ("/locus_tag=", '"YP_pPCP01"'),
                    (
                        "/note=",
                        '"Transposase and inactivated derivatives [DNA replication,'
                        ' recombination, and repair]; Region: COG4584"',
                    ),
                    ("/db_xref=", '"CDD:34222"'),
                ),
            ),
            (
                "misc_feature",
                "<111..209",
                (
                    ("/locus_tag=", '"YP_pPCP01"'),
                    (
                        "/note=",
                        '"Helix-turn-helix domain of Hin and related proteins, a family'
                        " of DNA-binding domains unique to bacteria and represented by"
                        " the Hin protein of Salmonella. The basic HTH domain is a"
                        " simple fold comprised of three core helices that form a"
                        ' right-handed...; Region: HTH_Hin_like; cl01116"',
                    ),
                    ("/db_xref=", '"CDD:186341"'),
                ),
            ),
            (
                "misc_feature",
                "438..812",
                (
                    ("/locus_tag=", '"YP_pPCP01"'),
                    ("/note=", '"Integrase core domain; Region: rve; cl01316"'),
                    ("/db_xref=", '"CDD:194099"'),
                ),
            ),
            (
                "gene",
                "1106..1888",
                (("/locus_tag=", '"YP_pPCP02"'), ("/db_xref=", '"GeneID:2767716"')),
            ),
            (
                "CDS",
                "1106..1888",
                (
                    ("/locus_tag=", '"YP_pPCP02"'),
                    (
                        "/note=",
                        '"similar to corresponding CDS form previously sequenced pPCP'
                        " plasmid of Yersinia pestis KIM (AF053945) and CO92"
                        " (AL109969), also many ATP-binding protein entries for"
                        " insertion sequence IS100 of Yersinia pestis. Contains"
                        " Chaperonin clpA/B (Interpro|IPR001270). Contains"
                        " ATP/GTP-binding site motif A (P-loop) (Interpro|IPR001687,"
                        " Molecular Function: ATP binding (GO:0005524)). Contains"
                        " Bacterial chromosomal replication initiator protein, DnaA"
                        " (Interpro|IPR001957, Molecular Function: DNA binding"
                        " (GO:0003677), Molecular Function: DNA replication origin"
                        " binding (GO:0003688), Molecular Function: ATP binding"
                        " (GO:0005524), Biological Process: DNA replication initiation"
                        " (GO:0006270), Biological Process: regulation of DNA"
                        " replication (GO:0006275)). Contains AAA ATPase"
                        " (Interpro|IPR003593, Molecular Function: nucleotide binding"
                        ' (GO:0000166))"',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"transposase/IS protein"'),
                    ("/protein_id=", '"NP_995568.1"'),
                    ("/db_xref=", '"GI:45478713"'),
                    ("/db_xref=", '"GeneID:2767716"'),
                    (
                        "/translation=",
                        '"MMMELQHQRLMALAGQLQLESLISAAPALSQQAVDQEWSYMDFLEHLLHEEKLARHQRKQA'
                        "MYTRMAAFPAVKTFEEYDFTFATGAPQKQLQSLRSLSFIERNENIVLLGPSGVGKTHLAIAM"
                        "GYEAVRAGIKVRFTTAADLLLQLSTAQRQGRYKTTLQRGVMAPRLLIIDEIGYLPFSQEEAK"
                        "LFFQVIAKRYEKSAMILTSNLPFGQWDQTFAGDAALTSAMLDRILHHSHVVQIKGESYRLRQ"
                        'KRKAGVIAEANPE"',
                    ),
                ),
            ),
            (
                "misc_feature",
                "1109..1885",
                (
                    ("/locus_tag=", '"YP_pPCP02"'),
                    (
                        "/note=",
                        '"transposase/IS protein; Provisional; Region: PRK09183"',
                    ),
                    ("/db_xref=", '"CDD:181681"'),
                ),
            ),
            (
                "misc_feature",
                "1367..>1669",
                (
                    ("/locus_tag=", '"YP_pPCP02"'),
                    (
                        "/note=",
                        '"The AAA+ (ATPases Associated with a wide variety of cellular'
                        " Activities) superfamily represents an ancient group of"
                        " ATPases belonging to the ASCE (for additional strand,"
                        " catalytic E) division of the P-loop NTPase fold. The ASCE"
                        ' division also includes...; Region: AAA; cd00009"',
                    ),
                    ("/db_xref=", '"CDD:99707"'),
                ),
            ),
            (
                "misc_feature",
                "1433..1456",
                (
                    ("/locus_tag=", '"YP_pPCP02"'),
                    ("/note=", '"Walker A motif; other site"'),
                    ("/db_xref=", '"CDD:99707"'),
                ),
            ),
            (
                "misc_feature",
                "order(1436..1459,1619..1621)",
                (
                    ("/locus_tag=", '"YP_pPCP02"'),
                    ("/note=", '"ATP binding site [chemical binding]; other site"'),
                    ("/db_xref=", '"CDD:99707"'),
                ),
            ),
            (
                "misc_feature",
                "1607..1624",
                (
                    ("/locus_tag=", '"YP_pPCP02"'),
                    ("/note=", '"Walker B motif; other site"'),
                    ("/db_xref=", '"CDD:99707"'),
                ),
            ),
            (
                "gene",
                "2925..3119",
                (
                    ("/gene=", '"rop"'),
                    ("/locus_tag=", '"YP_pPCP03"'),
                    ("/gene_synonym=", '"rom"'),
                    ("/db_xref=", '"GeneID:2767717"'),
                ),
            ),
            (
                "CDS",
                "2925..3119",
                (
                    ("/gene=", '"rop"'),
                    ("/locus_tag=", '"YP_pPCP03"'),
                    ("/gene_synonym=", '"rom"'),
                    (
                        "/note=",
                        '"Best Blastp hit =gi|16082682|ref|NP_395229.1| (NC_003132)'
                        " putative replication regulatory protein [Yersinia pestis],"
                        " gi|5763813|emb|CAB531 66.1| (AL109969) putative replication"
                        " regulatory protein [Yersinia pestis]; similar to"
                        " gb|AAK91579.1| (AY048853), RNAI modulator protein Rom"
                        " [Salmonella choleraesuis], Contains Regulatory protein Rop"
                        ' (Interpro|IPR000769)"',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"putative replication regulatory protein"'),
                    ("/protein_id=", '"NP_995569.1"'),
                    ("/db_xref=", '"GI:45478714"'),
                    ("/db_xref=", '"GeneID:2767717"'),
                    (
                        "/translation=",
                        '"MNKQQQTALNMARFIRSQSLILLEKLDALDADEQAAMCERLHELAEELQNSIQARFEAESE'
                        'TGT"',
                    ),
                ),
            ),
            (
                "misc_feature",
                "2925..3107",
                (
                    ("/gene=", '"rop"'),
                    ("/locus_tag=", '"YP_pPCP03"'),
                    ("/gene_synonym=", '"rom"'),
                    ("/note=", '"Rop protein; Region: Rop; pfam01815"'),
                    ("/db_xref=", '"CDD:145136"'),
                ),
            ),
            (
                "gene",
                "3486..3857",
                (("/locus_tag=", '"YP_pPCP04"'), ("/db_xref=", '"GeneID:2767720"')),
            ),
            (
                "CDS",
                "3486..3857",
                (
                    ("/locus_tag=", '"YP_pPCP04"'),
                    (
                        "/note=",
                        '"Best Blastp hit = gi|321919|pir||JQ1541 hypothetical 16.9K'
                        ' protein - Salmonella typhi murium plasmid NTP16."',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"hypothetical protein"'),
                    ("/protein_id=", '"NP_995570.1"'),
                    ("/db_xref=", '"GI:45478715"'),
                    ("/db_xref=", '"GeneID:2767720"'),
                    (
                        "/translation=",
                        '"MSKKRRPQKRPRRRRFFHRLRPPDEHHKNRRSSQRWRNPTGLKDTRRFPPEAPSCALLFRP'
                        'CRLPDTSPPFSLREAWRFLIAHAVGISVRCRSFAPSWAVCTNPPFSPTTAPYPVTIVLSPTR"',
                    ),
                ),
            ),
            (
                "misc_feature",
                "3498..3626",
                (
                    ("/locus_tag=", '"YP_pPCP04"'),
                    (
                        "/note=",
                        '"ProfileScan match to entry PS50323 ARG_RICH, E-value 8.981"',
                    ),
                ),
            ),
            (
                "gene",
                "4343..4780",
                (
                    ("/gene=", '"pim"'),
                    ("/locus_tag=", '"YP_pPCP05"'),
                    ("/db_xref=", '"GeneID:2767712"'),
                ),
            ),
            (
                "CDS",
                "4343..4780",
                (
                    ("/gene=", '"pim"'),
                    ("/locus_tag=", '"YP_pPCP05"'),
                    (
                        "/note=",
                        '"similar to many previously sequenced pesticin immunity'
                        " protein entries of Yersinia pestis plasmid pPCP, e.g. gi|"
                        " 16082683|,ref|NP_395230.1| (NC_003132) ,"
                        " gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655|"
                        " emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1|"
                        ' (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)"',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"pesticin immunity protein"'),
                    ("/protein_id=", '"NP_995571.1"'),
                    ("/db_xref=", '"GI:45478716"'),
                    ("/db_xref=", '"GeneID:2767712"'),
                    (
                        "/translation=",
                        '"MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTE'
                        "FTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVK"
                        'SGNFIVVKEIKKSIPGCTVYYH"',
                    ),
                ),
            ),
            (
                "gene",
                "complement(4815..5888)",
                (
                    ("/gene=", '"pst"'),
                    ("/locus_tag=", '"YP_pPCP06"'),
                    ("/db_xref=", '"GeneID:2767721"'),
                ),
            ),
            (
                "CDS",
                "complement(4815..5888)",
                (
                    ("/gene=", '"pst"'),
                    ("/locus_tag=", '"YP_pPCP06"'),
                    (
                        "/note=",
                        '"Best Blastp hit =|16082684|ref|NP_395231.1| (NC_003132)'
                        " pesticin [Yersinia pestis], gi|984824|gb|AAA75369.1| (U31974)"
                        " pesticin [Yersinia pestis], gi|1488654|emb|CAA63438.1|"
                        " (X92856) pesticin [Yersinia pestis],"
                        " gi|2996220|gb|AAC62544.1| (AF053945) pesticin [Yersinia"
                        " pestis], gi|5763815|emb|CAB53168.1| (AL1099 69) pesticin"
                        ' [Yersinia pestis]"',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"pesticin"'),
                    ("/protein_id=", '"NP_995572.1"'),
                    ("/db_xref=", '"GI:45478717"'),
                    ("/db_xref=", '"GeneID:2767721"'),
                    (
                        "/translation=",
                        '"MSDTMVVNGSGGVPAFLFSGSTLSSYRPNFEANSITIALPHYVDLPGRSNFKLMYIMGFPI'
                        "DTEMEKDSEYSNKIRQESKISKTEGTVSYEQKITVETGQEKDGVKVYRVMVLEGTIAESIEH"
                        "LDKKENEDILNNNRNRIVLADNTVINFDNISQLKEFLRRSVNIVDHDIFSSNGFEGFNPTSH"
                        "FPSNPSSDYFNSTGVTFGSGVDLGQRSKQDLLNDGVPQYIADRLDGYYMLRGKEAYDKVRTA"
                        "PLTLSDNEAHLLSNIYIDKFSHKIEGLFNDANIGLRFSDLPLRTRTALVSIGYQKGFKLSRT"
                        'APTVWNKVIAKDWNGLVNAFNNIVDGMSDRRKREGALVQKDIDSGLLK"',
                    ),
                ),
            ),
            (
                "variation",
                "5910..5911",
                (("/note=", '"compared to AF053945"'), ("/replace=", '""')),
            ),
            (
                "variation",
                "5933^5934",
                (("/note=", '"compared to AL109969"'), ("/replace=", '"a"')),
            ),
            (
                "variation",
                "5933^5934",
                (("/note=", '"compared to AF053945"'), ("/replace=", '"aa"')),
            ),
            (
                "variation",
                "5948",
                (("/note=", '"compared to AL109969"'), ("/replace=", '"c"')),
            ),
            (
                "gene",
                "6005..6421",
                (("/locus_tag=", '"YP_pPCP07"'), ("/db_xref=", '"GeneID:2767719"')),
            ),
            (
                "CDS",
                "6005..6421",
                (
                    ("/locus_tag=", '"YP_pPCP07"'),
                    (
                        "/note=",
                        '"Best Blastp hit = gi|16082685|ref|NP_395232.1| (NC_003132)'
                        " hypothetical protein [Yersinia pestis],"
                        " gi|5763816|emb|CAB53169.1| (AL109969) hypothetical protein"
                        ' [Yersinia pestis]"',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"hypothetical protein"'),
                    ("/protein_id=", '"NP_995573.1"'),
                    ("/db_xref=", '"GI:45478718"'),
                    ("/db_xref=", '"GeneID:2767719"'),
                    (
                        "/translation=",
                        '"MKFHFCDLNHSYKNQEGKIRSRKTAPGNIRKKQKGDNVSKTKSGRHRLSKTDKRLLAALVV'
                        "AGYEERTARDLIQKHVYTLTQADLRHLVSEISNGVGQSQAYDAIYQARRIRLARKYLSGKKP"
                        'EGVEPREGQEREDLP"',
                    ),
                ),
            ),
            (
                "variation",
                "6525",
                (
                    ("/note=", '"compared to AF053945 and AL109969"'),
                    ("/replace=", '"c"'),
                ),
            ),
            (
                "gene",
                "6664..7602",
                (
                    ("/gene=", '"pla"'),
                    ("/locus_tag=", '"YP_pPCP08"'),
                    ("/db_xref=", '"GeneID:2767715"'),
                ),
            ),
            (
                "CDS",
                "6664..7602",
                (
                    ("/gene=", '"pla"'),
                    ("/locus_tag=", '"YP_pPCP08"'),
                    ("/EC_number=", '"3.4.23.48"'),
                    (
                        "/note=",
                        '"outer membrane protease; involved in virulence in many'
                        " organisms; OmpT; IcsP; SopA; Pla; PgtE; omptin; in"
                        " Escherichia coli OmpT can degrade antimicrobial peptides; in"
                        " Yersinia Pla activates plasminogen during infection; in"
                        ' Shigella flexneria SopA cleaves the autotransporter IcsA"',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"outer membrane protease"'),
                    ("/protein_id=", '"NP_995574.1"'),
                    ("/db_xref=", '"GI:45478719"'),
                    ("/db_xref=", '"GeneID:2767715"'),
                    (
                        "/translation=",
                        '"MKKSSIVATIITILSGSANAASSQLIPNISPDSFTVAASTGMLSGKSHEMLYDAETGRKIS'
                        "QLDWKIKNVAILKGDISWDPYSFLTLNARGWTSLASGSGNMDDYDWMNENQSEWTDHSSHPA"
                        "TNVNHANEYDLNVKGWLLQDENYKAGITAGYQETRFSWTATGGSYSYNNGAYTGNFPKGVRV"
                        "IGYNQRFSMPYIGLAGQYRINDFELNALFKFSDWVRAHDNDEHYMRDLTFREKTSGSRYYGT"
                        "VINAGYYVTPNAKVFAEFTYSKYDEGKGGTQIIDKNSGDSVSIGGDAAGISNKNYTVTAGLQ"
                        'YRF"',
                    ),
                ),
            ),
            (
                "misc_feature",
                "6664..7599",
                (
                    ("/gene=", '"pla"'),
                    ("/locus_tag=", '"YP_pPCP08"'),
                    ("/note=", '"Omptin family; Region: Omptin; cl01886"'),
                    ("/db_xref=", '"CDD:186487"'),
                ),
            ),
            (
                "gene",
                "complement(7789..8088)",
                (("/locus_tag=", '"YP_pPCP09"'), ("/db_xref=", '"GeneID:2767713"')),
            ),
            (
                "CDS",
                "complement(7789..8088)",
                (
                    ("/locus_tag=", '"YP_pPCP09"'),
                    (
                        "/note=",
                        '"Best Blastp hit = gi|16082687|ref|NP_395234.1| (NC_003132)'
                        " putative transcriptional regulator [Yersinia pestis],"
                        " gi|5763818|emb|CAB53171.1| (AL109969) putative"
                        ' transcriptional regulator [Yersinia pestis]."',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"putative transcriptional regulator"'),
                    ("/protein_id=", '"NP_995575.1"'),
                    ("/db_xref=", '"GI:45478720"'),
                    ("/db_xref=", '"GeneID:2767713"'),
                    (
                        "/translation=",
                        '"MRTLDEVIASRSPESQTRIKEMADEMILEVGLQMMREELQLSQKQVAEAMGISQPAVTKLE'
                        'QRGNDLKLATLKRYVEAMGGKLSLDVELPTGRRVAFHV"',
                    ),
                ),
            ),
            (
                "misc_feature",
                "complement(7837..7995)",
                (
                    ("/locus_tag=", '"YP_pPCP09"'),
                    (
                        "/note=",
                        '"Helix-turn-helix XRE-family like proteins. Prokaryotic DNA'
                        " binding proteins belonging to the xenobiotic response element"
                        " family of transcriptional regulators; Region: HTH_XRE;"
                        ' cl09100"',
                    ),
                    ("/db_xref=", '"CDD:195788"'),
                ),
            ),
            (
                "gene",
                "complement(8088..8360)",
                (("/locus_tag=", '"YP_pPCP10"'), ("/db_xref=", '"GeneID:2767714"')),
            ),
            (
                "CDS",
                "complement(8088..8360)",
                (
                    ("/locus_tag=", '"YP_pPCP10"'),
                    (
                        "/note=",
                        '"Best Blastp hit = gi|16082688|ref|NP_395235.1| (NC_003132)'
                        " hypothetical protein [ Yersinia pestis],"
                        " gi|5763819|emb|CAB53172.1| (AL109969) hypothetical protein"
                        ' [Yersinia pestis]"',
                    ),
                    ("/codon_start=", "1"),
                    ("/transl_table=", "11"),
                    ("/product=", '"hypothetical protein"'),
                    ("/protein_id=", '"NP_995576.1"'),
                    ("/db_xref=", '"GI:45478721"'),
                    ("/db_xref=", '"GeneID:2767714"'),
                    (
                        "/translation=",
                        '"MADLKKLQVYGPELPRPYADTVKGSRYKNMKELRVQFSGRPIRAFYAFDPIRRAIVLCAGD'
                        'KSNDKRFYEKLVRIAEDEFTAHLNTLESK"',
                    ),
                ),
            ),
            (
                "misc_feature",
                "complement(8091..>8357)",
                (
                    ("/locus_tag=", '"YP_pPCP10"'),
                    (
                        "/note=",
                        '"Phage derived protein Gp49-like (DUF891); Region: Gp49;'
                        ' cl01470"',
                    ),
                    ("/db_xref=", '"CDD:194142"'),
                ),
            ),
            (
                "variation",
                "8529^8530",
                (("/note=", '"compared to AL109969"'), ("/replace=", '"tt"')),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_15(self):
        path = "GenBank/no_end_marker.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Premature end of file in sequence data
                record = next(records)
        length = 6497
        locus = "AB070938"
        definition = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        accession = ["AB070938"]
        titles = ()
        features = [
            (
                "source",
                "1..6497",
                (
                    ("/organism=", '"Streptomyces avermitilis"'),
                    ("/mol_type=", '"genomic DNA"'),
                    ("/db_xref=", '"taxon:33903"'),
                ),
            )
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_16(self):
        path = "GenBank/wrong_sequence_indent.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Invalid indentation for sequence line
                record = next(records)
        length = 6497
        locus = "AB070938"
        definition = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        accession = ["AB070938"]
        titles = ()
        features = [
            (
                "source",
                "1..6497",
                (
                    ("/organism=", '"Streptomyces avermitilis"'),
                    ("/mol_type=", '"genomic DNA"'),
                    ("/db_xref=", '"taxon:33903"'),
                ),
            )
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_17(self):
        path = "GenBank/invalid_locus_line_spacing.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Attempting to parse malformed locus line
                record = next(records)
        length = 6497
        locus = "AB070938"
        definition = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        accession = ["AB070938"]
        titles = ()
        features = [
            (
                "source",
                "1..6497",
                (
                    ("/organism=", '"Streptomyces avermitilis"'),
                    ("/mol_type=", '"genomic DNA"'),
                    ("/db_xref=", '"taxon:33903"'),
                ),
            )
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_18(self):
        path = "GenBank/empty_feature_qualifier.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 6497
        locus = "AB070938"
        definition = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        accession = ["AB070938"]
        titles = ()
        features = [
            (
                "source",
                "1..6497",
                (
                    ("/organism=", '"Streptomyces avermitilis"'),
                    ("/mol_type=", '"genomic DNA"'),
                    ("/db_xref=", '"taxon:33903"'),
                    ("/note=", '"This is a correct note, the following one isn\'t"'),
                    ("/note", ""),
                ),
            )
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_19(self):
        path = "GenBank/invalid_misc_feature.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: line too short to contain a feature
                record = next(records)
        length = 6497
        locus = "AB070938"
        definition = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        accession = ["AB070938"]
        titles = ()
        features = [
            (
                "source",
                "1..6497",
                (
                    ("/organism=", '"Streptomyces avermitilis"'),
                    ("/mol_type=", '"genomic DNA"'),
                    ("/db_xref=", '"taxon:33903"'),
                ),
            )
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_20(self):
        path = "GenBank/1MRR_A.gp"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 375
        locus = "1MRR_A"
        definition = (
            "Chain A, Substitution Of Manganese For Iron In Ribonucleotide Reductase"
            " From Escherichia Coli. Spectroscopic And Crystallographic"
            " Characterization"
        )
        accession = ["1MRR_A"]
        titles = (
            "Three-dimensional structure of the free radical protein of ribonucleotide"
            " reductase",
            "Substitution of manganese for iron in ribonucleotide reductase from"
            " Escherichia coli. Spectroscopic and crystallographic characterization",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..375",
                (("/organism=", '"Escherichia coli"'), ("/db_xref=", '"taxon:562"')),
            ),
            (
                "Region",
                "28..340",
                (
                    ("/region_name=", '"RNRR2"'),
                    (
                        "/note=",
                        '"Ribonucleotide Reductase, R2/beta subunit, ferritin-like'
                        ' diiron-binding domain; cd01049"',
                    ),
                    ("/db_xref=", '"CDD:153108"'),
                ),
            ),
            (
                "SecStr",
                "35..46",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 1"')),
            ),
            (
                "Site",
                "order(37,44,109..110,113,116..117,120,123,137..138,141)",
                (
                    ("/site_type=", '"other"'),
                    ("/note=", '"dimer interface [polypeptide binding]"'),
                    ("/db_xref=", '"CDD:153108"'),
                ),
            ),
            (
                "Site",
                "order(48,84,115,118,122,236..237,241)",
                (
                    ("/site_type=", '"other"'),
                    ("/note=", '"putative radical transfer pathway"'),
                    ("/db_xref=", '"CDD:153108"'),
                ),
            ),
            (
                "SecStr",
                "57..65",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 2"')),
            ),
            (
                "SecStr",
                "67..87",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 3"')),
            ),
            (
                "Site",
                "order(84,115,118,204,238,241)",
                (
                    ("/site_type=", '"other"'),
                    ("/note=", '"diiron center [ion binding]"'),
                    ("/db_xref=", '"CDD:153108"'),
                ),
            ),
            (
                "Het",
                "join(bond(84),bond(115),bond(118),bond(238))",
                (("/heterogen=", '"( MN,1000 )"'),),
            ),
            (
                "SecStr",
                "102..129",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 4"')),
            ),
            (
                "Het",
                "join(bond(115),bond(204),bond(238),bond(241))",
                (("/heterogen=", '"( MN,1001 )"'),),
            ),
            (
                "Site",
                "122",
                (
                    ("/site_type=", '"other"'),
                    ("/note=", '"tyrosyl radical"'),
                    ("/db_xref=", '"CDD:153108"'),
                ),
            ),
            (
                "SecStr",
                "133..140",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 5"')),
            ),
            (
                "SecStr",
                "143..151",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 6"')),
            ),
            (
                "SecStr",
                "153..169",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 7"')),
            ),
            (
                "SecStr",
                "172..177",
                (("/sec_str_type=", '"sheet"'), ("/note=", '"strand 1"')),
            ),
            (
                "SecStr",
                "180..185",
                (("/sec_str_type=", '"sheet"'), ("/note=", '"strand 2"')),
            ),
            (
                "SecStr",
                "186..216",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 8"')),
            ),
            ("Het", "join(bond(194),bond(272))", (("/heterogen=", '"( HG,1003 )"'),)),
            ("Het", "bond(196)", (("/heterogen=", '"( HG,1005 )"'),)),
            ("Het", "join(bond(196),bond(196))", (("/heterogen=", '"( HG,1002 )"'),)),
            (
                "Het",
                "join(bond(210),bond(214),bond(214))",
                (("/heterogen=", '"( HG,1004 )"'),),
            ),
            (
                "SecStr",
                "225..253",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 9"')),
            ),
            (
                "SecStr",
                "260..269",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 10"')),
            ),
            ("Bond", "bond(268,272)", (("/bond_type=", '"disulfide"'),)),
            (
                "SecStr",
                "270..285",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 11"')),
            ),
            (
                "Het",
                "join(bond(284),bond(305),bond(309),bond(305))",
                (("/heterogen=", '"( HG,1006 )"'),),
            ),
            (
                "SecStr",
                "301..319",
                (("/sec_str_type=", '"helix"'), ("/note=", '"helix 12"')),
            ),
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features
        )

    def test_record_parser_tsa(self):
        path = "GenBank/tsa_acropora.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 0
        locus = "GHGH01000000"
        definition = "TSA: Acropora millepora, transcriptome shotgun assembly"
        accession = ["GHGH00000000"]
        tsa = ["GHGH01000001", "GHGH01126539"]
        titles = (
            "Acropora millepora genome sequencing and assembly",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..126539",
                (
                    ("/organism=", '"Acropora millepora"'),
                    ("/mol_type=", '"transcribed RNA"'),
                    ("/db_xref=", '"taxon:45264"'),
                    ("/tissue_type=", '"late planula"'),
                    ("/country=", '"Australia: Queensland"'),
                    ("/collection_date=", '"2011"'),
                ),
            )
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features, tsa=tsa
        )

    def test_record_parser_tls(self):
        path = "GenBank/tls_KDHP01000000.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.rec_parser)
            record = next(records)
        length = 0
        locus = "KBUV01000000"
        definition = "TLS: soil metagenome 16S ribosomal RNA, targeted locus study"
        accession = ["KBUV00000000"]
        tls = ["KBUV01000001", "KBUV01003714"]
        titles = (
            "Spatio-temporal dynamics of soil bacterial communities in function of"
            " Amazon forest phenology",
            "Direct Submission",
        )
        features = [
            (
                "source",
                "1..3714",
                (
                    ("/organism=", '"soil metagenome"'),
                    ("/mol_type=", '"genomic DNA"'),
                    (
                        "/isolation_source=",
                        '"soil samples in a lowland tropical evergreen rain forest in'
                        ' Amazonia"',
                    ),
                    ("/db_xref=", '"taxon:410658"'),
                    ("/environmental_sample", ""),
                    ("/country=", '"Brazil: Manaus"'),
                    ("/lat_lon=", '"2.92 S 59.95 W"'),
                    ("/collection_date=", '"2013"'),
                    ("/note=", '"metagenomic"'),
                ),
            )
        ]
        self.perform_record_parser_test(
            record, length, locus, definition, accession, titles, features, tls=tls
        )


class TestFeatureParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.feat_parser = GenBank.FeatureParser(debug_level=0)

    def shorten(self, seq):
        try:
            s = str(seq)
        except UndefinedSequenceError:
            return None
        if len(s) <= 60:
            return s
        else:
            return s[:54] + "..." + s[-3:]

    def perform_feature_parser_test(
        self,
        record,
        seq,
        id,
        name,
        description,
        annotations,
        references,
        features,
        dbxrefs,
    ):
        self.assertEqual(self.shorten(record.seq), seq)
        self.assertEqual(record.id, id)
        self.assertEqual(record.name, name)
        self.assertEqual(record.description, description)
        references_found = []
        for key in record.annotations:
            if key == "references":
                for reference in record.annotations[key]:
                    references_found.append(str(reference))
            else:
                self.assertIn(key, annotations)
        for key in annotations:
            self.assertEqual(record.annotations[key], annotations[key])
        self.assertEqual(references_found, references)
        for feature1, (feature2, strand) in zip(record.features, features):
            self.assertEqual(str(feature1), feature2)
            self.assertEqual(feature1.strand, strand)
        self.assertEqual(record.dbxrefs, dbxrefs)

    def test_feature_parser_01(self):
        path = "GenBank/noref.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "GGCAAGATGGCGCCGGTGGGGGTGGAGAAGAAGCTGCTGCTAGGTCCCAACGGG...AAA"
        id = "NM_006141.1"
        name = "NM_006141"
        description = (
            "Homo sapiens dynein, cytoplasmic, light intermediate polypeptide 2"
            " (DNCLI2), mRNA"
        )
        annotations = {
            "accessions": ["NM_006141"],
            "comment": """\
PROVISIONAL REFSEQ: This record has not yet been subject to final
NCBI review. The reference sequence was derived from AF035812.1.""",
            "data_file_division": "PRI",
            "date": "01-NOV-2000",
            "gi": "5453633",
            "keywords": [""],
            "molecule_type": "mRNA",
            "organism": "Homo sapiens",
            "sequence_version": 1,
            "source": "human",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        }
        references = []
        features = (
            (
                """\
type: source
location: [0:1622](+)
qualifiers:
    Key: db_xref, Value: ['taxon:9606']
    Key: map, Value: ['16']
    Key: organism, Value: ['Homo sapiens']
""",
                1,
            ),
            (
                """\
type: gene
location: [0:1622](+)
qualifiers:
    Key: db_xref, Value: ['LocusID:1783']
    Key: gene, Value: ['DNCLI2']
    Key: note, Value: ['LIC2']
""",
                1,
            ),
            (
                """\
type: CDS
location: [6:1485](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['LocusID:1783', 'GI:5453634']
    Key: gene, Value: ['DNCLI2']
    Key: note, Value: ['similar to R. norvegicus and G. gallus dynein light intermediate chain 2, Swiss-Prot Accession Numbers Q62698 and Q90828, respectively']
    Key: product, Value: ['dynein, cytoplasmic, light intermediate polypeptide 2']
    Key: protein_id, Value: ['NP_006132.1']
    Key: translation, Value: ['MAPVGVEKKLLLGPNGPAVAAAGDLTSEEEEGQSLWSSILSEVSTRARSKLPSGKNILVFGEDGSGKTTLMTKLQGAEHGKKGRGLEYLYLSVHDEDRDDHTRCNVWILDGDLYHKGLLKFAVSAESLPETLVIFVADMSRPWTVMESLQKWASVLREHIDKMKIPPEKMRELERKFVKDFQDYMEPEEGCQGSPQRRGPLTSGSDEENVALPLGDNVLTHNLGIPVLVVCTKCDAVSVLEKEHDYRDEHLDFIQSHLRRFCLQYGAALIYTSVKEEKNLDLLYKYIVHKTYGFHFTTPALVVEKDAVFIPAGWDNEKKIAILHENFTTVKPEDAYEDFIVKPPVRKLVHDKELAAEDEQVFLMKQQSLLAKQPATPTRASESPARGPSGSPRTQGRGGPASVPSSSPGTSVKKPDPNIKNNAASEGVLASFFNSLLSKKTGSPGSPGAGGVQSTAKKSGQKTVLSNVQEELDRMTRKPDSMVTNSSTENEA']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_02(self):
        path = "GenBank/cor6_6.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
            seq = "AACAAAACACACATCAAAAACGATTTTACAAGAAAAAAATATCTGAAAAATGTC...AAA"
            id = "X55053.1"
            name = "ATCOR66M"
            description = "A.thaliana cor6.6 mRNA"
            annotations = {
                "accessions": ["X55053"],
                "comment": """\
Cor6.6 homologous to KIN1. KIN1 is a cold-regulated Arabidopsis
gene with suggested similarity to type I fish antifreeze proteins.""",
                "data_file_division": "PLN",
                "date": "02-MAR-1992",
                "gi": "16229",
                "keywords": [
                    "antifreeze protein homology",
                    "cold-regulated gene",
                    "cor6.6 gene",
                    "KIN1 homology",
                ],
                "molecule_type": "mRNA",
                "organism": "Arabidopsis thaliana",
                "sequence_version": 1,
                "source": "thale cress",
                "taxonomy": [
                    "Eukaryota",
                    "Viridiplantae",
                    "Streptophyta",
                    "Embryophyta",
                    "Tracheophyta",
                    "euphyllophytes",
                    "Spermatophyta",
                    "Magnoliophyta",
                    "eudicotyledons",
                    "Rosidae",
                    "Capparales",
                    "Brassicaceae",
                    "Arabidopsis",
                ],
            }
            references = [
                "location: [0:513]\nauthors: Thomashow,M.F.\ntitle: Direct"
                " Submission\njournal: Submitted (01-FEB-1991) M.F. Thomashow, Dept."
                " Crop and Soil Sciences, Dept. Microbiology, Michigan State"
                " University, East Lansing, Michigan 48824, USA\nmedline id: \npubmed"
                " id: \ncomment: \n",
                "location: [0:513]\nauthors: Gilmour,S.J., Artus,N.N. and"
                " Thomashow,M.F.\ntitle: cDNA sequence analysis and expression of two"
                " cold-regulated genes of Arabidopsis thaliana\njournal: Plant Mol."
                " Biol. 18 (1), 13-21 (1992)\nmedline id: 92119220\npubmed id:"
                " \ncomment: \n",
            ]
            features = (
                (
                    """\
type: source
location: [0:513](+)
qualifiers:
    Key: db_xref, Value: ['taxon:3702']
    Key: organism, Value: ['Arabidopsis thaliana']
    Key: strain, Value: ['Columbia']
""",
                    1,
                ),
                (
                    """\
type: gene
location: [49:250](+)
qualifiers:
    Key: gene, Value: ['cor6.6']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: [49:250](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:16230', 'SWISS-PROT:P31169']
    Key: gene, Value: ['cor6.6']
    Key: note, Value: ['cold regulated']
    Key: protein_id, Value: ['CAA38894.1']
    Key: translation, Value: ['MSETNKNAFQAGQAAGKAEEKSNVLLDKAKDAAAAAGASAQQAGKSISDAAVGGVNFVKDKTGLNK']
""",
                    1,
                ),
            )
            dbxrefs = []
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )
            record = next(records)
            seq = "ATTTGGCCTATAAATATAAACCCTTAAGCCCACATATCTTCTCAATCCATCACA...ATA"
            id = "X62281.1"
            name = "ATKIN2"
            description = "A.thaliana kin2 gene"
            annotations = {
                "accessions": ["X62281"],
                "data_file_division": "PLN",
                "date": "23-JUL-1992",
                "gi": "16353",
                "keywords": ["kin2 gene"],
                "molecule_type": "DNA",
                "organism": "Arabidopsis thaliana",
                "sequence_version": 1,
                "source": "thale cress",
                "taxonomy": [
                    "Eukaryota",
                    "Viridiplantae",
                    "Streptophyta",
                    "Embryophyta",
                    "Tracheophyta",
                    "euphyllophytes",
                    "Spermatophyta",
                    "Magnoliophyta",
                    "eudicotyledons",
                    "Rosidae",
                    "Capparales",
                    "Brassicaceae",
                    "Arabidopsis",
                ],
            }
            references = [
                "location: [0:880]\nauthors: Borg-Franck,M.E.\ntitle: Direct"
                " Submission\njournal: Submitted (27-SEP-1991) M.E. Borg-Franck, Inst"
                " of Biotechnology, University of Helsinki, Karvaamokuja 3, SF-00380"
                " Helsinki, FINLAND\nmedline id: \npubmed id: \ncomment: \n",
                "location: [0:880]\nauthors: Kurkela,S. and Borg-Franck,M.\ntitle:"
                " Structure and expression of kin2, one of two cold- and ABA-induced"
                " genes of Arabidopsis thaliana\njournal: Plant Mol. Biol. 19 (4),"
                " 689-692 (1992)\nmedline id: 92329728\npubmed id: \ncomment: \n",
            ]
            features = (
                (
                    """\
type: source
location: [0:880](+)
qualifiers:
    Key: db_xref, Value: ['taxon:3702']
    Key: organism, Value: ['Arabidopsis thaliana']
    Key: strain, Value: ['ssp. L. Heynh, Colombia']
""",
                    1,
                ),
                (
                    """\
type: TATA_signal
location: [8:20](+)
qualifiers:
""",
                    1,
                ),
                (
                    """\
type: exon
location: [43:160](+)
qualifiers:
    Key: gene, Value: ['kin2']
    Key: number, Value: ['1']
""",
                    1,
                ),
                (
                    """\
type: prim_transcript
location: [43:>579](+)
qualifiers:
    Key: gene, Value: ['kin2']
""",
                    1,
                ),
                (
                    """\
type: mRNA
location: join{[43:160](+), [319:390](+), [503:>579](+)}
qualifiers:
    Key: gene, Value: ['kin2']
""",
                    1,
                ),
                (
                    """\
type: gene
location: [43:579](+)
qualifiers:
    Key: gene, Value: ['kin2']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: join{[103:160](+), [319:390](+), [503:579](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:16354', 'SWISS-PROT:P31169']
    Key: gene, Value: ['kin2']
    Key: protein_id, Value: ['CAA44171.1']
    Key: translation, Value: ['MSETNKNAFQAGQAAGKAERRRAMFCWTRPRMLLLQLELPRNRAGKSISDAAVGGVNFVKDKTGLNK']
""",
                    1,
                ),
                (
                    """\
type: intron
location: [160:319](+)
qualifiers:
    Key: gene, Value: ['kin2']
    Key: number, Value: ['1']
""",
                    1,
                ),
                (
                    """\
type: exon
location: [319:390](+)
qualifiers:
    Key: gene, Value: ['kin2']
    Key: number, Value: ['2']
""",
                    1,
                ),
                (
                    """\
type: intron
location: [390:503](+)
qualifiers:
    Key: gene, Value: ['kin2']
    Key: number, Value: ['2']
""",
                    1,
                ),
                (
                    """\
type: exon
location: [503:>579](+)
qualifiers:
    Key: gene, Value: ['kin2']
    Key: number, Value: ['3']
""",
                    1,
                ),
                (
                    """\
type: polyA_signal
location: [619:625](+)
qualifiers:
""",
                    1,
                ),
                (
                    """\
type: polyA_signal
location: [640:646](+)
qualifiers:
""",
                    1,
                ),
                (
                    """\
type: polyA_site
location: [784:785](+)
qualifiers:
""",
                    1,
                ),
                (
                    """\
type: polyA_site
location: [799:800](+)
qualifiers:
""",
                    1,
                ),
            )
            dbxrefs = []
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )
            record = next(records)
            seq = "AAAAAAACACAACAAAACTCAATAAATAAACAAATGGCAGACAACAAGCAGAGC...TTC"
            id = "M81224.1"
            name = "BNAKINI"
            description = "Rapeseed Kin1 protein (kin1) mRNA, complete cds"
            annotations = {
                "accessions": ["M81224"],
                "data_file_division": "PLN",
                "date": "27-APR-1993",
                "gi": "167145",
                "keywords": [""],
                "molecule_type": "mRNA",
                "organism": "Brassica napus",
                "sequence_version": 1,
                "source": (
                    "Brassica napus (cultivar Jet neuf) cold induced leaf cDNA to mRNA"
                ),
                "taxonomy": [
                    "Eukaryota",
                    "Viridiplantae",
                    "Embryophyta",
                    "Tracheophyta",
                    "Spermatophyta",
                    "Magnoliophyta",
                    "eudicotyledons",
                    "core eudicots",
                    "Rosidae",
                    "eurosids II",
                    "Brassicales",
                    "Brassicaceae",
                    "Brassica",
                ],
            }
            references = [
                "location: [0:441]\nauthors: Orr,W., Iu,B., White,T., Robert,L.S. and"
                " Singh,J.\ntitle: Nucleotide sequence of a winter B. napus Kin 1"
                " cDNA\njournal: Plant Physiol. 98, 1532-1534 (1992)\nmedline id:"
                " \npubmed id: \ncomment: \n"
            ]
            features = (
                (
                    """\
type: source
location: [0:441](+)
qualifiers:
    Key: cultivar, Value: ['Jet neuf']
    Key: db_xref, Value: ['taxon:3708']
    Key: dev_stage, Value: ['cold induced']
    Key: organism, Value: ['Brassica napus']
    Key: tissue_type, Value: ['leaf']
""",
                    1,
                ),
                (
                    """\
type: gene
location: [33:300](+)
qualifiers:
    Key: gene, Value: ['kin1']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: [33:231](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:167146']
    Key: evidence, Value: ['experimental']
    Key: gene, Value: ['kin1']
    Key: protein_id, Value: ['AAA32993.1']
    Key: translation, Value: ['MADNKQSFQAGQASGRAEEKGNVLMDKVKDAATAAGASAQTAGQKITEAAGGAVNLVKEKTGMNK']
""",
                    1,
                ),
                (
                    """\
type: polyA_signal
location: [240:247](+)
qualifiers:
    Key: gene, Value: ['kin1']
    Key: note, Value: ['putative']
""",
                    1,
                ),
                (
                    """\
type: polyA_signal
location: [293:300](+)
qualifiers:
    Key: gene, Value: ['kin1']
    Key: note, Value: ['putative']
""",
                    1,
                ),
                (
                    """\
type: polyA_site
location: [440:441](+)
qualifiers:
    Key: gene, Value: ['kin1']
""",
                    1,
                ),
            )
            dbxrefs = []
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )
            record = next(records)
            seq = "GGACAAGGCCAAGGATGCTGCTGCTGCAGCTGGAGCTTCCGCGCAACAAGTAAA...GGC"
            id = "AJ237582.1"
            name = "ARU237582"
            description = "Armoracia rusticana csp14 gene (partial), exons 2-3"
            annotations = {
                "accessions": ["AJ237582"],
                "data_file_division": "PLN",
                "date": "24-MAR-1999",
                "gi": "4538892",
                "keywords": ["cold shock protein", "csp14 gene"],
                "molecule_type": "DNA",
                "organism": "Armoracia rusticana",
                "sequence_version": 1,
                "source": "horseradish",
                "taxonomy": [
                    "Eukaryota",
                    "Viridiplantae",
                    "Streptophyta",
                    "Embryophyta",
                    "Tracheophyta",
                    "euphyllophytes",
                    "Spermatophyta",
                    "Magnoliophyta",
                    "eudicotyledons",
                    "Rosidae",
                    "Capparales",
                    "Brassicaceae",
                    "Armoracia",
                ],
            }
            references = [
                "location: [0:206]\nauthors: Baymiev,A.K., Gimalov,F.R. and"
                " Vakhitov,V.A.\ntitle: \njournal: Unpublished\nmedline id: \npubmed"
                " id: \ncomment: \n",
                "location: [0:206]\nauthors: Baymiev,A.K.\ntitle: Direct"
                " Submission\njournal: Submitted (20-MAR-1999) Baymiev A.K.,"
                " Departament of Biochemistry and Cytochemistry, Ufa Scientific Centre,"
                " pr. Oktyabrya 69, Ufa, Bashkortostan, Russia, 450054, RUSSIA\nmedline"
                " id: \npubmed id: \ncomment: \n",
            ]
            features = (
                (
                    """\
type: source
location: [0:206](+)
qualifiers:
    Key: country, Value: ['Russia:Bashkortostan']
    Key: db_xref, Value: ['taxon:3704']
    Key: organism, Value: ['Armoracia rusticana']
""",
                    1,
                ),
                (
                    """\
type: mRNA
location: join{[<0:48](+), [142:>206](+)}
qualifiers:
    Key: gene, Value: ['csp14']
""",
                    1,
                ),
                (
                    """\
type: exon
location: [0:48](+)
qualifiers:
    Key: gene, Value: ['csp14']
    Key: number, Value: ['2']
""",
                    1,
                ),
                (
                    """\
type: gene
location: [0:206](+)
qualifiers:
    Key: gene, Value: ['csp14']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: join{[<0:48](+), [142:>206](+)}
qualifiers:
    Key: codon_start, Value: ['2']
    Key: db_xref, Value: ['GI:4538893']
    Key: gene, Value: ['csp14']
    Key: product, Value: ['cold shock protein']
    Key: protein_id, Value: ['CAB39890.1']
    Key: translation, Value: ['DKAKDAAAAAGASAQQAGKNISDAAAGGVNFVKEKTG']
""",
                    1,
                ),
                (
                    """\
type: intron
location: [48:142](+)
qualifiers:
    Key: gene, Value: ['csp14']
    Key: number, Value: ['2']
""",
                    1,
                ),
                (
                    """\
type: exon
location: [142:206](+)
qualifiers:
    Key: gene, Value: ['csp14']
    Key: number, Value: ['3']
""",
                    1,
                ),
            )
            dbxrefs = []
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )
            record = next(records)
            seq = "AACAAAACTCAATAAATAAACAAATGGCAGACAACAAGCAGAGCTTCCAAGCCG...TTT"
            id = "L31939.1"
            name = "BRRBIF72"
            description = "Brassica rapa (clone bif72) kin mRNA, complete cds"
            annotations = {
                "accessions": ["L31939"],
                "data_file_division": "PLN",
                "date": "01-MAR-1996",
                "gi": "1209261",
                "keywords": [""],
                "molecule_type": "mRNA",
                "organism": "Brassica rapa",
                "sequence_version": 1,
                "source": "Brassica rapa flower cDNA to mRNA",
                "taxonomy": [
                    "Eukaryota",
                    "Viridiplantae",
                    "Embryophyta",
                    "Tracheophyta",
                    "Spermatophyta",
                    "Magnoliophyta",
                    "eudicotyledons",
                    "core eudicots",
                    "Rosidae",
                    "eurosids II",
                    "Brassicales",
                    "Brassicaceae",
                    "Brassica",
                ],
            }
            references = [
                "location: [0:282]\nauthors: Kim,J.-B., Kim,H.-U., Park,B.-S.,"
                " Yun,C.-H., Cho,W.-S., Ryu,J.-C. and Chung,T.-Y.\ntitle: Nucleotide"
                " sequences of kin gene in chinese cabbage\njournal: Unpublished"
                " (1994)\nmedline id: \npubmed id: \ncomment: \n"
            ]
            features = (
                (
                    """\
type: source
location: [0:282](+)
qualifiers:
    Key: db_xref, Value: ['taxon:3711']
    Key: dev_stage, Value: ['flower']
    Key: organism, Value: ['Brassica rapa']
""",
                    1,
                ),
                (
                    """\
type: gene
location: [23:221](+)
qualifiers:
    Key: gene, Value: ['kin']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: [23:221](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:1209262']
    Key: gene, Value: ['kin']
    Key: protein_id, Value: ['AAA91051.1']
    Key: translation, Value: ['MADNKQSFQAGQAAGRAEEKGNVLLMDKVKDAATAAGALQTAGQKITEAAGGAVNLVKEKTGMNK']
""",
                    1,
                ),
            )
            dbxrefs = []
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )
            record = next(records)
            seq = "ATGGCAGACAACAAGCAGAGCTTCCAAGCCGGTCAAGCCGCTGGTCGTGCTGAG...TAG"
            id = "AF297471.1"
            name = "AF297471"
            description = "Brassica napus BN28a (BN28a) gene, complete cds"
            annotations = {
                "accessions": ["AF297471"],
                "data_file_division": "PLN",
                "date": "14-SEP-2000",
                "gi": "10121868",
                "keywords": [""],
                "molecule_type": "DNA",
                "organism": "Brassica napus",
                "sequence_version": 1,
                "source": "rape",
                "taxonomy": [
                    "Eukaryota",
                    "Viridiplantae",
                    "Embryophyta",
                    "Tracheophyta",
                    "Spermatophyta",
                    "Magnoliophyta",
                    "eudicotyledons",
                    "core eudicots",
                    "Rosidae",
                    "eurosids II",
                    "Brassicales",
                    "Brassicaceae",
                    "Brassica",
                ],
            }
            references = [
                "location: [0:497]\nauthors: Byass,L.J. and Flanagan,A.M.\ntitle:"
                " BN28a, a low temperature-induced gene of Brassica napus\njournal:"
                " Unpublished\nmedline id: \npubmed id: \ncomment: \n",
                "location: [0:497]\nauthors: Byass,L.J. and Flanagan,A.M.\ntitle:"
                " Direct Submission\njournal: Submitted (18-AUG-2000) AFNS, University"
                " of Alberta, 4-10 Agriculture/Forestry Centre, Edmonton, Alberta T6G"
                " 2P5, Canada\nmedline id: \npubmed id: \ncomment: \n",
            ]
            features = (
                (
                    """\
type: source
location: [0:497](+)
qualifiers:
    Key: cultivar, Value: ['Cascade']
    Key: db_xref, Value: ['taxon:3708']
    Key: organism, Value: ['Brassica napus']
""",
                    1,
                ),
                (
                    """\
type: mRNA
location: join{[<0:54](+), [240:309](+), [422:>497](+)}
qualifiers:
    Key: gene, Value: ['BN28a']
    Key: product, Value: ['BN28a']
""",
                    1,
                ),
                (
                    """\
type: gene
location: [<0:>497](+)
qualifiers:
    Key: gene, Value: ['BN28a']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: join{[0:54](+), [240:309](+), [422:497](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:10121869']
    Key: gene, Value: ['BN28a']
    Key: note, Value: ['low temperature-induced; similar to Brassica napus Kin1 in Accession Number M81224']
    Key: product, Value: ['BN28a']
    Key: protein_id, Value: ['AAG13407.1']
    Key: translation, Value: ['MADNKQSFQAGQAAGRAEEKGNVLMDKVKDAATAAGASAQTAGQKITEAAGGAVNLVKEKTGMNK']
""",
                    1,
                ),
            )
            dbxrefs = []
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )

    def test_feature_parser_03(self):
        path = "GenBank/iro.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "CACAGGCCCAGAGCCACTCCTGCCTACAGGTTCTGAGGGCTCAGGGGACCTCCT...AAA"
        id = "AL109817.1"
        name = "IRO125195"
        description = "Homo sapiens mRNA full length insert cDNA clone EUROIMAGE 125195"
        annotations = {
            "accessions": ["AL109817"],
            "comment": """\
EURO-IMAGE Consortium Contact: Auffray C
CNRS UPR 420 - Genetique Moleculaire et Biologie du Developement
IFR 1221 - Rue Guy Moquet 19, Batiment G - BP 8
94801 Villejuif  Cedex, FRANCE
Tel: ++33-1-49 58 34 98
Fax: ++33-1-49 58 35 09
e-mail: auffray@infobiogen.fr
This clone is available royalty-free through IMAGE Consortium
Distributors.
IMPORTANT: This sequence represents the full insert of this IMAGE
cDNA clone. No attempt has been made to verify whether this
corresponds to the full-length of the original mRNA from which it
was derived.""",
            "data_file_division": "PRI",
            "date": "11-AUG-1999",
            "gi": "5731880",
            "keywords": ["FLI_CDNA"],
            "molecule_type": "mRNA",
            "organism": "Homo sapiens",
            "sequence_version": 1,
            "source": "human",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        }
        references = [
            "location: [0:1326]\nauthors: Auffray,C., Ansorge,W., Ballabio,A.,"
            " Estivill,X., Gibson,K., Lehrach,H., Poustka,A. and Lundeberg,J.\ntitle:"
            " The European IMAGE consortium for integrated Molecular analysis of human"
            " gene transcripts\njournal: Unpublished\nmedline id: \npubmed id:"
            " \ncomment: \n",
            "location: [0:1326]\nauthors: Carim,L., Estivill,X., Sumoy,L. and"
            " Escarceller,M.\ntitle: Direct Submission\njournal: Submitted"
            " (11-AUG-1999) Dept. Genetica Molecular, Institut de Recerca Oncologica"
            " (IRO), Hospital Duran i Reynals, Autovia de Castelldefels Km 2,7"
            " L'Hospitalet de Llobregat, 08907 Barcelona, Catalunya, SPAIN. Tel:"
            " ++34-93-260-7775 Fax: ++34-93-260-7776 WWW site: http://www.iro.es e-mail"
            " enquiries: lsumoy@iro.es, mescarceller@iro.es\nmedline id: \npubmed id:"
            " \ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:1326](+)
qualifiers:
    Key: chromosome, Value: ['21']
    Key: clone, Value: ['IMAGE cDNA clone 125195']
    Key: clone_lib, Value: ['Soares fetal liver spleen 1NFLS']
    Key: db_xref, Value: ['taxon:9606']
    Key: note, Value: ['contains Alu repeat; likely to be be derived from unprocessed nuclear RNA or genomic DNA; encodes putative exons identical to FTCD; formimino transferase cyclodeaminase; formimino transferase (EC 2.1.2.5) /formimino tetrahydro folate cyclodeaminase (EC 4.3.1.4)']
    Key: organism, Value: ['Homo sapiens']
""",
                1,
            ),
            (
                """\
type: gene
location: [340:756](+)
qualifiers:
    Key: gene, Value: ['FTCD']
""",
                1,
            ),
            (
                """\
type: exon
location: [340:384](+)
qualifiers:
    Key: gene, Value: ['FTCD']
    Key: number, Value: ['1']
""",
                1,
            ),
            (
                """\
type: intron
location: [384:617](+)
qualifiers:
    Key: gene, Value: ['FTCD']
    Key: number, Value: ['1']
""",
                1,
            ),
            (
                """\
type: exon
location: [617:756](+)
qualifiers:
    Key: gene, Value: ['FTCD']
    Key: number, Value: ['2']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_04(self):
        path = "GenBank/pri1.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "GATCATGCATGCACTCCAGCCTGGGACAAGAGCGAAACTCCGTCTCAAAAAAAA...GCA"
        id = "U05344.1"
        name = "HUGLUT1"
        description = "Human fructose transporter (GLUT5) gene, promoter and exon 1"
        annotations = {
            "accessions": ["U05344"],
            "data_file_division": "PRI",
            "date": "16-NOV-1994",
            "gi": "452475",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Homo sapiens",
            "sequence_version": 1,
            "source": "human",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        }
        references = [
            "location: [0:741]\nauthors: Mahraoui,L., Takeda,J., Mesonero,J.,"
            " Chantret,I., Dussaulx,E., Bell,G.I. and Brot-Laroche,E.\ntitle:"
            " Regulation of expression of the human fructose transporter (GLUT5) by"
            " cyclic AMP\njournal: Biochem. J. 301 (Pt 1), 169-175 (1994)\nmedline id:"
            " 94311827\npubmed id: \ncomment: \n",
            "location: [0:741]\nauthors: Takeda,J.\ntitle: Direct Submission\njournal:"
            " Submitted (24-JAN-1994) Jun Takeda, Howard Hughes Medical Institute, The"
            " University of Chicago, 5841 S. Maryland Ave., Chicago, IL 60637,"
            " USA\nmedline id: \npubmed id: \ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:741](+)
qualifiers:
    Key: chromosome, Value: ['1']
    Key: clone, Value: ['lambda hGT5-157']
    Key: clone_lib, Value: ['partial Hae III/Alu I fetal human liver library in lambda Ch4A of Maniatis']
    Key: db_xref, Value: ['taxon:9606']
    Key: dev_stage, Value: ['fetal']
    Key: map, Value: ['1p31']
    Key: organism, Value: ['Homo sapiens']
    Key: tissue_type, Value: ['liver']
""",
                1,
            ),
            (
                """\
type: repeat_region
location: [0:73](+)
qualifiers:
    Key: rpt_family, Value: ['Alu']
""",
                1,
            ),
            (
                """\
type: promoter
location: [0:513](+)
qualifiers:
""",
                1,
            ),
            (
                """\
type: 5'UTR
location: [513:609](+)
qualifiers:
    Key: gene, Value: ['GLUT5']
""",
                1,
            ),
            (
                """\
type: exon
location: [513:642](+)
qualifiers:
    Key: gene, Value: ['GLUT5']
    Key: number, Value: ['1']
    Key: product, Value: ['fructose transporter']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_05(self):
        path = "GenBank/arab1.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "AAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTATTGTAACCTTAGG...CTT"
        id = "AC007323.5"
        name = "AC007323"
        description = (
            "Genomic sequence for Arabidopsis thaliana BAC T25K16 from chromosome I,"
            " complete sequence"
        )
        annotations = {
            "accessions": ["AC007323"],
            "comment": """\
On Dec 16, 1999 this sequence version replaced gi:5729683.""",
            "data_file_division": "PLN",
            "date": "19-JAN-2000",
            "gi": "6587720",
            "keywords": ["HTG"],
            "molecule_type": "DNA",
            "organism": "Arabidopsis thaliana",
            "sequence_version": 5,
            "source": "thale cress",
            "taxonomy": [
                "Eukaryota",
                "Viridiplantae",
                "Embryophyta",
                "Tracheophyta",
                "Spermatophyta",
                "Magnoliophyta",
                "eudicotyledons",
                "core eudicots",
                "Rosidae",
                "eurosids II",
                "Brassicales",
                "Brassicaceae",
                "Arabidopsis",
            ],
        }
        references = [
            "location: [0:86436]\nauthors: Dunn,P., Shinn,P., Brooks,S., Buehler,E.,"
            " Chao,Q., Johnson-Hopson,C., Khan,S., Kim,C., Altafi,H., Bei,Q., Chin,C.,"
            " Chiou,J., Choi,E., Conn,L., Conway,A., Gonzales,A., Hansen,N., Howing,B.,"
            " Koo,T., Lam,B., Lee,J., Lenz,C., Li,J., Liu,A., Liu,K., Liu,S.,"
            " Mukharsky,N., Nguyen,M., Palm,C., Pham,P., Sakano,H., Schwartz,J.,"
            " Southwick,A., Thaveri,A., Toriumi,M., Vaysberg,M., Yu,G.,"
            " Federspiel,N.A., Theologis,A. and Ecker,J.R.\ntitle: Genomic sequence for"
            " Arabidopsis thaliana BAC T25K16 from chromosome I\njournal:"
            " Unpublished\nmedline id: \npubmed id: \ncomment: \n",
            "location: [0:86436]\nauthors: Ecker,J.R.\ntitle: Direct"
            " Submission\njournal: Submitted (17-APR-1999) Arabidopsis thaliana Genome"
            " Center, Department of Biology, University of Pennsylvania, 38th Street"
            " and Hamilton Walk, Philadelphia, Pennsylvania 19104-6018, USA\nmedline"
            " id: \npubmed id: \ncomment: \n",
            "location: [0:86436]\nauthors: Ecker,J.R.\ntitle: Direct"
            " Submission\njournal: Submitted (11-AUG-1999) Arabidopsis thaliana Genome"
            " Center, Department of Biology, University of Pennsylvania, 38th Street"
            " and Hamilton Walk, Philadelphia, Pennsylvania 19104-6018, USA\nmedline"
            " id: \npubmed id: \ncomment: \n",
            "location: [0:86436]\nauthors: Ecker,J.R.\ntitle: Direct"
            " Submission\njournal: Submitted (16-DEC-1999) Arabidopsis thaliana Genome"
            " Center, Department of Biology, University of Pennsylvania, 38th Street"
            " and Hamilton Walk, Philadelphia, Pennsylvania 19104-6018, USA\nmedline"
            " id: \npubmed id: \ncomment: \n",
            "location: [0:86436]\nauthors: Chao,Q., Brooks,S., Buehler,E.,"
            " Johnson-Hopson,C., Khan,S., Kim,C., Shinn,P., Altafi,H., Bei,B., Chin,C.,"
            " Chiou,J., Choi,E., Conn,L., Conway,A., Gonzalez,A., Hansen,N., Howing,B.,"
            " Koo,T., Lam,B., Lee,J., Lenz,C., Li,J., Liu,A., Liu,J., Liu,S.,"
            " Mukharsky,N., Nguyen,M., Palm,C., Pham,P., Sakano,H., Schwartz,J.,"
            " Southwick,A., Thaveri,A., Toriumi,M., Vaysberg,M., Yu,G., Davis,R.,"
            " Federspiel,N., Theologis,A. and Ecker,J.\ntitle: Direct"
            " Submission\njournal: Submitted (19-JAN-2000) Arabidopsis thaliana Genome"
            " Center, Department of Biology, University of Pennsylvania, 38th and"
            " Hamilton Walk, Philadelphia, PA 19104-6018, USA\nmedline id: \npubmed id:"
            " \ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:86436](+)
qualifiers:
    Key: chromosome, Value: ['1']
    Key: clone, Value: ['T25K16']
    Key: db_xref, Value: ['taxon:3702']
    Key: organism, Value: ['Arabidopsis thaliana']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[3461:3615](+), [3697:3978](+), [4076:4307](+), [4407:4797](+), [4875:5028](+), [5140:5332](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715633']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['containing similarity to NAM-like proteins gi|3695378']
    Key: product, Value: ['T25K16.1']
    Key: protein_id, Value: ['AAF26460.1']
    Key: translation, Value: ['MEDQVGFGFRPNDEELVGHYLRNKIEGNTSRDVEVAISEVNICSYDPWNLRFQSKYKSRDAMWYFFSRRENNKGNRQSRTTVSGKWKLTGESVEVKDQWGFCSEGFRGKIGHKRVLVFLDGRYPDKTKSDWVIHEFHYDLLPEHQKLCNVTLFRFSSYFRLSLLSPMFYTDELMCLPPEILQRTYVICRLEYKGDDADILSAYAIDPTPAFVPNMTSSAGSVVNQSRQRNSGSYNTYSEYDSANHGQQFNENSNIMQQQPLQGSFNPLLEYDFANHGGQWLSDYIDLQQQVPYLAPYENESEMIWKHVIEENFEFLVDERTSMQQHYSDHRPKKPVSGVLPDDSSDTETGSMIFEDTSSSTDSVGSSDEPGHTRIDDIPSLNIIEPLHNYKAQEQPKQQSKEKVISSQKSECEWKMAEDSIKIPPSTNTVKQSWIVLENAQWNYLKNMIIGVLLFISVISWIILVG']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[8272:8368](-), [8086:8166](-), [7915:7998](-), [7463:7603](-), [7265:7351](-), [6616:6953](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715650']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['hypothetical protein']
    Key: product, Value: ['T25K16.2']
    Key: protein_id, Value: ['AAF26477.1']
    Key: translation, Value: ['MAASEHRCVGCGFRVKSLFIQYSPGNIRLMKCGNCKEVADEYIECERMVCFNHFLSLFGPKVYRHVLYNAINPATVNIQVKNYFNSTSRCVVGEIHRQTYLKSPELIIDRSLLLRKSDEESSFSDSPVLLSIKVLIGVLSANAAFIISFAIATKGLLNEVSRESLLLQVWEFPMSVIFFVDILLLTSNSMALKGQTFKMFSMQIVFCCCYFGISQCKFVFKPVMTESTMTRCIAVCLIAHLIRFLVGQIFEPTIFLIQIGSLLQYMSYFFRIV']
""",
                -1,
            ),
            (
                """\
type: CDS
location: [11565:12642](-)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715649']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['putative RAP2.8 protein gi|3695373']
    Key: product, Value: ['T25K16.3']
    Key: protein_id, Value: ['AAF26476.1']
    Key: translation, Value: ['MDLSLAPTTTTSSDQEQDRDQELTSNIGASSSSGPSGNNNNLPMMMIPPPEKEHMFDKVVTPSDVGKLNRLVIPKQHAERYFPLDSSNNQNGTLLNFQDRNGKMWRFRYSYWNSSQSYVMTKGWSRFVKEKKLDAGDIVSFQRGIGDESERSKLYIDWRHRPDMSLVQAHQFGNFGFNFNFPTTSQYSNRFHPLPEYNSVPIHRGLNIGNHQRSYYNTQRQEFVGYGYGNLAGRCYYTGSPLDHRNIVGSEPLVIDSVPVVPGRLTPVMLPPLPPPPSTAGKRLRLFGVNMECGNDYNQQEESWLVPRGEIGASSSSSSALRLNLSTDHDDDNDDGDDGDDDQFAKKGKSSLSLNFNP']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[23220:24174](+), [24243:24357](+), [24411:24664](+), [24742:25137](+), [25225:25445](+), [25526:25711](+), [25782:25905](+), [25993:26478](+), [26563:26730](+), [26813:26983](+), [27073:27235](+), [27319:27415](+), [27504:28133](+), [28313:28507](+), [28591:28782](+), [28861:30013](+), [30111:30518](+), [30603:30781](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715634']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['similar to UFD1 protein emb|CAB10321.1; similar to ESTs gb|H36434, gb|AI996152.1']
    Key: product, Value: ['T25K16.4']
    Key: protein_id, Value: ['AAF26461.1']
    Key: translation, Value: ['MVMEDEPREATIKPSYWLDACEDISCDLIDDLVSEFDPSSVAVNESTDENGVINDFFGGIDHILDSIKNGGGLPNNGVSDTNSQINEVTVTPQVIAKETVKENGLQKNGGKRDEFSKEEGDKDRKRARVCSYQSERSNLSGRGHVNNSREGDRFMNRKRTRNWDEAGNNKKKRECNNYRRDGRDREVRGYWERDKVGSNELVYRSGTWEADHERDVKKVSGGNRECDVKAEENKSKPEERKEKVVEEQARRYQLDVLEQAKAKNTIAFLETGAGKTLIAILLIKSVHKDLMSQNRKMLSVFLVPKVPLVYQVPPNKKHQAEVIRNQTCFQVGHYCGEMGQDFWDSRRWQREFESKQFLKLTSFFLFSSTQVLVMTAQILLNILRHSIIRMETIDLLILDECHHAVKKHPYSLVMSEFYHTTPKDKRPAIFGMTASPVNLKGVSSQVDCAIKIRNLETKLDSTVCTIKDRKELEKHVPMPSEIVVEYDKAATMWSLHETIKQMIAAVEEAAQASSRKSKWQFMGARDAGAKDELRQVYGVSERTESDGAANLIHKLRAINYTLAELGQWCAYKVGQSFLSALQSDERVNFQVDVKFQESYLSEVVSLLQCELLEGAAAEKVAAEVGKPENGNAHDEMEEGELPDDPVVSGGEHVDEVIGAAVADGKVTPKVQSLIKLLLKYQHTADFRAIVFVERVVAALVLPKVRIKVFAELPSLSFIRCASMIGHNNSQEMKSSQMQDTISKFRDGHVTLLVATSVAEEGLDIRQCNVVMRFDLAKTVLAYIQSRGRARKPGSDYILMVERYIKSFKNYILIFVTTGHQISTDMSTCVTCRGNVSHAAFLRNARNSEETLRKEAIERTDLSHLKDTSRLISIDAVPGTVYKVEATGAMVSLNSAVGLVHFYCSQLPGDRYAILRPEFSMEKHEKPGGHTEYSCRLQLPCNAPFEILEGPVCSSMRLAQQVDIIVSACKKLHEMGAFTDMLLPDKGSGQDAEKADQDDEGEPVPGTARHREFYPEGVADVLKGEWVSSGKEVCESSKLFHLYMYNVRCVDFGSSKDPFLSEVSEFAILFGNELDAEVLSMSMDLYVARAMITKASLAFKGSLDITENQLSSLKKFHVRLMSIVLDVDVEPSTTPWDPAKAYLFVPVTDNTSMEPIKGINWELVEKITKTTAWDNPLQRARPDVYLGTNERTLGGDRREYGFGKLRHNIVFGQKSHPTYGIRGAVASFDVVRASGLLPVRDAFEKEVEEDLSKGKLMMADGCMVAEDLIGKIVTAAHSGKRFYVDSICYDMSAETSFPRKEGYLGPLEYNTYADYYKQKIYVVQDRLFFYFLHNLRLLRLYKSSSIMLFIRYGVDLNCKQQPLIKGRGVSYCKNLLSPRFEQSGESETVLDKTYYVFLPPELCVVHPLSGSLIRGAQRLPSIMRRVESMLLAVQLKNLISYPIPTSKILEALTAASCQETFCYERAELLGDAYLKWVVSRFLFLKYPQKHEGQLTRMRQQMVSNMVLYQFALVKGLQSYIQADRFAPSRWSAPGVPPVFDEDTKDGGSSFFDEEQKPVSEENSDVFEDGEMEDGELEGDLSSYRVLSSKTLADVVEALIGVYYVEGGKIAANHLMKWIGIHVEDDPDEVDGTLKNVNVPESVLKSIDFVGLERALKYEFKEKGLLVEAITHASRPSSGVSCYQRLEFVGDAVLDHLITRHLFFTYTSLPPGRLTDLRAAAVNNENFARVAVKHKLHLYLRHGSSALEKQVNKIKKQSILFSKSFKCLTVWLLFVFQIREFVKEVQTESSKPGFNSFGLGDCKAPKVLGDIVESIAGAIFLDSGKDTTAAWKVFQPLLQPMVTPETLPMHPVRELQERCQQQAEGLEYKASRSGNTATVEVFIDGVQVGVAQNPQKKMAQKLAARNALAALKEKEIAESKEKHINNGNAGEDQGENENGNKKNGHQPFTRQTLNDICLRKNWPMPSYRCVKEGGPAHAKRFTFGVRVNTSDRGWTDECIGEPMPSVKKAKDSAAVLLLELLNKTFS']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[32248:32372](-), [32132:32161](-), [31983:32049](-), [31789:31897](-), [31634:31700](-), [31340:31515](-), [31222:31304](-), [31083:31126](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715648']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['putative inorganic pyrophosphatase gi|3510259; similar to ESTs gb|T42316, gb|AI994042.1, gb|AI994013.1, emb|Z29202']
    Key: product, Value: ['T25K16.5']
    Key: protein_id, Value: ['AAF26475.1']
    Key: translation, Value: ['MSEETKDNQRLQRPAPRLNERILSSLSRRSVAAHPWHDLEIGPGAPQIFNVVVEITKGSKVKYELDKKTGLIKVDRILYSSVVYPHNYGFVPRTLCEDNDPIDVLVIMQEPVLPGCFLRARAIGLMPMIDQGEKDDKIIAVCVDDPEYKHYTDIKELPPHRLSEIRRFFEDCILFLQCSSLFISIDLSTNKKNENKEVAVNDFLPSESAVEAIQYSMDLYAEYILHTLRR']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[36724:36763](-), [36511:36623](-), [36325:36387](-), [35431:35701](-), [35268:35349](-), [34102:35173](-), [33693:34029](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715647']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['putative late elongated hypocotyl emb|CAA07004; similar to ESTS gb|AI993521.1, gb|AA650979']
    Key: product, Value: ['T25K16.6']
    Key: protein_id, Value: ['AAF26474.1']
    Key: translation, Value: ['MDTNTSGEELLAKARKPYTITKQRERWTEDEHERFLEALRLYGRAWQRIEEHIGTKTAVQIRSHAQKFFTKFGKAHSFWFTFQLEKEAEVKGIPVCQALDIEIPPPRPKRKPNTPYPRKPGNNGTSSSQVSSAKDAKLVSSASSSQLNQAFLDLEKMPFSEKTSTGKENQDENCSGVSTVNKYPLPTKVSGDIETSKTSTVDNAVQDVPKKNKDKDGNDGTTVHSMQNYPWHFHADIVNGNIAKCPQNHPSGMVSQDFMFHPMREETHGHANLQATTASATTTASHQAFPACHSQDDYRSFLQISSTFSNLIMSTLLQNPAAHAAATFAASVWPYASVGNSGDSSTPMSSSPPSITAIAAATVAAATAWWASHGLLPVCAPAPITCVPFSTVAVPTPAMTEMDTVENTQPFEKQNTALQDQNLASKSPASSSDDSDETGVTKLNADSKTNDDKIEEVVVTAAVHDSNTAQKKNLVDRSSCGSNTPSGSDAETDALDKMEKDKEDVKETDENQPDVIELNNRKIKMRDNNSNNNATTDSWKEVSEEGRIAFQALFARERLPQSFSPPQVAENVNRKQSDTSMPLAPNFKSQDSCAADQEGVVMIGVGTCKSLKTRQTGFKPYKRCSMEVKESQVGNINNQSDEKVCKRLRLEGEAST']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[40376:40579](-), [39914:40031](-), [39110:39516](-), [38837:38989](-), [38599:38756](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715646']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['similar to Medicago truncatula MtN2 gi|3193308; similar to EST gb|H77065']
    Key: product, Value: ['T25K16.7']
    Key: protein_id, Value: ['AAF26473.1']
    Key: translation, Value: ['MAGDMQGVRVVEKYSPVIVMVMSNVAMGSVNALVKKALDVGVNHMVIGAYRMAISALILVPFAYVLERASLMQFFFLLGLSYTSATVSCALVSMLPAITFALALIFRTENVKILKTKAGMLKVIGTLICISGALFLTFYKGPQISNSHSHSHGGASHNNNDQDKANNWLLGCLYLTIGTVLLSLWMLFQGTLSIKYPCKYSSTCLMSIFAAFQCALLSLYKSRDVNDWIIDDRFVITVIIYAGVVGQAMTTVATTWGIKKLGAVFASAFFPLTLISATLFDFLILHTPLYLGSVIGSLVTITGLYMFLWGKNKETESSTALSSGMDNEAQYTTPNKDNDSKSPV']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[48637:48868](-), [47776:48554](-), [47447:47684](-), [46074:46313](-), [45718:45847](-), [45342:45656](-), [45149:45261](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715645']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['putative pyruvate dehydrogenase E1 alpha subunit gi|2454182; similar to ESTs emb|Z48417, gb|AW039459.1, gb|T15146, emb|Z48416, gb|AF066871, gb|T76832, gb|AI996061.1']
    Key: product, Value: ['T25K16.8']
    Key: protein_id, Value: ['AAF26472.1']
    Key: translation, Value: ['MATAFAPTKLTATVPLHGSHENRLLLPIRLAPPSSFLGSTRSLSLRRLNHSNATRRSPVVSVQEVVKEKQSTNNTSLLITKEEGLELYEDMILGRSFEDMCAQMYYRGKMFGFVHLYNGQEAVSTGFIKLLTKSDSVVSTYRDHVHALSKGVSARAVMSELFGKVTGCCRGQGGSMHMFSKEHNMLGGFAFIGEGIPVATGAAFSSKYRREVLKQDCDDVTVAFFGDGTCNNGQFFECLNMAALYKLPIIFVVENNLWAIGMSHLRATSDPEIWKKGPAFGMPGVHVDGMDVLKVREVAKEAVTRARRGEGPTLVECETYRFRGHSLADPDELRDAAEKAKYAARDPIAALKKYLIENKLAKEAELKSIEKKIDELVEEAVEFADASPQPGRSQLLENVFADPKGFGIGPDGRYRSQPLQIKVSSSELSVLDEEKEEEVVKGEAEPNKDSVVSKAEPVKKPRPCELYVCNIPRSYDIAQLLDMFQPFGTVISVEVVSRNPQTGESRGSGYVTMGSINSAKIAIASLDGTVRARETKKQEVGGREMRVRYSVDMNPGTRRNPEVLNSTPKKILMYESQHKVYVGNLPWFTQPDGLRNHFSKFGTIVSTRVLHDRKTGRNRVFAFLSFTSGEERDAALSFNGTVNNMKVAESSSEKVSRRVSRKPTVLLLLQRHLLDTNNV']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[50584:50656](-), [50120:50333](-), [49985:50039](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715644']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['similar to acidic ribosomal protein p1 gi|2252857; similar to ESTs gb|T42111, gb|AI099979, gb|AA728491']
    Key: product, Value: ['T25K16.9']
    Key: protein_id, Value: ['AAF26471.1']
    Key: translation, Value: ['MSTVGELACSYAVMILEDEGIAITADKIATLVKAAGVSIESYWPMLFAKMAEKRNVTDLIMNVGAGGGGGAPVAAAAPAAGGGAAAAPAAEEKKKDEPAEESDGDLGFGLFD']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[51940:52048](+), [52135:52432](+), [52639:52885](+), [53185:53326](+), [53404:54196](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715635']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['hypothetical protein']
    Key: product, Value: ['T25K16.10']
    Key: protein_id, Value: ['AAF26462.1']
    Key: translation, Value: ['MGKKNGSSSWLTAVKRAFRSPTKKDHSNDVEEDEEKKREKRRWFRKPATQESPVKSSGISPPAPQEDSLNVNSKPSPETAPSYATTTPPSNAGKPPSAVVPIATSASKTLAPRRIYYARENYAAVVIQTSFRGYLARRALRALKGLVKLQALVRGHNVRKQAKMTLRCMQALVRVQSRVLDQRKRLSHDGSRKSAFSDSHAVFESRYLQDLSDRQSMSREGSSAAEDWDDRPHTIDAVKVMLQRRRDTALRHDKTNLSQAFSQKMWRTVGNQSTEGHHEVELEEERPKWLDRWMATRPWDKRASSRASVDQRVSVKTVEIDTSQPYSRTGAGSPSRGQRPSSPSRTSHHYQSRNNFSATPSPAKSRPILIRSASPRCQRDPREDRDRAAYSYTSNTPSLRSNYSFTARSGCSISTTMVNNASLLPNYMASTESAKARIRSHSAPRQRPSTPERDRAGLVKKRLSYPVPPPAEYEDNNSLRSPSFKSVAGSHFGGMLEQQSNYSSCCTESNGVEISPASTSDFRNWLR']
""",
                1,
            ),
            (
                """\
type: CDS
location: [57093:58680](-)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715643']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['putative fatty acid elongase 3-ketoacyl-coA synthase 1 gi|4091810; similar to ESTs gb|T42377, gb|N96054, gb|T44368, gb|AI999379.1, emb|Z26005']
    Key: product, Value: ['T25K16.11']
    Key: protein_id, Value: ['AAF26470.1']
    Key: translation, Value: ['MERTNSIEMDRERLTAEMAFRDSSSAVIRIRRRLPDLLTSVKLKYVKLGLHNSCNVTTILFFLIILPLTGTVLVQLTGLTFDTFSELWSNQAVQLDTATRLTCLVFLSFVLTLYVANRSKPVYLVDFSCYKPEDERKISVDSFLTMTEENGSFTDDTVQFQQRISNRAGLGDETYLPRGITSTPPKLNMSEARAEAEAVMFGALDSLFEKTGIKPAEVGILIVNCSLFNPTPSLSAMIVNHYKMREDIKSYNLGGMGCSAGLISIDLANNLLKANPNSYAVVVSTENITLNWYFGNDRSMLLCNCIFRMGGAAILLSNRRQDRKKSKYSLVNVVRTHKGSDDKNYNCVYQKEDERGTIGVSLARELMSVAGDALKTNITTLGPMVLPLSEQLMFLISLVKRKMFKLKVKPYIPDFKLAFEHFCIHAGGRAVLDEVQKNLDLKDWHMEPSRMTLHRFGNTSSSSLWYEMAYTEAKGRVKAGDRLWQIAFGSGFKCNSAVWKALRPVSTEEMTGNAWAGSIDQYPVKVVQ']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[63132:63513](-), [61669:61826](-), [59507:59665](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715642']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['hypothetical protein']
    Key: product, Value: ['T25K16.12']
    Key: protein_id, Value: ['AAF26469.1']
    Key: translation, Value: ['MEKRSDSESVEILGDWDSPPPEERIVMVSVPTSPESDYARSNQPKEIESRVSDKETASASGEVAARRVLPPWMDPSYEWGGGKWKVDGRKNKNKKEKEKEKEEIIPFKEIIEALLGNSGDKVQQDNKVFEVAPSLHVVELRKTGDDTLEFHKVYFRFNLYQPVQLPLILFVVIRFSMLKIIHYHQFTMAHIKEFVCMWDTHLYKEITNLNIWDTLSSTLVLAIWTVNASHE']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[67025:67214](-), [66536:66599](-), [66379:66451](-), [66151:66259](-), [65963:66044](-), [65808:65862](-), [65434:65566](-), [65264:65354](-), [65032:65158](-), [64811:64919](-), [64602:64719](-), [64452:64509](-), [64271:64358](-), [64099:64177](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715641']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['similar to wpk4 protein kinase dbj|BAA34675; similar to ESTs dbj|AB015122, gb|AI997157.1']
    Key: product, Value: ['T25K16.13']
    Key: protein_id, Value: ['AAF26468.1']
    Key: translation, Value: ['MSGSRRKATPASRTRVGNYEMGRTLGEGSFAKVKYAKNTVTGDQAAIKILDREKVFRHKMVEQLKREISTMKLIKHPNVVEIIEVMASKTKIYIVLELVNGGELFDKIAQQGRLKEDEARRYFQQLINAVDYCHSRGVYHRDLKPENLILDANGVLKVSDFGLSAFSRQVREDGLLHTACGTPNYVAPEVLSDKGYDGAAADVWSCGVILFVLMAGYLPFDEPNLMTLYKRVRICKAEFSCPPWFSQGAKRVIKRILEPNPITRISIAELLEDEWFKKGYKPPSFDQDDEDITIDDVDAAFSNSKECLVTEKKEKPVSMNAFELISSSSEFSLENLFEKQAQLVKKETRFTSQRSASEIMSKMEETAKPLGFNVRKDNYKIKMKGDKSGRKGQLSVATEVFEVAPSLHVVELRKTGGDTLEFHKVCDSFYKNFSSGLKDVVWNTDAAAEEQKQ']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[71643:71700](-), [70742:71357](-), [70533:70670](-), [69830:69987](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715640']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['similar to ataxia-telangiectasia group D protein pir|A49618']
    Key: product, Value: ['T25K16.14']
    Key: protein_id, Value: ['AAF26467.1']
    Key: translation, Value: ['MVSDLPLDEDDIALLKSPYCDDGGDEDVNSAPNIFTYDNVPLKKRHYLGTSDTFRSFEPLNEHACIVCDIADDGVVPCSGNECPLAVHRKCVELDCEDPATFYCPYCWFKEQATRSTALRTRGVAAAKTLVQYGCSELRSGDIVMTRENSQLENGSDNSLPMQLHENLHQLQELVKHLKARNSQLDESTDQFIDMEKSCGEAYAVVNDQPKRVLWTVNEEKMLREGVEKFSDTINKNMPWKKILEMGKGIFHTTRNSSDLKDKWRNMVRIIILIWLRSRLTSSSSSQRSEIKMERERNAGVMKKMSPTGTIQRLEFVGWYL']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[72284:72371](+), [72788:72865](+), [72988:73097](+), [73189:73442](+), [73523:73585](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715636']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['similar to SYT gi|2252866; similar to ESTs emb|F14390, gb|H36066, emb|F14391']
    Key: product, Value: ['T25K16.15']
    Key: protein_id, Value: ['AAF26463.1']
    Key: translation, Value: ['MQQQQSPQMFPMVPSIPPANNITTEQIQKYLDENKKLIMAIMENQNLGKLAECAQYQALLQKNLMYLAAIADAQPPPPTPGPSPSTAVAAQMATPHSGMQPPSYFMQHPQASPAGIFAPRGPLQFGSPLQFQDPQQQQQIHQQAMQGHMGIRPMGMTNNGMQHAMQQPETGLGGNVGLRGGKQDGADGQGKDDGK']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[74035:74145](-), [73806:73990](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715639']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['similar to stress-induced protein OZI1 precursor pir|S59544; similar to EST gb|AI995719.1']
    Key: product, Value: ['T25K16.16']
    Key: protein_id, Value: ['AAF26466.1']
    Key: translation, Value: ['MASGGKAKYIIGALIGSFGISYIFDKVISDNKIFGGKDDLNGYLLVKISGTTPGTVSNKEWWAATDEKFQAWPRTAGPPVVMNPISRQNFIVKTRPE']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[75334:76249](+), [76515:76653](+), [76732:76982](+), [77014:77148](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715637']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['putative reverse transcriptase gb|AAD17395']
    Key: product, Value: ['T25K16.17']
    Key: protein_id, Value: ['AAF26464.1']
    Key: translation, Value: ['MKEDRRLPHKRDAFQFLKTKAAYVIVIVLTYAFGYFSAYHYHQPLQQQLPPSTTAVETTKPQVCSIDNFRVTTPCGNLVPPELIRQTVIDRIFNGTSPYIDFPPPHAKKFLRPKRIKGWGSYGAVFENLIRRVKPKTIVEVGSFLGASAIHMANLTRRLGLEETQILCVDDFRGWPGFRDRFKDMALVNGDVLLMYQFMQNVVISDFSGSILPVPFSTGSALEKLCEWGVTADLVEIDAGHDFNSAWADINRAVRILRPGGVIFGHDYFTAADNRGVRRAVNLFAEINRLKVKTDGQHWVIDSVKVINKGTRFAISKTVAKIKEDANQWFFAQVLENQDLVNEQAVHISVKVLRGFLRDEHGKVLIHARRSFASVHSKLDATFLCWQWAMESMKSLRVDKIIFASEDNDLIGAVTRLPSWPSYKFQIHFLLGELIRSSNLGAHLIAKSVTMEDRRQSYVATGFPFWLKHLFEKERSIA']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[83585:84581](-), [82750:83373](-), [82722:82738](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6715638']
    Key: evidence, Value: ['not_experimental']
    Key: note, Value: ['putative cytochrome P450 gi|3831440']
    Key: product, Value: ['T25K16.18']
    Key: protein_id, Value: ['AAF26465.1']
    Key: translation, Value: ['MFSLNMRTEIESLWVFALASKFNIYMQQHFASLLVAIAITWFTITIVFWSTPGGPAWGKYFFTRRFISLDYNRKYKNLIPGPRGFPLVGSMSLRSSHVAHQRIASVAEMSNAKRLMAFSLGDTKVVVTCHPAVAKEILNSSVFADRPVDETAYGLMFNRAMGFAPNGTYWRTLRRLGSNHLFNPKQIKQSEDQRRVIATQMVNAFARNPKSACAVRDLLKTASLCNMMGLVFGREYELESNNNLESECLKGLVEEGYDLLGTLNWTDHLPWLAGLDFQQIRFRCSQLVPKVNLLLSRIIHEQRAATGNFLDMLLSLQGSEKLSESDMVAVLWEMIFRGTDTVAVLVEWVLARIVMHPKVQLTVHDELDRVVGRSRTVDESDLPSLTYLTAMIKEVLRLHPPGPLLSWARLSITDTSVDGYHVPAGTTAMVNMWAIARDPHVWEDPLEFKPERFVAKEGEAEFSVFGSDLRLAPFGSGKRVCPGKNLGLTTVSFWVATLLHEFEWLPSVEANPPDLSEVLRLSCEMACPLIVNVSSRRKIIAWMF']
""",
                -1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_06(self):
        path = "GenBank/protein_refseq.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "MNNRWILHAAFLLCFSTTALSINYKQLQLQERTNIRKCQELLEQLNGKINLTYR...FQN"
        id = "NP_034640.1"
        name = "NP_034640"
        description = "interferon beta, fibroblast [Mus musculus]"
        annotations = {
            "accessions": ["NP_034640"],
            "comment": """\
PROVISIONAL REFSEQ: This record has not yet been subject to final
NCBI review. The reference sequence was derived from K00020.1.""",
            "data_file_division": "ROD",
            "date": "01-NOV-2000",
            "db_source": "REFSEQ: accession NM_010510.1",
            "gi": "6754304",
            "keywords": [""],
            "molecule_type": "protein",
            "organism": "Mus musculus",
            "pid": "g6754304",
            "sequence_version": 1,
            "source": "house mouse",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Rodentia",
                "Sciurognathi",
                "Muridae",
                "Murinae",
                "Mus",
            ],
        }
        references = [
            "location: [0:182]\nauthors: Higashi,Y., Sokawa,Y., Watanabe,Y., Kawade,Y.,"
            " Ohno,S., Takaoka,C. and Taniguchi,T.\ntitle: structure and expression of"
            " a cloned cdna for mouse interferon-beta\njournal: J. Biol. Chem. 258,"
            " 9522-9529 (1983)\nmedline id: 83265757\npubmed id: \ncomment: \n"
        ]
        features = (
            (
                """\
type: source
location: [0:182]
qualifiers:
    Key: chromosome, Value: ['4']
    Key: db_xref, Value: ['taxon:10090']
    Key: map, Value: ['4 42.6 cM']
    Key: organism, Value: ['Mus musculus']
""",
                None,
            ),
            (
                """\
type: Protein
location: [0:182]
qualifiers:
    Key: product, Value: ['interferon beta, fibroblast']
""",
                None,
            ),
            (
                """\
type: sig_peptide
location: [0:21]
qualifiers:
""",
                None,
            ),
            (
                """\
type: Region
location: [0:182]
qualifiers:
    Key: db_xref, Value: ['CDD:pfam00143']
    Key: note, Value: ['interferon']
    Key: region_name, Value: ['Interferon alpha/beta domain']
""",
                None,
            ),
            (
                """\
type: mat_peptide
location: [21:182]
qualifiers:
    Key: product, Value: ['ifn-beta']
""",
                None,
            ),
            (
                """\
type: Region
location: [55:170]
qualifiers:
    Key: db_xref, Value: ['CDD:IFabd']
    Key: note, Value: ['IFabd']
    Key: region_name, Value: ['Interferon alpha, beta and delta.']
""",
                None,
            ),
            (
                """\
type: CDS
location: [0:182]
qualifiers:
    Key: coded_by, Value: ['NM_010510.1:21..569']
    Key: db_xref, Value: ['LocusID:15977', 'MGD:MGI:107657']
    Key: gene, Value: ['Ifnb']
""",
                None,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_07(self):
        path = "GenBank/extra_keywords.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "TCCAGGGGATTCACGCGCAATATGTTTCCCTCGCTCGTCTGCAGGGTGTGGGAA...TTG"
        id = "AL138972.1"
        name = "DMBR25B3"
        description = "Drosophila melanogaster BAC clone BACR25B3"
        annotations = {
            "accessions": ["AL138972"],
            "comment": """\
Sequence submitted by Takis Benos, EMBL Outstation - The EBI,
Hinxton, Cambridge, CB10 1SD, U.K.
E-mail: benos@ebi.ac.uk on behalf of the European Drosophila Genome
Sequencing Consortium. For further information see the European
Drosophila Genome Sequencing Consortium's web site:
http://edgp.ebi.ac.uk/.
The syntax for the representation of annotation used in this record
is documented at:
ftp://ftp.ebi.ac.uk/pub/databases/edgp/sequence_annotation.README
Coding sequences are predicted from computer analysis, using both
gene and CDS prediction programs and matches to other sequences.
These predictions and matches have been evaluated by the annotators
and may have been refined by hand (in which case a Genefinder
prediction will have no score. The annotators have also used their
judgement on what matches to represent in this record. A far more
complete annotation record is available from FlyBase
(http://flybase.bio.indiana.edu/) through the FlyBase Annotation
Object linked by the db_xref qualifier in the Feature Table.
IMPORTANT:  This sequence is NOT necessarily the entire insert of
clone BACR25B3.  It may be shorter, since we are minimising the
overlap between clones to 100 bases, by trimming them. Sequence in
absolute orientation with respect to chromosome Clone=BACR25B3;
Contig ID=1; Length=154329; Status=Finished.""",
            "data_file_division": "INV",
            "date": "07-FEB-2000",
            "gi": "6946668",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Drosophila melanogaster",
            "sequence_version": 1,
            "source": "fruit fly",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Arthropoda",
                "Tracheata",
                "Hexapoda",
                "Insecta",
                "Pterygota",
                "Neoptera",
                "Endopterygota",
                "Diptera",
                "Brachycera",
                "Muscomorpha",
                "Ephydroidea",
                "Drosophilidae",
                "Drosophila",
            ],
        }
        references = [
            "location: [0:154329]\nauthors: Murphy,L., Harris,D. and Barrell,B.\ntitle:"
            " Sequencing the distal X chromosome of Drosophila melanogaster\njournal:"
            " Unpublished\nmedline id: \npubmed id: \ncomment: Sanger Centre, Wellcome"
            " Trust Genome Campus, Hinxton Hall, Hinxton, Cambridge CB10 1SA, U.K.\n",
            "location: [0:154329]\nauthors: Benos,P.\ntitle: Direct"
            " Submission\njournal: Submitted (06-FEB-2000) European Drosophila Genome"
            " Sequencing Consortium\nmedline id: \npubmed id: \ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:154329](+)
qualifiers:
    Key: clone, Value: ['BAC BACR25B3']
    Key: db_xref, Value: ['taxon:7227']
    Key: organism, Value: ['Drosophila melanogaster']
""",
                1,
            ),
            (
                """\
type: gene
location: [22147:27773](-)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.11']
    Key: note, Value: ['']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[27622:27773](-), [26676:27009](-), [25023:25178](-), [24615:24888](-), [23629:24555](-), [22859:23560](-), [22374:22791](-), [22147:22299](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946669']
    Key: gene, Value: ['EG:BACR25B3.11']
    Key: note, Value: ["/prediction=(method:''genefinder'', version:''084'', score:''105.71''); /prediction=(method:''genscan'', version:''1.0''); /match=(desc:''BASEMENT MEMBRANE-SPECIFIC HEPARAN SULFATE PROTEOGLYCAN CORE PROTEIN PRECURSOR (HSPG) (PERLECAN) (PLC)'', species:''Homo sapiens (Human)'', ranges:(query:24292..24549, target:SWISS-PROT::P98160:3713..3628, score:''201.00''), (query:24016..24291, target:SWISS-PROT::P98160:3815..3724, score:''139.00''), (query:23857..24006, target:SWISS-PROT::P98160:3866..3817, score:''99.00''), (query:24052..24327, target:SWISS-PROT::P98160:4059..3968, score:''143.00''), (query:24046..24312, target:SWISS-PROT::P98160:4341..4253, score:''116.00''), (query:23806..23901, target:SWISS-PROT::P98160:4177..4146, score:''76.00''), (query:23203..23382, target:SWISS-PROT::P98160:4062..4003, score:''116.00''), (query:22523..22777, target:SWISS-PROT::P98160:4288..4204, score:''112.00''), (query:22235..22300, target:SWISS-PROT::P98160:4358..4337, score:''64.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''GM03359.5prime GM Drosophila melanogaster ovary BlueScript Drosophila melanogaster cDNA clone GM03359 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:25024..25235, target:EMBL::AA801707:438..227, score:''1024.00''), (query:24851..24898, target:EMBL::AA801707:476..429, score:''204.00'')), method:''blastn'', version:''1.4.9''); /match=(desc:''LD08615.5prime LD Drosophila melanogaster embryo BlueScript Drosophila melanogaster cDNA clone LD08615 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:24629..24727, target:EMBL::AA264808:99..1, score:''495.00''), (query:24417..24566, target:EMBL::AA264808:250..101, score:''687.00''), (query:24048..24420, target:EMBL::AA264808:618..246, score:''1847.00''), (query:23986..24036, target:EMBL::AA264808:678..628, score:''237.00'')), method:''blastn'', version:''1.4.9''); /match=(desc:''HL02745.5prime HL Drosophila melanogaster head BlueScript Drosophila melanogaster cDNA clone HL02745 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:23944..24045, target:EMBL::AA697546:103..2, score:''510.00''), (query:23630..23943, target:EMBL::AA697546:416..103, score:''1570.00''), (query:23419..23561, target:EMBL::AA697546:558..416, score:''715.00''), (query:23306..23417, target:EMBL::AA697546:670..559, score:''524.00''), (query:23280..23316, target:EMBL::AA697546:695..659, score:''167.00'')), method:''blastn'', version:''1.4.9''); /match=(desc:''GM08137.5prime GM Drosophila melanogaster ovary BlueScript Drosophila melanogaster cDNA clone GM08137 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:23235..23278, target:EMBL::AA696682:44..1, score:''139.00''), (query:22986..23251, target:EMBL::AA696682:294..29, score:''1321.00'')), method:''blastn'', version:''1.4.9'')"]
    Key: protein_id, Value: ['CAB72284.1']
    Key: translation, Value: ['MACNCNQSMIYQSNERRDYNCPGAPQYPYNRFKGGVSLKDTPCMVLYICADFKSSKLSSAKPIISGPATTRAPAISYVCQPNDFKCVSHPHTCVRANMVCDGIYDCTDHSDEFNCIAGKGSGKSESNSGSGSFKRWKKSPEQGRRSLAKAVKNRKLRKRSFAKSRDYSLKLDDQSSNLRAGESTDVECYSSDDTYTDVVWERSDGAPLSNNVRQVGNRLVISNVSPSDAGNYVCKCKTDEGDLYTTSYKLEVEDQPHELKSSKIVYAKVGANADLQCGADESRQPTYRWSRQYGQLQAGRSLMNEKLSLDSVQANDAGTYICTAQYADGETADFPNILVVTGAIPQFRQEPRSYMSFPTLPNSSFKFNFELTFRPENGDGLLLFNGQTRGSGDYIALSLKDRYAEFRFDFGGKPMLVRAEEPLALNEWHTVRVSRFKRDGYIQVDEQHPVAFPTLQQIPQLDLIEDLYIGGVPNWELLPADAVSQQVGFVGCISRLTLQGRTVELIREAKYKEGITDCRPCAQGPCQNKGVCLESQTEQAYTCICQPGWTGRDCAIEGTQCTPGVCGAGRCENTENDMECLCPLNRSGDRCQYNEILNEHSLNFKGNSFAAYGTPKVTKVNITLSVRPASLEDSVILYTAESTLPSGDYLALVLRGGHAELLINTAARLDPVVVRSAEPLPLNRWTRIEIRRRLGEGILRVGDGPERKAKAPGSDRILSLKTHLYVGGYDRSTVKVNRDVNITKGFDGCISRLYNFQKPVNLLADIKDAANIQSCGETNMIGGDEDSDNEPPVPPPTPDVHENELQPYAMAPCASDPCENGGSCSEQEDVAVCSCPFGFSGKHCQEHLQLGFNASFRGDGYVELNRSHFQPALEQSYTSMGIVFTTNKPNGLLFWWGQEAGEEYTGQDFIAAAVVDGYVEYSMRLDGEEAVIRNSDIRVDNGERHIVIAKRDENTAILEVDRMLHSGETRPTSKKSMKLPGNVFVGGAPDLEVFTGFRYKHNLNGCIVVVEGETVGQINLSSAAVNGVNANVCPA']
""",
                -1,
            ),
            (
                """\
type: gene
location: [29925:33978](-)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.10']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[33816:33978](-), [33532:33713](-), [32685:33289](-), [32323:32634](-), [31658:31836](-), [31196:31591](-), [30616:31076](-), [30269:30519](-), [29925:30108](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946670']
    Key: gene, Value: ['EG:BACR25B3.10']
    Key: note, Value: ["/prediction=(method:''genefinder'', version:''084'', score:''98.50''); /prediction=(method:''genscan'', version:''1.0''); /match=(desc:''BASEMENT MEMBRANE-SPECIFIC HEPARAN SULFATE PROTEOGLYCAN CORE PROTEIN PRECURSOR (HSPG) (PERLECAN) (PLC)'', species:''Homo sapiens (Human)'', ranges:(query:33540..33716, target:SWISS-PROT::P98160:2716..2658, score:''113.00''), (query:32859..32963, target:SWISS-PROT::P98160:3341..3307, score:''63.00''), (query:33150..33215, target:SWISS-PROT::P98160:3530..3509, score:''73.00''), (query:32973..33089, target:SWISS-PROT::P98160:3588..3550, score:''71.00''), (query:32358..32567, target:SWISS-PROT::P98160:3650..3581, score:''107.00''), (query:31222..31323, target:SWISS-PROT::P98160:2620..2587, score:''80.00''), (query:31489..31572, target:SWISS-PROT::P98160:3387..3360, score:''72.00''), (query:31495..31593, target:SWISS-PROT::P98160:3575..3543, score:''60.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''GM02481.5prime GM Drosophila melanogaster ovary BlueScript Drosophila melanogaster cDNA clone GM02481 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:30008..30036, target:EMBL::AA695253:29..1, score:''145.00''), (query:29549..30004, target:EMBL::AA695253:487..32, score:''2262.00'')), method:''blastn'', version:''1.4.9'')"]
    Key: protein_id, Value: ['CAB72285.1']
    Key: translation, Value: ['MFLATLDTNDPTDIGTEDPVLTQIIVSIQKPEITIVPVGGSMTLSCSGRMRWSNSPVIVNWYKENSRLPENVEVQGGNLYLYDLQVSDSGVYICQAVNNETASVFKDTVSITITKKDQLSPAEIVNLPSHVTFEEYVNNEIICEVLGNPAPRVTWARVDGHADAQSTRTYDNRLIFDSPRKSDEGRYRCQAENDQNRDEKYVIVYVQSNPPQPPPQQDRLYITPEEINGLAGESFQLNCQFTSVASLRYDWSHNGRSLSSSPARNVEIRGNTLEVRDASESDSGVYTCVAYDVRTRRNFTESARVNIDRREEQPFGVLMRMMILTDSLINHSNKPIIESLEQNILIIQGEDYSITCEASGSPYPSIKWAKVHDFMPENVHISGNVLTIYGARFENRGVYSCVAENDHGSDLSSTSIDIEPRERPSVKIVSAPLQTFSVGAPASLYCTVEGIPDPTVEWVRVDGQPLSPRHKIQSPGYMVIDDIQLEDSGDYECRAKNIVGEATGVATITVQEPTLVQIIPDNRDLRLTEGDELSLTCVGSGVPNPEVEWVNEMALKRDLYSPPSNTAILKIYRVTKADAGIYTCHGKNEAGSDEAHVRVEVQERRGDIGGVDDDSDRDPINYNPPQQQNPGIHQPGSNQLLATDIGDNVTLTCDMFQPLNTRWERVDGAPLPRNAYTIKNRLEIVRVEQQNLGQYRCNGIGRDGNVKTYFVKELVLMPLPRIRFYPNIPLTVEAGQNLDVHCQVENVRPEDVHWSTDNNRPLPSSVRIVGSVLRFVSITQAAAGEYRCSAFNQYGNRSQIARVAVKKPADFHQVPQSQLQRHREGENIQLQCTVTDQYGVRAQDNVEFNWFRDDRRPLPNNARTDSQILVLTNLRPEDAGRYICNSYDVDRGQQLPEVSIDLQVLSE']
""",
                -1,
            ),
            (
                """\
type: gene
location: [36118:56153](-)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.1']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[56030:56153](-), [52848:52966](-), [50826:50937](-), [50318:50441](-), [49970:50102](-), [49762:49876](-), [49436:49562](-), [48410:48878](-), [47682:47831](-), [47221:47315](-), [46517:46688](-), [45975:46125](-), [45660:45793](-), [45147:45233](-), [44811:44928](-), [44240:44438](-), [43603:43837](-), [42750:42876](-), [42187:42415](-), [41854:42085](-), [41545:41620](-), [40680:40814](-), [40518:40612](-), [40344:40434](-), [39655:40042](-), [37280:39517](-), [36118:37213](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946671']
    Key: gene, Value: ['EG:BACR25B3.1']
    Key: note, Value: ["/prediction=(method:''genscan'', version:''1.0''); /prediction=(method:''genefinder'', version:''084''); /match=(desc:''LOW-DENSITY LIPOPROTEIN RECEPTOR-RELATED PROTEIN PRECURSOR (LRP)'', species:''Caenorhabditis elegans'', ranges:(query:50831..50941, target:SWISS-PROT::Q04833:1221..1185, score:''95.00''), (query:50840..51025, target:SWISS-PROT::Q04833:2865..2804, score:''102.00''), (query:50828..50935, target:SWISS-PROT::Q04833:3788..3753, score:''119.00''), (query:50323..50394, target:SWISS-PROT::Q04833:3706..3683, score:''77.00''), (query:50326..50433, target:SWISS-PROT::Q04833:1263..1228, score:''120.00''), (query:49948..50079, target:SWISS-PROT::Q04833:2917..2874, score:''88.00''), (query:49432..49587, target:SWISS-PROT::Q04833:4085..4034, score:''102.00''), (query:49429..49560, target:SWISS-PROT::Q04833:3915..3872, score:''97.00''), (query:48622..48720, target:SWISS-PROT::Q04833:1302..1270, score:''99.00''), (query:47698..47799, target:SWISS-PROT::Q04833:3996..3963, score:''88.00''), (query:47686..47775, target:SWISS-PROT::Q04833:3835..3806, score:''59.00''), (query:47692..47787, target:SWISS-PROT::Q04833:4041..4010, score:''83.00''), (query:47229..47315, target:SWISS-PROT::Q04833:3742..3714, score:''88.00''), (query:47220..47312, target:SWISS-PROT::Q04833:3829..3799, score:''67.00''), (query:47232..47318, target:SWISS-PROT::Q04833:3866..3838, score:''78.00''), (query:46552..46656, target:SWISS-PROT::Q04833:1344..1310, score:''95.00''), (query:46543..46650, target:SWISS-PROT::Q04833:3951..3916, score:''98.00''), (query:45983..46129, target:SWISS-PROT::Q04833:2870..2822, score:''82.00''), (query:45971..46096, target:SWISS-PROT::Q04833:4089..4048, score:''82.00''), (query:45678..45764, target:SWISS-PROT::Q04833:3666..3638, score:''80.00''), (query:45128..45238, target:SWISS-PROT::Q04833:94..58, score:''100.00''), (query:45158..45268, target:SWISS-PROT::Q04833:3990..3954, score:''80.00''), (query:44263..44379, target:SWISS-PROT::Q04833:85..47, score:''77.00''), (query:44251..44367, target:SWISS-PROT::Q04833:3995..3957, score:''100.00''), (query:43605..43688, target:SWISS-PROT::Q04833:2994..2967, score:''84.00''), (query:42764..42877, target:SWISS-PROT::Q04833:2951..2914, score:''77.00''), (query:42180..42377, target:SWISS-PROT::Q04833:260..195, score:''148.00''), (query:42234..42419, target:SWISS-PROT::Q04833:3199..3138, score:''106.00''), (query:39807..40013, target:SWISS-PROT::Q04833:2901..2833, score:''167.00''), (query:39645..39857, target:SWISS-PROT::Q04833:3138..3068, score:''151.00''), (query:39846..40046, target:SWISS-PROT::Q04833:3241..3175, score:''132.00''), (query:39654..39866, target:SWISS-PROT::Q04833:3913..3843, score:''201.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''LOW-DENSITY LIPOPROTEIN RECEPTOR-RELATED PROTEIN 2 PRECURSOR (MEGALIN) (GLYCOPROTEIN 330)'', species:''Homo sapiens (Human)'', ranges:(query:50834..50935, target:SWISS-PROT::P98164:2733..2700, score:''99.00''), (query:50840..50947, target:SWISS-PROT::P98164:3063..3028, score:''94.00''), (query:50831..50926, target:SWISS-PROT::P98164:3918..3887, score:''102.00''), (query:50326..50433, target:SWISS-PROT::P98164:1222..1187, score:''107.00''), (query:50302..50394, target:SWISS-PROT::P98164:3762..3732, score:''91.00''), (query:49773..49904, target:SWISS-PROT::P98164:2939..2896, score:''90.00''), (query:49438..49578, target:SWISS-PROT::P98164:217..171, score:''116.00''), (query:49429..49545, target:SWISS-PROT::P98164:3796..3758, score:''108.00''), (query:48622..48720, target:SWISS-PROT::P98164:3544..3512, score:''94.00''), (query:48595..48708, target:SWISS-PROT::P98164:3720..3683, score:''86.00''), (query:47701..47814, target:SWISS-PROT::P98164:2817..2780, score:''90.00''), (query:47692..47799, target:SWISS-PROT::P98164:3674..3639, score:''60.00''), (query:47217..47366, target:SWISS-PROT::P98164:3716..3667, score:''96.00''), (query:46543..46647, target:SWISS-PROT::P98164:1101..1067, score:''107.00''), (query:46552..46656, target:SWISS-PROT::P98164:3873..3839, score:''84.00''), (query:45989..46126, target:SWISS-PROT::P98164:3832..3787, score:''98.00''), (query:45149..45274, target:SWISS-PROT::P98164:2775..2734, score:''99.00''), (query:44780..44893, target:SWISS-PROT::P98164:268..231, score:''76.00''), (query:44813..44905, target:SWISS-PROT::P98164:1223..1193, score:''73.00''), (query:44251..44361, target:SWISS-PROT::P98164:3630..3594, score:''119.00''), (query:43602..43700, target:SWISS-PROT::P98164:179..147, score:''97.00''), (query:43674..43781, target:SWISS-PROT::P98164:191..156, score:''90.00''), (query:43584..43685, target:SWISS-PROT::P98164:1107..1074, score:''89.00''), (query:42758..42865, target:SWISS-PROT::P98164:1264..1229, score:''79.00''), (query:42204..42413, target:SWISS-PROT::P98164:2810..2741, score:''136.00''), (query:42189..42377, target:SWISS-PROT::P98164:3027..2965, score:''125.00''), (query:42186..42293, target:SWISS-PROT::P98164:3110..3075, score:''109.00''), (query:42198..42389, target:SWISS-PROT::P98164:3584..3521, score:''137.00''), (query:42309..42422, target:SWISS-PROT::P98164:3793..3756, score:''95.00''), (query:39654..39791, target:SWISS-PROT::P98164:63..18, score:''132.00''), (query:39786..40049, target:SWISS-PROT::P98164:1183..1096, score:''230.00''), (query:39657..39890, target:SWISS-PROT::P98164:3109..3032, score:''200.00''), (query:39780..39983, target:SWISS-PROT::P98164:3756..3689, score:''194.00''), (query:39618..39761, target:SWISS-PROT::P98164:3845..3798, score:''105.00''), (query:39651..39779, target:SWISS-PROT::P98164:3964..3922, score:''128.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''GM06086.5prime GM Drosophila melanogaster ovary BlueScript Drosophila melanogaster cDNA clone GM06086 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:50852..51290, target:EMBL::AA802674:672..234, score:''2195.00'')), method:''blastn'', version:''1.4.9''); /match=(desc:''SD04592.5prime SD Drosophila melanogaster Schneider L2 cell culture pOT2 Drosophila melanogaster cDNA clone SD04592 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:37280..37708, target:EMBL::AI532939:429..1, score:''2136.00''), (query:37097..37217, target:EMBL::AI532939:545..425, score:''569.00'')), method:''blastn'', version:''1.4.9''); /match=(desc:''GH03622.5prime GH Drosophila melanogaster head pOT2 Drosophila melanogaster cDNA clone GH03622 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:36446..37075, target:EMBL::AI063674:1..630, score:''3150.00'')), method:''blastn'', version:''1.4.9''); EST embl|AA802674|AA802674 comes from the 5' UTR"]
    Key: protein_id, Value: ['CAB72286.1']
    Key: translation, Value: ['MLLLQLLLQLLLLGKLLLGKTPPTVFGFRLLFAAFRFPLSLHFPHRMHDHFFVRGDTHSCGWKNSTTFTIRISAIYRYLNQCQANEFRCNNGDCIDARKRCNNVSDCSEGEDENEECPAACSGMEYQCRDGTRCISVSQQCDGHSDCSDGDDEEHCDGIVPKLRYTCPKGKFTCRDLSCISIVHRCDGRADCPNDRSDEEGCPCLYDKWQCDDGTCIAKELLCNGNIDCPEDISDERYCEGGYDSEECRFDEFHCGTGECIPMRQVCDNIYDCNDYSDEVNCVEGEEEDRVGIPIGHQPWRPASKHDDWLHEMDTSEYQVYQPSNVYEKANSQNPCASNQFRCTTSNVCIPLHLRCDGFYHCNDMSDEKSCEQYQRHTTTRRPLTLATPTSRITTQGPGLLERRNTTTATEASRWPWATKTTTIATTTSNPITTVGVANSPPQTCLENIEFACHNRDCISIESVCDGIPDCGRNEDEDDALCKCSGDKYKCQRGGGCIPKSQVCDGKPQCHDRSDESACHLHGRLNKTRLGVKCLESQYQCGDGSCISGYKRCNGIHDCADASDEYNCIYDYEDTYDTDPNNNPLNECDILEFECDYSQCLPLEKKCDGYADCEDMSDELECQSYTDHCLESEFECDSYCLPRDQLCNGIPNCQDGSDERNCTFCREDAYLCNTGECVADNQRCNGIADCADGSDERHCARIYCPPNKLACNGTCVSRRIKCDGIRDCLDGYDEMYCPETNNHYPTQNVNVIRPKLGPNPIPKSCRPHEWQCANLECIDSSLQCNEIKDCSDGSDEELSVCFGTATTRLKPSDCSPEQFYCDESCYNRSVRCNGHVDCSDGSDEVGCSLPCPQHQCPSGRCYTESERCDRHRHCEDGSDEANCTAILCKDNEFLCFDRQFCINATQQCDGYYDCRDFSDEQNCIGCYANQFRCNNGDCVSGSAPCNGYSECSDHSDELNCGGTQECLPNQFRCNSGQCVSSSVRCNGRTDCQDSSDEQNCGHRHTEVSQGLETTGVFTTSTTSTTAMTPLRIICPPTSFKCENGPCISLGLKCNGRVDCPYDGSDEADCGQISNDIDPADSNDRRPNQLNLKTYPDSQIIKESREVIFRCRDEGPARAKVKWSRPGGRPLPPGFTDRNGRLEIPNIRVEDAGTYVCEAVGYASYIPGQQVTVNLNVERSWGENKYEEIRSNRIRYGTVPHIDLEFFGLDNDVGSRPESACTEYQATCMNGECIDKSSICDGNPDCSDASDEQSCSLGLKCQPNQFMCSNSKCVDRTWRCDGENDCGDNSDETSCDPEPSGAPCRYNEFQCRSGHCIPKSFQCDNVPDCTDGTDEVGCMAPLPIRPPPQSVSLLEYEVLELTCVATGTPTPTIVWRLNWGHVPDKCESKSYGGTGTLRCPDMRPQDSGAYSCEIINTRGTHFVNPDTIVTVRPVRTDVCEAGFFNMLARKAEECVQCFCFGVAKACDSANLFTYAIHPPILSHRVVSVELSPLRQIVINEAAPGQDLLTLLHGVQFRATNVHFSGRETPYLALPADYMGNQLKSYGGNLRYEVNYRGSGRPVNGPDVIITGNRFTLTYRVRTQPGQNNRVSIPFVPGGWQKPDGRKASREEIMMILANVDNILIRLGYLDSTAREVDLINIALDSAGTADKGLGSASLVEKCQCPPGYVGDSCESCASGYVRQPGGPWLGHCVPFIPDSCPSGTYGDPRRGVPCKECPCPLTGSNNFASGCQQSPDGDVVCRCNEGYTGRRCEQCAAGYQGNPLAAGGICRRIPDTSCNVDGTYSVHSNGTCQCKDSVIGEQCDTCKSKSFHLNSFTYTGCIECFCSGVGLDCDSSTWYRDQVTSTFGRSRVDHGFVLVTNYMQPTPDTVPVSMAAEPNALSFIGSADQSGNTLYWSLPAAFLGNKLSSYGGKLTYTLSYSPLPNGIMSRNSAPDVVIKSGEDLRLIHYRKSQVVPSVANTYSVEIKESAWQRGDEVVANREHVLMALSDITAIYIKATYTTSTKEASLRQVTLDVATPTNLGTPRAVEVEQCRCPEGYLGLSCEQCAPGYARDPEGGIYLGLCRPCECNGHSKYCNSDTGDCEECSDNTEGPSCERCAAGYVGDATRGTIYDCQPDEGYPIPSPPAPGNQTLECTAYCQIEGIYDCRGNECLCKRNVIGDQCDQCRPGTYGLSAQNQDGCKECYCSGLASQCRSAALYRQLIPVDFILNAPLITDESGAVQDTENLIPDISRNMYTYTHTSYLPKYWSLRGSVLGNQLFSYGGRLSYSLIVESYGNYERGHDIVLIGNGLKLIWSRPDGNENQEEYNVRLHEDEQWTRQDRESARPASRSDFMTVLSDLQHILIRATPRVPTQSTSIGNVILESAVTTRTPGATHASDIELCQCPSGYVGTSCESCAPLHYRDASGSCSLCPCDVSNTESCDLVSGGYVECRCKARWKGDRCREIGE']
""",
                -1,
            ),
            (
                """\
type: gene
location: [70719:75241](-)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.2']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[75216:75241](-), [73085:73559](-), [72838:73016](-), [72604:72768](-), [71423:71621](-), [70719:70988](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946672']
    Key: gene, Value: ['EG:BACR25B3.2']
    Key: note, Value: ["/prediction=(method:''genefinder'', version:''084'', score:''41.82''); /prediction=(method:''genscan'', version:''1.0'')"]
    Key: protein_id, Value: ['CAB72287.1']
    Key: translation, Value: ['MANSKVVAHDESLQGINDSEWQLMGDDIDDGLLDDVDETLKPMETKSEEEDLPTGNWFSQSVHRVRRSINRLFGSDDNQERGRRQQRERSQRNRDAINRQKELRRRQKEDHNRWKQMRMERQLEKQRLVKRTNHVVFNRATDPRKRASDLYDENEASGYHEEDTTLYRTYFVVNEPYDNEYRDRESVQFQNLQKLLDDDLRNFFHSNYEGNDDEEQEIRSTLERVEPTNDNFKIRVQLRIELPTSVNDFGSKLQQQLNVYNRIENLSAATDGVFSFTESSDIEEEAIDVTLPQEEVEGSGSDDSSCRGDATFTCPRSGKTICDEMRCDREIQCPDGEDEEYCNYPNVCTEDQFKCDDKCLELKKRCDGSIDCLDQTDEAGCINAPEPEPEPEPEPEPEPESEPEAEPEPEPEPEPESEPEQEPEPQVPEANGKFY']
""",
                -1,
            ),
            (
                """\
type: gene
location: [121866:127124](+)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.3']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[121866:122046](+), [122173:122630](+), [123671:123823](+), [124062:124320](+), [124391:124688](+), [124754:125018](+), [125093:125254](+), [125316:125576](+), [126792:127124](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946673']
    Key: gene, Value: ['EG:BACR25B3.3']
    Key: note, Value: ["/prediction=(method:''genscan'', version:''1.0'', score:''174.91''); /prediction=(method:''genefinder'', version:''084''); /match=(desc:''PROBABLE G PROTEIN-COUPLED RECEPTOR C13B9.4 IN CHROMOSOME III'', species:''Caenorhabditis elegans'', ranges:(query:123671..123775, target:SWISS-PROT::Q09460:107..141, score:''80.00''), (query:123743..123829, target:SWISS-PROT::Q09460:235..263, score:''72.00''), (query:124072..124332, target:SWISS-PROT::Q09460:265..351, score:''161.00''), (query:124392..124691, target:SWISS-PROT::Q09460:349..448, score:''206.00''), (query:124755..124958, target:SWISS-PROT::Q09460:448..515, score:''123.00''), (query:124764..125027, target:SWISS-PROT::Q09460:454..541, score:''108.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''CALCITONIN RECEPTOR PRECURSOR (CT-R)'', species:''Sus scrofa (Pig)'', ranges:(query:124165..124236, target:SWISS-PROT::P25117:191..214, score:''54.00''), (query:124392..124580, target:SWISS-PROT::P25117:233..295, score:''118.00''), (query:124725..124886, target:SWISS-PROT::P25117:318..371, score:''127.00'')), method:''blastx'', version:''1.4.9'')"]
    Key: protein_id, Value: ['CAB72288.1']
    Key: translation, Value: ['MGAGNRKSETKTKTEAEIEIEMERDQFSIAANACMSMGPMLISKDKAPCSGGRVRHADSLHIYYAVDGKMTLLSNILDCGGCISAQRFTRLLRQSGSSGPSPSAPTAGTFESKSMLEPTSSHSLATGRVPLLHDFDASTTESPGTYVLDGVARVAQLALEPTVMDALPDSDTEQVLGNLNSSAPWNLTLASAAATNFENCSALFVNYTLPQTEFAIRKCELDGRWGSRPNATEVNPPGWTDYGPCYKPEIIRLMQQMGSKDFDAYIDIARRTRTLEIVGLCLSLFALIVSLLIFCTFRSLRNNRTKIHKNLFVAMVLQVIIRLTLYLDQFRRGNKEAATNTSLSVIENTPYLCEASYVLLEYARTAMFMWMFIEGLYLHNMVTVAVFQGSFPLKFFSRLGWCVPILMTTVWARCTVMYMDTSLGECLWNYNLTPYYWILEGPRLAVILLNFCFLVNIIRVLVMKLRQSQASDIEQTRKAVRAAIVLLPLLGITNLLHQLAPLKTATNFAVWSYGTHFLTSFQGFFIALIYCFLNGEVRAVLLKSLATQLSVRGHPEWAPKRASMYSGAYNTAPDTDAVQPAGDPSATGKRISPPNKRLNGRKPSSASIVMIHEPQQRQRLMPRLQNKAREKGKDRVEKTDAEAEPDPTISHIHSKEAGSARSRTRGSKWIMGICFRGQMCDAGLAKDAANIHDVANAADVDACSGSNNNYHNINNNNGSQNNNSIHCNHRDDDKVKGESQSDFKEPSNTNAESLVHLALFTAHTSNTQNNTHRNTIFTPIRRRNCS']
""",
                1,
            ),
            (
                """\
type: gene
location: [128488:129414](-)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.4']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[129373:129414](-), [129195:129313](-), [128776:129140](-), [128488:128715](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946674']
    Key: gene, Value: ['EG:BACR25B3.4']
    Key: note, Value: ["/prediction=(method:''genefinder'', version:''084'', score:''61.35''); /prediction=(method:''genscan'', version:''1.0''); /match=(desc:''VACUOLAR PROTON-ATPASE SUBUNIT D'', species:''Oryctolagus cuniculus (Rabbit)'', ranges:(query:129190..129324, target:SPTREMBL::O97755:55..11, score:''130.00''), (query:128778..129176, target:SPTREMBL::O97755:174..42, score:''472.00''), (query:128546..128716, target:SPTREMBL::O97755:231..175, score:''169.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''VACUOLAR ATP SYNTHASE SUBUNIT D (EC 3.6.1.34) (V-ATPASE D SUBUNIT) (V- ATPASE 28 KD ACCESSORY PROTEIN)'', species:''Bos taurus (Bovine)'', ranges:(query:129190..129324, target:SWISS-PROT::P39942:55..11, score:''130.00''), (query:128778..129176, target:SWISS-PROT::P39942:174..42, score:''471.00''), (query:128546..128716, target:SWISS-PROT::P39942:231..175, score:''173.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''GH28048.5prime GH Drosophila melanogaster head pOT2 Drosophila melanogaster cDNA clone GH28048 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:129196..129317, target:EMBL::AI517334:233..112, score:''412.00''), (query:128777..129145, target:EMBL::AI517334:597..229, score:''1251.00'')), method:''blastn'', version:''1.4.9''); /match=(desc:''GH07112.5prime GH Drosophila melanogaster head pOT2 Drosophila melanogaster cDNA clone GH07112 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:129196..129317, target:EMBL::AI108302:223..102, score:''412.00''), (query:128777..129145, target:EMBL::AI108302:587..219, score:''1251.00''), (query:128636..128716, target:EMBL::AI108302:667..587, score:''243.00'')), method:''blastn'', version:''1.4.9'')"]
    Key: protein_id, Value: ['CAB72289.1']
    Key: translation, Value: ['MAAKDRLPIFPSRGAQTLMKSRLAGATKGHGLLKKKADALQMRFRLILGKIIETKTLMGQVMKEAAFSLAEVKFTTGDINQIVLQNVTKAQIKIRTKKDNVAGVTLPIFEPYTDGVDTYELAGLARGGQQLAKLKKNYQSAVRLLVQLASLQTSFVTLDDVIKVTNRRVNAIEHVIIPRINRTIEYIISELDELEREEFYRLKKIQDKKREARKASDKLRAEQRLLGQMAEAQEVQNILDEDGDEDLLF']
""",
                -1,
            ),
            (
                """\
type: gene
location: [132239:132926](+)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.5']
""",
                1,
            ),
            (
                """\
type: CDS
location: [132239:132926](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946675']
    Key: gene, Value: ['EG:BACR25B3.5']
    Key: note, Value: ["/prediction=(method:''genefinder'', version:''084'', score:''48.06''); /prediction=(method:''genscan'', version:''1.0'', score:''132.90''); /match=(desc:''N-ACETYLTRANSFERASE'', species:''Drosophila melanogaster (Fruit fly)'', ranges:(query:132249..132326, target:SPTREMBL::Q94521:60..85, score:''64.00''), (query:132600..132842, target:SPTREMBL::Q94521:171..251, score:''105.00'')), method:''blastx'', version:''1.4.9''); EST embl|AI063093|AI063093 comes from the 3' UTR"]
    Key: protein_id, Value: ['CAB72290.1']
    Key: translation, Value: ['MEYKMIAPEHSEQVMEHLRRNFFADEPLNKAAGLCQNGSSCPALEAHCAEAIQHRMSVMAVDAKEKDTLKIVGVVLNGILKPGDTAKALSKLDCNDDADFRKIFDLLHRHNLKHNLFEHFDVDCMFDVRILSVDSCYRGQGIANELVKRSVAVAKKNGFRLLKADATGIFSQKIFRSHGFEVFSEQPYSKYTDENGKVILPVEAPHIKLQQLYKAICADDQDEKKQSL']
""",
                1,
            ),
            (
                """\
type: gene
location: [133491:134407](-)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.6']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[134197:134407](-), [133866:134135](-), [133662:133748](-), [133491:133595](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946676']
    Key: gene, Value: ['EG:BACR25B3.6']
    Key: note, Value: ["/prediction=(method:''genscan'', version:''1.0'', score:''119.22''); /prediction=(method:''genefinder'', version:''084''); /match=(desc:''LD41675.5prime LD Drosophila melanogaster embryo pOT2 Drosophila melanogaster cDNA clone LD41675 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:134192..134531, target:EMBL::AI515958:340..1, score:''1691.00''), (query:133879..134139, target:EMBL::AI515958:591..331, score:''1305.00'')), method:''blastn'', version:''1.4.9'')"]
    Key: protein_id, Value: ['CAB72291.1']
    Key: translation, Value: ['MNGLPPSKHYNLTHYQQRYNWDCGLSCIIMILSAQQREQLLGNFDAVCGEEGFGSSTWTIDLCYLLMRYQVRHEYFTQTLGIDPNYAQHTYYSKIIDKDERRVTRKFKDARAHGLRVEQRTVDMEVILRHLARHGPVILLTNASLLTCEVCKRNVLEKFGYAGHYVVLCGYDMAAQKLFYHNPEVHDGHICRCLIESMDTARRAYGTDEDIIFIYEKKETRE']
""",
                -1,
            ),
            (
                """\
type: gene
location: [135478:136829](+)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.7']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[135478:135749](+), [135960:136586](+), [136640:136829](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946677']
    Key: gene, Value: ['EG:BACR25B3.7']
    Key: note, Value: ["/prediction=(method:''genefinder'', version:''084'', score:''66.07''); /prediction=(method:''genscan'', version:''1.0'', score:''145.64''); /match=(desc:''HYPOTHETICAL 40.4 KD TRP-ASP REPEATS CONTAINING PROTEIN C14B1.4 IN CHROMOSOME III'', species:''Caenorhabditis elegans'', ranges:(query:135548..135748, target:SWISS-PROT::Q17963:39..105, score:''120.00''), (query:135957..136586, target:SWISS-PROT::Q17963:105..314, score:''899.00''), (query:136641..136823, target:SWISS-PROT::Q17963:315..375, score:''219.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''LD30385.5prime LD Drosophila melanogaster embryo pOT2 Drosophila melanogaster cDNA clone LD30385 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:135288..135749, target:EMBL::AA950546:102..563, score:''2301.00''), (query:135956..136047, target:EMBL::AA950546:559..650, score:''442.00'')), method:''blastn'', version:''1.4.9''); /match=(desc:''LD10938.5prime LD Drosophila melanogaster embryo BlueScript Drosophila melanogaster cDNA clone LD10938 5prime, mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:136108..136288, target:EMBL::AA392005:776..596, score:''212.00'')), method:''blastn'', version:''1.4.9'')"]
    Key: protein_id, Value: ['CAB72292.1']
    Key: translation, Value: ['MVPIGAVHGGHPGVVHPPQQPLPTAPSGPNSLQPNSVGQPGATTSSNSSASNKSSLSVKPNYTLKFTLAGHTKAVSAVKFSPNGEWLASSSADKLIKIWGAYDGKFEKTISGHKLGISDVAWSSDSRLLVSGSDDKTLKVWELSTGKSLKTLKGHSNYVFCCNFNPQSNLIVSGSFDESVRIWDVRTGKCLKTLPAHSDPVSAVHFNRDGSLIVSSSYDGLCRIWDTASGQCLKTLIDDDNPPVSFVKFSPNGKYILAATLDNTLKLWDYSKGKCLKTYTGHKNEKYCIFANFSVTGGKWIVSGSEDNMVYIWNLQSKEVVQKLQGHTDTVLCTACHPTENIIASAALENDKTIKLWKSDT']
""",
                1,
            ),
            (
                """\
type: gene
location: [145402:147087](+)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.8']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[145402:146203](+), [146514:147087](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946678']
    Key: gene, Value: ['EG:BACR25B3.8']
    Key: protein_id, Value: ['CAB72293.1']
    Key: translation, Value: ['MNSTTKHLLHCTLLITVIVTFEVFSGGIKIDENSFTLVDPWTEYGQLATVLLYLLRFLTLLTLPQVLFNFCGLVFYNAFPEKVVLKGSPLLAPFICIRVVTRGDFPDLVKTNVLRNMNTCLDTGLENFLIEVVTDKAVNLSQHRRIREIVVPKEYKTRTGALFKSRALQYCLEDNVNVLNDSDWIVHLDEETLLTENSVRGIINFVLDGKHPFGQGLITYANENVVNWLTTLADSFRVSDDMGKLRLQFKLFHKPLFSWKGSYVVTQVSAERSVSFDNGIDGSVAEDCFFAMRAFSQGYTFNFIEGEMYEKSPFTLLDFLQQRKRWLQGILLVVHSKMIPFKHKLLLGISVYSWVTMPLSTSNIIFAALYPIPCPNLVDFVCAFIAAINIYMYVFGVIKSFSLYRFGLFRFLACVLGAVCTIPVNVVIENVAVIWGLVGKKHKFYVVQKDVRVLETV']
""",
                1,
            ),
            (
                """\
type: gene
location: [148859:152785](-)
qualifiers:
    Key: gene, Value: ['EG:BACR25B3.9']
""",
                -1,
            ),
            (
                """\
type: CDS
location: join{[152105:152785](-), [151880:152032](-), [149545:151809](-), [148965:149462](-), [148859:148905](-)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946679']
    Key: gene, Value: ['EG:BACR25B3.9']
    Key: note, Value: ["/prediction=(method:''genscan'', version:''1.0''); /prediction=(method:''genefinder'', version:''084''); /match=(desc:''HYPOTHETICAL 135.8 KD PROTEIN'', species:''Drosophila melanogaster (Fruit fly)'', ranges:(query:152096..152785, target:SPTREMBL::Q9XZ29:230..1, score:''1147.00''), (query:151882..152043, target:SPTREMBL::Q9XZ29:277..224, score:''250.00''), (query:149546..151816, target:SPTREMBL::Q9XZ29:1032..276, score:''3735.00''), (query:148953..149465, target:SPTREMBL::Q9XZ29:1202..1032, score:''890.00''), (query:148863..148907, target:SPTREMBL::Q9XZ29:1212..1198, score:''76.00'')), method:''blastx'', version:''1.4.9''); /match=(desc:''LD21815.5prime LD Drosophila melanogaster embryo pOT2 Drosophila melanogaster cDNA clone LD21815 5prime similar to L19117: Drosophila melanogaster (chromosome X 3A6-8) kinesin-like protein of 3A (klp3A) mRNA sequence'', species:''Drosophila melanogaster (fruit fly)'', ranges:(query:152482..152787, target:EMBL::AA816942:460..155, score:''1485.00''), (query:152401..152483, target:EMBL::AA816942:540..458, score:''397.00'')), method:''blastn'', version:''1.4.9'')"]
    Key: protein_id, Value: ['CAB72294.1']
    Key: translation, Value: ['MSSEDPSCVAVALRVRPLVQSELDRGCRIAVERSADGAPQVTVNRNESYTYNYVFDIDDSQKDLFETCVQAKVKKLLNGYNVTILAYGQTGSGKTYTMGTAFNGVLDDHVGVIPRAVHDIFTAIAEMQSEFRFAVTCSFVELYQEQFYDLFSSKTRDKATVDIREVKNRIIMPGLTELVVTSAQQVTDHLIRGSAGRAVAATAMNETSSRSHAIFTLTLVATKLDGKQSVTTSRFNLVDLAGSERCSKTLASGDRFKEGVNINKGLLALGNVINALGSGQAAGYIPYRQSKLTRLLQDSLGGNSITLMIACVSPADYNVAETLSTLRYADRALQIKNKPVVNLDPHAAEVNMLKDVIQKLRVELLSGGKMSSSLISAVGAAGLGAIPCEESLAGSMANAAEIQRLKEQVRTLQDRNRKLQQELHQSLLDLTEKEMRAHIAEQAHDKLRSHVSELKNKLDQREQAQFGNENTNGDNEMRDFSLLVNRVHVELQRTQEELESQGHESRQRLSSRSHTEGGESGGDEVHEMLHSHSEEYTNKQMNFAGELRNINRQLDLKQELHERIMRNFSRLDSDDEDVKLRLCNQKIDDLEAERRDLMDQLRNIKSKDISAKLAEERRKRLQLLEQEISDLRRKLITQANLLKIRDKEREKIQNLSTEIRTMKESKVKLIRAMRGESEKFRQWKMVREKELTQLKSKDRKMQSEIVRQQTLHSKQRQVLKRKCEEALAANKRLKDALERQASAQAQRHKYKDNGGSAAGSSNANAKTDSWVDRELEIILSLIDAEHSLEQLMEDRAVINNHYHLLQQEKTSDPAEAAEQARILASLEEELEMRNAQISDLQQKVCPTDLDSRIRSLAEGVQSLGESRTVSKQLLKTLVQQRRLQASSLNEQRTTLDELRAQLLDAQQQEDAASKRLRLLQSQHEEQMLAQQRAYEEKVSVLIRTANQRWAEARSPAEDQQRNQILEELLSSREALQQELDKLRAKNKSKSKAVKSEPQDLDDSFQIVDGNETVVLSDVSDDPDWVPSTSKSKRIQSDSRNVISPPEKQDANVTSLGNSSIQSLNSTSATEDGKRCKGCKCRTKCTTKRCGCLSGNNACSETCVCKSNCRNPLNLKDHASQCGDGDGQKDETEDADKSDDDGDDEPQTSKENAVKFVTPEAPGKVVASPKQTLQEPKAAATPLMNSNVVEDINGPKLAKMSGLAFDTPKRKFF']
""",
                -1,
            ),
            (
                """\
type: gene
location: join{[153489:154269](+), AL121804.2[40:610](+), AL121804.2[671:1487](+)}
qualifiers:
    Key: gene, Value: ['EG:BACR7C10.3']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[153489:154269](+), AL121804.2[40:610](+), AL121804.2[671:1487](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:6946680']
    Key: gene, Value: ['EG:BACR7C10.3']
    Key: protein_id, Value: ['CAB72295.1']
    Key: translation, Value: ['MEEEAPRFNVLEEAFNGNGNGCANVEATQSAILKVLTRVNRFQMRVRKHIEDNYTEFLPNNTSPDIFLEESGSLNREIHDMLENLGSEGLDALDEANVKMAGNGRQLREILLGLGVSEHVLRIDELFQCVEEAKATKDYLVLLDLVGRLRAFIYGDDSVDGDAQVATPEVRRIFKALECYETIKVKYHVQAYMLQQSLQERFDRLVQLQCKSFPTSRCVTLQVSRDQTQLQDIVQALFQEPYNPARLAEFLLDNCIEPVIMRPVMADYSEEADGGTYVRLSLSYATKEPSSAHVRPNYKQVLENLRLLLHTLAGINCSVSRDQHVFGIIGDHVKDKMLKLLVDECLIPAVPESTEEYQTSTLCEDVAQLEQLLVDSFIINPEQDRALGQFVEKYETYYRNRMYRRVLETAREIIQRDLQDMVLVAPNNHSAEVANDPFLFPRCMISKSAQDFVKLMDRILRQPTDKLGDQEADPIAGVISIMLHTYINEVPKVHRKLLESIPQQAVLFHNNCMFFTHWVAQHANKGIESLAALAKTLQATGQQHFRVQVDYQSSILMGIMQEFEFESTHTLGSGPLKLVRQCLRQLELLKNVWANVLPETVYNATFCELINTFVAELIRRVFTLRHISAQMACELSDLIDVVLQRAPTLFREPNEVVQVLSWLKLQQLKAMLNASLMEITELWGDGVGPLTASYKSDEIKHLIRALFQDTDWRAKAITQIV']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_08(self):
        path = "GenBank/one_of.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "GAATTCAGATAGAATGTAGACAAGAGGGATGGTGAGGAAAACCTACGGCAAGCA...GGC"
        id = "U18266.1"
        name = "HSTMPO1"
        description = "Human thymopoietin (TMPO) gene, exon 1"
        annotations = {
            "accessions": ["U18266"],
            "data_file_division": "PRI",
            "date": "01-JUL-1995",
            "gi": "885676",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Homo sapiens",
            "segment": "1 of 6",
            "sequence_version": 1,
            "source": "human",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        }
        references = [
            "location: [0:2509]\nauthors: Harris,C.A., Andryuk,P.J., Cline,S.W.,"
            " Siekierka,J.J. and Goldstein,G.\ntitle: Structure and mapping of the"
            " human thymopoietin (TMPO) gene and relationship of TMPO beta to rat"
            " lamin-associated polypeptide 2\njournal: Unpublished\nmedline id:"
            " \npubmed id: \ncomment: \n",
            "location: [0:2509]\nauthors: Harris,C.A.\ntitle: Direct"
            " Submission\njournal: Submitted (07-DEC-1994) Crafford A. Harris,"
            " Immunobiology Research Institute, Route 22 East, Annandale, NJ"
            " 08801-0999, USA\nmedline id: \npubmed id: \ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:2509](+)
qualifiers:
    Key: chromosome, Value: ['12']
    Key: clone, Value: ['P1.516 (DMPC-HFFno.1B-0943F)']
    Key: clone_lib, Value: ['DuPont Merck Hum Fibroblast P1 Library no.1 Series B (compressed) (Genome Systems Inc)']
    Key: db_xref, Value: ['taxon:9606']
    Key: map, Value: ['12q22; 64% (% distance from centromere to telomere)']
    Key: organism, Value: ['Homo sapiens']
""",
                1,
            ),
            (
                """\
type: 5'UTR
location: [one-of(1887,1900):2200](+)
qualifiers:
    Key: gene, Value: ['TMPO']
""",
                1,
            ),
            (
                """\
type: gene
location: join{[1887:2509](+), U18267.1[0:270](+), U18268.1[0:309](+), U18270.1[0:6905](+), U18269.1[0:128](+), U18271.1[0:3234](+)}
qualifiers:
    Key: gene, Value: ['TMPO']
""",
                1,
            ),
            (
                """\
type: exon
location: [one-of(1887,1900):2479](+)
qualifiers:
    Key: gene, Value: ['TMPO']
    Key: number, Value: ['1']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[2200:2479](+), U18267.1[119:246](+), U18268.1[129:288](+), U18270.1[4690:4788](+), U18269.1[81:>128](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:885684']
    Key: gene, Value: ['TMPO']
    Key: product, Value: ['thymopoietin beta']
    Key: protein_id, Value: ['AAB60434.1']
    Key: translation, Value: ['MPEFLEDPSVLTKDKLKSELVANNVTLPAGEQRKDVYVQLYLQHLTARNRPPLPAGTNSKGPPDFSSDEEREPTPVLGSGAAAAGRSRAAVGRKATKKTDKPRQEDKDDLDVTELTNEDLLDQLVKYGVNPGPIVGTTRKLYEKKLLKLREQGTESRSSTPLPTISSSAENTRQNGSNDSDRYSDNEEDSKIELKLEKREPLKGRAKTPVTLKQRRVEHNQSYSQAGITETEWTSGS']
""",
                1,
            ),
            (
                """\
type: CDS
location: join{[2200:2479](+), U18267.1[119:246](+), U18268.1[129:288](+), U18270.1[38:1558](+)}
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:885683']
    Key: gene, Value: ['TMPO']
    Key: product, Value: ['thymopoietin alpha']
    Key: protein_id, Value: ['AAB60433.1']
    Key: translation, Value: ['MPEFLEDPSVLTKDKLKSELVANNVTLPAGEQRKDVYVQLYLQHLTARNRPPLPAGTNSKGPPDFSSDEEREPTPVLGSGAAAAGRSRAAVGRKATKKTDKPRQEDKDDLDVTELTNEDLLDQLVKYGVNPGPIVGTTRKLYEKKLLKLREQGTESRSSTPLPTISSSAENTRQNGSNDSDRYSDNEEGKKKEHKKVKSTRDIVPFSELGTTPSGGGFFQGISFPEISTRPPLGSTELQAAKKVHTSKGDLPREPLVATNLPGRGQLQKLASERNLFISCKSSHDRCLEKSSSSSSQPEHSAMLVSTAASPSLIKETTTGYYKDIVENICGREKSGIQPLCPERSHISDQSPLSSKRKALEESESSQLISPPLAQAIRDYVNSLLVQGGVGSLPGTSNSMPPLDVENIQKRIDQSKFQETEFLSPPRKVPRLSEKSVEERDSGSFVAFQNIPGSELMSSFAKTVVSHSLTTLGLEVAKQSQHDKIDASELSFPFHESILKVIEEEWQQVDRQLPSLACKYPVSSREATQILSVPKVDDEILGFISEATPLGGIQAASTESCNQQLDLALCRAYEAAASALQIATHTAFVAKAMQADISEAAQILSSDPSRTHQALGILSKTYDAASYICEAAFDEVKMAAHTMGNATVGRRYLWLKDCKINLASKNKLASTPFKGGTLFGGEVCKVIKKRGNKH']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_09(self):
        path = "GenBank/NT_019265.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = None
        id = "NT_019265.6"
        name = "NT_019265"
        description = "Homo sapiens chromosome 1 working draft sequence segment"
        annotations = {
            "accessions": ["NT_019265"],
            "comment": """\
GENOME ANNOTATION REFSEQ:  NCBI contigs are derived from assembled
genomic sequence data. They may include both draft and finished
sequence.
On Oct 16, 2001 this sequence version replaced gi:15294341.
COMPLETENESS: not full length.""",
            "contig": "join(AL391218.9:105173..108462,gap(100),complement(AL512330.12:1..182490),complement(AL590128.4:9034..81287),gap(100),AL591163.7:85799..94832,gap(100),AL591163.7:94933..113245,gap(100),AL591163.7:42173..44897,complement(AL590128.4:1..6208),AL591163.7:51307..52779,gap(100),AL591163.7:52880..85698,gap(100),AL591163.7:113346..126143,complement(AL159177.12:184729..186047),AL031447.4:1..112158,complement(AL159177.12:1..72671),complement(AL591866.12:23507..86371),AL031848.11:1..142965,AL031847.17:1..166418,AL035406.25:1..161651,complement(AL356261.20:94599..98345),complement(AC026968.3:54432..54579),gap(100),complement(AC062024.2:98529..107911),gap(100),AC062024.2:7713..11594,gap(100),complement(AL356261.20:1..94498),complement(AL356693.18:19988..70853),gap(100),AL356693.18:17351..19887,gap(100),complement(AL356693.18:3037..17250),gap(100),complement(AL356693.18:1..2936),gap(100),AC026968.3:675..2393,gap(100),AC026968.3:1..574,gap(100),AL356261.20:179029..182233)",
            "data_file_division": "CON",
            "date": "16-OCT-2001",
            "gi": "16156830",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Homo sapiens",
            "sequence_version": 6,
            "source": "human",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        }
        references = [
            "location: [0:1250660]\nauthors: NCBI Annotation Project.\ntitle: Direct"
            " Submission\njournal: Submitted (11-OCT-2001) National Center for"
            " Biotechnology Information, NIH, Bethesda, MD 20894, USA\nmedline id:"
            " \npubmed id: \ncomment: \n"
        ]
        features = (
            (
                """\
type: source
location: [0:1250660](+)
qualifiers:
    Key: chromosome, Value: ['1']
    Key: db_xref, Value: ['taxon:9606']
    Key: organism, Value: ['Homo sapiens']
""",
                1,
            ),
            (
                """\
type: source
location: [0:3290](+)
qualifiers:
    Key: clone, Value: ['RP11-13G5']
    Key: db_xref, Value: ['taxon:9606']
    Key: note, Value: ['Accession AL391218 sequenced by The Sanger Centre']
    Key: organism, Value: ['Homo sapiens']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [215901:365470](+)
qualifiers:
    Key: note, Value: ['FISH-mapped clone']
    Key: standard_name, Value: ['RP11-242F24']
""",
                1,
            ),
            (
                """\
type: variation
location: [217507:217508](+)
qualifiers:
    Key: allele, Value: ['T', 'C']
    Key: db_xref, Value: ['dbSNP:811400']
""",
                1,
            ),
            (
                """\
type: mRNA
location: join{[342429:342515](+), [363170:363300](+), [365740:365814](+), [376397:376499](+), [390168:390297](+), [391256:391379](+), [392605:392679](+), [398229:398419](+), [399081:399167](+), [399533:399650](+), [405843:405913](+), [406703:406761](+), [406867:407010](+), [407961:408091](+), [408507:409092](+)}
qualifiers:
    Key: db_xref, Value: ['LocusID:55735']
    Key: gene, Value: ['FLJ10737']
    Key: product, Value: ['hypothetical protein FLJ10737']
    Key: transcript_id, Value: ['XM_057697.1']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_10(self):
        path = "GenBank/origin_line.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "TTAATTAACTGTCTTCGATTGCGTTTAATTGACGGTTTTCGATTAAAAGCGGTA...CGC"
        id = "NC_002678.1"
        name = "NC_002678"
        description = "Mesorhizobium loti, complete genome (edited)"
        annotations = {
            "accessions": ["NC_002678"],
            "data_file_division": "BCT",
            "date": "28-MAR-2001",
            "gi": "13470324",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Mesorhizobium loti",
            "sequence_version": 1,
            "source": "Mesorhizobium loti",
            "taxonomy": [
                "Bacteria",
                "Proteobacteria",
                "alpha subdivision",
                "Rhizobiaceae group",
                "Phyllobacteriaceae",
                "Mesorhizobium",
            ],
            "topology": "circular",
        }
        references = [
            "authors: Kaneko,T., Nakamura,Y., Sato,S., Asamizu,E., Kato,T.,"
            " Sasamoto,S., Watanabe,A., Idesawa,K., Ishikawa,A., Kawashima,K.,"
            " Kimura,T., Kishida,Y., Kiyokawa,C., Kohara,M., Matsumoto,M., Matsuno,A.,"
            " Mochizuki,Y., Nakayama,S., Nakazaki,N., Shimpo,S., Sugimoto,M.,"
            " Takeuchi,C., Yamada,M. and Tabata,S.\ntitle: Complete genome structure of"
            " the nitrogen-fixing symbiotic bacterium Mesorhizobium loti\njournal: DNA"
            " Res. 7, 331-338 (2000)\nmedline id: \npubmed id: \ncomment: \n",
            "location: [0:180]\nauthors: Kaneko,T.\ntitle: Direct Submission\njournal:"
            " Submitted (05-DEC-2000) Takakazu Kaneko, Kazusa DNA Research Institute,"
            " The First Laboratory for Plant Gene Research; Yana 1532-3, Kisarazu,"
            " Chiba 292-0812, Japan (E-mail:kaneko@kazusa.or.jp,"
            " URL:http://www.kazusa.or.jp/rhizobase/, Tel:81-438-52-3935,"
            " Fax:81-438-52-3934)\nmedline id: \npubmed id: \ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:180](+)
qualifiers:
    Key: db_xref, Value: ['taxon:381']
    Key: organism, Value: ['Mesorhizobium loti']
    Key: strain, Value: ['MAFF303099']
""",
                1,
            ),
            (
                """\
type: gene
location: [19:120](+)
qualifiers:
    Key: gene, Value: ['fake']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_11(self):
        path = "GenBank/blank_seq.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "MEECWVTEIANGSKDGLDSNPMKDYMILSGPQKTAVAVLCTLLGLLSALENVAV...SDC"
        id = "NP_001832.1"
        name = "NP_001832"
        description = "cannabinoid receptor 2 (macrophage) [Homo sapiens]"
        annotations = {
            "accessions": ["NP_001832"],
            "comment": """\
REVIEWED REFSEQ: This record has been curated by NCBI staff. The
reference sequence was derived from X74328.1.
Summary: The cannabinoid delta-9-tetrahydrocannabinol is the
principal psychoactive ingredient of marijuana. The proteins
encoded by this gene and the cannabinoid receptor 1 (brain) (CNR1)
gene have the characteristics of a guanine nucleotide-binding
protein (G-protein)-coupled receptor for cannabinoids. They inhibit
adenylate cyclase activity in a dose-dependent, stereoselective,
and pertussis toxin-sensitive manner. These proteins have been
found to be involved in the cannabinoid-induced CNS effects
(including alterations in mood and cognition) experienced by users
of marijuana. The cannabinoid receptors are members of family 1 of
the G-protein-coupled receptors.""",
            "data_file_division": "PRI",
            "date": "18-DEC-2001",
            "db_source": "REFSEQ: accession NM_001841.1",
            "gi": "4502929",
            "keywords": [""],
            "molecule_type": "protein",
            "organism": "Homo sapiens",
            "pid": "g4502929",
            "sequence_version": 1,
            "source": "human",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
            "topology": "linear",
        }
        references = [
            "location: [0:360]\nauthors: Munro,S., Thomas,K.L. and Abu-Shaar,M.\ntitle:"
            " Molecular characterization of a peripheral receptor for"
            " cannabinoids\njournal: Nature 365 (6441), 61-65 (1993)\nmedline id:"
            " 93368659\npubmed id: 7689702\ncomment: \n",
            "location: [0:360]\nauthors: Galiegue,S., Mary,S., Marchand,J.,"
            " Dussossoy,D., Carriere,D., Carayon,P., Bouaboula,M., Shire,D., Le Fur,G."
            " and Casellas,P.\ntitle: Expression of central and peripheral cannabinoid"
            " receptors in human immune tissues and leukocyte subpopulations\njournal:"
            " Eur. J. Biochem. 232 (1), 54-61 (1995)\nmedline id: 96048028\npubmed id:"
            " 7556170\ncomment: \n",
            "location: [0:360]\nauthors: Shire,D., Calandra,B., Rinaldi-Carmona,M.,"
            " Oustric,D., Pessegue,B., Bonnin-Cabanne,O., Le Fur,G., Caput,D. and"
            " Ferrara,P.\ntitle: Molecular cloning, expression and function of the"
            " murine CB2 peripheral cannabinoid receptor\njournal: Biochim. Biophys."
            " Acta 1307 (2), 132-136 (1996)\nmedline id: 96283804\npubmed id:"
            " 8679694\ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:360]
qualifiers:
    Key: chromosome, Value: ['1']
    Key: db_xref, Value: ['taxon:9606']
    Key: map, Value: ['1p36.11']
    Key: organism, Value: ['Homo sapiens']
""",
                None,
            ),
            (
                """\
type: Protein
location: [0:360]
qualifiers:
    Key: product, Value: ['cannabinoid receptor 2 (macrophage)']
""",
                None,
            ),
            (
                """\
type: Region
location: [49:299]
qualifiers:
    Key: db_xref, Value: ['CDD:pfam00001']
    Key: note, Value: ['7tm_1']
    Key: region_name, Value: ['7 transmembrane receptor (rhodopsin family)']
""",
                None,
            ),
            (
                """\
type: CDS
location: [0:360]
qualifiers:
    Key: coded_by, Value: ['NM_001841.1:127..1209']
    Key: db_xref, Value: ['LocusID:1269', 'MIM:605051']
    Key: gene, Value: ['CNR2']
    Key: pseudo, Value: ['']
""",
                None,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_12(self):
        path = "GenBank/dbsource_wrap.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "VKDGYIVDDRNCTYFCGRNAYCNEECTKLKGESGYCQWASPYGNACYCYKVPDH...RCN"
        id = "P01485"
        name = "SCX3_BUTOC"
        description = "Neurotoxin III"
        annotations = {
            "accessions": ["P01485"],
            "comment": """\
[FUNCTION] BINDS TO SODIUM CHANNELS AND INHIBITS THE INACTIVATION
OF THE ACTIVATED CHANNELS, THEREBY BLOCKING NEURONAL TRANSMISSION.
[SUBCELLULAR LOCATION] SECRETED.
[SIMILARITY] BELONGS TO THE ALPHA/BETA-SCORPION TOXIN FAMILY.
ALPHA-TOXIN SUBFAMILY.""",
            "data_file_division": "INV",
            "date": "16-OCT-2001",
            "db_source": (
                "swissprot: locus SCX3_BUTOC, accession P01485; class: standard."
                " created: Jul 21, 1986. sequence updated: Jul 21, 1986. annotation"
                " updated: Oct 16, 2001. xrefs: gi: gi: 69530 xrefs (non-sequence"
                " databases): HSSP P01484, InterPro IPR003614, InterPro IPR002061,"
                " InterPro IPR001219, Pfam PF00537, PRINTS PR00284, ProDom PD000908,"
                " SMART SM00505"
            ),
            "gi": "134354",
            "keywords": ["Neurotoxin", "Sodium channel inhibitor", "Amidation"],
            "molecule_type": "protein",
            "organism": "Buthus occitanus tunetanus",
            "pid": "g134354",
            "source": "Tunisian scorpion",
            "taxonomy": [
                "Eukaryota",
                "Metazoa",
                "Arthropoda",
                "Chelicerata",
                "Arachnida",
                "Scorpiones",
                "Buthoidea",
                "Buthidae",
                "Buthus",
            ],
            "topology": "linear",
        }
        references = [
            "location: [0:64]\nauthors: Vargas,O., Gregoire,J., Martin,M.-F., Bechis,G."
            " and Rochat,H.\ntitle: Neurotoxins from the venoms of two scorpions:"
            " Buthus occitanus tunetanus and Buthus occitanus mardochei\njournal:"
            " Toxicon 20, 79-79 (1982)\nmedline id: \npubmed id: \ncomment: SEQUENCE.\n"
        ]
        features = (
            (
                """\
type: source
location: [0:64]
qualifiers:
    Key: db_xref, Value: ['taxon:6871']
    Key: organism, Value: ['Buthus occitanus tunetanus']
""",
                None,
            ),
            (
                """\
type: Protein
location: [0:64]
qualifiers:
    Key: product, Value: ['Neurotoxin III']
""",
                None,
            ),
            (
                """\
type: Bond
location: bond{[11:12], [62:63]}
qualifiers:
    Key: bond_type, Value: ['disulfide']
    Key: note, Value: ['BY SIMILARITY.']
""",
                None,
            ),
            (
                """\
type: Bond
location: bond{[15:16], [35:36]}
qualifiers:
    Key: bond_type, Value: ['disulfide']
    Key: note, Value: ['BY SIMILARITY.']
""",
                None,
            ),
            (
                """\
type: Bond
location: bond{[21:22], [45:46]}
qualifiers:
    Key: bond_type, Value: ['disulfide']
    Key: note, Value: ['BY SIMILARITY.']
""",
                None,
            ),
            (
                """\
type: Bond
location: bond{[25:26], [47:48]}
qualifiers:
    Key: bond_type, Value: ['disulfide']
    Key: note, Value: ['BY SIMILARITY.']
""",
                None,
            ),
            (
                """\
type: Site
location: [63:64]
qualifiers:
    Key: site_type, Value: ['amidation']
""",
                None,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_13(self):
        path = "GenBank/gbvrl1_start.seq"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
            seq = "ATGTCTGGCAACCAGTATACTGAGGAAGTTATGGAGGGAGTAAATTGGTTAAAG...TAA"
            id = "AB000048.1"
            name = "AB000048"
            description = (
                "Feline panleukopenia virus DNA for nonstructural protein 1,"
                " complete cds"
            )
            annotations = {
                "accessions": ["AB000048"],
                "data_file_division": "VRL",
                "date": "05-FEB-1999",
                "gi": "1769753",
                "keywords": ["nonstructural protein 1"],
                "molecule_type": "DNA",
                "organism": "Feline panleukopenia virus",
                "sequence_version": 1,
                "source": "Feline panleukopenia virus",
                "taxonomy": [
                    "Viruses",
                    "ssDNA viruses",
                    "Parvoviridae",
                    "Parvovirinae",
                    "Parvovirus",
                ],
                "topology": "linear",
            }
            references = [
                "location: [0:2007]\nauthors: Horiuchi,M.\ntitle: Evolutionary pattern"
                " of feline panleukopenia virus differs from that of canine"
                " parvovirus\njournal: Unpublished\nmedline id: \npubmed id:"
                " \ncomment: \n",
                "location: [0:2007]\nauthors: Horiuchi,M.\ntitle: Direct"
                " Submission\njournal: Submitted (22-DEC-1996) Motohiro Horiuchi,"
                " Obihiro University of Agriculture and Veterinary Medicine, Veterinary"
                " Public Health; Inada cho, Obihiro, Hokkaido 080, Japan"
                " (E-mail:horiuchi@obihiro.ac.jp, Tel:0155-49-5392,"
                " Fax:0155-49-5402)\nmedline id: \npubmed id: \ncomment: \n",
            ]
            features = (
                (
                    """\
type: source
location: [0:2007](+)
qualifiers:
    Key: db_xref, Value: ['taxon:10786']
    Key: isolate, Value: ['483']
    Key: lab_host, Value: ['Felis domesticus']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Feline panleukopenia virus']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: [0:2007](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:1769754']
    Key: product, Value: ['nonstructural protein 1']
    Key: protein_id, Value: ['BAA19009.1']
    Key: translation, Value: ['MSGNQYTEEVMEGVNWLKKHAEDEAFSFVFKCDNVQLNGKDVRWNNYTKPIQNEELTSLIRGAQTAMDQTEEEEMDWESEVDSLAKKQVQTFDALIKKCLFEVFVSKNIEPNECVWFIQHEWGKDQGWHCHVLLHSKNLQQATGKWLRRQMNMYWSRWLVTLCSINLTPTEKIKLREIAEDSEWVTILTYRHKQTKKDYVKMVHFGNMIAYYFLTKKKIVHMTKESGYFLSTDSGWKFNFMKYQDRHTVSTLYTEQMKPETVETTVTTAQETKRGRIQTKKEVSIKCTLRDLVSKRVTSPEDWMMLQPDSYIEMMAQPGGENLLKNTLEICTLTLARTKTAFELILEKADNTKLTNFDLANSRTCQIFRMHGWNWIKVCHAIACVLNRQGGKRNTVLFHGPASTGKSIIAQAIAQAVGNVGCYNAANVNFPFNDCTNKNLIWVEEAGNFGQQVNQFKAICSGQTIRIDQKGKGSKQIEPTPVIMTTNENITIVRIGCEERPEHTQPIRDRMLNIKLVCKLPGDFGLVDKEEWPLICAWLVKHGYQSTMANYTHHWGKVPEWDENWAEPKIQEGINSPGCKDLETQAASNPQSQDHVLTPLTPDVVDLALEPWSTPDTPIAETANQQSNQLGVTHKDVQASPTWSEIEADLRAIFTSEQLEEDFRDDLD']
""",
                    1,
                ),
            )
            dbxrefs = []
            self.maxDiff = None
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )
            record = next(records)
            seq = "ATGTCTGGCAACCAGTATACTGAGGAAGTTATGGAGGGAGTAAATTGGTTAAAG...TAA"
            id = "AB000049.1"
            name = "AB000049"
            description = (
                "Feline panleukopenia virus DNA for nonstructural protein 1,"
                " complete cds"
            )
            annotations = {
                "accessions": ["AB000049"],
                "data_file_division": "VRL",
                "date": "05-FEB-1999",
                "gi": "1769755",
                "keywords": ["nonstructural protein 1"],
                "molecule_type": "DNA",
                "organism": "Feline panleukopenia virus",
                "sequence_version": 1,
                "source": "Feline panleukopenia virus",
                "taxonomy": [
                    "Viruses",
                    "ssDNA viruses",
                    "Parvoviridae",
                    "Parvovirinae",
                    "Parvovirus",
                ],
                "topology": "linear",
            }
            references = [
                "location: [0:2007]\nauthors: Horiuchi,M.\ntitle: Evolutionary pattern"
                " of feline panleukopenia virus differs that of canine"
                " parvovirus\njournal: Unpublished\nmedline id: \npubmed id:"
                " \ncomment: \n",
                "location: [0:2007]\nauthors: Horiuchi,M.\ntitle: Direct"
                " Submission\njournal: Submitted (22-DEC-1996) Motohiro Horiuchi,"
                " Obihiro University of Agriculture and Veterinary Medicine, Veterinary"
                " Public Health; Inada cho, Obihiro, Hokkaido 080, Japan"
                " (E-mail:horiuchi@obihiro.ac.jp, Tel:0155-49-5392,"
                " Fax:0155-49-5402)\nmedline id: \npubmed id: \ncomment: \n",
            ]
            features = (
                (
                    """\
type: source
location: [0:2007](+)
qualifiers:
    Key: db_xref, Value: ['taxon:10786']
    Key: isolate, Value: ['94-1']
    Key: lab_host, Value: ['Felis domesticus']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Feline panleukopenia virus']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: [0:2007](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:1769756']
    Key: product, Value: ['nonstructural protein 1']
    Key: protein_id, Value: ['BAA19010.1']
    Key: translation, Value: ['MSGNQYTEEVMEGVNWLKKHAEDEAFSFVFKCDNVQLNGKDVRWNNYTKPIQNEELTSLIRGAQTAMDQTEEEEMDWESEVDSLAKKQVQTFDALIKKCLFEVFVSKNIEPNECVWFIQHEWGKDQGWHCHVLLHSKNLQQATGKWLRRQMNMYWSRWLVTLCSINLTPTEKIKLREIAEDSEWVTILTYRHKQTKKDYVKMVHFGNMIAYYFLTKKKIVHMTKESGYFLSTDSGWKFNFMKYQDRHTVSTLYTEQMKPETVETTVTTAQETKRGRIQTKKEVSIKCTLRDLVSKRVTSPEDWMMLQPDSYIEMMAQPGGENLLKNTLEICTLTLARTKTAFELILEKADNTKLTNFDLANSRTCQIFRMHGWNWIKVCHAIACVLNRQGGKRNTVLFHGPASTGKSIIAQAIAQAVGNVGCYNAANVNFPFNDCTNKNLIWVEEAGNFGQQVNQFKAICSGQTIRIDQKGKGSKQIEPTPVIMTTNENITIVRIGCEERPEHTQPIRDRMLNIKLVCKLPGDFGLVDKEEWPLICAWLVKHGYQSTMANYTHHWGKVPEWDENWAEPKIQEGINSPGCKDLETQAASNPQSQDHVLTPLTPDVVDLALEPWSTPDTPIAETANQQSNQLGVTHKDVQASPTWSEIEADLRAIFTSEQLEEDFRDDLD']
""",
                    1,
                ),
            )
            dbxrefs = []
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )
            record = next(records)
            seq = "ATGAGTGATGGAGCAGTTCAACCAGACGGTGGTCAACCTGCTGTCAGAAATGAA...TAA"
            id = "AB000050.1"
            name = "AB000050"
            description = (
                "Feline panleukopenia virus DNA for capsid protein 2, complete cds"
            )
            annotations = {
                "accessions": ["AB000050"],
                "data_file_division": "VRL",
                "date": "05-FEB-1999",
                "gi": "1769757",
                "keywords": ["capsid protein 2"],
                "molecule_type": "DNA",
                "organism": "Feline panleukopenia virus",
                "sequence_version": 1,
                "source": "Feline panleukopenia virus",
                "taxonomy": [
                    "Viruses",
                    "ssDNA viruses",
                    "Parvoviridae",
                    "Parvovirinae",
                    "Parvovirus",
                ],
                "topology": "linear",
            }
            references = [
                "location: [0:1755]\nauthors: Horiuchi,M.\ntitle: Evolutionary pattern"
                " of feline panleukopenia virus differs from that of canine"
                " parvovirus\njournal: Unpublished\nmedline id: \npubmed id:"
                " \ncomment: \n",
                "location: [0:1755]\nauthors: Horiuchi,M.\ntitle: Direct"
                " Submission\njournal: Submitted (22-DEC-1996) Motohiro Horiuchi,"
                " Obihiro University of Agriculture and Veterinary Medicine, Veterinary"
                " Public Health; Inada cho, Obihiro, Hokkaido 080, Japan"
                " (E-mail:horiuchi@obihiro.ac.jp, Tel:0155-49-5392,"
                " Fax:0155-49-5402)\nmedline id: \npubmed id: \ncomment: \n",
            ]
            features = (
                (
                    """\
type: source
location: [0:1755](+)
qualifiers:
    Key: db_xref, Value: ['taxon:10786']
    Key: isolate, Value: ['94-1']
    Key: lab_host, Value: ['Felis domesticus']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Feline panleukopenia virus']
""",
                    1,
                ),
                (
                    """\
type: CDS
location: [0:1755](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:1769758']
    Key: product, Value: ['capsid protein 2']
    Key: protein_id, Value: ['BAA19011.1']
    Key: translation, Value: ['MSDGAVQPDGGQPAVRNERATGSGNGSGGGGGGGSGGVGISTGTFNNQTEFKFLENGWVEITANSSRLVHLNMPESENYKRVVVNNMDKTAVKGNMALDDTHVQIVTPWSLVDANAWGVWFNPGDWQLIVNTMSELHLVSFEQEIFNVVLKTVSESATQPPTKVYNNDLTASLMVALDSNNTMPFTPAAMRSETLGFYPWKPTIPTPWRYYFQWDRTLIPSHTGTSGTPTNVYHGTDPDDVQFYTIENSVPVHLLRTGDEFATGTFFFDCKPCRLTHTWQTNRALGLPPFLNSLPQSEGATNFGDIGVQQDKRRGVTQMGNTDYITEATIMRPAEVGYSAPYYSFEASTQGPFKTPIAAGRGGAQTDENQAADGDPRYAFGRQHGQKTTTTGETPERFTYIAHQDTGRYPEGDWIQNINFNLPVTNDNVLLPTDPIGGKTGINYTNIFNTYGPLTALNNVPPVYPNGQIWDKEFDTDLKPRLHVNAPFVCQNNCPGQLFVKVAPNLTNEYDPDASANMSRIVTYSDFWWKGKLVFKAKLRASHTWNPIQQMSINVDNQFNYVPNNIGAMKIVYEKSQLAPRKLY']
""",
                    1,
                ),
            )
            dbxrefs = []
            self.perform_feature_parser_test(
                record,
                seq,
                id,
                name,
                description,
                annotations,
                references,
                features,
                dbxrefs,
            )

    def test_feature_parser_14(self):
        path = "GenBank/NC_005816.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Premature end of file in sequence data
                record = next(records)
        seq = "TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG"
        id = "NC_005816.1"
        name = "NC_005816"
        description = (
            "Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete"
            " sequence"
        )
        annotations = {
            "accessions": ["NC_005816"],
            "comment": """\
PROVISIONAL REFSEQ: This record has not yet been subject to final
NCBI review. The reference sequence was derived from AE017046.
COMPLETENESS: full length.""",
            "data_file_division": "BCT",
            "date": "21-JUL-2008",
            "gi": "45478711",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Yersinia pestis biovar Microtus str. 91001",
            "sequence_version": 1,
            "source": "Yersinia pestis biovar Microtus str. 91001",
            "taxonomy": [
                "Bacteria",
                "Proteobacteria",
                "Gammaproteobacteria",
                "Enterobacteriales",
                "Enterobacteriaceae",
                "Yersinia",
            ],
            "topology": "circular",
        }
        references = [
            "location: [0:9609]\nauthors: Zhou,D., Tong,Z., Song,Y., Han,Y., Pei,D.,"
            " Pang,X., Zhai,J., Li,M., Cui,B., Qi,Z., Jin,L., Dai,R., Du,Z., Wang,J.,"
            " Guo,Z., Wang,J., Huang,P. and Yang,R.\ntitle: Genetics of metabolic"
            " variations between Yersinia pestis biovars and the proposal of a new"
            " biovar, microtus\njournal: J. Bacteriol. 186 (15), 5147-5152"
            " (2004)\nmedline id: \npubmed id: 15262951\ncomment: \n",
            "location: [0:9609]\nauthors: Song,Y., Tong,Z., Wang,J., Wang,L., Guo,Z.,"
            " Han,Y., Zhang,J., Pei,D., Zhou,D., Qin,H., Pang,X., Han,Y., Zhai,J.,"
            " Li,M., Cui,B., Qi,Z., Jin,L., Dai,R., Chen,F., Li,S., Ye,C., Du,Z.,"
            " Lin,W., Wang,J., Yu,J., Yang,H., Wang,J., Huang,P. and Yang,R.\ntitle:"
            " Complete genome sequence of Yersinia pestis strain 91001, an isolate"
            " avirulent to humans\njournal: DNA Res. 11 (3), 179-197 (2004)\nmedline"
            " id: \npubmed id: 15368893\ncomment: \n",
            "location: [0:9609]\nauthors: \nconsrtm: NCBI Genome Project\ntitle: Direct"
            " Submission\njournal: Submitted (16-MAR-2004) National Center for"
            " Biotechnology Information, NIH, Bethesda, MD 20894, USA\nmedline id:"
            " \npubmed id: \ncomment: \n",
            "location: [0:9609]\nauthors: Song,Y., Tong,Z., Wang,L., Han,Y., Zhang,J.,"
            " Pei,D., Wang,J., Zhou,D., Han,Y., Pang,X., Zhai,J., Chen,F., Qin,H.,"
            " Wang,J., Li,S., Guo,Z., Ye,C., Du,Z., Lin,W., Wang,J., Yu,J., Yang,H.,"
            " Wang,J., Huang,P. and Yang,R.\ntitle: Direct Submission\njournal:"
            " Submitted (24-APR-2003) The Institute of Microbiology and Epidemiology,"
            " Academy of Military Medical Sciences, No. 20, Dongdajie Street, Fengtai"
            " District, Beijing 100071, People's Republic of China\nmedline id:"
            " \npubmed id: \ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:9609](+)
qualifiers:
    Key: biovar, Value: ['Microtus']
    Key: db_xref, Value: ['taxon:229193']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Yersinia pestis biovar Microtus str. 91001']
    Key: plasmid, Value: ['pPCP1']
    Key: strain, Value: ['91001']
""",
                1,
            ),
            (
                """\
type: repeat_region
location: [0:1954](+)
qualifiers:
""",
                1,
            ),
            (
                """\
type: gene
location: [86:1109](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767718']
    Key: locus_tag, Value: ['YP_pPCP01']
""",
                1,
            ),
            (
                """\
type: CDS
location: [86:1109](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478712', 'GeneID:2767718']
    Key: locus_tag, Value: ['YP_pPCP01']
    Key: note, Value: ['similar to corresponding CDS from previously sequenced pPCP plasmid of Yersinia pestis KIM (AF053945) and CO92 (AL109969), also many transposase entries for insertion sequence IS100 of Yersinia pestis. Contains IS21-like element transposase, HTH domain (Interpro|IPR007101)']
    Key: product, Value: ['putative transposase']
    Key: protein_id, Value: ['NP_995567.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MVTFETVMEIKILHKQGMSSRAIARELGISRNTVKRYLQAKSEPPKYTPRPAVASLLDEYRDYIRQRIADAHPYKIPATVIAREIRDQGYRGGMTILRAFIRSLSVPQEQEPAVRFETEPGRQMQVDWGTMRNGRSPLHVFVAVLGYSRMLYIEFTDNMRYDTLETCHRNAFRFFGGVPREVLYDNMKTVVLQRDAYQTGQHRFHPSLWQFGKEMGFSPRLCRPFRAQTKGKVERMVQYTRNSFYIPLMTRLRPMGITVDVETANRHGLRWLHDVANQRKHETIQARPCDRWLEEQQSMLALPPEKKEYDVHLDENLVNFDKHPLHHPLSIYDSFCRGVA']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [86:959](+)
qualifiers:
    Key: db_xref, Value: ['CDD:34222']
    Key: locus_tag, Value: ['YP_pPCP01']
    Key: note, Value: ['Transposase and inactivated derivatives [DNA replication, recombination, and repair]; Region: COG4584']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [<110:209](+)
qualifiers:
    Key: db_xref, Value: ['CDD:186341']
    Key: locus_tag, Value: ['YP_pPCP01']
    Key: note, Value: ['Helix-turn-helix domain of Hin and related proteins, a family of DNA-binding domains unique to bacteria and represented by the Hin protein of Salmonella. The basic HTH domain is a simple fold comprised of three core helices that form a right-handed...; Region: HTH_Hin_like; cl01116']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [437:812](+)
qualifiers:
    Key: db_xref, Value: ['CDD:194099']
    Key: locus_tag, Value: ['YP_pPCP01']
    Key: note, Value: ['Integrase core domain; Region: rve; cl01316']
""",
                1,
            ),
            (
                """\
type: gene
location: [1105:1888](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767716']
    Key: locus_tag, Value: ['YP_pPCP02']
""",
                1,
            ),
            (
                """\
type: CDS
location: [1105:1888](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478713', 'GeneID:2767716']
    Key: locus_tag, Value: ['YP_pPCP02']
    Key: note, Value: ['similar to corresponding CDS form previously sequenced pPCP plasmid of Yersinia pestis KIM (AF053945) and CO92 (AL109969), also many ATP-binding protein entries for insertion sequence IS100 of Yersinia pestis. Contains Chaperonin clpA/B (Interpro|IPR001270). Contains ATP/GTP-binding site motif A (P-loop) (Interpro|IPR001687, Molecular Function: ATP binding (GO:0005524)). Contains Bacterial chromosomal replication initiator protein, DnaA (Interpro|IPR001957, Molecular Function: DNA binding (GO:0003677), Molecular Function: DNA replication origin binding (GO:0003688), Molecular Function: ATP binding (GO:0005524), Biological Process: DNA replication initiation (GO:0006270), Biological Process: regulation of DNA replication (GO:0006275)). Contains AAA ATPase (Interpro|IPR003593, Molecular Function: nucleotide binding (GO:0000166))']
    Key: product, Value: ['transposase/IS protein']
    Key: protein_id, Value: ['NP_995568.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MMMELQHQRLMALAGQLQLESLISAAPALSQQAVDQEWSYMDFLEHLLHEEKLARHQRKQAMYTRMAAFPAVKTFEEYDFTFATGAPQKQLQSLRSLSFIERNENIVLLGPSGVGKTHLAIAMGYEAVRAGIKVRFTTAADLLLQLSTAQRQGRYKTTLQRGVMAPRLLIIDEIGYLPFSQEEAKLFFQVIAKRYEKSAMILTSNLPFGQWDQTFAGDAALTSAMLDRILHHSHVVQIKGESYRLRQKRKAGVIAEANPE']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [1108:1885](+)
qualifiers:
    Key: db_xref, Value: ['CDD:181681']
    Key: locus_tag, Value: ['YP_pPCP02']
    Key: note, Value: ['transposase/IS protein; Provisional; Region: PRK09183']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [1366:>1669](+)
qualifiers:
    Key: db_xref, Value: ['CDD:99707']
    Key: locus_tag, Value: ['YP_pPCP02']
    Key: note, Value: ['The AAA+ (ATPases Associated with a wide variety of cellular Activities) superfamily represents an ancient group of ATPases belonging to the ASCE (for additional strand, catalytic E) division of the P-loop NTPase fold. The ASCE division also includes...; Region: AAA; cd00009']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [1432:1456](+)
qualifiers:
    Key: db_xref, Value: ['CDD:99707']
    Key: locus_tag, Value: ['YP_pPCP02']
    Key: note, Value: ['Walker A motif; other site']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: order{[1435:1459](+), [1618:1621](+)}
qualifiers:
    Key: db_xref, Value: ['CDD:99707']
    Key: locus_tag, Value: ['YP_pPCP02']
    Key: note, Value: ['ATP binding site [chemical binding]; other site']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [1606:1624](+)
qualifiers:
    Key: db_xref, Value: ['CDD:99707']
    Key: locus_tag, Value: ['YP_pPCP02']
    Key: note, Value: ['Walker B motif; other site']
""",
                1,
            ),
            (
                """\
type: gene
location: [2924:3119](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767717']
    Key: gene, Value: ['rop']
    Key: gene_synonym, Value: ['rom']
    Key: locus_tag, Value: ['YP_pPCP03']
""",
                1,
            ),
            (
                """\
type: CDS
location: [2924:3119](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478714', 'GeneID:2767717']
    Key: gene, Value: ['rop']
    Key: gene_synonym, Value: ['rom']
    Key: locus_tag, Value: ['YP_pPCP03']
    Key: note, Value: ['Best Blastp hit =gi|16082682|ref|NP_395229.1| (NC_003132) putative replication regulatory protein [Yersinia pestis], gi|5763813|emb|CAB531 66.1| (AL109969) putative replication regulatory protein [Yersinia pestis]; similar to gb|AAK91579.1| (AY048853), RNAI modulator protein Rom [Salmonella choleraesuis], Contains Regulatory protein Rop (Interpro|IPR000769)']
    Key: product, Value: ['putative replication regulatory protein']
    Key: protein_id, Value: ['NP_995569.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MNKQQQTALNMARFIRSQSLILLEKLDALDADEQAAMCERLHELAEELQNSIQARFEAESETGT']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [2924:3107](+)
qualifiers:
    Key: db_xref, Value: ['CDD:145136']
    Key: gene, Value: ['rop']
    Key: gene_synonym, Value: ['rom']
    Key: locus_tag, Value: ['YP_pPCP03']
    Key: note, Value: ['Rop protein; Region: Rop; pfam01815']
""",
                1,
            ),
            (
                """\
type: gene
location: [3485:3857](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767720']
    Key: locus_tag, Value: ['YP_pPCP04']
""",
                1,
            ),
            (
                """\
type: CDS
location: [3485:3857](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478715', 'GeneID:2767720']
    Key: locus_tag, Value: ['YP_pPCP04']
    Key: note, Value: ['Best Blastp hit = gi|321919|pir||JQ1541 hypothetical 16.9K protein - Salmonella typhi murium plasmid NTP16.']
    Key: product, Value: ['hypothetical protein']
    Key: protein_id, Value: ['NP_995570.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MSKKRRPQKRPRRRRFFHRLRPPDEHHKNRRSSQRWRNPTGLKDTRRFPPEAPSCALLFRPCRLPDTSPPFSLREAWRFLIAHAVGISVRCRSFAPSWAVCTNPPFSPTTAPYPVTIVLSPTR']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [3497:3626](+)
qualifiers:
    Key: locus_tag, Value: ['YP_pPCP04']
    Key: note, Value: ['ProfileScan match to entry PS50323 ARG_RICH, E-value 8.981']
""",
                1,
            ),
            (
                """\
type: gene
location: [4342:4780](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767712']
    Key: gene, Value: ['pim']
    Key: locus_tag, Value: ['YP_pPCP05']
""",
                1,
            ),
            (
                """\
type: CDS
location: [4342:4780](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
    Key: gene, Value: ['pim']
    Key: locus_tag, Value: ['YP_pPCP05']
    Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
    Key: product, Value: ['pesticin immunity protein']
    Key: protein_id, Value: ['NP_995571.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
""",
                1,
            ),
            (
                """\
type: gene
location: [4814:5888](-)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767721']
    Key: gene, Value: ['pst']
    Key: locus_tag, Value: ['YP_pPCP06']
""",
                -1,
            ),
            (
                """\
type: CDS
location: [4814:5888](-)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478717', 'GeneID:2767721']
    Key: gene, Value: ['pst']
    Key: locus_tag, Value: ['YP_pPCP06']
    Key: note, Value: ['Best Blastp hit =|16082684|ref|NP_395231.1| (NC_003132) pesticin [Yersinia pestis], gi|984824|gb|AAA75369.1| (U31974) pesticin [Yersinia pestis], gi|1488654|emb|CAA63438.1| (X92856) pesticin [Yersinia pestis], gi|2996220|gb|AAC62544.1| (AF053945) pesticin [Yersinia pestis], gi|5763815|emb|CAB53168.1| (AL1099 69) pesticin [Yersinia pestis]']
    Key: product, Value: ['pesticin']
    Key: protein_id, Value: ['NP_995572.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MSDTMVVNGSGGVPAFLFSGSTLSSYRPNFEANSITIALPHYVDLPGRSNFKLMYIMGFPIDTEMEKDSEYSNKIRQESKISKTEGTVSYEQKITVETGQEKDGVKVYRVMVLEGTIAESIEHLDKKENEDILNNNRNRIVLADNTVINFDNISQLKEFLRRSVNIVDHDIFSSNGFEGFNPTSHFPSNPSSDYFNSTGVTFGSGVDLGQRSKQDLLNDGVPQYIADRLDGYYMLRGKEAYDKVRTAPLTLSDNEAHLLSNIYIDKFSHKIEGLFNDANIGLRFSDLPLRTRTALVSIGYQKGFKLSRTAPTVWNKVIAKDWNGLVNAFNNIVDGMSDRRKREGALVQKDIDSGLLK']
""",
                -1,
            ),
            (
                """\
type: variation
location: [5909:5911](+)
qualifiers:
    Key: note, Value: ['compared to AF053945']
    Key: replace, Value: ['']
""",
                1,
            ),
            (
                """\
type: variation
location: [5933:5933](+)
qualifiers:
    Key: note, Value: ['compared to AL109969']
    Key: replace, Value: ['a']
""",
                1,
            ),
            (
                """\
type: variation
location: [5933:5933](+)
qualifiers:
    Key: note, Value: ['compared to AF053945']
    Key: replace, Value: ['aa']
""",
                1,
            ),
            (
                """\
type: variation
location: [5947:5948](+)
qualifiers:
    Key: note, Value: ['compared to AL109969']
    Key: replace, Value: ['c']
""",
                1,
            ),
            (
                """\
type: gene
location: [6004:6421](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767719']
    Key: locus_tag, Value: ['YP_pPCP07']
""",
                1,
            ),
            (
                """\
type: CDS
location: [6004:6421](+)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478718', 'GeneID:2767719']
    Key: locus_tag, Value: ['YP_pPCP07']
    Key: note, Value: ['Best Blastp hit = gi|16082685|ref|NP_395232.1| (NC_003132) hypothetical protein [Yersinia pestis], gi|5763816|emb|CAB53169.1| (AL109969) hypothetical protein [Yersinia pestis]']
    Key: product, Value: ['hypothetical protein']
    Key: protein_id, Value: ['NP_995573.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MKFHFCDLNHSYKNQEGKIRSRKTAPGNIRKKQKGDNVSKTKSGRHRLSKTDKRLLAALVVAGYEERTARDLIQKHVYTLTQADLRHLVSEISNGVGQSQAYDAIYQARRIRLARKYLSGKKPEGVEPREGQEREDLP']
""",
                1,
            ),
            (
                """\
type: variation
location: [6524:6525](+)
qualifiers:
    Key: note, Value: ['compared to AF053945 and AL109969']
    Key: replace, Value: ['c']
""",
                1,
            ),
            (
                """\
type: gene
location: [6663:7602](+)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767715']
    Key: gene, Value: ['pla']
    Key: locus_tag, Value: ['YP_pPCP08']
""",
                1,
            ),
            (
                """\
type: CDS
location: [6663:7602](+)
qualifiers:
    Key: EC_number, Value: ['3.4.23.48']
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478719', 'GeneID:2767715']
    Key: gene, Value: ['pla']
    Key: locus_tag, Value: ['YP_pPCP08']
    Key: note, Value: ['outer membrane protease; involved in virulence in many organisms; OmpT; IcsP; SopA; Pla; PgtE; omptin; in Escherichia coli OmpT can degrade antimicrobial peptides; in Yersinia Pla activates plasminogen during infection; in Shigella flexneria SopA cleaves the autotransporter IcsA']
    Key: product, Value: ['outer membrane protease']
    Key: protein_id, Value: ['NP_995574.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MKKSSIVATIITILSGSANAASSQLIPNISPDSFTVAASTGMLSGKSHEMLYDAETGRKISQLDWKIKNVAILKGDISWDPYSFLTLNARGWTSLASGSGNMDDYDWMNENQSEWTDHSSHPATNVNHANEYDLNVKGWLLQDENYKAGITAGYQETRFSWTATGGSYSYNNGAYTGNFPKGVRVIGYNQRFSMPYIGLAGQYRINDFELNALFKFSDWVRAHDNDEHYMRDLTFREKTSGSRYYGTVINAGYYVTPNAKVFAEFTYSKYDEGKGGTQIIDKNSGDSVSIGGDAAGISNKNYTVTAGLQYRF']
""",
                1,
            ),
            (
                """\
type: misc_feature
location: [6663:7599](+)
qualifiers:
    Key: db_xref, Value: ['CDD:186487']
    Key: gene, Value: ['pla']
    Key: locus_tag, Value: ['YP_pPCP08']
    Key: note, Value: ['Omptin family; Region: Omptin; cl01886']
""",
                1,
            ),
            (
                """\
type: gene
location: [7788:8088](-)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767713']
    Key: locus_tag, Value: ['YP_pPCP09']
""",
                -1,
            ),
            (
                """\
type: CDS
location: [7788:8088](-)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478720', 'GeneID:2767713']
    Key: locus_tag, Value: ['YP_pPCP09']
    Key: note, Value: ['Best Blastp hit = gi|16082687|ref|NP_395234.1| (NC_003132) putative transcriptional regulator [Yersinia pestis], gi|5763818|emb|CAB53171.1| (AL109969) putative transcriptional regulator [Yersinia pestis].']
    Key: product, Value: ['putative transcriptional regulator']
    Key: protein_id, Value: ['NP_995575.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MRTLDEVIASRSPESQTRIKEMADEMILEVGLQMMREELQLSQKQVAEAMGISQPAVTKLEQRGNDLKLATLKRYVEAMGGKLSLDVELPTGRRVAFHV']
""",
                -1,
            ),
            (
                """\
type: misc_feature
location: [7836:7995](-)
qualifiers:
    Key: db_xref, Value: ['CDD:195788']
    Key: locus_tag, Value: ['YP_pPCP09']
    Key: note, Value: ['Helix-turn-helix XRE-family like proteins. Prokaryotic DNA binding proteins belonging to the xenobiotic response element family of transcriptional regulators; Region: HTH_XRE; cl09100']
""",
                -1,
            ),
            (
                """\
type: gene
location: [8087:8360](-)
qualifiers:
    Key: db_xref, Value: ['GeneID:2767714']
    Key: locus_tag, Value: ['YP_pPCP10']
""",
                -1,
            ),
            (
                """\
type: CDS
location: [8087:8360](-)
qualifiers:
    Key: codon_start, Value: ['1']
    Key: db_xref, Value: ['GI:45478721', 'GeneID:2767714']
    Key: locus_tag, Value: ['YP_pPCP10']
    Key: note, Value: ['Best Blastp hit = gi|16082688|ref|NP_395235.1| (NC_003132) hypothetical protein [ Yersinia pestis], gi|5763819|emb|CAB53172.1| (AL109969) hypothetical protein [Yersinia pestis]']
    Key: product, Value: ['hypothetical protein']
    Key: protein_id, Value: ['NP_995576.1']
    Key: transl_table, Value: ['11']
    Key: translation, Value: ['MADLKKLQVYGPELPRPYADTVKGSRYKNMKELRVQFSGRPIRAFYAFDPIRRAIVLCAGDKSNDKRFYEKLVRIAEDEFTAHLNTLESK']
""",
                -1,
            ),
            (
                """\
type: misc_feature
location: [8090:>8357](-)
qualifiers:
    Key: db_xref, Value: ['CDD:194142']
    Key: locus_tag, Value: ['YP_pPCP10']
    Key: note, Value: ['Phage derived protein Gp49-like (DUF891); Region: Gp49; cl01470']
""",
                -1,
            ),
            (
                """\
type: variation
location: [8529:8529](+)
qualifiers:
    Key: note, Value: ['compared to AL109969']
    Key: replace, Value: ['tt']
""",
                1,
            ),
        )
        dbxrefs = ["Project:58037"]
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_15(self):
        path = "GenBank/no_end_marker.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Premature end of file in sequence data
                record = next(records)
        seq = "CTAGCAGCCCGCATCGCCCTCGACGTTGGCGATCATCGTGCGCAGCACCTTGAG...TGA"
        id = "AB070938.1"
        name = "AB070938"
        description = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        annotations = {
            "accessions": ["AB070938"],
            "data_file_division": "BCT",
            "date": "11-OCT-2001",
            "gi": "15823953",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Streptomyces avermitilis",
            "sequence_version": 1,
            "source": "Streptomyces avermitilis",
            "taxonomy": [
                "Bacteria",
                "Actinobacteria",
                "Actinobacteridae",
                "Actinomycetales",
                "Streptomycineae",
                "Streptomycetaceae",
                "Streptomyces",
            ],
            "topology": "linear",
        }
        references = []
        features = (
            (
                """\
type: source
location: [0:6497](+)
qualifiers:
    Key: db_xref, Value: ['taxon:33903']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Streptomyces avermitilis']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_16(self):
        path = "GenBank/wrong_sequence_indent.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Invalid indentation for sequence line
                record = next(records)
        seq = "CTAGCAGCCCGCATCGCCCTCGACGTTGGCGATCATCGTGCGCAGCACCTTGAG...TGA"
        id = "AB070938.1"
        name = "AB070938"
        description = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        annotations = {
            "accessions": ["AB070938"],
            "data_file_division": "BCT",
            "date": "11-OCT-2001",
            "gi": "15823953",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Streptomyces avermitilis",
            "sequence_version": 1,
            "source": "Streptomyces avermitilis",
            "taxonomy": [
                "Bacteria",
                "Actinobacteria",
                "Actinobacteridae",
                "Actinomycetales",
                "Streptomycineae",
                "Streptomycetaceae",
                "Streptomyces",
            ],
            "topology": "linear",
        }
        references = []
        features = (
            (
                """\
type: source
location: [0:6497](+)
qualifiers:
    Key: db_xref, Value: ['taxon:33903']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Streptomyces avermitilis']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_17(self):
        path = "GenBank/invalid_locus_line_spacing.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Attempting to parse malformed locus line
                record = next(records)
        seq = "CTAGCAGCCCGCATCGCCCTCGACGTTGGCGATCATCGTGCGCAGCACCTTGAG...TGA"
        id = "AB070938.1"
        name = "AB070938"
        description = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        annotations = {
            "accessions": ["AB070938"],
            "data_file_division": "BCT",
            "date": "11-OCT-2001",
            "gi": "15823953",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Streptomyces avermitilis",
            "sequence_version": 1,
            "source": "Streptomyces avermitilis",
            "taxonomy": [
                "Bacteria",
                "Actinobacteria",
                "Actinobacteridae",
                "Actinomycetales",
                "Streptomycineae",
                "Streptomycetaceae",
                "Streptomyces",
            ],
            "topology": "linear",
        }
        references = []
        features = (
            (
                """\
type: source
location: [0:6497](+)
qualifiers:
    Key: db_xref, Value: ['taxon:33903']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Streptomyces avermitilis']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_18(self):
        path = "GenBank/empty_feature_qualifier.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            record = next(records)
        seq = "CTAGCAGCCCGCATCGCCCTCGACGTTGGCGATCATCGTGCGCAGCACCTTGAG...TGA"
        id = "AB070938.1"
        name = "AB070938"
        description = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        annotations = {
            "accessions": ["AB070938"],
            "data_file_division": "BCT",
            "date": "11-OCT-2001",
            "gi": "15823953",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Streptomyces avermitilis",
            "sequence_version": 1,
            "source": "Streptomyces avermitilis",
            "taxonomy": [
                "Bacteria",
                "Actinobacteria",
                "Actinobacteridae",
                "Actinomycetales",
                "Streptomycineae",
                "Streptomycetaceae",
                "Streptomyces",
            ],
            "topology": "linear",
        }
        references = []
        features = (
            (
                """\
type: source
location: [0:6497](+)
qualifiers:
    Key: db_xref, Value: ['taxon:33903']
    Key: mol_type, Value: ['genomic DNA']
    Key: note, Value: ["This is a correct note, the following one isn't"]
    Key: organism, Value: ['Streptomyces avermitilis']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_19(self):
        path = "GenBank/invalid_misc_feature.gb"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: line too short to contain a feature
                record = next(records)
        seq = "CTAGCAGCCCGCATCGCCCTCGACGTTGGCGATCATCGTGCGCAGCACCTTGAG...TGA"
        id = "AB070938.1"
        name = "AB070938"
        description = "Streptomyces avermitilis melanin biosynthetic gene cluster"
        annotations = {
            "accessions": ["AB070938"],
            "data_file_division": "BCT",
            "date": "11-OCT-2001",
            "gi": "15823953",
            "keywords": [""],
            "molecule_type": "DNA",
            "organism": "Streptomyces avermitilis",
            "sequence_version": 1,
            "source": "Streptomyces avermitilis",
            "taxonomy": [
                "Bacteria",
                "Actinobacteria",
                "Actinobacteridae",
                "Actinomycetales",
                "Streptomycineae",
                "Streptomycetaceae",
                "Streptomyces",
            ],
            "topology": "linear",
        }
        references = []
        features = (
            (
                """\
type: source
location: [0:6497](+)
qualifiers:
    Key: db_xref, Value: ['taxon:33903']
    Key: mol_type, Value: ['genomic DNA']
    Key: organism, Value: ['Streptomyces avermitilis']
""",
                1,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_feature_parser_20(self):
        path = "GenBank/1MRR_A.gp"
        with open(path) as handle:
            records = GenBank.Iterator(handle, self.feat_parser)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # BiopythonParserWarning: Dropping bond qualifier in feature location
                record = next(records)
        seq = "AYTTFSATKNDQLKEPMFFGQPVQVARYDQQKYDIFEKLIEKQLSFFWRPEEVD...FQL"
        id = "1MRR_A"
        name = "1MRR_A"
        description = (
            "Chain A, Substitution Of Manganese For Iron In Ribonucleotide Reductase"
            " From Escherichia Coli. Spectroscopic And Crystallographic"
            " Characterization"
        )
        annotations = {
            "accessions": ["1MRR_A"],
            "comment": """\
1 Ribonucleotide Reductase R1 Protein.""",
            "data_file_division": "BCT",
            "date": "10-OCT-2012",
            "db_source": (
                "pdb: molecule 1MRR, chain 65, release Aug 28, 2012; deposition: Jul"
                " 28, 1992; class: Reductase(Acting On Ch2); source: Mmdb_id: 50351,"
                " Pdb_id 1: 1MRR; Exp. method: X-Ray Diffraction."
            ),
            "gi": "494379",
            "keywords": [""],
            "molecule_type": "protein",
            "organism": "Escherichia coli",
            "source": "Escherichia coli",
            "taxonomy": [
                "Bacteria",
                "Proteobacteria",
                "Gammaproteobacteria",
                "Enterobacteriales",
                "Enterobacteriaceae",
                "Escherichia",
            ],
            "topology": "linear",
        }
        references = [
            "location: [0:375]\nauthors: Nordlund,P., Sjoberg,B.M. and"
            " Eklund,H.\ntitle: Three-dimensional structure of the free radical protein"
            " of ribonucleotide reductase\njournal: Nature 345 (6276), 593-598"
            " (1990)\nmedline id: \npubmed id: 2190093\ncomment: \n",
            "location: [0:375]\nauthors: Atta,M., Nordlund,P., Aberg,A., Eklund,H. and"
            " Fontecave,M.\ntitle: Substitution of manganese for iron in ribonucleotide"
            " reductase from Escherichia coli. Spectroscopic and crystallographic"
            " characterization\njournal: J. Biol. Chem. 267 (29), 20682-20688"
            " (1992)\nmedline id: \npubmed id: 1328209\ncomment: \n",
            "location: [0:375]\nauthors: Eklund,H. and Nordlund,P.\ntitle: Direct"
            " Submission\njournal: Submitted (28-JUL-1992)\nmedline id: \npubmed id:"
            " \ncomment: \n",
        ]
        features = (
            (
                """\
type: source
location: [0:375]
qualifiers:
    Key: db_xref, Value: ['taxon:562']
    Key: organism, Value: ['Escherichia coli']
""",
                None,
            ),
            (
                """\
type: Region
location: [27:340]
qualifiers:
    Key: db_xref, Value: ['CDD:153108']
    Key: note, Value: ['Ribonucleotide Reductase, R2/beta subunit, ferritin-like diiron-binding domain; cd01049']
    Key: region_name, Value: ['RNRR2']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [34:46]
qualifiers:
    Key: note, Value: ['helix 1']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: Site
location: order{[36:37], [43:44], [108:110], [112:113], [115:117], [119:120], [122:123], [136:138], [140:141]}
qualifiers:
    Key: db_xref, Value: ['CDD:153108']
    Key: note, Value: ['dimer interface [polypeptide binding]']
    Key: site_type, Value: ['other']
""",
                None,
            ),
            (
                """\
type: Site
location: order{[47:48], [83:84], [114:115], [117:118], [121:122], [235:237], [240:241]}
qualifiers:
    Key: db_xref, Value: ['CDD:153108']
    Key: note, Value: ['putative radical transfer pathway']
    Key: site_type, Value: ['other']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [56:65]
qualifiers:
    Key: note, Value: ['helix 2']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [66:87]
qualifiers:
    Key: note, Value: ['helix 3']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: Site
location: order{[83:84], [114:115], [117:118], [203:204], [237:238], [240:241]}
qualifiers:
    Key: db_xref, Value: ['CDD:153108']
    Key: note, Value: ['diiron center [ion binding]']
    Key: site_type, Value: ['other']
""",
                None,
            ),
            (
                """\
type: Het
location: join{[83:84], [114:115], [117:118], [237:238]}
qualifiers:
    Key: heterogen, Value: ['( MN,1000 )']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [101:129]
qualifiers:
    Key: note, Value: ['helix 4']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: Het
location: join{[114:115], [203:204], [237:238], [240:241]}
qualifiers:
    Key: heterogen, Value: ['( MN,1001 )']
""",
                None,
            ),
            (
                """\
type: Site
location: [121:122]
qualifiers:
    Key: db_xref, Value: ['CDD:153108']
    Key: note, Value: ['tyrosyl radical']
    Key: site_type, Value: ['other']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [132:140]
qualifiers:
    Key: note, Value: ['helix 5']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [142:151]
qualifiers:
    Key: note, Value: ['helix 6']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [152:169]
qualifiers:
    Key: note, Value: ['helix 7']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [171:177]
qualifiers:
    Key: note, Value: ['strand 1']
    Key: sec_str_type, Value: ['sheet']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [179:185]
qualifiers:
    Key: note, Value: ['strand 2']
    Key: sec_str_type, Value: ['sheet']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [185:216]
qualifiers:
    Key: note, Value: ['helix 8']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: Het
location: join{[193:194], [271:272]}
qualifiers:
    Key: heterogen, Value: ['( HG,1003 )']
""",
                None,
            ),
            (
                """\
type: Het
location: [195:196]
qualifiers:
    Key: heterogen, Value: ['( HG,1005 )']
""",
                None,
            ),
            (
                """\
type: Het
location: join{[195:196], [195:196]}
qualifiers:
    Key: heterogen, Value: ['( HG,1002 )']
""",
                None,
            ),
            (
                """\
type: Het
location: join{[209:210], [213:214], [213:214]}
qualifiers:
    Key: heterogen, Value: ['( HG,1004 )']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [224:253]
qualifiers:
    Key: note, Value: ['helix 9']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [259:269]
qualifiers:
    Key: note, Value: ['helix 10']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: Bond
location: bond{[267:268], [271:272]}
qualifiers:
    Key: bond_type, Value: ['disulfide']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [269:285]
qualifiers:
    Key: note, Value: ['helix 11']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
            (
                """\
type: Het
location: join{[283:284], [304:305], [308:309], [304:305]}
qualifiers:
    Key: heterogen, Value: ['( HG,1006 )']
""",
                None,
            ),
            (
                """\
type: SecStr
location: [300:319]
qualifiers:
    Key: note, Value: ['helix 12']
    Key: sec_str_type, Value: ['helix']
""",
                None,
            ),
        )
        dbxrefs = []
        self.perform_feature_parser_test(
            record,
            seq,
            id,
            name,
            description,
            annotations,
            references,
            features,
            dbxrefs,
        )

    def test_features_spanning_origin(self):
        """Test that features that span the origin on circular DNA are included correctly for different ways of specifying the topology."""
        # This first one should fail (location of the feature should be set to none), because
        # the file says the sequence is linear, but there is a feature that spans the origin.
        file_fails = "GenBank/addgene-plasmid-11664-sequence-180430.gbk"

        # The right warning is raised
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            record = SeqIO.read(file_fails, "gb")
            self.assertEqual(len(caught), 1)
            self.assertEqual(caught[0].category, BiopythonParserWarning)
            self.assertEqual(
                str(caught[0].message),
                "it appears that '8569..276' is a feature that spans the origin, but the sequence topology is undefined; setting feature location to None.",
            )

        # The last feature location is None
        self.assertIsNone(record.features[-1].location)

        # This one is circular and should include the features that span the origin
        file_succeeds = "GenBank/addgene-plasmid-39296-sequence-49545.gbk"

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            record = SeqIO.read(file_succeeds, "gb")
            # This gives the same error for features that span the origin
            self.assertEqual(len(caught), 2)
            for c in caught:
                self.assertEqual(c.category, BiopythonParserWarning)
                self.assertTrue("unintended behavior" in str(c.message))

        # The last two features should not be none
        self.assertIsNotNone(record.features[-1].location)
        self.assertIsNotNone(record.features[-2].location)


class GenBankTests(unittest.TestCase):
    """GenBank tests."""

    def test_invalid_product_line_raises_value_error(self):
        """Parsing invalid product line."""
        path = "GenBank/invalid_product.gb"
        self.assertRaises(ValueError, SeqIO.read, path, "genbank")

    def test_genbank_read(self):
        """GenBank.read(...) simple test."""
        path = "GenBank/NC_000932.gb"
        with open(path) as handle:
            record = GenBank.read(handle)
        self.assertEqual(["NC_000932"], record.accession)

    def test_genbank_read_multirecord(self):
        """GenBank.read(...) error on multiple record input."""
        path = "GenBank/cor6_6.gb"
        with open(path) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    def test_genbank_read_invalid(self):
        """GenBank.read(...) error on invalid file (e.g. FASTA file)."""
        path = "GenBank/NC_000932.faa"
        with open(path) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    def test_genbank_read_no_origin_no_end(self):
        """GenBank.read(...) error on malformed file."""
        path = "GenBank/no_origin_no_end.gb"
        with open(path) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    # Evil hack with 000 to manipulate sort order to ensure this is tested
    # first (otherwise something silences the warning)
    def test_000_genbank_bad_loc_wrap_warning(self):
        """Feature line wrapping warning."""
        path = "GenBank/bad_loc_wrap.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            with open(path) as handle:
                with self.assertRaises(BiopythonParserWarning) as cm:
                    GenBank.read(handle)
                self.assertEqual(
                    "Non-standard feature line wrapping (didn't break on comma)?",
                    str(cm.exception),
                )

    # Similar hack as we also want to catch that warning here
    def test_001_negative_location_warning(self):
        """Un-parsable feature location warning."""
        path = "GenBank/negative_location.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            with self.assertRaises(BiopythonParserWarning) as cm:
                record = SeqIO.read(path, "genbank")
            self.assertEqual(
                "negative starting position in feature location '-2..492'; setting feature location to None.",
                str(cm.exception),
            )

    def test_001_genbank_bad_origin_wrapping_location(self):
        """Bad origin wrapping."""
        path = "GenBank/bad_origin_wrap_linear.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            with self.assertRaises(BiopythonParserWarning) as cm:
                record = SeqIO.read(path, "genbank")
            self.assertEqual(
                "it appears that '6801..100' is a feature that spans the origin, but the sequence topology is undefined; setting feature location to None.",
                str(cm.exception),
            )

    def test_001_implicit_orign_wrap_fix(self):
        """Attempt to fix implied origin wrapping."""
        path = "GenBank/bad_origin_wrap.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            with self.assertRaises(BiopythonParserWarning) as cm:
                record = SeqIO.read(path, "genbank")
            self.assertEqual(
                str(cm.exception),
                "Attempting to fix invalid location '6801..100' "
                "as it looks like incorrect origin wrapping. "
                "Please fix input file, this could have "
                "unintended behavior.",
            )

    def test_compound_complex_origin_wrap(self):
        """Test the attempts to fix compound complex origin wrapping."""
        from Bio.SeqFeature import CompoundLocation

        path = "GenBank/bad_origin_wrap.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            record = SeqIO.read(path, "genbank")

            self.assertIsInstance(record.features[3].location, CompoundLocation)
            self.assertEqual(
                str(record.features[3].location),
                "join{[<5399:5600](+), [5699:6100](+), [6800:7000](+), [0:100](+)}",
            )

            self.assertIsInstance(record.features[4].location, CompoundLocation)
            self.assertEqual(
                str(record.features[4].location),
                "join{[5399:5600](+), [5699:6100](+), [<6800:7000](+), [0:100](+)}",
            )

            self.assertIsInstance(record.features[5].location, CompoundLocation)
            self.assertEqual(
                str(record.features[5].location),
                "join{[5399:5600](+), [5699:6100](+), [0:100](-), [<6800:7000](-)}",
            )

    def test_implicit_orign_wrap_extract_and_translate(self):
        """Test that features wrapped around origin give expected data."""
        path = "GenBank/bad_origin_wrap_CDS.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with open(path) as handle:
                seq_record = SeqIO.read(handle, "genbank")
        seq_features = seq_record.features
        self.assertEqual(
            seq_features[1].extract(seq_record).seq.lower(),
            "atgccctataaaacccagggctgccttggaaaaggcgcaaccccaaccccctcgagccgcggcatataa",
        )
        self.assertEqual(
            seq_features[2].extract(seq_record).seq.lower(),
            "atgccgcggctcgagggggttggggttgcgccttttccaaggcagccctgggttttatag",
        )
        self.assertEqual(
            seq_features[1].extract(seq_record).seq.translate(),
            "MPYKTQGCLGKGATPTPSSRGI*",
        )
        self.assertEqual(
            seq_features[2].extract(seq_record).seq.translate(), "MPRLEGVGVAPFPRQPWVL*"
        )

    def test_fuzzy_origin_wrap(self):
        """Test features that wrap an origin, and have fuzzy location."""
        path = "GenBank/bad_origin_wrap_fuzzy.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            with self.assertRaises(BiopythonParserWarning) as cm:
                record = SeqIO.read(path, "genbank")
            self.assertEqual(
                str(cm.exception),
                "Attempting to fix invalid location '<2644..159' "
                "as it looks like incorrect origin wrapping. "
                "Please fix input file, this could have "
                "unintended behavior.",
            )

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                with open(path) as handle:
                    seq_record = SeqIO.read(handle, "genbank")
                    self.assertEqual(
                        str(seq_record.features[3].location),
                        "join{[<2643:2686](+), [0:159](+)}",
                    )

    def test_genbank_bad_loc_wrap_parsing(self):
        """Bad location wrapping."""
        path = "GenBank/bad_loc_wrap.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with open(path) as handle:
                record = GenBank.read(handle)
        self.assertEqual(1, len(record.features))
        loc = record.features[0].location
        self.assertEqual(
            loc,
            "join(3462..3615,3698..3978,4077..4307,4408..4797,4876..5028,5141..5332)",
        )

    def test_negative_location(self):
        """Negative feature locations."""
        path = "GenBank/negative_location.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            record = SeqIO.read(path, "genbank")
            self.assertIsNone(record.features[-1].location)

    def test_dot_lineage(self):
        """Missing taxonomy lineage."""
        path = "GenBank/bad_loc_wrap.gb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            record = SeqIO.read(path, "genbank")
        self.assertEqual(record.annotations["organism"], ".")
        self.assertEqual(record.annotations["taxonomy"], [])

    def test_tsa(self):
        """Test TSA annotation parsing."""
        path = "GenBank/tsa_acropora.gb"
        record = SeqIO.read(path, "genbank")
        self.assertIn("tsa", record.annotations)
        self.assertEqual(record.annotations["tsa"], ["GHGH01000001", "GHGH01126539"])

    def test_dblink(self):
        """Parse GenBank record with old DBLINK project entry."""
        path = "GenBank/NC_005816.gb"
        record = SeqIO.read(path, "gb")
        self.assertEqual(record.dbxrefs, ["Project:58037"])
        gb = record.format("gb")
        self.assertIn("\nDBLINK      Project: 58037\n", gb)
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:58037;\nXX\n", embl)

    def test_dblink_two(self):
        """Parse GenBank record with old and new DBLINK project entries."""
        path = "GenBank/NP_416719.gbwithparts"
        record = SeqIO.read(path, "gb")
        self.assertEqual(record.dbxrefs, ["Project:57779", "BioProject:PRJNA57779"])
        gb = record.format("gb")
        self.assertIn(
            """
DBLINK      Project: 57779
            BioProject: PRJNA57779
KEYWORDS    """,
            gb,
        )
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:PRJNA57779;\nXX\n", embl)

    def test_dbline_gb_embl(self):
        """Parse GenBank/EMBL paired records with PR project entry: GenBank."""
        record = SeqIO.read("GenBank/DS830848.gb", "gb")
        self.assertIn("BioProject:PRJNA16232", record.dbxrefs)
        gb = record.format("gb")
        self.assertIn(
            """
DBLINK      BioProject: PRJNA16232
            BioSample: SAMN03004382
KEYWORDS    """,
            gb,
        )
        # Also check EMBL output
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:PRJNA16232;\nXX\n", embl)

    def test_dbline_embl_gb(self):
        """Parse GenBank/EMBL paired records with PR project entry: EMBL."""
        record = SeqIO.read("EMBL/DS830848.embl", "embl")
        # TODO: Should we map this to BioProject:PRJNA16232
        self.assertIn("Project:PRJNA16232", record.dbxrefs)
        gb = record.format("gb")
        self.assertIn(
            """
DBLINK      Project: PRJNA16232
            MD5: 387e72e4f7ae804780d06f875ab3bc41
            ENA: ABJB010000000
            ENA: ABJB000000000
            BioSample: SAMN03004382
KEYWORDS    """,
            gb,
        )
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:PRJNA16232;\nXX\n", embl)

    def test_structured_comment_parsing(self):
        """Structured comment parsing."""
        # GISAID_EpiFlu(TM)Data, HM138502.gbk has both
        # 'comment' and 'structured_comment'
        path = "GenBank/HM138502.gbk"
        record = SeqIO.read(path, "genbank")
        self.assertEqual(
            record.annotations["comment"],
            "Swine influenza A (H1N1) virus isolated during human swine flu\noutbreak"
            " of 2009.",
        )
        self.assertEqual(
            record.annotations["structured_comment"]["GISAID_EpiFlu(TM)Data"][
                "Lineage"
            ],
            "swl",
        )
        self.assertEqual(
            len(record.annotations["structured_comment"]["GISAID_EpiFlu(TM)Data"]), 3
        )
        path = "GenBank/HM138502_output.gbk"
        with open(path) as ifile:
            self.assertEqual(record.format("gb"), ifile.read())
        # FluData structured comment
        path = "GenBank/EU851978.gbk"
        record = SeqIO.read(path, "genbank")
        self.assertEqual(
            record.annotations["structured_comment"]["FluData"]["LabID"], "2008704957"
        )
        self.assertEqual(len(record.annotations["structured_comment"]["FluData"]), 5)
        path = "GenBank/EU851978_output.gbk"
        with open(path) as ifile:
            self.assertEqual(record.format("gb"), ifile.read())
        # Assembly-Data structured comment
        path = "GenBank/KF527485.gbk"
        record = SeqIO.read(path, "genbank")
        self.assertEqual(
            record.annotations["structured_comment"]["Assembly-Data"][
                "Assembly Method"
            ],
            "Lasergene v. 10",
        )
        self.assertEqual(
            len(record.annotations["structured_comment"]["Assembly-Data"]), 2
        )
        path = "GenBank/KF527485_output.gbk"
        with open(path) as ifile:
            self.assertEqual(record.format("gb"), ifile.read())
        # No structured comment in NC_000932.gb, just a regular comment
        path = "GenBank/NC_000932.gb"
        record = SeqIO.read(path, "genbank")
        self.assertNotIn("structured_comment", record.annotations)
        self.assertEqual(
            record.annotations["comment"],
            "REVIEWED REFSEQ: This record has been curated by NCBI staff. The\n"
            "reference sequence was derived from AP000423.\n"
            "COMPLETENESS: full length.",
        )

    def test_multiline_structured_comment_parsing(self):
        """Multiline structured comment parsing."""
        # GU949562.1, MIENS-Data, environment has value on multiple lines
        path = "GenBank/GU949562.1.gb"
        record = SeqIO.read(path, "genbank")
        self.assertEqual(
            record.annotations["structured_comment"]["MIENS-Data"]["environment"],
            "Temperate shelf and sea biome [ENVO:00000895], "
            "coastal water body [ENVO:02000049], "
            "coastal water [ENVO:00002150]",
        )

    def test_malformed_structured_comment_parsing(self):
        """Test malformed structured comment gives warning.

        The comment will be ignored if it is not read by the parser AYW00820.1;
        Malformed key-value delimiter used. Should be " :: ", but the record uses ": "
        """
        path = "GenBank/invalid_structured_comment.gb"

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            record = SeqIO.read(path, "genbank")
            self.assertNotIn("structured_comment", record.annotations)
            self.assertIn(
                "Structured comment not parsed for AYW00820.", str(caught[0].message)
            )

    def test_locus_line_topogoly(self):
        """Test if chromosome topology is conserved."""
        record = SeqIO.read("GenBank/DS830848.gb", "genbank")
        self.assertEqual(record.annotations["topology"], "linear")
        out_handle = StringIO()
        SeqIO.write([record], out_handle, "genbank")
        first_line = out_handle.getvalue().split("\n")[0]
        self.assertIn("linear", first_line)
        with open("GenBank/DS830848.gb") as fh:
            orig_first_line = fh.readline().strip()
        self.assertEqual(first_line, orig_first_line)

    def test_qualifier_order(self):
        """Check the qualifier order is preserved."""
        record = SeqIO.read("GenBank/DS830848.gb", "gb")
        f = record.features[0]
        self.assertEqual(
            list(f.qualifiers),
            ["organism", "mol_type", "strain", "db_xref", "dev_stage"],
        )

    def test_qualifier_escaping_read(self):
        """Check qualifier escaping is preserved when parsing."""
        # Make sure parsing improperly escaped qualifiers raises a warning
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            record = SeqIO.read("GenBank/qualifier_escaping_read.gb", "gb")
            self.assertEqual(len(caught), 4)
            self.assertEqual(caught[0].category, BiopythonParserWarning)
            self.assertEqual(
                str(caught[0].message),
                'The NCBI states double-quote characters like " should be escaped'
                ' as "" (two double - quotes), but here it was not: '
                "%r" % 'One missing ""quotation mark" here',
            )
        # Check records parsed as expected
        f1 = record.features[0]
        f2 = record.features[1]
        f3 = record.features[2]
        f4 = record.features[3]
        f5 = record.features[4]
        self.assertEqual(f1.qualifiers["note"][0], '"This" is "already" "escaped"')
        self.assertEqual(f2.qualifiers["note"][0], 'One missing "quotation mark" here')
        self.assertEqual(f3.qualifiers["note"][0], 'End not properly "escaped"')
        self.assertEqual(f4.qualifiers["note"][0], '"Start" not properly escaped')
        self.assertEqual(f5.qualifiers["note"][0], 'Middle not "properly" escaped')

    def test_qualifier_escaping_write(self):
        """Check qualifier escaping is preserved when writing."""
        # Write some properly escaped qualifiers and test
        genbank_out = "GenBank/qualifier_escaping_write.gb"
        record = SeqIO.read(genbank_out, "gb")
        f1 = record.features[0]
        f2 = record.features[1]
        f1.qualifiers["note"][0] = '"Should" now "be" escaped in "file"'
        f2.qualifiers["note"][0] = '"Should also be escaped in file"'
        SeqIO.write(record, genbank_out, "gb")
        # Read newly escaped qualifiers and test
        record = SeqIO.read(genbank_out, "gb")
        f1 = record.features[0]
        f2 = record.features[1]
        self.assertEqual(
            f1.qualifiers["note"][0], '"Should" now "be" escaped in "file"'
        )
        self.assertEqual(f2.qualifiers["note"][0], '"Should also be escaped in file"')

    def test_long_names(self):
        """Various GenBank names which push the column based LOCUS line."""
        original = SeqIO.read("GenBank/iro.gb", "gb")
        self.assertEqual(len(original), 1326)
        # Acceptability of LOCUS line with length > 80
        # invalidates some of these tests
        for name, seq_len, ok in [
            ("short", 1, True),
            ("max_length_of_16", 1000, True),
            ("overly_long_at_17", 1000, True),
            ("excessively_long_at_22", 99999, True),
            ("excessively_long_at_22", 100000, True),
            ("pushing_the_limits_at_24", 999, True),
            ("pushing_the_limits_at_24", 1000, True),
            ("old_max_name_length_was_26", 10, True),  # 2 digits
            ("old_max_name_length_was_26", 9, True),
        ]:  # 1 digit
            # Make the length match the desired target
            record = original[:]
            # TODO - Implement Seq * int
            record.seq = Seq("N" * seq_len)
            record.annotations["molecule_type"] = original.annotations["molecule_type"]
            # Set the identifier to the desired name
            record.id = record.name = name
            # Attempt to output the record...
            if not ok:
                # e.g. ValueError:
                # Locus identifier 'excessively_long_at_22' is too long
                self.assertRaises(ValueError, record.format, "gb")
                continue
            with warnings.catch_warnings():
                # e.g. BiopythonWarning: Stealing space from length
                # field to allow long name in LOCUS line
                warnings.simplefilter("ignore", BiopythonWarning)
                # output = record.format("gb")
                handle = StringIO()
                self.assertEqual(1, SeqIO.write(record, handle, "gb"))
            handle.seek(0)
            line = handle.readline()
            self.assertIn(f" {name} ", line)
            self.assertIn(" %i bp " % seq_len, line)
            # Splitting based on whitespace rather than position due to
            # updated GenBank specification
            name_and_length = line.split()[1:3]
            self.assertEqual(name_and_length, [name, str(seq_len)], line)
            handle.seek(0)
            with warnings.catch_warnings():
                # e.g. BiopythonParserWarning: GenBank LOCUS line
                # identifier over 16 characters
                warnings.simplefilter("ignore", BiopythonWarning)
                new = SeqIO.read(handle, "gb")
            self.assertEqual(name, new.name)
            self.assertEqual(seq_len, len(new))

    def test_genbank_date_default(self):
        """Check if default date is handled correctly."""
        sequence_object = Seq("ATGC")
        # check if default value is inserted correctly
        record = SeqRecord(
            sequence_object,
            id="123456789",
            name="UnitTest",
            description="Test case for date parsing",
            annotations={"molecule_type": "DNA"},
        )
        handle = StringIO()
        SeqIO.write(record, handle, "genbank")
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "01-JAN-1980")

    def test_genbank_date_correct(self):
        """Check if user provided date is inserted correctly."""
        sequence_object = Seq("ATGC")
        record = SeqRecord(
            sequence_object,
            id="123456789",
            name="UnitTest",
            description="Test case for date parsing",
            annotations={"molecule_type": "DNA"},
        )
        record.annotations["date"] = "24-DEC-2015"
        handle = StringIO()
        SeqIO.write(record, handle, "genbank")
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "24-DEC-2015")

    def test_genbank_date_list(self):
        """Check if date lists are handled correctly."""
        sequence_object = Seq("ATGC")
        record = SeqRecord(
            sequence_object,
            id="123456789",
            name="UnitTest",
            description="Test case for date parsing",
            annotations={"molecule_type": "DNA"},
        )
        record.annotations["date"] = ["24-DEC-2015"]
        handle = StringIO()
        SeqIO.write(record, handle, "genbank")
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "24-DEC-2015")

        record = SeqRecord(
            sequence_object,
            id="123456789",
            name="UnitTest",
            description="Test case for date parsing",
            annotations={"molecule_type": "DNA"},
        )
        record.annotations["date"] = ["24-DEC-2015", "25-JAN-2016"]
        handle = StringIO()
        SeqIO.write(record, handle, "genbank")
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "01-JAN-1980")

    def test_genbank_date_datetime(self):
        """Check if datetime objects are handled correctly."""
        sequence_object = Seq("ATGC")
        record = SeqRecord(
            sequence_object,
            id="123456789",
            name="UnitTest",
            description="Test case for date parsing",
            annotations={"molecule_type": "DNA"},
        )
        record.annotations["date"] = datetime(2000, 2, 2)
        handle = StringIO()
        SeqIO.write(record, handle, "genbank")
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "02-FEB-2000")

    def test_genbank_date_invalid(self):
        """Check if invalid dates are treated as default."""
        invalid_dates = ("invalid date", "29-2-1981", "35-1-2018", "1-1-80", "1-9-99")

        sequence_object = Seq("ATGC")
        for invalid_date in invalid_dates:
            record = SeqRecord(
                sequence_object,
                id="123456789",
                name="UnitTest",
                description="Test case for date parsing",
                annotations={"molecule_type": "DNA"},
            )

            record.annotations["date"] = invalid_date
            handle = StringIO()
            SeqIO.write(record, handle, "genbank")
            handle.seek(0)
            gb = SeqIO.read(handle, "gb")
            self.assertEqual(gb.annotations["date"], "01-JAN-1980")

    def test_longer_locus_line(self):
        """Check that we can read and write files with longer locus lines."""
        # Create example file from existing file
        path = "GenBank/DS830848.gb"
        with open(path) as inhandle:
            data = inhandle.readlines()
        data[0] = (
            "LOCUS       AZZZAA021234567891234 2147483647 bp    DNA     linear   PRI"
            " 15-OCT-2018\n"
        )

        # Create memory file from modified genbank file
        in_tmp = StringIO()
        in_tmp.writelines(data)
        in_tmp.seek(0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            in_tmp.seek(0)
            record = SeqIO.read(in_tmp, "genbank")

            # Create temporary output memory file
            out_tmp = StringIO()
            SeqIO.write(record, out_tmp, "genbank")

            # Check that the written file can be read back in
            out_tmp.seek(0)
            record_in = SeqIO.read(out_tmp, "genbank")
            self.assertEqual(record_in.id, "DS830848.1")
            self.assertEqual(record_in.name, "AZZZAA021234567891234")
            self.assertEqual(len(record_in.seq), 2147483647)

    if sys.maxsize > 2147483647:

        def test_extremely_long_sequence(self):
            """Tests if extremely long sequences can be read.

            This is only run if sys.maxsize > 2147483647.
            """
            # Create example file from existing file
            path = "GenBank/DS830848.gb"
            with open(path) as inhandle:
                data = inhandle.readlines()
            data[0] = (
                "LOCUS       AZZZAA02123456789 10000000000 bp    DNA     linear   PRI"
                " 15-OCT-2018\n"
            )

            # Create memory file from modified genbank file
            in_tmp = StringIO()
            in_tmp.writelines(data)
            in_tmp.seek(0)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                in_tmp.seek(0)
                record = SeqIO.read(in_tmp, "genbank")

                # Create temporary output memory file
                out_tmp = StringIO()
                SeqIO.write(record, out_tmp, "genbank")

                # Check that the written file can be read back in
                out_tmp.seek(0)
                record_in = SeqIO.read(out_tmp, "genbank")
                self.assertEqual(record_in.id, "DS830848.1")
                self.assertEqual(record_in.name, "AZZZAA02123456789")
                self.assertEqual(len(record_in.seq), 10000000000)

            def read_longer_than_maxsize():
                path = "GenBank/DS830848.gb"
                with open(path) as inhandle:
                    data2 = inhandle.readlines()
                data2[0] = (
                    "LOCUS       AZZZAA02123456789 "
                    + str(sys.maxsize + 1)
                    + " bp    DNA     linear   PRI 15-OCT-2018\n"
                )

                long_in_tmp = StringIO()
                long_in_tmp.writelines(data2)
                long_in_tmp.seek(0)
                record = SeqIO.read(long_in_tmp, "genbank")

            self.assertRaises(ValueError, read_longer_than_maxsize)


class LineOneTests(unittest.TestCase):
    """Check GenBank/EMBL topology / molecule_type parsing."""

    def test_topology_genbank(self):
        """Check GenBank LOCUS line parsing."""
        # This is a bit low level,
        # but can test parsing the LOCUS line only
        tests = [
            ("LOCUS       U00096", None, None, None, None),
            # This example is actually fungal,
            # accession U49845 from Saccharomyces cerevisiae:
            (
                "LOCUS       SCU49845     5028 bp    DNA             PLN      "
                " 21-JUN-1999",
                None,
                "DNA",
                "PLN",
                None,
            ),
            (
                "LOCUS       AB070938                6497 bp    DNA     linear   BCT"
                " 11-OCT-2001",
                "linear",
                "DNA",
                "BCT",
                None,
            ),
            (
                "LOCUS       NC_005816               9609 bp    DNA     circular BCT"
                " 21-JUL-2008",
                "circular",
                "DNA",
                "BCT",
                None,
            ),
            (
                "LOCUS       SCX3_BUTOC                64 aa            linear   INV"
                " 16-OCT-2001",
                "linear",
                None,
                "INV",
                None,
            ),
            (
                "LOCUS       pEH010                  5743 bp    DNA     circular",
                "circular",
                "DNA",
                None,
                [BiopythonParserWarning],
            ),
            # This is a test of the format > 80 chars long
            (
                "LOCUS       AZZZAA02123456789 1000000000 bp    DNA     linear   PRI"
                " 15-OCT-2018",
                "linear",
                "DNA",
                "PRI",
                None,
            ),
        ]
        for (line, topo, mol_type, div, warning_list) in tests:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                scanner = GenBank.Scanner.GenBankScanner()
                consumer = GenBank._FeatureConsumer(1, GenBank.FeatureValueCleaner)
                scanner._feed_first_line(consumer, line)
                t = consumer.data.annotations.get("topology", None)
                self.assertEqual(
                    t, topo, f"Wrong topology {t!r} not {topo!r} from {line!r}"
                )
                mt = consumer.data.annotations.get("molecule_type", None)
                self.assertEqual(
                    mt,
                    mol_type,
                    f"Wrong molecule_type {mt!r} not {mol_type!r} from {line!r}",
                )
                d = consumer.data.annotations.get("data_file_division", None)
                self.assertEqual(
                    d, div, f"Wrong division {d!r} not {div!r} from {line!r}"
                )
                if warning_list is None:
                    self.assertEqual(len(caught), 0)
                else:
                    self.assertEqual(len(caught), len(warning_list))
                    for i, warning_class in enumerate(warning_list):
                        self.assertEqual(caught[i].category, warning_class)

    def test_topology_embl(self):
        """Check EMBL ID line parsing."""
        # This is a bit low level, but can test parsing the ID line only
        tests = [
            # Modern examples with sequence version
            (
                "ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.",
                "linear",
                "mRNA",
                "PLN",
            ),
            (
                "ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.",
                "linear",
                "genomic DNA",
                "MAM",
            ),
            # Example to match GenBank example used above:
            (
                "ID   U49845; SV 1; linear; genomic DNA; STD; FUN; 5028 BP.",
                "linear",
                "genomic DNA",
                "FUN",
            ),
            # Old examples:
            (
                "ID   BSUB9999   standard; circular DNA; PRO; 4214630 BP.",
                "circular",
                "DNA",
                "PRO",
            ),
            ("ID   SC10H5 standard; DNA; PRO; 4870 BP.", None, "DNA", "PRO"),
            # Patent example from 2016-06-10
            # ftp://ftp.ebi.ac.uk/pub/databases/embl/patent/
            (
                "ID   A01679; SV 1; linear; unassigned DNA; PAT; MUS; 12 BP.",
                "linear",
                "unassigned DNA",
                "MUS",
            ),
            # Old patent examples
            ("ID   NRP_AX000635; PRT; NR1; 15 SQ", None, None, "NR1"),
            ("ID   NRP0000016E; PRT; NR2; 5 SQ", None, None, "NR2"),
            # KIPO patent examples
            ("ID   DI500001       STANDARD;      PRT;   111 AA.", None, None, None),
            ("ID   DI644510   standard; PRT;  1852 AA.", None, None, None),
        ]
        for (line, topo, mol_type, div) in tests:
            scanner = GenBank.Scanner.EmblScanner()
            consumer = GenBank._FeatureConsumer(1, GenBank.FeatureValueCleaner)
            scanner._feed_first_line(consumer, line)
            t = consumer.data.annotations.get("topology", None)
            self.assertEqual(
                t, topo, f"Wrong topology {t!r} not {topo!r} from {line!r}"
            )
            mt = consumer.data.annotations.get("molecule_type", None)
            self.assertEqual(
                mt,
                mol_type,
                f"Wrong molecule_type {mt!r} not {mol_type!r} from {line!r}",
            )
            d = consumer.data.annotations.get("data_file_division", None)
            self.assertEqual(d, div, f"Wrong division {d!r} not {div!r} from {line!r}")

    def test_first_line_imgt(self):
        """Check IMGT ID line parsing."""
        # This is a bit low level, but can test parsing the ID line only
        tests = [
            ("ID   HLA00001   standard; DNA; HUM; 3503 BP.", None, "DNA", "HUM"),
            ("ID   HLA00001; SV 1; standard; DNA; HUM; 3503 BP.", None, "DNA", "HUM"),
        ]
        for (line, topo, mol_type, div) in tests:
            scanner = GenBank.Scanner._ImgtScanner()
            consumer = GenBank._FeatureConsumer(1, GenBank.FeatureValueCleaner)
            scanner._feed_first_line(consumer, line)
            t = consumer.data.annotations.get("topology", None)
            self.assertEqual(
                t, topo, f"Wrong topology {t!r} not {topo!r} from {line!r}"
            )
            mt = consumer.data.annotations.get("molecule_type", None)
            self.assertEqual(
                mt,
                mol_type,
                f"Wrong molecule_type {mt!r} not {mol_type!r} from {line!r}",
            )
            d = consumer.data.annotations.get("data_file_division", None)
            self.assertEqual(d, div, f"Wrong division {d!r} not {div!r} from {line!r}")


class OutputTests(unittest.TestCase):
    """GenBank output tests."""

    def test_mad_dots(self):
        """Writing and reading back accesssion.version variants."""
        for identifier in ["example", "example.1a", "example.1.2", "example.1-2"]:
            old = SeqRecord(
                Seq("ACGT"),
                id=identifier,
                name=identifier,
                description="mad dots",
                annotations={"molecule_type": "DNA"},
            )
            new = SeqIO.read(StringIO(old.format("gb")), "gb")
            self.assertEqual(old.id, new.id)
            self.assertEqual(old.name, new.name)
            self.assertEqual(old.description, new.description)
            self.assertEqual(old.seq, new.seq)

    def test_seqrecord_default_description(self):
        """Read in file using SeqRecord default description."""
        old = SeqRecord(
            Seq("ACGT"),
            id="example",
            name="short",
            annotations={"molecule_type": "DNA"},
        )
        self.assertEqual(old.description, "<unknown description>")
        txt = old.format("gb")
        self.assertIn("DEFINITION  .\n", txt)
        new = SeqIO.read(StringIO(txt), "gb")
        self.assertEqual(old.id, new.id)
        self.assertEqual(old.name, new.name)
        self.assertEqual("", new.description)
        self.assertEqual(old.seq, new.seq)

    # Evil hack with 000 to manipulate sort order to ensure this is
    # tested first (otherwise something silences the warning)
    def test_000_write_invalid_but_parsed_locus_line(self):
        """Make sure we survive writing slightly invalid LOCUS lines we could parse."""
        # grab a valid file
        path = "GenBank/NC_005816.gb"
        with open(path) as handle:
            lines = handle.readlines()

        # futz with the molecule type to make it lower case
        invalid_line = (
            "LOCUS       NC_005816               9609 bp    dna     circular BCT"
            " 21-JUL-2008\n"
        )
        lines[0] = invalid_line
        fake_handle = StringIO("".join(lines))

        # Make sure parsing this actually raises a warning
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            rec = SeqIO.read(fake_handle, "genbank")
            self.assertEqual(len(caught), 1)
            self.assertEqual(caught[0].category, BiopythonParserWarning)
            self.assertEqual(
                str(caught[0].message),
                "Non-upper case molecule type in LOCUS line: dna",
            )

        out_handle = StringIO()

        ret = SeqIO.write([rec], out_handle, "genbank")
        self.assertEqual(ret, 1)

        out_handle.seek(0)
        out_lines = out_handle.readlines()
        self.assertEqual(out_lines[0], invalid_line)

    def test_write_tsa_data_division(self):
        """Make sure we don't kill the TSA data_file_division for TSA files."""
        with open("GenBank/tsa_acropora.gb") as infile:
            rec = SeqIO.read(infile, "genbank")
            infile.seek(0)
            first_line = infile.readline()

        outfile = StringIO()
        SeqIO.write([rec], outfile, "genbank")
        outfile.seek(0)
        first_line_written = outfile.readline()

        # ideally, we'd be able to compare these directly, but we also
        # break the "units" field at the moment, so use split instead
        original_division = first_line.split()[-2]
        written_division = first_line_written.split()[-2]

        self.assertEqual(original_division, written_division)


class GenBankScannerTests(unittest.TestCase):
    """GenBank Scanner tests, test parsing gbk and embl files."""

    gb_s = GenBank.Scanner.GenBankScanner()

    def gb_to_l_cds_f(self, filename, tags2id=None):
        """Gb file to Seq list parse CDS features."""
        with open(filename) as handle:
            if tags2id:
                l_cds_f = list(self.gb_s.parse_cds_features(handle, tags2id=tags2id))
            else:
                l_cds_f = list(self.gb_s.parse_cds_features(handle))
        return l_cds_f

    def gb_to_l_r(self, filename, do_features=False):
        """Gb file to Seq list parse records."""
        with open(filename) as handle:
            l_gb_r = list(self.gb_s.parse_records(handle, do_features=do_features))
        return l_gb_r

    def test_genbank_cds_interaction(self):
        """Test CDS interaction, parse CDS features on gb(k) files."""
        # Test parse CDS features on NC_000932.gb
        l_cds_f = self.gb_to_l_cds_f("GenBank/NC_000932.gb")
        # number of records, should be 85
        self.assertEqual(len(l_cds_f), 85)
        # Seq ID
        self.assertEqual(l_cds_f[0].id, "NP_051037.1")
        self.assertEqual(l_cds_f[84].id, "NP_051123.1")

        # Test parse CDS features on NC_005816.gb, Tag to ID
        l_cds_f = self.gb_to_l_cds_f(
            "GenBank/NC_005816.gb", tags2id=("gene", "locus_tag", "product")
        )
        # number of records, should be 10
        self.assertEqual(len(l_cds_f), 10)
        # Seq ID
        self.assertEqual(l_cds_f[0].id, "<unknown id>")
        self.assertEqual(l_cds_f[0].name, "YP_pPCP01")

        # Test parse CDS features on
        # NC_000932.gb and NC_005816.gb combined
        l_cds_f1 = self.gb_to_l_cds_f(
            "GenBank/NC_000932.gb", tags2id=("gene", "locus_tag", "product")
        )
        l_cds_f2 = self.gb_to_l_cds_f(
            "GenBank/NC_005816.gb", tags2id=("gene", "locus_tag", "product")
        )
        l_cds_combined = l_cds_f1 + l_cds_f2
        # number of records combined, should be 95
        self.assertEqual(len(l_cds_combined), 95)
        # Seq ID
        self.assertEqual(l_cds_combined[0].id, "rps12")
        self.assertEqual(l_cds_combined[0].description, "ribosomal protein S12")
        self.assertEqual(l_cds_combined[94].id, "<unknown id>")
        self.assertEqual(l_cds_combined[94].description, "hypothetical protein")

    def test_genbank_interaction(self):
        """Test GenBank records interaction on gbk files."""
        # Test parse records, on NC_005816, do_features False
        l_r = self.gb_to_l_r("GenBank/NC_005816.gb", do_features=False)
        # number of records, should be 1
        self.assertEqual(len(l_r), 1)
        self.assertEqual(l_r[0].id, "NC_005816.1")
        self.assertEqual(l_r[0].name, "NC_005816")
        self.assertEqual(
            l_r[0].description,
            "Yersinia pestis biovar "
            "Microtus str. 91001 plasmid "
            "pPCP1, complete sequence",
        )
        self.assertEqual(len(l_r[0].features), 0)

        # Test parse records on NC_005816, do_features True
        l_r = self.gb_to_l_r("GenBank/NC_005816.gb", do_features=True)
        # number of records, should be 1
        self.assertEqual(len(l_r), 1)
        self.assertEqual(l_r[0].id, "NC_005816.1")
        self.assertEqual(l_r[0].name, "NC_005816")
        self.assertEqual(
            l_r[0].description,
            "Yersinia pestis biovar "
            "Microtus str. 91001 plasmid "
            "pPCP1, complete sequence",
        )
        self.assertEqual(len(l_r[0].features), 41)

        # Test parse records on "GenBank/NC_000932.gb",
        # do_features False
        l_r = self.gb_to_l_r("GenBank/NC_000932.gb", do_features=False)
        # number of records, should be 1
        self.assertEqual(len(l_r), 1)
        self.assertEqual(l_r[0].id, "NC_000932.1")
        self.assertEqual(l_r[0].name, "NC_000932")
        self.assertEqual(
            l_r[0].description, "Arabidopsis thaliana chloroplast, complete genome"
        )
        self.assertEqual(len(l_r[0].features), 0)

        # Test parse records on NC_000932, do_features True
        l_r = self.gb_to_l_r("GenBank/NC_000932.gb", do_features=True)
        # number of records, should be 1
        self.assertEqual(len(l_r), 1)
        self.assertEqual(l_r[0].id, "NC_000932.1")
        self.assertEqual(l_r[0].name, "NC_000932")
        self.assertEqual(
            l_r[0].description, "Arabidopsis thaliana chloroplast, complete genome"
        )
        self.assertEqual(len(l_r[0].features), 259)

    def test_embl_cds_interaction(self):
        """Test EMBL CDS interaction, parse CDS features on embl files."""
        embl_s = GenBank.Scanner.EmblScanner()

        # Test parse CDS features on embl_file
        with open("EMBL/AE017046.embl") as handle_embl7046:
            l_cds_f = list(embl_s.parse_cds_features(handle_embl7046))
        # number of records, should be 10
        self.assertEqual(len(l_cds_f), 10)
        # Seq ID
        self.assertEqual(l_cds_f[0].id, "AAS58758.1")
        self.assertEqual(l_cds_f[0].description, "putative transposase")

    def test_embl_record_interaction(self):
        """Test EMBL Record interaction on embl files."""
        embl_s = GenBank.Scanner.EmblScanner()

        #  Test parse records on embl_file
        with open("EMBL/AE017046.embl") as handle_embl7046:
            l_embl_r = list(embl_s.parse_records(handle_embl7046, do_features=True))
        # number of records, should be 1
        self.assertEqual(len(l_embl_r), 1)
        self.assertEqual(l_embl_r[0].id, "AE017046.1")
        self.assertEqual(
            l_embl_r[0].description,
            "Yersinia pestis biovar Microtus "
            "str. 91001 plasmid pPCP1, complete "
            "sequence.",
        )
        self.assertEqual(len(l_embl_r[0].features), 29)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
