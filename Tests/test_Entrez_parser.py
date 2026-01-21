# Copyright 2008-2010 by Michiel de Hoon.  All rights reserved.
# Revisions copyright 2009-2016 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Testing code for Bio.Entrez parsers."""

import os
import pickle
import unittest
from io import BytesIO

from Bio import Entrez
from Bio import StreamModeError


class GeneralTests(unittest.TestCase):
    """General tests for Bio.Entrez."""

    def test_closed_file(self):
        """Test parsing closed file fails gracefully."""
        stream = open("Entrez/taxonomy.xml", "rb")
        stream.close()
        self.assertRaises(ValueError, Entrez.read, stream)

    def test_read_bytes_stream(self):
        """Test reading a file opened in binary mode."""
        with open("Entrez/taxonomy.xml", "rb") as stream:
            records = Entrez.read(stream)
        self.assertEqual(len(records), 2)
        for record in records:
            self.assertIn("TaxId", record)
            self.assertIn("ScientificName", record)

    def test_parse_bytes_stream(self):
        """Test parsing a file opened in binary mode."""
        with open("Entrez/taxonomy.xml", "rb") as stream:
            records = Entrez.parse(stream)
            n = 0
            for record in records:
                self.assertIn("TaxId", record)
                self.assertIn("ScientificName", record)
                n += 1
        self.assertEqual(n, 2)

    def test_read_text_file(self):
        """Test reading a file opened in text mode."""
        message = "^the XML file must be opened in binary mode.$"
        with open("Entrez/taxonomy.xml") as stream:
            with self.assertRaisesRegex(StreamModeError, message):
                Entrez.read(stream)

    def test_parse_text_file(self):
        """Test parsing a file opened in text mode."""
        message = "^the XML file must be opened in binary mode.$"
        with open("Entrez/taxonomy.xml") as stream:
            records = Entrez.parse(stream)
            with self.assertRaisesRegex(StreamModeError, message):
                next(records)

    def test_BytesIO(self):
        """Test parsing a BytesIO stream (bytes not string)."""
        with open("Entrez/taxonomy.xml", "rb") as stream:
            data = stream.read()
        stream = BytesIO(data)
        records = Entrez.read(stream)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["ScientificName"], "Canis lupus familiaris")
        self.assertEqual(records[1]["ScientificName"], "Felis catus")
        stream.close()

    def test_pickle(self):
        """Test if records created by the parser can be pickled."""
        directory = "Entrez"
        filenames = os.listdir(directory)
        for filename in sorted(filenames):
            basename, extension = os.path.splitext(filename)
            if extension != ".xml":
                continue
            if filename in (
                "biosample.xml",  # DTD not specified in XML file
                "einfo4.xml",  # XML corrupted
                "journals.xml",  # Missing XML declaration
            ):
                continue
            path = os.path.join(directory, filename)
            with open(path, "rb") as stream:
                if filename in ("epost2.xml", "epost3.xml", "esummary8.xml"):
                    # these include an ErrorElement
                    record = Entrez.read(stream, ignore_errors=True)
                else:
                    record = Entrez.read(stream)
            with BytesIO() as stream:
                pickle.dump(record, stream)
                stream.seek(0)
                pickled_record = pickle.load(stream)
            self.assertEqual(record, pickled_record)


class EInfoTest(unittest.TestCase):
    """Tests for parsing XML output returned by EInfo."""

    def test_list(self):
        """Test parsing database list returned by EInfo."""
        # To create the XML file, use
        # >>> Bio.Entrez.einfo()
        with open("Entrez/einfo1.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(
            record["DbList"],
            [
                "pubmed",
                "protein",
                "nuccore",
                "ipg",
                "nucleotide",
                "structure",
                "genome",
                "annotinfo",
                "assembly",
                "bioproject",
                "biosample",
                "blastdbinfo",
                "books",
                "cdd",
                "clinvar",
                "gap",
                "gapplus",
                "grasp",
                "dbvar",
                "gene",
                "gds",
                "geoprofiles",
                "medgen",
                "mesh",
                "nlmcatalog",
                "omim",
                "orgtrack",
                "pmc",
                "proteinclusters",
                "pcassay",
                "protfam",
                "pccompound",
                "pcsubstance",
                "seqannot",
                "snp",
                "sra",
                "taxonomy",
                "biocollections",
                "gtr",
            ],
        )

    def test_pubmed1(self):
        """Test parsing database info returned by EInfo."""
        # To create the XML file, use
        # >>> Bio.Entrez.einfo(db="pubmed")
        with open("Entrez/einfo2.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record["DbInfo"]), 8)
        self.assertEqual(record["DbInfo"]["DbName"], "pubmed")
        self.assertEqual(record["DbInfo"]["MenuName"], "PubMed")
        self.assertEqual(record["DbInfo"]["Description"], "PubMed bibliographic record")
        self.assertEqual(record["DbInfo"]["Count"], "39730388")
        self.assertEqual(record["DbInfo"]["LastUpdate"], "2025/11/27 06:33")

        self.assertEqual(len(record["DbInfo"]["FieldList"]), 50)

        self.assertEqual(record["DbInfo"]["FieldList"][0]["Name"], "ALL")
        self.assertEqual(record["DbInfo"]["FieldList"][0]["FullName"], "All Fields")
        self.assertEqual(
            record["DbInfo"]["FieldList"][0]["Description"],
            "All terms from all searchable fields",
        )

        self.assertEqual(record["DbInfo"]["FieldList"][0]["IsNumerical"], "N")
        self.assertEqual(record["DbInfo"]["FieldList"][0]["SingleToken"], "N")
        self.assertEqual(record["DbInfo"]["FieldList"][0]["Hierarchy"], "N")
        self.assertEqual(record["DbInfo"]["FieldList"][0]["IsHidden"], "N")

        self.assertEqual(len(record["DbInfo"]["LinkList"]), 57)

        self.assertEqual(record["DbInfo"]["LinkList"][0]["Name"], "pubmed_assembly")
        self.assertEqual(record["DbInfo"]["LinkList"][0]["Menu"], "Assembly")
        self.assertEqual(
            record["DbInfo"]["LinkList"][0]["Description"],
            "Assembly",
        )
        self.assertEqual(record["DbInfo"]["LinkList"][0]["DbTo"], "assembly")

        self.assertEqual(
            record["DbInfo"]["LinkList"][56]["Name"], "pubmed_taxonomy_entrez"
        )
        self.assertEqual(
            record["DbInfo"]["LinkList"][56]["Menu"], "Taxonomy via GenBank"
        )
        self.assertEqual(
            record["DbInfo"]["LinkList"][56]["Description"],
            "Related Taxonomy entry computed using other Entrez links",
        )
        self.assertEqual(record["DbInfo"]["LinkList"][56]["DbTo"], "taxonomy")

    def test_corrupted(self):
        """Test if corrupted XML is handled correctly."""
        # To create the XML file, use
        # >>> Bio.Entrez.einfo()
        # and manually delete the last couple of lines
        from Bio.Entrez import Parser

        with open("Entrez/einfo4.xml", "rb") as stream:
            self.assertRaises(Parser.CorruptedXMLError, Entrez.read, stream)


class ESearchTest(unittest.TestCase):
    """Tests for parsing XML output returned by ESearch."""

    def test_pubmed1(self):
        """Test parsing XML returned by ESearch from PubMed (first test)."""
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="biopython")
        with open("Entrez/esearch1.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["Count"], "63")
        self.assertEqual(record["RetMax"], "20")
        self.assertEqual(record["RetStart"], "0")
        self.assertEqual(len(record["IdList"]), 20)
        self.assertEqual(record["IdList"][0], "41282813")
        self.assertEqual(record["IdList"][1], "41148224")
        self.assertEqual(record["IdList"][2], "41011574")
        self.assertEqual(record["IdList"][3], "40959146")
        self.assertEqual(record["IdList"][4], "40937394")
        self.assertEqual(record["IdList"][5], "40657423")
        self.assertEqual(record["IdList"][6], "40651330")
        self.assertEqual(record["IdList"][7], "40572159")
        self.assertEqual(record["IdList"][8], "40160861")
        self.assertEqual(record["IdList"][9], "39883659")
        self.assertEqual(record["IdList"][10], "39882099")
        self.assertEqual(record["IdList"][11], "39717221")
        self.assertEqual(record["IdList"][12], "39546778")
        self.assertEqual(record["IdList"][13], "39507944")
        self.assertEqual(record["IdList"][14], "39445816")
        self.assertEqual(record["IdList"][15], "38808697")
        self.assertEqual(record["IdList"][16], "38650605")
        self.assertEqual(record["IdList"][17], "38365590")
        self.assertEqual(record["IdList"][18], "38235175")
        self.assertEqual(record["IdList"][19], "37810457")
        self.assertEqual(len(record["TranslationSet"]), 0)
        self.assertEqual(record["QueryTranslation"], '"biopython"[All Fields]')

    def test_pubmed2(self):
        """Test parsing XML returned by ESearch from PubMed (second test)."""
        # Search in PubMed for the term cancer for the entrez date from
        # the last 60 days and retrieve the first 100 IDs and translations
        # using the history parameter.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="cancer", reldate=60,
        #                        datetype="edat", retmax=100, usehistory="y")
        with open("Entrez/esearch2.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["Count"], "42249")
        self.assertEqual(record["RetMax"], "100")
        self.assertEqual(record["RetStart"], "0")
        self.assertEqual(record["QueryKey"], "1")
        self.assertEqual(record["WebEnv"], "MCID_6927d6e7fee3e90f880ec190")
        self.assertEqual(len(record["IdList"]), 100)
        self.assertEqual(record["IdList"][0], "41297076")
        self.assertEqual(record["IdList"][1], "41297074")
        self.assertEqual(record["IdList"][2], "41297068")
        self.assertEqual(record["IdList"][3], "41297055")
        self.assertEqual(record["IdList"][4], "41297046")
        self.assertEqual(record["IdList"][5], "41297043")
        self.assertEqual(record["IdList"][6], "41297032")
        self.assertEqual(record["IdList"][7], "41297028")
        self.assertEqual(record["IdList"][8], "41297023")
        self.assertEqual(record["IdList"][9], "41297015")
        self.assertEqual(record["IdList"][10], "41297014")
        self.assertEqual(record["IdList"][11], "41296995")
        self.assertEqual(record["IdList"][12], "41296994")
        self.assertEqual(record["IdList"][13], "41296991")
        self.assertEqual(record["IdList"][14], "41296990")
        self.assertEqual(record["IdList"][15], "41296985")
        self.assertEqual(record["IdList"][16], "41296927")
        self.assertEqual(record["IdList"][17], "41296920")
        self.assertEqual(record["IdList"][18], "41296918")
        self.assertEqual(record["IdList"][19], "41296886")
        self.assertEqual(record["IdList"][20], "41296881")
        self.assertEqual(record["IdList"][21], "41296866")
        self.assertEqual(record["IdList"][22], "41296862")
        self.assertEqual(record["IdList"][23], "41296854")
        self.assertEqual(record["IdList"][24], "41296842")
        self.assertEqual(record["IdList"][25], "41296829")
        self.assertEqual(record["IdList"][26], "41296826")
        self.assertEqual(record["IdList"][27], "41296824")
        self.assertEqual(record["IdList"][28], "41296822")
        self.assertEqual(record["IdList"][29], "41296821")
        self.assertEqual(record["IdList"][30], "41296811")
        self.assertEqual(record["IdList"][31], "41296783")
        self.assertEqual(record["IdList"][32], "41296773")
        self.assertEqual(record["IdList"][33], "41296763")
        self.assertEqual(record["IdList"][34], "41296753")
        self.assertEqual(record["IdList"][35], "41296747")
        self.assertEqual(record["IdList"][36], "41296732")
        self.assertEqual(record["IdList"][37], "41296731")
        self.assertEqual(record["IdList"][38], "41296714")
        self.assertEqual(record["IdList"][39], "41296711")
        self.assertEqual(record["IdList"][40], "41296710")
        self.assertEqual(record["IdList"][41], "41296689")
        self.assertEqual(record["IdList"][42], "41296685")
        self.assertEqual(record["IdList"][43], "41296683")
        self.assertEqual(record["IdList"][44], "41296675")
        self.assertEqual(record["IdList"][45], "41296667")
        self.assertEqual(record["IdList"][46], "41296665")
        self.assertEqual(record["IdList"][47], "41296663")
        self.assertEqual(record["IdList"][48], "41296658")
        self.assertEqual(record["IdList"][49], "41296656")
        self.assertEqual(record["IdList"][50], "41296616")
        self.assertEqual(record["IdList"][51], "41296609")
        self.assertEqual(record["IdList"][52], "41296608")
        self.assertEqual(record["IdList"][53], "41296606")
        self.assertEqual(record["IdList"][54], "41296594")
        self.assertEqual(record["IdList"][55], "41296590")
        self.assertEqual(record["IdList"][56], "41296586")
        self.assertEqual(record["IdList"][57], "41296583")
        self.assertEqual(record["IdList"][58], "41296580")
        self.assertEqual(record["IdList"][59], "41296569")
        self.assertEqual(record["IdList"][60], "41296567")
        self.assertEqual(record["IdList"][61], "41296563")
        self.assertEqual(record["IdList"][62], "41296561")
        self.assertEqual(record["IdList"][63], "41296555")
        self.assertEqual(record["IdList"][64], "41296544")
        self.assertEqual(record["IdList"][65], "41296517")
        self.assertEqual(record["IdList"][66], "41296514")
        self.assertEqual(record["IdList"][67], "41296508")
        self.assertEqual(record["IdList"][68], "41296505")
        self.assertEqual(record["IdList"][69], "41296501")
        self.assertEqual(record["IdList"][70], "41296494")
        self.assertEqual(record["IdList"][71], "41296485")
        self.assertEqual(record["IdList"][72], "41296462")
        self.assertEqual(record["IdList"][73], "41296458")
        self.assertEqual(record["IdList"][74], "41296454")
        self.assertEqual(record["IdList"][75], "41296452")
        self.assertEqual(record["IdList"][76], "41296443")
        self.assertEqual(record["IdList"][77], "41296441")
        self.assertEqual(record["IdList"][78], "41296440")
        self.assertEqual(record["IdList"][79], "41296437")
        self.assertEqual(record["IdList"][80], "41296436")
        self.assertEqual(record["IdList"][81], "41296434")
        self.assertEqual(record["IdList"][82], "41296433")
        self.assertEqual(record["IdList"][83], "41296432")
        self.assertEqual(record["IdList"][84], "41296430")
        self.assertEqual(record["IdList"][85], "41296429")
        self.assertEqual(record["IdList"][86], "41296425")
        self.assertEqual(record["IdList"][87], "41296417")
        self.assertEqual(record["IdList"][88], "41296413")
        self.assertEqual(record["IdList"][89], "41296412")
        self.assertEqual(record["IdList"][90], "41296408")
        self.assertEqual(record["IdList"][91], "41296406")
        self.assertEqual(record["IdList"][92], "41296402")
        self.assertEqual(record["IdList"][93], "41296396")
        self.assertEqual(record["IdList"][94], "41296387")
        self.assertEqual(record["IdList"][95], "41296385")
        self.assertEqual(record["IdList"][96], "41296382")
        self.assertEqual(record["IdList"][97], "41296379")
        self.assertEqual(record["IdList"][98], "41296369")
        self.assertEqual(record["IdList"][99], "41296368")
        self.assertEqual(len(record["TranslationSet"]), 1)
        self.assertEqual(record["TranslationSet"][0]["From"], "cancer")
        self.assertEqual(
            record["TranslationSet"][0]["To"],
            '"cancer\'s"[All Fields] OR "cancerated"[All Fields] OR "canceration"[All Fields] OR "cancerization"[All Fields] OR "cancerized"[All Fields] OR "cancerous"[All Fields] OR "neoplasms"[MeSH Terms] OR "neoplasms"[All Fields] OR "cancer"[All Fields] OR "cancers"[All Fields]',
        )
        self.assertEqual(
            record["QueryTranslation"],
            '("cancer s"[All Fields] OR "cancerated"[All Fields] OR "canceration"[All Fields] OR "cancerization"[All Fields] OR "cancerized"[All Fields] OR "cancerous"[All Fields] OR "neoplasms"[MeSH Terms] OR "neoplasms"[All Fields] OR "cancer"[All Fields] OR "cancers"[All Fields]) AND 2025/09/27:2025/11/26[Date - Entry]',
        )

    def test_pubmed3(self):
        """Test parsing XML returned by ESearch from PubMed (third test)."""
        # Search in PubMed for the journal PNAS Volume 97, and retrieve
        # 6 IDs starting at ID 7.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="PNAS[ta] AND 97[vi]",
        #                        retstart=6, retmax=6)
        with open("Entrez/esearch3.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["Count"], "2651")
        self.assertEqual(record["RetMax"], "6")
        self.assertEqual(record["RetStart"], "6")
        self.assertEqual(len(record["IdList"]), 6)
        self.assertEqual(record["IdList"][0], "11121077")
        self.assertEqual(record["IdList"][1], "11121076")
        self.assertEqual(record["IdList"][2], "11121075")
        self.assertEqual(record["IdList"][3], "11121074")
        self.assertEqual(record["IdList"][4], "11121073")
        self.assertEqual(record["IdList"][5], "11121072")
        self.assertEqual(len(record["TranslationSet"]), 1)
        self.assertEqual(record["TranslationSet"][0]["From"], "PNAS[ta]")
        self.assertEqual(
            record["TranslationSet"][0]["To"],
            '"Proc Natl Acad Sci U S A"[Journal:__jid7505876]',
        )
        self.assertEqual(
            record["QueryTranslation"],
            '"proc natl acad sci u s a"[Journal] AND "97"[Volume]',
        )

    def test_pmc(self):
        """Test parsing XML returned by ESearch from PubMed Central."""
        # Search in PubMed Central for stem cells in articles with an abstract.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pmc", term="stem cells AND hasabstract")
        with open("Entrez/esearch5.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["Count"], "5")
        self.assertEqual(record["RetMax"], "5")
        self.assertEqual(record["RetStart"], "0")
        self.assertEqual(len(record["IdList"]), 5)
        self.assertEqual(record["IdList"][0], "6540636")
        self.assertEqual(record["IdList"][1], "5830998")
        self.assertEqual(record["IdList"][2], "4460425")
        self.assertEqual(record["IdList"][3], "4466840")
        self.assertEqual(record["IdList"][4], "3443459")
        self.assertEqual(len(record["TranslationSet"]), 1)
        self.assertEqual(record["TranslationSet"][0]["From"], "stem cells")
        self.assertEqual(
            record["TranslationSet"][0]["To"],
            '"stem cells"[MeSH Terms] OR ("stem"[All Fields] AND "cells"[All Fields]) OR "stem cells"[All Fields]',
        )
        self.assertEqual(len(record["TranslationStack"]), 11)
        self.assertEqual(
            record["TranslationStack"][0]["Term"], '"stem cells"[MeSH Terms]'
        )
        self.assertEqual(record["TranslationStack"][0]["Field"], "MeSH Terms")
        self.assertEqual(record["TranslationStack"][0]["Count"], "109831")
        self.assertEqual(record["TranslationStack"][0]["Explode"], "Y")
        self.assertEqual(record["TranslationStack"][0].tag, "TermSet")
        self.assertEqual(record["TranslationStack"][1]["Term"], '"stem"[All Fields]')
        self.assertEqual(record["TranslationStack"][1]["Field"], "All Fields")
        self.assertEqual(record["TranslationStack"][1]["Count"], "1528859")
        self.assertEqual(record["TranslationStack"][1]["Explode"], "N")
        self.assertEqual(record["TranslationStack"][1].tag, "TermSet")
        self.assertEqual(
            record["TranslationStack"][2],
            {
                "Term": '"cells"[All Fields]',
                "Field": "All Fields",
                "Count": "5033502",
                "Explode": "N",
            },
        )
        self.assertEqual(record["TranslationStack"][2].tag, "TermSet")
        self.assertEqual(record["TranslationStack"][3], "AND")
        self.assertEqual(record["TranslationStack"][3].tag, "OP")
        self.assertEqual(record["TranslationStack"][4], "GROUP")
        self.assertEqual(record["TranslationStack"][4].tag, "OP")
        self.assertEqual(record["TranslationStack"][5], "OR")
        self.assertEqual(record["TranslationStack"][5].tag, "OP")
        self.assertEqual(
            record["TranslationStack"][6],
            {
                "Term": '"stem cells"[All Fields]',
                "Field": "All Fields",
                "Count": "715386",
                "Explode": "N",
            },
        )
        self.assertEqual(record["TranslationStack"][6].tag, "TermSet")
        self.assertEqual(record["TranslationStack"][7], "OR")
        self.assertEqual(record["TranslationStack"][7].tag, "OP")
        self.assertEqual(record["TranslationStack"][8], "GROUP")
        self.assertEqual(record["TranslationStack"][8].tag, "OP")
        self.assertEqual(
            record["TranslationStack"][9]["Term"], "hasabstract[All Fields]"
        )
        self.assertEqual(record["TranslationStack"][9]["Field"], "All Fields")
        self.assertEqual(record["TranslationStack"][9]["Count"], "83")
        self.assertEqual(record["TranslationStack"][9]["Explode"], "N")
        self.assertEqual(record["TranslationStack"][9].tag, "TermSet")
        self.assertEqual(record["TranslationStack"][10], "AND")
        self.assertEqual(record["TranslationStack"][10].tag, "OP")
        self.assertEqual(
            record["QueryTranslation"],
            '("stem cells"[MeSH Terms] OR ("stem"[All Fields] AND "cells"[All Fields]) OR "stem cells"[All Fields]) AND hasabstract[All Fields]',
        )

    def test_nucleotide(self):
        """Test parsing XML returned by ESearch from the Nucleotide database."""
        # Search in Nucleotide for a property of the sequence,
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="nucleotide", term="biomol trna[prop]")
        with open("Entrez/esearch6.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["Count"], "1018")
        self.assertEqual(record["RetMax"], "20")
        self.assertEqual(record["RetStart"], "0")
        self.assertEqual(len(record["IdList"]), 20)
        self.assertEqual(record["IdList"][0], "3069587804")
        self.assertEqual(record["IdList"][1], "3069587774")
        self.assertEqual(record["IdList"][2], "3069587771")
        self.assertEqual(record["IdList"][3], "3069587741")
        self.assertEqual(record["IdList"][4], "3069402178")
        self.assertEqual(record["IdList"][5], "3069218576")
        self.assertEqual(record["IdList"][6], "2737963026")
        self.assertEqual(record["IdList"][7], "2586967820")
        self.assertEqual(record["IdList"][8], "2274792564")
        self.assertEqual(record["IdList"][9], "2274792563")
        self.assertEqual(record["IdList"][10], "2274792562")
        self.assertEqual(record["IdList"][11], "2274792561")
        self.assertEqual(record["IdList"][12], "2274792560")
        self.assertEqual(record["IdList"][13], "2274792559")
        self.assertEqual(record["IdList"][14], "2274792558")
        self.assertEqual(record["IdList"][15], "2274792557")
        self.assertEqual(record["IdList"][16], "2274792556")
        self.assertEqual(record["IdList"][17], "2274792555")
        self.assertEqual(record["IdList"][18], "2274792554")
        self.assertEqual(record["IdList"][19], "2274792553")
        self.assertEqual(len(record["TranslationSet"]), 0)
        self.assertEqual(record["QueryTranslation"], "biomol trna[prop]")

    def test_protein(self):
        """Test parsing XML returned by ESearch from the Protein database."""
        # Search in Protein for a molecular weight
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="protein", term="200020[molecular weight]")
        with open("Entrez/esearch7.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["Count"], "330")
        self.assertEqual(record["RetMax"], "20")
        self.assertEqual(record["RetStart"], "0")
        self.assertEqual(len(record["IdList"]), 20)
        self.assertEqual(record["IdList"][0], "3108048616")
        self.assertEqual(record["IdList"][1], "3076096315")
        self.assertEqual(record["IdList"][2], "3070972985")
        self.assertEqual(record["IdList"][3], "3070963398")
        self.assertEqual(record["IdList"][4], "3042783758")
        self.assertEqual(record["IdList"][5], "3042778745")
        self.assertEqual(record["IdList"][6], "3042773739")
        self.assertEqual(record["IdList"][7], "3042768843")
        self.assertEqual(record["IdList"][8], "3042764246")
        self.assertEqual(record["IdList"][9], "3042754755")
        self.assertEqual(record["IdList"][10], "3042750053")
        self.assertEqual(record["IdList"][11], "3037813595")
        self.assertEqual(record["IdList"][12], "3032309824")
        self.assertEqual(record["IdList"][13], "1468867753")
        self.assertEqual(record["IdList"][14], "2994907531")
        self.assertEqual(record["IdList"][15], "2994891115")
        self.assertEqual(record["IdList"][16], "2993271659")
        self.assertEqual(record["IdList"][17], "2993266811")
        self.assertEqual(record["IdList"][18], "2993262038")
        self.assertEqual(record["IdList"][19], "2174081062")
        self.assertEqual(len(record["TranslationSet"]), 0)
        self.assertEqual(len(record["TranslationStack"]), 2)
        self.assertEqual(
            record["TranslationStack"][0]["Term"], "000200020[molecular weight]"
        )
        self.assertEqual(record["TranslationStack"][0]["Field"], "molecular weight")
        self.assertEqual(record["TranslationStack"][0]["Count"], "330")
        self.assertEqual(record["TranslationStack"][0]["Explode"], "N")
        self.assertEqual(record["TranslationStack"][0].tag, "TermSet")
        self.assertEqual(record["TranslationStack"][1], "GROUP")
        self.assertEqual(record["TranslationStack"][1].tag, "OP")
        self.assertEqual(record["QueryTranslation"], "000200020[molecular weight]")

    def test_notfound(self):
        """Test parsing XML returned by ESearch when no items were found."""
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="protein", term="abcXYZ")
        with open("Entrez/esearch8.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["Count"], "0")
        self.assertEqual(record["RetMax"], "0")
        self.assertEqual(record["RetStart"], "0")
        self.assertEqual(len(record["IdList"]), 0)
        self.assertEqual(len(record["TranslationSet"]), 0)
        self.assertEqual(record["QueryTranslation"], "(abcXYZ[All Fields])")
        self.assertEqual(len(record["ErrorList"]), 2)
        self.assertIn("PhraseNotFound", record["ErrorList"])
        self.assertIn("FieldNotFound", record["ErrorList"])
        self.assertEqual(len(record["ErrorList"]["PhraseNotFound"]), 1)
        self.assertEqual(len(record["ErrorList"]["FieldNotFound"]), 0)
        self.assertEqual(record["ErrorList"]["PhraseNotFound"][0], "abcXYZ")
        self.assertEqual(len(record["WarningList"]), 3)
        self.assertIn("PhraseIgnored", record["WarningList"])
        self.assertIn("QuotedPhraseNotFound", record["WarningList"])
        self.assertIn("OutputMessage", record["WarningList"])
        self.assertEqual(len(record["WarningList"]["PhraseIgnored"]), 0)
        self.assertEqual(len(record["WarningList"]["QuotedPhraseNotFound"]), 0)
        self.assertEqual(len(record["WarningList"]["OutputMessage"]), 1)
        self.assertEqual(record["WarningList"]["OutputMessage"][0], "No items found.")


class EPostTest(unittest.TestCase):
    """Tests for parsing XML output returned by EPost."""

    # Don't know how to get an InvalidIdList in the XML returned by EPost;
    # unable to test if we are parsing it correctly.
    def test_epost(self):
        """Test parsing XML returned by EPost."""
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="pubmed", id="11237011")
        with open("Entrez/epost1.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["QueryKey"], "1")
        self.assertEqual(record["WebEnv"], "MCID_692851c130aec64bed0a8c2a")

    def test_wrong(self):
        """Test parsing XML returned by EPost with incorrect arguments."""
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="nothing")
        with open("Entrez/epost2.xml", "rb") as stream:
            self.assertRaises(RuntimeError, Entrez.read, stream)
        with open("Entrez/epost2.xml", "rb") as stream:
            record = Entrez.read(stream, ignore_errors=True)
        self.assertEqual(len(record), 1)
        self.assertEqual(len(record.attributes), 0)
        self.assertEqual(record["ERROR"], "Invalid db name specified: nothing")
        self.assertEqual(record["ERROR"].tag, "ERROR")

    def test_invalid(self):
        """Test parsing XML returned by EPost with invalid id (overflow tag)."""
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="pubmed", id=99999999999999999999999999999999)
        with open("Entrez/epost3.xml", "rb") as stream:
            self.assertRaises(RuntimeError, Entrez.read, stream)
        with open("Entrez/epost3.xml", "rb") as stream:
            record = Entrez.read(stream, ignore_errors=True)
        self.assertEqual(len(record), 1)
        self.assertEqual(len(record.attributes), 0)
        self.assertEqual(
            record["ERROR"],
            "Some IDs have invalid value and were omitted. Maximum ID value 18446744073709551615",
        )
        self.assertEqual(record["ERROR"].tag, "ERROR")
        # Note that the first ERROR element is lost. Strictly speaking, the XML
        # is not consistent with the DTD, which allows only one ERROR element.


class ESummaryTest(unittest.TestCase):
    """Tests for parsing XML output returned by ESummary."""

    # Items have a type, which can be
    # (Integer|Date|String|Structure|List|Flags|Qualifier|Enumerator|Unknown)
    # I don't have an XML file where the type "Flags", "Qualifier",
    # "Enumerator", or "Unknown" is used, so they are not tested here.
    def test_pubmed(self):
        """Test parsing XML returned by ESummary from PubMed."""
        # In PubMed display records for PMIDs 11850928 and 11482001 in
        # xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="pubmed", id=["11850928","11482001"],
        #                         retmode="xml")
        with open("Entrez/esummary1.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record[0]["Id"], "11850928")
        self.assertEqual(record[0]["PubDate"], "1965 Aug")
        self.assertEqual(record[0]["EPubDate"], "")
        self.assertEqual(record[0]["Source"], "Arch Dermatol")
        self.assertEqual(len(record[0]["AuthorList"]), 2)
        self.assertEqual(record[0]["AuthorList"][0], "LoPresti PJ")
        self.assertEqual(record[0]["AuthorList"][1], "Hambrick GW Jr")
        self.assertEqual(record[0]["LastAuthor"], "Hambrick GW Jr")
        self.assertEqual(
            record[0]["Title"],
            "Zirconium granuloma following treatment of rhus dermatitis.",
        )
        self.assertEqual(record[0]["Volume"], "92")
        self.assertEqual(record[0]["Issue"], "2")
        self.assertEqual(record[0]["Pages"], "188-91")
        self.assertEqual(record[0]["LangList"], ["English"])
        self.assertEqual(record[0]["NlmUniqueID"], "0372433")
        self.assertEqual(record[0]["ISSN"], "0003-987X")
        self.assertEqual(record[0]["ESSN"], "")
        self.assertEqual(len(record[0]["PubTypeList"]), 1)
        self.assertEqual(record[0]["PubTypeList"][0], "Journal Article")
        self.assertEqual(record[0]["RecordStatus"], "PubMed - indexed for MEDLINE")
        self.assertEqual(record[0]["PubStatus"], "ppublish")
        self.assertEqual(len(record[0]["ArticleIds"]), 2)
        self.assertEqual(record[0]["ArticleIds"]["pubmed"], ["11850928"])
        self.assertEqual(record[0]["ArticleIds"]["medline"], [])
        self.assertEqual(len(record[0]["History"]), 3)
        self.assertEqual(record[0]["History"]["pubmed"], ["1965/08/01 00:00"])
        self.assertEqual(record[0]["History"]["medline"], ["2002/03/09 10:01"])
        self.assertEqual(record[0]["History"]["entrez"], "1965/08/01 00:00")
        self.assertEqual(len(record[0]["References"]), 0)
        self.assertEqual(record[0]["HasAbstract"], 1)
        self.assertEqual(record[0]["PmcRefCount"], 0)
        self.assertEqual(record[0]["FullJournalName"], "Archives of dermatology")
        self.assertEqual(record[0]["ELocationID"], "")
        self.assertEqual(record[0]["SO"], "1965 Aug;92(2):188-91")

        self.assertEqual(record[1]["Id"], "11482001")
        self.assertEqual(record[1]["PubDate"], "2001 Jun")
        self.assertEqual(record[1]["EPubDate"], "")
        self.assertEqual(record[1]["Source"], "Adverse Drug React Toxicol Rev")
        self.assertEqual(len(record[1]["AuthorList"]), 3)
        self.assertEqual(record[1]["AuthorList"][0], "Mantle D")
        self.assertEqual(record[1]["AuthorList"][1], "Gok MA")
        self.assertEqual(record[1]["AuthorList"][2], "Lennard TW")
        self.assertEqual(record[1]["LastAuthor"], "Lennard TW")
        self.assertEqual(
            record[1]["Title"],
            "Adverse and beneficial effects of plant extracts on skin and skin disorders.",
        )
        self.assertEqual(record[1]["Volume"], "20")
        self.assertEqual(record[1]["Issue"], "2")
        self.assertEqual(record[1]["Pages"], "89-103")
        self.assertEqual(len(record[1]["LangList"]), 1)
        self.assertEqual(record[1]["LangList"][0], "English")
        self.assertEqual(record[1]["NlmUniqueID"], "9109474")
        self.assertEqual(record[1]["ISSN"], "0964-198X")
        self.assertEqual(record[1]["ESSN"], "")
        self.assertEqual(len(record[1]["PubTypeList"]), 2)
        self.assertEqual(record[1]["PubTypeList"][0], "Journal Article")
        self.assertEqual(record[1]["PubTypeList"][1], "Review")
        self.assertEqual(record[1]["RecordStatus"], "PubMed - indexed for MEDLINE")
        self.assertEqual(record[1]["PubStatus"], "ppublish")
        self.assertEqual(len(record[1]["ArticleIds"]), 2)
        self.assertEqual(record[1]["ArticleIds"]["pubmed"], ["11482001"])
        self.assertEqual(record[1]["ArticleIds"]["medline"], [])
        self.assertEqual(len(record[1]["History"]), 3)
        self.assertEqual(record[1]["History"]["pubmed"], ["2001/08/03 10:00"])
        self.assertEqual(record[1]["History"]["medline"], ["2002/01/23 10:01"])
        self.assertEqual(record[1]["History"]["entrez"], "2001/08/03 10:00")
        self.assertEqual(len(record[1]["References"]), 0)
        self.assertEqual(record[1]["HasAbstract"], 1)
        self.assertEqual(record[1]["PmcRefCount"], 0)
        self.assertEqual(
            record[1]["FullJournalName"],
            "Adverse drug reactions and toxicological reviews",
        )
        self.assertEqual(record[1]["ELocationID"], "")
        self.assertEqual(record[1]["SO"], "2001 Jun;20(2):89-103")

    def test_protein(self):
        """Test parsing XML returned by ESummary from the Protein database."""
        # In Protein display records for GIs 28800982 and 28628843 in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="protein", id="28800982,28628843", retmode="xml")
        with open("Entrez/esummary3.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record[0]["Id"], "28800982")
        self.assertEqual(record[0]["Caption"], "AAO47091")
        self.assertEqual(record[0]["Title"], "hemochromatosis, partial [Homo sapiens]")
        self.assertEqual(record[0]["Extra"], "gi|28800982|gb|AAO47091.1|[28800982]")
        self.assertEqual(record[0]["Gi"], 28800982)
        self.assertEqual(record[0]["CreateDate"], "2003/03/03")
        self.assertEqual(record[0]["UpdateDate"], "2016/07/25")
        self.assertEqual(record[0]["Flags"], 0)
        self.assertEqual(record[0]["TaxId"], 9606)
        self.assertEqual(record[0]["Length"], 268)
        self.assertEqual(record[0]["Status"], "live")
        self.assertEqual(record[0]["ReplacedBy"], "")
        self.assertEqual(record[0]["Comment"], "  ")

        self.assertEqual(record[1]["Id"], "28628843")
        self.assertEqual(record[1]["Caption"], "AAO49381")
        self.assertEqual(
            record[1]["Title"], "erythroid associated factor [Homo sapiens]"
        )
        self.assertEqual(
            record[1]["Extra"], "gi|28628843|gb|AAO49381.1|AF485325_1[28628843]"
        )
        self.assertEqual(record[1]["Gi"], 28628843)
        self.assertEqual(record[1]["CreateDate"], "2003/03/02")
        self.assertEqual(record[1]["UpdateDate"], "2003/03/02")
        self.assertEqual(record[1]["Flags"], 0)
        self.assertEqual(record[1]["TaxId"], 9606)
        self.assertEqual(record[1]["Length"], 102)
        self.assertEqual(record[1]["Status"], "live")
        self.assertEqual(record[1]["ReplacedBy"], "")
        self.assertEqual(record[1]["Comment"], "  ")

    def test_nucleotide(self):
        """Test parsing XML returned by ESummary from the Nucleotide database."""
        # In Nucleotide display records for GIs 28864546 and 28800981
        # in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="nucleotide", id="28864546,28800981",
        #                         retmode="xml")
        with open("Entrez/esummary4.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record[0]["Id"], "28864546")
        self.assertEqual(record[0]["Caption"], "AY207443")
        self.assertEqual(
            record[0]["Title"],
            "Homo sapiens alpha hemoglobin (HBZP) pseudogene 3' UTR/AluJo repeat breakpoint junction",
        )
        self.assertEqual(record[0]["Extra"], "gi|28864546|gb|AY207443.1|[28864546]")
        self.assertEqual(record[0]["Gi"], 28864546)
        self.assertEqual(record[0]["CreateDate"], "2003/03/05")
        self.assertEqual(record[0]["UpdateDate"], "2003/03/05")
        self.assertEqual(record[0]["Flags"], 0)
        self.assertEqual(record[0]["TaxId"], 9606)
        self.assertEqual(record[0]["Length"], 491)
        self.assertEqual(record[0]["Status"], "live")
        self.assertEqual(record[0]["ReplacedBy"], "")
        self.assertEqual(record[0]["Comment"], "  ")

        self.assertEqual(record[1]["Id"], "28800981")
        self.assertEqual(record[1]["Caption"], "AY205604")
        self.assertEqual(
            record[1]["Title"], "Homo sapiens hemochromatosis (HFE) mRNA, partial cds"
        )
        self.assertEqual(record[1]["Extra"], "gi|28800981|gb|AY205604.1|[28800981]")
        self.assertEqual(record[1]["Gi"], 28800981)
        self.assertEqual(record[1]["CreateDate"], "2003/03/03")
        self.assertEqual(record[1]["UpdateDate"], "2016/07/25")
        self.assertEqual(record[1]["Flags"], 0)
        self.assertEqual(record[1]["TaxId"], 9606)
        self.assertEqual(record[1]["Length"], 860)
        self.assertEqual(record[1]["Status"], "live")
        self.assertEqual(record[1]["ReplacedBy"], "")
        self.assertEqual(record[1]["Comment"], "  ")

    def test_structure(self):
        """Test parsing XML returned by ESummary from the Structure database."""
        # In Nucleotide display records for GIs 28864546 and 28800981
        # in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="structure", id=["19923","12120"],
        #                         retmode="xml")
        with open("Entrez/esummary5.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record[0]["Id"], "19923")
        self.assertEqual(record[0]["PdbAcc"], "1L5J")
        self.assertEqual(
            record[0]["PdbDescr"], "CRYSTAL STRUCTURE OF E. COLI ACONITASE B"
        )
        self.assertEqual(record[0]["EC"], "")
        self.assertEqual(record[0]["Resolution"], "2.4")
        self.assertEqual(record[0]["ExpMethod"], "X-ray Diffraction")
        self.assertEqual(record[0]["PdbClass"], "LYASE")
        self.assertEqual(record[0]["PdbDepositDate"], "2002/03/07 00:00")
        self.assertEqual(record[0]["MMDBEntryDate"], "2002/07/11 00:00")
        self.assertEqual(record[0]["OrganismList"], ["Escherichia coli"])
        self.assertEqual(record[0]["LigCode"], "F3S|TRA")
        self.assertEqual(record[0]["LigCount"], "2")
        self.assertEqual(record[0]["ModProteinResCount"], "0")
        self.assertEqual(record[0]["ModDNAResCount"], "0")
        self.assertEqual(record[0]["ModRNAResCount"], "0")
        self.assertEqual(record[0]["ProteinChainCount"], "")
        self.assertEqual(record[0]["DNAChainCount"], "")
        self.assertEqual(record[0]["RNAChainCount"], "")

        self.assertEqual(record[1]["Id"], "12120")
        self.assertEqual(record[1]["PdbAcc"], "1B0K")
        self.assertEqual(
            record[1]["PdbDescr"], "S642A:FLUOROCITRATE COMPLEX OF ACONITASE"
        )
        self.assertEqual(record[1]["EC"], "")
        self.assertEqual(record[1]["Resolution"], "2.5")
        self.assertEqual(record[1]["ExpMethod"], "X-ray Diffraction")
        self.assertEqual(record[1]["PdbClass"], "LYASE")
        self.assertEqual(record[1]["PdbDepositDate"], "1998/11/11 00:00")
        self.assertEqual(record[1]["MMDBEntryDate"], "2000/01/24 00:00")
        self.assertEqual(record[1]["OrganismList"], ["Sus scrofa"])
        self.assertEqual(record[1]["LigCode"], "FLC|O|SF4")
        self.assertEqual(record[1]["LigCount"], "3")
        self.assertEqual(record[1]["ModProteinResCount"], "0")
        self.assertEqual(record[1]["ModDNAResCount"], "0")
        self.assertEqual(record[1]["ModRNAResCount"], "0")
        self.assertEqual(record[1]["ProteinChainCount"], "")
        self.assertEqual(record[1]["DNAChainCount"], "")
        self.assertEqual(record[1]["RNAChainCount"], "")

    def test_taxonomy(self):
        """Test parsing XML returned by ESummary from the Taxonomy database."""
        # In Taxonomy display records for TAXIDs 9913 and 30521 in
        # xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="taxonomy", id=["9913","30521"],
        #                         retmode="xml")
        with open("Entrez/esummary6.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record[0]["Id"], "9913")
        self.assertEqual(record[0]["Status"], "active")
        self.assertEqual(record[0]["Rank"], "species")
        self.assertEqual(record[0]["Division"], "even-toed ungulates & whales")
        self.assertEqual(record[0]["ScientificName"], "Bos taurus")
        self.assertEqual(record[0]["CommonName"], "domestic cattle")
        self.assertEqual(record[0]["TaxId"], 9913)
        self.assertEqual(record[0]["AkaTaxId"], 0)
        self.assertEqual(record[0]["Genus"], "")
        self.assertEqual(record[0]["Species"], "")
        self.assertEqual(record[0]["Subsp"], "")
        self.assertEqual(record[0]["ModificationDate"], "2024/08/09 00:00")

        self.assertEqual(record[1]["Id"], "30521")
        self.assertEqual(record[1]["Status"], "active")
        self.assertEqual(record[1]["Rank"], "species")
        self.assertEqual(record[1]["Division"], "even-toed ungulates & whales")
        self.assertEqual(record[1]["ScientificName"], "Bos grunniens")
        self.assertEqual(record[1]["CommonName"], "domestic yak")
        self.assertEqual(record[1]["TaxId"], 30521)
        self.assertEqual(record[1]["AkaTaxId"], 0)
        self.assertEqual(record[1]["Genus"], "")
        self.assertEqual(record[1]["Species"], "")
        self.assertEqual(record[1]["Subsp"], "")
        self.assertEqual(record[1]["ModificationDate"], "2020/04/29 00:00")

    def test_wrong(self):
        """Test parsing XML returned by ESummary with incorrect arguments."""
        # To create the XML file, use
        # >>> Bio.Entrez.esummary()
        with open("Entrez/esummary8.xml", "rb") as stream:
            self.assertRaises(RuntimeError, Entrez.read, stream)
        with open("Entrez/esummary8.xml", "rb") as stream:
            record = Entrez.read(stream, ignore_errors=True)
        self.assertEqual(len(record), 1)
        self.assertEqual(len(record.attributes), 0)
        self.assertEqual(record[0], "Neither query_key nor id specified")
        self.assertEqual(record[0].tag, "ERROR")

    def test_integer_none(self):
        """Test parsing ESummary XML where an Integer is not defined."""
        # To create the XML file, use
        # >>> Entrez.esummary(db='pccompound', id='7488')
        with open("Entrez/esummary9.xml", "rb") as stream:
            records = Entrez.read(stream)
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record["Id"], "7488")
        self.assertEqual(record["CID"], 7488)
        self.assertEqual(record["SourceNameList"], [])
        self.assertEqual(record["SourceIDList"], [])
        self.assertEqual(len(record["SourceCategoryList"]), 8)
        self.assertEqual(record["SourceCategoryList"][0], "Chemical Vendors")
        self.assertEqual(record["SourceCategoryList"][1], "Research and Development")
        self.assertEqual(record["SourceCategoryList"][2], "Curation Efforts")
        self.assertEqual(record["SourceCategoryList"][3], "Governmental Organizations")
        self.assertEqual(record["SourceCategoryList"][4], "Legacy Depositors")
        self.assertEqual(record["SourceCategoryList"][5], "Subscription Services")
        self.assertEqual(record["SourceCategoryList"][6], "Journal Publishers")
        self.assertEqual(record["SourceCategoryList"][7], "NIH Initiatives")
        self.assertEqual(record["CreateDate"], "2005/03/26 00:00")
        self.assertEqual(len(record["SynonymList"]), 77)
        self.assertEqual(record["SynonymList"][0], "Terephthaloyl chloride")
        self.assertEqual(record["SynonymList"][1], "100-20-9")
        self.assertEqual(record["SynonymList"][2], "Terephthaloyl dichloride")
        self.assertEqual(record["SynonymList"][3], "1,4-BENZENEDICARBONYL DICHLORIDE")
        self.assertEqual(record["SynonymList"][4], "Terephthalic acid dichloride")
        self.assertEqual(record["SynonymList"][5], "Terephthalic dichloride")
        self.assertEqual(record["SynonymList"][6], "p-Phthaloyl chloride")
        self.assertEqual(record["SynonymList"][7], "Terephthalic acid chloride")
        self.assertEqual(record["SynonymList"][8], "p-Phthalyl dichloride")
        self.assertEqual(record["SynonymList"][9], "p-Phthaloyl dichloride")
        self.assertEqual(record["SynonymList"][10], "Terephthalyl dichloride")
        self.assertEqual(record["SynonymList"][11], "1,4-Benzenedicarbonyl chloride")
        self.assertEqual(record["SynonymList"][12], "p-Phenylenedicarbonyl dichloride")
        self.assertEqual(record["SynonymList"][13], "benzene-1,4-dicarbonyl chloride")
        self.assertEqual(record["SynonymList"][14], "NSC 41885")
        self.assertEqual(record["SynonymList"][15], "terephthaloylchloride")
        self.assertEqual(record["SynonymList"][16], "UNII-G247CO9608")
        self.assertEqual(record["SynonymList"][17], "HSDB 5332")
        self.assertEqual(record["SynonymList"][18], "EINECS 202-829-5")
        self.assertEqual(record["SynonymList"][19], "BRN 0607796")
        self.assertEqual(record["SynonymList"][20], "LXEJRKJRKIFVNY-UHFFFAOYSA-N")
        self.assertEqual(record["SynonymList"][21], "MFCD00000693")
        self.assertEqual(record["SynonymList"][22], "G247CO9608")
        self.assertEqual(record["SynonymList"][23], "DSSTox_CID_6653")
        self.assertEqual(record["SynonymList"][24], "DSSTox_RID_78175")
        self.assertEqual(record["SynonymList"][25], "DSSTox_GSID_26653")
        self.assertEqual(record["SynonymList"][26], "Q-201791")
        self.assertEqual(record["SynonymList"][27], "Terephthaloyl chloride, 99+%")
        self.assertEqual(record["SynonymList"][28], "CAS-100-20-9")
        self.assertEqual(record["SynonymList"][29], "CCRIS 8626")
        self.assertEqual(record["SynonymList"][30], "p-Phthalyl chloride")
        self.assertEqual(record["SynonymList"][31], "terephthalic chloride")
        self.assertEqual(record["SynonymList"][32], "tere-phthaloyl chloride")
        self.assertEqual(record["SynonymList"][33], "AC1L1OVG")
        self.assertEqual(record["SynonymList"][34], "ACMC-2097nf")
        self.assertEqual(record["SynonymList"][35], "EC 202-829-5")
        self.assertEqual(record["SynonymList"][36], "1,4-Dichloroformyl benzene")
        self.assertEqual(record["SynonymList"][37], "SCHEMBL68148")
        self.assertEqual(
            record["SynonymList"][38], "4-09-00-03318 (Beilstein Handbook Reference)"
        )
        self.assertEqual(record["SynonymList"][39], "KSC174E9T")
        self.assertEqual(record["SynonymList"][40], "CHEMBL1893301")
        self.assertEqual(record["SynonymList"][41], "DTXSID7026653")
        self.assertEqual(record["SynonymList"][42], "KS-00000VAD")
        self.assertEqual(record["SynonymList"][43], "benzene-1,4-dicarbonyl dichloride")
        self.assertEqual(record["SynonymList"][44], "MolPort-003-926-079")
        self.assertEqual(record["SynonymList"][45], "BCP27385")
        self.assertEqual(record["SynonymList"][46], "NSC41885")
        self.assertEqual(record["SynonymList"][47], "Tox21_201899")
        self.assertEqual(record["SynonymList"][48], "Tox21_303166")
        self.assertEqual(record["SynonymList"][49], "ANW-14185")
        self.assertEqual(record["SynonymList"][50], "NSC-41885")
        self.assertEqual(record["SynonymList"][51], "ZINC38141445")
        self.assertEqual(record["SynonymList"][52], "AKOS015890038")
        self.assertEqual(record["SynonymList"][53], "FCH1319904")
        self.assertEqual(record["SynonymList"][54], "MCULE-9481285116")
        self.assertEqual(record["SynonymList"][55], "RP25985")
        self.assertEqual(
            record["SynonymList"][56], "Terephthaloyl chloride, >=99%, flakes"
        )
        self.assertEqual(record["SynonymList"][57], "NCGC00164045-01")
        self.assertEqual(record["SynonymList"][58], "NCGC00164045-02")
        self.assertEqual(record["SynonymList"][59], "NCGC00257127-01")
        self.assertEqual(record["SynonymList"][60], "NCGC00259448-01")
        self.assertEqual(record["SynonymList"][61], "AN-24545")
        self.assertEqual(record["SynonymList"][62], "I764")
        self.assertEqual(record["SynonymList"][63], "KB-10499")
        self.assertEqual(record["SynonymList"][64], "OR315758")
        self.assertEqual(record["SynonymList"][65], "SC-19185")
        self.assertEqual(record["SynonymList"][66], "LS-148753")
        self.assertEqual(record["SynonymList"][67], "RT-000669")
        self.assertEqual(record["SynonymList"][68], "ST51037908")
        self.assertEqual(record["SynonymList"][69], "6804-EP1441224A2")
        self.assertEqual(record["SynonymList"][70], "I01-5090")
        self.assertEqual(
            record["SynonymList"][71],
            "InChI=1/C8H4Cl2O2/c9-7(11)5-1-2-6(4-3-5)8(10)12/h1-4",
        )
        self.assertEqual(record["SynonymList"][72], "106158-15-0")
        self.assertEqual(record["SynonymList"][73], "108454-76-8")
        self.assertEqual(record["SynonymList"][74], "1640987-72-9")
        self.assertEqual(record["SynonymList"][75], "188665-55-6")
        self.assertEqual(record["SynonymList"][76], "1927884-58-9")
        self.assertEqual(len(record["MeSHHeadingList"]), 1)
        self.assertEqual(record["MeSHHeadingList"][0], "terephthaloyl chloride")
        self.assertEqual(len(record["MeSHTermList"]), 5)
        self.assertEqual(record["MeSHTermList"][0], "p-phthaloyl dichloride")
        self.assertEqual(record["MeSHTermList"][1], "terephthaloyl dichloride")
        self.assertEqual(record["MeSHTermList"][2], "1,4-benzenedicarbonyl dichloride")
        self.assertEqual(record["MeSHTermList"][3], "1,4-phthaloyl dichloride")
        self.assertEqual(record["MeSHTermList"][4], "terephthaloyl chloride")
        self.assertEqual(len(record["PharmActionList"]), 0)
        self.assertEqual(record["CommentList"], [])
        self.assertEqual(record["IUPACName"], "benzene-1,4-dicarbonyl chloride")
        self.assertEqual(record["CanonicalSmiles"], "C1=CC(=CC=C1C(=O)Cl)C(=O)Cl")
        self.assertEqual(record["IsomericSmiles"], "C1=CC(=CC=C1C(=O)Cl)C(=O)Cl")
        self.assertEqual(record["RotatableBondCount"], 2)
        self.assertEqual(record["MolecularFormula"], "C8H4Cl2O2")
        self.assertEqual(record["MolecularWeight"], "203.018")
        self.assertEqual(record["TotalFormalCharge"], 0)
        self.assertEqual(record["XLogP"], "4")
        self.assertEqual(record["HydrogenBondDonorCount"], 0)
        self.assertEqual(record["HydrogenBondAcceptorCount"], 2)
        self.assertEqual(record["Complexity"], "173.000")
        self.assertEqual(record["HeavyAtomCount"], 12)
        self.assertEqual(record["AtomChiralCount"], 0)
        self.assertEqual(record["AtomChiralDefCount"], 0)
        self.assertEqual(record["AtomChiralUndefCount"], 0)
        self.assertEqual(record["BondChiralCount"], 0)
        self.assertEqual(record["BondChiralDefCount"], 0)
        self.assertEqual(record["BondChiralUndefCount"], 0)
        self.assertEqual(record["IsotopeAtomCount"], 0)
        self.assertEqual(record["CovalentUnitCount"], 1)
        self.assertEqual(record["TautomerCount"], None)  # noqa: A502
        self.assertEqual(record["SubstanceIDList"], [])
        self.assertEqual(record["TPSA"], "34.1")
        self.assertEqual(record["AssaySourceNameList"], [])
        self.assertEqual(record["MinAC"], "")
        self.assertEqual(record["MaxAC"], "")
        self.assertEqual(record["MinTC"], "")
        self.assertEqual(record["MaxTC"], "")
        self.assertEqual(record["ActiveAidCount"], 1)
        self.assertEqual(record["InactiveAidCount"], None)
        self.assertEqual(record["TotalAidCount"], 243)
        self.assertEqual(record["InChIKey"], "LXEJRKJRKIFVNY-UHFFFAOYSA-N")
        self.assertEqual(
            record["InChI"], "InChI=1S/C8H4Cl2O2/c9-7(11)5-1-2-6(4-3-5)8(10)12/h1-4H"
        )


class ELinkTest(unittest.TestCase):
    """Tests for parsing XML output returned by ELink."""

    def test_pubmed1(self):
        """Test parsing pubmed links returned by ELink (first test)."""
        # Retrieve IDs from PubMed for PMID 9298984 to the PubMed database
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="9298984", cmd="neighbor")
        with open("Entrez/elink1.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 1)
        self.assertEqual(len(record[0]), 5)
        self.assertEqual(record[0]["DbFrom"], "pubmed")
        self.assertEqual(record[0]["IdList"], ["9298984"])
        self.assertEqual(len(record[0]["LinkSetDb"]), 7)
        self.assertEqual(record[0]["LinkSetDb"][0]["DbTo"], "pubmed")
        self.assertEqual(record[0]["LinkSetDb"][0]["LinkName"], "pubmed_pubmed")
        self.assertEqual(len(record[0]["LinkSetDb"][0]["Link"]), 101)
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][0]["Id"], "9298984")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][1]["Id"], "8794856")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][2]["Id"], "9700164")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][3]["Id"], "7914521")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][4]["Id"], "9914369")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][5]["Id"], "1339459")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][6]["Id"], "11590237")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][7]["Id"], "2211822")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][8]["Id"], "12686595")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][9]["Id"], "20980244")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][10]["Id"], "11146659")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][11]["Id"], "8978614")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][12]["Id"], "10893249")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][13]["Id"], "10402457")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][14]["Id"], "15371539")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][15]["Id"], "9074495")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][16]["Id"], "10806105")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][17]["Id"], "9490715")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][18]["Id"], "15915585")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][19]["Id"], "10545493")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][20]["Id"], "10523511")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][21]["Id"], "11483958")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][22]["Id"], "9869638")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][23]["Id"], "7690762")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][24]["Id"], "9425896")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][25]["Id"], "26892014")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][26]["Id"], "9378750")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][27]["Id"], "12515822")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][28]["Id"], "25919583")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][29]["Id"], "38830800")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][30]["Id"], "25081981")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][31]["Id"], "1691829")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][32]["Id"], "11146661")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][33]["Id"], "11685532")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][34]["Id"], "12080088")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][35]["Id"], "12034769")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][36]["Id"], "9852156")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][37]["Id"], "22733107")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][38]["Id"], "25194162")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][39]["Id"], "8923204")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][40]["Id"], "2022189")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][41]["Id"], "10985388")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][42]["Id"], "38402459")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][43]["Id"], "35609608")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][44]["Id"], "17222555")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][45]["Id"], "16741559")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][46]["Id"], "18936247")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][47]["Id"], "10749938")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][48]["Id"], "31150390")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][49]["Id"], "32794572")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][50]["Id"], "17895365")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][51]["Id"], "23300382")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][52]["Id"], "11914278")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][53]["Id"], "22563370")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][54]["Id"], "1541637")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][55]["Id"], "29706521")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][56]["Id"], "21118145")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][57]["Id"], "16732327")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][58]["Id"], "12388768")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][59]["Id"], "18202360")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][60]["Id"], "7585942")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][61]["Id"], "11179694")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][62]["Id"], "29158400")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][63]["Id"], "11352945")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][64]["Id"], "8056842")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][65]["Id"], "29046402")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][66]["Id"], "10398680")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][67]["Id"], "11267866")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][68]["Id"], "16516834")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][69]["Id"], "16839185")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][70]["Id"], "15616189")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][71]["Id"], "11266459")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][72]["Id"], "19641019")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][73]["Id"], "25976696")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][74]["Id"], "31250100")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][75]["Id"], "30784092")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][76]["Id"], "38019881")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][77]["Id"], "17182852")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][78]["Id"], "2211824")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][79]["Id"], "14522947")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][80]["Id"], "15268859")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][81]["Id"], "11252055")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][82]["Id"], "8175879")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][83]["Id"], "11102811")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][84]["Id"], "7904902")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][85]["Id"], "9606208")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][86]["Id"], "18460473")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][87]["Id"], "36920098")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][88]["Id"], "16510521")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][89]["Id"], "11092768")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][90]["Id"], "15824131")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][91]["Id"], "12235289")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][92]["Id"], "11266451")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][93]["Id"], "15485811")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][94]["Id"], "10898791")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][95]["Id"], "20729837")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][96]["Id"], "8548823")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][97]["Id"], "6807996")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][98]["Id"], "6791901")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][99]["Id"], "11715021")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][100]["Id"], "17333235")
        self.assertEqual(len(record[0]["LinkSetDb"][1]["Link"]), 39)
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][0]["Id"], "38830800")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][1]["Id"], "38188366")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][2]["Id"], "37424454")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][3]["Id"], "34205694")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][4]["Id"], "32052088")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][5]["Id"], "29915359")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][6]["Id"], "29475948")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][7]["Id"], "29423089")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][8]["Id"], "29192061")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][9]["Id"], "28320824")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][10]["Id"], "28125061")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][11]["Id"], "20439434")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][12]["Id"], "19273145")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][13]["Id"], "19177000")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][14]["Id"], "18936247")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][15]["Id"], "18268100")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][16]["Id"], "17699596")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][17]["Id"], "16563186")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][18]["Id"], "16505164")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][19]["Id"], "16107559")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][20]["Id"], "15824131")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][21]["Id"], "15289669")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][22]["Id"], "15029241")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][23]["Id"], "12906131")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][24]["Id"], "12686595")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][25]["Id"], "12498345")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][26]["Id"], "11756470")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][27]["Id"], "11553716")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][28]["Id"], "11500386")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][29]["Id"], "11402076")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][30]["Id"], "11331754")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][31]["Id"], "10780705")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][32]["Id"], "10545493")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][33]["Id"], "10402457")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][34]["Id"], "10402425")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][35]["Id"], "9914368")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][36]["Id"], "9763420")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][37]["Id"], "9700166")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][38]["Id"], "9700164")
        self.assertEqual(len(record[0]["LinkSetDb"][2]["Link"]), 5)
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][0]["Id"], "9298984")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][1]["Id"], "8794856")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][2]["Id"], "9700164")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][3]["Id"], "7914521")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][4]["Id"], "38830800")
        self.assertEqual(len(record[0]["LinkSetDb"][3]["Link"]), 5)
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][0]["Id"], "9298984")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][1]["Id"], "8794856")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][2]["Id"], "9700164")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][3]["Id"], "7914521")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][4]["Id"], "9914369")
        self.assertEqual(len(record[0]["LinkSetDb"][4]["Link"]), 56)
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][0]["Id"], "14732139")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][1]["Id"], "8909532")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][2]["Id"], "8898221")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][3]["Id"], "8824189")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][4]["Id"], "8824188")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][5]["Id"], "8794856")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][6]["Id"], "8763498")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][7]["Id"], "8706132")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][8]["Id"], "8706131")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][9]["Id"], "8647893")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][10]["Id"], "8617505")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][11]["Id"], "8560259")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][12]["Id"], "8521491")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][13]["Id"], "8505381")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][14]["Id"], "8485583")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][15]["Id"], "8416984")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][16]["Id"], "8267981")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][17]["Id"], "8143084")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][18]["Id"], "8023161")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][19]["Id"], "8005447")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][20]["Id"], "7914521")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][21]["Id"], "7906398")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][22]["Id"], "7860624")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][23]["Id"], "7854443")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][24]["Id"], "7854422")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][25]["Id"], "7846151")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][26]["Id"], "7821090")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][27]["Id"], "7758115")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][28]["Id"], "7739381")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][29]["Id"], "7704412")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][30]["Id"], "7698647")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][31]["Id"], "7664339")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][32]["Id"], "7642709")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][33]["Id"], "7642708")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][34]["Id"], "7579695")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][35]["Id"], "7542657")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][36]["Id"], "7502067")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][37]["Id"], "7172865")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][38]["Id"], "6966403")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][39]["Id"], "6793236")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][40]["Id"], "6684600")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][41]["Id"], "3928429")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][42]["Id"], "3670292")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][43]["Id"], "2686123")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][44]["Id"], "2683077")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][45]["Id"], "2512302")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][46]["Id"], "2498337")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][47]["Id"], "2195725")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][48]["Id"], "2185478")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][49]["Id"], "2139718")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][50]["Id"], "2139717")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][51]["Id"], "2022189")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][52]["Id"], "1999466")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][53]["Id"], "1684022")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][54]["Id"], "1406971")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][55]["Id"], "1339459")
        self.assertEqual(len(record[0]["LinkSetDb"][5]["Link"]), 1)
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][0]["Id"], "9298984")
        self.assertEqual(len(record[0]["LinkSetDb"][6]["Link"]), 5)
        self.assertEqual(record[0]["LinkSetDb"][6]["Link"][0]["Id"], "9298984")
        self.assertEqual(record[0]["LinkSetDb"][6]["Link"][1]["Id"], "8794856")
        self.assertEqual(record[0]["LinkSetDb"][6]["Link"][2]["Id"], "9700164")
        self.assertEqual(record[0]["LinkSetDb"][6]["Link"][3]["Id"], "7914521")
        self.assertEqual(record[0]["LinkSetDb"][6]["Link"][4]["Id"], "9914369")

    def test_nucleotide(self):
        """Test parsing Nucleotide to Protein links returned by ELink."""
        # Retrieve IDs from Nucleotide for GI  48819, 7140345 to Protein
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="nucleotide", db="protein",
        #                      id="48819,7140345")
        with open("Entrez/elink2.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 1)
        self.assertEqual(len(record[0]), 5)
        self.assertEqual(record[0]["DbFrom"], "nuccore")
        self.assertEqual(record[0]["IdList"], ["48819", "7140345"])
        self.assertEqual(len(record[0]["LinkSetDb"]), 1)
        self.assertEqual(len(record[0]["LinkSetDb"][0]), 3)
        self.assertEqual(record[0]["LinkSetDb"][0]["DbTo"], "protein")
        self.assertEqual(record[0]["LinkSetDb"][0]["LinkName"], "nuccore_protein")
        self.assertEqual(len(record[0]["LinkSetDb"][0]["Link"]), 1)
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][0]["Id"], "48820")

    def test_pubmed2(self):
        """Test parsing pubmed links returned by ELink (second test)."""
        # Retrieve PubMed related articles for PMIDs 11812492 11774222
        # with a publication date from 1995 to the present
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="11812492,11774222",
        #                      db="pubmed", mindate="1995", datetype="pdat")
        with open("Entrez/elink3.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 1)
        self.assertEqual(record[0]["DbFrom"], "pubmed")
        self.assertEqual(len(record[0]["IdList"]), 2)
        self.assertEqual(record[0]["IdList"][0], "11812492")
        self.assertEqual(record[0]["IdList"][1], "11774222")
        self.assertEqual(record[0]["LinkSetDb"][0]["DbTo"], "pubmed")
        self.assertEqual(record[0]["LinkSetDb"][0]["LinkName"], "pubmed_pubmed")
        self.assertEqual(len(record[0]["LinkSetDb"]), 6)
        self.assertEqual(len(record[0]["LinkSetDb"][0]["Link"]), 284)
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][0]["Id"], "39386366")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][1]["Id"], "39338906")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][2]["Id"], "39106844")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][281]["Id"], "10092480")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][282]["Id"], "9830540")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][283]["Id"], "9274032")
        self.assertEqual(len(record[0]["LinkSetDb"][1]["Link"]), 3)
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][0]["Id"], "18508935")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][1]["Id"], "15780005")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][2]["Id"], "15024419")
        self.assertEqual(len(record[0]["LinkSetDb"][2]["Link"]), 10)
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][0]["Id"], "16005284")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][1]["Id"], "15780005")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][2]["Id"], "15111095")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][7]["Id"], "11668631")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][8]["Id"], "10731564")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][9]["Id"], "10612825")
        self.assertEqual(len(record[0]["LinkSetDb"][3]["Link"]), 10)
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][0]["Id"], "28358880")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][1]["Id"], "24053607")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][2]["Id"], "15780005")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][7]["Id"], "11668631")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][8]["Id"], "10731564")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][9]["Id"], "10612825")
        self.assertEqual(len(record[0]["LinkSetDb"][4]["Link"]), 282)
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][0]["Id"], "40503031")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][1]["Id"], "39530225")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][2]["Id"], "39386366")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][279]["Id"], "10092480")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][280]["Id"], "9830540")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][281]["Id"], "9274032")
        self.assertEqual(len(record[0]["LinkSetDb"][5]["Link"]), 10)
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][0]["Id"], "28358880")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][1]["Id"], "24053607")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][2]["Id"], "15780005")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][7]["Id"], "11668631")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][8]["Id"], "10731564")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][9]["Id"], "10612825")

    def test_medline(self):
        """Test parsing medline indexed articles returned by ELink."""
        # Retrieve MEDLINE indexed only related articles for PMID 12242737
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12242737", db="pubmed",
        #                      term="medline[sb]")
        with open("Entrez/elink4.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 1)
        self.assertEqual(record[0]["DbFrom"], "pubmed")
        self.assertEqual(record[0]["IdList"], ["12242737"])
        self.assertEqual(len(record[0]["LinkSetDb"]), 6)
        self.assertEqual(record[0]["LinkSetDb"][0]["DbTo"], "pubmed")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][0]["Id"], "38997184")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][1]["Id"], "37474462")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][2]["Id"], "36642882")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][3]["Id"], "33097552")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][4]["Id"], "32874413")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][5]["Id"], "32345945")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][6]["Id"], "31558954")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][7]["Id"], "31415410")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][8]["Id"], "31011388")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][9]["Id"], "30255446")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][10]["Id"], "30255443")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][11]["Id"], "30115443")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][12]["Id"], "29340557")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][13]["Id"], "29331013")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][14]["Id"], "29233545")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][15]["Id"], "29115651")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][16]["Id"], "29022115")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][17]["Id"], "27597304")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][18]["Id"], "27315096")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][19]["Id"], "27304929")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][20]["Id"], "27142382")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][21]["Id"], "26965844")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][22]["Id"], "26481976")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][23]["Id"], "26427946")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][24]["Id"], "26331169")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][25]["Id"], "25690945")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][26]["Id"], "25669229")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][27]["Id"], "25646204")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][28]["Id"], "25584961")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][29]["Id"], "25555004")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][30]["Id"], "25470877")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][31]["Id"], "25223134")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][32]["Id"], "25167349")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][33]["Id"], "24934824")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][34]["Id"], "24879722")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][35]["Id"], "24836494")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][36]["Id"], "24309417")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][37]["Id"], "23978699")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][38]["Id"], "23759294")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][39]["Id"], "23570763")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][40]["Id"], "23255877")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][41]["Id"], "22688104")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][42]["Id"], "22661362")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][43]["Id"], "22648258")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][44]["Id"], "22521021")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][45]["Id"], "22424988")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][46]["Id"], "22369817")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][47]["Id"], "22368911")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][48]["Id"], "22194507")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][49]["Id"], "22156652")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][50]["Id"], "22109321")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][51]["Id"], "21984464")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][52]["Id"], "21944608")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][53]["Id"], "21908142")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][54]["Id"], "21715237")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][55]["Id"], "21153952")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][56]["Id"], "20860230")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][57]["Id"], "20718377")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][58]["Id"], "20674629")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][59]["Id"], "20558858")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][60]["Id"], "20533237")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][61]["Id"], "20016426")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][62]["Id"], "19843737")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][63]["Id"], "19616724")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][64]["Id"], "19520357")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][65]["Id"], "18783095")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][66]["Id"], "18582671")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][67]["Id"], "18554854")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][68]["Id"], "18053822")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][69]["Id"], "18021675")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][70]["Id"], "17875143")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][71]["Id"], "17875142")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][72]["Id"], "17879696")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][73]["Id"], "17602359")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][74]["Id"], "17601500")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][75]["Id"], "17376366")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][76]["Id"], "17354190")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][77]["Id"], "17325998")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][78]["Id"], "17243036")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][79]["Id"], "17205643")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][80]["Id"], "17193860")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][81]["Id"], "17174054")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][82]["Id"], "17040637")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][83]["Id"], "16999328")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][84]["Id"], "16988291")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][85]["Id"], "16580806")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][86]["Id"], "16566645")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][87]["Id"], "16552382")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][88]["Id"], "16357381")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][89]["Id"], "16284132")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][90]["Id"], "16133609")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][91]["Id"], "16096604")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][92]["Id"], "16046437")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][93]["Id"], "15835031")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][94]["Id"], "15788585")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][95]["Id"], "15788584")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][96]["Id"], "15505294")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][97]["Id"], "15278705")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][98]["Id"], "15236131")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][99]["Id"], "15143223")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][100]["Id"], "15141648")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][101]["Id"], "15136027")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][102]["Id"], "15094630")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][103]["Id"], "15022983")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][104]["Id"], "14661668")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][105]["Id"], "14661661")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][106]["Id"], "14650118")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][107]["Id"], "12878072")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][108]["Id"], "12846253")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][109]["Id"], "12822521")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][110]["Id"], "12733684")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][111]["Id"], "12719915")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][112]["Id"], "12563154")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][113]["Id"], "12242737")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][114]["Id"], "12226761")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][115]["Id"], "12164574")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][116]["Id"], "12069469")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][117]["Id"], "11973040")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][118]["Id"], "11895298")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][119]["Id"], "11781922")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][120]["Id"], "11775722")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][121]["Id"], "11762248")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][122]["Id"], "11702119")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][123]["Id"], "11368937")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][124]["Id"], "11329656")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][125]["Id"], "11329655")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][126]["Id"], "11329162")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][127]["Id"], "11274884")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][128]["Id"], "11218011")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][129]["Id"], "11125632")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][130]["Id"], "11016058")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][131]["Id"], "10688063")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][132]["Id"], "10499696")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][133]["Id"], "10222515")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][134]["Id"], "10222514")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][135]["Id"], "10024396")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][136]["Id"], "9793138")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][137]["Id"], "9757294")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][138]["Id"], "9575723")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][139]["Id"], "9510579")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][140]["Id"], "9456947")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][141]["Id"], "9314960")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][142]["Id"], "9314959")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][143]["Id"], "9269670")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][144]["Id"], "9193407")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][145]["Id"], "8872409")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][146]["Id"], "8756148")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][147]["Id"], "8903064")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][148]["Id"], "8599783")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][149]["Id"], "8153333")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][150]["Id"], "1343378")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][151]["Id"], "1535863")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][152]["Id"], "4818442")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][153]["Id"], "4808999")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][154]["Id"], "13969511")
        self.assertEqual(record[0]["LinkSetDb"][0]["Link"][155]["Id"], "13808134")
        self.assertEqual(record[0]["LinkSetDb"][0]["LinkName"], "pubmed_pubmed")
        self.assertEqual(record[0]["LinkSetDb"][1]["DbTo"], "pubmed")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][0]["Id"], "40900773")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][1]["Id"], "28779191")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][2]["Id"], "21102533")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][3]["Id"], "19517148")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][4]["Id"], "19132488")
        self.assertEqual(record[0]["LinkSetDb"][1]["Link"][5]["Id"], "15278705")
        self.assertEqual(record[0]["LinkSetDb"][1]["LinkName"], "pubmed_pubmed_citedin")
        self.assertEqual(record[0]["LinkSetDb"][2]["DbTo"], "pubmed")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][0]["Id"], "20718377")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][1]["Id"], "12242737")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][2]["Id"], "11329656")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][3]["Id"], "11218011")
        self.assertEqual(record[0]["LinkSetDb"][2]["Link"][4]["Id"], "9757294")
        self.assertEqual(
            record[0]["LinkSetDb"][2]["LinkName"], "pubmed_pubmed_combined"
        )
        self.assertEqual(record[0]["LinkSetDb"][3]["DbTo"], "pubmed")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][0]["Id"], "20718377")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][1]["Id"], "17193860")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][2]["Id"], "12242737")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][3]["Id"], "11218011")
        self.assertEqual(record[0]["LinkSetDb"][3]["Link"][4]["Id"], "9757294")
        self.assertEqual(record[0]["LinkSetDb"][3]["LinkName"], "pubmed_pubmed_five")
        self.assertEqual(record[0]["LinkSetDb"][4]["DbTo"], "pubmed")
        self.assertEqual(record[0]["LinkSetDb"][4]["Link"][0]["Id"], "12242737")
        self.assertEqual(record[0]["LinkSetDb"][4]["LinkName"], "pubmed_pubmed_reviews")
        self.assertEqual(record[0]["LinkSetDb"][5]["DbTo"], "pubmed")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][0]["Id"], "20718377")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][1]["Id"], "17193860")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][2]["Id"], "12242737")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][3]["Id"], "11218011")
        self.assertEqual(record[0]["LinkSetDb"][5]["Link"][4]["Id"], "9757294")
        self.assertEqual(
            record[0]["LinkSetDb"][5]["LinkName"], "pubmed_pubmed_reviews_five"
        )

    def test_pubmed3(self):
        """Test parsing pubmed link returned by ELink (third test)."""
        # Create a hyperlink to the first link available for PMID 10611131
        # in PubMed
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="10611131", cmd="prlinks")

        with open("Entrez/elink5.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 1)
        self.assertEqual(len(record[0]), 5)
        self.assertEqual(record[0]["DbFrom"], "pubmed")
        self.assertEqual(len(record[0]["LinkSetDb"]), 0)
        self.assertEqual(len(record[0]["LinkSetDbHistory"]), 0)
        self.assertEqual(len(record[0]["ERROR"]), 0)
        self.assertEqual(len(record[0]["IdUrlList"]), 2)
        self.assertEqual(len(record[0]["IdUrlList"]["FirstChars"]), 0)
        self.assertEqual(len(record[0]["IdUrlList"]["IdUrlSet"]), 1)

        self.assertEqual(record[0]["IdUrlList"]["IdUrlSet"][0]["Id"], "10611131")
        self.assertEqual(len(record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"]), 1)
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Url"],
            "https://academic.oup.com/brain/article-lookup/doi/10.1093/brain/123.1.171",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Url"].attributes,
            {"LNG": "EN"},
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["IconUrl"],
            "/corehtml/query/egifs/https:--academic.oup.com-images-oup_pubmed.png",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["IconUrl"].attributes,
            {"LNG": "EN"},
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["SubjectType"],
            [],
        )
        self.assertEqual(
            len(record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Attribute"]), 1
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Attribute"][0],
            "subscription/membership/fee required",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["Name"],
            "Silverchair Information Systems",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["NameAbbr"],
            "silverchair",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["Id"], "7898"
        )

    def test_pubmed4(self):
        """Test parsing pubmed links returned by ELink (fourth test)."""
        # List all available links in PubMed, except for libraries, for
        # PMIDs 12085856 and 12085853
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12085856,12085853", cmd="llinks")
        with open("Entrez/elink6.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record[0]["DbFrom"], "pubmed")
        self.assertEqual(len(record[0]["IdUrlList"]), 2)
        self.assertEqual(record[0]["IdUrlList"]["IdUrlSet"][0]["Id"], "12085856")
        self.assertEqual(len(record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"]), 1)
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Category"], ["Medical"]
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Url"],
            "https://medlineplus.gov/coronaryarterybypasssurgery.html",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Attribute"],
            ["free resource"],
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["SubjectType"],
            [],
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["IconUrl"],
            "/corehtml/query/egifs/https:--medlineplus.gov-images-linkout_sm.gif",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["Name"],
            "MedlinePlus Health Information",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["NameAbbr"],
            "MEDPLUS",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["Provider"]["Id"], "3162"
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][0]["ObjUrl"][0]["LinkName"],
            "Coronary Artery Bypass Surgery",
        )
        self.assertEqual(len(record[0]["IdUrlList"]["IdUrlSet"][1]), 2)
        self.assertEqual(record[0]["IdUrlList"]["IdUrlSet"][1]["Id"], "12085853")
        self.assertEqual(len(record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"]), 2)
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Category"], ["Medical"]
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Url"],
            "https://medlineplus.gov/exerciseandphysicalfitness.html",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["IconUrl"],
            "/corehtml/query/egifs/https:--medlineplus.gov-images-linkout_sm.gif",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Attribute"],
            ["free resource"],
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["SubjectType"],
            [],
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["LinkName"],
            "Exercise and Physical Fitness",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Provider"]["Name"],
            "MedlinePlus Health Information",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Provider"]["NameAbbr"],
            "MEDPLUS",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][0]["Provider"]["Id"], "3162"
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Category"], ["Medical"]
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Attribute"],
            ["free resource"],
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Url"],
            "https://medlineplus.gov/exerciseandphysicalfitness.html",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["IconUrl"],
            "/corehtml/query/egifs/https:--medlineplus.gov-images-linkout_sm.gif",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["LinkName"],
            "Exercise and Physical Fitness",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["SubjectType"],
            [],
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Provider"]["Name"],
            "MedlinePlus Consumer Health Information",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Provider"]["NameAbbr"],
            "medlineplus2",
        )
        self.assertEqual(
            record[0]["IdUrlList"]["IdUrlSet"][1]["ObjUrl"][1]["Provider"]["Id"],
            "10405",
        )

    def test_pubmed5(self):
        """Test parsing pubmed links returned by ELink (fifth test)."""
        # List Entrez database links for PubMed PMIDs 12169658 and 11748140
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12169658,11748140",
        #                      cmd="acheck")
        with open("Entrez/elink7.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 1)
        self.assertEqual(record[0]["DbFrom"], "pubmed")
        self.assertEqual(len(record[0]["IdCheckList"]), 2)
        self.assertEqual(record[0]["IdCheckList"]["IdLinkSet"][0]["Id"], "12169658")
        self.assertEqual(len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"]), 16)
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["DbTo"], "books"
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["LinkName"],
            "pubmed_books_refs",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["MenuTag"],
            "Cited in Books",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["HtmlTag"],
            "Cited in Books",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["Priority"], "185"
        )
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["DbTo"], "cdd"
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["LinkName"],
            "pubmed_cdd",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["MenuTag"],
            "Conserved Domain Links",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["HtmlTag"],
            "Conserved Domains",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["Priority"], "130"
        )
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["DbTo"], "gene"
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["LinkName"],
            "pubmed_gene",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["MenuTag"],
            "Gene Links",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["HtmlTag"],
            "Gene",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["Priority"], "128"
        )
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["DbTo"],
            "geoprofiles",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["LinkName"],
            "pubmed_geoprofiles",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["MenuTag"],
            "GEO Profile Links",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["HtmlTag"],
            "GEO Profiles",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["Priority"], "170"
        )
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["DbTo"], "nuccore"
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["LinkName"],
            "pubmed_nuccore",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["MenuTag"],
            "Nucleotide Links",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["HtmlTag"],
            "Nucleotide",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["Priority"], "128"
        )
        self.assertEqual(len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"]), 15)
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["DbTo"], "books"
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["LinkName"],
            "pubmed_books_refs",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["MenuTag"],
            "Cited in Books",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["HtmlTag"],
            "Cited in Books",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["Priority"], "185"
        )
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["DbTo"], "gene"
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["LinkName"],
            "pubmed_gene",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["MenuTag"],
            "Gene Links",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["HtmlTag"],
            "Gene",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["Priority"], "128"
        )
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["DbTo"],
            "geoprofiles",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["LinkName"],
            "pubmed_geoprofiles",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["MenuTag"],
            "GEO Profile Links",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["HtmlTag"],
            "GEO Profiles",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["Priority"], "170"
        )
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["DbTo"], "nuccore"
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["LinkName"],
            "pubmed_nuccore",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["MenuTag"],
            "Nucleotide Links",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["HtmlTag"],
            "Nucleotide",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["Priority"], "128"
        )
        self.assertEqual(
            len(record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]), 5
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["DbTo"], "nuccore"
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["LinkName"],
            "pubmed_nuccore_refseq",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["MenuTag"],
            "Nucleotide (RefSeq) Links",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["HtmlTag"],
            "Nucleotide (RefSeq)",
        )
        self.assertEqual(
            record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["Priority"], "128"
        )

    def test_pubmed6(self):
        """Test parsing pubmed links returned by ELink (sixth test)."""
        # Check for the existence of a Related Articles link for PMID
        # 12068369.
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12068369", cmd="ncheck")

        with open("Entrez/elink8.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 1)
        self.assertEqual(record[0]["DbFrom"], "pubmed")
        self.assertEqual(len(record[0]["IdCheckList"]), 2)
        self.assertEqual(len(record[0]["IdCheckList"]["Id"]), 1)
        self.assertEqual(record[0]["IdCheckList"]["Id"][0], "12068369")
        self.assertEqual(len(record[0]["IdCheckList"]["Id"][0].attributes), 1)
        self.assertEqual(
            record[0]["IdCheckList"]["Id"][0].attributes["HasNeighbor"], "Y"
        )
        self.assertEqual(len(record[0]["IdCheckList"]["IdLinkSet"]), 0)


class ESpellTest(unittest.TestCase):
    """Tests for parsing XML output returned by ESpell."""

    def test_espell(self):
        """Test parsing XML output returned by ESpell."""
        # Request suggestions for the PubMed search biopythooon
        # To create the XML file, use
        # >>> Bio.Entrez.espell(db="pubmed", term="biopythooon")
        with open("Entrez/espell.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record["Database"], "pubmed")
        self.assertEqual(record["Query"], "biopythooon")
        self.assertEqual(record["CorrectedQuery"], "biopython")
        self.assertEqual(len(record["SpelledQuery"]), 2)
        self.assertEqual(record["SpelledQuery"][0], "")
        self.assertEqual(record["SpelledQuery"][0].tag, "Original")
        self.assertEqual(record["SpelledQuery"][1], "biopython")
        self.assertEqual(record["SpelledQuery"][1].tag, "Replaced")


class EFetchTest(unittest.TestCase):
    """Tests for parsing XML output returned by EFetch."""

    def test_pubmed1(self):
        """Test parsing XML returned by EFetch, PubMed database (first test)."""
        # In PubMed display PMIDs 12091962 and 9997 in xml retrieval mode
        # and abstract retrieval type.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', id='12091962,9997',
        #                       retmode='xml', rettype='abstract')
        with open("Entrez/pubmed1.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        self.assertEqual(record["PubmedBookArticle"], [])
        record = record["PubmedArticle"]
        self.assertEqual(record[0]["MedlineCitation"].attributes["Owner"], "KIE")
        self.assertEqual(record[0]["MedlineCitation"].attributes["Status"], "MEDLINE")
        self.assertEqual(record[0]["MedlineCitation"]["PMID"], "12091962")
        self.assertEqual(record[0]["MedlineCitation"]["DateCompleted"]["Year"], "1991")
        self.assertEqual(record[0]["MedlineCitation"]["DateCompleted"]["Month"], "01")
        self.assertEqual(record[0]["MedlineCitation"]["DateCompleted"]["Day"], "22")
        self.assertEqual(record[0]["MedlineCitation"]["DateRevised"]["Year"], "2007")
        self.assertEqual(record[0]["MedlineCitation"]["DateRevised"]["Month"], "11")
        self.assertEqual(record[0]["MedlineCitation"]["DateRevised"]["Day"], "15")
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"].attributes["PubModel"], "Print"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"], "1043-1578"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes[
                "IssnType"
            ],
            "Print",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"],
            "Print",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Volume"
            ],
            "17",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"],
            "1",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Year"],
            "1990",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Season"],
            "Spring",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["Title"],
            "Social justice (San Francisco, Calif.)",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["ArticleTitle"],
            "The treatment of AIDS behind the walls of correctional facilities.",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"],
            "113-25",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"].attributes[
                "CompleteYN"
            ],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"],
            "Olivero",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"],
            "J Michael",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"], "JM"
        )
        self.assertEqual(record[0]["MedlineCitation"]["Article"]["Language"], ["eng"])
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["PublicationTypeList"],
            ["Journal Article", "Review"],
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MedlineJournalInfo"]["Country"],
            "United States",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"],
            "Soc Justice",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"], "9891830"
        )
        self.assertEqual(record[0]["MedlineCitation"]["CitationSubset"], [])
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"],
            "AIDS Serodiagnosis",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][0][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"],
            "Acquired Immunodeficiency Syndrome",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][1][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"],
            "Civil Rights",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][2][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"],
            "HIV Seropositivity",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][3][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"],
            "Humans",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][4][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"],
            "Jurisprudence",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][5][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"],
            "Law Enforcement",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][6][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"],
            "Mass Screening",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][7][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"],
            "Minority Groups",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][8][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"],
            "Organizational Policy",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][9][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"],
            "Patient Care",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][10][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"],
            "Prejudice",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][11][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"],
            "Prisoners",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][12][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][13]["DescriptorName"],
            "Public Policy",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][13][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][14]["DescriptorName"],
            "Quarantine",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][14][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][15]["DescriptorName"],
            "Social Control, Formal",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][15][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][16]["DescriptorName"],
            "Statistics as Topic",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][16][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][17]["DescriptorName"],
            "Stereotyping",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][17][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][18]["DescriptorName"],
            "United States",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][18][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(record[0]["MedlineCitation"]["NumberOfReferences"], "63")
        self.assertEqual(record[0]["MedlineCitation"]["OtherID"][0], "31840")
        self.assertEqual(
            record[0]["MedlineCitation"]["OtherID"][0].attributes["Source"], "KIE"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["KeywordList"][0].attributes["Owner"], "KIE"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["KeywordList"][0][0],
            "Health Care and Public Health",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["KeywordList"][0][0].attributes[
                "MajorTopicYN"
            ],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["KeywordList"][0][1], "Legal Approach"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["KeywordList"][0][1].attributes[
                "MajorTopicYN"
            ],
            "N",
        )
        self.assertEqual(record[0]["MedlineCitation"]["GeneralNote"][0], "14 fn.")
        self.assertEqual(
            record[0]["MedlineCitation"]["GeneralNote"][0].attributes["Owner"], "KIE"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["GeneralNote"][1],
            "KIE BoB Subject Heading: AIDS",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["GeneralNote"][1].attributes["Owner"], "KIE"
        )
        self.assertEqual(record[0]["MedlineCitation"]["GeneralNote"][2], "63 refs.")
        self.assertEqual(
            record[0]["MedlineCitation"]["GeneralNote"][2].attributes["Owner"], "KIE"
        )
        self.assertEqual(
            record[0]["PubmedData"]["History"][0].attributes["PubStatus"], "pubmed"
        )
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Year"], "1990")
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Month"], "4")
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Day"], "1")
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Hour"], "0")
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Minute"], "0")
        self.assertEqual(
            record[0]["PubmedData"]["History"][1].attributes["PubStatus"], "medline"
        )
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Year"], "2002")
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Month"], "7")
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Day"], "16")
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Hour"], "10")
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Minute"], "1")
        self.assertEqual(record[0]["PubmedData"]["PublicationStatus"], "ppublish")
        self.assertEqual(len(record[0]["PubmedData"]["ArticleIdList"]), 1)
        self.assertEqual(record[0]["PubmedData"]["ArticleIdList"][0], "12091962")
        self.assertEqual(
            record[0]["PubmedData"]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(record[1]["MedlineCitation"].attributes["Owner"], "NLM")
        self.assertEqual(record[1]["MedlineCitation"].attributes["Status"], "MEDLINE")
        self.assertEqual(record[1]["MedlineCitation"]["PMID"], "9997")
        self.assertEqual(record[1]["MedlineCitation"]["DateCompleted"]["Year"], "1976")
        self.assertEqual(record[1]["MedlineCitation"]["DateCompleted"]["Month"], "12")
        self.assertEqual(record[1]["MedlineCitation"]["DateCompleted"]["Day"], "30")
        self.assertEqual(record[1]["MedlineCitation"]["DateRevised"]["Year"], "2019")
        self.assertEqual(record[1]["MedlineCitation"]["DateRevised"]["Month"], "06")
        self.assertEqual(record[1]["MedlineCitation"]["DateRevised"]["Day"], "09")
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"].attributes["PubModel"], "Print"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"], "0006-3002"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes["IssnType"],
            "Print",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].attributes["CitedMedium"],
            "Print",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"],
            "446",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"],
            "1",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"],
            "1976",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Month"],
            "Sep",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Day"],
            "28",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["Title"],
            "Biochimica et biophysica acta",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"],
            "Biochim Biophys Acta",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["ArticleTitle"],
            "Magnetic studies of Chromatium flavocytochrome C552. A mechanism for heme-flavin interaction.",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"],
            "179-91",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"],
            ["Electron paramagnetic resonance and magnetic susceptibility studies of Chromatium flavocytochrome C552 and its diheme flavin-free subunit at temperatures below 45 degrees K are reported. The results show that in the intact protein and the subunit the two low-spin (S = 1/2) heme irons are distinguishable, giving rise to separate EPR signals. In the intact protein only, one of the heme irons exists in two different low spin environments in the pH range 5.5 to 10.5, while the other remains in a constant environment. Factors influencing the variable heme iron environment also influence flavin reactivity, indicating the existence of a mechanism for heme-flavin interaction."],
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"],
            "Strekas",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"], "T C"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"], "TC"
        )
        self.assertEqual(record[1]["MedlineCitation"]["Article"]["Language"], ["eng"])
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["PublicationTypeList"],
            ["Journal Article"],
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MedlineJournalInfo"]["Country"], "Netherlands"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"],
            "Biochim Biophys Acta",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"], "0217513"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["ChemicalList"][0]["RegistryNumber"], "0"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["ChemicalList"][0]["NameOfSubstance"],
            "Cytochrome c Group",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["ChemicalList"][1]["RegistryNumber"], "0"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["ChemicalList"][1]["NameOfSubstance"],
            "Flavins",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["ChemicalList"][2]["RegistryNumber"],
            "42VZT0U6YR",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["ChemicalList"][2]["NameOfSubstance"], "Heme"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["ChemicalList"][3]["RegistryNumber"],
            "E1UOL152H7",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["ChemicalList"][3]["NameOfSubstance"], "Iron"
        )
        self.assertEqual(record[1]["MedlineCitation"]["CitationSubset"], ["IM"])
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"],
            "Binding Sites",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"],
            "Chromatium",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0],
            "enzymology",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"],
            "Cytochrome c Group",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"],
            "Electron Spin Resonance Spectroscopy",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"],
            "Flavins",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"], "Heme"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"],
            "Hydrogen-Ion Concentration",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"], "Iron"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0],
            "analysis",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"],
            "Magnetics",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"],
            "Oxidation-Reduction",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"],
            "Protein Binding",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"],
            "Protein Conformation",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"],
            "Temperature",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(record[1]["PubmedData"]["History"][0].attributes["PubStatus"], "pubmed")
        self.assertEqual(record[1]["PubmedData"]["History"][0]["Year"], "1976")
        self.assertEqual(record[1]["PubmedData"]["History"][0]["Month"], "9")
        self.assertEqual(record[1]["PubmedData"]["History"][0]["Day"], "28")
        self.assertEqual(record[1]["PubmedData"]["History"][1].attributes["PubStatus"], "medline")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Year"], "1976")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Month"], "9")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Day"], "28")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Hour"], "0")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Minute"], "1")
        self.assertEqual(record[1]["PubmedData"]["PublicationStatus"], "ppublish")
        self.assertEqual(len(record[1]["PubmedData"]["ArticleIdList"]), 3)
        self.assertEqual(record[1]["PubmedData"]["ArticleIdList"][0], "9997")
        self.assertEqual(record[1]["PubmedData"]["ArticleIdList"][0].attributes["IdType"], "pubmed")
        self.assertEqual(record[1]["PubmedData"]["ArticleIdList"][1], "10.1016/0005-2795(76)90109-4")
        self.assertEqual(record[1]["PubmedData"]["ArticleIdList"][1].attributes["IdType"], "doi")
        self.assertEqual(record[1]["PubmedData"]["ArticleIdList"][2], "0005-2795(76)90109-4")
        self.assertEqual(record[1]["PubmedData"]["ArticleIdList"][2].attributes["IdType"], "pii")
        # fmt: on

    def test_pubmed2(self):
        """Test parsing XML returned by EFetch, PubMed database (second test)."""
        # In PubMed display PMIDs in xml retrieval mode.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', id="11748933,11700088",
        #                       retmode="xml")
        with open("Entrez/pubmed2.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        self.assertEqual(record["PubmedBookArticle"], [])
        record = record["PubmedArticle"]
        self.assertEqual(record[0]["MedlineCitation"].attributes["Owner"], "NLM")
        self.assertEqual(record[0]["MedlineCitation"].attributes["Status"], "MEDLINE")
        self.assertEqual(record[0]["MedlineCitation"]["PMID"], "11748933")
        self.assertEqual(record[0]["MedlineCitation"]["DateCompleted"]["Year"], "2002")
        self.assertEqual(record[0]["MedlineCitation"]["DateCompleted"]["Month"], "03")
        self.assertEqual(record[0]["MedlineCitation"]["DateCompleted"]["Day"], "04")
        self.assertEqual(record[0]["MedlineCitation"]["DateRevised"]["Year"], "2006")
        self.assertEqual(record[0]["MedlineCitation"]["DateRevised"]["Month"], "11")
        self.assertEqual(record[0]["MedlineCitation"]["DateRevised"]["Day"], "15")
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"].attributes["PubModel"], "Print"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"], "0011-2240"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes[
                "IssnType"
            ],
            "Print",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"],
            "Print",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Volume"
            ],
            "42",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"],
            "4",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Year"],
            "2001",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Month"],
            "Jun",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["Title"], "Cryobiology"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"],
            "Cryobiology",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["ArticleTitle"],
            "Is cryopreservation a homogeneous process? Ultrastructure and motility of untreated, prefreezing, and postthawed spermatozoa of Diplodus puntazzo (Cetti).",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"],
            "244-55",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"],
            ["This study subdivides the cryopreservation procedure for Diplodus puntazzo spermatozoa into three key phases, fresh, prefreezing (samples equilibrated in cryosolutions), and postthawed stages, and examines the ultrastructural anomalies and motility profiles of spermatozoa in each stage, with different cryodiluents. Two simple cryosolutions were evaluated: 0.17 M sodium chloride containing a final concentration of 15% dimethyl sulfoxide (Me(2)SO) (cryosolution A) and 0.1 M sodium citrate containing a final concentration of 10% Me(2)SO (cryosolution B). Ultrastructural anomalies of the plasmatic and nuclear membranes of the sperm head were common and the severity of the cryoinjury differed significantly between the pre- and the postfreezing phases and between the two cryosolutions. In spermatozoa diluted with cryosolution A, during the prefreezing phase, the plasmalemma of 61% of the cells was absent or damaged compared with 24% in the fresh sample (P < 0.001). In spermatozoa diluted with cryosolution B, there was a pronounced increase in the number of cells lacking the head plasmatic membrane from the prefreezing to the postthawed stages (from 32 to 52%, P < 0.01). In both cryosolutions, damages to nuclear membrane were significantly higher after freezing (cryosolution A: 8 to 23%, P < 0.01; cryosolution B: 5 to 38%, P < 0.001). With cryosolution A, the after-activation motility profile confirmed a consistent drop from fresh at the prefreezing stage, whereas freezing and thawing did not affect the motility much further and 50% of the cells were immotile by 60-90 s after activation. With cryosolution B, only the postthawing stage showed a sharp drop of motility profile. This study suggests that the different phases of the cryoprocess should be investigated to better understand the process of sperm damage.",]
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"],
            "Copyright 2001 Elsevier Science.",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"],
            "Taddei",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"], "A R"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"], "AR"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][0]["Affiliation"],
            "Dipartimento di Scienze Ambientali, Universit\xe0 degli Studi della Tuscia, 01100 Viterbo, Italy.",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][0]["Identifier"],
            [],
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][1].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"],
            "Barbato",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"], "F"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"], "F"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][2].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"],
            "Abelli",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"], "L"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"], "L"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][3].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"],
            "Canese",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"], "S"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"], "S"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][4].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"], "Moretti",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"], "F"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"], "F"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][5].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"], "Rana"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"], "K J"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"], "KJ"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][6].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["LastName"],
            "Fausto",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["ForeName"], "A M"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["Initials"], "AM"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][7].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["LastName"],
            "Mazzini",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["ForeName"], "M"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["Initials"], "M"
        )
        self.assertEqual(record[0]["MedlineCitation"]["Article"]["Language"], ["eng"])
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["PublicationTypeList"][0],
            "Journal Article",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["Article"]["PublicationTypeList"][1],
            "Research Support, Non-U.S. Gov't",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MedlineJournalInfo"]["Country"],
            "Netherlands",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"],
            "Cryobiology",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"], "0006252"
        )
        self.assertEqual(record[0]["MedlineCitation"]["CitationSubset"], ["IM"])
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"],
            "Animals",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][0][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"],
            "Cell Membrane",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][1][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0],
            "ultrastructure",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][
                0
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"],
            "Cryopreservation",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][2][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][2]["QualifierName"][0],
            "methods",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][2]["QualifierName"][
                0
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"], "Male"
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][3][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"],
            "Microscopy, Electron",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][4][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"],
            "Microscopy, Electron, Scanning",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][5][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"],
            "Nuclear Envelope",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][6][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][6]["QualifierName"][0],
            "ultrastructure",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][6]["QualifierName"][
                0
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"],
            "Sea Bream",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][7][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0],
            "anatomy & histology",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][
                0
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][1],
            "physiology",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][
                1
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"],
            "Semen Preservation",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][8][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][0],
            "adverse effects",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][
                0
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][1],
            "methods",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][
                1
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"],
            "Sperm Motility",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][9][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"],
            "Spermatozoa",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][10][
                "DescriptorName"
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][0],
            "physiology",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][
                0
            ].attributes["MajorTopicYN"],
            "N",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][1],
            "ultrastructure",
        )
        self.assertEqual(
            record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][
                1
            ].attributes["MajorTopicYN"],
            "Y",
        )
        self.assertEqual(
            record[0]["PubmedData"]["History"][0].attributes["PubStatus"], "pubmed"
        )
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Year"], "2001")
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Month"], "12")
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Day"], "26")
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Hour"], "10")
        self.assertEqual(record[0]["PubmedData"]["History"][0]["Minute"], "0")
        self.assertEqual(
            record[0]["PubmedData"]["History"][1].attributes["PubStatus"], "medline"
        )
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Year"], "2002")
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Month"], "3")
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Day"], "5")
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Hour"], "10")
        self.assertEqual(record[0]["PubmedData"]["History"][1]["Minute"], "1")
        self.assertEqual(record[0]["PubmedData"]["PublicationStatus"], "ppublish")
        self.assertEqual(record[0]["PubmedData"]["ArticleIdList"][0], "11748933")
        self.assertEqual(
            record[0]["PubmedData"]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            record[0]["PubmedData"]["ArticleIdList"][1], "10.1006/cryo.2001.2328"
        )
        self.assertEqual(
            record[0]["PubmedData"]["ArticleIdList"][1].attributes["IdType"], "doi"
        )
        self.assertEqual(
            record[0]["PubmedData"]["ArticleIdList"][2], "S0011-2240(01)92328-4"
        )
        self.assertEqual(
            record[0]["PubmedData"]["ArticleIdList"][2].attributes["IdType"], "pii"
        )

        self.assertEqual(record[1]["MedlineCitation"].attributes["Owner"], "NLM")
        self.assertEqual(
            record[1]["MedlineCitation"].attributes["Status"], "PubMed-not-MEDLINE"
        )
        self.assertEqual(record[1]["MedlineCitation"]["PMID"], "11700088")
        self.assertEqual(record[1]["MedlineCitation"]["DateCompleted"]["Year"], "2001")
        self.assertEqual(record[1]["MedlineCitation"]["DateCompleted"]["Month"], "12")
        self.assertEqual(record[1]["MedlineCitation"]["DateCompleted"]["Day"], "20")
        self.assertEqual(record[1]["MedlineCitation"]["DateRevised"]["Year"], "2003")
        self.assertEqual(record[1]["MedlineCitation"]["DateRevised"]["Month"], "10")
        self.assertEqual(record[1]["MedlineCitation"]["DateRevised"]["Day"], "31")
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"].attributes["PubModel"], "Print"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"], "1090-7807"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes[
                "IssnType"
            ],
            "Print",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"],
            "Print",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Volume"
            ],
            "153",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"],
            "1",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Year"],
            "2001",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Month"],
            "Nov",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["Title"],
            "Journal of magnetic resonance (San Diego, Calif. : 1997)",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"],
            "J Magn Reson",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["ArticleTitle"],
            "Proton MRI of (13)C distribution by J and chemical shift editing.",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"],
            "117-23",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"],
            ["The sensitivity of (13)C NMR imaging can be considerably favored by detecting the (1)H nuclei bound to (13)C nuclei via scalar J-interaction (X-filter). However, the J-editing approaches have difficulty in discriminating between compounds with similar J-constant as, for example, different glucose metabolites. In such cases, it is almost impossible to get J-edited images of a single-compound distribution, since the various molecules are distinguishable only via their chemical shift. In a recent application of J-editing to high-resolution spectroscopy, it has been shown that a more efficient chemical selectivity could be obtained by utilizing the larger chemical shift range of (13)C. This has been made by introducing frequency-selective (13)C pulses that allow a great capability of indirect chemical separation. Here a double-resonance imaging approach is proposed, based on both J-editing and (13)C chemical shift editing, which achieves a powerful chemical selectivity and is able to produce full maps of specific chemical compounds. Results are presented on a multicompartments sample containing solutions of glucose and lactic and glutamic acid in water."],
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"],
            "Copyright 2001 Academic Press.",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"],
            "Casieri",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"], "C"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"], "C"
        )
        self.assertEqual(record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["Identifier"], [])
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][0]["Affiliation"],
            "INFM and Department of Physics, University of L'Aquila, I-67100 L'Aquila, Italy.",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][1].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"],
            "Testa",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"], "C"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"], "C"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][2].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"],
            "Carpinelli",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"], "G"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"], "G"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][3].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"],
            "Canese",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"], "R"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"], "R"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][4].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"], "Podo"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"], "F"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"], "F"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][5].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"],
            "De Luca",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"], "F"
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"], "F"
        )
        self.assertEqual(record[1]["MedlineCitation"]["Article"]["Language"], ["eng"])
        self.assertEqual(
            record[1]["MedlineCitation"]["Article"]["PublicationTypeList"][0],
            "Journal Article",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MedlineJournalInfo"]["Country"],
            "United States",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"],
            "J Magn Reson",
        )
        self.assertEqual(
            record[1]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"], "9707935"
        )
        self.assertEqual(
            record[1]["PubmedData"]["History"][0].attributes["PubStatus"], "pubmed"
        )
        self.assertEqual(record[1]["PubmedData"]["History"][0]["Year"], "2001")
        self.assertEqual(record[1]["PubmedData"]["History"][0]["Month"], "11")
        self.assertEqual(record[1]["PubmedData"]["History"][0]["Day"], "9")
        self.assertEqual(record[1]["PubmedData"]["History"][0]["Hour"], "10")
        self.assertEqual(record[1]["PubmedData"]["History"][0]["Minute"], "0")
        self.assertEqual(
            record[1]["PubmedData"]["History"][1].attributes["PubStatus"], "medline"
        )
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Year"], "2001")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Month"], "11")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Day"], "9")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Hour"], "10")
        self.assertEqual(record[1]["PubmedData"]["History"][1]["Minute"], "1")
        self.assertEqual(record[1]["PubmedData"]["PublicationStatus"], "ppublish")
        self.assertEqual(record[1]["PubmedData"]["ArticleIdList"][0], "11700088")
        self.assertEqual(
            record[1]["PubmedData"]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            record[1]["PubmedData"]["ArticleIdList"][1], "10.1006/jmre.2001.2429"
        )
        self.assertEqual(
            record[1]["PubmedData"]["ArticleIdList"][1].attributes["IdType"], "doi"
        )
        self.assertEqual(
            record[1]["PubmedData"]["ArticleIdList"][2], "S1090-7807(01)92429-2"
        )
        self.assertEqual(
            record[1]["PubmedData"]["ArticleIdList"][2].attributes["IdType"], "pii"
        )
        # fmt: on

    def test_pubmed_html_tags(self):
        """Test parsing XML returned by EFetch, PubMed database with HTML tags."""
        # In PubMed display PMIDs in xml retrieval mode.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', retmode='xml', id='29106400')
        with open("Entrez/pubmed4.xml", "rb") as stream:
            records = Entrez.read(stream)
        # fmt: off
        self.assertEqual(len(records), 2)
        self.assertEqual(len(records["PubmedBookArticle"]), 0)
        self.assertEqual(len(records["PubmedArticle"]), 1)
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"].attributes["Status"],
            "MEDLINE",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"].attributes["Owner"], "NLM"
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["PMID"], "27797938"
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["PMID"].attributes[
                "Version"
            ],
            "1",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["DateCompleted"]["Year"],
            "2017",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["DateCompleted"]["Month"],
            "08",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["DateCompleted"]["Day"], "03"
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["DateRevised"]["Year"],
            "2018",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["DateRevised"]["Month"], "04"
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["DateRevised"]["Day"], "17"
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"].attributes[
                "PubModel"
            ],
            "Print-Electronic",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "ISSN"
            ],
            "1468-3288",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "ISSN"
            ].attributes["IssnType"],
            "Electronic",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"],
            "Internet",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ]["Volume"],
            "66",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ]["Issue"],
            "6",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ]["PubDate"]["Year"],
            "2017",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ]["PubDate"]["Month"],
            "06",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "Title"
            ],
            "Gut",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"][
                "ISOAbbreviation"
            ],
            "Gut",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"],
            "Leucocyte telomere length, genetic variants at the <i>TERT</i> gene region and risk of pancreatic cancer.",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Pagination"][
                "MedlinePgn"
            ],
            "1116-1122",
        )
        self.assertEqual(
            len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ELocationID"]
            ),
            1,
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ELocationID"][0],
            "10.1136/gutjnl-2016-312510",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ELocationID"][
                0
            ].attributes["EIdType"],
            "doi",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ELocationID"][
                0
            ].attributes["ValidYN"],
            "Y",
        )
        self.assertEqual(
            len(records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]),
            2,
        )
        self.assertEqual(
            len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                    "AbstractText"
                ]
            ),
            4,
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][0],
            "Telomere shortening occurs as an early event in pancreatic tumorigenesis, and genetic variants at the telomerase reverse transcriptase (<i>TERT</i>) gene region have been associated with pancreatic cancer risk. However, it is unknown whether prediagnostic leucocyte telomere length is associated with subsequent risk of pancreatic cancer.",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][0].attributes["Label"],
            "OBJECTIVE",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][1],
            "We measured prediagnostic leucocyte telomere length in 386 pancreatic cancer cases and 896 matched controls from five prospective US cohorts. ORs and 95% CIs were calculated using conditional logistic regression. Matching factors included year of birth, cohort (which also matches on sex), smoking status, fasting status and month/year of blood collection. We additionally examined single-nucleotide polymorphisms (SNPs) at the <i>TERT</i> region in relation to pancreatic cancer risk and leucocyte telomere length using logistic and linear regression, respectively.",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][1].attributes["Label"],
            "DESIGN",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][2],
            "Shorter prediagnostic leucocyte telomere length was associated with higher risk of pancreatic cancer (comparing extreme quintiles of telomere length, OR 1.72; 95% CI 1.07 to 2.78; p<sub>trend</sub>=0.048). Results remained unchanged after adjustment for diabetes, body mass index and physical activity. Three SNPs at <i>TERT</i> (linkage disequilibrium r<sup>2</sup><0.25) were associated with pancreatic cancer risk, including rs401681 (per minor allele OR 1.33; 95% CI 1.12 to 1.59; p=0.002), rs2736100 (per minor allele OR 1.36; 95% CI 1.13 to 1.63; p=0.001) and rs2736098 (per minor allele OR 0.75; 95% CI 0.63 to 0.90; p=0.002). The minor allele for rs401681 was associated with shorter telomere length (p=0.023).",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][2].attributes["Label"],
            "RESULTS",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][3],
            "Prediagnostic leucocyte telomere length and genetic variants at the <i>TERT</i> gene region were associated with risk of pancreatic cancer.",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"][
                "AbstractText"
            ][3].attributes["Label"],
            "CONCLUSIONS",
        )
        self.assertEqual(
            len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"]["AuthorList"]
            ),
            22,
        )
        self.assertEqual(
            len(records["PubmedArticle"][0]["MedlineCitation"]["Article"]["GrantList"]),
            35,
        )
        self.assertEqual(
            len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                    "PublicationTypeList"
                ]
            ),
            5,
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][0],
            "Journal Article",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][0].attributes["UI"],
            "D016428",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][1],
            "Observational Study",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][1].attributes["UI"],
            "D064888",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][2],
            "Research Support, N.I.H., Extramural",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][2].attributes["UI"],
            "D052061",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][3],
            "Research Support, U.S. Gov't, Non-P.H.S.",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][3].attributes["UI"],
            "D013486",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][4],
            "Research Support, Non-U.S. Gov't",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"][
                "PublicationTypeList"
            ][4].attributes["UI"],
            "D013485",
        )
        self.assertEqual(
            len(
                records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"]
            ),
            1,
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][
                0
            ].attributes["DateType"],
            "Electronic",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0][
                "Year"
            ],
            "2016",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0][
                "Month"
            ],
            "10",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0][
                "Day"
            ],
            "21",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["MedlineJournalInfo"][
                "Country"
            ],
            "England",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["MedlineJournalInfo"][
                "MedlineTA"
            ],
            "Gut",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["MedlineJournalInfo"][
                "NlmUniqueID"
            ],
            "2985108R",
        )
        self.assertEqual(
            records["PubmedArticle"][0]["MedlineCitation"]["MedlineJournalInfo"][
                "ISSNLinking"
            ],
            "0017-5749",
        )
        self.assertEqual(
            len(records["PubmedArticle"][0]["MedlineCitation"]["ChemicalList"]), 2
        )
        # fmt: on

    def test_pubmed_html_escaping(self):
        """Test parsing XML returned by EFetch, PubMed database with HTML tags and HTML escape characters."""
        # In PubMed display PMIDs in xml retrieval mode.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', retmode='xml', id='28775130')
        with open("Entrez/pubmed5.xml", "rb") as stream:
            record = Entrez.read(stream, escape=True)
        self.assertEqual(len(record), 2)
        self.assertEqual(len(record["PubmedArticle"]), 1)
        self.assertEqual(len(record["PubmedBookArticle"]), 0)
        article = record["PubmedArticle"][0]
        self.assertEqual(len(article), 2)
        self.assertEqual(len(article["PubmedData"]), 3)
        self.assertEqual(len(article["PubmedData"]["ArticleIdList"]), 5)
        self.assertEqual(article["PubmedData"]["ArticleIdList"][0], "28775130")
        self.assertEqual(
            article["PubmedData"]["ArticleIdList"][0].attributes, {"IdType": "pubmed"}
        )
        self.assertEqual(article["PubmedData"]["ArticleIdList"][1], "oemed-2017-104431")
        self.assertEqual(
            article["PubmedData"]["ArticleIdList"][1].attributes, {"IdType": "pii"}
        )
        self.assertEqual(
            article["PubmedData"]["ArticleIdList"][2], "10.1136/oemed-2017-104431"
        )
        self.assertEqual(
            article["PubmedData"]["ArticleIdList"][2].attributes, {"IdType": "doi"}
        )
        self.assertEqual(article["PubmedData"]["ArticleIdList"][3], "PMC5771820")
        self.assertEqual(
            article["PubmedData"]["ArticleIdList"][3].attributes, {"IdType": "pmc"}
        )
        self.assertEqual(article["PubmedData"]["ArticleIdList"][4], "NIHMS932407")
        self.assertEqual(
            article["PubmedData"]["ArticleIdList"][4].attributes, {"IdType": "mid"}
        )
        self.assertEqual(article["PubmedData"]["PublicationStatus"], "ppublish")
        self.assertEqual(len(article["PubmedData"]["History"]), 7)
        self.assertEqual(len(article["PubmedData"]["History"][0]), 3)
        self.assertEqual(article["PubmedData"]["History"][0]["Year"], "2017")
        self.assertEqual(article["PubmedData"]["History"][0]["Month"], "03")
        self.assertEqual(article["PubmedData"]["History"][0]["Day"], "10")
        self.assertEqual(
            article["PubmedData"]["History"][0].attributes, {"PubStatus": "received"}
        )
        self.assertEqual(len(article["PubmedData"]["History"][1]), 3)
        self.assertEqual(article["PubmedData"]["History"][1]["Year"], "2017")
        self.assertEqual(article["PubmedData"]["History"][1]["Month"], "06")
        self.assertEqual(article["PubmedData"]["History"][1]["Day"], "13")
        self.assertEqual(
            article["PubmedData"]["History"][1].attributes, {"PubStatus": "revised"}
        )
        self.assertEqual(len(article["PubmedData"]["History"][2]), 3)
        self.assertEqual(article["PubmedData"]["History"][2]["Year"], "2017")
        self.assertEqual(article["PubmedData"]["History"][2]["Month"], "06")
        self.assertEqual(article["PubmedData"]["History"][2]["Day"], "22")
        self.assertEqual(
            article["PubmedData"]["History"][2].attributes, {"PubStatus": "accepted"}
        )
        self.assertEqual(len(article["PubmedData"]["History"][3]), 3)
        self.assertEqual(article["PubmedData"]["History"][3]["Year"], "2019")
        self.assertEqual(article["PubmedData"]["History"][3]["Month"], "02")
        self.assertEqual(article["PubmedData"]["History"][3]["Day"], "01")
        self.assertEqual(
            article["PubmedData"]["History"][3].attributes, {"PubStatus": "pmc-release"}
        )
        self.assertEqual(len(article["PubmedData"]["History"][4]), 5)
        self.assertEqual(article["PubmedData"]["History"][4]["Year"], "2017")
        self.assertEqual(article["PubmedData"]["History"][4]["Month"], "8")
        self.assertEqual(article["PubmedData"]["History"][4]["Day"], "5")
        self.assertEqual(article["PubmedData"]["History"][4]["Hour"], "6")
        self.assertEqual(article["PubmedData"]["History"][4]["Minute"], "0")
        self.assertEqual(
            article["PubmedData"]["History"][4].attributes, {"PubStatus": "pubmed"}
        )
        self.assertEqual(len(article["PubmedData"]["History"][5]), 5)
        self.assertEqual(article["PubmedData"]["History"][5]["Year"], "2017")
        self.assertEqual(article["PubmedData"]["History"][5]["Month"], "8")
        self.assertEqual(article["PubmedData"]["History"][5]["Day"], "5")
        self.assertEqual(article["PubmedData"]["History"][5]["Hour"], "6")
        self.assertEqual(article["PubmedData"]["History"][5]["Minute"], "0")
        self.assertEqual(
            article["PubmedData"]["History"][5].attributes, {"PubStatus": "medline"}
        )
        self.assertEqual(len(article["PubmedData"]["History"][6]), 5)
        self.assertEqual(article["PubmedData"]["History"][6]["Year"], "2017")
        self.assertEqual(article["PubmedData"]["History"][6]["Month"], "8")
        self.assertEqual(article["PubmedData"]["History"][6]["Day"], "5")
        self.assertEqual(article["PubmedData"]["History"][6]["Hour"], "6")
        self.assertEqual(article["PubmedData"]["History"][6]["Minute"], "0")
        self.assertEqual(
            article["PubmedData"]["History"][6].attributes, {"PubStatus": "entrez"}
        )
        self.assertEqual(len(article["MedlineCitation"]), 12)
        self.assertEqual(len(article["MedlineCitation"]["CitationSubset"]), 0)
        self.assertEqual(
            article["MedlineCitation"]["CoiStatement"],
            "Competing interests: None declared.",
        )
        self.assertEqual(len(article["MedlineCitation"]["CommentsCorrectionsList"]), 40)
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][0]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][0]["RefSource"],
            "J Toxicol Environ Health A. 2003 Jun 13;66(11):965-86",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][0]["PMID"], "12775511"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][0].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][1]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][1]["RefSource"],
            "Ann Intern Med. 2015 May 5;162(9):641-50",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][1]["PMID"], "25798805"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][1].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][2]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][2]["RefSource"],
            "Cancer Causes Control. 1999 Dec;10(6):583-95",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][2]["PMID"], "10616827"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][2].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][3]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][3]["RefSource"],
            "Thyroid. 2010 Jul;20(7):755-61",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][3]["PMID"], "20578899"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][3].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][4]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][4]["RefSource"],
            "Environ Health Perspect. 1999 Mar;107(3):205-11",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][4]["PMID"], "10064550"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][4].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][5]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][5]["RefSource"],
            "J Clin Endocrinol Metab. 2006 Nov;91(11):4295-301",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][5]["PMID"], "16868053"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][5].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][6]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][6]["RefSource"],
            "Endocrinology. 1998 Oct;139(10):4252-63",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][6]["PMID"], "9751507"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][6].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][7]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][7]["RefSource"],
            "Eur J Endocrinol. 2016 Apr;174(4):409-14",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][7]["PMID"], "26863886"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][7].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][8]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][8]["RefSource"],
            "Eur J Endocrinol. 2000 Nov;143(5):639-47",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][8]["PMID"], "11078988"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][8].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][9]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][9]["RefSource"],
            "Environ Res. 2016 Nov;151:389-398",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][9]["PMID"], "27540871"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][9].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][10]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][10]["RefSource"],
            "Am J Epidemiol. 2010 Jan 15;171(2):242-52",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][10]["PMID"],
            "19951937",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][10].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][11]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][11]["RefSource"],
            "Thyroid. 1998 Sep;8(9):827-56",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][11]["PMID"], "9777756"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][11].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][12]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][12]["RefSource"],
            "Curr Opin Pharmacol. 2001 Dec;1(6):626-31",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][12]["PMID"],
            "11757819",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][12].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][13]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][13]["RefSource"],
            "Breast Cancer Res Treat. 2012 Jun;133(3):1169-77",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][13]["PMID"],
            "22434524",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][13].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][14]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][14]["RefSource"],
            "Int J Environ Res Public Health. 2011 Dec;8(12 ):4608-22",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][14]["PMID"],
            "22408592",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][14].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][15]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][15]["RefSource"],
            "Ann Oncol. 2014 Oct;25(10):2025-30",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][15]["PMID"],
            "25081899",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][15].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][16]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][16]["RefSource"],
            "Environ Health. 2006 Dec 06;5:32",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][16]["PMID"],
            "17147831",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][16].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][17]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][17]["RefSource"],
            "Environ Health Perspect. 1998 Aug;106(8):437-45",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][17]["PMID"], "9681970"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][17].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][18]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][18]["RefSource"],
            "Arch Intern Med. 2000 Feb 28;160(4):526-34",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][18]["PMID"],
            "10695693",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][18].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][19]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][19]["RefSource"],
            "Endocrine. 2011 Jun;39(3):259-65",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][19]["PMID"],
            "21161440",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][19].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][20]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][20]["RefSource"],
            "Cancer Epidemiol Biomarkers Prev. 2008 Aug;17(8):1880-3",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][20]["PMID"],
            "18708375",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][20].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][21]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][21]["RefSource"],
            "Am J Epidemiol. 2010 Feb 15;171(4):455-64",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][21]["PMID"],
            "20061368",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][21].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][22]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][22]["RefSource"],
            "J Clin Endocrinol Metab. 2002 Feb;87(2):489-99",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][22]["PMID"],
            "11836274",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][22].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][23]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][23]["RefSource"],
            "J Toxicol Environ Health A. 2015 ;78(21-22):1338-47",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][23]["PMID"],
            "26555155",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][23].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][24]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][24]["RefSource"],
            "Toxicol Sci. 2002 Jun;67(2):207-18",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][24]["PMID"],
            "12011480",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][24].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][25]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][25]["RefSource"],
            "Natl Cancer Inst Carcinog Tech Rep Ser. 1978;21:1-184",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][25]["PMID"],
            "12844187",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][25].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][26]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][26]["RefSource"],
            "Environ Res. 2013 Nov;127:7-15",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][26]["PMID"],
            "24183346",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][26].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][27]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][27]["RefSource"],
            "JAMA. 2004 Jan 14;291(2):228-38",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][27]["PMID"],
            "14722150",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][27].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][28]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][28]["RefSource"],
            "J Expo Sci Environ Epidemiol. 2010 Sep;20(6):559-69",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][28]["PMID"],
            "19888312",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][28].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][29]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][29]["RefSource"],
            "Environ Health Perspect. 1996 Apr;104(4):362-9",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][29]["PMID"], "8732939"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][29].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][30]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][30]["RefSource"],
            "Lancet. 2012 Mar 24;379(9821):1142-54",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][30]["PMID"],
            "22273398",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][30].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][31]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][31]["RefSource"],
            "JAMA. 1995 Mar 8;273(10):808-12",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][31]["PMID"], "7532241"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][31].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][32]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][32]["RefSource"],
            "Sci Total Environ. 2002 Aug 5;295(1-3):207-15",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][32]["PMID"],
            "12186288",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][32].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][33]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][33]["RefSource"],
            "Eur J Endocrinol. 2006 May;154(5):599-611",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][33]["PMID"],
            "16645005",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][33].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][34]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][34]["RefSource"],
            "J Occup Environ Med. 2013 Oct;55(10):1171-8",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][34]["PMID"],
            "24064777",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][34].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][35]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][35]["RefSource"],
            "Thyroid. 2007 Sep;17(9):811-7",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][35]["PMID"],
            "17956155",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][35].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][36]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][36]["RefSource"],
            "Rev Environ Contam Toxicol. 1991;120:1-82",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][36]["PMID"], "1899728"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][36].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][37]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][37]["RefSource"],
            "Environ Health Perspect. 1997 Oct;105(10):1126-30",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][37]["PMID"], "9349837"
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][37].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][38]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][38]["RefSource"],
            "J Biochem Mol Toxicol. 2005;19(3):175",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][38]["PMID"],
            "15977190",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][38].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(
            len(article["MedlineCitation"]["CommentsCorrectionsList"][39]), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][39]["RefSource"],
            "Immunogenetics. 2002 Jun;54(3):141-57",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][39]["PMID"],
            "12073143",
        )
        self.assertEqual(
            article["MedlineCitation"]["CommentsCorrectionsList"][39].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(article["MedlineCitation"]["DateRevised"]["Year"], "2018")
        self.assertEqual(article["MedlineCitation"]["DateRevised"]["Month"], "04")
        self.assertEqual(article["MedlineCitation"]["DateRevised"]["Day"], "25")
        self.assertEqual(len(article["MedlineCitation"]["DateRevised"].attributes), 0)
        self.assertEqual(len(article["MedlineCitation"]["GeneralNote"]), 0)
        self.assertEqual(len(article["MedlineCitation"]["KeywordList"]), 1)
        self.assertEqual(len(article["MedlineCitation"]["KeywordList"][0]), 5)
        self.assertEqual(article["MedlineCitation"]["KeywordList"][0][0], "agriculture")
        self.assertEqual(
            article["MedlineCitation"]["KeywordList"][0][0].attributes,
            {"MajorTopicYN": "N"},
        )
        self.assertEqual(
            article["MedlineCitation"]["KeywordList"][0][1], "hypothyroidism"
        )
        self.assertEqual(
            article["MedlineCitation"]["KeywordList"][0][1].attributes,
            {"MajorTopicYN": "N"},
        )
        self.assertEqual(article["MedlineCitation"]["KeywordList"][0][2], "pesticides")
        self.assertEqual(
            article["MedlineCitation"]["KeywordList"][0][2].attributes,
            {"MajorTopicYN": "N"},
        )
        self.assertEqual(
            article["MedlineCitation"]["KeywordList"][0][3], "thyroid disease"
        )
        self.assertEqual(
            article["MedlineCitation"]["KeywordList"][0][3].attributes,
            {"MajorTopicYN": "N"},
        )
        self.assertEqual(
            article["MedlineCitation"]["KeywordList"][0][4],
            "thyroid stimulating hormone",
        )
        self.assertEqual(
            article["MedlineCitation"]["KeywordList"][0][4].attributes,
            {"MajorTopicYN": "N"},
        )
        self.assertEqual(len(article["MedlineCitation"]["MedlineJournalInfo"]), 4)
        self.assertEqual(
            article["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"],
            "Occup Environ Med",
        )
        self.assertEqual(
            article["MedlineCitation"]["MedlineJournalInfo"]["Country"], "England"
        )
        self.assertEqual(
            article["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"], "9422759"
        )
        self.assertEqual(
            article["MedlineCitation"]["MedlineJournalInfo"]["ISSNLinking"], "1351-0711"
        )
        self.assertEqual(
            len(article["MedlineCitation"]["MedlineJournalInfo"].attributes), 0
        )
        self.assertEqual(len(article["MedlineCitation"]["OtherAbstract"]), 0)
        self.assertEqual(len(article["MedlineCitation"]["OtherID"]), 0)
        self.assertEqual(article["MedlineCitation"]["PMID"], "28775130")
        self.assertEqual(len(article["MedlineCitation"]["SpaceFlightMission"]), 0)
        self.assertEqual(len(article["MedlineCitation"]["Article"]["ArticleDate"]), 1)
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["ArticleDate"][0]), 3
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["ArticleDate"][0]["Month"], "08"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["ArticleDate"][0]["Day"], "03"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"], "2017"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["ArticleDate"][0].attributes,
            {"DateType": "Electronic"},
        )
        self.assertEqual(len(article["MedlineCitation"]["Article"]["Pagination"]), 1)
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"], "79-89"
        )
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["Pagination"].attributes), 0
        )
        self.assertEqual(len(article["MedlineCitation"]["Article"]["AuthorList"]), 12)
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"], "Lerro"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"], "CC"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][0]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][0][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][0]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][0][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"],
            "Catherine C",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][0].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"],
            "Beane Freeman",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"], "LE"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][1]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][1][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][1]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][1]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][1][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"],
            "Laura E",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][1].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"],
            "DellaValle",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"], "CT"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][2][
                    "AffiliationInfo"
                ]
            ),
            2,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][2][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2]["AffiliationInfo"][
                1
            ]["Affiliation"],
            "Environmental Working Group, Washington, DC, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2]["AffiliationInfo"][
                1
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][2][
                    "AffiliationInfo"
                ][1].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"], "Curt T"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][2].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"],
            "Kibriya",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"], "MG"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][3]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][3][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][3]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Department of Public Health Sciences, The University of Chicago, Chicago, Illinois, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][3]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][3][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"],
            "Muhammad G",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][3].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"],
            "Aschebrook-Kilfoy",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"], "B"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][4]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][4][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][4]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Department of Public Health Sciences, The University of Chicago, Chicago, Illinois, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][4]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][4][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"],
            "Briseis",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][4].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"],
            "Jasmine",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"], "F"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][5]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][5][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][5]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Department of Public Health Sciences, The University of Chicago, Chicago, Illinois, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][5]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][5][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"],
            "Farzana",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][5].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][6]["LastName"],
            "Koutros",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][6]["Initials"], "S"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][6]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][6][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][6]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][6]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][6][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][6]["ForeName"], "Stella"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][6].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][7]["LastName"], "Parks"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][7]["Initials"], "CG"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][7]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][7][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][7]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][7]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][7][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][7]["ForeName"],
            "Christine G",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][7].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][8]["LastName"],
            "Sandler",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][8]["Initials"], "DP"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][8]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][8][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][8]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "National Institute of Environmental Health Sciences, Research Triangle Park, North Carolina, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][8]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][8][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][8]["ForeName"], "Dale P"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][8].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9]["LastName"],
            "Alavanja",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9]["Initials"], "MCR"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][9][
                    "AffiliationInfo"
                ]
            ),
            2,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][9][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9]["AffiliationInfo"][
                1
            ]["Affiliation"],
            "Department of Biology, Hood College, Frederick, Maryland, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9]["AffiliationInfo"][
                1
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][9][
                    "AffiliationInfo"
                ][1].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9]["ForeName"],
            "Michael C R",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][9].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][10]["LastName"],
            "Hofmann",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][10]["Initials"], "JN"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][10]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][10][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][10]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][10]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][10][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][10]["ForeName"],
            "Jonathan N",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][10].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][11]["LastName"], "Ward"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][11]["Initials"], "MH"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][11]["Identifier"], []
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][11][
                    "AffiliationInfo"
                ]
            ),
            1,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][11]["AffiliationInfo"][
                0
            ]["Affiliation"],
            "Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, Maryland, USA.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][11]["AffiliationInfo"][
                0
            ]["Identifier"],
            [],
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["AuthorList"][11][
                    "AffiliationInfo"
                ][0].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][11]["ForeName"],
            "Mary H",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["AuthorList"][11].attributes,
            {"ValidYN": "Y"},
        )
        self.assertEqual(len(article["MedlineCitation"]["Article"]["Language"]), 1)
        self.assertEqual(article["MedlineCitation"]["Article"]["Language"][0], "eng")
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["PublicationTypeList"]), 1
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["PublicationTypeList"][0],
            "Journal Article",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["PublicationTypeList"][0].attributes,
            {"UI": "D016428"},
        )
        self.assertEqual(len(article["MedlineCitation"]["Article"]["Journal"]), 4)
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["ISSN"], "1470-7926"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes,
            {"IssnType": "Electronic"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"],
            "Occup Environ Med",
        )
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]), 3
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"],
            "75",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"],
            "2",
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                    "PubDate"
                ]
            ),
            2,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"][
                "Month"
            ],
            "Feb",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"][
                "Year"
            ],
            "2018",
        )
        self.assertEqual(
            len(
                article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                    "PubDate"
                ].attributes
            ),
            0,
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].attributes,
            {"CitedMedium": "Internet"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["Title"],
            "Occupational and environmental medicine",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes,
            {"IssnType": "Electronic"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["ArticleTitle"],
            "Occupational pesticide exposure and subclinical hypothyroidism among male pesticide applicators.",
        )
        self.assertEqual(len(article["MedlineCitation"]["Article"]["ELocationID"]), 1)
        self.assertEqual(
            article["MedlineCitation"]["Article"]["ELocationID"][0],
            "10.1136/oemed-2017-104431",
        )
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["ELocationID"][0].attributes), 2
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["ELocationID"][0].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["ELocationID"][0].attributes[
                "EIdType"
            ],
            "doi",
        )
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]), 4
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0],
            "Animal studies suggest that exposure to pesticides may alter thyroid function; however, few epidemiologic studies have examined this association. We evaluated the relationship between individual pesticides and thyroid function in 679 men enrolled in a substudy of the Agricultural Health Study, a cohort of licensed pesticide applicators.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][
                0
            ].attributes,
            {"NlmCategory": "OBJECTIVE", "Label": "OBJECTIVES"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][1],
            "Self-reported lifetime pesticide use was obtained at cohort enrolment (1993-1997). Intensity-weighted lifetime days were computed for 33 pesticides, which adjusts cumulative days of pesticide use for factors that modify exposure (eg, use of personal protective equipment). Thyroid-stimulating hormone (TSH), thyroxine (T4), triiodothyronine (T3) and antithyroid peroxidase (anti-TPO) autoantibodies were measured in serum collected in 2010-2013. We used multivariate logistic regression to estimate ORs and 95% CIs for subclinical hypothyroidism (TSH &gt;4.5 mIU/L) compared with normal TSH (0.4-<u>&lt;</u>4.5 mIU/L) and for anti-TPO positivity. We also examined pesticide associations with TSH, T4 and T3 in multivariate linear regression models.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][
                1
            ].attributes,
            {"NlmCategory": "METHODS", "Label": "METHODS"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][2],
            "Higher exposure to the insecticide aldrin (third and fourth quartiles of intensity-weighted days vs no exposure) was positively associated with subclinical hypothyroidism (OR<sub>Q3</sub>=4.15, 95% CI 1.56 to 11.01, OR<sub>Q4</sub>=4.76, 95% CI 1.53 to 14.82, p<sub>trend</sub> &lt;0.01), higher TSH (p<sub>trend</sub>=0.01) and lower T4 (p<sub>trend</sub>=0.04). Higher exposure to the herbicide pendimethalin was associated with subclinical hypothyroidism (fourth quartile vs no exposure: OR<sub>Q4</sub>=2.78, 95% CI 1.30 to 5.95, p<sub>trend</sub>=0.02), higher TSH (p<sub>trend</sub>=0.04) and anti-TPO positivity (p<sub>trend</sub>=0.01). The fumigant methyl bromide was inversely associated with TSH (p<sub>trend</sub>=0.02) and positively associated with T4 (p<sub>trend</sub>=0.01).",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][
                2
            ].attributes,
            {"NlmCategory": "RESULTS", "Label": "RESULTS"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][3],
            "Our results suggest that long-term exposure to aldrin, pendimethalin and methyl bromide may alter thyroid function among male pesticide applicators.",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][
                3
            ].attributes,
            {"NlmCategory": "CONCLUSIONS", "Label": "CONCLUSIONS"},
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"],
            "\xa9 Article author(s) (or their employer(s) unless otherwise stated in the text of the article) 2018. All rights reserved. No commercial use is permitted unless otherwise expressly granted.",
        )
        self.assertEqual(len(article["MedlineCitation"]["Article"]["GrantList"]), 3)
        self.assertEqual(len(article["MedlineCitation"]["Article"]["GrantList"][0]), 4)
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][0]["Acronym"], "CP"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][0]["Country"],
            "United States",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][0]["Agency"],
            "NCI NIH HHS",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][0]["GrantID"],
            "Z01 CP010119",
        )
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["GrantList"][0].attributes), 0
        )
        self.assertEqual(len(article["MedlineCitation"]["Article"]["GrantList"][1]), 4)
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][1]["Acronym"], "ES"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][1]["Country"],
            "United States",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][1]["Agency"],
            "NIEHS NIH HHS",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][1]["GrantID"],
            "Z01 ES049030",
        )
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["GrantList"][1].attributes), 0
        )
        self.assertEqual(len(article["MedlineCitation"]["Article"]["GrantList"][2]), 4)
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][2]["Acronym"], "NULL"
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][2]["Country"],
            "United States",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][2]["Agency"],
            "Intramural NIH HHS",
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"][2]["GrantID"],
            "Z99 CA999999",
        )
        self.assertEqual(
            len(article["MedlineCitation"]["Article"]["GrantList"][2].attributes), 0
        )
        self.assertEqual(
            article["MedlineCitation"]["Article"]["GrantList"].attributes,
            {"CompleteYN": "Y"},
        )

    def test_pubmed_html_mathml_tags(self):
        """Test parsing XML returned by EFetch, PubMed database, with both HTML and MathML tags."""
        # In PubMed display PMID 30108519 in xml retrieval mode, containing
        # both HTML and MathML tags in the abstract text.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="pubmed", id='30108519', rettype="null",
        #                       retmode="xml", parsed=True)
        with open("Entrez/pubmed6.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 2)
        self.assertEqual(record["PubmedBookArticle"], [])
        self.assertEqual(len(record["PubmedArticle"]), 1)
        pubmed_article = record["PubmedArticle"][0]
        self.assertEqual(len(pubmed_article), 2)
        self.assertEqual(len(pubmed_article["PubmedData"]), 3)
        self.assertEqual(len(pubmed_article["PubmedData"]["ArticleIdList"]), 3)
        self.assertEqual(pubmed_article["PubmedData"]["ArticleIdList"][0], "30108519")
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][0].attributes,
            {"IdType": "pubmed"},
        )
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][1], "10.3389/fphys.2018.01034"
        )
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][1].attributes,
            {"IdType": "doi"},
        )
        self.assertEqual(pubmed_article["PubmedData"]["ArticleIdList"][2], "PMC6079548")
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][2].attributes,
            {"IdType": "pmc"},
        )
        self.assertEqual(pubmed_article["PubmedData"]["PublicationStatus"], "epublish")
        self.assertEqual(len(pubmed_article["PubmedData"]["History"]), 5)
        self.assertEqual(len(pubmed_article["PubmedData"]["History"][0]), 3)
        self.assertEqual(pubmed_article["PubmedData"]["History"][0]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][0]["Month"], "05")
        self.assertEqual(pubmed_article["PubmedData"]["History"][0]["Day"], "22")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][0].attributes,
            {"PubStatus": "received"},
        )
        self.assertEqual(len(pubmed_article["PubmedData"]["History"][1]), 3)
        self.assertEqual(pubmed_article["PubmedData"]["History"][1]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][1]["Month"], "07")
        self.assertEqual(pubmed_article["PubmedData"]["History"][1]["Day"], "11")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][1].attributes,
            {"PubStatus": "accepted"},
        )
        self.assertEqual(len(pubmed_article["PubmedData"]["History"][2]), 5)
        self.assertEqual(pubmed_article["PubmedData"]["History"][2]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][2]["Month"], "8")
        self.assertEqual(pubmed_article["PubmedData"]["History"][2]["Day"], "16")
        self.assertEqual(pubmed_article["PubmedData"]["History"][2]["Hour"], "6")
        self.assertEqual(pubmed_article["PubmedData"]["History"][2]["Minute"], "0")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][2].attributes,
            {"PubStatus": "entrez"},
        )
        self.assertEqual(len(pubmed_article["PubmedData"]["History"][3]), 5)
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Month"], "8")
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Day"], "16")
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Hour"], "6")
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Minute"], "0")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][3].attributes,
            {"PubStatus": "pubmed"},
        )
        self.assertEqual(len(pubmed_article["PubmedData"]["History"][4]), 5)
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Month"], "8")
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Day"], "16")
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Hour"], "6")
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Minute"], "1")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][4].attributes,
            {"PubStatus": "medline"},
        )
        medline_citation = pubmed_article["MedlineCitation"]
        self.assertEqual(len(medline_citation), 11)
        self.assertEqual(medline_citation["GeneralNote"], [])
        self.assertEqual(len(medline_citation["KeywordList"]), 1)
        self.assertEqual(len(medline_citation["KeywordList"][0]), 8)
        self.assertEqual(medline_citation["KeywordList"][0][0], "Owles' point")
        self.assertEqual(
            medline_citation["KeywordList"][0][0].attributes, {"MajorTopicYN": "N"}
        )
        self.assertEqual(medline_citation["KeywordList"][0][1], "aerobic capacity")
        self.assertEqual(
            medline_citation["KeywordList"][0][1].attributes, {"MajorTopicYN": "N"}
        )
        self.assertEqual(medline_citation["KeywordList"][0][2], "aerobic threshold")
        self.assertEqual(
            medline_citation["KeywordList"][0][2].attributes, {"MajorTopicYN": "N"}
        )
        self.assertEqual(medline_citation["KeywordList"][0][3], "anaerobic threshold")
        self.assertEqual(
            medline_citation["KeywordList"][0][3].attributes, {"MajorTopicYN": "N"}
        )
        self.assertEqual(medline_citation["KeywordList"][0][4], "endurance assessment")
        self.assertEqual(
            medline_citation["KeywordList"][0][4].attributes, {"MajorTopicYN": "N"}
        )
        self.assertEqual(medline_citation["KeywordList"][0][5], "lactate threshold")
        self.assertEqual(
            medline_citation["KeywordList"][0][5].attributes, {"MajorTopicYN": "N"}
        )
        self.assertEqual(
            medline_citation["KeywordList"][0][6], "oxygen endurance performance limit"
        )
        self.assertEqual(
            medline_citation["KeywordList"][0][6].attributes, {"MajorTopicYN": "N"}
        )
        self.assertEqual(
            medline_citation["KeywordList"][0][7], "submaximal exercise testing"
        )
        self.assertEqual(
            medline_citation["KeywordList"][0][7].attributes, {"MajorTopicYN": "N"}
        )
        self.assertEqual(medline_citation["CitationSubset"], [])
        self.assertEqual(medline_citation["OtherAbstract"], [])
        self.assertEqual(medline_citation["OtherID"], [])
        self.assertEqual(medline_citation["SpaceFlightMission"], [])
        self.assertEqual(medline_citation["PMID"], "30108519")
        self.assertEqual(medline_citation["PMID"].attributes, {"Version": "1"})
        self.assertEqual(len(medline_citation["DateRevised"]), 3)
        self.assertEqual(medline_citation["DateRevised"]["Year"], "2018")
        self.assertEqual(medline_citation["DateRevised"]["Month"], "08")
        self.assertEqual(medline_citation["DateRevised"]["Day"], "17")
        self.assertEqual(medline_citation["DateRevised"].attributes, {})
        self.assertEqual(len(medline_citation["MedlineJournalInfo"]), 4)
        self.assertEqual(
            medline_citation["MedlineJournalInfo"]["Country"], "Switzerland"
        )
        self.assertEqual(
            medline_citation["MedlineJournalInfo"]["MedlineTA"], "Front Physiol"
        )
        self.assertEqual(
            medline_citation["MedlineJournalInfo"]["NlmUniqueID"], "101549006"
        )
        self.assertEqual(
            medline_citation["MedlineJournalInfo"]["ISSNLinking"], "1664-042X"
        )
        self.assertEqual(medline_citation["MedlineJournalInfo"].attributes, {})
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"]), 53)
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][0]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][0]["RefSource"],
            "Stat Med. 2008 Feb 28;27(5):778-80",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][0]["PMID"], "17907247"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][0]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][0].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][1]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][1]["RefSource"],
            "Int J Sports Med. 2009 Jan;30(1):40-45",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][1]["PMID"], "19202577"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][1]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][1].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][2]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][2]["RefSource"],
            "Med Sci Sports Exerc. 1995 Jun;27(6):863-7",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][2]["PMID"], "7658947"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][2]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][2].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][3]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][3]["RefSource"],
            "Eur J Appl Physiol. 2010 Apr;108(6):1153-67",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][3]["PMID"], "20033207"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][3]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][3].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][4]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][4]["RefSource"],
            "Med Sci Sports Exerc. 1999 Apr;31(4):578-82",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][4]["PMID"], "10211855"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][4]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][4].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][5]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][5]["RefSource"],
            "Br J Sports Med. 1988 Jun;22(2):51-4",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][5]["PMID"], "3167501"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][5]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][5].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][6]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][6]["RefSource"],
            "Front Physiol. 2017 Jun 08;8:389",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][6]["PMID"], "28642717"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][6]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][6].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][7]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][7]["RefSource"],
            "Med Sci Sports Exerc. 1999 Sep;31(9):1342-5",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][7]["PMID"], "10487378"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][7]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][7].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][8]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][8]["RefSource"],
            "Med Sci Sports Exerc. 1998 Aug;30(8):1304-13",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][8]["PMID"], "9710874"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][8]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][8].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][9]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][9]["RefSource"],
            "Med Sci Sports. 1979 Winter;11(4):338-44",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][9]["PMID"], "530025"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][9]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][9].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][10]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][10]["RefSource"],
            "J Strength Cond Res. 2005 May;19(2):364-8",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][10]["PMID"], "15903376"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][10]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][10].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][11]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][11]["RefSource"],
            "Eur J Appl Physiol Occup Physiol. 1984;53(3):196-9",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][11]["PMID"], "6542852"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][11]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][11].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][12]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][12]["RefSource"],
            "Eur J Appl Physiol Occup Physiol. 1978 Oct 20;39(4):219-27",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][12]["PMID"], "710387"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][12]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][12].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][13]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][13]["RefSource"],
            "J Appl Physiol Respir Environ Exerc Physiol. 1980 Mar;48(3):523-7",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][13]["PMID"], "7372524"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][13]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][13].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][14]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][14]["RefSource"],
            "Int J Sports Med. 2015 Dec;36(14):1142-8",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][14]["PMID"], "26332904"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][14]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][14].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][15]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][15]["RefSource"],
            "J Physiol. 1930 Apr 14;69(2):214-37",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][15]["PMID"], "16994099"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][15]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][15].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][16]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][16]["RefSource"],
            "J Strength Cond Res. 2015 Oct;29(10):2794-801",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][16]["PMID"], "25844867"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][16]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][16].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][17]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][17]["RefSource"],
            "PLoS One. 2018 Mar 13;13(3):e0194313",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][17]["PMID"], "29534108"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][17]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][17].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][18]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][18]["RefSource"],
            "J Cardiopulm Rehabil Prev. 2012 Nov-Dec;32(6):327-50",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][18]["PMID"], "23103476"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][18]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][18].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][19]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][19]["RefSource"],
            "Exerc Sport Sci Rev. 1982;10:49-83",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][19]["PMID"], "6811284"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][19]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][19].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][20]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][20]["RefSource"],
            "Int J Sports Physiol Perform. 2010 Sep;5(3):276-91",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][20]["PMID"], "20861519"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][20]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][20].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][21]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][21]["RefSource"],
            "Eur J Appl Physiol Occup Physiol. 1990;60(4):249-53",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][21]["PMID"], "2357979"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][21]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][21].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][22]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][22]["RefSource"],
            "Med Sci Sports Exerc. 2004 Oct;36(10):1737-42",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][22]["PMID"], "15595295"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][22]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][22].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][23]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][23]["RefSource"],
            "Int J Sports Med. 2016 Jun;37(7):539-46",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][23]["PMID"], "27116348"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][23]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][23].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][24]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][24]["RefSource"],
            "Scand J Med Sci Sports. 2017 May;27(5):462-473",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][24]["PMID"], "28181710"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][24]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][24].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][25]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][25]["RefSource"],
            "Int J Sports Med. 1983 Nov;4(4):226-30",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][25]["PMID"], "6654546"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][25]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][25].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][26]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][26]["RefSource"],
            "J Appl Physiol (1985). 1988 Jun;64(6):2622-30",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][26]["PMID"], "3403447"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][26]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][26].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][27]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][27]["RefSource"],
            "Med Sci Sports Exerc. 2009 Jan;41(1):3-13",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][27]["PMID"], "19092709"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][27]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][27].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][28]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][28]["RefSource"],
            "Int J Sports Med. 2009 Sep;30(9):643-6",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][28]["PMID"], "19569005"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][28]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][28].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][29]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][29]["RefSource"],
            "Eur J Appl Physiol Occup Physiol. 1988;57(4):420-4",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][29]["PMID"], "3396556"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][29]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][29].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][30]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][30]["RefSource"],
            "J Physiol. 2004 Jul 1;558(Pt 1):5-30",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][30]["PMID"], "15131240"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][30]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][30].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][31]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][31]["RefSource"],
            "Int J Sports Med. 1990 Feb;11(1):26-32",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][31]["PMID"], "2318561"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][31]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][31].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][32]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][32]["RefSource"],
            "J Appl Physiol. 1973 Aug;35(2):236-43",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][32]["PMID"], "4723033"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][32]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][32].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][33]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][33]["RefSource"],
            "Int J Sports Med. 1987 Dec;8(6):401-6",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][33]["PMID"], "3429086"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][33]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][33].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][34]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][34]["RefSource"],
            "J Sci Med Sport. 2008 Jun;11(3):280-6",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][34]["PMID"], "17553745"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][34]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][34].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][35]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][35]["RefSource"],
            "J Appl Physiol Respir Environ Exerc Physiol. 1984 May;56(5):1260-4",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][35]["PMID"], "6725086"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][35]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][35].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][36]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][36]["RefSource"],
            "Int J Sports Med. 2008 Jun;29(6):475-9",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][36]["PMID"], "18302077"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][36]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][36].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][37]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][37]["RefSource"],
            "Med Sci Sports Exerc. 1985 Feb;17(1):22-34",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][37]["PMID"], "3884959"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][37]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][37].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][38]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][38]["RefSource"],
            "Sports Med. 2009;39(6):469-90",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][38]["PMID"], "19453206"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][38]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][38].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][39]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][39]["RefSource"],
            "Int J Sports Med. 2004 Aug;25(6):403-8",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][39]["PMID"], "15346226"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][39]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][39].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][40]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][40]["RefSource"],
            "J Sports Med Phys Fitness. 2004 Jun;44(2):132-40",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][40]["PMID"], "15470310"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][40]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][40].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][41]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][41]["RefSource"],
            "Int J Sports Med. 1985 Jun;6(3):117-30",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][41]["PMID"], "4030186"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][41]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][41].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][42]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][42]["RefSource"],
            "Int J Sports Med. 1999 Feb;20(2):122-7",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][42]["PMID"], "10190774"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][42]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][42].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][43]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][43]["RefSource"],
            "Int J Sports Med. 2006 May;27(5):368-72",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][43]["PMID"], "16729378"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][43]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][43].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][44]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][44]["RefSource"],
            "Int J Sports Med. 1985 Jun;6(3):109-16",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][44]["PMID"], "3897079"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][44]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][44].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][45]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][45]["RefSource"],
            "Pneumologie. 1990 Jan;44(1):2-13",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][45]["PMID"], "2408033"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][45]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][45].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][46]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][46]["RefSource"],
            "Eur J Appl Physiol. 2018 Apr;118(4):691-728",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][46]["PMID"], "29322250"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][46]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][46].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][47]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][47]["RefSource"],
            "J Appl Physiol Respir Environ Exerc Physiol. 1983 Oct;55(4):1178-86",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][47]["PMID"], "6629951"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][47]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][47].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][48]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][48]["RefSource"],
            "Sports Med. 2014 Nov;44 Suppl 2:S139-47",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][48]["PMID"], "25200666"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][48]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][48].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][49]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][49]["RefSource"],
            "Front Physiol. 2015 Oct 30;6:308",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][49]["PMID"], "26578980"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][49]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][49].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][50]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][50]["RefSource"],
            "Int J Sports Med. 2013 Mar;34(3):196-9",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][50]["PMID"], "22972242"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][50]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][50].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][51]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][51]["RefSource"],
            "Int J Sports Med. 1992 Oct;13(7):518-22",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][51]["PMID"], "1459746"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][51]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][51].attributes,
            {"RefType": "Cites"},
        )
        self.assertEqual(len(medline_citation["CommentsCorrectionsList"][52]), 2)
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][52]["RefSource"],
            "Med Sci Sports Exerc. 1993 May;25(5):620-7",
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][52]["PMID"], "8492691"
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][52]["PMID"].attributes,
            {"Version": "1"},
        )
        self.assertEqual(
            medline_citation["CommentsCorrectionsList"][52].attributes,
            {"RefType": "Cites"},
        )
        article = medline_citation["Article"]
        self.assertEqual(len(article["ELocationID"]), 1)
        self.assertEqual(article["ELocationID"][0], "10.3389/fphys.2018.01034")
        self.assertEqual(
            article["ELocationID"][0].attributes, {"EIdType": "doi", "ValidYN": "Y"}
        )
        self.assertEqual(len(article["ArticleDate"]), 1)
        self.assertEqual(len(article["ArticleDate"][0]), 3)
        self.assertEqual(article["ArticleDate"][0]["Year"], "2018")
        self.assertEqual(article["ArticleDate"][0]["Month"], "07")
        self.assertEqual(article["ArticleDate"][0]["Day"], "31")
        self.assertEqual(
            article["ArticleDate"][0].attributes, {"DateType": "Electronic"}
        )
        self.assertEqual(article["Language"], ["eng"])
        self.assertEqual(len(article["Journal"]), 4)
        self.assertEqual(article["Journal"]["ISSN"], "1664-042X")
        self.assertEqual(article["Journal"]["ISSN"].attributes, {"IssnType": "Print"})
        self.assertEqual(article["Journal"]["JournalIssue"]["Volume"], "9")
        self.assertEqual(article["Journal"]["JournalIssue"]["PubDate"]["Year"], "2018")
        self.assertEqual(article["Journal"]["JournalIssue"]["PubDate"].attributes, {})
        self.assertEqual(
            article["Journal"]["JournalIssue"].attributes, {"CitedMedium": "Print"}
        )
        self.assertEqual(article["Journal"]["Title"], "Frontiers in physiology")
        self.assertEqual(article["Journal"]["ISOAbbreviation"], "Front Physiol")
        self.assertEqual(article["Journal"].attributes, {})
        self.assertEqual(len(article["PublicationTypeList"]), 1)
        self.assertEqual(article["PublicationTypeList"][0], "Journal Article")
        self.assertEqual(
            article["PublicationTypeList"][0].attributes, {"UI": "D016428"}
        )
        self.assertEqual(
            article["ArticleTitle"],
            'A "<i>Blood Relationship"</i> Between the Overlooked Minimum Lactate Equivalent and Maximal Lactate Steady State in Trained Runners. Back to the Old Days?',
        )
        self.assertEqual(len(article["Pagination"]), 1)
        self.assertEqual(article["Pagination"]["MedlinePgn"], "1034")
        self.assertEqual(article["Pagination"].attributes, {})
        self.assertEqual(len(article["AuthorList"]), 2)
        self.assertEqual(len(article["AuthorList"][0]), 5)
        self.assertEqual(article["AuthorList"][0]["Identifier"], [])
        self.assertEqual(len(article["AuthorList"][0]["AffiliationInfo"]), 1)
        self.assertEqual(len(article["AuthorList"][0]["AffiliationInfo"][0]), 2)
        self.assertEqual(
            article["AuthorList"][0]["AffiliationInfo"][0]["Identifier"], []
        )
        self.assertEqual(
            article["AuthorList"][0]["AffiliationInfo"][0]["Affiliation"],
            "Studies, Research and Sports Medicine Center, Government of Navarre, Pamplona, Spain.",
        )
        self.assertEqual(article["AuthorList"][0]["AffiliationInfo"][0].attributes, {})
        self.assertEqual(article["AuthorList"][0]["LastName"], "Garcia-Tabar")
        self.assertEqual(article["AuthorList"][0]["ForeName"], "Ibai")
        self.assertEqual(article["AuthorList"][0]["Initials"], "I")
        self.assertEqual(article["AuthorList"][0].attributes, {"ValidYN": "Y"})
        self.assertEqual(len(article["AuthorList"][1]), 5)
        self.assertEqual(article["AuthorList"][1]["Identifier"], [])
        self.assertEqual(len(article["AuthorList"][1]["AffiliationInfo"]), 1)
        self.assertEqual(len(article["AuthorList"][1]["AffiliationInfo"][0]), 2)
        self.assertEqual(
            article["AuthorList"][1]["AffiliationInfo"][0]["Identifier"], []
        )
        self.assertEqual(
            article["AuthorList"][1]["AffiliationInfo"][0]["Affiliation"],
            "Studies, Research and Sports Medicine Center, Government of Navarre, Pamplona, Spain.",
        )
        self.assertEqual(article["AuthorList"][1]["AffiliationInfo"][0].attributes, {})
        self.assertEqual(article["AuthorList"][1]["LastName"], "Gorostiaga")
        self.assertEqual(article["AuthorList"][1]["ForeName"], "Esteban M")
        self.assertEqual(article["AuthorList"][1]["Initials"], "EM")
        self.assertEqual(article["AuthorList"][1].attributes, {"ValidYN": "Y"})
        self.assertEqual(len(article["Abstract"]), 1)
        self.assertEqual(
            article["Abstract"]["AbstractText"][0],
            """\
Maximal Lactate Steady State (MLSS) and Lactate Threshold (LT) are physiologically-related and fundamental concepts within the sports and exercise sciences. Literature supporting their relationship, however, is scarce. Among the recognized LTs, we were particularly interested in the disused "Minimum Lactate Equivalent" (LE<sub>min</sub>), first described in the early 1980s. We hypothesized that velocity at LT, conceptually comprehended as in the old days (LE<sub>min</sub>), could predict velocity at MLSS (<sub>V</sub>MLSS) more accurate than some other blood lactate-related thresholds (BL<sub>R</sub>Ts) routinely used nowadays by many sport science practitioners. Thirteen male endurance-trained [<sub>V</sub>MLSS 15.0 ± 1.1 km·h<sup>-1</sup>; maximal oxygen uptake ( <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <msub>
                            <mrow>
                                <mover>
                                    <mrow>
                                        <mi>V</mi>
                                    </mrow>
                                    <mo>.</mo>
                                </mover>
                                <mi>O</mi>
                            </mrow>
                            <mrow>
                                <mn>2</mn>
                                <mi>m</mi>
                                <mi>a</mi>
                                <mi>x</mi>
                            </mrow>
                        </msub>
                    </math> ) 67.6 ± 4.1 ml·kg<sup>-1</sup>·min<sup>-1</sup>] homogeneous (coefficient of variation: ≈7%) runners conducted 1) a submaximal discontinuous incremental running test to determine several BL<sub>R</sub>Ts followed by a maximal ramp incremental running test for <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <msub>
                            <mrow>
                                <mover>
                                    <mrow>
                                        <mi>V</mi>
                                    </mrow>
                                    <mo>.</mo>
                                </mover>
                                <mi>O</mi>
                            </mrow>
                            <mrow>
                                <mn>2</mn>
                                <mi>m</mi>
                                <mi>a</mi>
                                <mi>x</mi>
                            </mrow>
                        </msub>
                        <mtext> </mtext>
                    </math> determination, and 2) several (4-5) constant velocity running tests to determine <sub>V</sub>MLSS with a precision of 0.20 km·h<sup>-1</sup>. Determined BL<sub>R</sub>Ts include LE<sub>min</sub> and LE<sub>min</sub>-related LE<sub>min</sub> plus 1 (LE<sub>min+1mM</sub>) and 1.5 mmol·L<sup>-1</sup> (LE<sub>min+1.5mM</sub>), along with well-established BL<sub>R</sub>Ts such as conventionally-calculated LT, D<sub>max</sub> and fixed blood lactate concentration thresholds. LE<sub>min</sub> did not differ from LT (<i>P</i> = 0.71; ES: 0.08) and was 27% lower than MLSS (<i>P</i> < 0.001; ES: 3.54). LE<sub>min+1mM</sub> was not different from MLSS (<i>P</i> = 0.47; ES: 0.09). LE<sub>min</sub> was the best predictor of <sub>V</sub>MLSS (<i>r</i> = 0.91; <i>P</i> < 0.001; SEE = 0.47 km·h<sup>-1</sup>), followed by LE<sub>min+1mM</sub> (<i>r</i> = 0.86; <i>P</i> < 0.001; SEE = 0.58 km·h<sup>-1</sup>) and LE<sub>min+1.5mM</sub> (<i>r</i> = 0.84; <i>P</i> < 0.001; SEE = 0.86 km·h<sup>-1</sup>). There was no statistical difference between MLSS and estimated MLSS using LE<sub>min</sub> prediction formula (<i>P</i> = 0.99; ES: 0.001). Mean bias and limits of agreement were 0.00 ± 0.45 km·h<sup>-1</sup> and ±0.89 km·h<sup>-1</sup>. Additionally, LE<sub>min</sub>, LE<sub>min+1mM</sub> and LE<sub>min+1.5mM</sub> were the best predictors of <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <msub>
                            <mrow>
                                <mover>
                                    <mrow>
                                        <mi>V</mi>
                                    </mrow>
                                    <mo>.</mo>
                                </mover>
                                <mi>O</mi>
                            </mrow>
                            <mrow>
                                <mn>2</mn>
                                <mi>m</mi>
                                <mi>a</mi>
                                <mi>x</mi>
                            </mrow>
                        </msub>
                    </math> (<i>r</i> = 0.72-0.79; <i>P</i> < 0.001). These results support LE<sub>min</sub>, an objective submaximal overlooked and underused BL<sub>R</sub>T, to be one of the best single MLSS predictors in endurance trained runners. Our study advocates factors controlling LE<sub>min</sub> to be shared, at least partly, with those controlling MLSS.""",
        )

    def test_pubmed_mathml_tags(self):
        """Test parsing XML returned by EFetch, PubMed database, with extensive MathML tags."""
        # In PubMed display PMID 29963580 in xml retrieval mode, containing
        # extensive MathML tags in the abstract text.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="pubmed", id="29963580", retmode="xml")
        with open("Entrez/pubmed7.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(len(record), 2)
        self.assertEqual(record["PubmedBookArticle"], [])
        self.assertEqual(len(record["PubmedArticle"]), 1)
        pubmed_article = record["PubmedArticle"][0]
        self.assertEqual(len(pubmed_article), 2)
        self.assertEqual(len(pubmed_article["MedlineCitation"].attributes), 2)
        self.assertEqual(
            pubmed_article["MedlineCitation"].attributes["Status"], "PubMed-not-MEDLINE"
        )
        self.assertEqual(pubmed_article["MedlineCitation"].attributes["Owner"], "NLM")
        self.assertEqual(pubmed_article["MedlineCitation"]["PMID"], "29963580")
        self.assertEqual(
            pubmed_article["MedlineCitation"]["PMID"].attributes["Version"], "1"
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["DateRevised"].attributes, {}
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["DateRevised"]["Year"], "2018"
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["DateRevised"]["Month"], "11"
        )
        self.assertEqual(pubmed_article["MedlineCitation"]["DateRevised"]["Day"], "14")
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["Article"].attributes), 1
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"].attributes["PubModel"],
            "Print-Electronic",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"].attributes, {}
        )
        self.assertEqual(
            len(
                pubmed_article["MedlineCitation"]["Article"]["Journal"][
                    "ISSN"
                ].attributes
            ),
            1,
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes[
                "IssnType"
            ],
            "Print",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"]["ISSN"], "2329-4302"
        )
        self.assertEqual(
            len(
                pubmed_article["MedlineCitation"]["Article"]["Journal"][
                    "JournalIssue"
                ].attributes
            ),
            1,
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"][
                "JournalIssue"
            ].attributes["CitedMedium"],
            "Print",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Volume"
            ],
            "5",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "Issue"
            ],
            "2",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Year"],
            "2018",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"][
                "PubDate"
            ]["Month"],
            "Apr",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"]["Title"],
            "Journal of medical imaging (Bellingham, Wash.)",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"],
            "J Med Imaging (Bellingham)",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["ArticleTitle"],
            "Development of a pulmonary imaging biomarker pipeline for phenotyping of chronic lung disease.",
        )
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["Article"]["Pagination"]), 1
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"],
            "026002",
        )
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["Article"]["ELocationID"]), 1
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["ELocationID"][0].attributes[
                "EIdType"
            ],
            "doi",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["ELocationID"][0].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["ELocationID"][0],
            "10.1117/1.JMI.5.2.026002",
        )
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["Article"]["Abstract"]), 1
        )
        self.assertEqual(
            len(
                pubmed_article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
            ),
            1,
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0],
            """\
We designed and generated pulmonary imaging biomarker pipelines to facilitate high-throughput research and point-of-care use in patients with chronic lung disease. Image processing modules and algorithm pipelines were embedded within a graphical user interface (based on the .NET framework) for pulmonary magnetic resonance imaging (MRI) and x-ray computed-tomography (CT) datasets. The software pipelines were generated using C++ and included: (1) inhaled <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <mrow>
                            <mmultiscripts>
                                <mrow>
                                    <mi>He</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>3</mn>
                                </mrow>
                            </mmultiscripts>
                            <mo>/</mo>
                            <mmultiscripts>
                                <mrow>
                                    <mi>Xe</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>129</mn>
                                </mrow>
                            </mmultiscripts>
                            <mtext> </mtext>
                            <mi>MRI</mi>
                        </mrow>
                    </math> ventilation and apparent diffusion coefficients, (2) CT-MRI coregistration for lobar and segmental ventilation and perfusion measurements, (3) ultrashort echo-time <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <mrow>
                            <mmultiscripts>
                                <mrow>
                                    <mi>H</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>1</mn>
                                </mrow>
                            </mmultiscripts>
                            <mtext> </mtext>
                            <mi>MRI</mi>
                        </mrow>
                    </math> proton density measurements, (4) free-breathing Fourier-decomposition <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <mrow>
                            <mmultiscripts>
                                <mrow>
                                    <mi>H</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>1</mn>
                                </mrow>
                            </mmultiscripts>
                            <mtext> </mtext>
                            <mi>MRI</mi>
                        </mrow>
                    </math> ventilation/perfusion and free-breathing <math xmlns="http://www.w3.org/1998/Math/MathML">
                        <mrow>
                            <mmultiscripts>
                                <mrow>
                                    <mi>H</mi>
                                </mrow>
                                <mprescripts></mprescripts>
                                <none></none>
                                <mrow>
                                    <mn>1</mn>
                                </mrow>
                            </mmultiscripts>
                            <mtext> </mtext>
                            <mi>MRI</mi>
                        </mrow>
                    </math> specific ventilation, (5)\u00a0multivolume CT and MRI parametric response maps, and (6)\u00a0MRI and CT texture analysis and radiomics. The image analysis framework was implemented on a desktop workstation/tablet to generate biomarkers of regional lung structure and function related to ventilation, perfusion, lung tissue texture, and integrity as well as multiparametric measures of gas trapping and airspace enlargement. All biomarkers were generated within 10 min with measurement reproducibility consistent with clinical and research requirements. The resultant pulmonary imaging biomarker pipeline provides real-time and automated lung imaging measurements for point-of-care and high-throughput research.""",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"].attributes[
                "CompleteYN"
            ],
            "Y",
        )
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["Article"]["AuthorList"]), 9
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"],
            "Guo",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"],
            "Fumin",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"],
            "F",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0][
                "AffiliationInfo"
            ][0]["Affiliation"],
            "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0][
                "AffiliationInfo"
            ][1]["Affiliation"],
            "University of Western Ontario, Graduate Program in Biomedical Engineering, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][0][
                "AffiliationInfo"
            ][2]["Affiliation"],
            "University of Toronto, Sunnybrook Research Institute, Toronto, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"],
            "Capaldi",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"],
            "Dante",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"],
            "D",
        )
        self.assertEqual(
            len(
                pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1][
                    "Identifier"
                ]
            ),
            1,
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["Identifier"][
                0
            ].attributes["Source"],
            "ORCID",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1]["Identifier"][
                0
            ],
            "https://orcid.org/0000-0002-4590-7461",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1][
                "AffiliationInfo"
            ][0]["Affiliation"],
            "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][1][
                "AffiliationInfo"
            ][1]["Affiliation"],
            "University of Western Ontario, Department of Medical Biophysics, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"],
            "Kirby",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"],
            "Miranda",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"],
            "M",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][2][
                "AffiliationInfo"
            ][0]["Affiliation"],
            "University of British Columbia, St. Paul's Hospital, Centre for Heart Lung Innovation, Vancouver, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"],
            "Sheikh",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"],
            "Khadija",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"],
            "K",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][3][
                "AffiliationInfo"
            ][0]["Affiliation"],
            "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"],
            "Svenningsen",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"],
            "Sarah",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"],
            "S",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][4][
                "AffiliationInfo"
            ][0]["Affiliation"],
            "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"],
            "McCormack",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"],
            "David G",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"],
            "DG",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][5][
                "AffiliationInfo"
            ][0]["Affiliation"],
            "University of Western Ontario, Division of Respirology, Department of Medicine, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["LastName"],
            "Fenster",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["ForeName"],
            "Aaron",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["Initials"],
            "A",
        )
        self.assertEqual(
            len(
                pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6][
                    "Identifier"
                ]
            ),
            1,
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["Identifier"][
                0
            ].attributes["Source"],
            "ORCID",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6]["Identifier"][
                0
            ],
            "https://orcid.org/0000-0003-3525-2788",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6][
                "AffiliationInfo"
            ][0]["Affiliation"],
            "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6][
                "AffiliationInfo"
            ][1]["Affiliation"],
            "University of Western Ontario, Graduate Program in Biomedical Engineering, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][6][
                "AffiliationInfo"
            ][2]["Affiliation"],
            "University of Western Ontario, Department of Medical Biophysics, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7]["LastName"],
            "Parraga",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7]["ForeName"],
            "Grace",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7]["Initials"],
            "G",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7][
                "AffiliationInfo"
            ][0]["Affiliation"],
            "University of Western Ontario, Robarts Research Institute, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7][
                "AffiliationInfo"
            ][1]["Affiliation"],
            "University of Western Ontario, Graduate Program in Biomedical Engineering, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][7][
                "AffiliationInfo"
            ][2]["Affiliation"],
            "University of Western Ontario, Department of Medical Biophysics, London, Ontario, Canada.",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][8].attributes[
                "ValidYN"
            ],
            "Y",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["AuthorList"][8][
                "CollectiveName"
            ],
            "Canadian Respiratory Research Network",
        )
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["Article"]["Language"]), 1
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["Language"][0], "eng"
        )
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["Article"]["PublicationTypeList"]), 1
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["PublicationTypeList"][
                0
            ].attributes["UI"],
            "D016428",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["PublicationTypeList"][0],
            "Journal Article",
        )
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["Article"]["ArticleDate"]), 1
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["ArticleDate"][0].attributes[
                "DateType"
            ],
            "Electronic",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"],
            "2018",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["ArticleDate"][0]["Month"],
            "06",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["Article"]["ArticleDate"][0]["Day"], "28"
        )
        self.assertEqual(
            len(pubmed_article["MedlineCitation"]["MedlineJournalInfo"]), 4
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["MedlineJournalInfo"]["Country"],
            "United States",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"],
            "J Med Imaging (Bellingham)",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"],
            "101643461",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["MedlineJournalInfo"]["ISSNLinking"],
            "2329-4302",
        )
        self.assertEqual(len(pubmed_article["MedlineCitation"]["KeywordList"]), 1)
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0].attributes["Owner"],
            "NOTNLM",
        )
        self.assertEqual(len(pubmed_article["MedlineCitation"]["KeywordList"][0]), 5)
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][0].attributes[
                "MajorTopicYN"
            ],
            "N",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][0], "asthma"
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][1].attributes[
                "MajorTopicYN"
            ],
            "N",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][1],
            "chronic obstructive lung disease",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][2].attributes[
                "MajorTopicYN"
            ],
            "N",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][2],
            "image processing, biomarkers",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][3].attributes[
                "MajorTopicYN"
            ],
            "N",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][3],
            "magnetic resonance imaging",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][4].attributes[
                "MajorTopicYN"
            ],
            "N",
        )
        self.assertEqual(
            pubmed_article["MedlineCitation"]["KeywordList"][0][4],
            "thoracic computed tomography",
        )
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][0].attributes["PubStatus"],
            "received",
        )
        self.assertEqual(pubmed_article["PubmedData"]["History"][0]["Year"], "2017")
        self.assertEqual(pubmed_article["PubmedData"]["History"][0]["Month"], "12")
        self.assertEqual(pubmed_article["PubmedData"]["History"][0]["Day"], "12")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][1].attributes["PubStatus"],
            "accepted",
        )
        self.assertEqual(pubmed_article["PubmedData"]["History"][1]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][1]["Month"], "06")
        self.assertEqual(pubmed_article["PubmedData"]["History"][1]["Day"], "14")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][2].attributes["PubStatus"],
            "pmc-release",
        )
        self.assertEqual(pubmed_article["PubmedData"]["History"][2]["Year"], "2019")
        self.assertEqual(pubmed_article["PubmedData"]["History"][2]["Month"], "06")
        self.assertEqual(pubmed_article["PubmedData"]["History"][2]["Day"], "28")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][3].attributes["PubStatus"], "entrez"
        )
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Month"], "7")
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Day"], "3")
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Hour"], "6")
        self.assertEqual(pubmed_article["PubmedData"]["History"][3]["Minute"], "0")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][4].attributes["PubStatus"], "pubmed"
        )
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Month"], "7")
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Day"], "3")
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Hour"], "6")
        self.assertEqual(pubmed_article["PubmedData"]["History"][4]["Minute"], "0")
        self.assertEqual(
            pubmed_article["PubmedData"]["History"][5].attributes["PubStatus"],
            "medline",
        )
        self.assertEqual(pubmed_article["PubmedData"]["History"][5]["Year"], "2018")
        self.assertEqual(pubmed_article["PubmedData"]["History"][5]["Month"], "7")
        self.assertEqual(pubmed_article["PubmedData"]["History"][5]["Day"], "3")
        self.assertEqual(pubmed_article["PubmedData"]["History"][5]["Hour"], "6")
        self.assertEqual(pubmed_article["PubmedData"]["History"][5]["Minute"], "1")
        self.assertEqual(pubmed_article["PubmedData"]["PublicationStatus"], "ppublish")
        self.assertEqual(len(pubmed_article["PubmedData"]["ArticleIdList"]), 4)
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][0].attributes["IdType"],
            "pubmed",
        )
        self.assertEqual(pubmed_article["PubmedData"]["ArticleIdList"][0], "29963580")
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][1].attributes["IdType"], "doi"
        )
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][1], "10.1117/1.JMI.5.2.026002"
        )
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][2].attributes["IdType"], "pii"
        )
        self.assertEqual(pubmed_article["PubmedData"]["ArticleIdList"][2], "17360RR")
        self.assertEqual(
            pubmed_article["PubmedData"]["ArticleIdList"][3].attributes["IdType"], "pmc"
        )
        self.assertEqual(pubmed_article["PubmedData"]["ArticleIdList"][3], "PMC6022861")
        self.assertEqual(len(pubmed_article["PubmedData"]["ReferenceList"]), 1)
        self.assertEqual(len(pubmed_article["PubmedData"]["ReferenceList"][0]), 2)
        self.assertEqual(
            len(pubmed_article["PubmedData"]["ReferenceList"][0]["ReferenceList"]), 0
        )
        references = pubmed_article["PubmedData"]["ReferenceList"][0]["Reference"]
        self.assertEqual(len(references), 49)
        self.assertEqual(references[0]["Citation"], "Radiology. 2015 Jan;274(1):250-9")
        self.assertEqual(len(references[0]["ArticleIdList"]), 1)
        self.assertEqual(references[0]["ArticleIdList"][0], "25144646")
        self.assertEqual(
            references[0]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[1]["Citation"], "Nature. 1994 Jul 21;370(6486):199-201"
        )
        self.assertEqual(len(references[1]["ArticleIdList"]), 1)
        self.assertEqual(references[1]["ArticleIdList"][0], "8028666")
        self.assertEqual(
            references[1]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[2]["Citation"], "Magn Reson Med. 2009 Sep;62(3):656-64"
        )
        self.assertEqual(len(references[2]["ArticleIdList"]), 1)
        self.assertEqual(references[2]["ArticleIdList"][0], "19585597")
        self.assertEqual(
            references[2]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[3]["Citation"], "Radiology. 2016 Feb;278(2):563-77")
        self.assertEqual(len(references[3]["ArticleIdList"]), 1)
        self.assertEqual(references[3]["ArticleIdList"][0], "26579733")
        self.assertEqual(
            references[3]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[4]["Citation"], "Radiology. 1991 Jun;179(3):777-81")
        self.assertEqual(len(references[4]["ArticleIdList"]), 1)
        self.assertEqual(references[4]["ArticleIdList"][0], "2027991")
        self.assertEqual(
            references[4]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[5]["Citation"], "Radiology. 2010 Jul;256(1):280-9")
        self.assertEqual(len(references[5]["ArticleIdList"]), 1)
        self.assertEqual(references[5]["ArticleIdList"][0], "20574101")
        self.assertEqual(
            references[5]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[6]["Citation"], "Phys Med Biol. 2001 May;46(5):R67-99"
        )
        self.assertEqual(len(references[6]["ArticleIdList"]), 1)
        self.assertEqual(references[6]["ArticleIdList"][0], "11384074")
        self.assertEqual(
            references[6]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[7]["Citation"], "IEEE Trans Med Imaging. 2011 Nov;30(11):1901-20"
        )
        self.assertEqual(len(references[7]["ArticleIdList"]), 1)
        self.assertEqual(references[7]["ArticleIdList"][0], "21632295")
        self.assertEqual(
            references[7]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[8]["Citation"], "Eur J Cancer. 2012 Mar;48(4):441-6"
        )
        self.assertEqual(len(references[8]["ArticleIdList"]), 1)
        self.assertEqual(references[8]["ArticleIdList"][0], "22257792")
        self.assertEqual(
            references[8]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[9]["Citation"], "Med Phys. 2017 May;44(5):1718-1733"
        )
        self.assertEqual(len(references[9]["ArticleIdList"]), 1)
        self.assertEqual(references[9]["ArticleIdList"][0], "28206676")
        self.assertEqual(
            references[9]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[10]["Citation"], "COPD. 2014 Apr;11(2):125-32")
        self.assertEqual(len(references[10]["ArticleIdList"]), 1)
        self.assertEqual(references[10]["ArticleIdList"][0], "22433011")
        self.assertEqual(
            references[10]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[11]["Citation"],
            "Am J Respir Crit Care Med. 2015 Nov 15;192(10):1215-22",
        )
        self.assertEqual(len(references[11]["ArticleIdList"]), 1)
        self.assertEqual(references[11]["ArticleIdList"][0], "26186608")
        self.assertEqual(
            references[11]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[12]["Citation"], "N Engl J Med. 2016 May 12;374(19):1811-21"
        )
        self.assertEqual(len(references[12]["ArticleIdList"]), 1)
        self.assertEqual(references[12]["ArticleIdList"][0], "27168432")
        self.assertEqual(
            references[12]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[13]["Citation"], "Am J Epidemiol. 2002 Nov 1;156(9):871-81"
        )
        self.assertEqual(len(references[13]["ArticleIdList"]), 1)
        self.assertEqual(references[13]["ArticleIdList"][0], "12397006")
        self.assertEqual(
            references[13]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[14]["Citation"], "BMC Cancer. 2014 Dec 11;14:934")
        self.assertEqual(len(references[14]["ArticleIdList"]), 1)
        self.assertEqual(references[14]["ArticleIdList"][0], "25496482")
        self.assertEqual(
            references[14]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[15]["Citation"], "Acad Radiol. 2015 Mar;22(3):320-9"
        )
        self.assertEqual(len(references[15]["ArticleIdList"]), 1)
        self.assertEqual(references[15]["ArticleIdList"][0], "25491735")
        self.assertEqual(
            references[15]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[16]["Citation"], "Chest. 1999 Dec;116(6):1750-61")
        self.assertEqual(len(references[16]["ArticleIdList"]), 1)
        self.assertEqual(references[16]["ArticleIdList"][0], "10593802")
        self.assertEqual(
            references[16]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[17]["Citation"], "Acad Radiol. 2012 Feb;19(2):141-52"
        )
        self.assertEqual(len(references[17]["ArticleIdList"]), 1)
        self.assertEqual(references[17]["ArticleIdList"][0], "22104288")
        self.assertEqual(
            references[17]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[18]["Citation"], "Med Phys. 2014 Mar;41(3):033502")
        self.assertEqual(len(references[18]["ArticleIdList"]), 1)
        self.assertEqual(references[18]["ArticleIdList"][0], "24593744")
        self.assertEqual(
            references[18]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[19]["Citation"], "Med Phys. 2008 Oct;35(10):4695-707"
        )
        self.assertEqual(len(references[19]["ArticleIdList"]), 1)
        self.assertEqual(references[19]["ArticleIdList"][0], "18975715")
        self.assertEqual(
            references[19]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[20]["Citation"], "Thorax. 2017 May;72(5):475-477")
        self.assertEqual(len(references[20]["ArticleIdList"]), 1)
        self.assertEqual(references[20]["ArticleIdList"][0], "28258250")
        self.assertEqual(
            references[20]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[21]["Citation"], "Nat Med. 1996 Nov;2(11):1236-9")
        self.assertEqual(len(references[21]["ArticleIdList"]), 1)
        self.assertEqual(references[21]["ArticleIdList"][0], "8898751")
        self.assertEqual(
            references[21]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[22]["Citation"], "J Magn Reson Imaging. 2015 May;41(5):1465-74"
        )
        self.assertEqual(len(references[22]["ArticleIdList"]), 1)
        self.assertEqual(references[22]["ArticleIdList"][0], "24965907")
        self.assertEqual(
            references[22]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[23]["Citation"], "Acad Radiol. 2008 Jun;15(6):776-85"
        )
        self.assertEqual(len(references[23]["ArticleIdList"]), 1)
        self.assertEqual(references[23]["ArticleIdList"][0], "18486013")
        self.assertEqual(
            references[23]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[24]["Citation"], "Magn Reson Med. 2000 Aug;44(2):174-9"
        )
        self.assertEqual(len(references[24]["ArticleIdList"]), 1)
        self.assertEqual(references[24]["ArticleIdList"][0], "10918314")
        self.assertEqual(
            references[24]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[25]["Citation"],
            "Am J Respir Crit Care Med. 2014 Jul 15;190(2):135-44",
        )
        self.assertEqual(len(references[25]["ArticleIdList"]), 1)
        self.assertEqual(references[25]["ArticleIdList"][0], "24873985")
        self.assertEqual(
            references[25]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[26]["Citation"], "Med Phys. 2016 Jun;43(6):2911-2926"
        )
        self.assertEqual(len(references[26]["ArticleIdList"]), 1)
        self.assertEqual(references[26]["ArticleIdList"][0], "27277040")
        self.assertEqual(
            references[26]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[27]["Citation"], "Nat Med. 2009 May;15(5):572-6")
        self.assertEqual(len(references[27]["ArticleIdList"]), 1)
        self.assertEqual(references[27]["ArticleIdList"][0], "19377487")
        self.assertEqual(
            references[27]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[28]["Citation"], "Eur J Radiol. 2014 Nov;83(11):2093-101"
        )
        self.assertEqual(len(references[28]["ArticleIdList"]), 1)
        self.assertEqual(references[28]["ArticleIdList"][0], "25176287")
        self.assertEqual(
            references[28]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[29]["Citation"], "Radiology. 2004 Sep;232(3):739-48"
        )
        self.assertEqual(len(references[29]["ArticleIdList"]), 1)
        self.assertEqual(references[29]["ArticleIdList"][0], "15333795")
        self.assertEqual(
            references[29]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[30]["Citation"], "Med Image Anal. 2015 Jul;23(1):43-55"
        )
        self.assertEqual(len(references[30]["ArticleIdList"]), 1)
        self.assertEqual(references[30]["ArticleIdList"][0], "25958028")
        self.assertEqual(
            references[30]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[31]["Citation"], "Radiology. 2015 Oct;277(1):192-205"
        )
        self.assertEqual(len(references[31]["ArticleIdList"]), 1)
        self.assertEqual(references[31]["ArticleIdList"][0], "25961632")
        self.assertEqual(
            references[31]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[32]["Citation"], "Med Image Anal. 2012 Oct;16(7):1423-35"
        )
        self.assertEqual(len(references[32]["ArticleIdList"]), 1)
        self.assertEqual(references[32]["ArticleIdList"][0], "22722056")
        self.assertEqual(
            references[32]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[33]["Citation"], "Radiology. 2016 May;279(2):597-608"
        )
        self.assertEqual(len(references[33]["ArticleIdList"]), 1)
        self.assertEqual(references[33]["ArticleIdList"][0], "26744928")
        self.assertEqual(
            references[33]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[34]["Citation"],
            "J Allergy Clin Immunol. 2003 Jun;111(6):1205-11",
        )
        self.assertEqual(len(references[34]["ArticleIdList"]), 1)
        self.assertEqual(references[34]["ArticleIdList"][0], "12789218")
        self.assertEqual(
            references[34]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[35]["Citation"], "J Magn Reson Imaging. 2016 Mar;43(3):544-57"
        )
        self.assertEqual(len(references[35]["ArticleIdList"]), 1)
        self.assertEqual(references[35]["ArticleIdList"][0], "26199216")
        self.assertEqual(
            references[35]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[36]["Citation"],
            "Am J Respir Crit Care Med. 2016 Oct 1;194(7):794-806",
        )
        self.assertEqual(len(references[36]["ArticleIdList"]), 1)
        self.assertEqual(references[36]["ArticleIdList"][0], "27482984")
        self.assertEqual(
            references[36]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[37]["Citation"], "Radiology. 1996 Nov;201(2):564-8")
        self.assertEqual(len(references[37]["ArticleIdList"]), 1)
        self.assertEqual(references[37]["ArticleIdList"][0], "8888259")
        self.assertEqual(
            references[37]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[38]["Citation"], "Thorax. 2014 May;69(5):491-4")
        self.assertEqual(len(references[38]["ArticleIdList"]), 1)
        self.assertEqual(references[38]["ArticleIdList"][0], "24029743")
        self.assertEqual(
            references[38]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[39]["Citation"], "J Magn Reson Imaging. 2017 Apr;45(4):1204-1215"
        )
        self.assertEqual(len(references[39]["ArticleIdList"]), 1)
        self.assertEqual(references[39]["ArticleIdList"][0], "27731948")
        self.assertEqual(
            references[39]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[40]["Citation"], "J Appl Physiol (1985). 2009 Oct;107(4):1258-65"
        )
        self.assertEqual(len(references[40]["ArticleIdList"]), 1)
        self.assertEqual(references[40]["ArticleIdList"][0], "19661452")
        self.assertEqual(
            references[40]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[41]["Citation"], "Acad Radiol. 2016 Feb;23(2):176-85"
        )
        self.assertEqual(len(references[41]["ArticleIdList"]), 1)
        self.assertEqual(references[41]["ArticleIdList"][0], "26601971")
        self.assertEqual(
            references[41]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[42]["Citation"], "Radiology. 2018 May;287(2):693-704"
        )
        self.assertEqual(len(references[42]["ArticleIdList"]), 1)
        self.assertEqual(references[42]["ArticleIdList"][0], "29470939")
        self.assertEqual(
            references[42]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[43]["Citation"], "Eur Respir J. 2016 Aug;48(2):370-9"
        )
        self.assertEqual(len(references[43]["ArticleIdList"]), 1)
        self.assertEqual(references[43]["ArticleIdList"][0], "27174885")
        self.assertEqual(
            references[43]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[44]["Citation"], "Radiology. 2011 Oct;261(1):283-92"
        )
        self.assertEqual(len(references[44]["ArticleIdList"]), 1)
        self.assertEqual(references[44]["ArticleIdList"][0], "21813741")
        self.assertEqual(
            references[44]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[45]["Citation"],
            "Am J Respir Crit Care Med. 2014 Mar 15;189(6):650-7",
        )
        self.assertEqual(len(references[45]["ArticleIdList"]), 1)
        self.assertEqual(references[45]["ArticleIdList"][0], "24401150")
        self.assertEqual(
            references[45]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[46]["Citation"],
            "Am J Respir Crit Care Med. 2012 Feb 15;185(4):356-62",
        )
        self.assertEqual(len(references[46]["ArticleIdList"]), 1)
        self.assertEqual(references[46]["ArticleIdList"][0], "22095547")
        self.assertEqual(
            references[46]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(references[47]["Citation"], "COPD. 2010 Feb;7(1):32-43")
        self.assertEqual(len(references[47]["ArticleIdList"]), 1)
        self.assertEqual(references[47]["ArticleIdList"][0], "20214461")
        self.assertEqual(
            references[47]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )
        self.assertEqual(
            references[48]["Citation"], "Eur Respir J. 2008 Apr;31(4):869-73"
        )
        self.assertEqual(len(references[48]["ArticleIdList"]), 1)
        self.assertEqual(references[48]["ArticleIdList"][0], "18216052")
        self.assertEqual(
            references[48]["ArticleIdList"][0].attributes["IdType"], "pubmed"
        )

    def test_pmc(self):
        """Test parsing XML returned by EFetch from PubMed Central."""
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pmc', id="2682512,3468381")
        with open("Entrez/efetch_pmc.xml", "rb") as stream:
            records = Entrez.parse(stream)
            records = list(records)
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(len(record), 2)
        self.assertEqual(len(record["front"]), 9)
        self.assertEqual(len(record["front"]["journal-meta"]), 10)
        self.assertEqual(len(record["front"]["journal-meta"]["journal-id"]), 4)
        self.assertEqual(
            record["front"]["journal-meta"]["journal-id"][0], "ERJ Open Res"
        )
        self.assertEqual(
            record["front"]["journal-meta"]["journal-id"][0].attributes,
            {"journal-id-type": "nlm-ta"},
        )
        self.assertEqual(
            record["front"]["journal-meta"]["journal-id"][1], "ERJ Open Res"
        )
        self.assertEqual(
            record["front"]["journal-meta"]["journal-id"][1].attributes,
            {"journal-id-type": "iso-abbrev"},
        )
        self.assertEqual(record["front"]["journal-meta"]["journal-id"][2], "ERJOR")
        self.assertEqual(
            record["front"]["journal-meta"]["journal-id"][2].attributes,
            {"journal-id-type": "publisher-id"},
        )
        self.assertEqual(record["front"]["journal-meta"]["journal-id"][3], "erjor")
        self.assertEqual(
            record["front"]["journal-meta"]["journal-id"][3].attributes,
            {"journal-id-type": "hwp"},
        )
        self.assertEqual(len(record["front"]["journal-meta"]["journal-title-group"]), 1)
        journal_title_group = record["front"]["journal-meta"]["journal-title-group"][0]
        self.assertEqual(len(journal_title_group), 4)
        self.assertEqual(journal_title_group["journal-title"], ["ERJ Open Research"])
        self.assertEqual(journal_title_group["journal-subtitle"], [])
        self.assertEqual(journal_title_group["abbrev-journal-title"], [])
        self.assertEqual(journal_title_group["trans-title-group"], [])
        self.assertEqual(len(record["front"]["journal-meta"]["issn"]), 1)
        self.assertEqual(record["front"]["journal-meta"]["issn"][0], "2312-0541")
        self.assertEqual(
            record["front"]["journal-meta"]["issn"][0].attributes, {"pub-type": "epub"}
        )
        self.assertEqual(len(record["front"]["journal-meta"]["publisher"]), 1)
        self.assertEqual(len(record["front"]["journal-meta"]["publisher"][0]), 1)
        self.assertEqual(
            record["front"]["journal-meta"]["publisher"][0][0],
            "European Respiratory Society",
        )
        self.assertEqual(
            record["front"]["journal-meta"]["publisher"][0][0].tag, "publisher-name"
        )
        self.assertEqual(record["front"]["journal-meta"]["contrib-group"], [])
        self.assertEqual(record["front"]["journal-meta"]["notes"], [])
        self.assertEqual(record["front"]["journal-meta"]["aff"], [])
        self.assertEqual(record["front"]["journal-meta"]["aff-alternatives"], [])
        self.assertEqual(record["front"]["journal-meta"]["self-uri"], [])
        self.assertEqual(record["front"]["journal-meta"]["isbn"], [])
        self.assertEqual(len(record["front"]["article-meta"]), 34)
        self.assertEqual(record["front"]["article-meta"]["abstract"], [])
        self.assertEqual(record["front"]["article-meta"]["funding-group"], [])
        self.assertEqual(record["front"]["article-meta"]["aff"], [])
        self.assertEqual(record["front"]["article-meta"]["issue-title"], [])
        self.assertEqual(len(record["front"]["article-meta"]["pub-date"]), 3)
        self.assertEqual(record["front"]["article-meta"]["pub-date"][0], ["7", "2021"])
        self.assertEqual(
            record["front"]["article-meta"]["pub-date"][0].attributes,
            {"pub-type": "collection"},
        )
        self.assertEqual(
            record["front"]["article-meta"]["pub-date"][1], ["13", "9", "2021"]
        )
        self.assertEqual(
            record["front"]["article-meta"]["pub-date"][1].attributes,
            {"pub-type": "epub"},
        )
        self.assertEqual(
            record["front"]["article-meta"]["pub-date"][2], ["13", "9", "2021"]
        )
        self.assertEqual(
            record["front"]["article-meta"]["pub-date"][2].attributes,
            {"pub-type": "pmc-release"},
        )
        self.assertEqual(record["front"]["article-meta"]["conference"], [])
        self.assertEqual(record["front"]["article-meta"]["supplementary-material"], [])
        self.assertEqual(len(record["front"]["article-meta"]["related-article"]), 1)
        self.assertEqual(record["front"]["article-meta"]["related-article"][0], "")
        self.assertEqual(
            record["front"]["article-meta"]["related-article"][0].attributes,
            {
                "related-article-type": "corrected-article",
                "id": "d31e52",
                "ext-link-type": "doi",
                "http://www.w3.org/1999/xlink href": "10.1183/23120541.00193-2021",
            },
        )
        self.assertEqual(record["front"]["article-meta"]["kwd-group"], [])
        self.assertEqual(record["front"]["article-meta"]["contrib-group"], [])
        self.assertEqual(record["front"]["article-meta"]["issue-sponsor"], [])
        self.assertEqual(record["front"]["article-meta"]["self-uri"], [])
        self.assertEqual(record["front"]["article-meta"]["product"], [])
        self.assertEqual(record["front"]["article-meta"]["issue"], ["3"])
        self.assertEqual(record["front"]["article-meta"]["ext-link"], [])
        self.assertEqual(record["front"]["article-meta"]["support-group"], [])
        self.assertEqual(len(record["front"]["article-meta"]["article-id"]), 4)
        self.assertEqual(record["front"]["article-meta"]["article-id"][0], "34527728")
        self.assertEqual(
            record["front"]["article-meta"]["article-id"][0].attributes,
            {"pub-id-type": "pmid"},
        )
        self.assertEqual(record["front"]["article-meta"]["article-id"][1], "8435807")
        self.assertEqual(
            record["front"]["article-meta"]["article-id"][1].attributes,
            {"pub-id-type": "pmc"},
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-id"][2],
            "10.1183/23120541.50193-2021",
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-id"][2].attributes,
            {"pub-id-type": "doi"},
        )
        self.assertEqual(record["front"]["article-meta"]["article-id"][3], "50193-2021")
        self.assertEqual(
            record["front"]["article-meta"]["article-id"][3].attributes,
            {"pub-id-type": "publisher-id"},
        )
        self.assertEqual(record["front"]["article-meta"]["issue-title-group"], [])
        self.assertEqual(record["front"]["article-meta"]["x"], [])
        self.assertEqual(record["front"]["article-meta"]["uri"], [])
        self.assertEqual(record["front"]["article-meta"]["email"], [])
        self.assertEqual(record["front"]["article-meta"]["volume-id"], [])
        self.assertEqual(record["front"]["article-meta"]["issue-id"], [])
        self.assertEqual(record["front"]["article-meta"]["trans-abstract"], [])
        self.assertEqual(record["front"]["article-meta"]["volume-issue-group"], [])
        self.assertEqual(record["front"]["article-meta"]["related-object"], [])
        self.assertEqual(record["front"]["article-meta"]["isbn"], [])
        self.assertEqual(record["front"]["article-meta"]["volume"], ["7"])
        self.assertEqual(record["front"]["article-meta"]["aff-alternatives"], [])
        self.assertEqual(
            record["front"]["article-meta"]["article-version"], "Version of Record"
        )
        self.assertEqual(
            len(record["front"]["article-meta"]["article-version"].attributes), 3
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-version"].attributes["vocab"],
            "JAV",
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-version"].attributes[
                "vocab-identifier"
            ],
            "http://www.niso.org/publications/rp/RP-8-2008.pdf",
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-version"].attributes[
                "article-version-type"
            ],
            "VoR",
        )
        self.assertEqual(len(record["front"]["article-meta"]["article-categories"]), 3)
        self.assertEqual(
            record["front"]["article-meta"]["article-categories"]["series-text"], []
        )
        self.assertEqual(
            len(record["front"]["article-meta"]["article-categories"]["subj-group"]), 1
        )
        self.assertEqual(
            len(record["front"]["article-meta"]["article-categories"]["subj-group"][0]),
            3,
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-categories"]["subj-group"][0][
                "subject"
            ],
            ["Author Correction"],
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-categories"]["subj-group"][0][
                "subj-group"
            ],
            [],
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-categories"]["subj-group"][0][
                "compound-subject"
            ],
            [],
        )
        self.assertEqual(
            record["front"]["article-meta"]["article-categories"]["subj-group"][
                0
            ].attributes,
            {"subj-group-type": "heading"},
        )
        self.assertEqual(len(record["front"]["article-meta"]["title-group"]), 4)
        self.assertEqual(
            record["front"]["article-meta"]["title-group"]["trans-title-group"], []
        )
        self.assertEqual(
            record["front"]["article-meta"]["title-group"]["alt-title"], []
        )
        self.assertEqual(record["front"]["article-meta"]["title-group"]["subtitle"], [])
        self.assertEqual(
            record["front"]["article-meta"]["title-group"]["article-title"],
            '“Lung diffusing capacity for nitric oxide measured by two commercial devices: a randomised crossover comparison in healthy adults”. Thomas Radtke, Quintin de Groot, Sarah R. Haile, Marion Maggi, Connie C.W. Hsia and Holger Dressel. <italic toggle="yes">ERJ Open Res</italic> 2021; 7: 00193-2021.',
        )
        self.assertEqual(record["front"]["article-meta"]["elocation-id"], "50193-2021")
        self.assertEqual(len(record["front"]["article-meta"]["permissions"]), 5)
        self.assertEqual(
            record["front"]["article-meta"]["permissions"]["copyright-year"], ["2021"]
        )
        self.assertEqual(
            record["front"]["article-meta"]["permissions"]["copyright-holder"], []
        )
        self.assertEqual(
            len(record["front"]["article-meta"]["permissions"]["license"]), 1
        )
        self.assertEqual(
            len(record["front"]["article-meta"]["permissions"]["license"][0]), 2
        )
        self.assertEqual(
            record["front"]["article-meta"]["permissions"]["license"][0][0],
            "https://creativecommons.org/licenses/by-nc/4.0/",
        )
        self.assertEqual(
            record["front"]["article-meta"]["permissions"]["license"][0][0].attributes,
            {"specific-use": "textmining", "content-type": "ccbynclicense"},
        )
        self.assertEqual(
            record["front"]["article-meta"]["permissions"]["license"][0][1],
            'This version is distributed under the terms of the Creative Commons Attribution Non-Commercial Licence 4.0. For commercial reproduction rights and permissions contact <ext-link ext-link-type="uri" http://www.w3.org/1999/xlink href="mailto:permissions@ersnet.org">permissions@ersnet.org</ext-link>',
        )
        self.assertEqual(
            record["front"]["article-meta"]["permissions"]["copyright-statement"],
            ["Copyright ©The authors 2021"],
        )
        self.assertEqual(
            record["front"]["article-meta"]["permissions"]["ali:free_to_read"], []
        )
        self.assertEqual(record["front"]["glossary"], [])
        self.assertEqual(record["front"]["fn-group"], [])
        self.assertEqual(record["front"]["notes"], [])
        self.assertEqual(record["front"]["bio"], [])
        self.assertEqual(record["front"]["list"], [])
        self.assertEqual(record["front"]["def-list"], [])
        self.assertEqual(record["front"]["ack"], [])
        self.assertEqual(len(record["body"]), 37)

        self.assertEqual(record["body"]["table-wrap-group"], [])
        self.assertEqual(record["body"]["disp-formula"], [])
        self.assertEqual(record["body"]["answer-set"], [])
        self.assertEqual(record["body"]["graphic"], [])
        self.assertEqual(record["body"]["statement"], [])
        self.assertEqual(record["body"]["fig-group"], [])
        self.assertEqual(record["body"]["verse-group"], [])
        self.assertEqual(record["body"]["supplementary-material"], [])
        self.assertEqual(record["body"]["related-article"], [])
        self.assertEqual(record["body"]["code"], [])
        self.assertEqual(record["body"]["question"], [])
        self.assertEqual(record["body"]["preformat"], [])
        self.assertEqual(record["body"]["tex-math"], [])
        self.assertEqual(record["body"]["mml:math"], [])
        self.assertEqual(record["body"]["speech"], [])
        self.assertEqual(record["body"]["block-alternatives"], [])
        self.assertEqual(record["body"]["explanation"], [])
        self.assertEqual(record["body"]["array"], [])
        self.assertEqual(record["body"]["question-wrap-group"], [])
        self.assertEqual(record["body"]["alternatives"], [])
        self.assertEqual(record["body"]["media"], [])
        self.assertEqual(record["body"]["x"], [])
        self.assertEqual(record["body"]["sec"], [])
        self.assertEqual(record["body"]["address"], [])
        self.assertEqual(record["body"]["disp-quote"], [])
        self.assertEqual(record["body"]["table-wrap"], [])
        self.assertEqual(record["body"]["ack"], [])
        self.assertEqual(record["body"]["chem-struct-wrap"], [])
        self.assertEqual(record["body"]["related-object"], [])
        self.assertEqual(record["body"]["list"], [])
        self.assertEqual(record["body"]["def-list"], [])
        self.assertEqual(
            record["body"]["p"],
            [
                "This article was originally published with an error in table 2. The upper 95% confidence limit of the per cent difference in the primary end-point (diffusing capacity of the lung for nitric oxide) was incorrectly given as 15.1% and has now been corrected to −15.1% in the published article.\n"
            ],
        )
        self.assertEqual(record["body"]["fig"], [])
        self.assertEqual(record["body"]["answer"], [])
        self.assertEqual(record["body"]["boxed-text"], [])
        self.assertEqual(record["body"]["disp-formula-group"], [])
        self.assertEqual(record["body"]["question-wrap"], [])

    def test_taxonomy(self):
        # Access the Taxonomy database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="taxonomy", id="9615,9685", retmode="xml")
        with open("Entrez/taxonomy.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        self.assertEqual(len(record), 2)
        self.assertEqual(record[0]["TaxId"], "9615")
        self.assertEqual(record[0]["ScientificName"], "Canis lupus familiaris")
        self.assertEqual(record[0]["OtherNames"]["GenbankCommonName"], "dog")
        self.assertEqual(len(record[0]["OtherNames"]["Synonym"]), 5)
        self.assertEqual(record[0]["OtherNames"]["Synonym"][0], "Canis borealis")
        self.assertEqual(record[0]["OtherNames"]["Synonym"][1], "Canis canis")
        self.assertEqual(record[0]["OtherNames"]["Synonym"][2], "Canis domesticus")
        self.assertEqual(record[0]["OtherNames"]["Synonym"][3], "Canis familiaris")
        self.assertEqual(record[0]["OtherNames"]["Synonym"][4], "Canis lupus borealis")
        self.assertEqual(len(record[0]["OtherNames"]["CommonName"]), 1)
        self.assertEqual(record[0]["OtherNames"]["CommonName"][0], "dogs")
        self.assertEqual(record[0]["OtherNames"]["Includes"][0], "beagle dog")
        self.assertEqual(record[0]["OtherNames"]["Includes"][1], "beagle dogs")
        self.assertEqual(record[0]["ParentTaxId"], "9612")
        self.assertEqual(record[0]["Rank"], "subspecies")
        self.assertEqual(record[0]["Division"], "Mammals")
        self.assertEqual(record[0]["GeneticCode"]["GCId"], "1")
        self.assertEqual(record[0]["GeneticCode"]["GCName"], "Standard")
        self.assertEqual(record[0]["MitoGeneticCode"]["MGCId"], "2")
        self.assertEqual(
            record[0]["MitoGeneticCode"]["MGCName"], "Vertebrate Mitochondrial"
        )
        self.assertEqual(
            record[0]["Lineage"],
            "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Laurasiatheria; Carnivora; Caniformia; Canidae; Canis; Canis lupus",
        )

        self.assertEqual(record[0]["LineageEx"][0]["TaxId"], "131567")
        self.assertEqual(
            record[0]["LineageEx"][0]["ScientificName"], "cellular organisms"
        )
        self.assertEqual(record[0]["LineageEx"][0]["Rank"], "cellular root")
        self.assertEqual(record[0]["LineageEx"][1]["TaxId"], "2759")
        self.assertEqual(record[0]["LineageEx"][1]["ScientificName"], "Eukaryota")
        self.assertEqual(record[0]["LineageEx"][1]["Rank"], "domain")
        self.assertEqual(record[0]["LineageEx"][2]["TaxId"], "33154")
        self.assertEqual(record[0]["LineageEx"][2]["ScientificName"], "Opisthokonta")
        self.assertEqual(record[0]["LineageEx"][2]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][3]["TaxId"], "33208")
        self.assertEqual(record[0]["LineageEx"][3]["ScientificName"], "Metazoa")
        self.assertEqual(record[0]["LineageEx"][3]["Rank"], "kingdom")
        self.assertEqual(record[0]["LineageEx"][4]["TaxId"], "6072")
        self.assertEqual(record[0]["LineageEx"][4]["ScientificName"], "Eumetazoa")
        self.assertEqual(record[0]["LineageEx"][4]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][5]["TaxId"], "33213")
        self.assertEqual(record[0]["LineageEx"][5]["ScientificName"], "Bilateria")
        self.assertEqual(record[0]["LineageEx"][5]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][6]["TaxId"], "33511")
        self.assertEqual(record[0]["LineageEx"][6]["ScientificName"], "Deuterostomia")
        self.assertEqual(record[0]["LineageEx"][6]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][7]["TaxId"], "7711")
        self.assertEqual(record[0]["LineageEx"][7]["ScientificName"], "Chordata")
        self.assertEqual(record[0]["LineageEx"][7]["Rank"], "phylum")
        self.assertEqual(record[0]["LineageEx"][8]["TaxId"], "89593")
        self.assertEqual(record[0]["LineageEx"][8]["ScientificName"], "Craniata")
        self.assertEqual(record[0]["LineageEx"][8]["Rank"], "subphylum")
        self.assertEqual(record[0]["LineageEx"][9]["TaxId"], "7742")
        self.assertEqual(record[0]["LineageEx"][9]["ScientificName"], "Vertebrata")
        self.assertEqual(record[0]["LineageEx"][9]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][10]["TaxId"], "7776")
        self.assertEqual(record[0]["LineageEx"][10]["ScientificName"], "Gnathostomata")
        self.assertEqual(record[0]["LineageEx"][10]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][11]["TaxId"], "117570")
        self.assertEqual(record[0]["LineageEx"][11]["ScientificName"], "Teleostomi")
        self.assertEqual(record[0]["LineageEx"][11]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][12]["TaxId"], "117571")
        self.assertEqual(record[0]["LineageEx"][12]["ScientificName"], "Euteleostomi")
        self.assertEqual(record[0]["LineageEx"][12]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][13]["TaxId"], "8287")
        self.assertEqual(record[0]["LineageEx"][13]["ScientificName"], "Sarcopterygii")
        self.assertEqual(record[0]["LineageEx"][13]["Rank"], "superclass")
        self.assertEqual(record[0]["LineageEx"][14]["TaxId"], "1338369")
        self.assertEqual(record[0]["LineageEx"][14]["ScientificName"], "Dipnotetrapodomorpha")
        self.assertEqual(record[0]["LineageEx"][14]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][15]["TaxId"], "32523")
        self.assertEqual(record[0]["LineageEx"][15]["ScientificName"], "Tetrapoda")
        self.assertEqual(record[0]["LineageEx"][15]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][16]["TaxId"], "32524")
        self.assertEqual(record[0]["LineageEx"][16]["ScientificName"], "Amniota")
        self.assertEqual(record[0]["LineageEx"][16]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][17]["TaxId"], "40674")
        self.assertEqual(record[0]["LineageEx"][17]["ScientificName"], "Mammalia")
        self.assertEqual(record[0]["LineageEx"][17]["Rank"], "class")
        self.assertEqual(record[0]["LineageEx"][18]["TaxId"], "32525")
        self.assertEqual(record[0]["LineageEx"][18]["ScientificName"], "Theria")
        self.assertEqual(record[0]["LineageEx"][18]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][19]["TaxId"], "9347")
        self.assertEqual(record[0]["LineageEx"][19]["ScientificName"], "Eutheria")
        self.assertEqual(record[0]["LineageEx"][19]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][20]["TaxId"], "1437010")
        self.assertEqual(record[0]["LineageEx"][20]["ScientificName"], "Boreoeutheria")
        self.assertEqual(record[0]["LineageEx"][20]["Rank"], "clade")
        self.assertEqual(record[0]["LineageEx"][21]["TaxId"], "314145")
        self.assertEqual(record[0]["LineageEx"][21]["ScientificName"], "Laurasiatheria")
        self.assertEqual(record[0]["LineageEx"][21]["Rank"], "superorder")
        self.assertEqual(record[0]["LineageEx"][22]["TaxId"], "33554")
        self.assertEqual(record[0]["LineageEx"][22]["ScientificName"], "Carnivora")
        self.assertEqual(record[0]["LineageEx"][22]["Rank"], "order")
        self.assertEqual(record[0]["LineageEx"][23]["TaxId"], "379584")
        self.assertEqual(record[0]["LineageEx"][23]["ScientificName"], "Caniformia")
        self.assertEqual(record[0]["LineageEx"][23]["Rank"], "suborder")
        self.assertEqual(record[0]["LineageEx"][24]["TaxId"], "9608")
        self.assertEqual(record[0]["LineageEx"][24]["ScientificName"], "Canidae")
        self.assertEqual(record[0]["LineageEx"][24]["Rank"], "family")
        self.assertEqual(record[0]["LineageEx"][25]["TaxId"], "9611")
        self.assertEqual(record[0]["LineageEx"][25]["ScientificName"], "Canis")
        self.assertEqual(record[0]["LineageEx"][25]["Rank"], "genus")
        self.assertEqual(record[0]["LineageEx"][26]["TaxId"], "9612")
        self.assertEqual(record[0]["LineageEx"][26]["ScientificName"], "Canis lupus")
        self.assertEqual(record[0]["LineageEx"][26]["Rank"], "species")
        self.assertEqual(record[0]["CreateDate"], "1995/02/27 09:24:00")
        self.assertEqual(record[0]["UpdateDate"], "2024/02/09 13:11:20")
        self.assertEqual(record[0]["PubDate"], "1993/04/27 01:00:00")
        self.assertEqual(record[1]["TaxId"], "9685")
        self.assertEqual(record[1]["ScientificName"], "Felis catus")
        self.assertEqual(record[1]["OtherNames"]["GenbankCommonName"], "domestic cat")
        self.assertEqual(len(record[1]["OtherNames"]["Synonym"]), 2)
        self.assertEqual(record[1]["OtherNames"]["Synonym"][0], "Felis domesticus")
        self.assertEqual(record[1]["OtherNames"]["Synonym"][1], "Felis silvestris catus")
        self.assertEqual(len(record[1]["OtherNames"]["CommonName"]), 2)
        self.assertEqual(record[1]["OtherNames"]["CommonName"][0], "cat")
        self.assertEqual(record[1]["OtherNames"]["CommonName"][1], "cats")
        self.assertEqual(record[1]["OtherNames"]["Includes"][0], "Korat cats")
        self.assertEqual(record[1]["ParentTaxId"], "9682")
        self.assertEqual(record[1]["Rank"], "species")
        self.assertEqual(record[1]["Division"], "Mammals")
        self.assertEqual(record[1]["GeneticCode"]["GCId"], "1")
        self.assertEqual(record[1]["GeneticCode"]["GCName"], "Standard")
        self.assertEqual(record[1]["MitoGeneticCode"]["MGCId"], "2")
        self.assertEqual(
            record[1]["MitoGeneticCode"]["MGCName"], "Vertebrate Mitochondrial"
        )
        self.assertEqual(
            record[1]["Lineage"],
            "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Laurasiatheria; Carnivora; Feliformia; Felidae; Felinae; Felis",
        )

        self.assertEqual(record[1]["LineageEx"][0]["TaxId"], "131567")
        self.assertEqual(
            record[1]["LineageEx"][0]["ScientificName"], "cellular organisms"
        )
        self.assertEqual(record[1]["LineageEx"][0]["Rank"], "cellular root")
        self.assertEqual(record[1]["LineageEx"][1]["TaxId"], "2759")
        self.assertEqual(record[1]["LineageEx"][1]["ScientificName"], "Eukaryota")
        self.assertEqual(record[1]["LineageEx"][1]["Rank"], "domain")
        self.assertEqual(record[1]["LineageEx"][2]["TaxId"], "33154")
        self.assertEqual(record[1]["LineageEx"][2]["ScientificName"], "Opisthokonta")
        self.assertEqual(record[1]["LineageEx"][2]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][3]["TaxId"], "33208")
        self.assertEqual(record[1]["LineageEx"][3]["ScientificName"], "Metazoa")
        self.assertEqual(record[1]["LineageEx"][3]["Rank"], "kingdom")
        self.assertEqual(record[1]["LineageEx"][4]["TaxId"], "6072")
        self.assertEqual(record[1]["LineageEx"][4]["ScientificName"], "Eumetazoa")
        self.assertEqual(record[1]["LineageEx"][4]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][5]["TaxId"], "33213")
        self.assertEqual(record[1]["LineageEx"][5]["ScientificName"], "Bilateria")
        self.assertEqual(record[1]["LineageEx"][5]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][6]["TaxId"], "33511")
        self.assertEqual(record[1]["LineageEx"][6]["ScientificName"], "Deuterostomia")
        self.assertEqual(record[1]["LineageEx"][6]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][7]["TaxId"], "7711")
        self.assertEqual(record[1]["LineageEx"][7]["ScientificName"], "Chordata")
        self.assertEqual(record[1]["LineageEx"][7]["Rank"], "phylum")
        self.assertEqual(record[1]["LineageEx"][8]["TaxId"], "89593")
        self.assertEqual(record[1]["LineageEx"][8]["ScientificName"], "Craniata")
        self.assertEqual(record[1]["LineageEx"][8]["Rank"], "subphylum")
        self.assertEqual(record[1]["LineageEx"][9]["TaxId"], "7742")
        self.assertEqual(record[1]["LineageEx"][9]["ScientificName"], "Vertebrata")
        self.assertEqual(record[1]["LineageEx"][9]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][10]["TaxId"], "7776")
        self.assertEqual(record[1]["LineageEx"][10]["ScientificName"], "Gnathostomata")
        self.assertEqual(record[1]["LineageEx"][10]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][11]["TaxId"], "117570")
        self.assertEqual(record[1]["LineageEx"][11]["ScientificName"], "Teleostomi")
        self.assertEqual(record[1]["LineageEx"][11]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][12]["TaxId"], "117571")
        self.assertEqual(record[1]["LineageEx"][12]["ScientificName"], "Euteleostomi")
        self.assertEqual(record[1]["LineageEx"][12]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][13]["TaxId"], "8287")
        self.assertEqual(record[1]["LineageEx"][13]["ScientificName"], "Sarcopterygii")
        self.assertEqual(record[1]["LineageEx"][13]["Rank"], "superclass")
        self.assertEqual(record[1]["LineageEx"][14]["TaxId"], "1338369")
        self.assertEqual(record[1]["LineageEx"][14]["ScientificName"], "Dipnotetrapodomorpha")
        self.assertEqual(record[1]["LineageEx"][14]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][15]["TaxId"], "32523")
        self.assertEqual(record[1]["LineageEx"][15]["ScientificName"], "Tetrapoda")
        self.assertEqual(record[1]["LineageEx"][15]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][16]["TaxId"], "32524")
        self.assertEqual(record[1]["LineageEx"][16]["ScientificName"], "Amniota")
        self.assertEqual(record[1]["LineageEx"][16]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][17]["TaxId"], "40674")
        self.assertEqual(record[1]["LineageEx"][17]["ScientificName"], "Mammalia")
        self.assertEqual(record[1]["LineageEx"][17]["Rank"], "class")
        self.assertEqual(record[1]["LineageEx"][18]["TaxId"], "32525")
        self.assertEqual(record[1]["LineageEx"][18]["ScientificName"], "Theria")
        self.assertEqual(record[1]["LineageEx"][18]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][19]["TaxId"], "9347")
        self.assertEqual(record[1]["LineageEx"][19]["ScientificName"], "Eutheria")
        self.assertEqual(record[1]["LineageEx"][19]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][20]["TaxId"], "1437010")
        self.assertEqual(record[1]["LineageEx"][20]["ScientificName"], "Boreoeutheria")
        self.assertEqual(record[1]["LineageEx"][20]["Rank"], "clade")
        self.assertEqual(record[1]["LineageEx"][21]["TaxId"], "314145")
        self.assertEqual(record[1]["LineageEx"][21]["ScientificName"], "Laurasiatheria")
        self.assertEqual(record[1]["LineageEx"][21]["Rank"], "superorder")
        self.assertEqual(record[1]["LineageEx"][22]["TaxId"], "33554")
        self.assertEqual(record[1]["LineageEx"][22]["ScientificName"], "Carnivora")
        self.assertEqual(record[1]["LineageEx"][22]["Rank"], "order")
        self.assertEqual(record[1]["LineageEx"][23]["TaxId"], "379583")
        self.assertEqual(record[1]["LineageEx"][23]["ScientificName"], "Feliformia")
        self.assertEqual(record[1]["LineageEx"][23]["Rank"], "suborder")
        self.assertEqual(record[1]["LineageEx"][24]["TaxId"], "9681")
        self.assertEqual(record[1]["LineageEx"][24]["ScientificName"], "Felidae")
        self.assertEqual(record[1]["LineageEx"][24]["Rank"], "family")
        self.assertEqual(record[1]["LineageEx"][25]["TaxId"], "338152")
        self.assertEqual(record[1]["LineageEx"][25]["ScientificName"], "Felinae")
        self.assertEqual(record[1]["LineageEx"][25]["Rank"], "subfamily")
        self.assertEqual(record[1]["LineageEx"][26]["TaxId"], "9682")
        self.assertEqual(record[1]["LineageEx"][26]["ScientificName"], "Felis")
        self.assertEqual(record[1]["LineageEx"][26]["Rank"], "genus")
        self.assertEqual(record[1]["CreateDate"], "1995/02/27 09:24:00")
        self.assertEqual(record[1]["UpdateDate"], "2024/03/03 11:27:08")
        self.assertEqual(record[1]["PubDate"], "1993/07/26 01:00:00")
        # fmt: on

    def test_nucleotide1(self):
        """Test parsing XML returned by EFetch, Nucleotide database (first test)."""
        # Access the nucleotide database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='nucleotide', id=5, retmode='xml')
        with open("Entrez/nucleotide1.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        self.assertEqual(record[0]["GBSeq_locus"], "X60065")
        self.assertEqual(record[0]["GBSeq_length"], "1136")
        self.assertEqual(record[0]["GBSeq_strandedness"], "single")
        self.assertEqual(record[0]["GBSeq_moltype"], "mRNA")
        self.assertEqual(record[0]["GBSeq_topology"], "linear")
        self.assertEqual(record[0]["GBSeq_division"], "MAM")
        self.assertEqual(record[0]["GBSeq_update-date"], "26-JUL-2016")
        self.assertEqual(record[0]["GBSeq_create-date"], "05-MAY-1992")
        self.assertEqual(
            record[0]["GBSeq_definition"],
            "B.bovis beta-2-gpI mRNA for beta-2-glycoprotein I",
        )
        self.assertEqual(record[0]["GBSeq_primary-accession"], "X60065")
        self.assertEqual(record[0]["GBSeq_accession-version"], "X60065.1")
        self.assertEqual(record[0]["GBSeq_other-seqids"][0], "emb|X60065.1|")
        self.assertEqual(record[0]["GBSeq_other-seqids"][1], "gi|5")
        self.assertEqual(record[0]["GBSeq_keywords"][0], "beta-2 glycoprotein I")
        self.assertEqual(record[0]["GBSeq_source"], "Bos taurus (domestic cattle)")
        self.assertEqual(record[0]["GBSeq_organism"], "Bos taurus")
        self.assertEqual(
            record[0]["GBSeq_taxonomy"],
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Artiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos",
        )
        self.assertEqual(record[0]["GBSeq_references"][0]["GBReference_reference"], "1")
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][0], "Bendixen,E."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][1], "Halkier,T."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][2], "Magnusson,S."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][3], "Sottrup-Jensen,L.",
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][4], "Kristensen,T."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_title"],
            "Complete primary structure of bovine beta 2-glycoprotein I: localization of the disulfide bridges",
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_journal"],
            "Biochemistry 31 (14), 3611-3617 (1992)",
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_pubmed"], "1567819"
        )
        self.assertEqual(record[0]["GBSeq_references"][1]["GBReference_reference"], "2")
        self.assertEqual(
            record[0]["GBSeq_references"][1]["GBReference_position"], "1..1136"
        )
        self.assertEqual(
            record[0]["GBSeq_references"][1]["GBReference_authors"][0], "Kristensen,T."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][1]["GBReference_title"], "Direct Submission"
        )
        self.assertEqual(
            record[0]["GBSeq_references"][1]["GBReference_journal"],
            "Submitted (11-JUN-1991) T. Kristensen, Dept of Mol Biology, University of Aarhus, C F Mollers Alle 130, DK-8000 Aarhus C, DENMARK",
        )
        self.assertEqual(len(record[0]["GBSeq_feature-table"]), 7)
        self.assertEqual(record[0]["GBSeq_feature-table"][0]["GBFeature_key"], "source")
        self.assertEqual(record[0]["GBSeq_feature-table"][0]["GBFeature_location"], "1..1136")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_from"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_to"],
            "1136",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_accession"],
            "X60065.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_name"],
            "organism",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_value"],
            "Bos taurus",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_name"],
            "mol_type",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_value"],
            "mRNA",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_value"],
            "taxon:9913",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_name"],
            "clone",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_value"],
            "pBB2I",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_name"],
            "tissue_type",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_value"],
            "liver",
        )
        self.assertEqual(record[0]["GBSeq_feature-table"][1]["GBFeature_key"], "gene")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_location"], "<1..1136"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_from"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_to"],
            "1136",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_accession"],
            "X60065.1",
        )
        self.assertEqual(record[0]["GBSeq_feature-table"][1]["GBFeature_partial5"], "")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_partial5"].attributes["value"],
            "true",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_name"],
            "gene",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_value"],
            "beta-2-gpI",
        )
        self.assertEqual(record[0]["GBSeq_feature-table"][2]["GBFeature_key"], "CDS")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_location"], "<1..1029"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_from"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_to"],
            "1029",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_accession"],
            "X60065.1",
        )
        self.assertEqual(record[0]["GBSeq_feature-table"][2]["GBFeature_partial5"], "")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_partial5"].attributes["value"],
            "true",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_name"],
            "gene",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_value"],
            "beta-2-gpI",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_name"],
            "codon_start",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_value"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_name"],
            "transl_table",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_value"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][3]["GBQualifier_name"],
            "product",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][3]["GBQualifier_value"],
            "beta-2-glycoprotein I",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][4]["GBQualifier_name"],
            "protein_id",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][4]["GBQualifier_value"],
            "CAA42669.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][5]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][5]["GBQualifier_value"],
            "GOA:P17690",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][6]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][6]["GBQualifier_value"],
            "InterPro:IPR000436",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][7]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][7]["GBQualifier_value"],
            "InterPro:IPR015104",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][8]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][8]["GBQualifier_value"],
            "InterPro:IPR016060",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_key"], "sig_peptide"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_location"], "<1..48"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_from"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_to"],
            "48",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_accession"],
            "X60065.1",
        )
        self.assertEqual(record[0]["GBSeq_feature-table"][3]["GBFeature_partial5"], "")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_partial5"].attributes["value"],
            "true",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_name"],
            "gene",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_value"],
            "beta-2-gpI",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_key"], "mat_peptide"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_location"], "49..1026"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_from"],
            "49",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_to"],
            "1026",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_accession"],
            "X60065.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][0]["GBQualifier_name"],
            "gene",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][0]["GBQualifier_value"],
            "beta-2-gpI",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][1]["GBQualifier_name"],
            "product",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][1]["GBQualifier_value"],
            "beta-2-glycoprotein I",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][2]["GBQualifier_name"],
            "peptide",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][2]["GBQualifier_value"],
            "GRTCPKPDELPFSTVVPLKRTYEPGEQIVFSCQPGYVSRGGIRRFTCPLTGLWPINTLKCMPRVCPFAGILENGTVRYTTFEYPNTISFSCHTGFYLKGASSAKCTEEGKWSPDLPVCAPITCPPPPIPKFASLSVYKPLAGNNSFYGSKAVFKCLPHHAMFGNDTVTCTEHGNWTQLPECREVRCPFPSRPDNGFVNHPANPVLYYKDTATFGCHETYSLDGPEEVECSKFGNWSAQPSCKASCKLSIKRATVIYEGERVAIQNKFKNGMLHGQKVSFFCKHKEKKCSYTEDAQCIDGTIEIPKCFKEHSSLAFWKTDASDVKPC",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_key"], "regulatory",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_location"], "1101..1106"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_from"],
            "1101",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_to"],
            "1106",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_accession"],
            "X60065.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][0]["GBQualifier_name"],
            "regulatory_class",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][0]["GBQualifier_value"],
            "polyA_signal_sequence",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][1]["GBQualifier_name"],
            "gene",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][1]["GBQualifier_value"],
            "beta-2-gpI",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][6]["GBFeature_key"], "polyA_site"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][6]["GBFeature_location"], "1130"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][6]["GBFeature_intervals"][0]["GBInterval_point"],
            "1130",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][6]["GBFeature_intervals"][0]["GBInterval_accession"],
            "X60065.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][6]["GBFeature_quals"][0]["GBQualifier_name"],
            "gene",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][6]["GBFeature_quals"][0]["GBQualifier_value"],
            "beta-2-gpI",
        )
        self.assertEqual(
            record[0]["GBSeq_sequence"],
            "ccagcgctcgtcttgctgttggggtttctctgccacgttgctatcgcaggacgaacctgccccaagccagatgagctaccgttttccacggtggttccactgaaacggacctatgagcccggggagcagatagtcttctcctgccagccgggctacgtgtcccggggagggatccggcggtttacatgcccgctcacaggactctggcccatcaacacgctgaaatgcatgcccagagtatgtccttttgctgggatcttagaaaacggaacggtacgctatacaacgtttgagtatcccaacaccatcagcttttcttgccacacggggttttatctgaaaggagctagttctgcaaaatgcactgaggaagggaagtggagcccagaccttcctgtctgtgcccctataacctgccctccaccacccatacccaagtttgcaagtctcagcgtttacaagccgttggctgggaacaactccttctatggcagcaaggcagtctttaagtgcttgccacaccacgcgatgtttggaaatgacaccgttacctgcacggaacatgggaactggacgcagttgccagaatgcagggaagtaagatgcccattcccatcaagaccagacaatgggtttgtgaaccatcctgcaaatccagtgctctactataaggacaccgccacctttggctgccatgaaacgtattccttggatggaccggaagaagtagaatgcagcaaattcggaaactggtctgcacagccaagctgtaaagcatcttgtaagttatctattaaaagagctactgtgatatatgaaggagagagagtagctatccagaacaaatttaagaatggaatgctgcatggccaaaaggtttctttcttctgcaagcataaggaaaagaagtgcagctacacagaagatgctcagtgcatagacggcaccatcgagattcccaaatgcttcaaggagcacagttctttagctttctggaaaacggatgcatctgacgtaaaaccatgctaagctggttttcacactgaaaattaaatgtcatgcttatatgtgtctgtctgagaatctgatggaaacggaaaaataaagagactgaatttaccgtgtcaagaaaaaaa",
        )
        # fmt: off

    def test_nucleotide2(self):
        """Test parsing XML returned by EFetch, Nucleotide database (second test)."""
        # Access the nucleotide database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='nucleotide', id=5,
        #                       rettype='fasta', complexity=0, retmode='xml')
        with open("Entrez/nucleotide2.xml", "rb") as stream:
            record = Entrez.read(stream)
        self.assertEqual(record[0]["TSeq_seqtype"], "")
        self.assertEqual(record[0]["TSeq_seqtype"].attributes["value"], "nucleotide")
        self.assertEqual(record[0]["TSeq_accver"], "X60065.1")
        self.assertEqual(record[0]["TSeq_taxid"], "9913")
        self.assertEqual(record[0]["TSeq_orgname"], "Bos taurus")
        self.assertEqual(
            record[0]["TSeq_defline"],
            "B.bovis beta-2-gpI mRNA for beta-2-glycoprotein I",
        )
        self.assertEqual(record[0]["TSeq_length"], "1136")
        self.assertEqual(
            record[0]["TSeq_sequence"],
            "CCAGCGCTCGTCTTGCTGTTGGGGTTTCTCTGCCACGTTGCTATCGCAGGACGAACCTGCCCCAAGCCAGATGAGCTACCGTTTTCCACGGTGGTTCCACTGAAACGGACCTATGAGCCCGGGGAGCAGATAGTCTTCTCCTGCCAGCCGGGCTACGTGTCCCGGGGAGGGATCCGGCGGTTTACATGCCCGCTCACAGGACTCTGGCCCATCAACACGCTGAAATGCATGCCCAGAGTATGTCCTTTTGCTGGGATCTTAGAAAACGGAACGGTACGCTATACAACGTTTGAGTATCCCAACACCATCAGCTTTTCTTGCCACACGGGGTTTTATCTGAAAGGAGCTAGTTCTGCAAAATGCACTGAGGAAGGGAAGTGGAGCCCAGACCTTCCTGTCTGTGCCCCTATAACCTGCCCTCCACCACCCATACCCAAGTTTGCAAGTCTCAGCGTTTACAAGCCGTTGGCTGGGAACAACTCCTTCTATGGCAGCAAGGCAGTCTTTAAGTGCTTGCCACACCACGCGATGTTTGGAAATGACACCGTTACCTGCACGGAACATGGGAACTGGACGCAGTTGCCAGAATGCAGGGAAGTAAGATGCCCATTCCCATCAAGACCAGACAATGGGTTTGTGAACCATCCTGCAAATCCAGTGCTCTACTATAAGGACACCGCCACCTTTGGCTGCCATGAAACGTATTCCTTGGATGGACCGGAAGAAGTAGAATGCAGCAAATTCGGAAACTGGTCTGCACAGCCAAGCTGTAAAGCATCTTGTAAGTTATCTATTAAAAGAGCTACTGTGATATATGAAGGAGAGAGAGTAGCTATCCAGAACAAATTTAAGAATGGAATGCTGCATGGCCAAAAGGTTTCTTTCTTCTGCAAGCATAAGGAAAAGAAGTGCAGCTACACAGAAGATGCTCAGTGCATAGACGGCACCATCGAGATTCCCAAATGCTTCAAGGAGCACAGTTCTTTAGCTTTCTGGAAAACGGATGCATCTGACGTAAAACCATGCTAAGCTGGTTTTCACACTGAAAATTAAATGTCATGCTTATATGTGTCTGTCTGAGAATCTGATGGAAACGGAAAAATAAAGAGACTGAATTTACCGTGTCAAGAAAAAAA",
        )
        self.assertEqual(record[1]["TSeq_seqtype"], "")
        self.assertEqual(record[1]["TSeq_seqtype"].attributes["value"], "protein")
        self.assertEqual(record[1]["TSeq_accver"], "CAA42669.1")
        self.assertEqual(record[1]["TSeq_taxid"], "9913")
        self.assertEqual(record[1]["TSeq_orgname"], "Bos taurus")
        self.assertEqual(
            record[1]["TSeq_defline"], "beta-2-glycoprotein I, partial [Bos taurus]"
        )
        self.assertEqual(record[1]["TSeq_length"], "342")
        self.assertEqual(
            record[1]["TSeq_sequence"],
            "PALVLLLGFLCHVAIAGRTCPKPDELPFSTVVPLKRTYEPGEQIVFSCQPGYVSRGGIRRFTCPLTGLWPINTLKCMPRVCPFAGILENGTVRYTTFEYPNTISFSCHTGFYLKGASSAKCTEEGKWSPDLPVCAPITCPPPPIPKFASLSVYKPLAGNNSFYGSKAVFKCLPHHAMFGNDTVTCTEHGNWTQLPECREVRCPFPSRPDNGFVNHPANPVLYYKDTATFGCHETYSLDGPEEVECSKFGNWSAQPSCKASCKLSIKRATVIYEGERVAIQNKFKNGMLHGQKVSFFCKHKEKKCSYTEDAQCIDGTIEIPKCFKEHSSLAFWKTDASDVKPC",
        )

    def test_protein(self):
        """Test parsing XML returned by EFetch, Protein database."""
        # Access the protein database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='protein', id=8, rettype='gp', retmode='xml')
        with open("Entrez/protein.xml", "rb") as stream:
            record = Entrez.read(stream)
        # fmt: off
        self.assertEqual(record[0]["GBSeq_locus"], "CAA35997")
        self.assertEqual(record[0]["GBSeq_length"], "100")
        self.assertEqual(record[0]["GBSeq_moltype"], "AA")
        self.assertEqual(record[0]["GBSeq_topology"], "linear")
        self.assertEqual(record[0]["GBSeq_division"], "MAM")
        self.assertEqual(record[0]["GBSeq_update-date"], "12-SEP-1993")
        self.assertEqual(record[0]["GBSeq_create-date"], "03-APR-1990")
        self.assertEqual(
            record[0]["GBSeq_definition"], "unnamed protein product [Bos taurus]"
        )
        self.assertEqual(record[0]["GBSeq_primary-accession"], "CAA35997")
        self.assertEqual(record[0]["GBSeq_accession-version"], "CAA35997.1")
        self.assertEqual(record[0]["GBSeq_other-seqids"][0], "emb|CAA35997.1|")
        self.assertEqual(record[0]["GBSeq_other-seqids"][1], "gi|8")
        self.assertEqual(record[0]["GBSeq_source"], "Bos taurus (domestic cattle)")
        self.assertEqual(record[0]["GBSeq_organism"], "Bos taurus")
        self.assertEqual(
            record[0]["GBSeq_taxonomy"],
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Artiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos",
        )
        self.assertEqual(record[0]["GBSeq_references"][0]["GBReference_reference"], "1")
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_position"], "1..100"
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][0], "Kiefer,M.C."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][1], "Saphire,A.C.S."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][2], "Bauer,D.M."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_authors"][3], "Barr,P.J."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][0]["GBReference_journal"], "Unpublished"
        )
        self.assertEqual(record[0]["GBSeq_references"][1]["GBReference_reference"], "2")
        self.assertEqual(
            record[0]["GBSeq_references"][1]["GBReference_position"], "1..100"
        )
        self.assertEqual(
            record[0]["GBSeq_references"][1]["GBReference_authors"][0], "Kiefer,M.C."
        )
        self.assertEqual(
            record[0]["GBSeq_references"][1]["GBReference_title"], "Direct Submission"
        )
        self.assertEqual(
            record[0]["GBSeq_references"][1]["GBReference_journal"],
            "Submitted (30-JAN-1990) Kiefer M.C., Chiron Corporation, 4560 Hortom St, Emeryville CA 94608-2916, U S A",
        )
        self.assertEqual(
            record[0]["GBSeq_comment"],
            "See <X15699> for Human sequence.~~Data kindly reviewed (08-MAY-1990) by Kiefer M.C.",
        )
        self.assertEqual(record[0]["GBSeq_source-db"], "embl accession X51700.1")
        self.assertEqual(record[0]["GBSeq_feature-table"][0]["GBFeature_key"], "source")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_location"], "1..100"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_from"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_to"],
            "100",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_accession"],
            "CAA35997.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_name"],
            "organism",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_value"],
            "Bos taurus",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_value"],
            "taxon:9913",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_name"],
            "clone",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_value"],
            "bBGP-3",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_name"],
            "tissue_type",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_value"],
            "bone matrix",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_name"],
            "clone_lib",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_value"],
            "Zap-bb",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_key"], "Protein"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_location"], "1..100"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_from"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_to"],
            "100",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_accession"],
            "CAA35997.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_name"],
            "name",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_value"],
            "unnamed protein product",
        )
        self.assertEqual(record[0]["GBSeq_feature-table"][2]["GBFeature_key"], "Region")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_location"], "33..97"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_from"],
            "33",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_to"],
            "97",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_accession"],
            "CAA35997.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_name"],
            "region_name",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_value"],
            "GLA",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_name"],
            "note",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_value"],
            "Domain containing Gla (gamma-carboxyglutamate) residues; smart00069",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_value"],
            "CDD:214503",
        )
        self.assertEqual(record[0]["GBSeq_feature-table"][3]["GBFeature_key"], "CDS")
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_location"], "1..100"
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_from"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_to"],
            "100",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_accession"],
            "CAA35997.1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_name"],
            "coded_by",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_value"],
            "X51700.1:28..330",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][1]["GBQualifier_name"],
            "note",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][1]["GBQualifier_value"],
            "bone Gla precursor (100 AA)",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][2]["GBQualifier_name"],
            "transl_table",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][2]["GBQualifier_value"],
            "1",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][3]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][3]["GBQualifier_value"],
            "GOA:P02820",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][4]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][4]["GBQualifier_value"],
            "InterPro:IPR000294",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][5]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][5]["GBQualifier_value"],
            "InterPro:IPR002384",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][6]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][6]["GBQualifier_value"],
            "PDB:1Q3M",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][7]["GBQualifier_name"],
            "db_xref",
        )
        self.assertEqual(
            record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][7]["GBQualifier_value"],
            "UniProtKB/Swiss-Prot:P02820",
        )
        self.assertEqual(
            record[0]["GBSeq_sequence"],
            "mrtpmllallalatlclagradakpgdaesgkgaafvskqegsevvkrlrryldhwlgapapypdplepkrevcelnpdcdeladhigfqeayrrfygpv",
        )
        # fmt: on

    def test_efetch_schemas(self):
        """Test parsing XML using Schemas."""
        # To create the XML file,use
        # >>> Bio.Entrez.efetch("protein", id="783730874", rettype="ipg", retmode="xml")
        with open("Entrez/efetch_schemas.xml", "rb") as stream:
            records = Entrez.read(stream)
        self.assertEqual(len(records), 1)
        record = records["IPGReport"]
        self.assertEqual(len(record.attributes), 2)
        self.assertEqual(record.attributes["product_acc"], "KJV04014.1")
        self.assertEqual(record.attributes["ipg"], "79092155")
        self.assertEqual(len(record), 3)
        self.assertEqual(record["Product"], "")
        self.assertEqual(record["Product"].attributes["kingdom"], "Bacteria")
        self.assertEqual(record["Product"].attributes["slen"], "513")
        self.assertEqual(
            record["Product"].attributes["name"],
            "methylmalonate-semialdehyde dehydrogenase (CoA acylating)",
        )
        self.assertEqual(record["Product"].attributes["org"], "Rhodococcus sp. PML026")
        self.assertEqual(record["Product"].attributes["kingdom_taxid"], "2")
        self.assertEqual(record["Product"].attributes["accver"], "WP_045840896.1")
        self.assertEqual(record["Product"].attributes["taxid"], "1356405")
        self.assertEqual(len(record["ProteinList"]), 2)
        protein = record["ProteinList"][0]
        self.assertEqual(protein.tag, "Protein")
        self.assertEqual(protein.attributes["accver"], "KJV04014.1")
        self.assertEqual(protein.attributes["kingdom"], "Bacteria")
        self.assertEqual(protein.attributes["kingdom_taxid"], "2")
        self.assertEqual(
            protein.attributes["name"],
            "methylmalonic acid semialdehyde dehydrogenase mmsa",
        )
        self.assertEqual(protein.attributes["org"], "Rhodococcus sp. PML026")
        self.assertEqual(protein.attributes["priority"], "0")
        self.assertEqual(protein.attributes["source"], "INSDC")
        self.assertEqual(protein.attributes["taxid"], "1356405")
        self.assertEqual(len(protein), 1)
        self.assertEqual(protein["CDSList"].tag, "CDSList")
        self.assertEqual(protein["CDSList"].attributes, {})
        self.assertEqual(len(protein["CDSList"]), 2)
        self.assertEqual(protein["CDSList"][0], "")
        self.assertEqual(protein["CDSList"][0].attributes["kingdom"], "Bacteria")
        self.assertEqual(
            protein["CDSList"][0].attributes["assembly"], "GCA_000963615.1"
        )
        self.assertEqual(protein["CDSList"][0].attributes["start"], "264437")
        self.assertEqual(protein["CDSList"][0].attributes["stop"], "265978")
        self.assertEqual(protein["CDSList"][0].attributes["taxid"], "1356405")
        self.assertEqual(protein["CDSList"][0].attributes["strain"], "PML026")
        self.assertEqual(
            protein["CDSList"][0].attributes["org"], "Rhodococcus sp. PML026"
        )
        self.assertEqual(protein["CDSList"][0].attributes["kingdom_taxid"], "2")
        self.assertEqual(protein["CDSList"][0].attributes["accver"], "JZIS01000004.1")
        self.assertEqual(protein["CDSList"][0].attributes["strand"], "-")
        self.assertEqual(protein["CDSList"][1], "")
        self.assertEqual(protein["CDSList"][1].attributes["kingdom"], "Bacteria")
        self.assertEqual(
            protein["CDSList"][1].attributes["assembly"], "GCA_000963615.1"
        )
        self.assertEqual(protein["CDSList"][1].attributes["start"], "264437")
        self.assertEqual(protein["CDSList"][1].attributes["stop"], "265978")
        self.assertEqual(protein["CDSList"][1].attributes["taxid"], "1356405")
        self.assertEqual(protein["CDSList"][1].attributes["strain"], "PML026")
        self.assertEqual(
            protein["CDSList"][1].attributes["org"], "Rhodococcus sp. PML026"
        )
        self.assertEqual(protein["CDSList"][1].attributes["kingdom_taxid"], "2")
        self.assertEqual(protein["CDSList"][1].attributes["accver"], "KQ031368.1")
        self.assertEqual(protein["CDSList"][1].attributes["strand"], "-")
        protein = record["ProteinList"][1]
        self.assertEqual(protein.attributes["accver"], "WP_045840896.1")
        self.assertEqual(protein.attributes["source"], "RefSeq")
        self.assertEqual(
            protein.attributes["name"],
            "methylmalonate-semialdehyde dehydrogenase (CoA acylating)",
        )
        self.assertEqual(protein.attributes["taxid"], "1356405")
        self.assertEqual(protein.attributes["org"], "Rhodococcus sp. PML026")
        self.assertEqual(protein.attributes["kingdom_taxid"], "2")
        self.assertEqual(protein.attributes["kingdom"], "Bacteria")
        self.assertEqual(protein.attributes["priority"], "1")
        self.assertEqual(len(protein), 1)
        self.assertEqual(protein["CDSList"].tag, "CDSList")
        self.assertEqual(protein["CDSList"].attributes, {})
        self.assertEqual(len(protein["CDSList"]), 1)
        self.assertEqual(
            protein["CDSList"][0].attributes["assembly"], "GCF_000963615.1"
        )
        self.assertEqual(protein["CDSList"][0].attributes["start"], "264437")
        self.assertEqual(protein["CDSList"][0].attributes["stop"], "265978")
        self.assertEqual(protein["CDSList"][0].attributes["taxid"], "1356405")
        self.assertEqual(protein["CDSList"][0].attributes["strain"], "PML026")
        self.assertEqual(
            protein["CDSList"][0].attributes["org"], "Rhodococcus sp. PML026"
        )
        self.assertEqual(protein["CDSList"][0].attributes["kingdom_taxid"], "2")
        self.assertEqual(protein["CDSList"][0].attributes["accver"], "NZ_KQ031368.1")
        self.assertEqual(protein["CDSList"][0].attributes["strand"], "-")
        self.assertEqual(record["Statistics"], "")
        self.assertEqual(record["Statistics"].attributes["assmb_count"], "2")
        self.assertEqual(record["Statistics"].attributes["nuc_count"], "3")
        self.assertEqual(record["Statistics"].attributes["prot_count"], "2")

    def test_genbank(self):
        """Test error handling when presented with GenBank non-XML data."""
        # Access the nucleotide database using efetch, but return the data
        # in GenBank format.
        # To create the GenBank file, use
        # >>> Bio.Entrez.efetch(db='nucleotide', id='NT_019265', rettype='gb')
        from Bio.Entrez import Parser

        with open("GenBank/NT_019265.gb", "rb") as stream:
            self.assertRaises(Parser.NotXMLError, Entrez.read, stream)
        with open("GenBank/NT_019265.gb", "rb") as stream:
            iterator = Entrez.parse(stream)
            self.assertRaises(Parser.NotXMLError, next, iterator)

    def test_fasta(self):
        """Test error handling when presented with Fasta non-XML data."""
        from Bio.Entrez import Parser

        with open("Fasta/wisteria.nu", "rb") as stream:
            self.assertRaises(Parser.NotXMLError, Entrez.read, stream)
        with open("Fasta/wisteria.nu", "rb") as stream:
            iterator = Entrez.parse(stream)
            self.assertRaises(Parser.NotXMLError, next, iterator)

    def test_pubmed_html(self):
        """Test error handling when presented with HTML (so XML-like) data."""
        # To create the HTML file, use
        # >>> Bio.Entrez.efetch(db="pubmed", id="19304878")
        from Bio.Entrez import Parser

        with open("Entrez/pubmed3.html", "rb") as stream:
            self.assertRaises(Parser.NotXMLError, Entrez.read, stream)
        # Test if the error is also raised with Entrez.parse
        with open("Entrez/pubmed3.html", "rb") as stream:
            records = Entrez.parse(stream)
            self.assertRaises(Parser.NotXMLError, next, records)

    def test_xml_without_declaration(self):
        """Test error handling for a missing XML declaration."""
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="journals",id="2830,6011,7473",retmode='xml')
        from Bio.Entrez import Parser

        with open("Entrez/journals.xml", "rb") as stream:
            self.assertRaises(Parser.NotXMLError, Entrez.read, stream)
        # Test if the error is also raised with Entrez.parse
        with open("Entrez/journals.xml", "rb") as stream:
            records = Entrez.parse(stream)
            self.assertRaises(Parser.NotXMLError, next, records)

    def test_xml_without_definition(self):
        """Test error handling for a missing DTD or XML Schema."""
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="biosample", id="3502652", rettype="xml")
        with open("Entrez/biosample.xml", "rb") as stream:
            self.assertRaises(ValueError, Entrez.read, stream)
        # Test if the error is also raised with Entrez.parse
        with open("Entrez/biosample.xml", "rb") as stream:
            records = Entrez.parse(stream)
            self.assertRaises(ValueError, next, records)

    def test_truncated_xml(self):
        """Test error handling for a truncated XML declaration."""
        from io import BytesIO

        from Bio.Entrez.Parser import CorruptedXMLError

        truncated_xml = b"""<?xml version="1.0" encoding="UTF-8"  ?>
<!DOCTYPE GBSet PUBLIC "-//NCBI//NCBI GBSeq/EN" "https://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd">
<GBSet>
  <GBSeq>

    <GBSeq_locus>
        """
        stream = BytesIO()
        stream.write(truncated_xml)
        stream.seek(0)
        records = Entrez.parse(stream)
        self.assertRaises(CorruptedXMLError, next, records)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
