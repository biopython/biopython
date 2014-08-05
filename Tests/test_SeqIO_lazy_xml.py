import unittest
import os
import sys
from io import BytesIO

from Bio.SeqIO import _lazy, UniprotIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio import SeqIO
from Bio._py3k import _string_to_bytes

#when possible, use 'zip_longest'
if sys.version[0] == "3":
    from itertools import zip_longest
    compatible_zip = zip_longest
elif sys.version[0] == "2":
    compatible_zip = zip

SIMPLEXML = ['<begin>\n  <acount n="3" />\n"',
             '<a>\n    <b number="1">\n      ',
             '<c>\n        text for c\n      </c>\n    </b>\n  </a>',
             '\n  ',
             '<a>\n    <b number="2">\n      <c>\n        textc1\n      ',
             '</c>\n      <c>\n        textc2\n      </c>\n    </b>\n  </a>',
             '\n</begin>\n']

# No newlines
CONDENSEDXML = ['<begin><acount n="3">',
                '<a><b number="1">',
                '<c>text for c</c></b></a>',
                '',
                '<a><b number="2"><c>textc1',
                '</c><c>textc2</c></b></a>',
                '</begin>']
# Erratic newline use and somewhat nonstandard whitespace
WEIRDXML = ['<begin>\n <acount   n="3"/>\n ',
            '<a>\n <b number="1">\n ',
            '<c>\n text for c\n </c>\n </b>\n </a>',
            '\n ',
            '<a>\n <b number="2">\n<c>\n textc1\n ',
            '</c>\n <c >\n textc2\n </c>          </b>\n </a>',
            '\n</begin>\n']

def xml_parser_iter(ioobject, targetfield, tagstoparse):
    """Use ExpatHandler to iterate through all 'a' records

    In addition to allowing simplified tests, this function provides
    a nice tutorial on using the ExpatHandler to find all records
    of a given type (in this case, 'a') in an XML file.
    """
    position = 0
    parser = _lazy.ExpatHandler(ioobject, targetfield, tagstoparse)
    while True:
        root = parser.parse_from_position(position)
        yield root
        if root.lastrecord is True:
            break
        position = root.nextelementoffset


class XmlIndexerTests(unittest.TestCase):

    xmllist = SIMPLEXML

    def setUp(self):
        self.ioobject = BytesIO(_string_to_bytes(''.join(self.xmllist)))
        targetfield = "a"
        tagstoparse = ["a", "b", "c"]
        self.resultparser = xml_parser_iter(self.ioobject, \
                                targetfield, tagstoparse)

    def test_iter_finds_2_roots(self):
        count = 0
        for root in self.resultparser:
            count += 1
            self.assertEqual(root.tag, "ROOT")
        self.assertEqual(2, count)

    def test_each_root_contains_correct_tree(self):
        for root in self.resultparser:
            #test children of root

            a = root.children[0]
            rootchildrenlen = len(root.children)
            self.assertEqual(a.tag, 'a')
            self.assertEqual(rootchildrenlen, 1)
            #test children of a; expected a single 'b' with
            b = a.children[0]
            achildrenlen = len(root.children)
            self.assertEqual(b.tag, 'b')
            #testing b attributes
            self.assertTrue("number" in b.attributes)
            self.assertEqual(len(b.attributes.keys()), 1)
            self.assertEqual(achildrenlen, 1)
            #test children of b, count variabl
            ctags = b.children
            bchildrenlen = len(b.children)
            correctblen = int(b.attributes["number"])
            self.assertTrue(any(map(lambda c: c.tag == "c", ctags)))
            self.assertEqual(bchildrenlen, correctblen)

    def test_root_extracts_tag(self):
        root1 = next(self.resultparser)
        extractedtext = repr(root1.extract_from_handle(self.ioobject))
        knowntext = repr("".join(self.xmllist[1:3]))
        self.assertEqual(extractedtext, knowntext)
        root2 = next(self.resultparser)
        extractedtext = repr(root2.extract_from_handle(self.ioobject))
        knowntext = repr("".join(self.xmllist[4:6]))
        self.assertEqual(extractedtext, knowntext)

class XmlIndexerTestsCONDENSEDXML(XmlIndexerTests):
    xmllist = CONDENSEDXML

class XmlIndexerTestsErraticFormat(XmlIndexerTests):
    xmllist = WEIRDXML

class UniProtXmlTest(unittest.TestCase):

    ioobject = open("SwissProt/uni001", 'rb')

    def setUp(self):
        targetfield = "entry"
        tagstoparse = ["entry", "accession", "feature", "sequence"]
        self.resultparser = xml_parser_iter(self.ioobject, \
                                targetfield, tagstoparse)

    def test_iter_finds_roots(self):
        count = 0
        for root in self.resultparser:
            count += 1
            self.assertEqual(root.tag, "ROOT")
            self.assertTrue(len(root.children) == 1)
            self.assertTrue(len(root.children[0].children) >= 2)
        self.assertTrue(count >= 1)

class TestUniprot_copy_delete_this(unittest.TestCase):

    def test_uni001(self):
        "Parsing Uniprot file uni001"
        filename = 'uni001'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "uniprot-xml")

        self.assertTrue(isinstance(seq_record, SeqRecord))

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(seq_record.id, "Q91G55")
        self.assertEqual(seq_record.name, "043L_IIV6")
        self.assertEqual(seq_record.description, "Uncharacterized protein 043L")
        self.assertEqual(repr(seq_record.seq), "Seq('MDLINNKLNIEIQKFCLDLEKKYNINYNNLIDLWFNKESTERLIKCEVNLENKI...IPI', ProteinAlphabet())")

        # self.assertEqual(seq_record.accessions, ['Q91G55']) #seq_record.accessions does not exist
        # self.assertEqual(seq_record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        # self.assertEqual(record.seqinfo, (348, 39676, '75818910'))

        self.assertEqual(len(seq_record.features), 1)
        self.assertEqual(repr(seq_record.features[0]), "SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(116)), type='chain', id='PRO_0000377969')")

        self.assertEqual(len(seq_record.annotations['references']), 2)
        self.assertEqual(seq_record.annotations['references'][0].authors, 'Jakob N.J., Mueller K., Bahr U., Darai G.')
        self.assertEqual(seq_record.annotations['references'][0].title, 'Analysis of the first complete DNA sequence of an invertebrate iridovirus: coding strategy of the genome of Chilo iridescent virus.')
        self.assertEqual(seq_record.annotations['references'][0].journal, 'Virology 286:182-196(2001)')
        self.assertEqual(seq_record.annotations['references'][0].comment, 'journal article | 2001 | Scope: NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA] | ')

        self.assertEqual(len(seq_record.dbxrefs), 11)
        self.assertEqual(seq_record.dbxrefs[0], 'DOI:10.1006/viro.2001.0963')

        self.assertEqual(seq_record.annotations['sequence_length'], 116)
        self.assertEqual(seq_record.annotations['sequence_checksum'], '4A29B35FB716523C')
        self.assertEqual(seq_record.annotations['modified'], '2009-07-07')
        self.assertEqual(seq_record.annotations['accessions'], ['Q91G55'])
        self.assertEqual(seq_record.annotations['taxonomy'], ['Viruses', 'dsDNA viruses, no RNA stage', 'Iridoviridae', 'Iridovirus'])
        self.assertEqual(seq_record.annotations['sequence_mass'], 13673)
        self.assertEqual(seq_record.annotations['dataset'], 'Swiss-Prot')
        self.assertEqual(seq_record.annotations['gene_name_ORF'], ['IIV6-043L'])
        self.assertEqual(seq_record.annotations['version'], 21)
        self.assertEqual(seq_record.annotations['sequence_modified'], '2001-12-01')
        self.assertEqual(seq_record.annotations['keywords'], ['Complete proteome', 'Virus reference strain'])
        self.assertEqual(seq_record.annotations['organism_host'], ['Acheta domesticus', 'House cricket', 'Chilo suppressalis', 'striped riceborer', 'Gryllus bimaculatus', 'Two-spotted cricket', 'Gryllus campestris', 'Spodoptera frugiperda', 'Fall armyworm'])
        self.assertEqual(seq_record.annotations['created'], '2009-06-16')
        self.assertEqual(seq_record.annotations['organism_name'], ['Chilo iridescent virus'])
        self.assertEqual(seq_record.annotations['organism'], 'Invertebrate iridescent virus 6 (IIV-6)')
        self.assertEqual(seq_record.annotations['recommendedName_fullName'], ['Uncharacterized protein 043L'])
        self.assertEqual(seq_record.annotations['sequence_version'], 1)
        self.assertEqual(seq_record.annotations['proteinExistence'], ['Predicted'])

class TestUniprotLazyComparison(unittest.TestCase):
    recordfile = "uni001"

    def setUp(self):
        self.filename = os.path.join('SwissProt', self.recordfile)
        self.handle = open(self.filename, 'rb')
        returncls = UniprotIO.UniprotXMLSeqRecProxy
        self.parser = SeqIO._lazy.LazyIterator(self.handle, returncls)
        self.oldparser = SeqIO.parse(self.filename, "uniprot-xml")
        self.handle.seek(0)

    def tearDown(self):
        self.handle.close()

    def test_parser_init(self):
        record = next(iter(self.parser))

    def test_id_name(self):
        for old, lazy in compatible_zip(self.oldparser, self.parser):
            self.assertEqual(old.id, lazy.id)
            self.assertEqual(old.name, lazy.name)

    def test_record_description_annotations(self):
        for old, lazy in compatible_zip(self.oldparser, self.parser):
            self.assertEqual(repr(old.description), repr(lazy.description))
            self.assertEqual(repr(old.annotations), repr(lazy.annotations))

    def test_id_seq(self):
        for old, lazy in compatible_zip(self.oldparser, self.parser):
            self.assertEqual(str(old[0:5].seq), str(lazy[0:5].seq))
            self.assertEqual(str(old[50:60].seq), str(lazy[50:60].seq))
            self.assertEqual(str(old.seq), str(lazy.seq))
            self.assertEqual(str(old[50:60].seq), str(lazy[50:60].seq))

    def test_lazy_record_has_features(self):
        for old, lazy in compatible_zip(self.oldparser, self.parser):
            #if the old parser has features, then we can check the lazy
            if len(old.features) > 0:
                self.assertTrue(isinstance(old.features[0], SeqFeature))
                self.assertTrue(isinstance(lazy.features[0], SeqFeature))

    def make_feature_tuple_list(self, featurelist):
        returnlist = [(f.location.nofuzzy_start,
                       f.location.nofuzzy_end,
                       f.type) for f in featurelist]
        returnlist.sort(key = lambda t: t[0]+(t[1]/100.0) + \
                          sum(ord(c) for c in t[2])/100000)
        return returnlist

    def test_record_has_same_set_of_features(self):
        for old, lazy in compatible_zip(self.oldparser, self.parser):
            self.assertEqual(len(old.features), len(lazy.features))
            #default features don't allow comparison
            oldset = self.make_feature_tuple_list(old.features)
            newset = self.make_feature_tuple_list(lazy.features)
            self.assertEqual(oldset, newset)

class TestUniprotLazyComparison002(TestUniprotLazyComparison):
    recordfile = "uni002"

class TestUniprotLazyComparison003(TestUniprotLazyComparison):
    recordfile = "uni003"

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
