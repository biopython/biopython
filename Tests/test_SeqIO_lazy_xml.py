import unittest
from io import BytesIO

from Bio.SeqIO import _lazy
from Bio._py3k import _string_to_bytes

simplexml = '<begin>\n  <acount n="3" />\n  <a>\n    <b number="1">\n      ' +\
            '<c>\n        text for c\n      </c>\n    </b>\n  </a>\n  ' +\
            '<a>\n    <b number="2">\n      <c>\n        textc1\n      ' +\
            '</c>\n      <c>\n        textc2\n      </c>\n    </b>\n  ' +\
            '</a>\n</begin>\n'
simplexmlio = BytesIO(_string_to_bytes(simplexml))
condensed = '<begin><acount n="3"> <a> <b number="1"> ' +\
            '<c> text for c </c> </b> </a> ' +\
            '<a> <b number="2"><c> textc1 ' +\
            '</c> <c> textc2 </c> </b> ' +\
            '</a></begin>'
condensedio = BytesIO(_string_to_bytes(condensed))
weirdxml = '<begin>\n <acount   n="3"/>\n <a>\n <b number="1">\n ' +\
           '<c>\n text for c\n </c>\n </b>\n </a>\n ' +\
           '<a>\n <b number="2">\n<c>\n textc1\n ' +\
           '</c>\n <c >\n textc2\n </c>          </b>\n ' +\
           '</a>\n</begin>\n'
weirdxmlio = BytesIO(_string_to_bytes(weirdxml))


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
        position = root.indexend


class XmlIndexerTests(unittest.TestCase):
    ioobject = simplexmlio
    def setUp(self):
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

class XmlIndexerTestsCondensed(XmlIndexerTests):
    ioobject = condensedio

class XmlIndexerTestsErraticFormat(XmlIndexerTests):
    ioobject = weirdxmlio

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

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)


