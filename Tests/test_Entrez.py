"""Testing code for Bio.Entrez parsers.
"""

import sys
import unittest

from Bio import Entrez

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [EInfoTest, ]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class EInfoTest(unittest.TestCase):
    """Tests for dealing with basic enzymes using the Restriction package.
    """
    def t_list(self):
        """Test parsing database list returned by EInfo
        """
        input = open("Entrez/einfo1.xml")
        record = Entrez.read(input)
        assert record==["pubmed",
                        "protein",
                        "nucleotide",
                        "nuccore",
                        "nucgss",
                        "nucest",
                        "structure",
                        "genome",
                        "books",
                        "cancerchromosomes",
                        "cdd",
                        "gap",
                        "domains",
                        "gene",
                        "genomeprj",
                        "gensat",
                        "geo",
                        "gds",
                        "homologene",
                        "journals",
                        "mesh",
                        "ncbisearch",
                        "nlmcatalog",
                        "omia",
                        "omim",
                        "pmc",
                        "popset",
                        "probe",
                        "proteinclusters",
                        "pcassay",
                        "pccompound",
                        "pcsubstance",
                        "snp",
                        "taxonomy",
                        "toolkit",
                        "unigene",
                        "unists"
                       ]
    def t_pubmed(self):
        """Test parsing database info returned by EInfo
        """
        input = open("Entrez/einfo2.xml")
        record = Entrez.read(input)
        assert record["DbName"]=="pubmed"
        assert record["MenuName"]=="PubMed"
        assert record["Description"]=="PubMed bibliographic record"
        assert record["Count"]==17905967
        assert record["LastUpdate"]=="2008/04/15 06:42"

        assert len(record["FieldList"])==40

        assert record["FieldList"][0]["Name"]=="ALL"
        assert record["FieldList"][0]["FullName"]=="All Fields"
        assert record["FieldList"][0]["Description"]=="All terms from all searchable fields"
        assert record["FieldList"][0]["TermCount"]==70792830
        assert record["FieldList"][0]["IsDate"]==False
        assert record["FieldList"][0]["IsNumerical"]==False
        assert record["FieldList"][0]["SingleToken"]==False
        assert record["FieldList"][0]["Hierarchy"]==False
        assert record["FieldList"][0]["IsHidden"]==False

        assert len(record["LinkList"])==46

        assert record["LinkList"][0]["Name"]=="pubmed_books_refs"
        assert record["LinkList"][0]["Menu"]=="Cited in Books"
        assert record["LinkList"][0]["Description"]=="PubMed links associated with Books"
        assert record["LinkList"][0]["DbTo"]=="books"

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
