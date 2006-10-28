import unittest, sys, urllib2

from Bio import EUtils
from Bio.EUtils import ThinClient, HistoryClient, Datatypes, MultiDict
import support_test

picklestore = None
PICKLESTORE_FILENAME = "test_HistoryClient.pickle"
    
class TestBasicHistoryClient(unittest.TestCase):
    # Do a few queries and make sure things come out right
    def testHistory(self):
        entrez = HistoryClient.HistoryClient(
            eutils = picklestore.client())
        
        results1 = entrez.search("Dalke", field = "au",
                                 daterange = EUtils.DateRange("1995", "1998"))
        self.assertEquals(len(results1), 10)
        sizes = []
        expression = results1.metadata.expression
        for x in expression:
            if isinstance(x, Datatypes.Term):
                n = x.count
                assert n, n  # cannot be 0 or None
                sizes.append(n)
        self.assertEquals(len(sizes), 3)
        if sizes[0] < 30:
            raise AssertionError(sizes)
        self.assertEquals(sizes[1], -1)
        self.assertEquals(sizes[2], -1)

        self.assertEquals(expression.left.term, "Dalke[Author]")
        self.assertEquals(expression.left.field, "Author")

        self.assertEquals(expression.right.left.term, "1995[EDAT]")
        self.assertEquals(expression.right.left.field, "EDAT")
        self.assertEquals(expression.right.right.term, "1998[EDAT]")
        self.assertEquals(expression.right.right.field, "EDAT")


        expected_dbids = Datatypes.DBIds("pubmed", [
            "9454215", "9454196", "9454186", "9390282", "9303476",
            "9300720", "8763495", "8744570", "8566008", "7648552"])
        self.assertEquals(results1.dbids, expected_dbids)

        # this is a no-no, since EDAT isn't a searchable field
        self.failUnlessRaises(EUtils.EUtilsSearchError,
                              entrez.search,
                              "poliovirus AND 1995:1998[EDAT]",
                              db = "nucleotide")

        results2 = entrez.search("poliovirus AND 1995:1998[PDAT]",
                                 db = "pubmed")

        if len(results2) < 1160:
            raise AssertionError(len(results2))

        all_ids = results2.dbids
        self.assertEquals(len(all_ids), len(results2))
        self.assertEquals(all_ids[:20], results2[:20].dbids)
        self.assertEquals(all_ids[5:20], results2[5:20].dbids)
        self.assertEquals(all_ids[-5:], results2[-5:].dbids)
        self.assertEquals(all_ids[-5:-1], results2[-5:-1].dbids)
        self.assertEquals(all_ids[10:-14], results2[10:-14].dbids)

        # This is illegal because pubmed isn't a sequence database
        self.failUnlessRaises(TypeError, results2.efetch, seq_start = 0)

        # Try a different database
        results3 = entrez.search("poliovirus AND 1995:1998[PDAT]",
                                 db = "nucleotide")

        # Make sure I can still access fields from the first database
        self.assertEquals(results1.dbids, expected_dbids)

        # This is illegal because it mixes databases
        self.failUnlessRaises(EUtils.EUtilsSearchError,
                              entrez.search,
                              "#%s OR #%s" % (results1.query_key,
                                              results3.query_key))


        # However, this should yield the same as results3
        results4 = entrez.search("poliovirus", db = "nucleotide")
        results5 = entrez.search("#%s AND 1995:1998[PDAT]" % results4.query_key,
                                 db = results4.db)
        self.assertEquals(len(results3), len(results5))
        results3_dbids = results3.dbids
        self.assertEquals(results3_dbids, results5.dbids)

        # Get the sequence as FASTA one way
        s = results3[0].efetch(retmode = 'text', rettype = 'fasta').read()
        # And another way
        t = entrez.eutils.efetch_using_dbids(results3_dbids[:1], retmode = 'text',
                                             rettype = 'fasta').read()
        self.assertEquals(s, t)

class TestPostHistoryClient(unittest.TestCase):
    def testPostHistory(self):
        entrez = HistoryClient.HistoryClient(
            eutils = picklestore.client())
        # Upload a couple nucleotide record identifiers
        dbids = Datatypes.DBIds("nucleotide",
                                ["26285514", "16445404"])
        results = entrez.post(dbids)

        self.assertEquals(dbids, results.dbids)

        # Restrict to those which also mention "dopamine"
        results2 = entrez.search("#%s AND dopamine" % results.query_key,
                                 db = results.db)
        # Should be one hit
        self.assertEquals(Datatypes.DBIds("nucleotide", ["16445404"]),
                          results2.dbids)

        summary = results2[0].summary()
        self.assertEquals(summary.id, "16445404")
        gold_data = {
            "16445404": (
              ("Caption", "NM_000794"),
              ("Title", "Homo sapiens dopamine receptor D1 (DRD1), mRNA"),
              ("Extra", "gi|16445404|ref|NM_000794.2|"),
              ("Gi", 16445404),
              ("CreateDate", "1999/03/19"),
              ("UpdateDate", "2002/11/05"),
              ("Flags", 0),
              ("TaxId", 9606)),
            "26285514": (
              ("Caption", "AY165035"),
              ("Title", "Saccharomyces cerevisiae scavenger mRNA decapping enzyme (DCS1) mRNA, complete cds"),
              ("Extra", "gi|26285514|gb|AY165035.1|"),
              ("Gi", 26285514),
              ("CreateDate", "2002/12/10"),
              ("UpdateDate", "2002/12/10"),
              ("Flags", 0),
              ("TaxId", 4932)),
            }
        for fieldname, value in gold_data[summary.id]:
            self.assertEquals(summary.dataitems[fieldname], value)

        for x in results2.summary():
            for fieldname, value in gold_data[x.id]:
                self.assertEquals(x.dataitems[fieldname], value)

def _trim_lines(lines):
    while lines[-1] in ("", "\n"):
        del lines[-1]

class TestSequenceFetch(unittest.TestCase):
    def testSeqFetch(self):
        entrez = HistoryClient.HistoryClient(
            eutils = picklestore.client())
        dbids = EUtils.DBIds("nucleotide", "8886082")
        result = entrez.post(dbids)
        lines = result.efetch(retmode = "text", rettype = "gb").readlines()
        _trim_lines(lines)
        self.assertEquals(lines[-2], """\
     2761 ttcctagata aacac
""")
        lines = result.efetch(retmode = "text", rettype = "gb",
                              seq_start = 1, seq_stop = 60).readlines()
        _trim_lines(lines)
        self.assertEquals(lines[-2], """\
        1 ttcagtttat ttggacggaa gaatggtggc tcaattattg acgacttgtg gcctaaattt
""")
        lines = result.efetch(retmode = "text", rettype = "gb",
                              seq_start = 3, seq_stop = 40,
                              strand = EUtils.PLUS_STRAND).readlines()
        _trim_lines(lines)
        self.assertEquals(lines[-2], """\
        1 cagtttattt ggacggaaga atggtggctc aattattg
""")
        lines = result.efetch(retmode = "text", rettype = "gb",
                              seq_start = 3, seq_stop = 40,
                              strand = EUtils.MINUS_STRAND).readlines()
        _trim_lines(lines)
        self.assertEquals(lines[-2], """\
        1 caataattga gccaccattc ttccgtccaa ataaactg
""")

class TestLinks(unittest.TestCase):
    def setup(self):
        self.entrez = HistoryClient.HistoryClient(
            eutils = picklestore.client(1))
        self.result1 = self.entrez.post(EUtils.DBIds("protein", "4579714"))
        self.result2 = self.entrez.post(EUtils.DBIds("nucleotide",
                                        ["18250303", "18250301", "18250299"]))
        
    def testNeighborLinks(self):
        self.setup()
        # Get the document neighbors in PubMed
        neighbors = self.result1.neighbor_links(EUtils.databases.PUBMED)
        self.assertEquals(neighbors.dbids, self.result1.dbids)

        expected = MultiDict.OrderedMultiDict(
            [("protein_pubmed", Datatypes.LinkSetDb("pubmed", "protein_pubmed",
                                                    [Datatypes.Link("9878396")]))])
        self.assertEquals(neighbors.linksetdbs, expected)
                                              
        # Get protein neighbors
        neighbors = self.result2.neighbor_links(EUtils.databases.PROTEIN)
        self.assertEquals(neighbors.dbids, self.result2.dbids)
        self.assertEquals(neighbors.linksetdbs,
                          MultiDict.OrderedMultiDict(
            [("nucleotide_protein", Datatypes.LinkSetDb("protein", "nucleotide_protein",
                                                        [Datatypes.Link("18250304")]))]
            ))
        
    def testNCheck(self):
        self.setup()
        ncheck = self.result1.ncheck(EUtils.databases.PUBMED)
        self.assertEquals(ncheck.dbfrom, "protein")
        self.assertEquals(ncheck.idchecks,
                          [Datatypes.IdCheck("4579714", 0, 1)])

        
        ncheck = self.result2.ncheck(EUtils.databases.PUBMED)
        self.assertEquals(ncheck.idchecks,
                          [Datatypes.IdCheck("18250303", 0, 1)])
        
    def testLCheck(self):
        self.setup()
        lcheck = self.result1.lcheck(EUtils.databases.PUBMED)
        self.assertEquals(lcheck.dbfrom, "protein")
        self.assertEquals(lcheck.idchecks,
                          [Datatypes.IdCheck("4579714", 1, 0)])

        lcheck = self.result2.lcheck(EUtils.databases.PUBMED)
        self.assertEquals(lcheck.idchecks,
                          [Datatypes.IdCheck("18250303", 1, 0)])

    def testLLinks(self):
        self.setup()
        self.assertEquals(self.result1.llinks(),
              Datatypes.LinksLinkSet("protein",
                  [Datatypes.IdUrlSet("4579714",
        [Datatypes.ObjUrl(["protein identification/characterization"],
                          Datatypes.Provider(
                                   "Domain Architecture Retrieval Tool",
                                   "DART", "3240", None, None),
                          "Domain Neighbors",
                          "http://www.ncbi.nlm.nih.gov/Structure/lexington/lexington.cgi?cmd=seq&FILTER=on&EXPECT=0.01&DATALIB=oasis_sap&fasta=4579714",
                          [])]
                                      )]))

        self.assertEquals(self.result2.llinks(),
          Datatypes.LinksLinkSet("nucleotide",
            [Datatypes.IdUrlSet("18250303",
              [Datatypes.ObjUrl(
                  [],
                  Datatypes.Provider(
                    "NCBI LocusLink",
                    "LocusLink", "3086",
                    "http://www.ncbi.nlm.nih.gov/LocusLink",
                    "http://www.ncbi.nlm.nih.gov/LocusLink/IMG/LLogo.gif"),
                    None,
                    "http://www.ncbi.nlm.nih.gov/LocusLink/list.cgi?V=0&Q=18250303[ngi]",
                    []),
               Datatypes.ObjUrl(
                  ["sequence screening/similarity/alignment"],
                  Datatypes.Provider(
                    "UCSC Genome Browser", "UCSCgb", "3715",
                    "http://www.cse.ucsc.edu/centers/cbe/Genome/",
                    "http://www.cse.ucsc.edu/cbe/images/ucsc.gif"),
                  None,
                 "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg11&position=NM_080831",
                  []),
               Datatypes.ObjUrl(
                  ["DNA/protein sequence"],
                  Datatypes.Provider(
                    "UniGene", "UniGene", "3170",
                    "http://www.ncbi.nlm.nih.gov/UniGene/",
                    "http://www.ncbi.nlm.nih.gov/UniGene/IMG/button.gif"),
                  None,
                  "http://www.ncbi.nlm.nih.gov/UniGene/query.cgi?TEXT=@gi(18250303)&ORG=Hs",
                  [])
               ]
                                )
            ]))
                          
    
    def testPRLinks(self):
        self.setup()
        self.assertEquals(self.result1.prlinks(),
                          Datatypes.LinksLinkSet('protein',
                              [Datatypes.IdUrlSet('4579714', [])]))

        self.assertEquals(self.result2.prlinks(),
                          Datatypes.LinksLinkSet('nucleotide',
                              [Datatypes.IdUrlSet('18250303', [])]))

        

def main():
    global picklestore

    picklestore = support_test.UsePickleStore(PICKLESTORE_FILENAME)
    
    if len(sys.argv) == 2:
        if sys.argv[1] == "--use-live":
            picklestore = support_test.NoPickleStore()
            del sys.argv[1]
        elif sys.argv[1] == "--create-pickle":
            picklestore = support_test.CreatePickleStore(PICKLESTORE_FILENAME)
            del sys.argv[1]
        elif sys.argv[1] == "--use-pickle":
            picklestore = support_test.UsePickleStore(PICKLESTORE_FILENAME)
            del sys.argv[1]
        elif sys.argv[1] == "--help":
            print """\
Usage: test_HistoryClient.py [--use-live|--create-pickle|--use-pickle|--help]

If you specify nothing, this runs the unittest against the golden pickle file.
"""
            sys.exit(0)
    try:
        unittest.main()
    finally:
        picklestore.done()
        
if __name__ == "__main__":
    #ThinClient.DUMP_URL = 1
    main()
