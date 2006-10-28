import re, sys, urllib2, unittest, urllib

from Bio import EUtils
from Bio.EUtils import ThinClient

# This test uses a local file instead of going over the net.
# This helps because:
#  - the network/server doesn't need to be up
#  - it is a lot faster
#  - some fields on the server might change
#
# NOTE: this requires that the data always be accessed
#  in the same order.  You may need to regenerate the
#  reference ("golden") pickle data set if things change.
#
# To generate the golden pickle, use '--create-pickle'
# To test from the golden pickle, use '--use-picke'
# To test the server directly, use '--use-live'
# To run the standard unit test, don't use any parameters


import support_test

picklestore = None
PICKLESTORE_FILENAME = "test_ThinClient.pickle"



####################


# Not unittest friendly
class TestThinClient(unittest.TestCase):
    def testClient(self):
        eutils = picklestore.client()
        infile = eutils.esearch("Dalke", field="au",
                                daterange = EUtils.DateRange("1995", "1998"),
                                retstart = 1, retmax = 5,
                                usehistory = 1, webenv = None)
        s = infile.read()
        try:
            counts = map(int, re.findall(r"<Count>(-?\d+)</Count>", s))
            assert len(counts) == 4, counts
            assert counts[0] == 10
            assert counts[1] >= 30
            assert counts[2] == -1
            assert counts[3] == -1

            assert s.find("<RetMax>5</RetMax>") != -1
            assert s.find("<RetStart>1</RetStart>") != -1
            ids = re.findall(r"<Id>(\d+)</Id>", s)
            assert len(ids) == 5, ids

            terms = re.findall(r"<Term>([^<]+)</Term>", s)
            assert len(terms) == 3
            assert terms[0] == "Dalke[Author]"
            assert terms[1] == "1995[EDAT]"
            assert terms[2] == "1998[EDAT]"

            query_key1 = re.findall(r"<QueryKey>(\d+)</QueryKey>", s)[0]
            assert query_key1 == "1", query_key1  # always true?

            quoted_webenv = re.findall("<WebEnv>([^>]+)</WebEnv>", s)[0]
            webenv = urllib.unquote(quoted_webenv)

        except:
            print "ERROR!"
            print s
            raise

        try:
            # Can I refetch those same Ids using the history?
            t = ""
            t = eutils.efetch_using_history(
                db = "pubmed", webenv = webenv, query_key = query_key1,
                retstart = 1, retmax = 5,
                retmode = "text", rettype = "uilist").read()
            new_ids = t.split()
            assert ids == new_ids, (ids, new_ids)  # Must be in same order too!
        except:
            print "ERROR!"
            print s
            print " -- and --"
            print t
            raise

        # Make sure I'm getting the same XML summary through history and id
        sum1 = sum2 = None
        try:
            sum1 = eutils.esummary_using_history(
                db = "pubmed", webenv = webenv, query_key = query_key1,
                retstart = 1, retmax = 1).read()
            sum2 = eutils.esummary_using_dbids(
                dbids = EUtils.DBIds("pubmed", [ids[0]])).read()
            assert sum1 == sum2
        except:
            print "Summary 1"
            print sum1
            print "-----------------"
            print "Summary 2"
            print sum2
            raise

        # Make sure I'm getting the same XML version of the records
        rec1 = rec2 = None
        try:
            rec1 = eutils.efetch_using_history(
                db = "pubmed", webenv = webenv, query_key = query_key1,
                retmode = "xml", retstart = 1, retmax = 1).read()
            rec2 = eutils.efetch_using_dbids(
                dbids = EUtils.DBIds("pubmed", [ids[0]]),
                retmode = "xml").read()
            assert rec1 == rec2
        except:
            print "Record 1"
            print rec1
            print "-----------------"
            print "Record 2"
            print rec2
            raise


        # Post a few GIs (from the protein database) to the server
        # This appends to the existing history so should be query_key #2.
        post_ids = ["914034", "5263173", "1769808", "1060883"]
        infile = eutils.epost(EUtils.DBIds("protein", post_ids),
                              webenv = webenv)
        post_results = infile.read()
        try:
            query_key2 =  re.findall(r"<QueryKey>(\d+)</QueryKey>",
                                     post_results)[0]
            assert query_key2 == "2"

            quoted_webenv = re.findall("<WebEnv>([^>]+)</WebEnv>",
                                       post_results)[0]
            webenv = urllib.unquote(quoted_webenv)
        except:
            print "ERROR"
            print post_results
            raise

        # Verify that the posted ids are correct
        posted_ids = eutils.efetch_using_history(
            db = "pubmed", webenv = webenv, query_key = query_key2,
            retstart = 0, retmax = len(post_ids),
            retmode = "text", rettype = "uilist").read().split()
        x1 = posted_ids[:]  # Make copies since I need the correct
        x1.sort()           # order for getting the FASTA version, below
        x2 = post_ids[:]
        x2.sort()
        assert x1 == x2, (post_ids, posted_ids)

        # Now fetch them as FASTA format
        fasta1 = fasta2 = None
        try:
            fasta1 = eutils.efetch_using_history(
                db = "protein", webenv = webenv, query_key = query_key2,
                retstart = 0, retmax = len(post_ids),
                retmode = "text", rettype = "fasta").read()
            fasta2 = eutils.efetch_using_dbids(
                dbids = EUtils.DBIds("protein", posted_ids),
                retmode = "text", rettype = "fasta").read()
            assert fasta1 == fasta2
        except:
            print "ERROR FASTA1"
            print fasta1
            print "ERROR FASTA2"
            print fasta2
            raise


        # It's much harder to test the ELink capabilities.

        # Get the VMD paper
        results = None
        try:
            results = eutils.esearch(
                "Humphrey W. AND Dalke A. AND Schulten K. AND VMD[Title]",
                field = "au").read()
            # There should only be one match
            ids = re.findall(r"<Id>(\d+)</Id>", results)
            assert ids == ["8744570"]
        except:
            print "Error"
            print results
            raise

        # Look at the related publications and we should find
        # my Tcl paper, which is 9390282
        links = None
        try:
            links = eutils.elink_using_dbids(EUtils.DBIds("pubmed", ids),
                                             cmd = "neighbor").read()

            # remember, the first id comes from the <Id> in <IdList>
            related_ids = re.findall(r"<Id>(\d+)</Id>", links)[1:]
            assert "9390282" in related_ids
        except:
            print "Error"
            print links
            raise

        # Get the taxonomy record for the "posted_ids".
        # NOTE: This test original compared the 2nd element from
        # that list, but elink_using_history doesn't support the
        # retstart/retmax parameters.
        #
        # This comes from query_key2 in the history.  Do it both ways
        # to compare results.
        link1 = link2 = None
        try:
            link1 = eutils.elink_using_dbids(EUtils.DBIds("protein", posted_ids),
                                             db = "taxonomy",
                                             cmd = "neighbor").read()
            link2 = eutils.elink_using_history(
                dbfrom = "protein", webenv = webenv, query_key = query_key2,
                db = "taxonomy",
                cmd = "neighbor").read()
            assert link1 == link2
            taxids = re.findall(r"<Id>(\d+)</Id>", link1)[len(posted_ids):]
            assert taxids == ["43776", "29282", "28442", "2237"], taxids
        except:
            print "Error",
            print link1
            print "----------------"
            print link2
            raise

        # See if there are linkouts 914034
        # Should be at least one, to DART 3240.
        llinks = None
        try:
            llinks = eutils.elink_using_dbids(EUtils.DBIds("protein", [posted_ids[1]]),
                                              cmd = "llinks").read()
            assert llinks.find("<ObjUrl>") != -1
            assert "3240" in re.findall("<Id>(\d+)</Id>", llinks)
        except:
            print "ERROR"
            print llinks
            raise

        # Finally, check that I can limit the seach to an Entrez query string
        # I'm using the example
        #   "retrieve MEDLINE indexed only related articles for PMID 12242737"
        #   elink.fcgi?dbfrom=pubmed&id=12242737&db=pubmed&term=medline[sb]
        full = restricted = None
        try:
            full = eutils.elink_using_dbids(EUtils.DBIds("pubmed", ["12242737"])).read()
            restricted =  eutils.elink_using_dbids(EUtils.DBIds("pubmed", ["12242737"]),
                                                   term = "medline[sb]").read()
            counts1 = full.count("<Link>")
            counts2 = restricted.count("<Link>")
            assert counts1 > counts2
        except:
            print "ERROR"
            print full
            print "---------"
            print restricted
            raise

####################################3

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
Usage: %s [--use-live|--create-pickle|--use-pickle|--help]

If you specify nothing, this runs the unittest against the golden pickle file.
""" % (sys.argv[0],)
            sys.exit(0)
    try:
        unittest.main()
    finally:
        picklestore.done()
        
if __name__ == "__main__":
    #ThinClient.DUMP_URL = 1
    main()
