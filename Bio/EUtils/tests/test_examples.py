# These are burn tests -- I mostly just check if everything works

import unittest, sys, urllib2, time

from Bio import EUtils
from Bio.EUtils import ThinClient, HistoryClient, DBIdsClient, Datatypes
import support_test

picklestore = None
PICKLESTORE_FILENAME = "test_examples.pickle"

class TestForNewHits(unittest.TestCase):
    def testForNewItems(self):
        # See if there are Halobacterium salinarum hits

        known_dbids = EUtils.DBIds('protein', [
#            '461608', '133739', '121859', '121858', '114808',
            '26399547', '21362599', '20807968', '18310162',
            '10954592', '7674074', '3913880', '3913879', '3913877',
            '1350915', '462366', '462365', '461784', '119171',
            '21617825', '18202987', '14423949', '13878670',
            '13431821', '13124532', '12230579', '7531150', '6094165',
            '1352354', '1350833', '1169376', '140149', '133039',
            '133038', '132873', '132841', '114547', '114516',
            '18138454', '11385314', '11385312', '11385311',
            '11383975', '11383843', '11383705', '11383685',
            '11383631', '11383600', '11383544', '11383521',
            '11383449', '11383413', '11362758', '11362757',
            '11362756', '11362633', '11361208', '11361204',
            '11361203', '11361202', '11360968', '11360864',
            '11360863', '11360531', '11291672', '11290031',
            '11283208', '11281612', '11280741', '11280740',
            '11280737', '11280563', '11280281', '11280223',
            '11279456', '11277771', '11277752', '11277727',
            '11277613', '11277612', '11277611', '11277609',
            '11277608', '11277607', '11277606', '11277605',
            '11277604', '11277603', '11277602', '11277601',
            '11277597', '11277160', '11277109', '11277100',
            '11271301', '11270360', '11261072', '11260624',
            '11255573', '9972746', '7443396', '7440307', '7427602',
            '2129406', '2129405', '2129404', '2117914', '1363464',
            '1363463', '1084274', '1076149', '1076148', '1076147',
            '629399', '629398', '629396', '629394', '487075',
            '487074', '486684', '477959', '477726', '477027',
            '421715', '421714', '282664', '282663', '282662',
            '282661', '282660', '282659', '282658', '282657',
            '282656', '282655', '282654', '282653', '282652',
            '281164', '281162', '281161', '281160', '281159',
            '281158', '281157', '280348', '99214', '99212', '99211',
            '99210', '99209', '99207', '99206', '99205', '99203',
            '99202', '99201', '99200', '99198', '99195', '99192',
            '99191', '99190', '99188', '99187', '99185', '99183',
            '81076', '81071', '81070', '81069', '81067', '81065',
            '81063', '81062', '81061', '81060', '81058', '81057',
            '81055', '81054', '81052', '81051', '81047', '81045',
            '81043', '81042', '81041', '81040', '81039', '81038',
            '81037', '81036', '81035', '81034', '81032', '81031',
            '81029', '81026', '81023', '81022', '81020', '81019',
            '81018', '81017', '81016', '81015', '81013', '81012',
            '81011', '76372', '76329', '76324', '72632', '71284',
            '71211', '71179', '71146', '71063', '66382', '21264507',
            '18000395', '7531067', '132751', '117568', '14916722',
            '133992', '3122880', '114811', '24158915', '24158914',
            '24158913', '10640268', '3913878', '3023997', '544287',
            '19909691', '19909603', '3183098', '141254', '2492920',
            '809698', '43513', '18144841', '4467437', '4467436',
            '6226495', '3122803', '3122663', '2499389', '2494605',
            '1350834', '548760', '141356', '141353', '140804',
            '139955', '134663', '134662', '133462', '133078',
            '132785', '132646', '120245', '120232', '114812',
            '20873485', '20873480', '17942995', '17942994',
            '17942993', '20516540', '20151159', '20150922',
            '20150921', '12644298', '12644014', '120249', '10120917',
            '13124533', '133398', '11132078', '6225334', '1168858',
            '549772', '140827', '140825', '120239', '114838',
            '16974947', '16974946', '3219755', '347255', '8569313',
            '4929918', '899266', '10121032', '1621047', '14278674',
            '14278519', '11071687', '11071686', '11071685',
            '11071684', '6729723', '6729722', '6435624', '6435623',
            '3745775', '12655906', '12655904', '12655903', '6435626',
            '6435625', '11992133', '8918496', '10121037', '10121036',
            '10121033', '10120892', '10120851', '10120850', '2851428',
            '1710568', '1710498', '548772', '548769', '133031',
            '132995', '132985', '132879', '132784', '132750',
            '132728', '132708', '132639', '120251', '7327959',
            '7107278', '6435593', '455306', '455305', '455304',
            '455303', '455302', '455301', '455300', '455299',
            '455298', '455297', '150418', '150417', '6435594',
            '6172231', '6172230', '6172229', '6172228', '6172227',
            '6137453', '130153', '118549', '2499384', '2499383',
            '461612', '461611', '461610', '5822280', '285806',
            '4930169', '2425186', '2425185', '2425184', '2425183',
            '2425182', '2425181', '2425180', '2425179', '2425178',
            '2425177', '2425176', '2425175', '2425174', '3659953',
            '3659944', '4469246', '4001706', '4001704', '4388967',
            '4378986', '4377598', '809702', '809701', '43668',
            '43667', '43666', '43665', '43663', '43662', '43661',
            '43660', '43658', '43491', '4322492', '4104487',
            '4104485', '4104483', '3928158', '43508', '1633466',
            '2072795', '1527138', '1527137', '1353676', '3015619',
            '598123', '515085', '493889', '493888', '229726',
            '2351849', '2351848', '1154790', '1154789', '1154788',
            '1154787', '1154786', '1154785', '1154784', '1154783',
            '1154782', '1154781', '1154780', '1154779', '1154778',
            '1154777', '2190417', '2190416', '984742', '984741',
            '984740', '984739', '2760612', '2648028', '1070346',
            '1070345', '1070344', '2209068', '1487875', '1199752',
            '1199750', '509675', '509674', '285817', '285810',
            '285808', '216709', '43455', '43454', '43453', '43452',
            '43451', '43450', '807110', '1654427', '1654425',
            '1654423', '1654421', '1654419', '225948', '1583108',
            '1583107', '1094422', '448273', '226716', '226715',
            '350080', '226309', '225904', '225428', '1483625',
            '1235894', '1435134', '1435132', '1435131', '1435129',
            '994803', '994802', '226310', '223370', '223077',
            '223076', '223063', '1333716', '671101', '671100',
            '550341', '311841', '297410', '49046', '43656', '43655',
            '43654', '43653', '43641', '43640', '43559', '43557',
            '43552', '43551', '43550', '43548', '43546', '43545',
            '43544', '43543', '43542', '43541', '43540', '43539',
            '43537', '43535', '43534', '43533', '43531', '43530',
            '43526', '43525', '43524', '43523', '43522', '43521',
            '43520', '43519', '43518', '43517', '43511', '43510',
            '43505', '43504', '43503', '43501', '43499', '43498',
            '43496', '43495', '43493', '517390', '305353', '305352',
            '148816', '148814', '148812', '148794', '148793',
            '148792', '148768', '148767', '148766', '148764',
            '148763', '148759', '148757', '148753', '148749',
            '148747', '148745', '305350'])

        # Upload to the server
        client = HistoryClient.HistoryClient(
            eutils = picklestore.client())

        old_records = client.post(known_dbids)

        # Now see if there's anything new
        new_records = client.search("Halobacterium salinarum BUTNOT #%s" %
                                             (old_records.query_key,),
                                    db = "protein")
        print
        print "There are", len(new_records), "new Halobacterium salinarum records"
        new_dbids = new_records.dbids
        print "The record identifiers are:", ", ".join(map(str, new_dbids))
        print new_records.efetch(retmode = "text", rettype = "summary").read()

        # These should exist, since I commented them out above :)
        for x in ('461608', '133739', '121859', '121858', '114808'):
            assert x in new_dbids, "Cannot find expected record %r" % x


class TestRecentItems(unittest.TestCase):
    def testRecentItems(self):
        # Get publications mentioning 'rhodopsin' which were added to
        # PubMed in the last 10 days
        client = DBIdsClient.DBIdsClient(
            eutils = picklestore.client())
        results = client.search("rhodopsin", daterange = EUtils.WithinNDays(10))

        print
        print "There are", len(results), "rhodopsin records from the last 10 days"
        print "The record identifiers are:", ", ".join(map(str, results.dbids))
        print results.efetch(retmode = "text", rettype = "docsum").read()

class TestRelatedItems(unittest.TestCase):
    def testRelatedItems(self):
        # Get proteins similar to 4579714 (bacteriorhodopsin) which
        # were published in 2002
        client = DBIdsClient.DBIdsClient(
            eutils = picklestore.client())
        results = client.from_dbids(EUtils.DBIds("protein", "4579714"))

        neighbors = results.neighbor_links("protein",
                 daterange = EUtils.DateRange("2002/01/01", "2002/12/31", "pdat"))
        
        dbids = neighbors.linksetdbs["protein_protein"].dbids

        print
        print len(dbids), "sequences similar to GI:4579714 were published in 2002"
        print "The record identifiers are:", ", ".join(map(str, dbids))
        print client.from_dbids(dbids).efetch(
                     retmode = "text", rettype = "summary").read()

    def testRelatedSequences(self):
        # Get all protein sequences similar to 4579714 (bacteriorhodopsin)
        # in FASTA format
        client = HistoryClient.HistoryClient(
            eutils = picklestore.client())

        result = client.post(EUtils.DBIds("protein", "4579714"))
        related = result.neighbor_links("protein")
        related_dbids = related.linksetdbs["protein_protein"].dbids

        proteins = client.post(related_dbids)
        fasta_infile = proteins.efetch(retmode = "text", rettype = "fasta")
        fasta = fasta_infile.read()
        print fasta[:1000]

# http://www.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html
class TestEUtilsESearchExamples(unittest.TestCase):
    def setup(self):
        eutils = picklestore.client(1, "ESearch")
        self.dbids_client = DBIdsClient.DBIdsClient(eutils)
        self.history_client = DBIdsClient.DBIdsClient(eutils)
            
    def test1(self):
        self.setup()
        # Search in PubMed for the term cancer for the entrez date from the
        # last 60 days and retrieve the first 100 IDs and translations using
        # the history parameter:
        # http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?
        #   db=pubmed&term=cancer&reldate=60&datetype=edat&retmax=100&usehistory=y
        result = self.history_client.search("cancer",
                                            db = EUtils.databases.PUBMED,
                                            daterange = EUtils.WithinNDays(60, "edat"))
        dbids = result[:100].dbids

    def test2(self):
        self.setup()
        # Search in PubMed for the journal PNAS Volume 97, and retrieve 6 IDs
        # starting at ID 7
        # http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?
        #  db=pubmed&term=PNAS[ta]+AND+97[vi]&retstart=6&retmax=6
        result = self.dbids_client.search("PNAS[ta] AND 97[vi]",
                                          db = "pubmed",
                                          retstart = 6,
                                          retmax = 6)
        dbids = result.dbids

    def test3(self):
        self.setup()
        # Search in Journals for the term obstetrics:
        # http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?
        #   db=journals&term=obstetrics
        result = self.history_client.search("obstetrics", db = "journals")
        n = len(result)


# http://www.ncbi.nlm.nih.gov/entrez/query/static/elink_help.html
class TestEUtilsELinkExamples(unittest.TestCase):
    def setup(self):
        self.client = DBIdsClient.DBIdsClient(
            picklestore.client(1, "ELink"))

    def test1(self):
        self.setup()
        # To retrieve IDs and relvancy scores from pubmed for PMID 9298984
        # to the pubmed database:
        #   dbfrom=pubmed&id=9298984&cmd=neighbor
        results = self.client.from_dbids(EUtils.DBIds(
                                    "pubmed", "9298984")).neighbor_links()

    def test2(self):
        self.setup()
        # To retrieve IDs from nucleotide for GI 18250303, 18250301, 18250299
        # to protein:
        #    dbfrom=nucleotide&db=protein&id=18250303,18250307
        records = self.client.from_dbids(EUtils.DBIds("nucleotide",
                        ["18250303", "18250301", "18250299"]))
        neighbors = records.neighbor_links("protein")

    def test3(self):
        self.setup()
        # To retrieve pubmed related articles for PMIDs 11812492 11774222
        # with a publication date from 1995 to the present:
        #   dbfrom=pubmed&id=11812492,11774222&db=pubmed&mindate=1995&datetype=pdat
        # NOTE: this example is wrong -- it need the current date too!
        records = self.client.from_dbids(EUtils.DBIds("pubmed",
                                                      ["11812492", "11774222"]))
        # Can't actually use this since that would cause the regression
        # tests to fail.
        today = time.strftime("%Y/%m/%d")
        today = "2003/01/12"
        related = records.neighbor_links(
            daterange = EUtils.DateRange("1995", today, "pdat"))
        
    def test4(self):
        self.setup()
        # To retrieve MEDLINE indexed only related articles for PMID 12242737:
        #  dbfrom=pubmed&id=12242737&db=pubmed&term=medline[sb]
        records = self.client.from_dbids(EUtils.DBIds("pubmed", "12242737"))
        related = records.neighbor_links(term = "medline[sb]")

    def test5(self):
        self.setup()
        # To create a hyperlink to the first link available for PMID 10611131
        # in pubmed (NOTE: the 'prlinks' is only useful for a web browser
        # so I'm ignoring that part.)
        #  dbfrom=pubmed&id=10611131cmd=prlinks
        prlinks = self.client.from_dbids(EUtils.DBIds("pubmed", "10611131")).prlinks()

    def test6(self):
        self.setup()
        # To list all available links in pubmed for PMIDs 12085856 and 12085853:
        #   dbfrom=pubmed&id=12085856,12085853&cmd=llinks
        records = self.client.from_dbids(EUtils.DBIds("pubmed",
                                                      ["10611131", "12085853"]))
        llinks = records.llinks()

    def test7(self):
        self.setup()
        # To check for the existence of a Related Articles link for PMIDs
        # 0611131, 111645 and 12068369:
        #  dbfrom=pubmed&id=10611131+111645&id=12068369&cmd=ncheck
        records = self.client.from_dbids(EUtils.DBIds("pubmed",
                                             ["0611131", "111645", "12068369"]))
        check = records.ncheck()
 


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
