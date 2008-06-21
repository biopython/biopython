'''Testing code for Bio.Entrez parsers.
'''

import sys
import unittest

from Bio import Entrez

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    '''Generate the suite of tests.
    '''
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [EInfoTest, ESearchTest, EPostTest, ESummaryTest, ELinkTest, EFetchTest, EGQueryTest, ESpellTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class EInfoTest(unittest.TestCase):
    '''Tests for parsing XML output returned by EInfo
    '''
    def t_list(self):
        '''Test parsing database list returned by EInfo
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.einfo()
        input = open('Entrez/einfo1.xml')
        record = Entrez.read(input)
        assert record["DbList"]==['pubmed',
                                  'protein',
                                  'nucleotide',
                                  'nuccore',
                                  'nucgss',
                                  'nucest',
                                  'structure',
                                  'genome',
                                  'books',
                                  'cancerchromosomes',
                                  'cdd',
                                  'gap',
                                  'domains',
                                  'gene',
                                  'genomeprj',
                                  'gensat',
                                  'geo',
                                  'gds',
                                  'homologene',
                                  'journals',
                                  'mesh',
                                  'ncbisearch',
                                  'nlmcatalog',
                                  'omia',
                                  'omim',
                                  'pmc',
                                  'popset',
                                  'probe',
                                  'proteinclusters',
                                  'pcassay',
                                  'pccompound',
                                  'pcsubstance',
                                  'snp',
                                  'taxonomy',
                                  'toolkit',
                                  'unigene',
                                  'unists'
                                 ]
    def t_pubmed(self):
        '''Test parsing database info returned by EInfo
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.einfo(db="pubmed")
        input = open('Entrez/einfo2.xml')
        record = Entrez.read(input)
        assert record["DbInfo"]['DbName']=='pubmed'
        assert record["DbInfo"]['MenuName']=='PubMed'
        assert record["DbInfo"]['Description']=='PubMed bibliographic record'
        assert record["DbInfo"]['Count']=="17905967"
        assert record["DbInfo"]['LastUpdate']=='2008/04/15 06:42'

        assert len(record["DbInfo"]['FieldList'])==40

        assert record["DbInfo"]['FieldList'][0]['Name']=='ALL'
        assert record["DbInfo"]['FieldList'][0]['FullName']=='All Fields'
        assert record["DbInfo"]['FieldList'][0]['Description']=='All terms from all searchable fields'
        assert record["DbInfo"]['FieldList'][0]['TermCount']=="70792830"
        assert record["DbInfo"]['FieldList'][0]['IsDate']=='N'
        assert record["DbInfo"]['FieldList'][0]['IsNumerical']=='N'
        assert record["DbInfo"]['FieldList'][0]['SingleToken']=='N'
        assert record["DbInfo"]['FieldList'][0]['Hierarchy']=='N'
        assert record["DbInfo"]['FieldList'][0]['IsHidden']=='N'

        assert len(record["DbInfo"]['LinkList'])==46

        assert record["DbInfo"]['LinkList'][0]['Name']=='pubmed_books_refs'
        assert record["DbInfo"]['LinkList'][0]['Menu']=='Cited in Books'
        assert record["DbInfo"]['LinkList'][0]['Description']=='PubMed links associated with Books'
        assert record["DbInfo"]['LinkList'][0]['DbTo']=='books'

class ESearchTest(unittest.TestCase):
    '''Tests for parsing XML output returned by ESearch
    '''
    def t_pubmed1(self):
        '''Test parsing XML returned by ESearch from PubMed (first test)
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="biopython")
        input = open('Entrez/esearch1.xml')
        record = Entrez.read(input)
        assert record['Count']=='5'
        assert record['RetMax']=='5'
        assert record['RetStart']=='0'
        assert len(record['IdList'])==5
        assert record['IdList'][0]=='16403221'
        assert record['IdList'][1]=='16377612'
        assert record['IdList'][2]=='14871861'
        assert record['IdList'][3]=='14630660'
        assert record['IdList'][4]=='12230038'
        assert len(record['TranslationSet'])==0
        assert len(record['TranslationStack'])==2
        assert record['TranslationStack'][0]['Term']=='biopython[All Fields]'
        assert record['TranslationStack'][0]['Field']=='All Fields'
        assert record['TranslationStack'][0]['Count']=='5'
        assert record['TranslationStack'][0]['Explode']=='Y'
        assert record['TranslationStack'][1]=='GROUP'
        assert record['QueryTranslation']=='biopython[All Fields]'

    def t_pubmed2(self):
        '''Test parsing XML returned by ESearch from PubMed (second test)
        '''
        # Search in PubMed for the term cancer for the entrez date from
        # the last 60 days and retrieve the first 100 IDs and translations
        # using the history parameter. 
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="cancer", reldate=60,
        #                        datetype="edat", retmax=100, usehistory="y")
        input = open('Entrez/esearch2.xml')
        record = Entrez.read(input)
        assert record['Count']=="10238"
        assert record['RetMax']=="100"
        assert record['RetStart']=="0"
        assert record['QueryKey']=='12'
        assert record['WebEnv']=='0rYFb69LfbTFXfG7-0HPo2BU-ZFWF1s_51WtYR5e0fAzThQCR0WIW12inPQRRIj1xUzSfGgG9ovT9-@263F6CC86FF8F760_0173SID'
        assert len(record['IdList'])==100
        assert record['IdList'][0]=='18411453'
        assert record['IdList'][1]=='18411431'
        assert record['IdList'][2]=='18411430'
        assert record['IdList'][3]=='18411429'
        assert record['IdList'][4]=='18411428'
        assert record['IdList'][5]=='18411402'
        assert record['IdList'][6]=='18411381'
        assert record['IdList'][7]=='18411373'
        assert record['IdList'][8]=='18411372'
        assert record['IdList'][9]=='18411371'
        assert record['IdList'][10]=='18411370'
        assert record['IdList'][11]=='18411367'
        assert record['IdList'][12]=='18411306'
        assert record['IdList'][13]=='18411292'
        assert record['IdList'][14]=='18411277'
        assert record['IdList'][15]=='18411260'
        assert record['IdList'][16]=='18411234'
        assert record['IdList'][17]=='18411200'
        assert record['IdList'][18]=='18411199'
        assert record['IdList'][19]=='18411198'
        assert record['IdList'][20]=='18411197'
        assert record['IdList'][21]=='18411195'
        assert record['IdList'][22]=='18411194'
        assert record['IdList'][23]=='18411193'
        assert record['IdList'][24]=='18411192'
        assert record['IdList'][25]=='18411191'
        assert record['IdList'][26]=='18411052'
        assert record['IdList'][27]=='18411048'
        assert record['IdList'][28]=='18411046'
        assert record['IdList'][29]=='18411019'
        assert record['IdList'][30]=='18411018'
        assert record['IdList'][31]=='18411017'
        assert record['IdList'][32]=='18411015'
        assert record['IdList'][33]=='18411014'
        assert record['IdList'][34]=='18411011'
        assert record['IdList'][35]=='18411010'
        assert record['IdList'][36]=='18411005'
        assert record['IdList'][37]=='18411003'
        assert record['IdList'][38]=='18411001'
        assert record['IdList'][39]=='18411000'
        assert record['IdList'][40]=='18410999'
        assert record['IdList'][41]=='18410998'
        assert record['IdList'][42]=='18410997'
        assert record['IdList'][43]=='18410995'
        assert record['IdList'][44]=='18410977'
        assert record['IdList'][45]=='18410975'
        assert record['IdList'][46]=='18410966'
        assert record['IdList'][47]=='18410954'
        assert record['IdList'][48]=='18410953'
        assert record['IdList'][49]=='18410934'
        assert record['IdList'][50]=='18410925'
        assert record['IdList'][51]=='18410903'
        assert record['IdList'][52]=='18410826'
        assert record['IdList'][53]=='18410739'
        assert record['IdList'][54]=='18410720'
        assert record['IdList'][55]=='18410716'
        assert record['IdList'][56]=='18410709'
        assert record['IdList'][57]=='18410705'
        assert record['IdList'][58]=='18410692'
        assert record['IdList'][59]=='18410690'
        assert record['IdList'][60]=='18410634'
        assert record['IdList'][61]=='18410618'
        assert record['IdList'][62]=='18410610'
        assert record['IdList'][63]=='18410593'
        assert record['IdList'][64]=='18410587'
        assert record['IdList'][65]=='18410567'
        assert record['IdList'][66]=='18410539'
        assert record['IdList'][67]=='18410530'
        assert record['IdList'][68]=='18410528'
        assert record['IdList'][69]=='18410461'
        assert record['IdList'][70]=='18410455'
        assert record['IdList'][71]=='18410444'
        assert record['IdList'][72]=='18410443'
        assert record['IdList'][73]=='18410442'
        assert record['IdList'][74]=='18410441'
        assert record['IdList'][75]=='18410440'
        assert record['IdList'][76]=='18410439'
        assert record['IdList'][77]=='18410437'
        assert record['IdList'][78]=='18410436'
        assert record['IdList'][79]=='18410435'
        assert record['IdList'][80]=='18410431'
        assert record['IdList'][81]=='18410430'
        assert record['IdList'][82]=='18410428'
        assert record['IdList'][83]=='18410427'
        assert record['IdList'][84]=='18410405'
        assert record['IdList'][85]=='18410404'
        assert record['IdList'][86]=='18410355'
        assert record['IdList'][87]=='18410327'
        assert record['IdList'][88]=='18410312'
        assert record['IdList'][89]=='18410311'
        assert record['IdList'][90]=='18410307'
        assert record['IdList'][91]=='18410259'
        assert record['IdList'][92]=='18410249'
        assert record['IdList'][93]=='18410245'
        assert record['IdList'][94]=='18410243'
        assert record['IdList'][95]=='18410242'
        assert record['IdList'][96]=='18410060'
        assert record['IdList'][97]=='18410013'
        assert record['IdList'][98]=='18409992'
        assert record['IdList'][99]=='18409991'
        assert len(record['TranslationSet'])==1
        assert record['TranslationSet'][0]['From']=='cancer'
        assert record['TranslationSet'][0]['To']=='(("neoplasms"[TIAB] NOT Medline[SB]) OR "neoplasms"[MeSH Terms] OR cancer[Text Word])'
        assert len(record['TranslationStack'])==13
        assert record['TranslationStack'][0]['Term']=='"neoplasms"[TIAB]'
        assert record['TranslationStack'][0]['Field']=='TIAB'
        assert record['TranslationStack'][0]['Count']=="52104"
        assert record['TranslationStack'][0]['Explode']=='Y'
        assert record['TranslationStack'][1]['Term']=='Medline[SB]'
        assert record['TranslationStack'][1]['Field']=='SB'
        assert record['TranslationStack'][1]['Count']=="16509514"
        assert record['TranslationStack'][1]['Explode']=='Y'
        assert record['TranslationStack'][2]=='NOT'
        assert record['TranslationStack'][3]=='GROUP'
        assert record['TranslationStack'][4]['Term']=='"neoplasms"[MeSH Terms]'
        assert record['TranslationStack'][4]['Field']=='MeSH Terms'
        assert record['TranslationStack'][4]['Count']=="1918010"
        assert record['TranslationStack'][4]['Explode']=='Y'
        assert record['TranslationStack'][5]=='OR'
        assert record['TranslationStack'][6]['Term']=='cancer[Text Word]'
        assert record['TranslationStack'][6]['Field']=='Text Word'
        assert record['TranslationStack'][6]['Count']=="638849"
        assert record['TranslationStack'][6]['Explode']=='Y'
        assert record['TranslationStack'][7]=='OR'
        assert record['TranslationStack'][8]=='GROUP'
        assert record['TranslationStack'][9]['Term']=='2008/02/16[EDAT]'
        assert record['TranslationStack'][9]['Field']=='EDAT'
        assert record['TranslationStack'][9]['Count']=="-1"
        assert record['TranslationStack'][9]['Explode']=='Y'
        assert record['TranslationStack'][10]['Term']=='2008/04/16[EDAT]'
        assert record['TranslationStack'][10]['Field']=='EDAT'
        assert record['TranslationStack'][10]['Count']=="-1"
        assert record['TranslationStack'][10]['Explode']=='Y'
        assert record['TranslationStack'][11]=='RANGE'
        assert record['TranslationStack'][12]=='AND'
        assert record['QueryTranslation']=='(("neoplasms"[TIAB] NOT Medline[SB]) OR "neoplasms"[MeSH Terms] OR cancer[Text Word]) AND 2008/02/16[EDAT] : 2008/04/16[EDAT]'

    def t_pubmed3(self):
        '''Test parsing XML returned by ESearch from PubMed (third test)
        '''
        # Search in PubMed for the journal PNAS Volume 97, and retrieve
        # 6 IDs starting at ID 7.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pubmed", term="PNAS[ta] AND 97[vi]",
        #                        retstart=6, retmax=6)
        input = open('Entrez/esearch3.xml')
        record = Entrez.read(input)
        assert record['Count']=='2652'
        assert record['RetMax']=='6'
        assert record['RetStart']=='6'
        assert len(record['IdList'])==6
        assert record['IdList'][0]=='11121077'
        assert record['IdList'][1]=='11121076'
        assert record['IdList'][2]=='11121075'
        assert record['IdList'][3]=='11121074'
        assert record['IdList'][4]=='11121073'
        assert record['IdList'][5]=='11121072'
        assert len(record['TranslationSet'])==1
        assert record['TranslationSet'][0]['From']=='PNAS[ta]'
        assert record['TranslationSet'][0]['To']=='"Proc Natl Acad Sci U S A"[Journal:__jrid6653]'
        assert len(record['TranslationStack'])==3
        assert record['TranslationStack'][0]['Term']=='"Proc Natl Acad Sci U S A"[Journal]'
        assert record['TranslationStack'][0]['Field']=='Journal'
        assert record['TranslationStack'][0]['Count']=='91806'
        assert record['TranslationStack'][0]['Explode']=='Y'
        assert record['TranslationStack'][1]['Term']=='97[vi]'
        assert record['TranslationStack'][1]['Field']=='vi'
        assert record['TranslationStack'][1]['Count']=='58681'
        assert record['TranslationStack'][1]['Explode']=='Y'
        assert record['TranslationStack'][2]=='AND'
        assert record['QueryTranslation']=='"Proc Natl Acad Sci U S A"[Journal] AND 97[vi]'

    def t_journals(self):
        '''Test parsing XML returned by ESearch from the Journals database
        '''
        # Search in Journals for the term obstetrics.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="journals", term="obstetrics")
        input = open('Entrez/esearch4.xml')
        record = Entrez.read(input)
        assert record['Count']=='177'
        assert record['RetMax']=='20'
        assert record['RetStart']=='0'
        assert len(record['IdList'])==20
        assert record['IdList'][0]=='75'
        assert record['IdList'][1]=='138'
        assert record['IdList'][2]=='136'
        assert record['IdList'][3]=='137'
        assert record['IdList'][4]=='139'
        assert record['IdList'][5]=='140'
        assert record['IdList'][6]=='355'
        assert record['IdList'][7]=='354'
        assert record['IdList'][8]=='27731'
        assert record['IdList'][9]=='439'
        assert record['IdList'][10]=='564'
        assert record['IdList'][11]=='617'
        assert record['IdList'][12]=='749'
        assert record['IdList'][13]=='735'
        assert record['IdList'][14]=='815'
        assert record['IdList'][15]=='905'
        assert record['IdList'][16]=='903'
        assert record['IdList'][17]=='932'
        assert record['IdList'][18]=='933'
        assert record['IdList'][19]=='875'
        assert len(record['TranslationSet'])==0
        assert len(record['TranslationStack'])==2
        assert record['TranslationStack'][0]['Term']=='obstetrics[All Fields]'
        assert record['TranslationStack'][0]['Field']=='All Fields'
        assert record['TranslationStack'][0]['Count']=='177'
        assert record['TranslationStack'][0]['Explode']=='Y'
        assert record['TranslationStack'][0].tag=="TermSet"
        assert record['TranslationStack'][1]=='GROUP'
        assert record['TranslationStack'][1].tag=="OP"
        assert record['QueryTranslation']=='obstetrics[All Fields]'

    def t_pmc(self):
        '''Test parsing XML returned by ESearch from PubMed Central
        '''
        # Search in PubMed Central for stem cells in free fulltext articles.
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="pmc",
        #                        term="stem cells AND free fulltext[filter]")
        input = open('Entrez/esearch5.xml')
        record = Entrez.read(input)
        assert record['Count']=='23492'
        assert record['RetMax']=='20'
        assert record['RetStart']=='0'
        assert len(record['IdList'])==20
        assert record['IdList'][0]=='1894783'
        assert record['IdList'][1]=='2064507'
        assert record['IdList'][2]=='520747'
        assert record['IdList'][3]=='2043120'
        assert record['IdList'][4]=='2118723'
        assert record['IdList'][5]=='1815228'
        assert record['IdList'][6]=='1253596'
        assert record['IdList'][7]=='2077853'
        assert record['IdList'][8]=='1308908'
        assert record['IdList'][9]=='2233634'
        assert record['IdList'][10]=='556262'
        assert record['IdList'][11]=='1925137'
        assert record['IdList'][12]=='1860068'
        assert record['IdList'][13]=='1626529'
        assert record['IdList'][14]=='2217616'
        assert record['IdList'][15]=='1584276'
        assert record['IdList'][16]=='2000702'
        assert record['IdList'][17]=='186324'
        assert record['IdList'][18]=='1959362'
        assert record['IdList'][19]=='1413911'
        assert len(record['TranslationSet'])==1
        assert record['TranslationSet'][0]['From']=='stem cells'
        assert record['TranslationSet'][0]['To']=='("stem cells"[MeSH Terms] OR stem cells[Acknowledgments] OR stem cells[Figure/Table Caption] OR stem cells[Section Title] OR stem cells[Body - All Words] OR stem cells[Title] OR stem cells[Abstract])'
        assert len(record['TranslationStack'])==16
        assert record['TranslationStack'][0]['Term']=='"stem cells"[MeSH Terms]'
        assert record['TranslationStack'][0]['Field']=='MeSH Terms'
        assert record['TranslationStack'][0]['Count']=='12224'
        assert record['TranslationStack'][0]['Explode']=='Y'
        assert record['TranslationStack'][1]['Term']=='stem cells[Acknowledgments]'
        assert record['TranslationStack'][1]['Field']=='Acknowledgments'
        assert record['TranslationStack'][1]['Count']=='79'
        assert record['TranslationStack'][1]['Explode']=='Y'
        assert record['TranslationStack'][2]=='OR'
        assert record['TranslationStack'][3]['Term']=='stem cells[Figure/Table Caption]'
        assert record['TranslationStack'][3]['Field']=='Figure/Table Caption'
        assert record['TranslationStack'][3]['Count']=='806'
        assert record['TranslationStack'][3]['Explode']=='Y'
        assert record['TranslationStack'][4]=='OR'
        assert record['TranslationStack'][5]['Term']=='stem cells[Section Title]'
        assert record['TranslationStack'][5]['Field']=='Section Title'
        assert record['TranslationStack'][5]['Count']=='522'
        assert record['TranslationStack'][5]['Explode']=='Y'
        assert record['TranslationStack'][6]=='OR'
        assert record['TranslationStack'][7]['Term']=='stem cells[Body - All Words]'
        assert record['TranslationStack'][7]['Field']=='Body - All Words'
        assert record['TranslationStack'][7]['Count']=='13936'
        assert record['TranslationStack'][7]['Explode']=='Y'
        assert record['TranslationStack'][8]=='OR'
        assert record['TranslationStack'][9]['Term']=='stem cells[Title]'
        assert record['TranslationStack'][9]['Field']=='Title'
        assert record['TranslationStack'][9]['Count']=='1005'
        assert record['TranslationStack'][9]['Explode']=='Y'
        assert record['TranslationStack'][10]=='OR'
        assert record['TranslationStack'][11]['Term']=='stem cells[Abstract]'
        assert record['TranslationStack'][11]['Field']=='Abstract'
        assert record['TranslationStack'][11]['Count']=='2503'
        assert record['TranslationStack'][11]['Explode']=='Y'
        assert record['TranslationStack'][12]=='OR'
        assert record['TranslationStack'][13]=='GROUP'
        assert record['TranslationStack'][14]['Term']=='free fulltext[filter]'
        assert record['TranslationStack'][14]['Field']=='filter'
        assert record['TranslationStack'][14]['Count']=='1412839'
        assert record['TranslationStack'][14]['Explode']=='Y'
        assert record['TranslationStack'][15]=='AND'
        assert record['QueryTranslation']=='("stem cells"[MeSH Terms] OR stem cells[Acknowledgments] OR stem cells[Figure/Table Caption] OR stem cells[Section Title] OR stem cells[Body - All Words] OR stem cells[Title] OR stem cells[Abstract]) AND free fulltext[filter]'

    def t_nucleotide(self):
        '''Test parsing XML returned by ESearch from the Nucleotide database
        '''
        # Search in Nucleotide for a property of the sequence,
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="nucleotide", term="biomol trna[prop]")
        input = open('Entrez/esearch6.xml')
        record = Entrez.read(input)
        assert record['Count']=="699"
        assert record['RetMax']=="20"
        assert record['RetStart']=="0"
        assert len(record['IdList'])==20
        assert record['IdList'][0]=='220161'
        assert record['IdList'][1]=='220160'
        assert record['IdList'][2]=='220159'
        assert record['IdList'][3]=='220263'
        assert record['IdList'][4]=='220162'
        assert record['IdList'][5]=='159885659'
        assert record['IdList'][6]=='156572228'
        assert record['IdList'][7]=='2648075'
        assert record['IdList'][8]=='287595'
        assert record['IdList'][9]=='402544'
        assert record['IdList'][10]=='402506'
        assert record['IdList'][11]=='402505'
        assert record['IdList'][12]=='287594'
        assert record['IdList'][13]=='287593'
        assert record['IdList'][14]=='287592'
        assert record['IdList'][15]=='287591'
        assert record['IdList'][16]=='287590'
        assert record['IdList'][17]=='287589'
        assert record['IdList'][18]=='287588'
        assert record['IdList'][19]=='287587'
        assert len(record['TranslationSet'])==0
        assert record['QueryTranslation']==''

    def t_protein(self):
        '''Test parsing XML returned by ESearch from the Protein database
        '''
        # Search in Protein for a molecular weight
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="protein", term="200020[molecular weight]")
        input = open('Entrez/esearch7.xml')
        record = Entrez.read(input)
        assert record['Count']=='3'
        assert record['RetMax']=='3'
        assert record['RetStart']=='0'
        assert len(record['IdList'])==3
        assert record['IdList'][0]=='16766766'
        assert record['IdList'][1]=='16422035'
        assert record['IdList'][2]=='4104812'
        assert len(record['TranslationSet'])==0
        assert len(record['TranslationStack'])==2
        assert record['TranslationStack'][0]['Term']=='000200020[molecular weight]'
        assert record['TranslationStack'][0]['Field']=='molecular weight'
        assert record['TranslationStack'][0]['Count']=='3'
        assert record['TranslationStack'][0]['Explode']=='Y'
        assert record['TranslationStack'][1]=='GROUP'
        assert record['QueryTranslation']=='000200020[molecular weight]'

    def t_notfound(self):
        '''Test parsing XML returned by ESearch when no items were found
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="protein", term="abcXYZ")
        input = open('Entrez/esearch8.xml')
        record = Entrez.read(input)
        assert record['Count']=="0"
        assert record['RetMax']=="0"
        assert record['RetStart']=="0"
        assert len(record['IdList'])==0
        assert len(record['TranslationSet'])==0
        assert record['QueryTranslation']==''
        assert len(record['ErrorList'])==2
        assert "PhraseNotFound" in record['ErrorList']
        assert "FieldNotFound" in record['ErrorList']
        assert len(record['ErrorList']["PhraseNotFound"])==1
        assert len(record['ErrorList']["FieldNotFound"])==0
        assert record['ErrorList']["PhraseNotFound"][0]=="abcXYZ"
        assert len(record['WarningList'])==3
        assert "PhraseIgnored" in record['WarningList']
        assert "QuotedPhraseNotFound" in record['WarningList']
        assert "OutputMessage" in record['WarningList']
        assert len(record['WarningList']["PhraseIgnored"])==0
        assert len(record['WarningList']["QuotedPhraseNotFound"])==0
        assert len(record['WarningList']["OutputMessage"])==1
        assert record['WarningList']["OutputMessage"][0]=="No items found."

class EPostTest(unittest.TestCase):
    '''Tests for parsing XML output returned by EPost
    '''
    # Don't know how to get an InvalidIdList in the XML returned by EPost;
    # unable to test if we are parsing it correctly.
    def t_epost(self):
        '''Test parsing XML returned by EPost
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="pubmed", id="11237011")
        input = open('Entrez/epost1.xml')
        record = Entrez.read(input)
        assert record["QueryKey"]== '1'
        assert record["WebEnv"]=='0zYsuLk3zG_lRMkblPBEqnT8nIENUGw4HAy8xXChTnoVm7GEnWY71jv3nz@1FC077F3806DE010_0042SID'

    def t_wrong(self):
        '''Test parsing XML returned by EPost with incorrect arguments
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="nothing")
        input = open('Entrez/epost2.xml')
        exception_triggered = False
        try:
            record = Entrez.read(input)
        except RuntimeError, exception:
            assert str(exception)=="Wrong DB name"
            exception_triggered = True
        assert exception_triggered


    def t_invalid(self):
        '''Test parsing XML returned by EPost with an invalid id (overflow tag)
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.epost(db="pubmed", id=99999999999999999999999999999999)
        input = open('Entrez/epost3.xml')
        record = Entrez.read(input)

        assert record["InvalidIdList"]==["-1"]
        assert record["QueryKey"]=="1"
        assert record["WebEnv"]=="08AIUeBsfIk6BfdzKnd3GM2RtCudczC9jm5aeb4US0o7azCTQCeCsr-xg0@1EDE54E680D03C40_0011SID"


class ESummaryTest(unittest.TestCase):
    '''Tests for parsing XML output returned by ESummary
    '''
    # Items have a type, which can be
    # (Integer|Date|String|Structure|List|Flags|Qualifier|Enumerator|Unknown)
    # I don't have an XML file where the type "Flags", "Qualifier",
    # "Enumerator", or "Unknown" is used, so they are not tested here.
    def t_pubmed(self):
        '''Test parsing XML returned by ESummary from PubMed
        '''
        # In PubMed display records for PMIDs 11850928 and 11482001 in
        # xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="pubmed", id=["11850928","11482001"],
        #                         retmode="xml")
        input = open('Entrez/esummary1.xml')
        record = Entrez.read(input)
        assert record[0]["Id"]=="11850928"
        assert record[0]["PubDate"]=="1965 Aug"
        assert record[0]["EPubDate"]==""
        assert record[0]["Source"]=="Arch Dermatol"
        assert len(record[0]["AuthorList"])==2
        assert record[0]["AuthorList"][0]=="LoPresti PJ"
        assert record[0]["AuthorList"][1]=="Hambrick GW Jr"
        assert record[0]["LastAuthor"]=="Hambrick GW Jr"
        assert record[0]["Title"]=="Zirconium granuloma following treatment of rhus dermatitis."
        assert record[0]["Volume"]=="92"
        assert record[0]["Issue"]=="2"
        assert record[0]["Pages"]=="188-91"
        assert record[0]["LangList"]==["English"]
        assert record[0]["NlmUniqueID"]=="0372433"
        assert record[0]["ISSN"]=="0003-987X"
        assert record[0]["ESSN"]=="1538-3652"
        assert len(record[0]["PubTypeList"])==1
        assert record[0]["PubTypeList"][0]=="Journal Article"
        assert record[0]["RecordStatus"]=="PubMed - indexed for MEDLINE"
        assert record[0]["PubStatus"]=="ppublish"
        assert len(record[0]["ArticleIds"])==2
        assert record[0]["ArticleIds"]["pubmed"]==["11850928"]
        assert record[0]["ArticleIds"]["medline"]==[]
        assert len(record[0]["History"])==2
        assert record[0]["History"]["pubmed"]==["1965/08/01 00:00"]
        assert record[0]["History"]["medline"]==["2002/03/09 10:01"]
        assert len(record[0]["References"])==0
        assert record[0]["HasAbstract"]==1
        assert record[0]["PmcRefCount"]==0
        assert record[0]["FullJournalName"]=="Archives of dermatology"
        assert record[0]["ELocationID"]==""
        assert record[0]["SO"]=="1965 Aug;92(2):188-91"

        assert record[1]["Id"]=="11482001"
        assert record[1]["PubDate"]=="2001 Jun"
        assert record[1]["EPubDate"]==""
        assert record[1]["Source"]=="Adverse Drug React Toxicol Rev"
        assert len(record[1]["AuthorList"])==3
        assert record[1]["AuthorList"][0]=="Mantle D"
        assert record[1]["AuthorList"][1]=="Gok MA"
        assert record[1]["AuthorList"][2]=="Lennard TW"
        assert record[1]["LastAuthor"]=="Lennard TW"
        assert record[1]["Title"]=="Adverse and beneficial effects of plant extracts on skin and skin disorders."
        assert record[1]["Volume"]=="20"
        assert record[1]["Issue"]=="2"
        assert record[1]["Pages"]=="89-103"
        assert len(record[1]["LangList"])==1
        assert record[1]["LangList"][0]=="English"
        assert record[1]["NlmUniqueID"]=="9109474"
        assert record[1]["ISSN"]=="0964-198X"
        assert record[1]["ESSN"]==""
        assert len(record[1]["PubTypeList"])==2
        assert record[1]["PubTypeList"][0]=="Journal Article"
        assert record[1]["PubTypeList"][1]=="Review"
        assert record[1]["RecordStatus"]=="PubMed - indexed for MEDLINE"
        assert record[1]["PubStatus"]=="ppublish"
        assert len(record[1]["ArticleIds"])==2
        assert record[1]["ArticleIds"]["pubmed"]==["11482001"]
        assert record[1]["ArticleIds"]["medline"]==[]
        assert len(record[1]["History"])==2
        assert record[1]["History"]["pubmed"]==["2001/08/03 10:00"]
        assert record[1]["History"]["medline"]==["2002/01/23 10:01"]
        assert len(record[1]["References"])==0
        assert record[1]["HasAbstract"]==1
        assert record[1]["PmcRefCount"]==0
        assert record[1]["FullJournalName"]=="Adverse drug reactions and toxicological reviews"
        assert record[1]["ELocationID"]==""
        assert record[1]["SO"]=="2001 Jun;20(2):89-103"

    def t_journals(self):
        '''Test parsing XML returned by ESummary from the Journals database
        '''
        # In Journals display records for journal IDs 27731,439,735,905 
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="journals", id="27731,439,735,905")
        input = open('Entrez/esummary2.xml')
        record = Entrez.read(input)
        assert record[0]["Id"]=="27731"
        assert record[0]["Title"]=="The American journal of obstetrics and diseases of women and children"
        assert record[0]["MedAbbr"]=="Am J Obstet Dis Women Child"
        assert record[0]["IsoAbbr"]==""
        assert record[0]["NlmId"]=="14820330R"
        assert record[0]["pISSN"]=="0894-5543"
        assert record[0]["eISSN"]==""
        assert record[0]["PublicationStartYear"]=="1868"
        assert record[0]["PublicationEndYear"]=="1919"
        assert record[0]["Publisher"]=="W.A. Townsend & Adams, $c [1868-1919]"
        assert record[0]["Language"]=="eng"
        assert record[0]["Country"]=="United States"
        assert len(record[0]["BroadHeading"])==0
        assert record[0]["ContinuationNotes"]==""

        assert record[1]["Id"]=="439"
        assert record[1]["Title"]=="American journal of obstetrics and gynecology"
        assert record[1]["MedAbbr"]=="Am J Obstet Gynecol"
        assert record[1]["IsoAbbr"]=="Am. J. Obstet. Gynecol."
        assert record[1]["NlmId"]=="0370476"
        assert record[1]["pISSN"]=="0002-9378"
        assert record[1]["eISSN"]=="1097-6868"
        assert record[1]["PublicationStartYear"]=="1920"
        assert record[1]["PublicationEndYear"]==""
        assert record[1]["Publisher"]=="Elsevier,"
        assert record[1]["Language"]=="eng"
        assert record[1]["Country"]=="United States"
        assert len(record[1]["BroadHeading"])==2
        assert record[1]["BroadHeading"][0]=="Gynecology"
        assert record[1]["BroadHeading"][1]=="Obstetrics"
        assert record[1]["ContinuationNotes"]=="Continues: American journal of obstetrics and diseases of women and children. "

        assert record[2]["Id"]=="735"
        assert record[2]["Title"]=="Archives of gynecology and obstetrics"
        assert record[2]["MedAbbr"]=="Arch Gynecol Obstet"
        assert record[2]["IsoAbbr"]=="Arch. Gynecol. Obstet."
        assert record[2]["NlmId"]=="8710213"
        assert record[2]["pISSN"]=="0932-0067"
        assert record[2]["eISSN"]=="1432-0711"
        assert record[2]["PublicationStartYear"]=="1987"
        assert record[2]["PublicationEndYear"]==""
        assert record[2]["Publisher"]=="Springer Verlag"
        assert record[2]["Language"]=="eng"
        assert record[2]["Country"]=="Germany"
        assert len(record[2]["BroadHeading"])==2
        assert record[2]["BroadHeading"][0]=="Gynecology"
        assert record[2]["BroadHeading"][1]=="Obstetrics"
        assert record[2]["ContinuationNotes"]=="Continues: Archives of gynecology. "

        assert record[3]["Id"]=="905"
        assert record[3]["Title"]=="Asia-Oceania journal of obstetrics and gynaecology / AOFOG"
        assert record[3]["MedAbbr"]=="Asia Oceania J Obstet Gynaecol"
        assert record[3]["IsoAbbr"]==""
        assert record[3]["NlmId"]=="8102781"
        assert record[3]["pISSN"]=="0389-2328"
        assert record[3]["eISSN"]==""
        assert record[3]["PublicationStartYear"]=="1980"
        assert record[3]["PublicationEndYear"]=="1994"
        assert record[3]["Publisher"]=="University Of Tokyo Press"
        assert record[3]["Language"]=="eng"
        assert record[3]["Country"]=="Japan"
        assert len(record[3]["BroadHeading"])==2
        assert record[3]["BroadHeading"][0]=="Gynecology"
        assert record[3]["BroadHeading"][1]=="Obstetrics"
        assert record[3]["ContinuationNotes"]=="Continues: Journal of the Asian Federation of Obstetrics and Gynaecology. Continued by: Journal of obstetrics and gynaecology (Tokyo, Japan). "

    def t_protein(self):
        '''Test parsing XML returned by ESummary from the Protein database
        '''
        # In Protein display records for GIs 28800982 and 28628843 in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="protein", id="28800982,28628843", retmode="xml")
        input = open('Entrez/esummary3.xml')
        record = Entrez.read(input)

        assert record[0]["Id"]=="28800982"
        assert record[0]["Caption"]=="AAO47091"
        assert record[0]["Title"]=="hemochromatosis [Homo sapiens]"
        assert record[0]["Extra"]=="gi|28800982|gb|AAO47091.1|[28800982]"
        assert record[0]["Gi"]==28800982
        assert record[0]["CreateDate"]=="2003/03/03"
        assert record[0]["UpdateDate"]=="2003/03/03"
        assert record[0]["Flags"]==0
        assert record[0]["TaxId"]==9606
        assert record[0]["Length"]==268
        assert record[0]["Status"]=="live"
        assert record[0]["ReplacedBy"]==""
        assert record[0]["Comment"]=="  "

        assert record[1]["Id"]=="28628843"
        assert record[1]["Caption"]=="AAO49381"
        assert record[1]["Title"]=="erythroid associated factor [Homo sapiens]"
        assert record[1]["Extra"]=="gi|28628843|gb|AAO49381.1|AF485325_1[28628843]"
        assert record[1]["Gi"]== 28628843
        assert record[1]["CreateDate"]=="2003/03/02"
        assert record[1]["UpdateDate"]=="2003/03/02"
        assert record[1]["Flags"]==0
        assert record[1]["TaxId"]==9606
        assert record[1]["Length"]==102
        assert record[1]["Status"]=="live"
        assert record[1]["ReplacedBy"]==""
        assert record[1]["Comment"]=="  "

    def t_nucleotide(self):
        '''Test parsing XML returned by ESummary from the Nucleotide database
        '''
        # In Nucleotide display records for GIs 28864546 and 28800981
        # in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="nucleotide", id="28864546,28800981",
        #                         retmode="xml")
        input = open('Entrez/esummary4.xml')
        record = Entrez.read(input)

        assert record[0]["Id"]=="28864546"
        assert record[0]["Caption"]=="AY207443"
        assert record[0]["Title"]=="Homo sapiens alpha hemoglobin (HBZP) pseudogene 3' UTR/AluJo repeat breakpoint junction"
        assert record[0]["Extra"]=="gi|28864546|gb|AY207443.1|[28864546]"
        assert record[0]["Gi"]==28864546
        assert record[0]["CreateDate"]=="2003/03/05"
        assert record[0]["UpdateDate"]=="2003/03/05"
        assert record[0]["Flags"]==0
        assert record[0]["TaxId"]==9606
        assert record[0]["Length"]==491
        assert record[0]["Status"]=="live"
        assert record[0]["ReplacedBy"]==""
        assert record[0]["Comment"]=="  "

        assert record[1]["Id"]=="28800981"
        assert record[1]["Caption"]=="AY205604"
        assert record[1]["Title"]=="Homo sapiens hemochromatosis (HFE) mRNA, partial cds"
        assert record[1]["Extra"]=="gi|28800981|gb|AY205604.1|[28800981]"
        assert record[1]["Gi"]==28800981
        assert record[1]["CreateDate"]=="2003/03/03"
        assert record[1]["UpdateDate"]=="2003/03/03"
        assert record[1]["Flags"]==0
        assert record[1]["TaxId"]==9606
        assert record[1]["Length"]==860
        assert record[1]["Status"]=="live"
        assert record[1]["ReplacedBy"]==""
        assert record[1]["Comment"]=="  "

    def t_structure(self):
        '''Test parsing XML returned by ESummary from the Structure database
        '''
        # In Nucleotide display records for GIs 28864546 and 28800981
        # in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="structure", id=["19923","12120"],
        #                         retmode="xml")
        input = open('Entrez/esummary5.xml')
        record = Entrez.read(input)
        assert record[0]["Id"]=="19923"
        assert record[0]["PdbAcc"]=="1L5J"
        assert record[0]["PdbDescr"]=="Crystal Structure Of E. Coli Aconitase B"
        assert record[0]["EC"]=="4.2.1.3"
        assert record[0]["Resolution"]=="2.4"
        assert record[0]["ExpMethod"]=="X-Ray Diffraction"
        assert record[0]["PdbClass"]=="Lyase"
        assert record[0]["PdbReleaseDate"]=="2007/8/27"
        assert record[0]["PdbDepositDate"]=="2002/3/7"
        assert record[0]["DepositDate"]=="2007/10/25"
        assert record[0]["ModifyDate"]=="2007/10/25"
        assert record[0]["LigCode"]=="F3S|TRA"
        assert record[0]["LigCount"]=="2"
        assert record[0]["ModProteinResCount"]=="0"
        assert record[0]["ModDNAResCount"]=="0"
        assert record[0]["ModRNAResCount"]=="0"
        assert record[0]["ProteinChainCount"]=="2"
        assert record[0]["DNAChainCount"]=="0"
        assert record[0]["RNAChainCount"]=="0"

        assert record[1]["Id"]=="12120"
        assert record[1]["PdbAcc"]=="1B0K"
        assert record[1]["PdbDescr"]=="S642a:fluorocitrate Complex Of Aconitase"
        assert record[1]["EC"]=="4.2.1.3"
        assert record[1]["Resolution"]=="2.5"
        assert record[1]["ExpMethod"]=="X-Ray Diffraction"
        assert record[1]["PdbClass"]=="Lyase"
        assert record[1]["PdbReleaseDate"]=="2007/8/27"
        assert record[1]["PdbDepositDate"]=="1998/11/11"
        assert record[1]["DepositDate"]=="2007/10/07"
        assert record[1]["ModifyDate"]=="2007/10/07"
        assert record[1]["LigCode"]=="FLC|O|SF4"
        assert record[1]["LigCount"]=="3"
        assert record[1]["ModProteinResCount"]=="0"
        assert record[1]["ModDNAResCount"]=="0"
        assert record[1]["ModRNAResCount"]=="0"
        assert record[1]["ProteinChainCount"]=="1"
        assert record[1]["DNAChainCount"]=="0"
        assert record[1]["RNAChainCount"]=="0"

    def t_taxonomy(self):
        '''Test parsing XML returned by ESummary from the Taxonomy database
        '''
        # In Taxonomy display records for TAXIDs 9913 and 30521 in
        # xml retrieval mode 
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="taxonomy", id=["9913","30521"],
        #                         retmode="xml")
        input = open('Entrez/esummary6.xml')
        record = Entrez.read(input)
        assert record[0]["Id"]=="9913"
        assert record[0]["Rank"]=="species"
        assert record[0]["Division"]=="even-toed ungulates"
        assert record[0]["ScientificName"]=="Bos taurus"
        assert record[0]["CommonName"]=="cattle"
        assert record[0]["TaxId"]==9913
        assert record[0]["NucNumber"]==2264214
        assert record[0]["ProtNumber"]==55850
        assert record[0]["StructNumber"]==1517
        assert record[0]["GenNumber"]==31
        assert record[0]["GeneNumber"]==29651
        assert record[0]["Genus"]==""
        assert record[0]["Species"]==""
        assert record[0]["Subsp"]==""

        assert record[1]["Id"]=="30521"
        assert record[1]["Rank"]=="species"
        assert record[1]["Division"]=="even-toed ungulates"
        assert record[1]["ScientificName"]=="Bos grunniens"
        assert record[1]["CommonName"]=="domestic yak"
        assert record[1]["TaxId"]==30521
        assert record[1]["NucNumber"]==560
        assert record[1]["ProtNumber"]==254
        assert record[1]["StructNumber"]==0
        assert record[1]["GenNumber"]==1
        assert record[1]["GeneNumber"]==13
        assert record[1]["Genus"]==""
        assert record[1]["Species"]==""
        assert record[1]["Subsp"]==""

    def t_unists(self):
        '''Test parsing XML returned by ESummary from the UniSTS database
        '''
        # In UniSTS display records for IDs 254085 and 254086 in xml
        # retrieval mode 
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="unists", id=["254085","254086"],
        #                         retmode="xml")
        input = open('Entrez/esummary7.xml')
        record = Entrez.read(input)

        assert record[0]["Id"]=="254085"
        assert record[0]["Marker_Name"]=="SE234324"
        assert len(record[0]["Map_Gene_Summary_List"])==1
        assert record[0]["Map_Gene_Summary_List"][0]["Org"]=="Sus scrofa"
        assert record[0]["Map_Gene_Summary_List"][0]["Chr"]==" chromosome 7"
        assert record[0]["Map_Gene_Summary_List"][0]["Locus"]==""
        assert record[0]["EPCR_Summary"]=="Found by e-PCR in sequences from Sus scrofa."
        assert record[0]["LocusId"]==""

        assert record[1]["Id"]=="254086"
        assert record[1]["Marker_Name"]=="SE259162"
        assert len(record[1]["Map_Gene_Summary_List"])==1
        assert record[1]["Map_Gene_Summary_List"][0]["Org"]=="Sus scrofa"
        assert record[1]["Map_Gene_Summary_List"][0]["Chr"]==" chromosome 12"
        assert record[1]["Map_Gene_Summary_List"][0]["Locus"]==""
        assert record[1]["EPCR_Summary"]=="Found by e-PCR in sequences from Sus scrofa."
        assert record[1]["LocusId"]==""

    def t_wrong(self):
        '''Test parsing XML returned by ESummary with incorrect arguments
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.esummary()
        input = open('Entrez/esummary8.xml')
        exception_triggered = False
        try:
            record = Entrez.read(input)
        except RuntimeError, exception:
            assert str(exception)=="Neither query_key nor id specified"
            exception_triggered = True
        assert exception_triggered


class ELinkTest(unittest.TestCase):
    '''Tests for parsing XML output returned by ELink
    '''
    def t_pubmed1(self):
        '''Test parsing pubmed links returned by ELink (first test)
        '''
        # Retrieve IDs from PubMed for PMID 9298984 to the PubMed database
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="9298984", cmd="neighbor")
        input = open('Entrez/elink1.xml')
        record = Entrez.read(input)

        assert len(record)==1
        assert record[0]["DbFrom"]=="pubmed"
        assert record[0]["IdList"]==["9298984"]
        assert record[0]["LinkSetDb"][0]["DbTo"]=="pubmed"
        assert record[0]["LinkSetDb"][0]["LinkName"]=="pubmed_pubmed"
        assert len(record[0]["LinkSetDb"][0]["Link"])==144
        assert record[0]["LinkSetDb"][0]["Link"][0]["Id"]=="9298984"
        assert record[0]["LinkSetDb"][0]["Link"][0]["Score"]=="2147483647"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Id"]=="8794856"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Score"]=="65259341"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Id"]=="9700164"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Score"]=="60347327"
        assert record[0]["LinkSetDb"][0]["Link"][3]["Id"]=="7914521"
        assert record[0]["LinkSetDb"][0]["Link"][3]["Score"]=="54343405"
        assert record[0]["LinkSetDb"][0]["Link"][4]["Id"]=="1339459"
        assert record[0]["LinkSetDb"][0]["Link"][4]["Score"]=="53014422"
        assert record[0]["LinkSetDb"][0]["Link"][5]["Id"]=="9914369"
        assert record[0]["LinkSetDb"][0]["Link"][5]["Score"]=="52741538"
        assert record[0]["LinkSetDb"][0]["Link"][6]["Id"]=="11590237"
        assert record[0]["LinkSetDb"][0]["Link"][6]["Score"]=="52493903"
        assert record[0]["LinkSetDb"][0]["Link"][7]["Id"]=="12686595"
        assert record[0]["LinkSetDb"][0]["Link"][7]["Score"]=="48734007"
        assert record[0]["LinkSetDb"][0]["Link"][8]["Id"]=="9074495"
        assert record[0]["LinkSetDb"][0]["Link"][8]["Score"]=="48220447"
        assert record[0]["LinkSetDb"][0]["Link"][9]["Id"]=="11146659"
        assert record[0]["LinkSetDb"][0]["Link"][9]["Score"]=="46604626"
        assert record[0]["LinkSetDb"][0]["Link"][10]["Id"]=="10893249"
        assert record[0]["LinkSetDb"][0]["Link"][10]["Score"]=="46254167"
        assert record[0]["LinkSetDb"][0]["Link"][11]["Id"]=="8978614"
        assert record[0]["LinkSetDb"][0]["Link"][11]["Score"]=="46166362"
        assert record[0]["LinkSetDb"][0]["Link"][12]["Id"]=="15371539"
        assert record[0]["LinkSetDb"][0]["Link"][12]["Score"]=="45060488"
        assert record[0]["LinkSetDb"][0]["Link"][13]["Id"]=="10806105"
        assert record[0]["LinkSetDb"][0]["Link"][13]["Score"]=="44825774"
        assert record[0]["LinkSetDb"][0]["Link"][14]["Id"]=="10402457"
        assert record[0]["LinkSetDb"][0]["Link"][14]["Score"]=="44338522"
        assert record[0]["LinkSetDb"][0]["Link"][15]["Id"]=="10545493"
        assert record[0]["LinkSetDb"][0]["Link"][15]["Score"]=="43860609"
        assert record[0]["LinkSetDb"][0]["Link"][16]["Id"]=="10523511"
        assert record[0]["LinkSetDb"][0]["Link"][16]["Score"]=="43268800"
        assert record[0]["LinkSetDb"][0]["Link"][17]["Id"]=="12515822"
        assert record[0]["LinkSetDb"][0]["Link"][17]["Score"]=="43215343"
        assert record[0]["LinkSetDb"][0]["Link"][18]["Id"]=="15915585"
        assert record[0]["LinkSetDb"][0]["Link"][18]["Score"]=="43029760"
        assert record[0]["LinkSetDb"][0]["Link"][19]["Id"]=="11483958"
        assert record[0]["LinkSetDb"][0]["Link"][19]["Score"]=="42348877"
        assert record[0]["LinkSetDb"][0]["Link"][20]["Id"]=="11685532"
        assert record[0]["LinkSetDb"][0]["Link"][20]["Score"]=="42262104"
        assert record[0]["LinkSetDb"][0]["Link"][21]["Id"]=="9869638"
        assert record[0]["LinkSetDb"][0]["Link"][21]["Score"]=="41143593"
        assert record[0]["LinkSetDb"][0]["Link"][22]["Id"]=="12080088"
        assert record[0]["LinkSetDb"][0]["Link"][22]["Score"]=="40849490"
        assert record[0]["LinkSetDb"][0]["Link"][23]["Id"]=="12034769"
        assert record[0]["LinkSetDb"][0]["Link"][23]["Score"]=="40841328"
        assert record[0]["LinkSetDb"][0]["Link"][24]["Id"]=="9852156"
        assert record[0]["LinkSetDb"][0]["Link"][24]["Score"]=="40793501"
        assert record[0]["LinkSetDb"][0]["Link"][25]["Id"]=="9735366"
        assert record[0]["LinkSetDb"][0]["Link"][25]["Score"]=="40661605"
        assert record[0]["LinkSetDb"][0]["Link"][26]["Id"]=="10749938"
        assert record[0]["LinkSetDb"][0]["Link"][26]["Score"]=="40486739"
        assert record[0]["LinkSetDb"][0]["Link"][27]["Id"]=="9490715"
        assert record[0]["LinkSetDb"][0]["Link"][27]["Score"]=="40311339"
        assert record[0]["LinkSetDb"][0]["Link"][28]["Id"]=="9425896"
        assert record[0]["LinkSetDb"][0]["Link"][28]["Score"]=="40056298"
        assert record[0]["LinkSetDb"][0]["Link"][29]["Id"]=="11266459"
        assert record[0]["LinkSetDb"][0]["Link"][29]["Score"]=="39883140"
        assert record[0]["LinkSetDb"][0]["Link"][30]["Id"]=="14522947"
        assert record[0]["LinkSetDb"][0]["Link"][30]["Score"]=="39683976"
        assert record[0]["LinkSetDb"][0]["Link"][31]["Id"]=="15616189"
        assert record[0]["LinkSetDb"][0]["Link"][31]["Score"]=="39518630"
        assert record[0]["LinkSetDb"][0]["Link"][32]["Id"]=="16732327"
        assert record[0]["LinkSetDb"][0]["Link"][32]["Score"]=="39425668"
        assert record[0]["LinkSetDb"][0]["Link"][33]["Id"]=="11179694"
        assert record[0]["LinkSetDb"][0]["Link"][33]["Score"]=="39183565"
        assert record[0]["LinkSetDb"][0]["Link"][34]["Id"]=="10898791"
        assert record[0]["LinkSetDb"][0]["Link"][34]["Score"]=="39159761"
        assert record[0]["LinkSetDb"][0]["Link"][35]["Id"]=="11146661"
        assert record[0]["LinkSetDb"][0]["Link"][35]["Score"]=="39116609"
        assert record[0]["LinkSetDb"][0]["Link"][36]["Id"]=="11914278"
        assert record[0]["LinkSetDb"][0]["Link"][36]["Score"]=="39028004"
        assert record[0]["LinkSetDb"][0]["Link"][37]["Id"]=="10985388"
        assert record[0]["LinkSetDb"][0]["Link"][37]["Score"]=="39002572"
        assert record[0]["LinkSetDb"][0]["Link"][38]["Id"]=="16839185"
        assert record[0]["LinkSetDb"][0]["Link"][38]["Score"]=="38916726"
        assert record[0]["LinkSetDb"][0]["Link"][39]["Id"]=="7585942"
        assert record[0]["LinkSetDb"][0]["Link"][39]["Score"]=="38747288"
        assert record[0]["LinkSetDb"][0]["Link"][40]["Id"]=="2022189"
        assert record[0]["LinkSetDb"][0]["Link"][40]["Score"]=="38717145"
        assert record[0]["LinkSetDb"][0]["Link"][41]["Id"]=="7690762"
        assert record[0]["LinkSetDb"][0]["Link"][41]["Score"]=="38647275"
        assert record[0]["LinkSetDb"][0]["Link"][42]["Id"]=="7904902"
        assert record[0]["LinkSetDb"][0]["Link"][42]["Score"]=="38557343"
        assert record[0]["LinkSetDb"][0]["Link"][43]["Id"]=="9378750"
        assert record[0]["LinkSetDb"][0]["Link"][43]["Score"]=="38484849"
        assert record[0]["LinkSetDb"][0]["Link"][44]["Id"]=="12388768"
        assert record[0]["LinkSetDb"][0]["Link"][44]["Score"]=="38454422"
        assert record[0]["LinkSetDb"][0]["Link"][45]["Id"]=="11352945"
        assert record[0]["LinkSetDb"][0]["Link"][45]["Score"]=="38449836"
        assert record[0]["LinkSetDb"][0]["Link"][46]["Id"]=="11267866"
        assert record[0]["LinkSetDb"][0]["Link"][46]["Score"]=="38419058"
        assert record[0]["LinkSetDb"][0]["Link"][47]["Id"]=="17222555"
        assert record[0]["LinkSetDb"][0]["Link"][47]["Score"]=="38368546"
        assert record[0]["LinkSetDb"][0]["Link"][48]["Id"]=="11252055"
        assert record[0]["LinkSetDb"][0]["Link"][48]["Score"]=="38257516"
        assert record[0]["LinkSetDb"][0]["Link"][49]["Id"]=="16585270"
        assert record[0]["LinkSetDb"][0]["Link"][49]["Score"]=="37800856"
        assert record[0]["LinkSetDb"][0]["Link"][50]["Id"]=="9606208"
        assert record[0]["LinkSetDb"][0]["Link"][50]["Score"]=="37669054"
        assert record[0]["LinkSetDb"][0]["Link"][51]["Id"]=="17182852"
        assert record[0]["LinkSetDb"][0]["Link"][51]["Score"]=="37621285"
        assert record[0]["LinkSetDb"][0]["Link"][52]["Id"]=="9933569"
        assert record[0]["LinkSetDb"][0]["Link"][52]["Score"]=="37398470"
        assert record[0]["LinkSetDb"][0]["Link"][53]["Id"]=="15268859"
        assert record[0]["LinkSetDb"][0]["Link"][53]["Score"]=="37340582"
        assert record[0]["LinkSetDb"][0]["Link"][54]["Id"]=="12235289"
        assert record[0]["LinkSetDb"][0]["Link"][54]["Score"]=="37247450"
        assert record[0]["LinkSetDb"][0]["Link"][55]["Id"]=="16741559"
        assert record[0]["LinkSetDb"][0]["Link"][55]["Score"]=="37198716"
        assert record[0]["LinkSetDb"][0]["Link"][56]["Id"]=="11266451"
        assert record[0]["LinkSetDb"][0]["Link"][56]["Score"]=="37142542"
        assert record[0]["LinkSetDb"][0]["Link"][57]["Id"]=="15075237"
        assert record[0]["LinkSetDb"][0]["Link"][57]["Score"]=="36897578"
        assert record[0]["LinkSetDb"][0]["Link"][58]["Id"]=="15485811"
        assert record[0]["LinkSetDb"][0]["Link"][58]["Score"]=="36804297"
        assert record[0]["LinkSetDb"][0]["Link"][59]["Id"]=="14699129"
        assert record[0]["LinkSetDb"][0]["Link"][59]["Score"]=="36782062"
        assert record[0]["LinkSetDb"][0]["Link"][60]["Id"]=="16510521"
        assert record[0]["LinkSetDb"][0]["Link"][60]["Score"]=="36724370"
        assert record[0]["LinkSetDb"][0]["Link"][61]["Id"]=="15824131"
        assert record[0]["LinkSetDb"][0]["Link"][61]["Score"]=="36695341"
        assert record[0]["LinkSetDb"][0]["Link"][62]["Id"]=="15371340"
        assert record[0]["LinkSetDb"][0]["Link"][62]["Score"]=="36685694"
        assert record[0]["LinkSetDb"][0]["Link"][63]["Id"]=="9878245"
        assert record[0]["LinkSetDb"][0]["Link"][63]["Score"]=="36684230"
        assert record[0]["LinkSetDb"][0]["Link"][64]["Id"]=="10398680"
        assert record[0]["LinkSetDb"][0]["Link"][64]["Score"]=="36573411"
        assert record[0]["LinkSetDb"][0]["Link"][65]["Id"]=="16516834"
        assert record[0]["LinkSetDb"][0]["Link"][65]["Score"]=="36525654"
        assert record[0]["LinkSetDb"][0]["Link"][66]["Id"]=="11715021"
        assert record[0]["LinkSetDb"][0]["Link"][66]["Score"]=="36518129"
        assert record[0]["LinkSetDb"][0]["Link"][67]["Id"]=="14622138"
        assert record[0]["LinkSetDb"][0]["Link"][67]["Score"]=="36496009"
        assert record[0]["LinkSetDb"][0]["Link"][68]["Id"]=="11092768"
        assert record[0]["LinkSetDb"][0]["Link"][68]["Score"]=="36457186"
        assert record[0]["LinkSetDb"][0]["Link"][69]["Id"]=="12514103"
        assert record[0]["LinkSetDb"][0]["Link"][69]["Score"]=="36385909"
        assert record[0]["LinkSetDb"][0]["Link"][70]["Id"]=="17525528"
        assert record[0]["LinkSetDb"][0]["Link"][70]["Score"]=="36316439"
        assert record[0]["LinkSetDb"][0]["Link"][71]["Id"]=="11402064"
        assert record[0]["LinkSetDb"][0]["Link"][71]["Score"]=="36172957"
        assert record[0]["LinkSetDb"][0]["Link"][72]["Id"]=="9258677"
        assert record[0]["LinkSetDb"][0]["Link"][72]["Score"]=="35989143"
        assert record[0]["LinkSetDb"][0]["Link"][73]["Id"]=="14499625"
        assert record[0]["LinkSetDb"][0]["Link"][73]["Score"]=="35978627"
        assert record[0]["LinkSetDb"][0]["Link"][74]["Id"]=="10428958"
        assert record[0]["LinkSetDb"][0]["Link"][74]["Score"]=="35924800"
        assert record[0]["LinkSetDb"][0]["Link"][75]["Id"]=="14972679"
        assert record[0]["LinkSetDb"][0]["Link"][75]["Score"]=="35915578"
        assert record[0]["LinkSetDb"][0]["Link"][76]["Id"]=="9396743"
        assert record[0]["LinkSetDb"][0]["Link"][76]["Score"]=="35883749"
        assert record[0]["LinkSetDb"][0]["Link"][77]["Id"]=="16219694"
        assert record[0]["LinkSetDb"][0]["Link"][77]["Score"]=="35870689"
        assert record[0]["LinkSetDb"][0]["Link"][78]["Id"]=="11369198"
        assert record[0]["LinkSetDb"][0]["Link"][78]["Score"]=="35838048"
        assert record[0]["LinkSetDb"][0]["Link"][79]["Id"]=="17333235"
        assert record[0]["LinkSetDb"][0]["Link"][79]["Score"]=="35815282"
        assert record[0]["LinkSetDb"][0]["Link"][80]["Id"]=="11102811"
        assert record[0]["LinkSetDb"][0]["Link"][80]["Score"]=="35783566"
        assert record[0]["LinkSetDb"][0]["Link"][81]["Id"]=="10207147"
        assert record[0]["LinkSetDb"][0]["Link"][81]["Score"]=="35594009"
        assert record[0]["LinkSetDb"][0]["Link"][82]["Id"]=="10477755"
        assert record[0]["LinkSetDb"][0]["Link"][82]["Score"]=="35589601"
        assert record[0]["LinkSetDb"][0]["Link"][83]["Id"]=="10747094"
        assert record[0]["LinkSetDb"][0]["Link"][83]["Score"]=="35548072"
        assert record[0]["LinkSetDb"][0]["Link"][84]["Id"]=="15215209"
        assert record[0]["LinkSetDb"][0]["Link"][84]["Score"]=="35526869"
        assert record[0]["LinkSetDb"][0]["Link"][85]["Id"]=="11157774"
        assert record[0]["LinkSetDb"][0]["Link"][85]["Score"]=="35510607"
        assert record[0]["LinkSetDb"][0]["Link"][86]["Id"]=="10669599"
        assert record[0]["LinkSetDb"][0]["Link"][86]["Score"]=="35462246"
        assert record[0]["LinkSetDb"][0]["Link"][87]["Id"]=="17448445"
        assert record[0]["LinkSetDb"][0]["Link"][87]["Score"]=="35398470"
        assert record[0]["LinkSetDb"][0]["Link"][88]["Id"]=="17878237"
        assert record[0]["LinkSetDb"][0]["Link"][88]["Score"]=="35231311"
        assert record[0]["LinkSetDb"][0]["Link"][89]["Id"]=="10411903"
        assert record[0]["LinkSetDb"][0]["Link"][89]["Score"]=="35202708"
        assert record[0]["LinkSetDb"][0]["Link"][90]["Id"]=="12773390"
        assert record[0]["LinkSetDb"][0]["Link"][90]["Score"]=="35171743"
        assert record[0]["LinkSetDb"][0]["Link"][91]["Id"]=="12498686"
        assert record[0]["LinkSetDb"][0]["Link"][91]["Score"]=="35131906"
        assert record[0]["LinkSetDb"][0]["Link"][92]["Id"]=="9009204"
        assert record[0]["LinkSetDb"][0]["Link"][92]["Score"]=="34993776"
        assert record[0]["LinkSetDb"][0]["Link"][93]["Id"]=="17576797"
        assert record[0]["LinkSetDb"][0]["Link"][93]["Score"]=="34988639"
        assert record[0]["LinkSetDb"][0]["Link"][94]["Id"]=="10225945"
        assert record[0]["LinkSetDb"][0]["Link"][94]["Score"]=="34950419"
        assert record[0]["LinkSetDb"][0]["Link"][95]["Id"]=="11161560"
        assert record[0]["LinkSetDb"][0]["Link"][95]["Score"]=="34912466"
        assert record[0]["LinkSetDb"][0]["Link"][96]["Id"]=="11967147"
        assert record[0]["LinkSetDb"][0]["Link"][96]["Score"]=="34900540"
        assert record[0]["LinkSetDb"][0]["Link"][97]["Id"]=="14711415"
        assert record[0]["LinkSetDb"][0]["Link"][97]["Score"]=="34883714"
        assert record[0]["LinkSetDb"][0]["Link"][98]["Id"]=="2211824"
        assert record[0]["LinkSetDb"][0]["Link"][98]["Score"]=="34843507"
        assert record[0]["LinkSetDb"][0]["Link"][99]["Id"]=="15737064"
        assert record[0]["LinkSetDb"][0]["Link"][99]["Score"]=="34828187"
        assert record[0]["LinkSetDb"][0]["Link"][100]["Id"]=="7720068"
        assert record[0]["LinkSetDb"][0]["Link"][100]["Score"]=="34811182"
        assert record[0]["LinkSetDb"][0]["Link"][101]["Id"]=="9472001"
        assert record[0]["LinkSetDb"][0]["Link"][101]["Score"]=="34799321"
        assert record[0]["LinkSetDb"][0]["Link"][102]["Id"]=="11792803"
        assert record[0]["LinkSetDb"][0]["Link"][102]["Score"]=="34697393"
        assert record[0]["LinkSetDb"][0]["Link"][103]["Id"]=="11386760"
        assert record[0]["LinkSetDb"][0]["Link"][103]["Score"]=="34684610"
        assert record[0]["LinkSetDb"][0]["Link"][104]["Id"]=="15094189"
        assert record[0]["LinkSetDb"][0]["Link"][104]["Score"]=="34684021"
        assert record[0]["LinkSetDb"][0]["Link"][105]["Id"]=="9763420"
        assert record[0]["LinkSetDb"][0]["Link"][105]["Score"]=="34666950"
        assert record[0]["LinkSetDb"][0]["Link"][106]["Id"]=="10629219"
        assert record[0]["LinkSetDb"][0]["Link"][106]["Score"]=="34422925"
        assert record[0]["LinkSetDb"][0]["Link"][107]["Id"]=="11238410"
        assert record[0]["LinkSetDb"][0]["Link"][107]["Score"]=="34318521"
        assert record[0]["LinkSetDb"][0]["Link"][108]["Id"]=="17199038"
        assert record[0]["LinkSetDb"][0]["Link"][108]["Score"]=="34255594"
        assert record[0]["LinkSetDb"][0]["Link"][109]["Id"]=="12944469"
        assert record[0]["LinkSetDb"][0]["Link"][109]["Score"]=="34249507"
        assert record[0]["LinkSetDb"][0]["Link"][110]["Id"]=="15616192"
        assert record[0]["LinkSetDb"][0]["Link"][110]["Score"]=="34110517"
        assert record[0]["LinkSetDb"][0]["Link"][111]["Id"]=="11146660"
        assert record[0]["LinkSetDb"][0]["Link"][111]["Score"]=="34063257"
        assert record[0]["LinkSetDb"][0]["Link"][112]["Id"]=="11402066"
        assert record[0]["LinkSetDb"][0]["Link"][112]["Score"]=="34012520"
        assert record[0]["LinkSetDb"][0]["Link"][113]["Id"]=="6791901"
        assert record[0]["LinkSetDb"][0]["Link"][113]["Score"]=="33311119"
        assert record[0]["LinkSetDb"][0]["Link"][114]["Id"]=="7172865"
        assert record[0]["LinkSetDb"][0]["Link"][114]["Score"]=="32934223"
        assert record[0]["LinkSetDb"][0]["Link"][115]["Id"]=="8270646"
        assert record[0]["LinkSetDb"][0]["Link"][115]["Score"]=="32898701"
        assert record[0]["LinkSetDb"][0]["Link"][116]["Id"]=="1916263"
        assert record[0]["LinkSetDb"][0]["Link"][116]["Score"]=="32707765"
        assert record[0]["LinkSetDb"][0]["Link"][117]["Id"]=="7588080"
        assert record[0]["LinkSetDb"][0]["Link"][117]["Score"]=="32503526"
        assert record[0]["LinkSetDb"][0]["Link"][118]["Id"]=="7391142"
        assert record[0]["LinkSetDb"][0]["Link"][118]["Score"]=="31632645"
        assert record[0]["LinkSetDb"][0]["Link"][119]["Id"]=="6793236"
        assert record[0]["LinkSetDb"][0]["Link"][119]["Score"]=="31522175"
        assert record[0]["LinkSetDb"][0]["Link"][120]["Id"]=="2512302"
        assert record[0]["LinkSetDb"][0]["Link"][120]["Score"]=="30339372"
        assert record[0]["LinkSetDb"][0]["Link"][121]["Id"]=="7720069"
        assert record[0]["LinkSetDb"][0]["Link"][121]["Score"]=="30024525"
        assert record[0]["LinkSetDb"][0]["Link"][122]["Id"]=="8257792"
        assert record[0]["LinkSetDb"][0]["Link"][122]["Score"]=="29834355"
        assert record[0]["LinkSetDb"][0]["Link"][123]["Id"]=="3417141"
        assert record[0]["LinkSetDb"][0]["Link"][123]["Score"]=="27920818"
        assert record[0]["LinkSetDb"][0]["Link"][124]["Id"]=="3315496"
        assert record[0]["LinkSetDb"][0]["Link"][124]["Score"]=="27422009"
        assert record[0]["LinkSetDb"][0]["Link"][125]["Id"]=="1993311"
        assert record[0]["LinkSetDb"][0]["Link"][125]["Score"]=="26763828"
        assert record[0]["LinkSetDb"][0]["Link"][126]["Id"]=="6185450"
        assert record[0]["LinkSetDb"][0]["Link"][126]["Score"]=="26100420"
        assert record[0]["LinkSetDb"][0]["Link"][127]["Id"]=="1819515"
        assert record[0]["LinkSetDb"][0]["Link"][127]["Score"]=="26036804"
        assert record[0]["LinkSetDb"][0]["Link"][128]["Id"]=="7250964"
        assert record[0]["LinkSetDb"][0]["Link"][128]["Score"]=="25738652"
        assert record[0]["LinkSetDb"][0]["Link"][129]["Id"]=="8489280"
        assert record[0]["LinkSetDb"][0]["Link"][129]["Score"]=="25587858"
        assert record[0]["LinkSetDb"][0]["Link"][130]["Id"]=="7096444"
        assert record[0]["LinkSetDb"][0]["Link"][130]["Score"]=="24642544"
        assert record[0]["LinkSetDb"][0]["Link"][131]["Id"]=="348629"
        assert record[0]["LinkSetDb"][0]["Link"][131]["Score"]=="24432498"
        assert record[0]["LinkSetDb"][0]["Link"][132]["Id"]=="2275018"
        assert record[0]["LinkSetDb"][0]["Link"][132]["Score"]=="23077593"
        assert record[0]["LinkSetDb"][0]["Link"][133]["Id"]=="1747872"
        assert record[0]["LinkSetDb"][0]["Link"][133]["Score"]=="22933494"
        assert record[0]["LinkSetDb"][0]["Link"][134]["Id"]=="3547036"
        assert record[0]["LinkSetDb"][0]["Link"][134]["Score"]=="22925639"
        assert record[0]["LinkSetDb"][0]["Link"][135]["Id"]=="18291669"
        assert record[0]["LinkSetDb"][0]["Link"][135]["Score"]=="22762310"
        assert record[0]["LinkSetDb"][0]["Link"][136]["Id"]=="1576878"
        assert record[0]["LinkSetDb"][0]["Link"][136]["Score"]=="20846041"
        assert record[0]["LinkSetDb"][0]["Link"][137]["Id"]=="6230555"
        assert record[0]["LinkSetDb"][0]["Link"][137]["Score"]=="19354488"
        assert record[0]["LinkSetDb"][0]["Link"][138]["Id"]=="7627547"
        assert record[0]["LinkSetDb"][0]["Link"][138]["Score"]=="18940607"
        assert record[0]["LinkSetDb"][0]["Link"][139]["Id"]=="17678444"
        assert record[0]["LinkSetDb"][0]["Link"][139]["Score"]=="18834135"
        assert record[0]["LinkSetDb"][0]["Link"][140]["Id"]=="3366468"
        assert record[0]["LinkSetDb"][0]["Link"][140]["Score"]=="14831756"
        assert record[0]["LinkSetDb"][0]["Link"][141]["Id"]=="1959920"
        assert record[0]["LinkSetDb"][0]["Link"][141]["Score"]=="14156414"
        assert record[0]["LinkSetDb"][0]["Link"][142]["Id"]=="13242628"
        assert record[0]["LinkSetDb"][0]["Link"][142]["Score"]=="12584732"
        assert record[0]["LinkSetDb"][0]["Link"][143]["Id"]=="17248312"
        assert record[0]["LinkSetDb"][0]["Link"][143]["Score"]=="7610436"

    def t_nucleotide(self):
        '''Test parsing Nucleotide to Protein links returned by ELink
        '''
        # Retrieve IDs from Nucleotide for GI  48819, 7140345 to Protein
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="nucleotide", db="protein",
        #                      id="48819,7140345")
        input = open('Entrez/elink2.xml')
        record = Entrez.read(input)

        assert record[0]["DbFrom"]=="nucleotide"
        assert record[0]["IdList"]==["48819", "7140345"]

    def t_pubmed2(self):
        '''Test parsing pubmed links returned by ELink (second test)
        '''
        # Retrieve PubMed related articles for PMIDs 11812492 11774222
        # with a publication date from 1995 to the present
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="11812492,11774222",
        #                      db="pubmed", mindate="1995", datetype="pdat")
        input = open('Entrez/elink3.xml')
        record = Entrez.read(input)

        assert len(record)==1
        assert record[0]["DbFrom"]=="pubmed"
        assert len(record[0]['IdList'])==2
        assert record[0]['IdList'][0]=="11812492"
        assert record[0]['IdList'][1]=="11774222"
        assert record[0]["LinkSetDb"][0]["DbTo"]=="pubmed"
        assert record[0]["LinkSetDb"][0]["LinkName"]=="pubmed_pubmed"
        assert record[0]["LinkSetDb"][0]["Link"][0]["Id"]=="11812492"
        assert record[0]["LinkSetDb"][0]["Link"][0]["Score"]=="2147483647"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Id"]=="11774222"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Score"]=="2147483647"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Id"]=="11668631"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Score"]=="86345306"
        assert record[0]["LinkSetDb"][0]["Link"][3]["Id"]=="15111095"
        assert record[0]["LinkSetDb"][0]["Link"][3]["Score"]=="81604359"
        assert record[0]["LinkSetDb"][0]["Link"][4]["Id"]=="10731564"
        assert record[0]["LinkSetDb"][0]["Link"][4]["Score"]=="65665112"
        assert record[0]["LinkSetDb"][0]["Link"][5]["Id"]=="15780005"
        assert record[0]["LinkSetDb"][0]["Link"][5]["Score"]=="62251079"
        assert record[0]["LinkSetDb"][0]["Link"][6]["Id"]=="17885136"
        assert record[0]["LinkSetDb"][0]["Link"][6]["Score"]=="50322134"
        assert record[0]["LinkSetDb"][0]["Link"][7]["Id"]=="17470297"
        assert record[0]["LinkSetDb"][0]["Link"][7]["Score"]=="49148434"
        assert record[0]["LinkSetDb"][0]["Link"][8]["Id"]=="16005284"
        assert record[0]["LinkSetDb"][0]["Link"][8]["Score"]=="49035508"
        assert record[0]["LinkSetDb"][0]["Link"][9]["Id"]=="10856373"
        assert record[0]["LinkSetDb"][0]["Link"][9]["Score"]=="48363137"
        assert record[0]["LinkSetDb"][0]["Link"][10]["Id"]=="15383292"
        assert record[0]["LinkSetDb"][0]["Link"][10]["Score"]=="48347159"
        assert record[0]["LinkSetDb"][0]["Link"][11]["Id"]=="17040125"
        assert record[0]["LinkSetDb"][0]["Link"][11]["Score"]=="48301243"
        assert record[0]["LinkSetDb"][0]["Link"][12]["Id"]=="10770808"
        assert record[0]["LinkSetDb"][0]["Link"][12]["Score"]=="47696325"
        assert record[0]["LinkSetDb"][0]["Link"][13]["Id"]=="11125122"
        assert record[0]["LinkSetDb"][0]["Link"][13]["Score"]=="45889695"
        assert record[0]["LinkSetDb"][0]["Link"][14]["Id"]=="15287587"
        assert record[0]["LinkSetDb"][0]["Link"][14]["Score"]=="45599733"
        assert record[0]["LinkSetDb"][0]["Link"][15]["Id"]=="15839745"
        assert record[0]["LinkSetDb"][0]["Link"][15]["Score"]=="44650620"
        assert record[0]["LinkSetDb"][0]["Link"][16]["Id"]=="10612825"
        assert record[0]["LinkSetDb"][0]["Link"][16]["Score"]=="44445812"
        assert record[0]["LinkSetDb"][0]["Link"][17]["Id"]=="15024419"
        assert record[0]["LinkSetDb"][0]["Link"][17]["Score"]=="44075047"
        assert record[0]["LinkSetDb"][0]["Link"][18]["Id"]=="12743802"
        assert record[0]["LinkSetDb"][0]["Link"][18]["Score"]=="43873158"
        assert record[0]["LinkSetDb"][0]["Link"][19]["Id"]=="15238684"
        assert record[0]["LinkSetDb"][0]["Link"][19]["Score"]=="43856864"
        assert record[0]["LinkSetDb"][0]["Link"][20]["Id"]=="12386340"
        assert record[0]["LinkSetDb"][0]["Link"][20]["Score"]=="43770229"
        assert record[0]["LinkSetDb"][0]["Link"][21]["Id"]=="16269725"
        assert record[0]["LinkSetDb"][0]["Link"][21]["Score"]=="43712594"
        assert record[0]["LinkSetDb"][0]["Link"][22]["Id"]=="10592273"
        assert record[0]["LinkSetDb"][0]["Link"][22]["Score"]=="43640108"
        assert record[0]["LinkSetDb"][0]["Link"][23]["Id"]=="15383308"
        assert record[0]["LinkSetDb"][0]["Link"][23]["Score"]=="42835474"
        assert record[0]["LinkSetDb"][0]["Link"][24]["Id"]=="15676075"
        assert record[0]["LinkSetDb"][0]["Link"][24]["Score"]=="42272663"
        assert record[0]["LinkSetDb"][0]["Link"][25]["Id"]=="11774221"
        assert record[0]["LinkSetDb"][0]["Link"][25]["Score"]=="42058380"
        assert record[0]["LinkSetDb"][0]["Link"][26]["Id"]=="10592272"
        assert record[0]["LinkSetDb"][0]["Link"][26]["Score"]=="41719917"
        assert record[0]["LinkSetDb"][0]["Link"][27]["Id"]=="15997407"
        assert record[0]["LinkSetDb"][0]["Link"][27]["Score"]=="41535461"
        assert record[0]["LinkSetDb"][0]["Link"][28]["Id"]=="15774024"
        assert record[0]["LinkSetDb"][0]["Link"][28]["Score"]=="41351079"
        assert record[0]["LinkSetDb"][0]["Link"][29]["Id"]=="11233160"
        assert record[0]["LinkSetDb"][0]["Link"][29]["Score"]=="41268965"
        assert record[0]["LinkSetDb"][0]["Link"][30]["Id"]=="14702162"
        assert record[0]["LinkSetDb"][0]["Link"][30]["Score"]=="41147661"
        assert record[0]["LinkSetDb"][0]["Link"][31]["Id"]=="16616613"
        assert record[0]["LinkSetDb"][0]["Link"][31]["Score"]=="41073100"
        assert record[0]["LinkSetDb"][0]["Link"][32]["Id"]=="17202370"
        assert record[0]["LinkSetDb"][0]["Link"][32]["Score"]=="40819600"
        assert record[0]["LinkSetDb"][0]["Link"][33]["Id"]=="15478601"
        assert record[0]["LinkSetDb"][0]["Link"][33]["Score"]=="40578911"
        assert record[0]["LinkSetDb"][0]["Link"][34]["Id"]=="15322925"
        assert record[0]["LinkSetDb"][0]["Link"][34]["Score"]=="40548101"
        assert record[0]["LinkSetDb"][0]["Link"][35]["Id"]=="11472559"
        assert record[0]["LinkSetDb"][0]["Link"][35]["Score"]=="40508356"
        assert record[0]["LinkSetDb"][0]["Link"][36]["Id"]=="11925998"
        assert record[0]["LinkSetDb"][0]["Link"][36]["Score"]=="39844751"
        assert record[0]["LinkSetDb"][0]["Link"][37]["Id"]=="12372145"
        assert record[0]["LinkSetDb"][0]["Link"][37]["Score"]=="39809277"
        assert record[0]["LinkSetDb"][0]["Link"][38]["Id"]=="17562224"
        assert record[0]["LinkSetDb"][0]["Link"][38]["Score"]=="38850094"
        assert record[0]["LinkSetDb"][0]["Link"][39]["Id"]=="15037105"
        assert record[0]["LinkSetDb"][0]["Link"][39]["Score"]=="38758229"
        assert record[0]["LinkSetDb"][0]["Link"][40]["Id"]=="14998511"
        assert record[0]["LinkSetDb"][0]["Link"][40]["Score"]=="38608049"
        assert record[0]["LinkSetDb"][0]["Link"][41]["Id"]=="10092480"
        assert record[0]["LinkSetDb"][0]["Link"][41]["Score"]=="38410463"
        assert record[0]["LinkSetDb"][0]["Link"][42]["Id"]=="7729881"
        assert record[0]["LinkSetDb"][0]["Link"][42]["Score"]=="38329800"
        assert record[0]["LinkSetDb"][0]["Link"][43]["Id"]=="12933853"
        assert record[0]["LinkSetDb"][0]["Link"][43]["Score"]=="37881850"
        assert record[0]["LinkSetDb"][0]["Link"][44]["Id"]=="16818783"
        assert record[0]["LinkSetDb"][0]["Link"][44]["Score"]=="37835096"
        assert record[0]["LinkSetDb"][0]["Link"][45]["Id"]=="16406333"
        assert record[0]["LinkSetDb"][0]["Link"][45]["Score"]=="37775136"
        assert record[0]["LinkSetDb"][0]["Link"][46]["Id"]=="11472553"
        assert record[0]["LinkSetDb"][0]["Link"][46]["Score"]=="37750745"
        assert record[0]["LinkSetDb"][0]["Link"][47]["Id"]=="11403387"
        assert record[0]["LinkSetDb"][0]["Link"][47]["Score"]=="37707525"
        assert record[0]["LinkSetDb"][0]["Link"][48]["Id"]=="17306254"
        assert record[0]["LinkSetDb"][0]["Link"][48]["Score"]=="37685833"
        assert record[0]["LinkSetDb"][0]["Link"][49]["Id"]=="11516587"
        assert record[0]["LinkSetDb"][0]["Link"][49]["Score"]=="37620966"
        assert record[0]["LinkSetDb"][0]["Link"][50]["Id"]=="9274032"
        assert record[0]["LinkSetDb"][0]["Link"][50]["Score"]=="37528832"
        assert record[0]["LinkSetDb"][0]["Link"][51]["Id"]=="12856318"
        assert record[0]["LinkSetDb"][0]["Link"][51]["Score"]=="37484650"
        assert record[0]["LinkSetDb"][0]["Link"][52]["Id"]=="14695526"
        assert record[0]["LinkSetDb"][0]["Link"][52]["Score"]=="37429895"
        assert record[0]["LinkSetDb"][0]["Link"][53]["Id"]=="12481045"
        assert record[0]["LinkSetDb"][0]["Link"][53]["Score"]=="37051674"
        assert record[0]["LinkSetDb"][0]["Link"][54]["Id"]=="11752345"
        assert record[0]["LinkSetDb"][0]["Link"][54]["Score"]=="36875760"
        assert record[0]["LinkSetDb"][0]["Link"][55]["Id"]=="12467974"
        assert record[0]["LinkSetDb"][0]["Link"][55]["Score"]=="36787103"
        assert record[0]["LinkSetDb"][0]["Link"][56]["Id"]=="11214099"
        assert record[0]["LinkSetDb"][0]["Link"][56]["Score"]=="36710749"
        assert record[0]["LinkSetDb"][0]["Link"][57]["Id"]=="14638788"
        assert record[0]["LinkSetDb"][0]["Link"][57]["Score"]=="36667774"
        assert record[0]["LinkSetDb"][0]["Link"][58]["Id"]=="16278157"
        assert record[0]["LinkSetDb"][0]["Link"][58]["Score"]=="36598908"
        assert record[0]["LinkSetDb"][0]["Link"][59]["Id"]=="11752242"
        assert record[0]["LinkSetDb"][0]["Link"][59]["Score"]=="36555638"
        assert record[0]["LinkSetDb"][0]["Link"][60]["Id"]=="14681474"
        assert record[0]["LinkSetDb"][0]["Link"][60]["Score"]=="36317853"
        assert record[0]["LinkSetDb"][0]["Link"][61]["Id"]=="15944077"
        assert record[0]["LinkSetDb"][0]["Link"][61]["Score"]=="36264027"
        assert record[0]["LinkSetDb"][0]["Link"][62]["Id"]=="12625936"
        assert record[0]["LinkSetDb"][0]["Link"][62]["Score"]=="36088314"
        assert record[0]["LinkSetDb"][0]["Link"][63]["Id"]=="16672453"
        assert record[0]["LinkSetDb"][0]["Link"][63]["Score"]=="35985060"
        assert record[0]["LinkSetDb"][0]["Link"][64]["Id"]=="14695451"
        assert record[0]["LinkSetDb"][0]["Link"][64]["Score"]=="35971708"
        assert record[0]["LinkSetDb"][0]["Link"][65]["Id"]=="12402526"
        assert record[0]["LinkSetDb"][0]["Link"][65]["Score"]=="35942170"
        assert record[0]["LinkSetDb"][0]["Link"][66]["Id"]=="10592200"
        assert record[0]["LinkSetDb"][0]["Link"][66]["Score"]=="35932875"
        assert record[0]["LinkSetDb"][0]["Link"][67]["Id"]=="17584494"
        assert record[0]["LinkSetDb"][0]["Link"][67]["Score"]=="35869907"
        assert record[0]["LinkSetDb"][0]["Link"][68]["Id"]=="17761848"
        assert record[0]["LinkSetDb"][0]["Link"][68]["Score"]=="35868206"
        assert record[0]["LinkSetDb"][0]["Link"][69]["Id"]=="16697384"
        assert record[0]["LinkSetDb"][0]["Link"][69]["Score"]=="35792791"
        assert record[0]["LinkSetDb"][0]["Link"][70]["Id"]=="8784774"
        assert record[0]["LinkSetDb"][0]["Link"][70]["Score"]=="35787497"
        assert record[0]["LinkSetDb"][0]["Link"][71]["Id"]=="18000556"
        assert record[0]["LinkSetDb"][0]["Link"][71]["Score"]=="35701408"
        assert record[0]["LinkSetDb"][0]["Link"][72]["Id"]=="15828434"
        assert record[0]["LinkSetDb"][0]["Link"][72]["Score"]=="35604052"
        assert record[0]["LinkSetDb"][0]["Link"][73]["Id"]=="10511685"
        assert record[0]["LinkSetDb"][0]["Link"][73]["Score"]=="35598319"
        assert record[0]["LinkSetDb"][0]["Link"][74]["Id"]=="15608284"
        assert record[0]["LinkSetDb"][0]["Link"][74]["Score"]=="35439627"
        assert record[0]["LinkSetDb"][0]["Link"][75]["Id"]=="11125071"
        assert record[0]["LinkSetDb"][0]["Link"][75]["Score"]=="35414962"
        assert record[0]["LinkSetDb"][0]["Link"][76]["Id"]=="11791238"
        assert record[0]["LinkSetDb"][0]["Link"][76]["Score"]=="35411948"
        assert record[0]["LinkSetDb"][0]["Link"][77]["Id"]=="15710433"
        assert record[0]["LinkSetDb"][0]["Link"][77]["Score"]=="35197152"
        assert record[0]["LinkSetDb"][0]["Link"][78]["Id"]=="16164550"
        assert record[0]["LinkSetDb"][0]["Link"][78]["Score"]=="35172458"
        assert record[0]["LinkSetDb"][0]["Link"][79]["Id"]=="17697334"
        assert record[0]["LinkSetDb"][0]["Link"][79]["Score"]=="35121478"
        assert record[0]["LinkSetDb"][0]["Link"][80]["Id"]=="12537121"
        assert record[0]["LinkSetDb"][0]["Link"][80]["Score"]=="35054632"
        assert record[0]["LinkSetDb"][0]["Link"][81]["Id"]=="12860672"
        assert record[0]["LinkSetDb"][0]["Link"][81]["Score"]=="35046651"
        assert record[0]["LinkSetDb"][0]["Link"][82]["Id"]=="15630619"
        assert record[0]["LinkSetDb"][0]["Link"][82]["Score"]=="35034076"
        assert record[0]["LinkSetDb"][0]["Link"][83]["Id"]=="15125639"
        assert record[0]["LinkSetDb"][0]["Link"][83]["Score"]=="35007338"
        assert record[0]["LinkSetDb"][0]["Link"][84]["Id"]=="11443570"
        assert record[0]["LinkSetDb"][0]["Link"][84]["Score"]=="34935553"
        assert record[0]["LinkSetDb"][0]["Link"][85]["Id"]=="12208043"
        assert record[0]["LinkSetDb"][0]["Link"][85]["Score"]=="34923107"
        assert record[0]["LinkSetDb"][0]["Link"][86]["Id"]=="11731507"
        assert record[0]["LinkSetDb"][0]["Link"][86]["Score"]=="34875290"
        assert record[0]["LinkSetDb"][0]["Link"][87]["Id"]=="11988510"
        assert record[0]["LinkSetDb"][0]["Link"][87]["Score"]=="34773036"
        assert record[0]["LinkSetDb"][0]["Link"][88]["Id"]=="11125038"
        assert record[0]["LinkSetDb"][0]["Link"][88]["Score"]=="34754724"
        assert record[0]["LinkSetDb"][0]["Link"][89]["Id"]=="16381944"
        assert record[0]["LinkSetDb"][0]["Link"][89]["Score"]=="34747225"
        assert record[0]["LinkSetDb"][0]["Link"][90]["Id"]=="17135206"
        assert record[0]["LinkSetDb"][0]["Link"][90]["Score"]=="34735015"
        assert record[0]["LinkSetDb"][0]["Link"][91]["Id"]=="17099226"
        assert record[0]["LinkSetDb"][0]["Link"][91]["Score"]=="34698054"
        assert record[0]["LinkSetDb"][0]["Link"][92]["Id"]=="15608233"
        assert record[0]["LinkSetDb"][0]["Link"][92]["Score"]=="34588400"
        assert record[0]["LinkSetDb"][0]["Link"][93]["Id"]=="16672057"
        assert record[0]["LinkSetDb"][0]["Link"][93]["Score"]=="34583177"
        assert record[0]["LinkSetDb"][0]["Link"][94]["Id"]=="15687015"
        assert record[0]["LinkSetDb"][0]["Link"][94]["Score"]=="34357840"
        assert record[0]["LinkSetDb"][0]["Link"][95]["Id"]=="10782070"
        assert record[0]["LinkSetDb"][0]["Link"][95]["Score"]=="34326746"
        assert record[0]["LinkSetDb"][0]["Link"][96]["Id"]=="14970722"
        assert record[0]["LinkSetDb"][0]["Link"][96]["Score"]=="34217911"
        assert record[0]["LinkSetDb"][0]["Link"][97]["Id"]=="18027007"
        assert record[0]["LinkSetDb"][0]["Link"][97]["Score"]=="34185436"
        assert record[0]["LinkSetDb"][0]["Link"][98]["Id"]=="12387845"
        assert record[0]["LinkSetDb"][0]["Link"][98]["Score"]=="34083368"
        assert record[0]["LinkSetDb"][0]["Link"][99]["Id"]=="16237012"
        assert record[0]["LinkSetDb"][0]["Link"][99]["Score"]=="34070163"
        assert record[0]["LinkSetDb"][0]["Link"][100]["Id"]=="16351742"
        assert record[0]["LinkSetDb"][0]["Link"][100]["Score"]=="33775198"
        assert record[0]["LinkSetDb"][0]["Link"][101]["Id"]=="12203989"
        assert record[0]["LinkSetDb"][0]["Link"][101]["Score"]=="33759170"
        assert record[0]["LinkSetDb"][0]["Link"][102]["Id"]=="15474306"
        assert record[0]["LinkSetDb"][0]["Link"][102]["Score"]=="33737675"
        assert record[0]["LinkSetDb"][0]["Link"][103]["Id"]=="15270538"
        assert record[0]["LinkSetDb"][0]["Link"][103]["Score"]=="33697306"
        assert record[0]["LinkSetDb"][0]["Link"][104]["Id"]=="17518759"
        assert record[0]["LinkSetDb"][0]["Link"][104]["Score"]=="33695140"
        assert record[0]["LinkSetDb"][0]["Link"][105]["Id"]=="16085497"
        assert record[0]["LinkSetDb"][0]["Link"][105]["Score"]=="33652537"
        assert record[0]["LinkSetDb"][0]["Link"][106]["Id"]=="16423288"
        assert record[0]["LinkSetDb"][0]["Link"][106]["Score"]=="33564554"
        assert record[0]["LinkSetDb"][0]["Link"][107]["Id"]=="16251775"
        assert record[0]["LinkSetDb"][0]["Link"][107]["Score"]=="33547325"
        assert record[0]["LinkSetDb"][0]["Link"][108]["Id"]=="12632152"
        assert record[0]["LinkSetDb"][0]["Link"][108]["Score"]=="33497998"
        assert record[0]["LinkSetDb"][0]["Link"][109]["Id"]=="11269648"
        assert record[0]["LinkSetDb"][0]["Link"][109]["Score"]=="33493800"
        assert record[0]["LinkSetDb"][0]["Link"][110]["Id"]=="16103603"
        assert record[0]["LinkSetDb"][0]["Link"][110]["Score"]=="33378796"
        assert record[0]["LinkSetDb"][0]["Link"][111]["Id"]=="12816546"
        assert record[0]["LinkSetDb"][0]["Link"][111]["Score"]=="33316167"
        assert record[0]["LinkSetDb"][0]["Link"][112]["Id"]=="10221636"
        assert record[0]["LinkSetDb"][0]["Link"][112]["Score"]=="33310814"
        assert record[0]["LinkSetDb"][0]["Link"][113]["Id"]=="16381973"
        assert record[0]["LinkSetDb"][0]["Link"][113]["Score"]=="33236048"
        assert record[0]["LinkSetDb"][0]["Link"][114]["Id"]=="15977173"
        assert record[0]["LinkSetDb"][0]["Link"][114]["Score"]=="33222497"
        assert record[0]["LinkSetDb"][0]["Link"][115]["Id"]=="16351753"
        assert record[0]["LinkSetDb"][0]["Link"][115]["Score"]=="33205084"
        assert record[0]["LinkSetDb"][0]["Link"][116]["Id"]=="15317790"
        assert record[0]["LinkSetDb"][0]["Link"][116]["Score"]=="33195439"
        assert record[0]["LinkSetDb"][0]["Link"][117]["Id"]=="17135198"
        assert record[0]["LinkSetDb"][0]["Link"][117]["Score"]=="33189951"
        assert record[0]["LinkSetDb"][0]["Link"][118]["Id"]=="12701381"
        assert record[0]["LinkSetDb"][0]["Link"][118]["Score"]=="33172200"
        assert record[0]["LinkSetDb"][0]["Link"][119]["Id"]=="12203988"
        assert record[0]["LinkSetDb"][0]["Link"][119]["Score"]=="33172077"
        assert record[0]["LinkSetDb"][0]["Link"][120]["Id"]=="11456466"
        assert record[0]["LinkSetDb"][0]["Link"][120]["Score"]=="33124900"
        assert record[0]["LinkSetDb"][0]["Link"][121]["Id"]=="16936055"
        assert record[0]["LinkSetDb"][0]["Link"][121]["Score"]=="33081742"
        assert record[0]["LinkSetDb"][0]["Link"][122]["Id"]=="17183477"
        assert record[0]["LinkSetDb"][0]["Link"][122]["Score"]=="33005068"
        assert record[0]["LinkSetDb"][0]["Link"][123]["Id"]=="9455480"
        assert record[0]["LinkSetDb"][0]["Link"][123]["Score"]=="32997067"
        assert record[0]["LinkSetDb"][0]["Link"][124]["Id"]=="12490454"
        assert record[0]["LinkSetDb"][0]["Link"][124]["Score"]=="32995041"
        assert record[0]["LinkSetDb"][0]["Link"][125]["Id"]=="12435493"
        assert record[0]["LinkSetDb"][0]["Link"][125]["Score"]=="32990122"
        assert record[0]["LinkSetDb"][0]["Link"][126]["Id"]=="11038309"
        assert record[0]["LinkSetDb"][0]["Link"][126]["Score"]=="32977663"
        assert record[0]["LinkSetDb"][0]["Link"][127]["Id"]=="10366827"
        assert record[0]["LinkSetDb"][0]["Link"][127]["Score"]=="32903347"
        assert record[0]["LinkSetDb"][0]["Link"][128]["Id"]=="10466136"
        assert record[0]["LinkSetDb"][0]["Link"][128]["Score"]=="32869387"
        assert record[0]["LinkSetDb"][0]["Link"][129]["Id"]=="16381840"
        assert record[0]["LinkSetDb"][0]["Link"][129]["Score"]=="32816923"
        assert record[0]["LinkSetDb"][0]["Link"][130]["Id"]=="11825250"
        assert record[0]["LinkSetDb"][0]["Link"][130]["Score"]=="32776183"
        assert record[0]["LinkSetDb"][0]["Link"][131]["Id"]=="12234534"
        assert record[0]["LinkSetDb"][0]["Link"][131]["Score"]=="32708547"
        assert record[0]["LinkSetDb"][0]["Link"][132]["Id"]=="14624247"
        assert record[0]["LinkSetDb"][0]["Link"][132]["Score"]=="32708542"
        assert record[0]["LinkSetDb"][0]["Link"][133]["Id"]=="12886019"
        assert record[0]["LinkSetDb"][0]["Link"][133]["Score"]=="32653276"
        assert record[0]["LinkSetDb"][0]["Link"][134]["Id"]=="12041732"
        assert record[0]["LinkSetDb"][0]["Link"][134]["Score"]=="32607185"
        assert record[0]["LinkSetDb"][0]["Link"][135]["Id"]=="15336912"
        assert record[0]["LinkSetDb"][0]["Link"][135]["Score"]=="32596453"
        assert record[0]["LinkSetDb"][0]["Link"][136]["Id"]=="12652910"
        assert record[0]["LinkSetDb"][0]["Link"][136]["Score"]=="32567397"
        assert record[0]["LinkSetDb"][0]["Link"][137]["Id"]=="14681353"
        assert record[0]["LinkSetDb"][0]["Link"][137]["Score"]=="32549157"
        assert record[0]["LinkSetDb"][0]["Link"][138]["Id"]=="12586873"
        assert record[0]["LinkSetDb"][0]["Link"][138]["Score"]=="32504063"
        assert record[0]["LinkSetDb"][0]["Link"][139]["Id"]=="11481430"
        assert record[0]["LinkSetDb"][0]["Link"][139]["Score"]=="32462602"
        assert record[0]["LinkSetDb"][0]["Link"][140]["Id"]=="15254259"
        assert record[0]["LinkSetDb"][0]["Link"][140]["Score"]=="32441737"
        assert record[0]["LinkSetDb"][0]["Link"][141]["Id"]=="16873516"
        assert record[0]["LinkSetDb"][0]["Link"][141]["Score"]=="32433603"
        assert record[0]["LinkSetDb"][0]["Link"][142]["Id"]=="17170002"
        assert record[0]["LinkSetDb"][0]["Link"][142]["Score"]=="32425626"
        assert record[0]["LinkSetDb"][0]["Link"][143]["Id"]=="12519941"
        assert record[0]["LinkSetDb"][0]["Link"][143]["Score"]=="32367760"
        assert record[0]["LinkSetDb"][0]["Link"][144]["Id"]=="11197770"
        assert record[0]["LinkSetDb"][0]["Link"][144]["Score"]=="32362623"
        assert record[0]["LinkSetDb"][0]["Link"][145]["Id"]=="11240843"
        assert record[0]["LinkSetDb"][0]["Link"][145]["Score"]=="32347064"
        assert record[0]["LinkSetDb"][0]["Link"][146]["Id"]=="11328780"
        assert record[0]["LinkSetDb"][0]["Link"][146]["Score"]=="32333807"
        assert record[0]["LinkSetDb"][0]["Link"][147]["Id"]=="11875041"
        assert record[0]["LinkSetDb"][0]["Link"][147]["Score"]=="32312036"
        assert record[0]["LinkSetDb"][0]["Link"][148]["Id"]=="11752243"
        assert record[0]["LinkSetDb"][0]["Link"][148]["Score"]=="32268199"
        assert record[0]["LinkSetDb"][0]["Link"][149]["Id"]=="16907992"
        assert record[0]["LinkSetDb"][0]["Link"][149]["Score"]=="32247019"
        assert record[0]["LinkSetDb"][0]["Link"][150]["Id"]=="15046636"
        assert record[0]["LinkSetDb"][0]["Link"][150]["Score"]=="32214942"
        assert record[0]["LinkSetDb"][0]["Link"][151]["Id"]=="10592169"
        assert record[0]["LinkSetDb"][0]["Link"][151]["Score"]=="32137798"
        assert record[0]["LinkSetDb"][0]["Link"][152]["Id"]=="17919582"
        assert record[0]["LinkSetDb"][0]["Link"][152]["Score"]=="32137767"
        assert record[0]["LinkSetDb"][0]["Link"][153]["Id"]=="18025705"
        assert record[0]["LinkSetDb"][0]["Link"][153]["Score"]=="32131322"
        assert record[0]["LinkSetDb"][0]["Link"][154]["Id"]=="11029673"
        assert record[0]["LinkSetDb"][0]["Link"][154]["Score"]=="32126363"
        assert record[0]["LinkSetDb"][0]["Link"][155]["Id"]=="9047337"
        assert record[0]["LinkSetDb"][0]["Link"][155]["Score"]=="32090163"
        assert record[0]["LinkSetDb"][0]["Link"][156]["Id"]=="11080372"
        assert record[0]["LinkSetDb"][0]["Link"][156]["Score"]=="31924475"
        assert record[0]["LinkSetDb"][0]["Link"][157]["Id"]=="18045790"
        assert record[0]["LinkSetDb"][0]["Link"][157]["Score"]=="31834367"
        assert record[0]["LinkSetDb"][0]["Link"][158]["Id"]=="10215019"
        assert record[0]["LinkSetDb"][0]["Link"][158]["Score"]=="31823989"
        assert record[0]["LinkSetDb"][0]["Link"][159]["Id"]=="14706096"
        assert record[0]["LinkSetDb"][0]["Link"][159]["Score"]=="31781977"
        assert record[0]["LinkSetDb"][0]["Link"][160]["Id"]=="17537593"
        assert record[0]["LinkSetDb"][0]["Link"][160]["Score"]=="31771566"
        assert record[0]["LinkSetDb"][0]["Link"][161]["Id"]=="12819149"
        assert record[0]["LinkSetDb"][0]["Link"][161]["Score"]=="31683943"
        assert record[0]["LinkSetDb"][0]["Link"][162]["Id"]=="17880721"
        assert record[0]["LinkSetDb"][0]["Link"][162]["Score"]=="31630816"
        assert record[0]["LinkSetDb"][0]["Link"][163]["Id"]=="14681478"
        assert record[0]["LinkSetDb"][0]["Link"][163]["Score"]=="31620457"
        assert record[0]["LinkSetDb"][0]["Link"][164]["Id"]=="11985867"
        assert record[0]["LinkSetDb"][0]["Link"][164]["Score"]=="31544318"
        assert record[0]["LinkSetDb"][0]["Link"][165]["Id"]=="15608248"
        assert record[0]["LinkSetDb"][0]["Link"][165]["Score"]=="31542256"
        assert record[0]["LinkSetDb"][0]["Link"][166]["Id"]=="17401150"
        assert record[0]["LinkSetDb"][0]["Link"][166]["Score"]=="31497289"
        assert record[0]["LinkSetDb"][0]["Link"][167]["Id"]=="10359795"
        assert record[0]["LinkSetDb"][0]["Link"][167]["Score"]=="31460779"
        assert record[0]["LinkSetDb"][0]["Link"][168]["Id"]=="15608286"
        assert record[0]["LinkSetDb"][0]["Link"][168]["Score"]=="31435112"
        assert record[0]["LinkSetDb"][0]["Link"][169]["Id"]=="15774022"
        assert record[0]["LinkSetDb"][0]["Link"][169]["Score"]=="31425851"
        assert record[0]["LinkSetDb"][0]["Link"][170]["Id"]=="9921679"
        assert record[0]["LinkSetDb"][0]["Link"][170]["Score"]=="31396086"
        assert record[0]["LinkSetDb"][0]["Link"][171]["Id"]=="17038195"
        assert record[0]["LinkSetDb"][0]["Link"][171]["Score"]=="31380822"
        assert record[0]["LinkSetDb"][0]["Link"][172]["Id"]=="15491544"
        assert record[0]["LinkSetDb"][0]["Link"][172]["Score"]=="31294370"
        assert record[0]["LinkSetDb"][0]["Link"][173]["Id"]=="10469257"
        assert record[0]["LinkSetDb"][0]["Link"][173]["Score"]=="31291548"
        assert record[0]["LinkSetDb"][0]["Link"][174]["Id"]=="15487498"
        assert record[0]["LinkSetDb"][0]["Link"][174]["Score"]=="31268351"
        assert record[0]["LinkSetDb"][0]["Link"][175]["Id"]=="15383303"
        assert record[0]["LinkSetDb"][0]["Link"][175]["Score"]=="31264596"
        assert record[0]["LinkSetDb"][0]["Link"][176]["Id"]=="15643605"
        assert record[0]["LinkSetDb"][0]["Link"][176]["Score"]=="31259953"
        assert record[0]["LinkSetDb"][0]["Link"][177]["Id"]=="16418238"
        assert record[0]["LinkSetDb"][0]["Link"][177]["Score"]=="31259003"
        assert record[0]["LinkSetDb"][0]["Link"][178]["Id"]=="15500248"
        assert record[0]["LinkSetDb"][0]["Link"][178]["Score"]=="31252080"
        assert record[0]["LinkSetDb"][0]["Link"][179]["Id"]=="15479945"
        assert record[0]["LinkSetDb"][0]["Link"][179]["Score"]=="31249988"
        assert record[0]["LinkSetDb"][0]["Link"][180]["Id"]=="16962738"
        assert record[0]["LinkSetDb"][0]["Link"][180]["Score"]=="31249405"
        assert record[0]["LinkSetDb"][0]["Link"][181]["Id"]=="15094394"
        assert record[0]["LinkSetDb"][0]["Link"][181]["Score"]=="31200337"
        assert record[0]["LinkSetDb"][0]["Link"][182]["Id"]=="11758285"
        assert record[0]["LinkSetDb"][0]["Link"][182]["Score"]=="31180435"
        assert record[0]["LinkSetDb"][0]["Link"][183]["Id"]=="15723693"
        assert record[0]["LinkSetDb"][0]["Link"][183]["Score"]=="31083464"
        assert record[0]["LinkSetDb"][0]["Link"][184]["Id"]=="16710453"
        assert record[0]["LinkSetDb"][0]["Link"][184]["Score"]=="31083136"
        assert record[0]["LinkSetDb"][0]["Link"][185]["Id"]=="15311460"
        assert record[0]["LinkSetDb"][0]["Link"][185]["Score"]=="31068402"
        assert record[0]["LinkSetDb"][0]["Link"][186]["Id"]=="16549670"
        assert record[0]["LinkSetDb"][0]["Link"][186]["Score"]=="30995148"
        assert record[0]["LinkSetDb"][0]["Link"][187]["Id"]=="18180957"
        assert record[0]["LinkSetDb"][0]["Link"][187]["Score"]=="30973190"
        assert record[0]["LinkSetDb"][0]["Link"][188]["Id"]=="14681351"
        assert record[0]["LinkSetDb"][0]["Link"][188]["Score"]=="30968930"
        assert record[0]["LinkSetDb"][0]["Link"][189]["Id"]=="10902212"
        assert record[0]["LinkSetDb"][0]["Link"][189]["Score"]=="30960861"
        assert record[0]["LinkSetDb"][0]["Link"][190]["Id"]=="15357877"
        assert record[0]["LinkSetDb"][0]["Link"][190]["Score"]=="30947680"
        assert record[0]["LinkSetDb"][0]["Link"][191]["Id"]=="12356773"
        assert record[0]["LinkSetDb"][0]["Link"][191]["Score"]=="30910321"
        assert record[0]["LinkSetDb"][0]["Link"][192]["Id"]=="17537669"
        assert record[0]["LinkSetDb"][0]["Link"][192]["Score"]=="30893205"
        assert record[0]["LinkSetDb"][0]["Link"][193]["Id"]=="16551372"
        assert record[0]["LinkSetDb"][0]["Link"][193]["Score"]=="30889080"
        assert record[0]["LinkSetDb"][0]["Link"][194]["Id"]=="15231810"
        assert record[0]["LinkSetDb"][0]["Link"][194]["Score"]=="30863616"
        assert record[0]["LinkSetDb"][0]["Link"][195]["Id"]=="12819150"
        assert record[0]["LinkSetDb"][0]["Link"][195]["Score"]=="30847027"
        assert record[0]["LinkSetDb"][0]["Link"][196]["Id"]=="15608257"
        assert record[0]["LinkSetDb"][0]["Link"][196]["Score"]=="30840234"
        assert record[0]["LinkSetDb"][0]["Link"][197]["Id"]=="17384426"
        assert record[0]["LinkSetDb"][0]["Link"][197]["Score"]=="30827754"
        assert record[0]["LinkSetDb"][0]["Link"][198]["Id"]=="15811532"
        assert record[0]["LinkSetDb"][0]["Link"][198]["Score"]=="30823185"
        assert record[0]["LinkSetDb"][0]["Link"][199]["Id"]=="10612821"
        assert record[0]["LinkSetDb"][0]["Link"][199]["Score"]=="30822187"
        assert record[0]["LinkSetDb"][0]["Link"][200]["Id"]=="17062145"
        assert record[0]["LinkSetDb"][0]["Link"][200]["Score"]=="30813605"
        assert record[0]["LinkSetDb"][0]["Link"][201]["Id"]=="11355885"
        assert record[0]["LinkSetDb"][0]["Link"][201]["Score"]=="30810648"
        assert record[0]["LinkSetDb"][0]["Link"][202]["Id"]=="15746365"
        assert record[0]["LinkSetDb"][0]["Link"][202]["Score"]=="30784209"
        assert record[0]["LinkSetDb"][0]["Link"][203]["Id"]=="16282300"
        assert record[0]["LinkSetDb"][0]["Link"][203]["Score"]=="30782807"
        assert record[0]["LinkSetDb"][0]["Link"][204]["Id"]=="15546336"
        assert record[0]["LinkSetDb"][0]["Link"][204]["Score"]=="30773578"
        assert record[0]["LinkSetDb"][0]["Link"][205]["Id"]=="11741630"
        assert record[0]["LinkSetDb"][0]["Link"][205]["Score"]=="30764995"
        assert record[0]["LinkSetDb"][0]["Link"][206]["Id"]=="15980532"
        assert record[0]["LinkSetDb"][0]["Link"][206]["Score"]=="30735790"
        assert record[0]["LinkSetDb"][0]["Link"][207]["Id"]=="12519977"
        assert record[0]["LinkSetDb"][0]["Link"][207]["Score"]=="30707395"
        assert record[0]["LinkSetDb"][0]["Link"][208]["Id"]=="12436197"
        assert record[0]["LinkSetDb"][0]["Link"][208]["Score"]=="30705501"
        assert record[0]["LinkSetDb"][0]["Link"][209]["Id"]=="11125059"
        assert record[0]["LinkSetDb"][0]["Link"][209]["Score"]=="30614888"
        assert record[0]["LinkSetDb"][0]["Link"][210]["Id"]=="11163442"
        assert record[0]["LinkSetDb"][0]["Link"][210]["Score"]=="30550965"
        assert record[0]["LinkSetDb"][0]["Link"][211]["Id"]=="12519964"
        assert record[0]["LinkSetDb"][0]["Link"][211]["Score"]=="30518025"
        assert record[0]["LinkSetDb"][0]["Link"][212]["Id"]=="12083398"
        assert record[0]["LinkSetDb"][0]["Link"][212]["Score"]=="30466595"
        assert record[0]["LinkSetDb"][0]["Link"][213]["Id"]=="11908756"
        assert record[0]["LinkSetDb"][0]["Link"][213]["Score"]=="30462080"
        assert record[0]["LinkSetDb"][0]["Link"][214]["Id"]=="15608226"
        assert record[0]["LinkSetDb"][0]["Link"][214]["Score"]=="30335152"
        assert record[0]["LinkSetDb"][0]["Link"][215]["Id"]=="16845091"
        assert record[0]["LinkSetDb"][0]["Link"][215]["Score"]=="30277120"
        assert record[0]["LinkSetDb"][0]["Link"][216]["Id"]=="17338820"
        assert record[0]["LinkSetDb"][0]["Link"][216]["Score"]=="30208452"
        assert record[0]["LinkSetDb"][0]["Link"][217]["Id"]=="10407783"
        assert record[0]["LinkSetDb"][0]["Link"][217]["Score"]=="30171504"
        assert record[0]["LinkSetDb"][0]["Link"][218]["Id"]=="17130148"
        assert record[0]["LinkSetDb"][0]["Link"][218]["Score"]=="30160136"
        assert record[0]["LinkSetDb"][0]["Link"][219]["Id"]=="14681471"
        assert record[0]["LinkSetDb"][0]["Link"][219]["Score"]=="30155757"
        assert record[0]["LinkSetDb"][0]["Link"][220]["Id"]=="17445272"
        assert record[0]["LinkSetDb"][0]["Link"][220]["Score"]=="30015229"
        assert record[0]["LinkSetDb"][0]["Link"][221]["Id"]=="11279516"
        assert record[0]["LinkSetDb"][0]["Link"][221]["Score"]=="29947199"
        assert record[0]["LinkSetDb"][0]["Link"][222]["Id"]=="17221864"
        assert record[0]["LinkSetDb"][0]["Link"][222]["Score"]=="29893674"
        assert record[0]["LinkSetDb"][0]["Link"][223]["Id"]=="15827081"
        assert record[0]["LinkSetDb"][0]["Link"][223]["Score"]=="29891924"
        assert record[0]["LinkSetDb"][0]["Link"][224]["Id"]=="11222582"
        assert record[0]["LinkSetDb"][0]["Link"][224]["Score"]=="29878915"
        assert record[0]["LinkSetDb"][0]["Link"][225]["Id"]=="11384164"
        assert record[0]["LinkSetDb"][0]["Link"][225]["Score"]=="29871698"
        assert record[0]["LinkSetDb"][0]["Link"][226]["Id"]=="17877839"
        assert record[0]["LinkSetDb"][0]["Link"][226]["Score"]=="29843765"
        assert record[0]["LinkSetDb"][0]["Link"][227]["Id"]=="17151077"
        assert record[0]["LinkSetDb"][0]["Link"][227]["Score"]=="29841695"
        assert record[0]["LinkSetDb"][0]["Link"][228]["Id"]=="16381974"
        assert record[0]["LinkSetDb"][0]["Link"][228]["Score"]=="29740312"
        assert record[0]["LinkSetDb"][0]["Link"][229]["Id"]=="10592263"
        assert record[0]["LinkSetDb"][0]["Link"][229]["Score"]=="29633946"
        assert record[0]["LinkSetDb"][0]["Link"][230]["Id"]=="15608212"
        assert record[0]["LinkSetDb"][0]["Link"][230]["Score"]=="29621479"
        assert record[0]["LinkSetDb"][0]["Link"][231]["Id"]=="9847217"
        assert record[0]["LinkSetDb"][0]["Link"][231]["Score"]=="29618439"
        assert record[0]["LinkSetDb"][0]["Link"][232]["Id"]=="17142236"
        assert record[0]["LinkSetDb"][0]["Link"][232]["Score"]=="29577611"
        assert record[0]["LinkSetDb"][0]["Link"][233]["Id"]=="17059604"
        assert record[0]["LinkSetDb"][0]["Link"][233]["Score"]=="29569767"
        assert record[0]["LinkSetDb"][0]["Link"][234]["Id"]=="16845079"
        assert record[0]["LinkSetDb"][0]["Link"][234]["Score"]=="29506663"
        assert record[0]["LinkSetDb"][0]["Link"][235]["Id"]=="14727153"
        assert record[0]["LinkSetDb"][0]["Link"][235]["Score"]=="29368276"
        assert record[0]["LinkSetDb"][0]["Link"][236]["Id"]=="18045498"
        assert record[0]["LinkSetDb"][0]["Link"][236]["Score"]=="29364312"
        assert record[0]["LinkSetDb"][0]["Link"][237]["Id"]=="17185755"
        assert record[0]["LinkSetDb"][0]["Link"][237]["Score"]=="29331905"
        assert record[0]["LinkSetDb"][0]["Link"][238]["Id"]=="18025704"
        assert record[0]["LinkSetDb"][0]["Link"][238]["Score"]=="29323161"
        assert record[0]["LinkSetDb"][0]["Link"][239]["Id"]=="15215374"
        assert record[0]["LinkSetDb"][0]["Link"][239]["Score"]=="29306559"
        assert record[0]["LinkSetDb"][0]["Link"][240]["Id"]=="17135185"
        assert record[0]["LinkSetDb"][0]["Link"][240]["Score"]=="29236297"
        assert record[0]["LinkSetDb"][0]["Link"][241]["Id"]=="10466135"
        assert record[0]["LinkSetDb"][0]["Link"][241]["Score"]=="29231855"
        assert record[0]["LinkSetDb"][0]["Link"][242]["Id"]=="17148475"
        assert record[0]["LinkSetDb"][0]["Link"][242]["Score"]=="29229044"
        assert record[0]["LinkSetDb"][0]["Link"][243]["Id"]=="15657101"
        assert record[0]["LinkSetDb"][0]["Link"][243]["Score"]=="29209567"
        assert record[0]["LinkSetDb"][0]["Link"][244]["Id"]=="14681490"
        assert record[0]["LinkSetDb"][0]["Link"][244]["Score"]=="29189708"
        assert record[0]["LinkSetDb"][0]["Link"][245]["Id"]=="15714328"
        assert record[0]["LinkSetDb"][0]["Link"][245]["Score"]=="29183488"
        assert record[0]["LinkSetDb"][0]["Link"][246]["Id"]=="14960477"
        assert record[0]["LinkSetDb"][0]["Link"][246]["Score"]=="29040531"
        assert record[0]["LinkSetDb"][0]["Link"][247]["Id"]=="11015564"
        assert record[0]["LinkSetDb"][0]["Link"][247]["Score"]=="29011368"
        assert record[0]["LinkSetDb"][0]["Link"][248]["Id"]=="18064491"
        assert record[0]["LinkSetDb"][0]["Link"][248]["Score"]=="28956740"
        assert record[0]["LinkSetDb"][0]["Link"][249]["Id"]=="12734009"
        assert record[0]["LinkSetDb"][0]["Link"][249]["Score"]=="28950064"
        assert record[0]["LinkSetDb"][0]["Link"][250]["Id"]=="17094804"
        assert record[0]["LinkSetDb"][0]["Link"][250]["Score"]=="28906953"
        assert record[0]["LinkSetDb"][0]["Link"][251]["Id"]=="17908294"
        assert record[0]["LinkSetDb"][0]["Link"][251]["Score"]=="28897717"
        assert record[0]["LinkSetDb"][0]["Link"][252]["Id"]=="16176584"
        assert record[0]["LinkSetDb"][0]["Link"][252]["Score"]=="28874470"
        assert record[0]["LinkSetDb"][0]["Link"][253]["Id"]=="14715089"
        assert record[0]["LinkSetDb"][0]["Link"][253]["Score"]=="28763886"
        assert record[0]["LinkSetDb"][0]["Link"][254]["Id"]=="14681408"
        assert record[0]["LinkSetDb"][0]["Link"][254]["Score"]=="28697827"
        assert record[0]["LinkSetDb"][0]["Link"][255]["Id"]=="14594716"
        assert record[0]["LinkSetDb"][0]["Link"][255]["Score"]=="28686075"
        assert record[0]["LinkSetDb"][0]["Link"][256]["Id"]=="16528802"
        assert record[0]["LinkSetDb"][0]["Link"][256]["Score"]=="28644452"
        assert record[0]["LinkSetDb"][0]["Link"][257]["Id"]=="16010002"
        assert record[0]["LinkSetDb"][0]["Link"][257]["Score"]=="28637570"
        assert record[0]["LinkSetDb"][0]["Link"][258]["Id"]=="17430565"
        assert record[0]["LinkSetDb"][0]["Link"][258]["Score"]=="28635513"
        assert record[0]["LinkSetDb"][0]["Link"][259]["Id"]=="16452787"
        assert record[0]["LinkSetDb"][0]["Link"][259]["Score"]=="28631832"
        assert record[0]["LinkSetDb"][0]["Link"][260]["Id"]=="11197127"
        assert record[0]["LinkSetDb"][0]["Link"][260]["Score"]=="28619225"
        assert record[0]["LinkSetDb"][0]["Link"][261]["Id"]=="8682188"
        assert record[0]["LinkSetDb"][0]["Link"][261]["Score"]=="28592521"
        assert record[0]["LinkSetDb"][0]["Link"][262]["Id"]=="12519940"
        assert record[0]["LinkSetDb"][0]["Link"][262]["Score"]=="28573991"
        assert record[0]["LinkSetDb"][0]["Link"][263]["Id"]=="17121775"
        assert record[0]["LinkSetDb"][0]["Link"][263]["Score"]=="28448726"
        assert record[0]["LinkSetDb"][0]["Link"][264]["Id"]=="16371163"
        assert record[0]["LinkSetDb"][0]["Link"][264]["Score"]=="28373394"
        assert record[0]["LinkSetDb"][0]["Link"][265]["Id"]=="15300845"
        assert record[0]["LinkSetDb"][0]["Link"][265]["Score"]=="28338477"
        assert record[0]["LinkSetDb"][0]["Link"][266]["Id"]=="15248903"
        assert record[0]["LinkSetDb"][0]["Link"][266]["Score"]=="28323328"
        assert record[0]["LinkSetDb"][0]["Link"][267]["Id"]=="11319266"
        assert record[0]["LinkSetDb"][0]["Link"][267]["Score"]=="28293166"
        assert record[0]["LinkSetDb"][0]["Link"][268]["Id"]=="16336665"
        assert record[0]["LinkSetDb"][0]["Link"][268]["Score"]=="28231249"
        assert record[0]["LinkSetDb"][0]["Link"][269]["Id"]=="14681350"
        assert record[0]["LinkSetDb"][0]["Link"][269]["Score"]=="28227327"
        assert record[0]["LinkSetDb"][0]["Link"][270]["Id"]=="16216831"
        assert record[0]["LinkSetDb"][0]["Link"][270]["Score"]=="28224610"
        assert record[0]["LinkSetDb"][0]["Link"][271]["Id"]=="15494741"
        assert record[0]["LinkSetDb"][0]["Link"][271]["Score"]=="28190925"
        assert record[0]["LinkSetDb"][0]["Link"][272]["Id"]=="17088289"
        assert record[0]["LinkSetDb"][0]["Link"][272]["Score"]=="28168901"
        assert record[0]["LinkSetDb"][0]["Link"][273]["Id"]=="17099235"
        assert record[0]["LinkSetDb"][0]["Link"][273]["Score"]=="28159766"
        assert record[0]["LinkSetDb"][0]["Link"][274]["Id"]=="15215372"
        assert record[0]["LinkSetDb"][0]["Link"][274]["Score"]=="28129693"
        assert record[0]["LinkSetDb"][0]["Link"][275]["Id"]=="9169870"
        assert record[0]["LinkSetDb"][0]["Link"][275]["Score"]=="28117392"
        assert record[0]["LinkSetDb"][0]["Link"][276]["Id"]=="10077537"
        assert record[0]["LinkSetDb"][0]["Link"][276]["Score"]=="27911205"
        assert record[0]["LinkSetDb"][0]["Link"][277]["Id"]=="18172929"
        assert record[0]["LinkSetDb"][0]["Link"][277]["Score"]=="27885172"
        assert record[0]["LinkSetDb"][0]["Link"][278]["Id"]=="9571806"
        assert record[0]["LinkSetDb"][0]["Link"][278]["Score"]=="27841468"
        assert record[0]["LinkSetDb"][0]["Link"][279]["Id"]=="11752280"
        assert record[0]["LinkSetDb"][0]["Link"][279]["Score"]=="27795833"
        assert record[0]["LinkSetDb"][0]["Link"][280]["Id"]=="11414208"
        assert record[0]["LinkSetDb"][0]["Link"][280]["Score"]=="27725996"
        assert record[0]["LinkSetDb"][0]["Link"][281]["Id"]=="9298642"
        assert record[0]["LinkSetDb"][0]["Link"][281]["Score"]=="27716027"
        assert record[0]["LinkSetDb"][0]["Link"][282]["Id"]=="18073380"
        assert record[0]["LinkSetDb"][0]["Link"][282]["Score"]=="27437383"
        assert record[0]["LinkSetDb"][0]["Link"][283]["Id"]=="14527308"
        assert record[0]["LinkSetDb"][0]["Link"][283]["Score"]=="27332641"
        assert record[0]["LinkSetDb"][0]["Link"][284]["Id"]=="9847220"
        assert record[0]["LinkSetDb"][0]["Link"][284]["Score"]=="27083894"
        assert record[0]["LinkSetDb"][0]["Link"][285]["Id"]=="10413661"
        assert record[0]["LinkSetDb"][0]["Link"][285]["Score"]=="27073030"
        assert record[0]["LinkSetDb"][0]["Link"][286]["Id"]=="10407677"
        assert record[0]["LinkSetDb"][0]["Link"][286]["Score"]=="26907635"
        assert record[0]["LinkSetDb"][0]["Link"][287]["Id"]=="11244060"
        assert record[0]["LinkSetDb"][0]["Link"][287]["Score"]=="26897688"
        assert record[0]["LinkSetDb"][0]["Link"][288]["Id"]=="10227170"
        assert record[0]["LinkSetDb"][0]["Link"][288]["Score"]=="26766431"
        assert record[0]["LinkSetDb"][0]["Link"][289]["Id"]=="8719164"
        assert record[0]["LinkSetDb"][0]["Link"][289]["Score"]=="26515360"
        assert record[0]["LinkSetDb"][0]["Link"][290]["Id"]=="18359019"
        assert record[0]["LinkSetDb"][0]["Link"][290]["Score"]=="26225983"
        assert record[0]["LinkSetDb"][0]["Link"][291]["Id"]=="10511680"
        assert record[0]["LinkSetDb"][0]["Link"][291]["Score"]=="26031196"
        assert record[0]["LinkSetDb"][0]["Link"][292]["Id"]=="9884329"
        assert record[0]["LinkSetDb"][0]["Link"][292]["Score"]=="25992564"
        assert record[0]["LinkSetDb"][0]["Link"][293]["Id"]=="17827295"
        assert record[0]["LinkSetDb"][0]["Link"][293]["Score"]=="25989152"
        assert record[0]["LinkSetDb"][0]["Link"][294]["Id"]=="10899154"
        assert record[0]["LinkSetDb"][0]["Link"][294]["Score"]=="25843128"
        assert record[0]["LinkSetDb"][0]["Link"][295]["Id"]=="11668619"
        assert record[0]["LinkSetDb"][0]["Link"][295]["Score"]=="25822950"
        assert record[0]["LinkSetDb"][0]["Link"][296]["Id"]=="18386064"
        assert record[0]["LinkSetDb"][0]["Link"][296]["Score"]=="25702942"
        assert record[0]["LinkSetDb"][0]["Link"][297]["Id"]=="11092731"
        assert record[0]["LinkSetDb"][0]["Link"][297]["Score"]=="25618899"
        assert record[0]["LinkSetDb"][0]["Link"][298]["Id"]=="9520376"
        assert record[0]["LinkSetDb"][0]["Link"][298]["Score"]=="25549761"
        assert record[0]["LinkSetDb"][0]["Link"][299]["Id"]=="11756688"
        assert record[0]["LinkSetDb"][0]["Link"][299]["Score"]=="25440634"
        assert record[0]["LinkSetDb"][0]["Link"][300]["Id"]=="10737802"
        assert record[0]["LinkSetDb"][0]["Link"][300]["Score"]=="25362744"
        assert record[0]["LinkSetDb"][0]["Link"][301]["Id"]=="9879937"
        assert record[0]["LinkSetDb"][0]["Link"][301]["Score"]=="25277089"
        assert record[0]["LinkSetDb"][0]["Link"][302]["Id"]=="17822801"
        assert record[0]["LinkSetDb"][0]["Link"][302]["Score"]=="25252984"
        assert record[0]["LinkSetDb"][0]["Link"][303]["Id"]=="10965872"
        assert record[0]["LinkSetDb"][0]["Link"][303]["Score"]=="25208185"
        assert record[0]["LinkSetDb"][0]["Link"][304]["Id"]=="10511682"
        assert record[0]["LinkSetDb"][0]["Link"][304]["Score"]=="25183443"
        assert record[0]["LinkSetDb"][0]["Link"][305]["Id"]=="10851186"
        assert record[0]["LinkSetDb"][0]["Link"][305]["Score"]=="25092764"
        assert record[0]["LinkSetDb"][0]["Link"][306]["Id"]=="9775388"
        assert record[0]["LinkSetDb"][0]["Link"][306]["Score"]=="25026910"
        assert record[0]["LinkSetDb"][0]["Link"][307]["Id"]=="10810023"
        assert record[0]["LinkSetDb"][0]["Link"][307]["Score"]=="24904718"
        assert record[0]["LinkSetDb"][0]["Link"][308]["Id"]=="18032438"
        assert record[0]["LinkSetDb"][0]["Link"][308]["Score"]=="24509777"
        assert record[0]["LinkSetDb"][0]["Link"][309]["Id"]=="18377816"
        assert record[0]["LinkSetDb"][0]["Link"][309]["Score"]=="24373788"
        assert record[0]["LinkSetDb"][0]["Link"][310]["Id"]=="11774190"
        assert record[0]["LinkSetDb"][0]["Link"][310]["Score"]=="24185658"
        assert record[0]["LinkSetDb"][0]["Link"][311]["Id"]=="10484179"
        assert record[0]["LinkSetDb"][0]["Link"][311]["Score"]=="24122767"
        assert record[0]["LinkSetDb"][0]["Link"][312]["Id"]=="9625791"
        assert record[0]["LinkSetDb"][0]["Link"][312]["Score"]=="24049917"
        assert record[0]["LinkSetDb"][0]["Link"][313]["Id"]=="11446511"
        assert record[0]["LinkSetDb"][0]["Link"][313]["Score"]=="24048253"
        assert record[0]["LinkSetDb"][0]["Link"][314]["Id"]=="10066467"
        assert record[0]["LinkSetDb"][0]["Link"][314]["Score"]=="23968405"
        assert record[0]["LinkSetDb"][0]["Link"][315]["Id"]=="11783003"
        assert record[0]["LinkSetDb"][0]["Link"][315]["Score"]=="23393870"
        assert record[0]["LinkSetDb"][0]["Link"][316]["Id"]=="10611059"
        assert record[0]["LinkSetDb"][0]["Link"][316]["Score"]=="23255298"
        assert record[0]["LinkSetDb"][0]["Link"][317]["Id"]=="10587943"
        assert record[0]["LinkSetDb"][0]["Link"][317]["Score"]=="23014503"
        assert record[0]["LinkSetDb"][0]["Link"][318]["Id"]=="10612820"
        assert record[0]["LinkSetDb"][0]["Link"][318]["Score"]=="22990878"
        assert record[0]["LinkSetDb"][0]["Link"][319]["Id"]=="9685316"
        assert record[0]["LinkSetDb"][0]["Link"][319]["Score"]=="22771348"
        assert record[0]["LinkSetDb"][0]["Link"][320]["Id"]=="11125121"
        assert record[0]["LinkSetDb"][0]["Link"][320]["Score"]=="22732820"
        assert record[0]["LinkSetDb"][0]["Link"][321]["Id"]=="10075567"
        assert record[0]["LinkSetDb"][0]["Link"][321]["Score"]=="22670427"
        assert record[0]["LinkSetDb"][0]["Link"][322]["Id"]=="11084929"
        assert record[0]["LinkSetDb"][0]["Link"][322]["Score"]=="22397665"
        assert record[0]["LinkSetDb"][0]["Link"][323]["Id"]=="11357826"
        assert record[0]["LinkSetDb"][0]["Link"][323]["Score"]=="22362882"
        assert record[0]["LinkSetDb"][0]["Link"][324]["Id"]=="17983575"
        assert record[0]["LinkSetDb"][0]["Link"][324]["Score"]=="22305320"
        assert record[0]["LinkSetDb"][0]["Link"][325]["Id"]=="11038308"
        assert record[0]["LinkSetDb"][0]["Link"][325]["Score"]=="22115670"
        assert record[0]["LinkSetDb"][0]["Link"][326]["Id"]=="18257289"
        assert record[0]["LinkSetDb"][0]["Link"][326]["Score"]=="22053176"
        assert record[0]["LinkSetDb"][0]["Link"][327]["Id"]=="10419978"
        assert record[0]["LinkSetDb"][0]["Link"][327]["Score"]=="22016184"
        assert record[0]["LinkSetDb"][0]["Link"][328]["Id"]=="9421619"
        assert record[0]["LinkSetDb"][0]["Link"][328]["Score"]=="21957407"
        assert record[0]["LinkSetDb"][0]["Link"][329]["Id"]=="10592198"
        assert record[0]["LinkSetDb"][0]["Link"][329]["Score"]=="21803908"
        assert record[0]["LinkSetDb"][0]["Link"][330]["Id"]=="11483982"
        assert record[0]["LinkSetDb"][0]["Link"][330]["Score"]=="20783817"
        assert record[0]["LinkSetDb"][0]["Link"][331]["Id"]=="11329386"
        assert record[0]["LinkSetDb"][0]["Link"][331]["Score"]=="20223493"
        assert record[0]["LinkSetDb"][0]["Link"][332]["Id"]=="10587942"
        assert record[0]["LinkSetDb"][0]["Link"][332]["Score"]=="20208799"
        assert record[0]["LinkSetDb"][0]["Link"][333]["Id"]=="10810024"
        assert record[0]["LinkSetDb"][0]["Link"][333]["Score"]=="19989188"
        assert record[0]["LinkSetDb"][0]["Link"][334]["Id"]=="11480780"
        assert record[0]["LinkSetDb"][0]["Link"][334]["Score"]=="19974101"
        assert record[0]["LinkSetDb"][0]["Link"][335]["Id"]=="11802378"
        assert record[0]["LinkSetDb"][0]["Link"][335]["Score"]=="19738532"
        assert record[0]["LinkSetDb"][0]["Link"][336]["Id"]=="10610803"
        assert record[0]["LinkSetDb"][0]["Link"][336]["Score"]=="19359100"
        assert record[0]["LinkSetDb"][0]["Link"][337]["Id"]=="10407668"
        assert record[0]["LinkSetDb"][0]["Link"][337]["Score"]=="19070525"
        assert record[0]["LinkSetDb"][0]["Link"][338]["Id"]=="18287701"
        assert record[0]["LinkSetDb"][0]["Link"][338]["Score"]=="19065945"
        assert record[0]["LinkSetDb"][0]["Link"][339]["Id"]=="10963611"
        assert record[0]["LinkSetDb"][0]["Link"][339]["Score"]=="18962273"
        assert record[0]["LinkSetDb"][0]["Link"][340]["Id"]=="10447503"
        assert record[0]["LinkSetDb"][0]["Link"][340]["Score"]=="17406980"
        assert record[0]["LinkSetDb"][0]["Link"][341]["Id"]=="9830540"
        assert record[0]["LinkSetDb"][0]["Link"][341]["Score"]=="17143709"
        assert record[0]["LinkSetDb"][0]["Link"][342]["Id"]=="11462837"
        assert record[0]["LinkSetDb"][0]["Link"][342]["Score"]=="16819799"
        assert record[0]["LinkSetDb"][0]["Link"][343]["Id"]=="10637631"
        assert record[0]["LinkSetDb"][0]["Link"][343]["Score"]=="16390796"
        assert record[0]["LinkSetDb"][0]["Link"][344]["Id"]=="11387032"
        assert record[0]["LinkSetDb"][0]["Link"][344]["Score"]=="15698695"
        assert record[0]["LinkSetDb"][0]["Link"][345]["Id"]=="18365535"
        assert record[0]["LinkSetDb"][0]["Link"][345]["Score"]=="15494816"
        assert record[0]["LinkSetDb"][0]["Link"][346]["Id"]=="15181901"
        assert record[0]["LinkSetDb"][0]["Link"][346]["Score"]=="14385628"

    def t_medline(self):
        '''Test parsing medline indexed articles returned by ELink
        '''
        # Retrieve MEDLINE indexed only related articles for PMID 12242737
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12242737", db="pubmed",
        #                      term="medline[sb]")
        input = open('Entrez/elink4.xml')
        record = Entrez.read(input)
        
        assert len(record)==1
        assert record[0]["DbFrom"]=="pubmed"
        assert record[0]["IdList"]==["12242737"]
        assert record[0]["LinkSetDb"][0]["DbTo"]=="pubmed"
        assert record[0]["LinkSetDb"][0]["LinkName"]=="pubmed_pubmed"
        assert record[0]["LinkSetDb"][0]["Link"][0]["Id"]=="12242737"
        assert record[0]["LinkSetDb"][0]["Link"][0]["Score"]=="2147483647"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Id"]=="11218011"
        assert record[0]["LinkSetDb"][0]["Link"][1]["Score"]=="50825961"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Id"]=="11329656"
        assert record[0]["LinkSetDb"][0]["Link"][2]["Score"]=="49822043"
        assert record[0]["LinkSetDb"][0]["Link"][3]["Id"]=="9757294"
        assert record[0]["LinkSetDb"][0]["Link"][3]["Score"]=="42645380"
        assert record[0]["LinkSetDb"][0]["Link"][4]["Id"]=="9456947"
        assert record[0]["LinkSetDb"][0]["Link"][4]["Score"]=="39871666"
        assert record[0]["LinkSetDb"][0]["Link"][5]["Id"]=="17193860"
        assert record[0]["LinkSetDb"][0]["Link"][5]["Score"]=="39717388"
        assert record[0]["LinkSetDb"][0]["Link"][6]["Id"]=="11274884"
        assert record[0]["LinkSetDb"][0]["Link"][6]["Score"]=="39233276"
        assert record[0]["LinkSetDb"][0]["Link"][7]["Id"]=="12878072"
        assert record[0]["LinkSetDb"][0]["Link"][7]["Score"]=="37748327"
        assert record[0]["LinkSetDb"][0]["Link"][8]["Id"]=="11125632"
        assert record[0]["LinkSetDb"][0]["Link"][8]["Score"]=="36227857"
        assert record[0]["LinkSetDb"][0]["Link"][9]["Id"]=="12822521"
        assert record[0]["LinkSetDb"][0]["Link"][9]["Score"]=="36170366"
        assert record[0]["LinkSetDb"][0]["Link"][10]["Id"]=="16999328"
        assert record[0]["LinkSetDb"][0]["Link"][10]["Score"]=="36107139"
        assert record[0]["LinkSetDb"][0]["Link"][11]["Id"]=="17875142"
        assert record[0]["LinkSetDb"][0]["Link"][11]["Score"]=="35736802"
        assert record[0]["LinkSetDb"][0]["Link"][12]["Id"]=="9510579"
        assert record[0]["LinkSetDb"][0]["Link"][12]["Score"]=="35206779"
        assert record[0]["LinkSetDb"][0]["Link"][13]["Id"]=="17354190"
        assert record[0]["LinkSetDb"][0]["Link"][13]["Score"]=="34792954"
        assert record[0]["LinkSetDb"][0]["Link"][14]["Id"]=="11702119"
        assert record[0]["LinkSetDb"][0]["Link"][14]["Score"]=="34618984"
        assert record[0]["LinkSetDb"][0]["Link"][15]["Id"]=="10024396"
        assert record[0]["LinkSetDb"][0]["Link"][15]["Score"]=="33877753"
        assert record[0]["LinkSetDb"][0]["Link"][16]["Id"]=="14650118"
        assert record[0]["LinkSetDb"][0]["Link"][16]["Score"]=="33746160"
        assert record[0]["LinkSetDb"][0]["Link"][17]["Id"]=="17243036"
        assert record[0]["LinkSetDb"][0]["Link"][17]["Score"]=="33198930"
        assert record[0]["LinkSetDb"][0]["Link"][18]["Id"]=="16580806"
        assert record[0]["LinkSetDb"][0]["Link"][18]["Score"]=="33117197"
        assert record[0]["LinkSetDb"][0]["Link"][19]["Id"]=="15278705"
        assert record[0]["LinkSetDb"][0]["Link"][19]["Score"]=="33002826"
        assert record[0]["LinkSetDb"][0]["Link"][20]["Id"]=="15236131"
        assert record[0]["LinkSetDb"][0]["Link"][20]["Score"]=="32808406"
        assert record[0]["LinkSetDb"][0]["Link"][21]["Id"]=="11368937"
        assert record[0]["LinkSetDb"][0]["Link"][21]["Score"]=="32277701"
        assert record[0]["LinkSetDb"][0]["Link"][22]["Id"]=="10688065"
        assert record[0]["LinkSetDb"][0]["Link"][22]["Score"]=="32052850"
        assert record[0]["LinkSetDb"][0]["Link"][23]["Id"]=="15635471"
        assert record[0]["LinkSetDb"][0]["Link"][23]["Score"]=="31938251"
        assert record[0]["LinkSetDb"][0]["Link"][24]["Id"]=="16357381"
        assert record[0]["LinkSetDb"][0]["Link"][24]["Score"]=="31780147"
        assert record[0]["LinkSetDb"][0]["Link"][25]["Id"]=="8153333"
        assert record[0]["LinkSetDb"][0]["Link"][25]["Score"]=="31542202"
        assert record[0]["LinkSetDb"][0]["Link"][26]["Id"]=="16284132"
        assert record[0]["LinkSetDb"][0]["Link"][26]["Score"]=="31290577"
        assert record[0]["LinkSetDb"][0]["Link"][27]["Id"]=="11329162"
        assert record[0]["LinkSetDb"][0]["Link"][27]["Score"]=="31163088"
        assert record[0]["LinkSetDb"][0]["Link"][28]["Id"]=="11973040"
        assert record[0]["LinkSetDb"][0]["Link"][28]["Score"]=="31156707"
        assert record[0]["LinkSetDb"][0]["Link"][29]["Id"]=="15143223"
        assert record[0]["LinkSetDb"][0]["Link"][29]["Score"]=="31025329"
        assert record[0]["LinkSetDb"][0]["Link"][30]["Id"]=="17040637"
        assert record[0]["LinkSetDb"][0]["Link"][30]["Score"]=="30990506"
        assert record[0]["LinkSetDb"][0]["Link"][31]["Id"]=="11016058"
        assert record[0]["LinkSetDb"][0]["Link"][31]["Score"]=="30966482"
        assert record[0]["LinkSetDb"][0]["Link"][32]["Id"]=="9317094"
        assert record[0]["LinkSetDb"][0]["Link"][32]["Score"]=="30935529"
        assert record[0]["LinkSetDb"][0]["Link"][33]["Id"]=="16133609"
        assert record[0]["LinkSetDb"][0]["Link"][33]["Score"]=="30580027"
        assert record[0]["LinkSetDb"][0]["Link"][34]["Id"]=="17325998"
        assert record[0]["LinkSetDb"][0]["Link"][34]["Score"]=="30130533"
        assert record[0]["LinkSetDb"][0]["Link"][35]["Id"]=="15505294"
        assert record[0]["LinkSetDb"][0]["Link"][35]["Score"]=="29430378"
        assert record[0]["LinkSetDb"][0]["Link"][36]["Id"]=="17268692"
        assert record[0]["LinkSetDb"][0]["Link"][36]["Score"]=="29166153"
        assert record[0]["LinkSetDb"][0]["Link"][37]["Id"]=="11329655"
        assert record[0]["LinkSetDb"][0]["Link"][37]["Score"]=="29112282"
        assert record[0]["LinkSetDb"][0]["Link"][38]["Id"]=="11775722"
        assert record[0]["LinkSetDb"][0]["Link"][38]["Score"]=="28940754"
        assert record[0]["LinkSetDb"][0]["Link"][39]["Id"]=="11907356"
        assert record[0]["LinkSetDb"][0]["Link"][39]["Score"]=="28860163"
        assert record[0]["LinkSetDb"][0]["Link"][40]["Id"]=="10222515"
        assert record[0]["LinkSetDb"][0]["Link"][40]["Score"]=="28807143"
        assert record[0]["LinkSetDb"][0]["Link"][41]["Id"]=="17174054"
        assert record[0]["LinkSetDb"][0]["Link"][41]["Score"]=="28790302"
        assert record[0]["LinkSetDb"][0]["Link"][42]["Id"]=="9314960"
        assert record[0]["LinkSetDb"][0]["Link"][42]["Score"]=="28750160"
        assert record[0]["LinkSetDb"][0]["Link"][43]["Id"]=="14661661"
        assert record[0]["LinkSetDb"][0]["Link"][43]["Score"]=="28361423"
        assert record[0]["LinkSetDb"][0]["Link"][44]["Id"]=="17879696"
        assert record[0]["LinkSetDb"][0]["Link"][44]["Score"]=="28120568"
        assert record[0]["LinkSetDb"][0]["Link"][45]["Id"]=="4818442"
        assert record[0]["LinkSetDb"][0]["Link"][45]["Score"]=="28058957"
        assert record[0]["LinkSetDb"][0]["Link"][46]["Id"]=="15141648"
        assert record[0]["LinkSetDb"][0]["Link"][46]["Score"]=="28011681"
        assert record[0]["LinkSetDb"][0]["Link"][47]["Id"]=="8855688"
        assert record[0]["LinkSetDb"][0]["Link"][47]["Score"]=="27711822"
        assert record[0]["LinkSetDb"][0]["Link"][48]["Id"]=="17875143"
        assert record[0]["LinkSetDb"][0]["Link"][48]["Score"]=="27711025"
        assert record[0]["LinkSetDb"][0]["Link"][49]["Id"]=="1481295"
        assert record[0]["LinkSetDb"][0]["Link"][49]["Score"]=="27707751"
        assert record[0]["LinkSetDb"][0]["Link"][50]["Id"]=="8599783"
        assert record[0]["LinkSetDb"][0]["Link"][50]["Score"]=="27683273"
        assert record[0]["LinkSetDb"][0]["Link"][51]["Id"]=="10499696"
        assert record[0]["LinkSetDb"][0]["Link"][51]["Score"]=="27623848"
        assert record[0]["LinkSetDb"][0]["Link"][52]["Id"]=="12733684"
        assert record[0]["LinkSetDb"][0]["Link"][52]["Score"]=="27527242"
        assert record[0]["LinkSetDb"][0]["Link"][53]["Id"]=="18021675"
        assert record[0]["LinkSetDb"][0]["Link"][53]["Score"]=="27495074"
        assert record[0]["LinkSetDb"][0]["Link"][54]["Id"]=="12226761"
        assert record[0]["LinkSetDb"][0]["Link"][54]["Score"]=="27366064"
        assert record[0]["LinkSetDb"][0]["Link"][55]["Id"]=="4808999"
        assert record[0]["LinkSetDb"][0]["Link"][55]["Score"]=="27304472"
        assert record[0]["LinkSetDb"][0]["Link"][56]["Id"]=="16988291"
        assert record[0]["LinkSetDb"][0]["Link"][56]["Score"]=="27295295"
        assert record[0]["LinkSetDb"][0]["Link"][57]["Id"]=="10575758"
        assert record[0]["LinkSetDb"][0]["Link"][57]["Score"]=="27243181"
        assert record[0]["LinkSetDb"][0]["Link"][58]["Id"]=="8903064"
        assert record[0]["LinkSetDb"][0]["Link"][58]["Score"]=="27206664"
        assert record[0]["LinkSetDb"][0]["Link"][59]["Id"]=="10811354"
        assert record[0]["LinkSetDb"][0]["Link"][59]["Score"]=="27088219"
        assert record[0]["LinkSetDb"][0]["Link"][60]["Id"]=="16096604"
        assert record[0]["LinkSetDb"][0]["Link"][60]["Score"]=="26862979"
        assert record[0]["LinkSetDb"][0]["Link"][61]["Id"]=="15788584"
        assert record[0]["LinkSetDb"][0]["Link"][61]["Score"]=="26759584"
        assert record[0]["LinkSetDb"][0]["Link"][62]["Id"]=="17376366"
        assert record[0]["LinkSetDb"][0]["Link"][62]["Score"]=="26743241"
        assert record[0]["LinkSetDb"][0]["Link"][63]["Id"]=="16566645"
        assert record[0]["LinkSetDb"][0]["Link"][63]["Score"]=="26725076"
        assert record[0]["LinkSetDb"][0]["Link"][64]["Id"]=="17259035"
        assert record[0]["LinkSetDb"][0]["Link"][64]["Score"]=="26595433"
        assert record[0]["LinkSetDb"][0]["Link"][65]["Id"]=="9314959"
        assert record[0]["LinkSetDb"][0]["Link"][65]["Score"]=="26445900"
        assert record[0]["LinkSetDb"][0]["Link"][66]["Id"]=="11895298"
        assert record[0]["LinkSetDb"][0]["Link"][66]["Score"]=="26256774"
        assert record[0]["LinkSetDb"][0]["Link"][67]["Id"]=="11740602"
        assert record[0]["LinkSetDb"][0]["Link"][67]["Score"]=="26158189"
        assert record[0]["LinkSetDb"][0]["Link"][68]["Id"]=="15022983"
        assert record[0]["LinkSetDb"][0]["Link"][68]["Score"]=="25889186"
        assert record[0]["LinkSetDb"][0]["Link"][69]["Id"]=="15300544"
        assert record[0]["LinkSetDb"][0]["Link"][69]["Score"]=="25837458"
        assert record[0]["LinkSetDb"][0]["Link"][70]["Id"]=="12719915"
        assert record[0]["LinkSetDb"][0]["Link"][70]["Score"]=="25831232"
        assert record[0]["LinkSetDb"][0]["Link"][71]["Id"]=="14661306"
        assert record[0]["LinkSetDb"][0]["Link"][71]["Score"]=="25788023"
        assert record[0]["LinkSetDb"][0]["Link"][72]["Id"]=="16362812"
        assert record[0]["LinkSetDb"][0]["Link"][72]["Score"]=="25565076"
        assert record[0]["LinkSetDb"][0]["Link"][73]["Id"]=="17320773"
        assert record[0]["LinkSetDb"][0]["Link"][73]["Score"]=="25504305"
        assert record[0]["LinkSetDb"][0]["Link"][74]["Id"]=="11762248"
        assert record[0]["LinkSetDb"][0]["Link"][74]["Score"]=="25504002"
        assert record[0]["LinkSetDb"][0]["Link"][75]["Id"]=="10665303"
        assert record[0]["LinkSetDb"][0]["Link"][75]["Score"]=="25384388"
        assert record[0]["LinkSetDb"][0]["Link"][76]["Id"]=="17453494"
        assert record[0]["LinkSetDb"][0]["Link"][76]["Score"]=="25226372"
        assert record[0]["LinkSetDb"][0]["Link"][77]["Id"]=="9575723"
        assert record[0]["LinkSetDb"][0]["Link"][77]["Score"]=="25174136"
        assert record[0]["LinkSetDb"][0]["Link"][78]["Id"]=="12744498"
        assert record[0]["LinkSetDb"][0]["Link"][78]["Score"]=="24971179"
        assert record[0]["LinkSetDb"][0]["Link"][79]["Id"]=="12352163"
        assert record[0]["LinkSetDb"][0]["Link"][79]["Score"]=="24915990"
        assert record[0]["LinkSetDb"][0]["Link"][80]["Id"]=="8290724"
        assert record[0]["LinkSetDb"][0]["Link"][80]["Score"]=="24909462"
        assert record[0]["LinkSetDb"][0]["Link"][81]["Id"]=="11973504"
        assert record[0]["LinkSetDb"][0]["Link"][81]["Score"]=="24878058"
        assert record[0]["LinkSetDb"][0]["Link"][82]["Id"]=="14661668"
        assert record[0]["LinkSetDb"][0]["Link"][82]["Score"]=="24779779"
        assert record[0]["LinkSetDb"][0]["Link"][83]["Id"]=="16552382"
        assert record[0]["LinkSetDb"][0]["Link"][83]["Score"]=="24760919"
        assert record[0]["LinkSetDb"][0]["Link"][84]["Id"]=="17709829"
        assert record[0]["LinkSetDb"][0]["Link"][84]["Score"]=="24743292"
        assert record[0]["LinkSetDb"][0]["Link"][85]["Id"]=="14528718"
        assert record[0]["LinkSetDb"][0]["Link"][85]["Score"]=="24686212"
        assert record[0]["LinkSetDb"][0]["Link"][86]["Id"]=="15008163"
        assert record[0]["LinkSetDb"][0]["Link"][86]["Score"]=="24612994"
        assert record[0]["LinkSetDb"][0]["Link"][87]["Id"]=="10051883"
        assert record[0]["LinkSetDb"][0]["Link"][87]["Score"]=="24492331"
        assert record[0]["LinkSetDb"][0]["Link"][88]["Id"]=="11027076"
        assert record[0]["LinkSetDb"][0]["Link"][88]["Score"]=="24410525"
        assert record[0]["LinkSetDb"][0]["Link"][89]["Id"]=="17543650"
        assert record[0]["LinkSetDb"][0]["Link"][89]["Score"]=="24371825"
        assert record[0]["LinkSetDb"][0]["Link"][90]["Id"]=="17658095"
        assert record[0]["LinkSetDb"][0]["Link"][90]["Score"]=="24331965"
        assert record[0]["LinkSetDb"][0]["Link"][91]["Id"]=="9193407"
        assert record[0]["LinkSetDb"][0]["Link"][91]["Score"]=="24240252"
        assert record[0]["LinkSetDb"][0]["Link"][92]["Id"]=="10578418"
        assert record[0]["LinkSetDb"][0]["Link"][92]["Score"]=="24091226"
        assert record[0]["LinkSetDb"][0]["Link"][93]["Id"]=="12592155"
        assert record[0]["LinkSetDb"][0]["Link"][93]["Score"]=="24001341"
        assert record[0]["LinkSetDb"][0]["Link"][94]["Id"]=="17157468"
        assert record[0]["LinkSetDb"][0]["Link"][94]["Score"]=="23984321"
        assert record[0]["LinkSetDb"][0]["Link"][95]["Id"]=="15094630"
        assert record[0]["LinkSetDb"][0]["Link"][95]["Score"]=="23912874"
        assert record[0]["LinkSetDb"][0]["Link"][96]["Id"]=="8794574"
        assert record[0]["LinkSetDb"][0]["Link"][96]["Score"]=="23900764"
        assert record[0]["LinkSetDb"][0]["Link"][97]["Id"]=="9125660"
        assert record[0]["LinkSetDb"][0]["Link"][97]["Score"]=="23884352"
        assert record[0]["LinkSetDb"][0]["Link"][98]["Id"]=="8819381"
        assert record[0]["LinkSetDb"][0]["Link"][98]["Score"]=="23839719"
        assert record[0]["LinkSetDb"][0]["Link"][99]["Id"]=="14661666"
        assert record[0]["LinkSetDb"][0]["Link"][99]["Score"]=="23748510"
        assert record[0]["LinkSetDb"][0]["Link"][100]["Id"]=="9658901"
        assert record[0]["LinkSetDb"][0]["Link"][100]["Score"]=="23667126"
        assert record[0]["LinkSetDb"][0]["Link"][101]["Id"]=="12744499"
        assert record[0]["LinkSetDb"][0]["Link"][101]["Score"]=="23647189"
        assert record[0]["LinkSetDb"][0]["Link"][102]["Id"]=="12164574"
        assert record[0]["LinkSetDb"][0]["Link"][102]["Score"]=="23623853"
        assert record[0]["LinkSetDb"][0]["Link"][103]["Id"]=="15136027"
        assert record[0]["LinkSetDb"][0]["Link"][103]["Score"]=="23572558"
        assert record[0]["LinkSetDb"][0]["Link"][104]["Id"]=="14872380"
        assert record[0]["LinkSetDb"][0]["Link"][104]["Score"]=="23460906"
        assert record[0]["LinkSetDb"][0]["Link"][105]["Id"]=="3905087"
        assert record[0]["LinkSetDb"][0]["Link"][105]["Score"]=="23305022"
        assert record[0]["LinkSetDb"][0]["Link"][106]["Id"]=="15642291"
        assert record[0]["LinkSetDb"][0]["Link"][106]["Score"]=="23234831"
        assert record[0]["LinkSetDb"][0]["Link"][107]["Id"]=="16928974"
        assert record[0]["LinkSetDb"][0]["Link"][107]["Score"]=="23223298"
        assert record[0]["LinkSetDb"][0]["Link"][108]["Id"]=="6072516"
        assert record[0]["LinkSetDb"][0]["Link"][108]["Score"]=="23042548"
        assert record[0]["LinkSetDb"][0]["Link"][109]["Id"]=="12949462"
        assert record[0]["LinkSetDb"][0]["Link"][109]["Score"]=="23001441"
        assert record[0]["LinkSetDb"][0]["Link"][110]["Id"]=="10761553"
        assert record[0]["LinkSetDb"][0]["Link"][110]["Score"]=="22995991"
        assert record[0]["LinkSetDb"][0]["Link"][111]["Id"]=="14661663"
        assert record[0]["LinkSetDb"][0]["Link"][111]["Score"]=="22986720"
        assert record[0]["LinkSetDb"][0]["Link"][112]["Id"]=="16338316"
        assert record[0]["LinkSetDb"][0]["Link"][112]["Score"]=="22933288"
        assert record[0]["LinkSetDb"][0]["Link"][113]["Id"]=="17464254"
        assert record[0]["LinkSetDb"][0]["Link"][113]["Score"]=="22912253"
        assert record[0]["LinkSetDb"][0]["Link"][114]["Id"]=="15529836"
        assert record[0]["LinkSetDb"][0]["Link"][114]["Score"]=="22892154"
        assert record[0]["LinkSetDb"][0]["Link"][115]["Id"]=="12361530"
        assert record[0]["LinkSetDb"][0]["Link"][115]["Score"]=="22871698"
        assert record[0]["LinkSetDb"][0]["Link"][116]["Id"]=="12876813"
        assert record[0]["LinkSetDb"][0]["Link"][116]["Score"]=="22822515"
        assert record[0]["LinkSetDb"][0]["Link"][117]["Id"]=="10749221"
        assert record[0]["LinkSetDb"][0]["Link"][117]["Score"]=="22794373"
        assert record[0]["LinkSetDb"][0]["Link"][118]["Id"]=="6482054"
        assert record[0]["LinkSetDb"][0]["Link"][118]["Score"]=="22791927"
        assert record[0]["LinkSetDb"][0]["Link"][119]["Id"]=="9016217"
        assert record[0]["LinkSetDb"][0]["Link"][119]["Score"]=="22738432"
        assert record[0]["LinkSetDb"][0]["Link"][120]["Id"]=="14702442"
        assert record[0]["LinkSetDb"][0]["Link"][120]["Score"]=="22722123"
        assert record[0]["LinkSetDb"][0]["Link"][121]["Id"]=="15279747"
        assert record[0]["LinkSetDb"][0]["Link"][121]["Score"]=="22698787"
        assert record[0]["LinkSetDb"][0]["Link"][122]["Id"]=="7892443"
        assert record[0]["LinkSetDb"][0]["Link"][122]["Score"]=="22642038"
        assert record[0]["LinkSetDb"][0]["Link"][123]["Id"]=="616459"
        assert record[0]["LinkSetDb"][0]["Link"][123]["Score"]=="22591277"
        assert record[0]["LinkSetDb"][0]["Link"][124]["Id"]=="8886718"
        assert record[0]["LinkSetDb"][0]["Link"][124]["Score"]=="22542938"
        assert record[0]["LinkSetDb"][0]["Link"][125]["Id"]=="17245521"
        assert record[0]["LinkSetDb"][0]["Link"][125]["Score"]=="22538649"
        assert record[0]["LinkSetDb"][0]["Link"][126]["Id"]=="1535863"
        assert record[0]["LinkSetDb"][0]["Link"][126]["Score"]=="22468774"
        assert record[0]["LinkSetDb"][0]["Link"][127]["Id"]=="15537403"
        assert record[0]["LinkSetDb"][0]["Link"][127]["Score"]=="22458002"
        assert record[0]["LinkSetDb"][0]["Link"][128]["Id"]=="16040910"
        assert record[0]["LinkSetDb"][0]["Link"][128]["Score"]=="22452119"
        assert record[0]["LinkSetDb"][0]["Link"][129]["Id"]=="16929028"
        assert record[0]["LinkSetDb"][0]["Link"][129]["Score"]=="22433988"
        assert record[0]["LinkSetDb"][0]["Link"][130]["Id"]=="16697589"
        assert record[0]["LinkSetDb"][0]["Link"][130]["Score"]=="22366606"
        assert record[0]["LinkSetDb"][0]["Link"][131]["Id"]=="531835"
        assert record[0]["LinkSetDb"][0]["Link"][131]["Score"]=="22366454"
        assert record[0]["LinkSetDb"][0]["Link"][132]["Id"]=="2308313"
        assert record[0]["LinkSetDb"][0]["Link"][132]["Score"]=="22330898"
        assert record[0]["LinkSetDb"][0]["Link"][133]["Id"]=="12522920"
        assert record[0]["LinkSetDb"][0]["Link"][133]["Score"]=="22178764"
        assert record[0]["LinkSetDb"][0]["Link"][134]["Id"]=="10222521"
        assert record[0]["LinkSetDb"][0]["Link"][134]["Score"]=="22135023"
        assert record[0]["LinkSetDb"][0]["Link"][135]["Id"]=="10499697"
        assert record[0]["LinkSetDb"][0]["Link"][135]["Score"]=="22130302"
        assert record[0]["LinkSetDb"][0]["Link"][136]["Id"]=="8903058"
        assert record[0]["LinkSetDb"][0]["Link"][136]["Score"]=="22113132"
        assert record[0]["LinkSetDb"][0]["Link"][137]["Id"]=="17441569"
        assert record[0]["LinkSetDb"][0]["Link"][137]["Score"]=="22085858"
        assert record[0]["LinkSetDb"][0]["Link"][138]["Id"]=="15284932"
        assert record[0]["LinkSetDb"][0]["Link"][138]["Score"]=="22075791"
        assert record[0]["LinkSetDb"][0]["Link"][139]["Id"]=="15466771"
        assert record[0]["LinkSetDb"][0]["Link"][139]["Score"]=="22075418"
        assert record[0]["LinkSetDb"][0]["Link"][140]["Id"]=="17145267"
        assert record[0]["LinkSetDb"][0]["Link"][140]["Score"]=="22033864"
        assert record[0]["LinkSetDb"][0]["Link"][141]["Id"]=="11329662"
        assert record[0]["LinkSetDb"][0]["Link"][141]["Score"]=="22012948"
        assert record[0]["LinkSetDb"][0]["Link"][142]["Id"]=="10222514"
        assert record[0]["LinkSetDb"][0]["Link"][142]["Score"]=="22009777"
        assert record[0]["LinkSetDb"][0]["Link"][143]["Id"]=="17383530"
        assert record[0]["LinkSetDb"][0]["Link"][143]["Score"]=="22003600"
        assert record[0]["LinkSetDb"][0]["Link"][144]["Id"]=="12455800"
        assert record[0]["LinkSetDb"][0]["Link"][144]["Score"]=="21992674"
        assert record[0]["LinkSetDb"][0]["Link"][145]["Id"]=="15845051"
        assert record[0]["LinkSetDb"][0]["Link"][145]["Score"]=="21946257"
        assert record[0]["LinkSetDb"][0]["Link"][146]["Id"]=="11443295"
        assert record[0]["LinkSetDb"][0]["Link"][146]["Score"]=="21908841"
        assert record[0]["LinkSetDb"][0]["Link"][147]["Id"]=="15162233"
        assert record[0]["LinkSetDb"][0]["Link"][147]["Score"]=="21903624"
        assert record[0]["LinkSetDb"][0]["Link"][148]["Id"]=="16133610"
        assert record[0]["LinkSetDb"][0]["Link"][148]["Score"]=="21872203"
        assert record[0]["LinkSetDb"][0]["Link"][149]["Id"]=="12845461"
        assert record[0]["LinkSetDb"][0]["Link"][149]["Score"]=="21864314"
        assert record[0]["LinkSetDb"][0]["Link"][150]["Id"]=="16947073"
        assert record[0]["LinkSetDb"][0]["Link"][150]["Score"]=="21832153"
        assert record[0]["LinkSetDb"][0]["Link"][151]["Id"]=="7415301"
        assert record[0]["LinkSetDb"][0]["Link"][151]["Score"]=="21822396"
        assert record[0]["LinkSetDb"][0]["Link"][152]["Id"]=="16416239"
        assert record[0]["LinkSetDb"][0]["Link"][152]["Score"]=="21820165"
        assert record[0]["LinkSetDb"][0]["Link"][153]["Id"]=="4848922"
        assert record[0]["LinkSetDb"][0]["Link"][153]["Score"]=="21786194"
        assert record[0]["LinkSetDb"][0]["Link"][154]["Id"]=="12720164"
        assert record[0]["LinkSetDb"][0]["Link"][154]["Score"]=="21785319"
        assert record[0]["LinkSetDb"][0]["Link"][155]["Id"]=="17093987"
        assert record[0]["LinkSetDb"][0]["Link"][155]["Score"]=="21750370"
        assert record[0]["LinkSetDb"][0]["Link"][156]["Id"]=="16769006"
        assert record[0]["LinkSetDb"][0]["Link"][156]["Score"]=="21735873"
        assert record[0]["LinkSetDb"][0]["Link"][157]["Id"]=="17954835"
        assert record[0]["LinkSetDb"][0]["Link"][157]["Score"]=="21733933"
        assert record[0]["LinkSetDb"][0]["Link"][158]["Id"]=="15236134"
        assert record[0]["LinkSetDb"][0]["Link"][158]["Score"]=="21640099"
        assert record[0]["LinkSetDb"][0]["Link"][159]["Id"]=="12524603"
        assert record[0]["LinkSetDb"][0]["Link"][159]["Score"]=="21636724"
        assert record[0]["LinkSetDb"][0]["Link"][160]["Id"]=="16749985"
        assert record[0]["LinkSetDb"][0]["Link"][160]["Score"]=="21628926"
        assert record[0]["LinkSetDb"][0]["Link"][161]["Id"]=="3213296"
        assert record[0]["LinkSetDb"][0]["Link"][161]["Score"]=="21490232"
        assert record[0]["LinkSetDb"][0]["Link"][162]["Id"]=="11409026"
        assert record[0]["LinkSetDb"][0]["Link"][162]["Score"]=="21061296"
        assert record[0]["LinkSetDb"][0]["Link"][163]["Id"]=="9725288"
        assert record[0]["LinkSetDb"][0]["Link"][163]["Score"]=="21053585"
        assert record[0]["LinkSetDb"][0]["Link"][164]["Id"]=="6217136"
        assert record[0]["LinkSetDb"][0]["Link"][164]["Score"]=="21042914"
        assert record[0]["LinkSetDb"][0]["Link"][165]["Id"]=="663071"
        assert record[0]["LinkSetDb"][0]["Link"][165]["Score"]=="20926141"
        assert record[0]["LinkSetDb"][0]["Link"][166]["Id"]=="10341802"
        assert record[0]["LinkSetDb"][0]["Link"][166]["Score"]=="20797282"
        assert record[0]["LinkSetDb"][0]["Link"][167]["Id"]=="6473764"
        assert record[0]["LinkSetDb"][0]["Link"][167]["Score"]=="20757680"
        assert record[0]["LinkSetDb"][0]["Link"][168]["Id"]=="2584497"
        assert record[0]["LinkSetDb"][0]["Link"][168]["Score"]=="20521350"
        assert record[0]["LinkSetDb"][0]["Link"][169]["Id"]=="8338105"
        assert record[0]["LinkSetDb"][0]["Link"][169]["Score"]=="20501334"
        assert record[0]["LinkSetDb"][0]["Link"][170]["Id"]=="18053822"
        assert record[0]["LinkSetDb"][0]["Link"][170]["Score"]=="20275078"
        assert record[0]["LinkSetDb"][0]["Link"][171]["Id"]=="4058411"
        assert record[0]["LinkSetDb"][0]["Link"][171]["Score"]=="20161667"
        assert record[0]["LinkSetDb"][0]["Link"][172]["Id"]=="11669077"
        assert record[0]["LinkSetDb"][0]["Link"][172]["Score"]=="19993282"
        assert record[0]["LinkSetDb"][0]["Link"][173]["Id"]=="11781922"
        assert record[0]["LinkSetDb"][0]["Link"][173]["Score"]=="19969425"
        assert record[0]["LinkSetDb"][0]["Link"][174]["Id"]=="9793138"
        assert record[0]["LinkSetDb"][0]["Link"][174]["Score"]=="19952972"
        assert record[0]["LinkSetDb"][0]["Link"][175]["Id"]=="9391495"
        assert record[0]["LinkSetDb"][0]["Link"][175]["Score"]=="19815538"
        assert record[0]["LinkSetDb"][0]["Link"][176]["Id"]=="10803203"
        assert record[0]["LinkSetDb"][0]["Link"][176]["Score"]=="19495693"
        assert record[0]["LinkSetDb"][0]["Link"][177]["Id"]=="7326186"
        assert record[0]["LinkSetDb"][0]["Link"][177]["Score"]=="19273989"
        assert record[0]["LinkSetDb"][0]["Link"][178]["Id"]=="11868066"
        assert record[0]["LinkSetDb"][0]["Link"][178]["Score"]=="19220137"
        assert record[0]["LinkSetDb"][0]["Link"][179]["Id"]=="10904988"
        assert record[0]["LinkSetDb"][0]["Link"][179]["Score"]=="19203510"
        assert record[0]["LinkSetDb"][0]["Link"][180]["Id"]=="3288780"
        assert record[0]["LinkSetDb"][0]["Link"][180]["Score"]=="18958114"
        assert record[0]["LinkSetDb"][0]["Link"][181]["Id"]=="2047316"
        assert record[0]["LinkSetDb"][0]["Link"][181]["Score"]=="18907473"
        assert record[0]["LinkSetDb"][0]["Link"][182]["Id"]=="12237004"
        assert record[0]["LinkSetDb"][0]["Link"][182]["Score"]=="18751474"
        assert record[0]["LinkSetDb"][0]["Link"][183]["Id"]=="5627987"
        assert record[0]["LinkSetDb"][0]["Link"][183]["Score"]=="18741903"
        assert record[0]["LinkSetDb"][0]["Link"][184]["Id"]=="9269670"
        assert record[0]["LinkSetDb"][0]["Link"][184]["Score"]=="18666426"
        assert record[0]["LinkSetDb"][0]["Link"][185]["Id"]=="8903059"
        assert record[0]["LinkSetDb"][0]["Link"][185]["Score"]=="18653874"
        assert record[0]["LinkSetDb"][0]["Link"][186]["Id"]=="5594242"
        assert record[0]["LinkSetDb"][0]["Link"][186]["Score"]=="18548780"
        assert record[0]["LinkSetDb"][0]["Link"][187]["Id"]=="7068417"
        assert record[0]["LinkSetDb"][0]["Link"][187]["Score"]=="18390022"
        assert record[0]["LinkSetDb"][0]["Link"][188]["Id"]=="7330196"
        assert record[0]["LinkSetDb"][0]["Link"][188]["Score"]=="18371587"
        assert record[0]["LinkSetDb"][0]["Link"][189]["Id"]=="7408592"
        assert record[0]["LinkSetDb"][0]["Link"][189]["Score"]=="18275541"
        assert record[0]["LinkSetDb"][0]["Link"][190]["Id"]=="8835983"
        assert record[0]["LinkSetDb"][0]["Link"][190]["Score"]=="18176923"
        assert record[0]["LinkSetDb"][0]["Link"][191]["Id"]=="6940010"
        assert record[0]["LinkSetDb"][0]["Link"][191]["Score"]=="18011066"
        assert record[0]["LinkSetDb"][0]["Link"][192]["Id"]=="10499712"
        assert record[0]["LinkSetDb"][0]["Link"][192]["Score"]=="17943586"
        assert record[0]["LinkSetDb"][0]["Link"][193]["Id"]=="4539876"
        assert record[0]["LinkSetDb"][0]["Link"][193]["Score"]=="17915154"
        assert record[0]["LinkSetDb"][0]["Link"][194]["Id"]=="1943587"
        assert record[0]["LinkSetDb"][0]["Link"][194]["Score"]=="17752606"
        assert record[0]["LinkSetDb"][0]["Link"][195]["Id"]=="9847909"
        assert record[0]["LinkSetDb"][0]["Link"][195]["Score"]=="17568386"
        assert record[0]["LinkSetDb"][0]["Link"][196]["Id"]=="11578071"
        assert record[0]["LinkSetDb"][0]["Link"][196]["Score"]=="17561413"
        assert record[0]["LinkSetDb"][0]["Link"][197]["Id"]=="11789473"
        assert record[0]["LinkSetDb"][0]["Link"][197]["Score"]=="17435433"
        assert record[0]["LinkSetDb"][0]["Link"][198]["Id"]=="9885599"
        assert record[0]["LinkSetDb"][0]["Link"][198]["Score"]=="17383598"
        assert record[0]["LinkSetDb"][0]["Link"][199]["Id"]=="7423836"
        assert record[0]["LinkSetDb"][0]["Link"][199]["Score"]=="17196872"
        assert record[0]["LinkSetDb"][0]["Link"][200]["Id"]=="10688063"
        assert record[0]["LinkSetDb"][0]["Link"][200]["Score"]=="16453112"
        assert record[0]["LinkSetDb"][0]["Link"][201]["Id"]=="11695100"
        assert record[0]["LinkSetDb"][0]["Link"][201]["Score"]=="16352760"
        assert record[0]["LinkSetDb"][0]["Link"][202]["Id"]=="11329658"
        assert record[0]["LinkSetDb"][0]["Link"][202]["Score"]=="16089885"
        assert record[0]["LinkSetDb"][0]["Link"][203]["Id"]=="11939665"
        assert record[0]["LinkSetDb"][0]["Link"][203]["Score"]=="15947974"
        assert record[0]["LinkSetDb"][0]["Link"][204]["Id"]=="5512349"
        assert record[0]["LinkSetDb"][0]["Link"][204]["Score"]=="15647685"
        assert record[0]["LinkSetDb"][0]["Link"][205]["Id"]=="2222794"
        assert record[0]["LinkSetDb"][0]["Link"][205]["Score"]=="14981157"
        assert record[0]["LinkSetDb"][0]["Link"][206]["Id"]=="5998281"
        assert record[0]["LinkSetDb"][0]["Link"][206]["Score"]=="14226588"
        assert record[0]["LinkSetDb"][0]["Link"][207]["Id"]=="10475937"
        assert record[0]["LinkSetDb"][0]["Link"][207]["Score"]=="13934390"
        assert record[0]["LinkSetDb"][0]["Link"][208]["Id"]=="5046513"
        assert record[0]["LinkSetDb"][0]["Link"][208]["Score"]=="12769605"
        assert record[0]["LinkSetDb"][0]["Link"][209]["Id"]=="1539132"
        assert record[0]["LinkSetDb"][0]["Link"][209]["Score"]=="12395064"
        assert record[0]["LinkSetDb"][0]["Link"][210]["Id"]=="4414214"
        assert record[0]["LinkSetDb"][0]["Link"][210]["Score"]=="10113539"

    def t_pubmed3(self):
        '''Test parsing pubmed link returned by ELink (third test)
        '''
        # Create a hyperlink to the first link available for PMID 10611131
        # in PubMed
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="10611131", cmd="prlinks")

        input = open('Entrez/elink5.xml')
        record = Entrez.read(input)

        assert record[0]["DbFrom"]=="pubmed"
        assert len(record[0]["IdUrlList"])==1
        assert record[0]["IdUrlList"][0]
        assert record[0]["IdUrlList"][0]["Id"]=="10611131"
        assert len(record[0]["IdUrlList"][0]["ObjUrl"])==1
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Url"]=="http://brain.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=10611131"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["IconUrl"]=="http://www.ncbi.nlm.nih.gov/entrez/query/egifs/http:--highwire.stanford.edu-icons-externalservices-pubmed-custom-oxfordjournals_final_free.gif"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["SubjectType"]==["publishers/providers"]
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Attribute"][0]=="publisher of information in URL"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Attribute"][1]=="full-text online"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["Name"]=="HighWire Press"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["NameAbbr"]=="HighWire"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["Id"]=="3051"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["Url"]=="http://highwire.stanford.edu"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["IconUrl"]=="http://highwire.stanford.edu/icons/externalservices/pubmed/highwirepress.jpg"


    def t_pubmed4(self):
        '''Test parsing pubmed links returned by ELink (fourth test)
        '''
        # List all available links in PubMed, except for libraries, for
        # PMIDs 12085856 and 12085853
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12085856,12085853", cmd="llinks")
        input = open('Entrez/elink6.xml')
        record = Entrez.read(input)
        
        assert record[0]["DbFrom"]=="pubmed"
        assert len(record[0]["IdUrlList"])==2
        assert record[0]["IdUrlList"][0]["Id"]=="12085856"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Url"]=="http://symptomresearch.nih.gov/chapter_1/index.htm"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["SubjectType"]==["online tutorials/courses"]
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["Name"]=="New England Research Institutes Inc."
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["NameAbbr"]=="NERI"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["Id"]=="3291"
        assert record[0]["IdUrlList"][0]["ObjUrl"][0]["Provider"]["Url"]=="http://www.symptomresearch.com"
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["Url"]=="http://www.nlm.nih.gov/medlineplus/coronaryarterybypasssurgery.html"
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["IconUrl"]=="http://www.ncbi.nlm.nih.gov/entrez/query/egifs/http:--www.nlm.nih.gov-medlineplus-images-linkout_sm.gif"
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["LinkName"]=="Coronary Artery Bypass Surgery"
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["SubjectType"]==["consumer health"]
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["Provider"]["Name"]=="MedlinePlus Health Information"
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["Provider"]["NameAbbr"]=="MEDPLUS"
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["Provider"]["Id"]=="3162"
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["Provider"]["Url"]=="http://medlineplus.gov/"
        assert record[0]["IdUrlList"][0]["ObjUrl"][1]["Provider"]["IconUrl"]=="http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif"

        assert record[0]["IdUrlList"][1]["Id"]=="12085853"
        assert len(record[0]["IdUrlList"][1]["ObjUrl"])==4
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["Url"]=="http://www.nlm.nih.gov/medlineplus/arrhythmia.html"
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["IconUrl"]=="http://www.ncbi.nlm.nih.gov/entrez/query/egifs/http:--www.nlm.nih.gov-medlineplus-images-linkout_sm.gif"
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["LinkName"]=="Arrhythmia"
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["SubjectType"]==["consumer health"]
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["Provider"]["Name"]=="MedlinePlus Health Information"
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["Provider"]["NameAbbr"]=="MEDPLUS"
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["Provider"]["Id"]=="3162"
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["Provider"]["Url"]=="http://medlineplus.gov/"
        assert record[0]["IdUrlList"][1]["ObjUrl"][0]["Provider"]["IconUrl"]=="http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif"
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["Url"]=="http://www.nlm.nih.gov/medlineplus/exerciseandphysicalfitness.html"
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["IconUrl"]=="http://www.ncbi.nlm.nih.gov/entrez/query/egifs/http:--www.nlm.nih.gov-medlineplus-images-linkout_sm.gif"
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["LinkName"]=="Exercise and Physical Fitness"
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["SubjectType"]==["consumer health"]
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["Provider"]["Name"]=="MedlinePlus Health Information"
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["Provider"]["NameAbbr"]=="MEDPLUS"
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["Provider"]["Id"]=="3162"
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["Provider"]["Url"]=="http://medlineplus.gov/"
        assert record[0]["IdUrlList"][1]["ObjUrl"][1]["Provider"]["IconUrl"]=="http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif"
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["Url"]=="http://www.nlm.nih.gov/medlineplus/exerciseforchildren.html"
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["IconUrl"]=="http://www.ncbi.nlm.nih.gov/entrez/query/egifs/http:--www.nlm.nih.gov-medlineplus-images-linkout_sm.gif"
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["LinkName"]=="Exercise for Children"
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["SubjectType"]==["consumer health"]
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["Provider"]["Name"]=="MedlinePlus Health Information"
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["Provider"]["NameAbbr"]=="MEDPLUS"
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["Provider"]["Id"]=="3162"
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["Provider"]["Url"]=="http://medlineplus.gov/"
        assert record[0]["IdUrlList"][1]["ObjUrl"][2]["Provider"]["IconUrl"]=="http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif"
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["Url"]=="http://www.nlm.nih.gov/medlineplus/pacemakersandimplantabledefibrillators.html"
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["IconUrl"]=="http://www.ncbi.nlm.nih.gov/entrez/query/egifs/http:--www.nlm.nih.gov-medlineplus-images-linkout_sm.gif"
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["LinkName"]=="Pacemakers and Implantable Defibrillators"
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["SubjectType"]==["consumer health"]
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["Provider"]["Name"]=="MedlinePlus Health Information"
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["Provider"]["NameAbbr"]=="MEDPLUS"
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["Provider"]["Id"]=="3162"
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["Provider"]["Url"]=="http://medlineplus.gov/"
        assert record[0]["IdUrlList"][1]["ObjUrl"][3]["Provider"]["IconUrl"]=="http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif"

    def t_pubmed5(self):
        '''Test parsing pubmed links returned by ELink (fifth test)
        '''
        # List Entrez database links for PubMed PMIDs 12169658 and 11748140
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12169658,11748140",
        #                      cmd="acheck")
        input = open('Entrez/elink7.xml')
        record = Entrez.read(input)
        
        assert len(record)==1
        assert record[0]["DbFrom"]=="pubmed"
        assert len(record[0]["IdCheckList"])==2
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["Id"]=="12169658"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["DbTo"]=="books"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["LinkName"]=="pubmed_books_refs"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["MenuTag"]=="Cited in Books"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["HtmlTag"]=="Cited in Books"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][0]["Priority"]=="185"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["DbTo"]=="gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["LinkName"]=="pubmed_gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["MenuTag"]=="Gene Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["HtmlTag"]=="Gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][1]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["DbTo"]=="geo"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["LinkName"]=="pubmed_geo"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["MenuTag"]=="GEO Profile Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["HtmlTag"]=="GEO Profiles"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][2]["Priority"]=="170"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["DbTo"]=="homologene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["LinkName"]=="pubmed_homologene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["MenuTag"]=="HomoloGene Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["HtmlTag"]=="HomoloGene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][3]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["DbTo"]=="nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["LinkName"]=="pubmed_nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["MenuTag"]=="CoreNucleotide Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["HtmlTag"]=="CoreNucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][4]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][5]["DbTo"]=="nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][5]["LinkName"]=="pubmed_nuccore_refseq"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][5]["MenuTag"]=="CoreNucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][5]["HtmlTag"]=="CoreNucleotide (RefSeq)"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][5]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][6]["DbTo"]=="nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][6]["LinkName"]=="pubmed_nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][6]["MenuTag"]=="Nucleotide Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][6]["HtmlTag"]=="Nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][6]["Priority"]=="135"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][7]["DbTo"]=="nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][7]["LinkName"]=="pubmed_nucleotide_refseq"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][7]["MenuTag"]=="Nucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][7]["HtmlTag"]=="Nucleotide (RefSeq)"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][7]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][8]["DbTo"]=="pcsubstance"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][8]["LinkName"]=="pubmed_pcsubstance_mesh"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][8]["MenuTag"]=="Substance (MeSH Keyword)"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][8]["HtmlTag"]=="Substance (MeSH Keyword)"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][8]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][9]["DbTo"]=="pmc"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][9]["LinkName"]=="pubmed_pmc_refs"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][9]["MenuTag"]=="Cited in PMC"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][9]["HtmlTag"]=="Cited in PMC"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][9]["Url"]=="http://www.pubmedcentral.gov/tocrender.fcgi?action=cited&tool=pubmed&pubmedid=<@UID@>"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][9]["Priority"]=="180"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][10]["DbTo"]=="protein"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][10]["LinkName"]=="pubmed_protein"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][10]["MenuTag"]=="Protein Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][10]["HtmlTag"]=="Protein"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][10]["Priority"]=="140"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][11]["DbTo"]=="protein"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][11]["LinkName"]=="pubmed_protein_refseq"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][11]["MenuTag"]=="Protein (RefSeq) Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][11]["HtmlTag"]=="Protein (RefSeq)"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][11]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][12]["DbTo"]=="pubmed"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][12]["LinkName"]=="pubmed_pubmed"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][12]["MenuTag"]=="Related Articles"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][12]["HtmlTag"]=="Related Articles"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][12]["Priority"]=="1"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][13]["DbTo"]=="taxonomy"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][13]["LinkName"]=="pubmed_taxonomy_entrez"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][13]["MenuTag"]=="Taxonomy via GenBank"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][13]["HtmlTag"]=="Taxonomy via GenBank"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][13]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][14]["DbTo"]=="unigene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][14]["LinkName"]=="pubmed_unigene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][14]["MenuTag"]=="UniGene Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][14]["HtmlTag"]=="UniGene"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][14]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][15]["DbTo"]=="LinkOut"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][15]["LinkName"]=="ExternalLink"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][15]["MenuTag"]=="LinkOut"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][15]["HtmlTag"]=="LinkOut"
        assert record[0]["IdCheckList"]["IdLinkSet"][0]["LinkInfo"][15]["Priority"]=="255"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["Id"]=="11748140"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["DbTo"]=="books"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["LinkName"]=="pubmed_books_refs"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["MenuTag"]=="Cited in Books"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["HtmlTag"]=="Cited in Books"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][0]["Priority"]=="185"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["DbTo"]=="gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["LinkName"]=="pubmed_gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["MenuTag"]=="Gene Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["HtmlTag"]=="Gene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][1]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["DbTo"]=="geo"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["LinkName"]=="pubmed_geo"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["MenuTag"]=="GEO Profile Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["HtmlTag"]=="GEO Profiles"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][2]["Priority"]=="170"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["DbTo"]=="nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["LinkName"]=="pubmed_nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["MenuTag"]=="CoreNucleotide Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["HtmlTag"]=="CoreNucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][3]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["DbTo"]=="nuccore"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["LinkName"]=="pubmed_nuccore_refseq"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["MenuTag"]=="CoreNucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["HtmlTag"]=="CoreNucleotide (RefSeq)"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][4]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][5]["DbTo"]=="nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][5]["LinkName"]=="pubmed_nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][5]["MenuTag"]=="Nucleotide Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][5]["HtmlTag"]=="Nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][5]["Priority"]=="135"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][6]["DbTo"]=="nucleotide"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][6]["LinkName"]=="pubmed_nucleotide_refseq"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][6]["MenuTag"]=="Nucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][6]["HtmlTag"]=="Nucleotide (RefSeq)"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][6]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][7]["DbTo"]=="pmc"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][7]["LinkName"]=="pubmed_pmc_refs"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][7]["MenuTag"]=="Cited in PMC"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][7]["HtmlTag"]=="Cited in PMC"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][7]["Url"]=="http://www.pubmedcentral.gov/tocrender.fcgi?action=cited&tool=pubmed&pubmedid=<@UID@>"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][7]["Priority"]=="180"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][8]["DbTo"]=="protein"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][8]["LinkName"]=="pubmed_protein"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][8]["MenuTag"]=="Protein Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][8]["HtmlTag"]=="Protein"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][8]["Priority"]=="140"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][9]["DbTo"]=="protein"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][9]["LinkName"]=="pubmed_protein_refseq"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][9]["MenuTag"]=="Protein (RefSeq) Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][9]["HtmlTag"]=="Protein (RefSeq)"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][9]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][10]["DbTo"]=="pubmed"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][10]["LinkName"]=="pubmed_pubmed"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][10]["MenuTag"]=="Related Articles"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][10]["HtmlTag"]=="Related Articles"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][10]["Priority"]=="1"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][11]["DbTo"]=="taxonomy"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][11]["LinkName"]=="pubmed_taxonomy_entrez"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][11]["MenuTag"]=="Taxonomy via GenBank"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][11]["HtmlTag"]=="Taxonomy via GenBank"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][11]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][12]["DbTo"]=="unigene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][12]["LinkName"]=="pubmed_unigene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][12]["MenuTag"]=="UniGene Links"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][12]["HtmlTag"]=="UniGene"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][12]["Priority"]=="128"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][13]["DbTo"]=="LinkOut"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][13]["LinkName"]=="ExternalLink"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][13]["MenuTag"]=="LinkOut"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][13]["HtmlTag"]=="LinkOut"
        assert record[0]["IdCheckList"]["IdLinkSet"][1]["LinkInfo"][13]["Priority"]=="255"


    def t_pubmed6(self):
        '''Test parsing pubmed links returned by ELink (sixth test)
        '''
        # Check for the existence of a Related Articles link for PMIDs
        # 0611131, 111645 and 12068369
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="10611131,111645",
        #                      id="12068369", cmd="ncheck")
        input = open('Entrez/elink8.xml')
        record = Entrez.read(input)

        assert len(record)==1
        assert record[0]["DbFrom"]=="pubmed"
        assert len(record[0]["IdCheckList"])==2
        assert len(record[0]["IdCheckList"]["Id"])==1
        assert record[0]["IdCheckList"]["Id"][0]=="12068369"
        assert len(record[0]["IdCheckList"]["Id"][0].attributes)==1
        assert record[0]["IdCheckList"]["Id"][0].attributes["HasNeighbor"]=="Y"

    def t_cancerchromosomes(self):
        '''Test parsing cancerchromosomes links returned by ELink
        '''
        # Retrieve neighbors for Cancer Chromosomes ID 2662 to the link
        # subset related by cytogenetics
        # To create the XML file, use
        # >>> Bio.Entrez.elink(dbfrom="cancerchromosomes",
        #                       db="cancerchromosomes", id="2662",
        #                       cmd="neighbor",
        #                       linkname="cancerchromosomes_cancerchromosomes_cyto")
        input = open('Entrez/elink9.xml')
        record = Entrez.read(input)
        
        assert record[0]["DbFrom"]=="cancerchromosomes"
	assert record[0]["IdList"]==["2662"]


class EGQueryTest(unittest.TestCase):
    '''Tests for parsing XML output returned by EGQuery
    '''
    def t_egquery1(self):
        '''Test parsing XML output returned by EGQuery (first test)
        '''
        # Display counts in XML for stem cells in each Entrez database
        # To create the XML file, use
        # >>> Bio.Entrez.egquery(term="stem cells")
        input = open('Entrez/egquery1.xml')
        record = Entrez.read(input)

        assert record["Term"]=="stem cells"

        assert record["eGQueryResult"][0]["DbName"]=="pubmed"
        assert record["eGQueryResult"][0]["MenuName"]=="PubMed"
        assert record["eGQueryResult"][0]["Count"]=="392"
        assert record["eGQueryResult"][0]["Status"]=="Ok"
        assert record["eGQueryResult"][1]["DbName"]=="pmc"
        assert record["eGQueryResult"][1]["MenuName"]=="PMC"
        assert record["eGQueryResult"][1]["Count"]=="173"
        assert record["eGQueryResult"][1]["Status"]=="Ok"
        assert record["eGQueryResult"][2]["DbName"]=="journals"
        assert record["eGQueryResult"][2]["MenuName"]=="Journals"
        assert record["eGQueryResult"][2]["Count"]=="0"
        assert record["eGQueryResult"][2]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][3]["DbName"]=="mesh"
        assert record["eGQueryResult"][3]["MenuName"]=="MeSH"
        assert record["eGQueryResult"][3]["Count"]=="0"
        assert record["eGQueryResult"][3]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][4]["DbName"]=="books"
        assert record["eGQueryResult"][4]["MenuName"]=="Books"
        assert record["eGQueryResult"][4]["Count"]=="10"
        assert record["eGQueryResult"][4]["Status"]=="Ok"
        assert record["eGQueryResult"][5]["DbName"]=="omim"
        assert record["eGQueryResult"][5]["MenuName"]=="OMIM"
        assert record["eGQueryResult"][5]["Count"]=="0"
        assert record["eGQueryResult"][5]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][6]["DbName"]=="omia"
        assert record["eGQueryResult"][6]["MenuName"]=="OMIA"
        assert record["eGQueryResult"][6]["Count"]=="0"
        assert record["eGQueryResult"][6]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][7]["DbName"]=="ncbisearch"
        assert record["eGQueryResult"][7]["MenuName"]=="NCBI Web Site"
        assert record["eGQueryResult"][7]["Count"]=="0"
        assert record["eGQueryResult"][7]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][8]["DbName"]=="nuccore"
        assert record["eGQueryResult"][8]["MenuName"]=="CoreNucleotide"
        assert record["eGQueryResult"][8]["Count"]=="0"
        assert record["eGQueryResult"][8]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][9]["DbName"]=="nucgss"
        assert record["eGQueryResult"][9]["MenuName"]=="GSS"
        assert record["eGQueryResult"][9]["Count"]=="0"
        assert record["eGQueryResult"][9]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][10]["DbName"]=="nucest"
        assert record["eGQueryResult"][10]["MenuName"]=="EST"
        assert record["eGQueryResult"][10]["Count"]=="0"
        assert record["eGQueryResult"][10]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][11]["DbName"]=="protein"
        assert record["eGQueryResult"][11]["MenuName"]=="Protein"
        assert record["eGQueryResult"][11]["Count"]=="0"
        assert record["eGQueryResult"][11]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][12]["DbName"]=="genome"
        assert record["eGQueryResult"][12]["MenuName"]=="Genome"
        assert record["eGQueryResult"][12]["Count"]=="0"
        assert record["eGQueryResult"][12]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][13]["DbName"]=="structure"
        assert record["eGQueryResult"][13]["MenuName"]=="Structure"
        assert record["eGQueryResult"][13]["Count"]=="0"
        assert record["eGQueryResult"][13]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][14]["DbName"]=="taxonomy"
        assert record["eGQueryResult"][14]["MenuName"]=="Taxonomy"
        assert record["eGQueryResult"][14]["Count"]=="0"
        assert record["eGQueryResult"][14]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][15]["DbName"]=="snp"
        assert record["eGQueryResult"][15]["MenuName"]=="SNP"
        assert record["eGQueryResult"][15]["Count"]=="0"
        assert record["eGQueryResult"][15]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][16]["DbName"]=="gene"
        assert record["eGQueryResult"][16]["MenuName"]=="Gene"
        assert record["eGQueryResult"][16]["Count"]=="0"
        assert record["eGQueryResult"][16]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][17]["DbName"]=="unigene"
        assert record["eGQueryResult"][17]["MenuName"]=="UniGene"
        assert record["eGQueryResult"][17]["Count"]=="0"
        assert record["eGQueryResult"][17]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][18]["DbName"]=="cdd"
        assert record["eGQueryResult"][18]["MenuName"]=="Conserved Domains"
        assert record["eGQueryResult"][18]["Count"]=="0"
        assert record["eGQueryResult"][18]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][19]["DbName"]=="domains"
        assert record["eGQueryResult"][19]["MenuName"]=="3D Domains"
        assert record["eGQueryResult"][19]["Count"]=="0"
        assert record["eGQueryResult"][19]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][20]["DbName"]=="unists"
        assert record["eGQueryResult"][20]["MenuName"]=="UniSTS"
        assert record["eGQueryResult"][20]["Count"]=="0"
        assert record["eGQueryResult"][20]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][21]["DbName"]=="popset"
        assert record["eGQueryResult"][21]["MenuName"]=="PopSet"
        assert record["eGQueryResult"][21]["Count"]=="0"
        assert record["eGQueryResult"][21]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][22]["DbName"]=="geo"
        assert record["eGQueryResult"][22]["MenuName"]=="GEO Profiles"
        assert record["eGQueryResult"][22]["Count"]=="0"
        assert record["eGQueryResult"][22]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][23]["DbName"]=="gds"
        assert record["eGQueryResult"][23]["MenuName"]=="GEO DataSets"
        assert record["eGQueryResult"][23]["Count"]=="0"
        assert record["eGQueryResult"][23]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][24]["DbName"]=="homologene"
        assert record["eGQueryResult"][24]["MenuName"]=="HomoloGene"
        assert record["eGQueryResult"][24]["Count"]=="0"
        assert record["eGQueryResult"][24]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][25]["DbName"]=="cancerchromosomes"
        assert record["eGQueryResult"][25]["MenuName"]=="CancerChromosomes"
        assert record["eGQueryResult"][25]["Count"]=="0"
        assert record["eGQueryResult"][25]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][26]["DbName"]=="pccompound"
        assert record["eGQueryResult"][26]["MenuName"]=="PubChem Compound"
        assert record["eGQueryResult"][26]["Count"]=="0"
        assert record["eGQueryResult"][26]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][27]["DbName"]=="pcsubstance"
        assert record["eGQueryResult"][27]["MenuName"]=="PubChem Substance"
        assert record["eGQueryResult"][27]["Count"]=="0"
        assert record["eGQueryResult"][27]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][28]["DbName"]=="pcassay"
        assert record["eGQueryResult"][28]["MenuName"]=="PubChem BioAssay"
        assert record["eGQueryResult"][28]["Count"]=="0"
        assert record["eGQueryResult"][28]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][29]["DbName"]=="nlmcatalog"
        assert record["eGQueryResult"][29]["MenuName"]=="NLM Catalog"
        assert record["eGQueryResult"][29]["Count"]=="2"
        assert record["eGQueryResult"][29]["Status"]=="Ok"
        assert record["eGQueryResult"][30]["DbName"]=="gensat"
        assert record["eGQueryResult"][30]["MenuName"]=="GENSAT"
        assert record["eGQueryResult"][30]["Count"]=="0"
        assert record["eGQueryResult"][30]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][31]["DbName"]=="probe"
        assert record["eGQueryResult"][31]["MenuName"]=="Probe"
        assert record["eGQueryResult"][31]["Count"]=="0"
        assert record["eGQueryResult"][31]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][32]["DbName"]=="genomeprj"
        assert record["eGQueryResult"][32]["MenuName"]=="Genome Project"
        assert record["eGQueryResult"][32]["Count"]=="0"
        assert record["eGQueryResult"][32]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][33]["DbName"]=="gap"
        assert record["eGQueryResult"][33]["MenuName"]=="dbGaP"
        assert record["eGQueryResult"][33]["Count"]=="0"
        assert record["eGQueryResult"][33]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][34]["DbName"]=="proteinclusters"
        assert record["eGQueryResult"][34]["MenuName"]=="Protein Clusters"
        assert record["eGQueryResult"][34]["Count"]=="0"
        assert record["eGQueryResult"][34]["Status"]=="Term or Database is not found"

    def t_egquery2(self):
        '''Test parsing XML output returned by EGQuery (second test)
        '''
        # Display counts in XML for brca1 or brca2 for each Entrez database
        # To create the XML file, use
        # >>> Bio.Entrez.egquery(term="brca1 OR brca2")
        input = open('Entrez/egquery2.xml')
        record = Entrez.read(input)

        assert record["Term"]=="brca1 OR brca2"

        assert record["eGQueryResult"][0]["DbName"]=="pubmed"
        assert record["eGQueryResult"][0]["MenuName"]=="PubMed"
        assert record["eGQueryResult"][0]["Count"]=="0"
        assert record["eGQueryResult"][0]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][1]["DbName"]=="pmc"
        assert record["eGQueryResult"][1]["MenuName"]=="PMC"
        assert record["eGQueryResult"][1]["Count"]=="2739"
        assert record["eGQueryResult"][1]["Status"]=="Ok"
        assert record["eGQueryResult"][2]["DbName"]=="journals"
        assert record["eGQueryResult"][2]["MenuName"]=="Journals"
        assert record["eGQueryResult"][2]["Count"]=="0"
        assert record["eGQueryResult"][2]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][3]["DbName"]=="mesh"
        assert record["eGQueryResult"][3]["MenuName"]=="MeSH"
        assert record["eGQueryResult"][3]["Count"]=="29"
        assert record["eGQueryResult"][3]["Status"]=="Ok"
        assert record["eGQueryResult"][4]["DbName"]=="books"
        assert record["eGQueryResult"][4]["MenuName"]=="Books"
        assert record["eGQueryResult"][4]["Count"]=="392"
        assert record["eGQueryResult"][4]["Status"]=="Ok"
        assert record["eGQueryResult"][5]["DbName"]=="omim"
        assert record["eGQueryResult"][5]["MenuName"]=="OMIM"
        assert record["eGQueryResult"][5]["Count"]=="149"
        assert record["eGQueryResult"][5]["Status"]=="Ok"
        assert record["eGQueryResult"][6]["DbName"]=="omia"
        assert record["eGQueryResult"][6]["MenuName"]=="OMIA"
        assert record["eGQueryResult"][6]["Count"]=="0"
        assert record["eGQueryResult"][6]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][7]["DbName"]=="ncbisearch"
        assert record["eGQueryResult"][7]["MenuName"]=="NCBI Web Site"
        assert record["eGQueryResult"][7]["Count"]=="13"
        assert record["eGQueryResult"][7]["Status"]=="Ok"
        assert record["eGQueryResult"][8]["DbName"]=="nuccore"
        assert record["eGQueryResult"][8]["MenuName"]=="CoreNucleotide"
        assert record["eGQueryResult"][8]["Count"]=="4917"
        assert record["eGQueryResult"][8]["Status"]=="Ok"
        assert record["eGQueryResult"][9]["DbName"]=="nucgss"
        assert record["eGQueryResult"][9]["MenuName"]=="GSS"
        assert record["eGQueryResult"][9]["Count"]=="184"
        assert record["eGQueryResult"][9]["Status"]=="Ok"
        assert record["eGQueryResult"][10]["DbName"]=="nucest"
        assert record["eGQueryResult"][10]["MenuName"]=="EST"
        assert record["eGQueryResult"][10]["Count"]=="600"
        assert record["eGQueryResult"][10]["Status"]=="Ok"
        assert record["eGQueryResult"][11]["DbName"]=="protein"
        assert record["eGQueryResult"][11]["MenuName"]=="Protein"
        assert record["eGQueryResult"][11]["Count"]=="6779"
        assert record["eGQueryResult"][11]["Status"]=="Ok"
        assert record["eGQueryResult"][12]["DbName"]=="genome"
        assert record["eGQueryResult"][12]["MenuName"]=="Genome"
        assert record["eGQueryResult"][12]["Count"]=="44"
        assert record["eGQueryResult"][12]["Status"]=="Ok"
        assert record["eGQueryResult"][13]["DbName"]=="structure"
        assert record["eGQueryResult"][13]["MenuName"]=="Structure"
        assert record["eGQueryResult"][13]["Count"]=="29"
        assert record["eGQueryResult"][13]["Status"]=="Ok"
        assert record["eGQueryResult"][14]["DbName"]=="taxonomy"
        assert record["eGQueryResult"][14]["MenuName"]=="Taxonomy"
        assert record["eGQueryResult"][14]["Count"]=="0"
        assert record["eGQueryResult"][14]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][15]["DbName"]=="snp"
        assert record["eGQueryResult"][15]["MenuName"]=="SNP"
        assert record["eGQueryResult"][15]["Count"]=="2013"
        assert record["eGQueryResult"][15]["Status"]=="Ok"
        assert record["eGQueryResult"][16]["DbName"]=="gene"
        assert record["eGQueryResult"][16]["MenuName"]=="Gene"
        assert record["eGQueryResult"][16]["Count"]=="1775"
        assert record["eGQueryResult"][16]["Status"]=="Ok"
        assert record["eGQueryResult"][17]["DbName"]=="unigene"
        assert record["eGQueryResult"][17]["MenuName"]=="UniGene"
        assert record["eGQueryResult"][17]["Count"]=="207"
        assert record["eGQueryResult"][17]["Status"]=="Ok"
        assert record["eGQueryResult"][18]["DbName"]=="cdd"
        assert record["eGQueryResult"][18]["MenuName"]=="Conserved Domains"
        assert record["eGQueryResult"][18]["Count"]=="17"
        assert record["eGQueryResult"][18]["Status"]=="Ok"
        assert record["eGQueryResult"][19]["DbName"]=="domains"
        assert record["eGQueryResult"][19]["MenuName"]=="3D Domains"
        assert record["eGQueryResult"][19]["Count"]=="131"
        assert record["eGQueryResult"][19]["Status"]=="Ok"
        assert record["eGQueryResult"][20]["DbName"]=="unists"
        assert record["eGQueryResult"][20]["MenuName"]=="UniSTS"
        assert record["eGQueryResult"][20]["Count"]=="198"
        assert record["eGQueryResult"][20]["Status"]=="Ok"
        assert record["eGQueryResult"][21]["DbName"]=="popset"
        assert record["eGQueryResult"][21]["MenuName"]=="PopSet"
        assert record["eGQueryResult"][21]["Count"]=="43"
        assert record["eGQueryResult"][21]["Status"]=="Ok"
        assert record["eGQueryResult"][22]["DbName"]=="geo"
        assert record["eGQueryResult"][22]["MenuName"]=="GEO Profiles"
        assert record["eGQueryResult"][22]["Count"]=="128692"
        assert record["eGQueryResult"][22]["Status"]=="Ok"
        assert record["eGQueryResult"][23]["DbName"]=="gds"
        assert record["eGQueryResult"][23]["MenuName"]=="GEO DataSets"
        assert record["eGQueryResult"][23]["Count"]=="21"
        assert record["eGQueryResult"][23]["Status"]=="Ok"
        assert record["eGQueryResult"][24]["DbName"]=="homologene"
        assert record["eGQueryResult"][24]["MenuName"]=="HomoloGene"
        assert record["eGQueryResult"][24]["Count"]=="50"
        assert record["eGQueryResult"][24]["Status"]=="Ok"
        assert record["eGQueryResult"][25]["DbName"]=="cancerchromosomes"
        assert record["eGQueryResult"][25]["MenuName"]=="CancerChromosomes"
        assert record["eGQueryResult"][25]["Count"]=="18"
        assert record["eGQueryResult"][25]["Status"]=="Ok"
        assert record["eGQueryResult"][26]["DbName"]=="pccompound"
        assert record["eGQueryResult"][26]["MenuName"]=="PubChem Compound"
        assert record["eGQueryResult"][26]["Count"]=="0"
        assert record["eGQueryResult"][26]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][27]["DbName"]=="pcsubstance"
        assert record["eGQueryResult"][27]["MenuName"]=="PubChem Substance"
        assert record["eGQueryResult"][27]["Count"]=="26"
        assert record["eGQueryResult"][27]["Status"]=="Ok"
        assert record["eGQueryResult"][28]["DbName"]=="pcassay"
        assert record["eGQueryResult"][28]["MenuName"]=="PubChem BioAssay"
        assert record["eGQueryResult"][28]["Count"]=="0"
        assert record["eGQueryResult"][28]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][29]["DbName"]=="nlmcatalog"
        assert record["eGQueryResult"][29]["MenuName"]=="NLM Catalog"
        assert record["eGQueryResult"][29]["Count"]=="31"
        assert record["eGQueryResult"][29]["Status"]=="Ok"
        assert record["eGQueryResult"][30]["DbName"]=="gensat"
        assert record["eGQueryResult"][30]["MenuName"]=="GENSAT"
        assert record["eGQueryResult"][30]["Count"]=="0"
        assert record["eGQueryResult"][30]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][31]["DbName"]=="probe"
        assert record["eGQueryResult"][31]["MenuName"]=="Probe"
        assert record["eGQueryResult"][31]["Count"]=="1410"
        assert record["eGQueryResult"][31]["Status"]=="Ok"
        assert record["eGQueryResult"][32]["DbName"]=="genomeprj"
        assert record["eGQueryResult"][32]["MenuName"]=="Genome Project"
        assert record["eGQueryResult"][32]["Count"]=="0"
        assert record["eGQueryResult"][32]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][33]["DbName"]=="gap"
        assert record["eGQueryResult"][33]["MenuName"]=="dbGaP"
        assert record["eGQueryResult"][33]["Count"]=="0"
        assert record["eGQueryResult"][33]["Status"]=="Term or Database is not found"
        assert record["eGQueryResult"][34]["DbName"]=="proteinclusters"
        assert record["eGQueryResult"][34]["MenuName"]=="Protein Clusters"
        assert record["eGQueryResult"][34]["Count"]=="0"
        assert record["eGQueryResult"][34]["Status"]=="Term or Database is not found"

class ESpellTest(unittest.TestCase):
    '''Tests for parsing XML output returned by ESpell
    '''
    def t_espell(self):
        '''Test parsing XML output returned by ESpell
        '''
        # Request suggestions for the PubMed search biopythooon
        # To create the XML file, use
        # >>> Bio.Entrez.espell(db="pubmed", term="biopythooon")
        input = open('Entrez/espell.xml')
        record = Entrez.read(input)
        assert record["Database"]=="pubmed"
        assert record["Query"]=="biopythooon"
        assert record["CorrectedQuery"]=="biopython"
        assert len(record["SpelledQuery"])==1
        assert record["SpelledQuery"][0]=="biopython"
        assert record["SpelledQuery"][0].tag=="Replaced"


class EFetchTest(unittest.TestCase):
    '''Tests for parsing XML output returned by EFetch
    '''
    def t_pubmed1(self):
        '''Test parsing XML returned by EFetch, PubMed database (first test)
        '''
        # In PubMed display PMIDs 12091962 and 9997 in xml retrieval mode
        # and abstract retrieval type.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', id='12091962,9997',
        #                       retmode='xml', rettype='abstract')
        input = open('Entrez/pubmed1.xml')
        record = Entrez.read(input)

        assert record[0]["MedlineCitation"].attributes["Owner"]=="KIE"
        assert record[0]["MedlineCitation"].attributes["Status"]=="MEDLINE"
        assert record[0]["MedlineCitation"]["PMID"]=="12091962"
        assert record[0]["MedlineCitation"]["DateCreated"]["Year"]=="1991"
        assert record[0]["MedlineCitation"]["DateCreated"]["Month"]=="01"
        assert record[0]["MedlineCitation"]["DateCreated"]["Day"]=="22"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Year"]=="1991"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Month"]=="01"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Day"]=="22"
        assert record[0]["MedlineCitation"]["DateRevised"]["Year"]=="2007"
        assert record[0]["MedlineCitation"]["DateRevised"]["Month"]=="11"
        assert record[0]["MedlineCitation"]["DateRevised"]["Day"]=="15"
        assert record[0]["MedlineCitation"]["Article"].attributes["PubModel"]=="Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"]=="1043-1578"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes["IssnType"]=="Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].attributes["CitedMedium"]=="Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"]=="17"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"]=="1"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]=="1990"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Season"]=="Spring"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["Title"]=="Social justice (San Francisco, Calif.)"
        assert record[0]["MedlineCitation"]["Article"]["ArticleTitle"]=="The treatment of AIDS behind the walls of correctional facilities."
        assert record[0]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"]=="113-25"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"]=='Y'
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"]=="Olivero"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"]=="J Michael"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"]=="JM"
        assert record[0]["MedlineCitation"]["Article"]["Language"]==["eng"]
        assert record[0]["MedlineCitation"]["Article"]["PublicationTypeList"]==["Journal Article", "Review"]
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["Country"]=="United States"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"]=="Soc Justice"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"]=="9891830"
        assert record[0]["MedlineCitation"]["CitationSubset"]==["E"]
        assert record[0]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"]=="AIDS Serodiagnosis"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"]=="Acquired Immunodeficiency Syndrome"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"]=="Civil Rights"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"]=="HIV Seropositivity"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"]=="Humans"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"]=="Jurisprudence"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"]=="Law Enforcement"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"]=="Mass Screening"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"]=="Minority Groups"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"]=="Organizational Policy"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"]=="Patient Care"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"]=="Prejudice"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"]=="Prisoners"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][13]["DescriptorName"]=="Public Policy"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][13]["DescriptorName"].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][14]["DescriptorName"]=="Quarantine"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][14]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][15]["DescriptorName"]=="Social Control, Formal"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][15]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][16]["DescriptorName"]=="Statistics as Topic"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][16]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][17]["DescriptorName"]=="Stereotyping"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][17]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][18]["DescriptorName"]=="United States"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][18]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["NumberOfReferences"]=="63"
        assert record[0]["MedlineCitation"]["OtherID"][0]=="31840"
        assert record[0]["MedlineCitation"]["OtherID"][0].attributes["Source"]=="KIE"
        assert record[0]["MedlineCitation"]["KeywordList"][0].attributes["Owner"]=="KIE"
        assert record[0]["MedlineCitation"]["KeywordList"][0][0]=="Health Care and Public Health"
        assert record[0]["MedlineCitation"]["KeywordList"][0][0].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["KeywordList"][0][1]=="Legal Approach"
        assert record[0]["MedlineCitation"]["KeywordList"][0][1].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["GeneralNote"][0]=="14 fn."
        assert record[0]["MedlineCitation"]["GeneralNote"][0].attributes["Owner"]=="KIE"
        assert record[0]["MedlineCitation"]["GeneralNote"][1]=="KIE BoB Subject Heading: AIDS"
        assert record[0]["MedlineCitation"]["GeneralNote"][1].attributes["Owner"]=="KIE"
        assert record[0]["MedlineCitation"]["GeneralNote"][2]=="63 refs."
        assert record[0]["MedlineCitation"]["GeneralNote"][2].attributes["Owner"]=="KIE"
        assert record[0]["PubmedData"]["History"][0][0].attributes["PubStatus"]=="pubmed"
        assert record[0]["PubmedData"]["History"][0][0]["Year"]=="1990"
        assert record[0]["PubmedData"]["History"][0][0]["Month"]=="4"
        assert record[0]["PubmedData"]["History"][0][0]["Day"]=="1"
        assert record[0]["PubmedData"]["History"][0][0]["Hour"]=="0"
        assert record[0]["PubmedData"]["History"][0][0]["Minute"]=="0"
        assert record[0]["PubmedData"]["History"][0][1].attributes["PubStatus"]=="medline"
        assert record[0]["PubmedData"]["History"][0][1]["Year"]=="2002"
        assert record[0]["PubmedData"]["History"][0][1]["Month"]=="7"
        assert record[0]["PubmedData"]["History"][0][1]["Day"]=="16"
        assert record[0]["PubmedData"]["History"][0][1]["Hour"]=="10"
        assert record[0]["PubmedData"]["History"][0][1]["Minute"]=="1"
        assert record[0]["PubmedData"]["PublicationStatus"]=="ppublish"
        assert len(record[0]["PubmedData"]["ArticleIdList"])==1
        assert record[0]["PubmedData"]["ArticleIdList"][0]=="12091962"
        assert record[0]["PubmedData"]["ArticleIdList"][0].attributes["IdType"]=="pubmed"
        assert record[1]["MedlineCitation"].attributes["Owner"]=="NLM"
        assert record[1]["MedlineCitation"].attributes["Status"]=="MEDLINE"
        assert record[1]["MedlineCitation"]["PMID"]=="9997"
        assert record[1]["MedlineCitation"]["DateCreated"]["Year"]=="1976"
        assert record[1]["MedlineCitation"]["DateCreated"]["Month"]=="12"
        assert record[1]["MedlineCitation"]["DateCreated"]["Day"]=="30"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Year"]=="1976"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Month"]=="12"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Day"]=="30"
        assert record[1]["MedlineCitation"]["DateRevised"]["Year"]=="2003"
        assert record[1]["MedlineCitation"]["DateRevised"]["Month"]=="11"
        assert record[1]["MedlineCitation"]["DateRevised"]["Day"]=="14"
        assert record[1]["MedlineCitation"]["Article"].attributes["PubModel"]=="Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"]=="0006-3002"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes["IssnType"]=="Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].attributes["CitedMedium"]=="Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"]=="446"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"]=="1"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]=="1976"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Month"]=="Sep"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Day"]=="28"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["Title"]=="Biochimica et biophysica acta"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"]=="Biochim. Biophys. Acta"
        assert record[1]["MedlineCitation"]["Article"]["ArticleTitle"]=="Magnetic studies of Chromatium flavocytochrome C552. A mechanism for heme-flavin interaction."
        assert record[1]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"]=="179-91"
        assert record[1]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]=="Electron paramagnetic resonance and magnetic susceptibility studies of Chromatium flavocytochrome C552 and its diheme flavin-free subunit at temperatures below 45 degrees K are reported. The results show that in the intact protein and the subunit the two low-spin (S = 1/2) heme irons are distinguishable, giving rise to separate EPR signals. In the intact protein only, one of the heme irons exists in two different low spin environments in the pH range 5.5 to 10.5, while the other remains in a constant environment. Factors influencing the variable heme iron environment also influence flavin reactivity, indicating the existence of a mechanism for heme-flavin interaction."
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"]=="Strekas"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"]=="T C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"]=="TC"
        assert record[1]["MedlineCitation"]["Article"]["Language"]==["eng"]
        assert record[1]["MedlineCitation"]["Article"]["PublicationTypeList"]==["Journal Article"]
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["Country"]=="NETHERLANDS"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"]=="Biochim Biophys Acta"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"]=="0217513"
        assert record[1]["MedlineCitation"]["ChemicalList"][0]["RegistryNumber"]=="0"
        assert record[1]["MedlineCitation"]["ChemicalList"][0]["NameOfSubstance"]=="Cytochrome c Group"
        assert record[1]["MedlineCitation"]["ChemicalList"][1]["RegistryNumber"]=="0"
        assert record[1]["MedlineCitation"]["ChemicalList"][1]["NameOfSubstance"]=="Flavins"
        assert record[1]["MedlineCitation"]["ChemicalList"][2]["RegistryNumber"]=="14875-96-8"
        assert record[1]["MedlineCitation"]["ChemicalList"][2]["NameOfSubstance"]=="Heme"
        assert record[1]["MedlineCitation"]["ChemicalList"][3]["RegistryNumber"]=="7439-89-6"
        assert record[1]["MedlineCitation"]["ChemicalList"][3]["NameOfSubstance"]=="Iron"
        assert record[1]["MedlineCitation"]["CitationSubset"]==["IM"]
        assert record[1]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"]=="Binding Sites"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"]=="Chromatium"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0]=="enzymology"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0].attributes["MajorTopicYN"]=="Y"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"]=="Cytochrome c Group"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"].attributes["MajorTopicYN"]=="Y"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"]=="Electron Spin Resonance Spectroscopy"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"]=="Flavins"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"]=="Heme"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"]=="Hydrogen-Ion Concentration"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"]=="Iron"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0]=="analysis"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"]=="Magnetics"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"]=="Oxidation-Reduction"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"]=="Protein Binding"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"]=="Protein Conformation"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][11]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"]=="Temperature"
        assert record[1]["MedlineCitation"]["MeshHeadingList"][12]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[1]["PubmedData"]["History"][0][0].attributes["PubStatus"]=="pubmed"
        assert record[1]["PubmedData"]["History"][0][0]["Year"]=="1976"
        assert record[1]["PubmedData"]["History"][0][0]["Month"]=="9"
        assert record[1]["PubmedData"]["History"][0][0]["Day"]=="28"
        assert record[1]["PubmedData"]["History"][0][1].attributes["PubStatus"]=="medline"
        assert record[1]["PubmedData"]["History"][0][1]["Year"]=="1976"
        assert record[1]["PubmedData"]["History"][0][1]["Month"]=="9"
        assert record[1]["PubmedData"]["History"][0][1]["Day"]=="28"
        assert record[1]["PubmedData"]["History"][0][1]["Hour"]=="0"
        assert record[1]["PubmedData"]["History"][0][1]["Minute"]=="1"
        assert record[1]["PubmedData"]["PublicationStatus"]=="ppublish"
        assert len(record[1]["PubmedData"]["ArticleIdList"])==1
        assert record[1]["PubmedData"]["ArticleIdList"][0]=="9997"
        assert record[1]["PubmedData"]["ArticleIdList"][0].attributes["IdType"]=="pubmed"


    def t_pubmed2(self):
        '''Test parsing XML returned by EFetch, PubMed database (second test)
        '''
        # In PubMed display PMIDs in xml retrieval mode.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='pubmed', id="11748933,11700088",
        #                       retmode="xml")
        input = open('Entrez/pubmed2.xml')
        record = Entrez.read(input)

        assert record[0]["MedlineCitation"].attributes["Owner"]=="NLM"
        assert record[0]["MedlineCitation"].attributes["Status"]=="MEDLINE"
        assert record[0]["MedlineCitation"]["PMID"]=="11748933"
        assert record[0]["MedlineCitation"]["DateCreated"]["Year"]=="2001"
        assert record[0]["MedlineCitation"]["DateCreated"]["Month"]=="12"
        assert record[0]["MedlineCitation"]["DateCreated"]["Day"]=="25"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Year"]=="2002"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Month"]=="03"
        assert record[0]["MedlineCitation"]["DateCompleted"]["Day"]=="04"
        assert record[0]["MedlineCitation"]["DateRevised"]["Year"]=="2006"
        assert record[0]["MedlineCitation"]["DateRevised"]["Month"]=="11"
        assert record[0]["MedlineCitation"]["DateRevised"]["Day"]=="15"
        assert record[0]["MedlineCitation"]["Article"].attributes["PubModel"]=="Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"]=="0011-2240"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes["IssnType"]=="Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].attributes["CitedMedium"]=="Print"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"]=="42"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"]=="4"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]=="2001"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Month"]=="Jun"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["Title"]=="Cryobiology"
        assert record[0]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"]=="Cryobiology"
        assert record[0]["MedlineCitation"]["Article"]["ArticleTitle"]=="Is cryopreservation a homogeneous process? Ultrastructure and motility of untreated, prefreezing, and postthawed spermatozoa of Diplodus puntazzo (Cetti)."
        assert record[0]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"]=="244-55"
        assert record[0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]=="This study subdivides the cryopreservation procedure for Diplodus puntazzo spermatozoa into three key phases, fresh, prefreezing (samples equilibrated in cryosolutions), and postthawed stages, and examines the ultrastructural anomalies and motility profiles of spermatozoa in each stage, with different cryodiluents. Two simple cryosolutions were evaluated: 0.17 M sodium chloride containing a final concentration of 15% dimethyl sulfoxide (Me(2)SO) (cryosolution A) and 0.1 M sodium citrate containing a final concentration of 10% Me(2)SO (cryosolution B). Ultrastructural anomalies of the plasmatic and nuclear membranes of the sperm head were common and the severity of the cryoinjury differed significantly between the pre- and the postfreezing phases and between the two cryosolutions. In spermatozoa diluted with cryosolution A, during the prefreezing phase, the plasmalemma of 61% of the cells was absent or damaged compared with 24% in the fresh sample (P < 0.001). In spermatozoa diluted with cryosolution B, there was a pronounced increase in the number of cells lacking the head plasmatic membrane from the prefreezing to the postthawed stages (from 32 to 52%, P < 0.01). In both cryosolutions, damages to nuclear membrane were significantly higher after freezing (cryosolution A: 8 to 23%, P < 0.01; cryosolution B: 5 to 38%, P < 0.001). With cryosolution A, the after-activation motility profile confirmed a consistent drop from fresh at the prefreezing stage, whereas freezing and thawing did not affect the motility much further and 50% of the cells were immotile by 60-90 s after activation. With cryosolution B, only the postthawing stage showed a sharp drop of motility profile. This study suggests that the different phases of the cryoprocess should be investigated to better understand the process of sperm damage."
        assert record[0]["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"]=="Copyright 2001 Elsevier Science."
        assert record[0]["MedlineCitation"]["Article"]["Affiliation"]==u'Dipartimento di Scienze Ambientali, Universit\xe0 degli Studi della Tuscia, 01100 Viterbo, Italy.'
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"]=="Taddei"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"]=="A R"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"]=="AR"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][1].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"]=="Barbato"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"]=="F"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"]=="F"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][2].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"]=="Abelli"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"]=="L"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"]=="L"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][3].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"]=="Canese"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"]=="S"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"]=="S"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][4].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"]=="Moretti"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"]=="F"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"]=="F"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][5].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"]=="Rana"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"]=="K J"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"]=="KJ"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][6].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["LastName"]=="Fausto"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["ForeName"]=="A M"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][6]["Initials"]=="AM"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][7].attributes["ValidYN"]=="Y"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["LastName"]=="Mazzini"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["ForeName"]=="M"
        assert record[0]["MedlineCitation"]["Article"]["AuthorList"][7]["Initials"]=="M"
        assert record[0]["MedlineCitation"]["Article"]["Language"]==["eng"]
        assert record[0]["MedlineCitation"]["Article"]["PublicationTypeList"][0]=="Journal Article"
        assert record[0]["MedlineCitation"]["Article"]["PublicationTypeList"][1]=="Research Support, Non-U.S. Gov't"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["Country"]=="United States"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"]=="Cryobiology"
        assert record[0]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"]=="0006252"
        assert record[0]["MedlineCitation"]["CitationSubset"]==["IM"]
        assert record[0]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"]=="Animals"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][0]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"]=="Cell Membrane"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0]=="ultrastructure"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][1]["QualifierName"][0].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"]=="Cryopreservation"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["QualifierName"][0]=="methods"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][2]["QualifierName"][0].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"]=="Male"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][3]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"]=="Microscopy, Electron"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][4]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"]=="Microscopy, Electron, Scanning"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][5]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"]=="Nuclear Envelope"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["QualifierName"][0]=="ultrastructure"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][6]["QualifierName"][0].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"]=="Sea Bream"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0]=="anatomy & histology"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][0].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][1]=="physiology"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][7]["QualifierName"][1].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"]=="Semen Preservation"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][0]=="adverse effects"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][0].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][1]=="methods"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][8]["QualifierName"][1].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"]=="Sperm Motility"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][9]["DescriptorName"].attributes["MajorTopicYN"]=="Y"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"]=="Spermatozoa"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["DescriptorName"].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][0]=="physiology"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][0].attributes["MajorTopicYN"]=="N"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][1]=="ultrastructure"
        assert record[0]["MedlineCitation"]["MeshHeadingList"][10]["QualifierName"][1].attributes["MajorTopicYN"]=="Y"
        assert record[0]["PubmedData"]["History"][0][0].attributes["PubStatus"]=="pubmed"
        assert record[0]["PubmedData"]["History"][0][0]["Year"]=="2001"
        assert record[0]["PubmedData"]["History"][0][0]["Month"]=="12"
        assert record[0]["PubmedData"]["History"][0][0]["Day"]=="26"
        assert record[0]["PubmedData"]["History"][0][0]["Hour"]=="10"
        assert record[0]["PubmedData"]["History"][0][0]["Minute"]=="0"
        assert record[0]["PubmedData"]["History"][0][1].attributes["PubStatus"]=="medline"
        assert record[0]["PubmedData"]["History"][0][1]["Year"]=="2002"
        assert record[0]["PubmedData"]["History"][0][1]["Month"]=="3"
        assert record[0]["PubmedData"]["History"][0][1]["Day"]=="5"
        assert record[0]["PubmedData"]["History"][0][1]["Hour"]=="10"
        assert record[0]["PubmedData"]["History"][0][1]["Minute"]=="1"
        assert record[0]["PubmedData"]["PublicationStatus"]=="ppublish"
        assert record[0]["PubmedData"]["ArticleIdList"][0]=="11748933"
        assert record[0]["PubmedData"]["ArticleIdList"][0].attributes["IdType"]=="pubmed"
        assert record[0]["PubmedData"]["ArticleIdList"][1]=="10.1006/cryo.2001.2328"
        assert record[0]["PubmedData"]["ArticleIdList"][1].attributes["IdType"]=="doi"
        assert record[0]["PubmedData"]["ArticleIdList"][2]=="S0011-2240(01)92328-4"
        assert record[0]["PubmedData"]["ArticleIdList"][2].attributes["IdType"]=="pii"

        assert record[1]["MedlineCitation"].attributes["Owner"]=="NLM"
        assert record[1]["MedlineCitation"].attributes["Status"]=="PubMed-not-MEDLINE"
        assert record[1]["MedlineCitation"]["PMID"]=="11700088"
        assert record[1]["MedlineCitation"]["DateCreated"]["Year"]=="2001"
        assert record[1]["MedlineCitation"]["DateCreated"]["Month"]=="11"
        assert record[1]["MedlineCitation"]["DateCreated"]["Day"]=="08"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Year"]=="2001"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Month"]=="12"
        assert record[1]["MedlineCitation"]["DateCompleted"]["Day"]=="20"
        assert record[1]["MedlineCitation"]["DateRevised"]["Year"]=="2003"
        assert record[1]["MedlineCitation"]["DateRevised"]["Month"]=="10"
        assert record[1]["MedlineCitation"]["DateRevised"]["Day"]=="31"
        assert record[1]["MedlineCitation"]["Article"].attributes["PubModel"]=="Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"]=="1090-7807"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISSN"].attributes["IssnType"]=="Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"].attributes["CitedMedium"]=="Print"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"]=="153"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"]=="1"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]=="2001"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Month"]=="Nov"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["Title"]=="Journal of magnetic resonance (San Diego, Calif. : 1997)"
        assert record[1]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"]=="J. Magn. Reson."
        assert record[1]["MedlineCitation"]["Article"]["ArticleTitle"]=="Proton MRI of (13)C distribution by J and chemical shift editing."
        assert record[1]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"]=="117-23"
        assert record[1]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]=="The sensitivity of (13)C NMR imaging can be considerably favored by detecting the (1)H nuclei bound to (13)C nuclei via scalar J-interaction (X-filter). However, the J-editing approaches have difficulty in discriminating between compounds with similar J-constant as, for example, different glucose metabolites. In such cases, it is almost impossible to get J-edited images of a single-compound distribution, since the various molecules are distinguishable only via their chemical shift. In a recent application of J-editing to high-resolution spectroscopy, it has been shown that a more efficient chemical selectivity could be obtained by utilizing the larger chemical shift range of (13)C. This has been made by introducing frequency-selective (13)C pulses that allow a great capability of indirect chemical separation. Here a double-resonance imaging approach is proposed, based on both J-editing and (13)C chemical shift editing, which achieves a powerful chemical selectivity and is able to produce full maps of specific chemical compounds. Results are presented on a multicompartments sample containing solutions of glucose and lactic and glutamic acid in water."
        assert record[1]["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"]=="Copyright 2001 Academic Press."
        assert record[1]["MedlineCitation"]["Article"]["Affiliation"]=="INFM and Department of Physics, University of L'Aquila, I-67100 L'Aquila, Italy."
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"].attributes["CompleteYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0].attributes["ValidYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"]=="Casieri"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["ForeName"]=="C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"]=="C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][1].attributes["ValidYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["LastName"]=="Testa"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["ForeName"]=="C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][1]["Initials"]=="C"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][2].attributes["ValidYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["LastName"]=="Carpinelli"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["ForeName"]=="G"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][2]["Initials"]=="G"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][3].attributes["ValidYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["LastName"]=="Canese"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["ForeName"]=="R"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][3]["Initials"]=="R"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][4].attributes["ValidYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["LastName"]=="Podo"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["ForeName"]=="F"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][4]["Initials"]=="F"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][5].attributes["ValidYN"]=="Y"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["LastName"]=="De Luca"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["ForeName"]=="F"
        assert record[1]["MedlineCitation"]["Article"]["AuthorList"][5]["Initials"]=="F"
        assert record[1]["MedlineCitation"]["Article"]["Language"]==["eng"]
        assert record[1]["MedlineCitation"]["Article"]["PublicationTypeList"][0]=="Journal Article"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["Country"]=="United States"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"]=="J Magn Reson"
        assert record[1]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"]=="9707935"
        assert record[1]["PubmedData"]["History"][0][0].attributes["PubStatus"]=="pubmed"
        assert record[1]["PubmedData"]["History"][0][0]["Year"]=="2001"
        assert record[1]["PubmedData"]["History"][0][0]["Month"]=="11"
        assert record[1]["PubmedData"]["History"][0][0]["Day"]=="9"
        assert record[1]["PubmedData"]["History"][0][0]["Hour"]=="10"
        assert record[1]["PubmedData"]["History"][0][0]["Minute"]=="0"
        assert record[1]["PubmedData"]["History"][0][1].attributes["PubStatus"]=="medline"
        assert record[1]["PubmedData"]["History"][0][1]["Year"]=="2001"
        assert record[1]["PubmedData"]["History"][0][1]["Month"]=="11"
        assert record[1]["PubmedData"]["History"][0][1]["Day"]=="9"
        assert record[1]["PubmedData"]["History"][0][1]["Hour"]=="10"
        assert record[1]["PubmedData"]["History"][0][1]["Minute"]=="1"
        assert record[1]["PubmedData"]["PublicationStatus"]=="ppublish"
        assert record[1]["PubmedData"]["ArticleIdList"][0]=="11700088"
        assert record[1]["PubmedData"]["ArticleIdList"][0].attributes["IdType"]=="pubmed"
        assert record[1]["PubmedData"]["ArticleIdList"][1]=="10.1006/jmre.2001.2429"
        assert record[1]["PubmedData"]["ArticleIdList"][1].attributes["IdType"]=="doi"
        assert record[1]["PubmedData"]["ArticleIdList"][2]=="S1090-7807(01)92429-2"
        assert record[1]["PubmedData"]["ArticleIdList"][2].attributes["IdType"]=="pii"

    def t_journals(self):
        '''Test parsing XML returned by EFetch, Journals database
        '''
        # In Journals display records for journal IDs 22682,21698,1490
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="journals", id=["22682","21698","1490"],
        #                       rettype="full", retmode='xml')
        input = open('Entrez/serialset.xml')
        record = Entrez.read(input)

        assert record[0]["NlmUniqueID"]=="100971611"
        assert record[0]["Title"]=="21st century science & technology"
        assert record[0]["MedlineTA"]=="21st Century Sci Technol"
        assert len(record[0]["PublicationInfo"])==5
        assert record[0]["PublicationInfo"]["Country"]=="United States"
        assert record[0]["PublicationInfo"]["Place"]=="[Washington, D.C. :"
        assert record[0]["PublicationInfo"]["Publisher"]=="21st Century Science Associates,"
        assert record[0]["PublicationInfo"]["PublicationFirstYear"]=="1988"
        assert record[0]["PublicationInfo"]["Frequency"]=="Quarterly,"
        assert record[0]["ISSN"]=="0895-6820"
        assert record[0]["ISSN"].attributes["IssnType"]=="Print"
        assert record[0]["Language"]==["eng"]
        assert record[0]["AcidFreeYN"]=="N"
        assert record[0]["MinorTitleChangeYN"]=="N"
        assert record[0]["CurrentlyIndexedYN"]=="N"
        assert record[0]["IndexOnlineYN"]=="N"
        assert record[0]["IndexingSubset"]=="S"
        assert len(record[0]["CrossReferenceList"])==5
        assert record[0]["CrossReferenceList"][0]["XrTitle"]=="21 century"
        assert record[0]["CrossReferenceList"][0].attributes["XrType"]=="X"
        assert record[0]["CrossReferenceList"][1]["XrTitle"]=="21st century science & technology."
        assert record[0]["CrossReferenceList"][1].attributes["XrType"]=="A"
        assert record[0]["CrossReferenceList"][2]["XrTitle"]=="21st century science and technology"
        assert record[0]["CrossReferenceList"][2].attributes["XrType"]=="X"
        assert record[0]["CrossReferenceList"][3]["XrTitle"]=="Twenty-first century science & technology"
        assert record[0]["CrossReferenceList"][3].attributes["XrType"]=="X"
        assert record[0]["CrossReferenceList"][4]["XrTitle"]=="Twenty-first century science and technology"
        assert record[0]["CrossReferenceList"][4].attributes["XrType"]=="X"
        assert record[0]["SortSerialName"]=="21ST CENTURY SCIENCE & TECHNOLOGY"
        assert record[0]["IlsCreatedTimestamp"]["Year"]=="2000"
        assert record[0]["IlsCreatedTimestamp"]["Month"]=="11"
        assert record[0]["IlsCreatedTimestamp"]["Day"]=="22"
        assert record[0]["IlsUpdatedTimestamp"]["Year"]=="2006"
        assert record[0]["IlsUpdatedTimestamp"]["Month"]=="10"
        assert record[0]["IlsUpdatedTimestamp"]["Day"]=="21"

        assert record[1].attributes["DataCreationMethod"]=="P"
        assert record[1]["NlmUniqueID"]=="100939625"
        assert record[1]["Title"]=="AIHAJ : a journal for the science of occupational and environmental\nhealth and safety"
        assert record[1]["MedlineTA"]=="AIHAJ"
        assert len(record[1]["PublicationInfo"])==6
        assert record[1]["PublicationInfo"]["Country"]=="United States"
        assert record[1]["PublicationInfo"]["Place"]=="Fairfax, VA :"
        assert record[1]["PublicationInfo"]["Publisher"]=="American Industrial Hygiene Association,"
        assert record[1]["PublicationInfo"]["PublicationFirstYear"]=="2000"
        assert record[1]["PublicationInfo"]["PublicationEndYear"]=="2001"
        assert record[1]["PublicationInfo"]["Frequency"]=="Bimonthly"
        assert record[1]["ISSN"]=="1529-8663"
        assert record[1]["ISSN"].attributes["IssnType"]=="Print"
        assert record[1]["Language"]==["eng"]
        assert record[1]["AcidFreeYN"]=="N"
        assert record[1]["ContinuationNotes"]=="Continues: American Industrial Hygiene Association\njournal. Continued by: AIHA journal. "
        assert record[1]["MinorTitleChangeYN"]=="N"
        assert len(record[1]["IndexingHistoryList"])==2
        assert record[1]["IndexingHistoryList"][0].attributes["CitationSubset"]=="IM"
        assert record[1]["IndexingHistoryList"][0].attributes["IndexingTreatment"]=="Full"
        assert record[1]["IndexingHistoryList"][0].attributes["IndexingStatus"]=="Currently-indexed-Title-changed"
        assert record[1]["IndexingHistoryList"][0]["DateOfAction"]["Year"]=="2000"
        assert record[1]["IndexingHistoryList"][0]["DateOfAction"]["Month"]=="03"
        assert record[1]["IndexingHistoryList"][0]["DateOfAction"]["Day"]=="24"
        assert record[1]["IndexingHistoryList"][0]["Coverage"]=="v61n1,Jan./Feb. 2000-v62n6,Nov./Dec. 2001"
        assert record[1]["IndexingHistoryList"][1].attributes["CitationSubset"]=="IM"
        assert record[1]["IndexingHistoryList"][1].attributes["IndexingTreatment"]=="Full"
        assert record[1]["IndexingHistoryList"][1].attributes["IndexingStatus"]=="Continued-by-another-indexed-title"
        assert record[1]["IndexingHistoryList"][1]["DateOfAction"]["Year"]=="2002"
        assert record[1]["IndexingHistoryList"][1]["DateOfAction"]["Month"]=="06"
        assert record[1]["IndexingHistoryList"][1]["DateOfAction"]["Day"]=="03"
        assert record[1]["CurrentlyIndexedYN"]=="N"
        assert record[1]["IndexOnlineYN"]=="N"
        assert record[1]["IndexingSubset"]=="IM"
        assert record[1]["BroadJournalHeadingList"][0]=="Occupational Medicine"
        assert len(record[1]["CrossReferenceList"])==2
        assert record[1]["CrossReferenceList"][0]["XrTitle"]=="AIHAJ :"
        assert record[1]["CrossReferenceList"][0].attributes["XrType"]=="A"
        assert record[1]["CrossReferenceList"][1]["XrTitle"]=="American Industrial Hygiene Association journal"
        assert record[1]["CrossReferenceList"][1].attributes["XrType"]=="X"
        assert record[1]["SortSerialName"]=="AIHAJ : A JOURNAL FOR THE SCIENCE OF OCCUPATIONAL AND\nENVIRONMENTAL HEALTH AND SAFETY"
        assert record[1]["IlsCreatedTimestamp"]["Year"]=="2000"
        assert record[1]["IlsCreatedTimestamp"]["Month"]=="03"
        assert record[1]["IlsCreatedTimestamp"]["Day"]=="22"
        assert record[1]["IlsUpdatedTimestamp"]["Year"]=="2005"
        assert record[1]["IlsUpdatedTimestamp"]["Month"]=="11"
        assert record[1]["IlsUpdatedTimestamp"]["Day"]=="20"

        assert record[2].attributes["DataCreationMethod"]=="P"
        assert record[2]["NlmUniqueID"]=="8403252"
        assert record[2]["Title"]=="Acta crystallographica. Section B, Structural science"
        assert record[2]["MedlineTA"]=="Acta Crystallogr B"
        assert len(record[2]["PublicationInfo"])==5
        assert record[2]["PublicationInfo"]["Country"]=="Denmark"
        assert record[2]["PublicationInfo"]["Place"]=="Copenhagen"
        assert record[2]["PublicationInfo"]["Publisher"]=="Munksgaard International Publishers For The International\nUnion Of Crystallography"
        assert record[2]["PublicationInfo"]["PublicationFirstYear"]=="1983"
        assert record[2]["PublicationInfo"]["Frequency"]=="Bimonthly"
        assert record[2]["ISSN"]=="0108-7681"
        assert record[2]["ISSN"].attributes["IssnType"]=="Print"
        assert record[2]["ISOAbbreviation"]=="Acta Crystallogr., B"
        assert record[2]["Language"]==["eng", "fre", "ger"]
        assert record[2]["AcidFreeYN"]=="N"
        assert record[2]["Coden"]=="ASBSDK"
        assert record[2]["ContinuationNotes"]=="Continues: Acta crystallographica. Section B, Structural\ncrystallography and crystal chemistry. "
        assert record[2]["MinorTitleChangeYN"]=="N"
        assert len(record[2]["IndexingHistoryList"])==1
        assert record[2]["IndexingHistoryList"][0].attributes["CitationSubset"]=="IM"
        assert record[2]["IndexingHistoryList"][0].attributes["IndexingTreatment"]=="Selective"
        assert record[2]["IndexingHistoryList"][0].attributes["IndexingStatus"]=="Currently-indexed"
        assert record[2]["IndexingHistoryList"][0]["DateOfAction"]["Year"]=="1989"
        assert record[2]["IndexingHistoryList"][0]["DateOfAction"]["Month"]=="11"
        assert record[2]["IndexingHistoryList"][0]["DateOfAction"]["Day"]=="06"
        assert record[2]["IndexingHistoryList"][0]["Coverage"]=="v44n1, 1988-"
        assert record[2]["CurrentlyIndexedYN"]=="Y"
        assert record[2]["CurrentlyIndexedForSubset"]==""
        assert record[2]["CurrentlyIndexedForSubset"].attributes["CurrentSubset"]=="IM"
        assert record[2]["CurrentlyIndexedForSubset"].attributes["CurrentIndexingTreatment"]=="Selective"
        assert record[2]["IndexOnlineYN"]=="N"
        assert record[2]["IndexingSubset"]=="IM"
        assert record[2]["BroadJournalHeadingList"][0]=="Chemistry, Analytical"
        assert len(record[2]["CrossReferenceList"])==4
        assert record[2]["CrossReferenceList"][0]["XrTitle"]=="ACTA CRYSTALLOGR B"
        assert record[2]["CrossReferenceList"][0].attributes["XrType"]=="A"
        assert record[2]["CrossReferenceList"][1]["XrTitle"]=="Acta Crystallogr.,Sect.B"
        assert record[2]["CrossReferenceList"][1].attributes["XrType"]=="A"
        assert record[2]["CrossReferenceList"][2]["XrTitle"]=="Acta crystallographica. Section B, Structural\nscience."
        assert record[2]["CrossReferenceList"][2].attributes["XrType"]=="A"
        assert record[2]["CrossReferenceList"][3]["XrTitle"]=="Structural science"
        assert record[2]["CrossReferenceList"][3].attributes["XrType"]=="X"
        assert record[2]["SortSerialName"]=="ACTA CRYSTALLOGRAPHICA. SECTION B, STRUCTURAL\nSCIENCE"
        assert record[2]["IlsCreatedTimestamp"]["Year"]=="1998"
        assert record[2]["IlsCreatedTimestamp"]["Month"]=="11"
        assert record[2]["IlsCreatedTimestamp"]["Day"]=="05"
        assert record[2]["IlsUpdatedTimestamp"]["Year"]=="2008"
        assert record[2]["IlsUpdatedTimestamp"]["Month"]=="04"
        assert record[2]["IlsUpdatedTimestamp"]["Day"]=="04"

    def t_omim(self):
        '''Test parsing XML returned by EFetch, OMIM database
        '''
        # In OMIM show the full record for MIM number 601100 as XML
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="omim", id="601100", retmode='xml',
        #                       rettype='full')
        input = open('Entrez/ncbi_mim.xml')
        record = Entrez.read(input)

        assert len(record)==1
        assert record[0]["Mim-entry_mimNumber"]=="601100"
        assert record[0]["Mim-entry_mimType"]=="1"
        assert record[0]["Mim-entry_mimType"].attributes["value"]=="star"
        assert record[0]["Mim-entry_title"]=="STRESS 70 PROTEIN CHAPERONE, MICROSOME-ASSOCIATED, 60-KD; STCH"
        assert record[0]["Mim-entry_copyright"]=="Copyright (c) 1966-2008 Johns Hopkins University"
        assert record[0]["Mim-entry_symbol"]=="STCH"
        assert record[0]["Mim-entry_locus"]=="21q11.1"
        assert len(record[0]["Mim-entry_text"])==2
        assert record[0]["Mim-entry_text"][0]["Mim-text_label"]=="TEXT"
        assert record[0]["Mim-entry_text"][0]["Mim-text_text"]=="The stress-70 chaperone family consists of proteins that bind to denatured or incorrectly folded polypeptides and play a major role in the processing of cytosolic and secretory proteins. {2:Otterson et al. (1994)} cloned a human cDNA encoding a predicted 471-amino acid protein (60 kD) which they designated STCH. {1:Brodsky et al. (1995)} stated that the protein sequence is very similar to that of HSP70 ({140550}) and BiP ({138120}). As with other members of the family, the STCH protein contains an ATPase domain at the amino terminus whose activity was shown to be independent of peptide stimulation. The protein was found to be microsome-associated and constitutively expressed in all cell types examined."
        assert len(record[0]["Mim-entry_text"][0]["Mim-text_neighbors"])==1
        assert record[0]["Mim-entry_text"][0]["Mim-text_neighbors"]["Mim-link"]["Mim-link_num"]=="30"
        assert record[0]["Mim-entry_text"][0]["Mim-text_neighbors"]["Mim-link"]["Mim-link_uids"]=="8131751,9358068,10675567,9488737,8757872,11048651,2559088,10982831,2105497,16572726,9083109,17181539,14508011,15028727,10651811,9108392,11599566,2661019,11836248,7594475,12406544,8536694,12389629,10430932,9177027,9837933,8522346,2928112,12834280,8702658"
        assert record[0]["Mim-entry_text"][0]["Mim-text_neighbors"]["Mim-link"]["Mim-link_numRelevant"]=="0"
        assert record[0]["Mim-entry_text"][1]["Mim-text_label"]=="TEXT"
        assert record[0]["Mim-entry_text"][1]["Mim-text_text"]=="{1:Brodsky et al. (1995)} mapped the STCH gene to chromosome 21q11.1 with a high-resolution somatic cell hybrid panel for chromosome 21 and by fluorescence in situ hybridization with a YAC containing the gene. By interspecific backcross analysis, {3:Reeves et al. (1998)} mapped the mouse Stch gene to chromosome 16."
        assert len(record[0]["Mim-entry_text"][1]["Mim-text_neighbors"])==1
        assert record[0]["Mim-entry_text"][1]["Mim-text_neighbors"]["Mim-link"]["Mim-link_num"]=="30"
        assert record[0]["Mim-entry_text"][1]["Mim-text_neighbors"]["Mim-link"]["Mim-link_uids"]=="1354597,8244375,8597637,8838809,9143508,1427875,7806216,9852683,7835904,11060461,10083745,7789175,7806232,7513297,8020937,12014109,1769649,2045096,9747039,8034329,8088815,1783375,8275716,8020959,7956352,8020952,10198174,7655454,8750197,11272792"
        assert record[0]["Mim-entry_text"][1]["Mim-text_neighbors"]["Mim-link"]["Mim-link_numRelevant"]=="0"
        assert record[0]["Mim-entry_hasSummary"]==""
        assert record[0]["Mim-entry_hasSummary"].attributes["value"]=="false"
        assert record[0]["Mim-entry_hasSynopsis"]==""
        assert record[0]["Mim-entry_hasSynopsis"].attributes["value"]=="false"
        assert len(record[0]["Mim-entry_editHistory"])==6
        assert record[0]["Mim-entry_editHistory"][0]["Mim-edit-item_author"]=="terry"
        assert record[0]["Mim-entry_editHistory"][0]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"]=="1999"
        assert record[0]["Mim-entry_editHistory"][0]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"]=="3"
        assert record[0]["Mim-entry_editHistory"][0]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"]=="9"
        assert record[0]["Mim-entry_editHistory"][1]["Mim-edit-item_author"]=="carol"
        assert record[0]["Mim-entry_editHistory"][1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"]=="1999"
        assert record[0]["Mim-entry_editHistory"][1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"]=="3"
        assert record[0]["Mim-entry_editHistory"][1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"]=="7"
        assert record[0]["Mim-entry_editHistory"][2]["Mim-edit-item_author"]=="carol"
        assert record[0]["Mim-entry_editHistory"][2]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"]=="1998"
        assert record[0]["Mim-entry_editHistory"][2]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"]=="7"
        assert record[0]["Mim-entry_editHistory"][2]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"]=="8"
        assert record[0]["Mim-entry_editHistory"][3]["Mim-edit-item_author"]=="terry"
        assert record[0]["Mim-entry_editHistory"][3]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"]=="1996"
        assert record[0]["Mim-entry_editHistory"][3]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"]=="5"
        assert record[0]["Mim-entry_editHistory"][3]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"]=="24"
        assert record[0]["Mim-entry_editHistory"][4]["Mim-edit-item_author"]=="mark"
        assert record[0]["Mim-entry_editHistory"][4]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"]=="1996"
        assert record[0]["Mim-entry_editHistory"][4]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"]=="3"
        assert record[0]["Mim-entry_editHistory"][4]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"]=="1"
        assert record[0]["Mim-entry_editHistory"][5]["Mim-edit-item_author"]=="mark"
        assert record[0]["Mim-entry_editHistory"][5]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"]=="1996"
        assert record[0]["Mim-entry_editHistory"][5]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"]=="3"
        assert record[0]["Mim-entry_editHistory"][5]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"]=="1"
        assert record[0]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_author"]=="Alan F. Scott"
        assert record[0]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"]=="1996"
        assert record[0]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"]=="3"
        assert record[0]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"]=="1"
        assert len(record[0]["Mim-entry_references"])==3
        assert record[0]["Mim-entry_references"][0]["Mim-reference_number"]=="1"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_origNumber"]=="1"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_type"]==""
        assert record[0]["Mim-entry_references"][0]["Mim-reference_type"].attributes["value"]=="citation"
        assert len(record[0]["Mim-entry_references"][0]["Mim-reference_authors"])==6
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][0]["Mim-author_name"]=="Brodsky, G."
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][0]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][1]["Mim-author_name"]=="Otterson, G. A."
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][1]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][2]["Mim-author_name"]=="Parry, B. B."
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][2]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][3]["Mim-author_name"]=="Hart, I."
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][3]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][4]["Mim-author_name"]=="Patterson, D."
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][4]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][5]["Mim-author_name"]=="Kaye, F. J."
        assert record[0]["Mim-entry_references"][0]["Mim-reference_authors"][5]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_primaryAuthor"]=="Brodsky"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_otherAuthors"]=="et al."
        assert record[0]["Mim-entry_references"][0]["Mim-reference_citationTitle"]=="Localization of STCH to human chromosome 21q11.1."
        assert record[0]["Mim-entry_references"][0]["Mim-reference_citationType"]=="0"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_volume"]=="30"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_journal"]=="Genomics"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_year"]=="1995"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_month"]=="0"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_day"]=="0"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_pages"][0]["Mim-page_from"]=="627"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_pages"][0]["Mim-page_to"]=="628"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_pubmedUID"]=="8825657"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_ambiguous"]==""
        assert record[0]["Mim-entry_references"][0]["Mim-reference_ambiguous"].attributes["value"]=="false"
        assert record[0]["Mim-entry_references"][0]["Mim-reference_noLink"]==""
        assert record[0]["Mim-entry_references"][0]["Mim-reference_noLink"].attributes["value"]=="false"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_number"]=="2"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_origNumber"]=="2"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_type"]==""
        assert record[0]["Mim-entry_references"][1]["Mim-reference_type"].attributes["value"]=="citation"
        assert len(record[0]["Mim-entry_references"][1]["Mim-reference_authors"])==6
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][0]["Mim-author_name"]=="Otterson, G. A."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][0]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][1]["Mim-author_name"]=="Flynn, G. C."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][1]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][2]["Mim-author_name"]=="Kratzke, R. A."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][2]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][3]["Mim-author_name"]=="Coxon, A."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][3]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][4]["Mim-author_name"]=="Johnston, P. G."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][4]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][5]["Mim-author_name"]=="Kaye, F. J."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_authors"][5]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_primaryAuthor"]=="Otterson"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_otherAuthors"]=="et al."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_citationTitle"]=="Stch encodes the 'ATPase core' of a microsomal stress70 protein."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_citationType"]=="0"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_volume"]=="13"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_journal"]=="EMBO J."
        assert record[0]["Mim-entry_references"][1]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_year"]=="1994"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_month"]=="0"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_day"]=="0"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_pages"][0]["Mim-page_from"]=="1216"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_pages"][0]["Mim-page_to"]=="1225"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_pubmedUID"]=="8131751"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_ambiguous"]==""
        assert record[0]["Mim-entry_references"][1]["Mim-reference_ambiguous"].attributes["value"]=="false"
        assert record[0]["Mim-entry_references"][1]["Mim-reference_noLink"]==""
        assert record[0]["Mim-entry_references"][1]["Mim-reference_noLink"].attributes["value"]=="false"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_number"]=="3"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_origNumber"]=="3"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_type"]==""
        assert record[0]["Mim-entry_references"][2]["Mim-reference_type"].attributes["value"]=="citation"
        assert len(record[0]["Mim-entry_references"][2]["Mim-reference_authors"])==4
        assert record[0]["Mim-entry_references"][2]["Mim-reference_authors"][0]["Mim-author_name"]=="Reeves, R. H."
        assert record[0]["Mim-entry_references"][2]["Mim-reference_authors"][0]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_authors"][1]["Mim-author_name"]=="Rue, E."
        assert record[0]["Mim-entry_references"][2]["Mim-reference_authors"][1]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_authors"][2]["Mim-author_name"]=="Yu, J."
        assert record[0]["Mim-entry_references"][2]["Mim-reference_authors"][2]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_authors"][3]["Mim-author_name"]=="Kao, F.-T."
        assert record[0]["Mim-entry_references"][2]["Mim-reference_authors"][3]["Mim-author_index"]=="1"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_primaryAuthor"]=="Reeves"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_otherAuthors"]=="et al."
        assert record[0]["Mim-entry_references"][2]["Mim-reference_citationTitle"]=="Stch maps to mouse chromosome 16, extending the conserved synteny with human chromosome 21."
        assert record[0]["Mim-entry_references"][2]["Mim-reference_citationType"]=="0"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_volume"]=="49"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_journal"]=="Genomics"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_year"]=="1998"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_month"]=="0"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_day"]=="0"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_pages"][0]["Mim-page_from"]=="156"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_pages"][0]["Mim-page_to"]=="157"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_pubmedUID"]=="9570963"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_ambiguous"]==""
        assert record[0]["Mim-entry_references"][2]["Mim-reference_ambiguous"].attributes["value"]=="false"
        assert record[0]["Mim-entry_references"][2]["Mim-reference_noLink"]==""
        assert record[0]["Mim-entry_references"][2]["Mim-reference_noLink"].attributes["value"]=="false"
        assert record[0]["Mim-entry_attribution"][0]["Mim-edit-item_author"]=="Carol A. Bocchini - updated"
        assert record[0]["Mim-entry_attribution"][0]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"]=="1999"
        assert record[0]["Mim-entry_attribution"][0]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"]=="3"
        assert record[0]["Mim-entry_attribution"][0]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"]=="7"
        assert record[0]["Mim-entry_numGeneMaps"]=="1"
        assert len(record[0]["Mim-entry_medlineLinks"])==1
        assert record[0]["Mim-entry_medlineLinks"]["Mim-link"]["Mim-link_num"]=="3"
        assert record[0]["Mim-entry_medlineLinks"]["Mim-link"]["Mim-link_uids"]=="8825657,8131751,9570963"
        assert record[0]["Mim-entry_medlineLinks"]["Mim-link"]["Mim-link_numRelevant"]=="0"
        assert len(record[0]["Mim-entry_proteinLinks"])==1
        assert record[0]["Mim-entry_proteinLinks"]["Mim-link"]["Mim-link_num"]=="7"
        assert record[0]["Mim-entry_proteinLinks"]["Mim-link"]["Mim-link_uids"]=="148747550,67461586,48928056,30089677,2352621,1351125,460148"
        assert record[0]["Mim-entry_proteinLinks"]["Mim-link"]["Mim-link_numRelevant"]=="0"
        assert len(record[0]["Mim-entry_nucleotideLinks"])==1
        assert record[0]["Mim-entry_nucleotideLinks"]["Mim-link"]["Mim-link_num"]=="5"
        assert record[0]["Mim-entry_nucleotideLinks"]["Mim-link"]["Mim-link_uids"]=="148747549,55741785,48928055,2352620,460147"
        assert record[0]["Mim-entry_nucleotideLinks"]["Mim-link"]["Mim-link_numRelevant"]=="0"

    def t_taxonomy(self):
        '''Test parsing XML returned by EFetch, Taxonomy database
        '''
        # Access the Taxonomy database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db="taxonomy", id="9685", retmode="xml")
        input = open('Entrez/taxonomy.xml')
        record = Entrez.read(input)

        assert len(record)==1
        assert record[0]["TaxId"]=="9685"
        assert record[0]["ScientificName"]=="Felis catus"
        assert record[0]["OtherNames"]["GenbankCommonName"]=="domestic cat"
        assert record[0]["OtherNames"]["Synonym"][0]=="Felis silvestris catus"
        assert record[0]["OtherNames"]["Synonym"][1]=="Felis domesticus"
        assert record[0]["OtherNames"]["CommonName"][0]=="cat"
        assert record[0]["OtherNames"]["CommonName"][1]=="cats"
        assert record[0]["OtherNames"]["Includes"][0]=="Korat cats"
        assert record[0]["ParentTaxId"]=="9682"
        assert record[0]["Rank"]=="species"
        assert record[0]["Division"]=="Mammals"
        assert record[0]["GeneticCode"]["GCId"]=="1"
        assert record[0]["GeneticCode"]["GCName"]=="Standard"
        assert record[0]["MitoGeneticCode"]["MGCId"]=="2"
        assert record[0]["MitoGeneticCode"]["MGCName"]=="Vertebrate Mitochondrial"
        assert record[0]["Lineage"]=="cellular organisms; Eukaryota; Fungi/Metazoa group; Metazoa; Eumetazoa; Bilateria; Coelomata; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Laurasiatheria; Carnivora; Feliformia; Felidae; Felinae; Felis"

        assert record[0]["LineageEx"][0]["TaxId"]=="131567"
        assert record[0]["LineageEx"][0]["ScientificName"]=="cellular organisms"
        assert record[0]["LineageEx"][0]["Rank"]=="no rank"
        assert record[0]["LineageEx"][1]["TaxId"]=="2759"
        assert record[0]["LineageEx"][1]["ScientificName"]=="Eukaryota"
        assert record[0]["LineageEx"][1]["Rank"]=="superkingdom"
        assert record[0]["LineageEx"][2]["TaxId"]=="33154"
        assert record[0]["LineageEx"][2]["ScientificName"]=="Fungi/Metazoa group"
        assert record[0]["LineageEx"][2]["Rank"]=="no rank"
        assert record[0]["LineageEx"][3]["TaxId"]=="33208"
        assert record[0]["LineageEx"][3]["ScientificName"]=="Metazoa"
        assert record[0]["LineageEx"][3]["Rank"]=="kingdom"
        assert record[0]["LineageEx"][4]["TaxId"]=="6072"
        assert record[0]["LineageEx"][4]["ScientificName"]=="Eumetazoa"
        assert record[0]["LineageEx"][4]["Rank"]=="no rank"
        assert record[0]["LineageEx"][5]["TaxId"]=="33213"
        assert record[0]["LineageEx"][5]["ScientificName"]=="Bilateria"
        assert record[0]["LineageEx"][5]["Rank"]=="no rank"
        assert record[0]["LineageEx"][6]["TaxId"]=="33316"
        assert record[0]["LineageEx"][6]["ScientificName"]=="Coelomata"
        assert record[0]["LineageEx"][6]["Rank"]=="no rank"
        assert record[0]["LineageEx"][7]["TaxId"]=="33511"
        assert record[0]["LineageEx"][7]["ScientificName"]=="Deuterostomia"
        assert record[0]["LineageEx"][7]["Rank"]=="no rank"
        assert record[0]["LineageEx"][8]["TaxId"]=="7711"
        assert record[0]["LineageEx"][8]["ScientificName"]=="Chordata"
        assert record[0]["LineageEx"][8]["Rank"]=="phylum"
        assert record[0]["LineageEx"][9]["TaxId"]=="89593"
        assert record[0]["LineageEx"][9]["ScientificName"]=="Craniata"
        assert record[0]["LineageEx"][9]["Rank"]=="subphylum"
        assert record[0]["LineageEx"][10]["TaxId"]=="7742"
        assert record[0]["LineageEx"][10]["ScientificName"]=="Vertebrata"
        assert record[0]["LineageEx"][10]["Rank"]=="no rank"
        assert record[0]["LineageEx"][11]["TaxId"]=="7776"
        assert record[0]["LineageEx"][11]["ScientificName"]=="Gnathostomata"
        assert record[0]["LineageEx"][11]["Rank"]=="superclass"
        assert record[0]["LineageEx"][12]["TaxId"]=="117570"
        assert record[0]["LineageEx"][12]["ScientificName"]=="Teleostomi"
        assert record[0]["LineageEx"][12]["Rank"]=="no rank"
        assert record[0]["LineageEx"][13]["TaxId"]=="117571"
        assert record[0]["LineageEx"][13]["ScientificName"]=="Euteleostomi"
        assert record[0]["LineageEx"][13]["Rank"]=="no rank"
        assert record[0]["LineageEx"][14]["TaxId"]=="8287"
        assert record[0]["LineageEx"][14]["ScientificName"]=="Sarcopterygii"
        assert record[0]["LineageEx"][14]["Rank"]=="no rank"
        assert record[0]["LineageEx"][15]["TaxId"]=="32523"
        assert record[0]["LineageEx"][15]["ScientificName"]=="Tetrapoda"
        assert record[0]["LineageEx"][15]["Rank"]=="no rank"
        assert record[0]["LineageEx"][16]["TaxId"]=="32524"
        assert record[0]["LineageEx"][16]["ScientificName"]=="Amniota"
        assert record[0]["LineageEx"][16]["Rank"]=="no rank"
        assert record[0]["LineageEx"][17]["TaxId"]=="40674"
        assert record[0]["LineageEx"][17]["ScientificName"]=="Mammalia"
        assert record[0]["LineageEx"][17]["Rank"]=="class"
        assert record[0]["LineageEx"][18]["TaxId"]=="32525"
        assert record[0]["LineageEx"][18]["ScientificName"]=="Theria"
        assert record[0]["LineageEx"][18]["Rank"]=="no rank"
        assert record[0]["LineageEx"][19]["TaxId"]=="9347"
        assert record[0]["LineageEx"][19]["ScientificName"]=="Eutheria"
        assert record[0]["LineageEx"][19]["Rank"]=="no rank"
        assert record[0]["LineageEx"][20]["TaxId"]=="314145"
        assert record[0]["LineageEx"][20]["ScientificName"]=="Laurasiatheria"
        assert record[0]["LineageEx"][20]["Rank"]=="superorder"
        assert record[0]["LineageEx"][21]["TaxId"]=="33554"
        assert record[0]["LineageEx"][21]["ScientificName"]=="Carnivora"
        assert record[0]["LineageEx"][21]["Rank"]=="order"
        assert record[0]["LineageEx"][22]["TaxId"]=="379583"
        assert record[0]["LineageEx"][22]["ScientificName"]=="Feliformia"
        assert record[0]["LineageEx"][22]["Rank"]=="suborder"
        assert record[0]["LineageEx"][23]["TaxId"]=="9681"
        assert record[0]["LineageEx"][23]["ScientificName"]=="Felidae"
        assert record[0]["LineageEx"][23]["Rank"]=="family"
        assert record[0]["LineageEx"][24]["TaxId"]=="338152"
        assert record[0]["LineageEx"][24]["ScientificName"]=="Felinae"
        assert record[0]["LineageEx"][24]["Rank"]=="subfamily"
        assert record[0]["LineageEx"][25]["TaxId"]=="9682"
        assert record[0]["LineageEx"][25]["ScientificName"]=="Felis"
        assert record[0]["LineageEx"][25]["Rank"]=="genus"
        assert record[0]["CreateDate"]=="1995/02/27"
        assert record[0]["UpdateDate"]=="2007/09/04"
        assert record[0]["PubDate"]=="1993/07/26"

    def t_nucleotide1(self):
        '''Test parsing XML returned by EFetch, Nucleotide database (first test)
        '''
        # Access the nucleotide database using efetch.
        # To create the XML file, use
        # >>> Bio.Entrez.efetch(db='nucleotide', id=5, retmode='xml')
        input = open('Entrez/nucleotide1.xml')
        record = Entrez.read(input)

        assert record[0]["GBSeq_locus"]=="X60065"
        assert record[0]["GBSeq_length"]=="1136"
        assert record[0]["GBSeq_strandedness"]=="single"
        assert record[0]["GBSeq_moltype"]=="mRNA"
        assert record[0]["GBSeq_topology"]=="linear"
        assert record[0]["GBSeq_division"]=="MAM"
        assert record[0]["GBSeq_update-date"]=="14-NOV-2006"
        assert record[0]["GBSeq_create-date"]=="05-MAY-1992"
        assert record[0]["GBSeq_definition"]=="B.bovis beta-2-gpI mRNA for beta-2-glycoprotein I"
        assert record[0]["GBSeq_primary-accession"]=="X60065"
        assert record[0]["GBSeq_accession-version"]=="X60065.1"
        assert record[0]["GBSeq_other-seqids"][0]=="emb|X60065.1|"
        assert record[0]["GBSeq_other-seqids"][1]=="gi|5"
        assert record[0]["GBSeq_keywords"][0]=="beta-2 glycoprotein I"
        assert record[0]["GBSeq_source"]=="Bos taurus (cattle)"
        assert record[0]["GBSeq_organism"]=="Bos taurus"
        assert record[0]["GBSeq_taxonomy"]=="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Cetartiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos"
        assert record[0]["GBSeq_references"][0]["GBReference_reference"]=="1"
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][0]=="Bendixen,E."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][1]=="Halkier,T."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][2]=="Magnusson,S."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][3]=="Sottrup-Jensen,L."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][4]=="Kristensen,T."
        assert record[0]["GBSeq_references"][0]["GBReference_title"]=="Complete primary structure of bovine beta 2-glycoprotein I: localization of the disulfide bridges"
        assert record[0]["GBSeq_references"][0]["GBReference_journal"]=="Biochemistry 31 (14), 3611-3617 (1992)"
        assert record[0]["GBSeq_references"][0]["GBReference_pubmed"]=="1567819"
        assert record[0]["GBSeq_references"][1]["GBReference_reference"]=="2"
        assert record[0]["GBSeq_references"][1]["GBReference_position"]=="1..1136"
        assert record[0]["GBSeq_references"][1]["GBReference_authors"][0]=="Kristensen,T."
        assert record[0]["GBSeq_references"][1]["GBReference_title"]=="Direct Submission"
        assert record[0]["GBSeq_references"][1]["GBReference_journal"]=="Submitted (11-JUN-1991) T. Kristensen, Dept of Mol Biology, University of Aarhus, C F Mollers Alle 130, DK-8000 Aarhus C, DENMARK"
        assert len(record[0]["GBSeq_feature-table"])==7
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_key"]=="source"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_location"]=="1..1136"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_from"]=="1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_to"]=="1136"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_accession"]=="X60065.1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_name"]=="organism"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_value"]=="Bos taurus"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_name"]=="mol_type"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_value"]=="mRNA"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_value"]=="taxon:9913"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_name"]=="clone"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_value"]=="pBB2I"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_name"]=="tissue_type"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_value"]=="liver"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_key"]=="gene"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_location"]=="<1..1136"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_from"]=="1"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_to"]=="1136"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_accession"]=="X60065.1"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_partial5"]==""
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_partial5"].attributes["value"]=="true"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_name"]=="gene"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_value"]=="beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_key"]=="CDS"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_location"]=="<1..1029"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_from"]=="1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_to"]=="1029"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_accession"]=="X60065.1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_partial5"]==""
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_partial5"].attributes["value"]=="true"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_name"]=="gene"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_value"]=="beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_name"]=="codon_start"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_value"]=="1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_name"]=="transl_table"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_value"]=="1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][3]["GBQualifier_name"]=="product"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][3]["GBQualifier_value"]=="beta-2-glycoprotein I"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][4]["GBQualifier_name"]=="protein_id"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][4]["GBQualifier_value"]=="CAA42669.1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][5]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][5]["GBQualifier_value"]=="GI:6"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][6]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][6]["GBQualifier_value"]=="GOA:P17690"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][7]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][7]["GBQualifier_value"]=="UniProtKB/Swiss-Prot:P17690"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][8]["GBQualifier_name"]=="translation"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][8]["GBQualifier_value"]=="PALVLLLGFLCHVAIAGRTCPKPDELPFSTVVPLKRTYEPGEQIVFSCQPGYVSRGGIRRFTCPLTGLWPINTLKCMPRVCPFAGILENGTVRYTTFEYPNTISFSCHTGFYLKGASSAKCTEEGKWSPDLPVCAPITCPPPPIPKFASLSVYKPLAGNNSFYGSKAVFKCLPHHAMFGNDTVTCTEHGNWTQLPECREVRCPFPSRPDNGFVNHPANPVLYYKDTATFGCHETYSLDGPEEVECSKFGNWSAQPSCKASCKLSIKRATVIYEGERVAIQNKFKNGMLHGQKVSFFCKHKEKKCSYTEDAQCIDGTIEIPKCFKEHSSLAFWKTDASDVKPC"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_key"]=="sig_peptide"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_location"]=="<1..48"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_from"]=="1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_to"]=="48"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_accession"]=="X60065.1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_partial5"]==""
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_partial5"].attributes["value"]=="true"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_name"]=="gene"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_value"]=="beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_key"]=="mat_peptide"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_location"]=="49..1026"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_from"]=="49"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_to"]=="1026"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_intervals"][0]["GBInterval_accession"]=="X60065.1"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][0]["GBQualifier_name"]=="gene"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][0]["GBQualifier_value"]=="beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][1]["GBQualifier_name"]=="product"
        assert record[0]["GBSeq_feature-table"][4]["GBFeature_quals"][1]["GBQualifier_value"]=="beta-2-glycoprotein I"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_key"]=="polyA_signal"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_location"]=="1101..1106"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_from"]=="1101"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_to"]=="1106"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_intervals"][0]["GBInterval_accession"]=="X60065.1"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][0]["GBQualifier_name"]=="gene"
        assert record[0]["GBSeq_feature-table"][5]["GBFeature_quals"][0]["GBQualifier_value"]=="beta-2-gpI"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_key"]=="polyA_site"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_location"]=="1130"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_intervals"][0]["GBInterval_point"]=="1130"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_intervals"][0]["GBInterval_accession"]=="X60065.1"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_quals"][0]["GBQualifier_name"]=="gene"
        assert record[0]["GBSeq_feature-table"][6]["GBFeature_quals"][0]["GBQualifier_value"]=="beta-2-gpI"
        assert record[0]["GBSeq_sequence"]=="ccagcgctcgtcttgctgttggggtttctctgccacgttgctatcgcaggacgaacctgccccaagccagatgagctaccgttttccacggtggttccactgaaacggacctatgagcccggggagcagatagtcttctcctgccagccgggctacgtgtcccggggagggatccggcggtttacatgcccgctcacaggactctggcccatcaacacgctgaaatgcatgcccagagtatgtccttttgctgggatcttagaaaacggaacggtacgctatacaacgtttgagtatcccaacaccatcagcttttcttgccacacggggttttatctgaaaggagctagttctgcaaaatgcactgaggaagggaagtggagcccagaccttcctgtctgtgcccctataacctgccctccaccacccatacccaagtttgcaagtctcagcgtttacaagccgttggctgggaacaactccttctatggcagcaaggcagtctttaagtgcttgccacaccacgcgatgtttggaaatgacaccgttacctgcacggaacatgggaactggacgcagttgccagaatgcagggaagtaagatgcccattcccatcaagaccagacaatgggtttgtgaaccatcctgcaaatccagtgctctactataaggacaccgccacctttggctgccatgaaacgtattccttggatggaccggaagaagtagaatgcagcaaattcggaaactggtctgcacagccaagctgtaaagcatcttgtaagttatctattaaaagagctactgtgatatatgaaggagagagagtagctatccagaacaaatttaagaatggaatgctgcatggccaaaaggtttctttcttctgcaagcataaggaaaagaagtgcagctacacagaagatgctcagtgcatagacggcaccatcgagattcccaaatgcttcaaggagcacagttctttagctttctggaaaacggatgcatctgacgtaaaaccatgctaagctggttttcacactgaaaattaaatgtcatgcttatatgtgtctgtctgagaatctgatggaaacggaaaaataaagagactgaatttaccgtgtcaagaaaaaaa"

    def t_nucleotide2(self):
        '''Test parsing XML returned by EFetch, Nucleotide database (second test)
        '''
        # Access the nucleotide database using efetch.
        # To create the XML file, use
        # >>> Bio.handle = Entrez.efetch(db='nucleotide', id=5,
        #                                rettype='fasta', complexity=0,
        #                                retmode='xml')
        input = open('Entrez/nucleotide2.xml')
        record = Entrez.read(input)

        assert record[0]["TSeq_seqtype"]==""
        assert record[0]["TSeq_seqtype"].attributes["value"]=="nucleotide"
        assert record[0]["TSeq_gi"]=="5"
        assert record[0]["TSeq_accver"]=="X60065.1"
        assert record[0]["TSeq_taxid"]=="9913"
        assert record[0]["TSeq_orgname"]=="Bos taurus"
        assert record[0]["TSeq_defline"]=="B.bovis beta-2-gpI mRNA for beta-2-glycoprotein I"
        assert record[0]["TSeq_length"]=="1136"
        assert record[0]["TSeq_sequence"]=="CCAGCGCTCGTCTTGCTGTTGGGGTTTCTCTGCCACGTTGCTATCGCAGGACGAACCTGCCCCAAGCCAGATGAGCTACCGTTTTCCACGGTGGTTCCACTGAAACGGACCTATGAGCCCGGGGAGCAGATAGTCTTCTCCTGCCAGCCGGGCTACGTGTCCCGGGGAGGGATCCGGCGGTTTACATGCCCGCTCACAGGACTCTGGCCCATCAACACGCTGAAATGCATGCCCAGAGTATGTCCTTTTGCTGGGATCTTAGAAAACGGAACGGTACGCTATACAACGTTTGAGTATCCCAACACCATCAGCTTTTCTTGCCACACGGGGTTTTATCTGAAAGGAGCTAGTTCTGCAAAATGCACTGAGGAAGGGAAGTGGAGCCCAGACCTTCCTGTCTGTGCCCCTATAACCTGCCCTCCACCACCCATACCCAAGTTTGCAAGTCTCAGCGTTTACAAGCCGTTGGCTGGGAACAACTCCTTCTATGGCAGCAAGGCAGTCTTTAAGTGCTTGCCACACCACGCGATGTTTGGAAATGACACCGTTACCTGCACGGAACATGGGAACTGGACGCAGTTGCCAGAATGCAGGGAAGTAAGATGCCCATTCCCATCAAGACCAGACAATGGGTTTGTGAACCATCCTGCAAATCCAGTGCTCTACTATAAGGACACCGCCACCTTTGGCTGCCATGAAACGTATTCCTTGGATGGACCGGAAGAAGTAGAATGCAGCAAATTCGGAAACTGGTCTGCACAGCCAAGCTGTAAAGCATCTTGTAAGTTATCTATTAAAAGAGCTACTGTGATATATGAAGGAGAGAGAGTAGCTATCCAGAACAAATTTAAGAATGGAATGCTGCATGGCCAAAAGGTTTCTTTCTTCTGCAAGCATAAGGAAAAGAAGTGCAGCTACACAGAAGATGCTCAGTGCATAGACGGCACCATCGAGATTCCCAAATGCTTCAAGGAGCACAGTTCTTTAGCTTTCTGGAAAACGGATGCATCTGACGTAAAACCATGCTAAGCTGGTTTTCACACTGAAAATTAAATGTCATGCTTATATGTGTCTGTCTGAGAATCTGATGGAAACGGAAAAATAAAGAGACTGAATTTACCGTGTCAAGAAAAAAA"
        assert record[1]["TSeq_seqtype"]==""
        assert record[1]["TSeq_seqtype"].attributes["value"]=="protein"
        assert record[1]["TSeq_gi"]=="6"
        assert record[1]["TSeq_accver"]=="CAA42669.1"
        assert record[1]["TSeq_taxid"]=="9913"
        assert record[1]["TSeq_orgname"]=="Bos taurus"
        assert record[1]["TSeq_defline"]=="beta-2-glycoprotein I [Bos taurus]"
        assert record[1]["TSeq_length"]=="342"
        assert record[1]["TSeq_sequence"]=="PALVLLLGFLCHVAIAGRTCPKPDELPFSTVVPLKRTYEPGEQIVFSCQPGYVSRGGIRRFTCPLTGLWPINTLKCMPRVCPFAGILENGTVRYTTFEYPNTISFSCHTGFYLKGASSAKCTEEGKWSPDLPVCAPITCPPPPIPKFASLSVYKPLAGNNSFYGSKAVFKCLPHHAMFGNDTVTCTEHGNWTQLPECREVRCPFPSRPDNGFVNHPANPVLYYKDTATFGCHETYSLDGPEEVECSKFGNWSAQPSCKASCKLSIKRATVIYEGERVAIQNKFKNGMLHGQKVSFFCKHKEKKCSYTEDAQCIDGTIEIPKCFKEHSSLAFWKTDASDVKPC"

    def t_nucleotide2(self):
        '''Test parsing XML returned by EFetch, Protein database
        '''
        # Access the protein database using efetch.
        # To create the XML file, use
        # >>> Bio.handle = Entrez.efetch(db='protein', id=8,
        #                                rettype='gp', retmode='xml')
        input = open('Entrez/protein.xml')
        record = Entrez.read(input)

        assert record[0]["GBSeq_locus"]=="CAA35997"
        assert record[0]["GBSeq_length"]=="100"
        assert record[0]["GBSeq_moltype"]=="AA"
        assert record[0]["GBSeq_topology"]=="linear"
        assert record[0]["GBSeq_division"]=="MAM"
        assert record[0]["GBSeq_update-date"]=="12-SEP-1993"
        assert record[0]["GBSeq_create-date"]=="03-APR-1990"
        assert record[0]["GBSeq_definition"]=="unnamed protein product [Bos taurus]"
        assert record[0]["GBSeq_primary-accession"]=="CAA35997"
        assert record[0]["GBSeq_accession-version"]=="CAA35997.1"
        assert record[0]["GBSeq_other-seqids"][0]=="emb|CAA35997.1|"
        assert record[0]["GBSeq_other-seqids"][1]=="gi|8"
        assert record[0]["GBSeq_source"]=="Bos taurus (cattle)"
        assert record[0]["GBSeq_organism"]=="Bos taurus"
        assert record[0]["GBSeq_taxonomy"]=="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Cetartiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos"
        assert record[0]["GBSeq_references"][0]["GBReference_reference"]=="1"
        assert record[0]["GBSeq_references"][0]["GBReference_position"]=="1..100"
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][0]=="Kiefer,M.C."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][1]=="Saphire,A.C.S."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][2]=="Bauer,D.M."
        assert record[0]["GBSeq_references"][0]["GBReference_authors"][3]=="Barr,P.J."
        assert record[0]["GBSeq_references"][0]["GBReference_journal"]=="Unpublished"
        assert record[0]["GBSeq_references"][1]["GBReference_reference"]=="2"
        assert record[0]["GBSeq_references"][1]["GBReference_position"]=="1..100"
        assert record[0]["GBSeq_references"][1]["GBReference_authors"][0]=="Kiefer,M.C."
        assert record[0]["GBSeq_references"][1]["GBReference_title"]=="Direct Submission"
        assert record[0]["GBSeq_references"][1]["GBReference_journal"]=="Submitted (30-JAN-1990) Kiefer M.C., Chiron Corporation, 4560 Hortom St, Emeryville CA 94608-2916, U S A"
        assert record[0]["GBSeq_comment"]=="See <X15699> for Human sequence.~Data kindly reviewed (08-MAY-1990) by Kiefer M.C."
        assert record[0]["GBSeq_source-db"]=="embl accession X51700.1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_key"]=="source"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_location"]=="1..100"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_from"]=="1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_to"]=="100"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_intervals"][0]["GBInterval_accession"]=="CAA35997.1"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_name"]=="organism"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][0]["GBQualifier_value"]=="Bos taurus"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][1]["GBQualifier_value"]=="taxon:9913"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_name"]=="clone"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][2]["GBQualifier_value"]=="bBGP-3"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_name"]=="tissue_type"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][3]["GBQualifier_value"]=="bone matrix"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_name"]=="clone_lib"
        assert record[0]["GBSeq_feature-table"][0]["GBFeature_quals"][4]["GBQualifier_value"]=="Zap-bb"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_key"]=="Protein"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_location"]=="1..100"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_from"]=="1"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_to"]=="100"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_intervals"][0]["GBInterval_accession"]=="CAA35997.1"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_name"]=="name"
        assert record[0]["GBSeq_feature-table"][1]["GBFeature_quals"][0]["GBQualifier_value"]=="unnamed protein product"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_key"]=="Region"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_location"]=="33..97"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_from"]=="33"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_to"]=="97"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_intervals"][0]["GBInterval_accession"]=="CAA35997.1"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_name"]=="region_name"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][0]["GBQualifier_value"]=="Gla"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_name"]=="note"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][1]["GBQualifier_value"]=="Vitamin K-dependent carboxylation/gamma-carboxyglutamic (GLA) domain. This domain is responsible for the high-affinity binding of calcium ions. This domain contains post-translational modifications of many glutamate residues by Vitamin K-dependent...; cl02449"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][2]["GBFeature_quals"][2]["GBQualifier_value"]=="CDD:92835"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_key"]=="CDS"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_location"]=="1..100"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_from"]=="1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_to"]=="100"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_intervals"][0]["GBInterval_accession"]=="CAA35997.1"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_name"]=="coded_by"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][0]["GBQualifier_value"]=="X51700.1:28..330"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][1]["GBQualifier_name"]=="note"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][1]["GBQualifier_value"]=="bone Gla precursor (100 AA)"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][2]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][2]["GBQualifier_value"]=="GOA:P02820"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][3]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][3]["GBQualifier_value"]=="InterPro:IPR000294"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][4]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][4]["GBQualifier_value"]=="InterPro:IPR002384"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][5]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][5]["GBQualifier_value"]=="PDB:1Q3M"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][6]["GBQualifier_name"]=="db_xref"
        assert record[0]["GBSeq_feature-table"][3]["GBFeature_quals"][6]["GBQualifier_value"]=="UniProtKB/Swiss-Prot:P02820"
        assert record[0]["GBSeq_sequence"]=="mrtpmllallalatlclagradakpgdaesgkgaafvskqegsevvkrlrryldhwlgapapypdplepkrevcelnpdcdeladhigfqeayrrfygpv"



if __name__ == '__main__':
    sys.exit(run_tests(sys.argv))
