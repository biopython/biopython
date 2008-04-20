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
    tests = [EInfoTest, ESearchTest, EPostTest, ESummaryTest, ELinkTest]
    
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
        assert record==['pubmed',
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
        assert record['DbName']=='pubmed'
        assert record['MenuName']=='PubMed'
        assert record['Description']=='PubMed bibliographic record'
        assert record['Count']==17905967
        assert record['LastUpdate']=='2008/04/15 06:42'

        assert len(record['FieldList'])==40

        assert record['FieldList'][0]['Name']=='ALL'
        assert record['FieldList'][0]['FullName']=='All Fields'
        assert record['FieldList'][0]['Description']=='All terms from all searchable fields'
        assert record['FieldList'][0]['TermCount']==70792830
        assert record['FieldList'][0]['IsDate']==False
        assert record['FieldList'][0]['IsNumerical']==False
        assert record['FieldList'][0]['SingleToken']==False
        assert record['FieldList'][0]['Hierarchy']==False
        assert record['FieldList'][0]['IsHidden']==False

        assert len(record['LinkList'])==46

        assert record['LinkList'][0]['Name']=='pubmed_books_refs'
        assert record['LinkList'][0]['Menu']=='Cited in Books'
        assert record['LinkList'][0]['Description']=='PubMed links associated with Books'
        assert record['LinkList'][0]['DbTo']=='books'

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
        assert record['Count']==5
        assert record['RetMax']==5
        assert record['RetStart']==0
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
        assert record['TranslationStack'][0]['Count']==5
        assert record['TranslationStack'][0]['Explode']==True
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
        assert record['Count']==10238
        assert record['RetMax']==100
        assert record['RetStart']==0
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
        assert record['TranslationStack'][0]['Count']==52104
        assert record['TranslationStack'][0]['Explode']==True
        assert record['TranslationStack'][1]['Term']=='Medline[SB]'
        assert record['TranslationStack'][1]['Field']=='SB'
        assert record['TranslationStack'][1]['Count']==16509514
        assert record['TranslationStack'][1]['Explode']==True
        assert record['TranslationStack'][2]=='NOT'
        assert record['TranslationStack'][3]=='GROUP'
        assert record['TranslationStack'][4]['Term']=='"neoplasms"[MeSH Terms]'
        assert record['TranslationStack'][4]['Field']=='MeSH Terms'
        assert record['TranslationStack'][4]['Count']==1918010
        assert record['TranslationStack'][4]['Explode']==True
        assert record['TranslationStack'][5]=='OR'
        assert record['TranslationStack'][6]['Term']=='cancer[Text Word]'
        assert record['TranslationStack'][6]['Field']=='Text Word'
        assert record['TranslationStack'][6]['Count']==638849
        assert record['TranslationStack'][6]['Explode']==True
        assert record['TranslationStack'][7]=='OR'
        assert record['TranslationStack'][8]=='GROUP'
        assert record['TranslationStack'][9]['Term']=='2008/02/16[EDAT]'
        assert record['TranslationStack'][9]['Field']=='EDAT'
        assert record['TranslationStack'][9]['Count']==-1
        assert record['TranslationStack'][9]['Explode']==True
        assert record['TranslationStack'][10]['Term']=='2008/04/16[EDAT]'
        assert record['TranslationStack'][10]['Field']=='EDAT'
        assert record['TranslationStack'][10]['Count']==-1
        assert record['TranslationStack'][10]['Explode']==True
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
        assert record['Count']==2652
        assert record['RetMax']==6
        assert record['RetStart']==6
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
        assert record['TranslationStack'][0]['Count']==91806
        assert record['TranslationStack'][0]['Explode']==True
        assert record['TranslationStack'][1]['Term']=='97[vi]'
        assert record['TranslationStack'][1]['Field']=='vi'
        assert record['TranslationStack'][1]['Count']==58681
        assert record['TranslationStack'][1]['Explode']==True
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
        assert record['Count']==177
        assert record['RetMax']==20
        assert record['RetStart']==0
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
        assert record['TranslationStack'][0]['Count']==177
        assert record['TranslationStack'][0]['Explode']==True
        assert record['TranslationStack'][1]=='GROUP'
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
        assert record['Count']==23492
        assert record['RetMax']==20
        assert record['RetStart']==0
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
        assert record['TranslationStack'][0]['Count']==12224
        assert record['TranslationStack'][0]['Explode']==True
        assert record['TranslationStack'][1]['Term']=='stem cells[Acknowledgments]'
        assert record['TranslationStack'][1]['Field']=='Acknowledgments'
        assert record['TranslationStack'][1]['Count']==79
        assert record['TranslationStack'][1]['Explode']==True
        assert record['TranslationStack'][2]=='OR'
        assert record['TranslationStack'][3]['Term']=='stem cells[Figure/Table Caption]'
        assert record['TranslationStack'][3]['Field']=='Figure/Table Caption'
        assert record['TranslationStack'][3]['Count']==806
        assert record['TranslationStack'][3]['Explode']==True
        assert record['TranslationStack'][4]=='OR'
        assert record['TranslationStack'][5]['Term']=='stem cells[Section Title]'
        assert record['TranslationStack'][5]['Field']=='Section Title'
        assert record['TranslationStack'][5]['Count']==522
        assert record['TranslationStack'][5]['Explode']==True
        assert record['TranslationStack'][6]=='OR'
        assert record['TranslationStack'][7]['Term']=='stem cells[Body - All Words]'
        assert record['TranslationStack'][7]['Field']=='Body - All Words'
        assert record['TranslationStack'][7]['Count']==13936
        assert record['TranslationStack'][7]['Explode']==True
        assert record['TranslationStack'][8]=='OR'
        assert record['TranslationStack'][9]['Term']=='stem cells[Title]'
        assert record['TranslationStack'][9]['Field']=='Title'
        assert record['TranslationStack'][9]['Count']==1005
        assert record['TranslationStack'][9]['Explode']==True
        assert record['TranslationStack'][10]=='OR'
        assert record['TranslationStack'][11]['Term']=='stem cells[Abstract]'
        assert record['TranslationStack'][11]['Field']=='Abstract'
        assert record['TranslationStack'][11]['Count']==2503
        assert record['TranslationStack'][11]['Explode']==True
        assert record['TranslationStack'][12]=='OR'
        assert record['TranslationStack'][13]=='GROUP'
        assert record['TranslationStack'][14]['Term']=='free fulltext[filter]'
        assert record['TranslationStack'][14]['Field']=='filter'
        assert record['TranslationStack'][14]['Count']==1412839
        assert record['TranslationStack'][14]['Explode']==True
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
        assert record['Count']==699
        assert record['RetMax']==20
        assert record['RetStart']==0
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
        assert record['Count']==3
        assert record['RetMax']==3
        assert record['RetStart']==0
        assert len(record['IdList'])==3
        assert record['IdList'][0]=='16766766'
        assert record['IdList'][1]=='16422035'
        assert record['IdList'][2]=='4104812'
        assert len(record['TranslationSet'])==0
        assert len(record['TranslationStack'])==2
        assert record['TranslationStack'][0]['Term']=='000200020[molecular weight]'
        assert record['TranslationStack'][0]['Field']=='molecular weight'
        assert record['TranslationStack'][0]['Count']==3
        assert record['TranslationStack'][0]['Explode']==True
        assert record['TranslationStack'][1]=='GROUP'
        assert record['QueryTranslation']=='000200020[molecular weight]'

    def t_notfound(self):
        '''Test parsing XML returned by ESearch when no items were found
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.esearch(db="protein", term="abcXYZ")
        input = open('Entrez/esearch8.xml')
        record = Entrez.read(input)
        assert record['Count']==0
        assert record['RetMax']==0
        assert record['RetStart']==0
        assert len(record['IdList'])==0
        assert len(record['TranslationSet'])==0
        assert record['QueryTranslation']==''
        assert len(record['ErrorList'])==2
        assert "PhraseNotFound" in record['ErrorList']
        assert "FieldNotFound" in record['ErrorList']
        assert len(record['ErrorList']["PhraseNotFound"])==1
        assert record['ErrorList']["PhraseNotFound"][0]=="abcXYZ"
        assert len(record['ErrorList']["FieldNotFound"])==0
        assert len(record['WarningList'])==3
        assert "PhraseIgnored" in record['WarningList']
        assert "QuotedPhraseNotFound" in record['WarningList']
        assert "OutputMessage" in record['WarningList']
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
        except ValueError, exception:
            assert exception.message=="Wrong DB name"
            exception_triggered = True
        assert exception_triggered


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
        assert record[0]["Item"][0]==["PubDate", "Date", "1965 Aug"]
        assert record[0]["Item"][1]==["EPubDate", "Date", ""]
        assert record[0]["Item"][2]==["Source", "String", "Arch Dermatol"]
        assert record[0]["Item"][3][:2]==["AuthorList", "List"]
        assert len(record[0]["Item"][3][2])==2
        assert record[0]["Item"][3][2][0]==["Author", "String", "LoPresti PJ"]
        assert record[0]["Item"][3][2][1]==["Author", "String", "Hambrick GW Jr"]
        assert record[0]["Item"][4]==["LastAuthor", "String", "Hambrick GW Jr"]
        assert record[0]["Item"][5]==["Title", "String", "Zirconium granuloma following treatment of rhus dermatitis."]
        assert record[0]["Item"][6]==["Volume", "String", "92"]
        assert record[0]["Item"][7]==["Issue", "String", "2"]
        assert record[0]["Item"][8]==["Pages", "String", "188-91"]
        assert record[0]["Item"][9][:2]==["LangList", "List"]
        assert len(record[0]["Item"][9][2])==1
        assert record[0]["Item"][9][2][0]==["Lang", "String", "English"]
        assert record[0]["Item"][10]==["NlmUniqueID", "String", "0372433"]
        assert record[0]["Item"][11]==["ISSN", "String", "0003-987X"]
        assert record[0]["Item"][12]==["ESSN", "String", "1538-3652"]
        assert record[0]["Item"][13][:2]==["PubTypeList", "List"]
        assert len(record[0]["Item"][13][2])==1
        assert record[0]["Item"][13][2][0]==["PubType", "String", "Journal Article"]
        assert record[0]["Item"][14]==["RecordStatus", "String", "PubMed - indexed for MEDLINE"]
        assert record[0]["Item"][15]==["PubStatus", "String", "ppublish"]
        assert record[0]["Item"][16][:2]==["ArticleIds", "List"]
        assert len(record[0]["Item"][16][2])==1
        assert record[0]["Item"][16][2][0]==["pubmed", "String", "11850928"]
        assert record[0]["Item"][17][:2]==["History", "List"]
        assert len(record[0]["Item"][17][2])==2
        assert record[0]["Item"][17][2][0]==["pubmed", "Date", "1965/08/01 00:00"]
        assert record[0]["Item"][17][2][1]==["medline", "Date", "2002/03/09 10:01"]
        assert record[0]["Item"][18][:2]==["References", "List"]
        assert len(record[0]["Item"][18][2])==0
        assert record[0]["Item"][19]==["HasAbstract", "Integer", 1]
        assert record[0]["Item"][20]==["PmcRefCount", "Integer", 0]
        assert record[0]["Item"][21]==["FullJournalName", "String", "Archives of dermatology"]
        assert record[0]["Item"][22]==["ELocationID", "String", ""]
        assert record[0]["Item"][23]==["SO", "String", "1965 Aug;92(2):188-91"]

        assert record[1]["Id"]=="11482001"
        assert record[1]["Item"][0]==["PubDate", "Date", "2001 Jun"]
        assert record[1]["Item"][1]==["EPubDate", "Date", ""]
        assert record[1]["Item"][2]==["Source", "String", "Adverse Drug React Toxicol Rev"]
        assert record[1]["Item"][3][:2]==["AuthorList", "List"]
        assert record[1]["Item"][3][2][0]==["Author", "String", "Mantle D"]
        assert record[1]["Item"][3][2][1]==["Author", "String", "Gok MA"]
        assert record[1]["Item"][3][2][2]==["Author", "String", "Lennard TW"]
        assert record[1]["Item"][4]==["LastAuthor", "String", "Lennard TW"]
        assert record[1]["Item"][5]==["Title", "String", "Adverse and beneficial effects of plant extracts on skin and skin disorders."]
        assert record[1]["Item"][6]==["Volume", "String", "20"]
        assert record[1]["Item"][7]==["Issue", "String", "2"]
        assert record[1]["Item"][8]==["Pages", "String", "89-103"]
        assert record[1]["Item"][9][:2]==["LangList", "List"]
        assert len(record[1]["Item"][9][2])==1
        assert record[1]["Item"][9][2][0]==["Lang", "String", "English"]
        assert record[1]["Item"][10]==["NlmUniqueID", "String", "9109474"]
        assert record[1]["Item"][11]==["ISSN", "String", "0964-198X"]
        assert record[1]["Item"][12]==["ESSN", "String", ""]
        assert record[1]["Item"][13][:2]==["PubTypeList", "List"]
        assert len(record[1]["Item"][13][2])==2
        assert record[1]["Item"][13][2][0]==["PubType", "String", "Journal Article"]
        assert record[1]["Item"][13][2][1]==["PubType", "String", "Review"]
        assert record[1]["Item"][14]==["RecordStatus", "String", "PubMed - indexed for MEDLINE"]
        assert record[1]["Item"][15]==["PubStatus", "String", "ppublish"]
        assert record[1]["Item"][16][:2]==["ArticleIds", "List"]
        assert len(record[1]["Item"][16][2])==1
        assert record[1]["Item"][16][2][0]==["pubmed", "String", "11482001"]
        assert record[1]["Item"][17][:2]==["History", "List"]
        assert len(record[1]["Item"][17][2])==2
        assert record[1]["Item"][17][2][0]==["pubmed", "Date", "2001/08/03 10:00"]
        assert record[1]["Item"][17][2][1]==["medline", "Date", "2002/01/23 10:01"]
        assert record[1]["Item"][18][:2]==["References", "List"]
        assert len(record[1]["Item"][18][2])==0
        assert record[1]["Item"][19]==["HasAbstract", "Integer", 1]
        assert record[1]["Item"][20]==["PmcRefCount", "Integer", 0]
        assert record[1]["Item"][21]==["FullJournalName", "String", "Adverse drug reactions and toxicological reviews"]
        assert record[1]["Item"][22]==["ELocationID", "String", ""]
        assert record[1]["Item"][23]==["SO", "String", "2001 Jun;20(2):89-103"]

    def t_journals(self):
        '''Test parsing XML returned by ESummary from the Journals database
        '''
        # In Journals display records for journal IDs 27731,439,735,905 
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="journals", id="27731,439,735,905")
        input = open('Entrez/esummary2.xml')
        record = Entrez.read(input)
        assert record[0]["Id"]=="27731"
        assert record[0]["Item"][0]==["Title", "String", "The American journal of obstetrics and diseases of women and children"]
        assert record[0]["Item"][1]==["MedAbbr", "String", "Am J Obstet Dis Women Child"]
        assert record[0]["Item"][2]==["IsoAbbr", "String", ""]
        assert record[0]["Item"][3]==["NlmId", "String", "14820330R"]
        assert record[0]["Item"][4]==["pISSN", "String", "0894-5543"]
        assert record[0]["Item"][5]==["eISSN", "String", ""]
        assert record[0]["Item"][6]==["PublicationStartYear", "String", "1868"]
        assert record[0]["Item"][7]==["PublicationEndYear", "String", "1919"]
        assert record[0]["Item"][8]==["Publisher", "String", "W.A. Townsend & Adams, $c [1868-1919]"]
        assert record[0]["Item"][9]==["Language", "String", "eng"]
        assert record[0]["Item"][10]==["Country", "String", "United States"]
        assert record[0]["Item"][11][:2]==["BroadHeading", "List"]
        assert len(record[0]["Item"][11][2])==0
        assert record[0]["Item"][12]==["ContinuationNotes", "String", ""]

        assert record[1]["Id"]=="439"
        assert record[1]["Item"][0]==["Title", "String", "American journal of obstetrics and gynecology"]
        assert record[1]["Item"][1]==["MedAbbr", "String", "Am J Obstet Gynecol"]
        assert record[1]["Item"][2]==["IsoAbbr", "String", "Am. J. Obstet. Gynecol."]
        assert record[1]["Item"][3]==["NlmId", "String", "0370476"]
        assert record[1]["Item"][4]==["pISSN", "String", "0002-9378"]
        assert record[1]["Item"][5]==["eISSN", "String", "1097-6868"]
        assert record[1]["Item"][6]==["PublicationStartYear", "String", "1920"]
        assert record[1]["Item"][7]==["PublicationEndYear", "String", ""]
        assert record[1]["Item"][8]==["Publisher", "String", "Elsevier,"]
        assert record[1]["Item"][9]==["Language", "String", "eng"]
        assert record[1]["Item"][10]==["Country", "String", "United States"]
        assert record[1]["Item"][11][:2]==["BroadHeading", "List"]
        assert len(record[1]["Item"][11][2])==2
        assert record[1]["Item"][11][2][0]==["string", "String", "Gynecology"]
        assert record[1]["Item"][11][2][1]==["string", "String", "Obstetrics"]
        assert record[1]["Item"][12]==["ContinuationNotes", "String", "Continues: American journal of obstetrics and diseases of women and children. "]

        assert record[2]["Id"]=="735"
        assert record[2]["Item"][0]==["Title", "String", "Archives of gynecology and obstetrics"]
        assert record[2]["Item"][1]==["MedAbbr", "String", "Arch Gynecol Obstet"]
        assert record[2]["Item"][2]==["IsoAbbr", "String", "Arch. Gynecol. Obstet."]
        assert record[2]["Item"][3]==["NlmId", "String", "8710213"]
        assert record[2]["Item"][4]==["pISSN", "String", "0932-0067"]
        assert record[2]["Item"][5]==["eISSN", "String", "1432-0711"]
        assert record[2]["Item"][6]==["PublicationStartYear", "String", "1987"]
        assert record[2]["Item"][7]==["PublicationEndYear", "String", ""]
        assert record[2]["Item"][8]==["Publisher", "String", "Springer Verlag"]
        assert record[2]["Item"][9]==["Language", "String", "eng"]
        assert record[2]["Item"][10]==["Country", "String", "Germany"]
        assert record[2]["Item"][11][:2]==["BroadHeading", "List"]
        assert len(record[2]["Item"][11][2])==2
        assert record[2]["Item"][11][2][0]==["string", "String", "Gynecology"]
        assert record[2]["Item"][11][2][1]==["string", "String", "Obstetrics"]
        assert record[2]["Item"][12]==["ContinuationNotes", "String", "Continues: Archives of gynecology. "]

        assert record[3]["Id"]=="905"
        assert record[3]["Item"][0]==["Title", "String", "Asia-Oceania journal of obstetrics and gynaecology / AOFOG"]
        assert record[3]["Item"][1]==["MedAbbr", "String", "Asia Oceania J Obstet Gynaecol"]
        assert record[3]["Item"][2]==["IsoAbbr", "String", ""]
        assert record[3]["Item"][3]==["NlmId", "String", "8102781"]
        assert record[3]["Item"][4]==["pISSN", "String", "0389-2328"]
        assert record[3]["Item"][5]==["eISSN", "String", ""]
        assert record[3]["Item"][6]==["PublicationStartYear", "String", "1980"]
        assert record[3]["Item"][7]==["PublicationEndYear", "String", "1994"]
        assert record[3]["Item"][8]==["Publisher", "String", "University Of Tokyo Press"]
        assert record[3]["Item"][9]==["Language", "String", "eng"]
        assert record[3]["Item"][10]==["Country", "String", "Japan"]
        assert record[3]["Item"][11][:2]==["BroadHeading", "List"]
        assert len(record[3]["Item"][11][2])==2
        assert record[3]["Item"][11][2][0]==["string", "String", "Gynecology"]
        assert record[3]["Item"][11][2][1]==["string", "String", "Obstetrics"]
        assert record[3]["Item"][12]==["ContinuationNotes", "String", "Continues: Journal of the Asian Federation of Obstetrics and Gynaecology. Continued by: Journal of obstetrics and gynaecology (Tokyo, Japan). "]

    def t_protein(self):
        '''Test parsing XML returned by ESummary from the Protein database
        '''
        # In Protein display records for GIs 28800982 and 28628843 in xml retrieval mode
        # To create the XML file, use
        # >>> Bio.Entrez.esummary(db="protein", id="28800982,28628843", retmode="xml")
        input = open('Entrez/esummary3.xml')
        record = Entrez.read(input)

        assert record[0]["Id"]=="28800982"
        assert record[0]["Item"][0]==["Caption", "String", "AAO47091"]
        assert record[0]["Item"][1]==["Title", "String", "hemochromatosis [Homo sapiens]"]
        assert record[0]["Item"][2]==["Extra", "String", "gi|28800982|gb|AAO47091.1|[28800982]"]
        assert record[0]["Item"][3]==["Gi", "Integer", 28800982]
        assert record[0]["Item"][4]==["CreateDate", "String", "2003/03/03"]
        assert record[0]["Item"][5]==["UpdateDate", "String", "2003/03/03"]
        assert record[0]["Item"][6]==["Flags", "Integer", 0]
        assert record[0]["Item"][7]==["TaxId", "Integer", 9606]
        assert record[0]["Item"][8]==["Length", "Integer", 268]
        assert record[0]["Item"][9]==["Status", "String", "live"]
        assert record[0]["Item"][10]==["ReplacedBy", "String", ""]
        assert record[0]["Item"][11]==["Comment", "String", "  "]

        assert record[1]["Id"]=="28628843"
        assert record[1]["Item"][0]==["Caption", "String", "AAO49381"]
        assert record[1]["Item"][1]==["Title", "String", "erythroid associated factor [Homo sapiens]"]
        assert record[1]["Item"][2]==["Extra", "String", "gi|28628843|gb|AAO49381.1|AF485325_1[28628843]"]
        assert record[1]["Item"][3]==["Gi", "Integer", 28628843]
        assert record[1]["Item"][4]==["CreateDate", "String", "2003/03/02"]
        assert record[1]["Item"][5]==["UpdateDate", "String", "2003/03/02"]
        assert record[1]["Item"][6]==["Flags", "Integer", 0]
        assert record[1]["Item"][7]==["TaxId", "Integer", 9606]
        assert record[1]["Item"][8]==["Length", "Integer", 102]
        assert record[1]["Item"][9]==["Status", "String", "live"]
        assert record[1]["Item"][10]==["ReplacedBy", "String", ""]
        assert record[1]["Item"][11]==["Comment", "String", "  "]

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
        assert record[0]["Item"][0]==["Caption", "String", "AY207443"]
        assert record[0]["Item"][1]==["Title", "String", "Homo sapiens alpha hemoglobin (HBZP) pseudogene 3' UTR/AluJo repeat breakpoint junction"]
        assert record[0]["Item"][2]==["Extra", "String", "gi|28864546|gb|AY207443.1|[28864546]"]
        assert record[0]["Item"][3]==["Gi", "Integer", 28864546]
        assert record[0]["Item"][4]==["CreateDate", "String", "2003/03/05"]
        assert record[0]["Item"][5]==["UpdateDate", "String", "2003/03/05"]
        assert record[0]["Item"][6]==["Flags", "Integer", 0]
        assert record[0]["Item"][7]==["TaxId", "Integer", 9606]
        assert record[0]["Item"][8]==["Length", "Integer", 491]
        assert record[0]["Item"][9]==["Status", "String", "live"]
        assert record[0]["Item"][10]==["ReplacedBy", "String", ""]
        assert record[0]["Item"][11]==["Comment", "String", "  "]

        assert record[1]["Id"]=="28800981"
        assert record[1]["Item"][0]==["Caption", "String", "AY205604"]
        assert record[1]["Item"][1]==["Title", "String", "Homo sapiens hemochromatosis (HFE) mRNA, partial cds"]
        assert record[1]["Item"][2]==["Extra", "String", "gi|28800981|gb|AY205604.1|[28800981]"]
        assert record[1]["Item"][3]==["Gi", "Integer", 28800981]
        assert record[1]["Item"][4]==["CreateDate", "String", "2003/03/03"]
        assert record[1]["Item"][5]==["UpdateDate", "String", "2003/03/03"]
        assert record[1]["Item"][6]==["Flags", "Integer", 0]
        assert record[1]["Item"][7]==["TaxId", "Integer", 9606]
        assert record[1]["Item"][8]==["Length", "Integer", 860]
        assert record[1]["Item"][9]==["Status", "String", "live"]
        assert record[1]["Item"][10]==["ReplacedBy", "String", ""]
        assert record[1]["Item"][11]==["Comment", "String", "  "]

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
        assert record[0]["Item"][0]==["PdbAcc", "String", "1L5J"]
        assert record[0]["Item"][1]==["PdbDescr", "String", "Crystal Structure Of E. Coli Aconitase B"]
        assert record[0]["Item"][2]==["EC", "String", "4.2.1.3"]
        assert record[0]["Item"][3]==["Resolution", "String", "2.4"]
        assert record[0]["Item"][4]==["ExpMethod", "String", "X-Ray Diffraction"]
        assert record[0]["Item"][5]==["PdbClass", "String", "Lyase"]
        assert record[0]["Item"][6]==["PdbReleaseDate", "String", "2007/8/27"]
        assert record[0]["Item"][7]==["PdbDepositDate", "String", "2002/3/7"]
        assert record[0]["Item"][8]==["DepositDate", "String", "2007/10/25"]
        assert record[0]["Item"][9]==["ModifyDate", "String", "2007/10/25"]
        assert record[0]["Item"][10]==["LigCode", "String", "F3S|TRA"]
        assert record[0]["Item"][11]==["LigCount", "String", "2"]
        assert record[0]["Item"][12]==["ModProteinResCount", "String", "0"]
        assert record[0]["Item"][13]==["ModDNAResCount", "String", "0"]
        assert record[0]["Item"][14]==["ModRNAResCount", "String", "0"]
        assert record[0]["Item"][15]==["ProteinChainCount", "String", "2"]
        assert record[0]["Item"][16]==["DNAChainCount", "String", "0"]
        assert record[0]["Item"][17]==["RNAChainCount", "String", "0"]

        assert record[1]["Id"]=="12120"
        assert record[1]["Item"][0]==["PdbAcc", "String", "1B0K"]
        assert record[1]["Item"][1]==["PdbDescr", "String", "S642a:fluorocitrate Complex Of Aconitase"]
        assert record[1]["Item"][2]==["EC", "String", "4.2.1.3"]
        assert record[1]["Item"][3]==["Resolution", "String", "2.5"]
        assert record[1]["Item"][4]==["ExpMethod", "String", "X-Ray Diffraction"]
        assert record[1]["Item"][5]==["PdbClass", "String", "Lyase"]
        assert record[1]["Item"][6]==["PdbReleaseDate", "String", "2007/8/27"]
        assert record[1]["Item"][7]==["PdbDepositDate", "String", "1998/11/11"]
        assert record[1]["Item"][8]==["DepositDate", "String", "2007/10/07"]
        assert record[1]["Item"][9]==["ModifyDate", "String", "2007/10/07"]
        assert record[1]["Item"][10]==["LigCode", "String", "FLC|O|SF4"]
        assert record[1]["Item"][11]==["LigCount", "String", "3"]
        assert record[1]["Item"][12]==["ModProteinResCount", "String", "0"]
        assert record[1]["Item"][13]==["ModDNAResCount", "String", "0"]
        assert record[1]["Item"][14]==["ModRNAResCount", "String", "0"]
        assert record[1]["Item"][15]==["ProteinChainCount", "String", "1"]
        assert record[1]["Item"][16]==["DNAChainCount", "String", "0"]
        assert record[1]["Item"][17]==["RNAChainCount", "String", "0"]

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
        assert record[0]["Item"][0]==["Rank", "String", "species"]
        assert record[0]["Item"][1]==["Division", "String", "even-toed ungulates"]
        assert record[0]["Item"][2]==["ScientificName", "String", "Bos taurus"]
        assert record[0]["Item"][3]==["CommonName", "String", "cattle"]
        assert record[0]["Item"][4]==["TaxId", "Integer", 9913]
        assert record[0]["Item"][5]==["NucNumber", "Integer", 2264214]
        assert record[0]["Item"][6]==["ProtNumber", "Integer", 55850]
        assert record[0]["Item"][7]==["StructNumber", "Integer", 1517]
        assert record[0]["Item"][8]==["GenNumber", "Integer", 31]
        assert record[0]["Item"][9]==["GeneNumber", "Integer", 29651]
        assert record[0]["Item"][10]==["Genus", "String", ""]
        assert record[0]["Item"][11]==["Species", "String", ""]
        assert record[0]["Item"][12]==["Subsp", "String", ""]

        assert record[1]["Id"]=="30521"
        assert record[1]["Item"][0]==["Rank", "String", "species"]
        assert record[1]["Item"][1]==["Division", "String", "even-toed ungulates"]
        assert record[1]["Item"][2]==["ScientificName", "String", "Bos grunniens"]
        assert record[1]["Item"][3]==["CommonName", "String", "domestic yak"]
        assert record[1]["Item"][4]==["TaxId", "Integer", 30521]
        assert record[1]["Item"][5]==["NucNumber", "Integer", 560]
        assert record[1]["Item"][6]==["ProtNumber", "Integer", 254]
        assert record[1]["Item"][7]==["StructNumber", "Integer", 0]
        assert record[1]["Item"][8]==["GenNumber", "Integer", 1]
        assert record[1]["Item"][9]==["GeneNumber", "Integer", 13]
        assert record[1]["Item"][10]==["Genus", "String", ""]
        assert record[1]["Item"][11]==["Species", "String", ""]
        assert record[1]["Item"][12]==["Subsp", "String", ""]

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
        assert record[0]["Item"][0]==["Marker_Name", "String", "SE234324"]
        assert record[0]["Item"][1][:2]==["Map_Gene_Summary_List", "List"]
        assert len(record[0]["Item"][1][2])==1
        assert record[0]["Item"][1][2][0][:2]==["Map_Gene_Summary", "Structure"]
        assert len(record[0]["Item"][1][2][0][2])==3
        assert record[0]["Item"][1][2][0][2][0]==["Org", "String", "Sus scrofa"]
        assert record[0]["Item"][1][2][0][2][1]==["Chr", "String", " chromosome 7"]
        assert record[0]["Item"][1][2][0][2][2]==["Locus", "String", ""]
        assert record[0]["Item"][2]==["EPCR_Summary", "String", "Found by e-PCR in sequences from Sus scrofa."]
        assert record[0]["Item"][3]==["LocusId", "String", ""]

        assert record[1]["Id"]=="254086"
        assert record[1]["Item"][0]==["Marker_Name", "String", "SE259162"]
        assert record[1]["Item"][1][:2]==["Map_Gene_Summary_List", "List"]
        assert len(record[1]["Item"][1][2])==1
        assert record[1]["Item"][1][2][0][:2]==["Map_Gene_Summary", "Structure"]
        assert len(record[1]["Item"][1][2][0][2])==3
        assert record[1]["Item"][1][2][0][2][0]==["Org", "String", "Sus scrofa"]
        assert record[1]["Item"][1][2][0][2][1]==["Chr", "String", " chromosome 12"]
        assert record[1]["Item"][1][2][0][2][2]==["Locus", "String", ""]
        assert record[1]["Item"][2]==["EPCR_Summary", "String", "Found by e-PCR in sequences from Sus scrofa."]
        assert record[1]["Item"][3]==["LocusId", "String", ""]

    def t_wrong(self):
        '''Test parsing XML returned by ESummary with incorrect arguments
        '''
        # To create the XML file, use
        # >>> Bio.Entrez.esummary()
        input = open('Entrez/esummary8.xml')
        exception_triggered = False
        try:
            record = Entrez.read(input)
        except ValueError, exception:
            assert exception.message=="Neither query_key nor id specified"
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
        assert record[0]["LinkSetDb"]["DbTo"]=="pubmed"
        assert record[0]["LinkSetDb"]["LinkName"]=="pubmed_pubmed"
        assert len(record[0]["LinkSetDb"]["Link"])==144
        assert record[0]["LinkSetDb"]["Link"][0]["Id"]=="9298984"
        assert record[0]["LinkSetDb"]["Link"][0]["Score"]==2147483647
        assert record[0]["LinkSetDb"]["Link"][1]["Id"]=="8794856"
        assert record[0]["LinkSetDb"]["Link"][1]["Score"]==65259341
        assert record[0]["LinkSetDb"]["Link"][2]["Id"]=="9700164"
        assert record[0]["LinkSetDb"]["Link"][2]["Score"]==60347327
        assert record[0]["LinkSetDb"]["Link"][3]["Id"]=="7914521"
        assert record[0]["LinkSetDb"]["Link"][3]["Score"]==54343405
        assert record[0]["LinkSetDb"]["Link"][4]["Id"]=="1339459"
        assert record[0]["LinkSetDb"]["Link"][4]["Score"]==53014422
        assert record[0]["LinkSetDb"]["Link"][5]["Id"]=="9914369"
        assert record[0]["LinkSetDb"]["Link"][5]["Score"]==52741538
        assert record[0]["LinkSetDb"]["Link"][6]["Id"]=="11590237"
        assert record[0]["LinkSetDb"]["Link"][6]["Score"]==52493903
        assert record[0]["LinkSetDb"]["Link"][7]["Id"]=="12686595"
        assert record[0]["LinkSetDb"]["Link"][7]["Score"]==48734007
        assert record[0]["LinkSetDb"]["Link"][8]["Id"]=="9074495"
        assert record[0]["LinkSetDb"]["Link"][8]["Score"]==48220447
        assert record[0]["LinkSetDb"]["Link"][9]["Id"]=="11146659"
        assert record[0]["LinkSetDb"]["Link"][9]["Score"]==46604626
        assert record[0]["LinkSetDb"]["Link"][10]["Id"]=="10893249"
        assert record[0]["LinkSetDb"]["Link"][10]["Score"]==46254167
        assert record[0]["LinkSetDb"]["Link"][11]["Id"]=="8978614"
        assert record[0]["LinkSetDb"]["Link"][11]["Score"]==46166362
        assert record[0]["LinkSetDb"]["Link"][12]["Id"]=="15371539"
        assert record[0]["LinkSetDb"]["Link"][12]["Score"]==45060488
        assert record[0]["LinkSetDb"]["Link"][13]["Id"]=="10806105"
        assert record[0]["LinkSetDb"]["Link"][13]["Score"]==44825774
        assert record[0]["LinkSetDb"]["Link"][14]["Id"]=="10402457"
        assert record[0]["LinkSetDb"]["Link"][14]["Score"]==44338522
        assert record[0]["LinkSetDb"]["Link"][15]["Id"]=="10545493"
        assert record[0]["LinkSetDb"]["Link"][15]["Score"]==43860609
        assert record[0]["LinkSetDb"]["Link"][16]["Id"]=="10523511"
        assert record[0]["LinkSetDb"]["Link"][16]["Score"]==43268800
        assert record[0]["LinkSetDb"]["Link"][17]["Id"]=="12515822"
        assert record[0]["LinkSetDb"]["Link"][17]["Score"]==43215343
        assert record[0]["LinkSetDb"]["Link"][18]["Id"]=="15915585"
        assert record[0]["LinkSetDb"]["Link"][18]["Score"]==43029760
        assert record[0]["LinkSetDb"]["Link"][19]["Id"]=="11483958"
        assert record[0]["LinkSetDb"]["Link"][19]["Score"]==42348877
        assert record[0]["LinkSetDb"]["Link"][20]["Id"]=="11685532"
        assert record[0]["LinkSetDb"]["Link"][20]["Score"]==42262104
        assert record[0]["LinkSetDb"]["Link"][21]["Id"]=="9869638"
        assert record[0]["LinkSetDb"]["Link"][21]["Score"]==41143593
        assert record[0]["LinkSetDb"]["Link"][22]["Id"]=="12080088"
        assert record[0]["LinkSetDb"]["Link"][22]["Score"]==40849490
        assert record[0]["LinkSetDb"]["Link"][23]["Id"]=="12034769"
        assert record[0]["LinkSetDb"]["Link"][23]["Score"]==40841328
        assert record[0]["LinkSetDb"]["Link"][24]["Id"]=="9852156"
        assert record[0]["LinkSetDb"]["Link"][24]["Score"]==40793501
        assert record[0]["LinkSetDb"]["Link"][25]["Id"]=="9735366"
        assert record[0]["LinkSetDb"]["Link"][25]["Score"]==40661605
        assert record[0]["LinkSetDb"]["Link"][26]["Id"]=="10749938"
        assert record[0]["LinkSetDb"]["Link"][26]["Score"]==40486739
        assert record[0]["LinkSetDb"]["Link"][27]["Id"]=="9490715"
        assert record[0]["LinkSetDb"]["Link"][27]["Score"]==40311339
        assert record[0]["LinkSetDb"]["Link"][28]["Id"]=="9425896"
        assert record[0]["LinkSetDb"]["Link"][28]["Score"]==40056298
        assert record[0]["LinkSetDb"]["Link"][29]["Id"]=="11266459"
        assert record[0]["LinkSetDb"]["Link"][29]["Score"]==39883140
        assert record[0]["LinkSetDb"]["Link"][30]["Id"]=="14522947"
        assert record[0]["LinkSetDb"]["Link"][30]["Score"]==39683976
        assert record[0]["LinkSetDb"]["Link"][31]["Id"]=="15616189"
        assert record[0]["LinkSetDb"]["Link"][31]["Score"]==39518630
        assert record[0]["LinkSetDb"]["Link"][32]["Id"]=="16732327"
        assert record[0]["LinkSetDb"]["Link"][32]["Score"]==39425668
        assert record[0]["LinkSetDb"]["Link"][33]["Id"]=="11179694"
        assert record[0]["LinkSetDb"]["Link"][33]["Score"]==39183565
        assert record[0]["LinkSetDb"]["Link"][34]["Id"]=="10898791"
        assert record[0]["LinkSetDb"]["Link"][34]["Score"]==39159761
        assert record[0]["LinkSetDb"]["Link"][35]["Id"]=="11146661"
        assert record[0]["LinkSetDb"]["Link"][35]["Score"]==39116609
        assert record[0]["LinkSetDb"]["Link"][36]["Id"]=="11914278"
        assert record[0]["LinkSetDb"]["Link"][36]["Score"]==39028004
        assert record[0]["LinkSetDb"]["Link"][37]["Id"]=="10985388"
        assert record[0]["LinkSetDb"]["Link"][37]["Score"]==39002572
        assert record[0]["LinkSetDb"]["Link"][38]["Id"]=="16839185"
        assert record[0]["LinkSetDb"]["Link"][38]["Score"]==38916726
        assert record[0]["LinkSetDb"]["Link"][39]["Id"]=="7585942"
        assert record[0]["LinkSetDb"]["Link"][39]["Score"]==38747288
        assert record[0]["LinkSetDb"]["Link"][40]["Id"]=="2022189"
        assert record[0]["LinkSetDb"]["Link"][40]["Score"]==38717145
        assert record[0]["LinkSetDb"]["Link"][41]["Id"]=="7690762"
        assert record[0]["LinkSetDb"]["Link"][41]["Score"]==38647275
        assert record[0]["LinkSetDb"]["Link"][42]["Id"]=="7904902"
        assert record[0]["LinkSetDb"]["Link"][42]["Score"]==38557343
        assert record[0]["LinkSetDb"]["Link"][43]["Id"]=="9378750"
        assert record[0]["LinkSetDb"]["Link"][43]["Score"]==38484849
        assert record[0]["LinkSetDb"]["Link"][44]["Id"]=="12388768"
        assert record[0]["LinkSetDb"]["Link"][44]["Score"]==38454422
        assert record[0]["LinkSetDb"]["Link"][45]["Id"]=="11352945"
        assert record[0]["LinkSetDb"]["Link"][45]["Score"]==38449836
        assert record[0]["LinkSetDb"]["Link"][46]["Id"]=="11267866"
        assert record[0]["LinkSetDb"]["Link"][46]["Score"]==38419058
        assert record[0]["LinkSetDb"]["Link"][47]["Id"]=="17222555"
        assert record[0]["LinkSetDb"]["Link"][47]["Score"]==38368546
        assert record[0]["LinkSetDb"]["Link"][48]["Id"]=="11252055"
        assert record[0]["LinkSetDb"]["Link"][48]["Score"]==38257516
        assert record[0]["LinkSetDb"]["Link"][49]["Id"]=="16585270"
        assert record[0]["LinkSetDb"]["Link"][49]["Score"]==37800856
        assert record[0]["LinkSetDb"]["Link"][50]["Id"]=="9606208"
        assert record[0]["LinkSetDb"]["Link"][50]["Score"]==37669054
        assert record[0]["LinkSetDb"]["Link"][51]["Id"]=="17182852"
        assert record[0]["LinkSetDb"]["Link"][51]["Score"]==37621285
        assert record[0]["LinkSetDb"]["Link"][52]["Id"]=="9933569"
        assert record[0]["LinkSetDb"]["Link"][52]["Score"]==37398470
        assert record[0]["LinkSetDb"]["Link"][53]["Id"]=="15268859"
        assert record[0]["LinkSetDb"]["Link"][53]["Score"]==37340582
        assert record[0]["LinkSetDb"]["Link"][54]["Id"]=="12235289"
        assert record[0]["LinkSetDb"]["Link"][54]["Score"]==37247450
        assert record[0]["LinkSetDb"]["Link"][55]["Id"]=="16741559"
        assert record[0]["LinkSetDb"]["Link"][55]["Score"]==37198716
        assert record[0]["LinkSetDb"]["Link"][56]["Id"]=="11266451"
        assert record[0]["LinkSetDb"]["Link"][56]["Score"]==37142542
        assert record[0]["LinkSetDb"]["Link"][57]["Id"]=="15075237"
        assert record[0]["LinkSetDb"]["Link"][57]["Score"]==36897578
        assert record[0]["LinkSetDb"]["Link"][58]["Id"]=="15485811"
        assert record[0]["LinkSetDb"]["Link"][58]["Score"]==36804297
        assert record[0]["LinkSetDb"]["Link"][59]["Id"]=="14699129"
        assert record[0]["LinkSetDb"]["Link"][59]["Score"]==36782062
        assert record[0]["LinkSetDb"]["Link"][60]["Id"]=="16510521"
        assert record[0]["LinkSetDb"]["Link"][60]["Score"]==36724370
        assert record[0]["LinkSetDb"]["Link"][61]["Id"]=="15824131"
        assert record[0]["LinkSetDb"]["Link"][61]["Score"]==36695341
        assert record[0]["LinkSetDb"]["Link"][62]["Id"]=="15371340"
        assert record[0]["LinkSetDb"]["Link"][62]["Score"]==36685694
        assert record[0]["LinkSetDb"]["Link"][63]["Id"]=="9878245"
        assert record[0]["LinkSetDb"]["Link"][63]["Score"]==36684230
        assert record[0]["LinkSetDb"]["Link"][64]["Id"]=="10398680"
        assert record[0]["LinkSetDb"]["Link"][64]["Score"]==36573411
        assert record[0]["LinkSetDb"]["Link"][65]["Id"]=="16516834"
        assert record[0]["LinkSetDb"]["Link"][65]["Score"]==36525654
        assert record[0]["LinkSetDb"]["Link"][66]["Id"]=="11715021"
        assert record[0]["LinkSetDb"]["Link"][66]["Score"]==36518129
        assert record[0]["LinkSetDb"]["Link"][67]["Id"]=="14622138"
        assert record[0]["LinkSetDb"]["Link"][67]["Score"]==36496009
        assert record[0]["LinkSetDb"]["Link"][68]["Id"]=="11092768"
        assert record[0]["LinkSetDb"]["Link"][68]["Score"]==36457186
        assert record[0]["LinkSetDb"]["Link"][69]["Id"]=="12514103"
        assert record[0]["LinkSetDb"]["Link"][69]["Score"]==36385909
        assert record[0]["LinkSetDb"]["Link"][70]["Id"]=="17525528"
        assert record[0]["LinkSetDb"]["Link"][70]["Score"]==36316439
        assert record[0]["LinkSetDb"]["Link"][71]["Id"]=="11402064"
        assert record[0]["LinkSetDb"]["Link"][71]["Score"]==36172957
        assert record[0]["LinkSetDb"]["Link"][72]["Id"]=="9258677"
        assert record[0]["LinkSetDb"]["Link"][72]["Score"]==35989143
        assert record[0]["LinkSetDb"]["Link"][73]["Id"]=="14499625"
        assert record[0]["LinkSetDb"]["Link"][73]["Score"]==35978627
        assert record[0]["LinkSetDb"]["Link"][74]["Id"]=="10428958"
        assert record[0]["LinkSetDb"]["Link"][74]["Score"]==35924800
        assert record[0]["LinkSetDb"]["Link"][75]["Id"]=="14972679"
        assert record[0]["LinkSetDb"]["Link"][75]["Score"]==35915578
        assert record[0]["LinkSetDb"]["Link"][76]["Id"]=="9396743"
        assert record[0]["LinkSetDb"]["Link"][76]["Score"]==35883749
        assert record[0]["LinkSetDb"]["Link"][77]["Id"]=="16219694"
        assert record[0]["LinkSetDb"]["Link"][77]["Score"]==35870689
        assert record[0]["LinkSetDb"]["Link"][78]["Id"]=="11369198"
        assert record[0]["LinkSetDb"]["Link"][78]["Score"]==35838048
        assert record[0]["LinkSetDb"]["Link"][79]["Id"]=="17333235"
        assert record[0]["LinkSetDb"]["Link"][79]["Score"]==35815282
        assert record[0]["LinkSetDb"]["Link"][80]["Id"]=="11102811"
        assert record[0]["LinkSetDb"]["Link"][80]["Score"]==35783566
        assert record[0]["LinkSetDb"]["Link"][81]["Id"]=="10207147"
        assert record[0]["LinkSetDb"]["Link"][81]["Score"]==35594009
        assert record[0]["LinkSetDb"]["Link"][82]["Id"]=="10477755"
        assert record[0]["LinkSetDb"]["Link"][82]["Score"]==35589601
        assert record[0]["LinkSetDb"]["Link"][83]["Id"]=="10747094"
        assert record[0]["LinkSetDb"]["Link"][83]["Score"]==35548072
        assert record[0]["LinkSetDb"]["Link"][84]["Id"]=="15215209"
        assert record[0]["LinkSetDb"]["Link"][84]["Score"]==35526869
        assert record[0]["LinkSetDb"]["Link"][85]["Id"]=="11157774"
        assert record[0]["LinkSetDb"]["Link"][85]["Score"]==35510607
        assert record[0]["LinkSetDb"]["Link"][86]["Id"]=="10669599"
        assert record[0]["LinkSetDb"]["Link"][86]["Score"]==35462246
        assert record[0]["LinkSetDb"]["Link"][87]["Id"]=="17448445"
        assert record[0]["LinkSetDb"]["Link"][87]["Score"]==35398470
        assert record[0]["LinkSetDb"]["Link"][88]["Id"]=="17878237"
        assert record[0]["LinkSetDb"]["Link"][88]["Score"]==35231311
        assert record[0]["LinkSetDb"]["Link"][89]["Id"]=="10411903"
        assert record[0]["LinkSetDb"]["Link"][89]["Score"]==35202708
        assert record[0]["LinkSetDb"]["Link"][90]["Id"]=="12773390"
        assert record[0]["LinkSetDb"]["Link"][90]["Score"]==35171743
        assert record[0]["LinkSetDb"]["Link"][91]["Id"]=="12498686"
        assert record[0]["LinkSetDb"]["Link"][91]["Score"]==35131906
        assert record[0]["LinkSetDb"]["Link"][92]["Id"]=="9009204"
        assert record[0]["LinkSetDb"]["Link"][92]["Score"]==34993776
        assert record[0]["LinkSetDb"]["Link"][93]["Id"]=="17576797"
        assert record[0]["LinkSetDb"]["Link"][93]["Score"]==34988639
        assert record[0]["LinkSetDb"]["Link"][94]["Id"]=="10225945"
        assert record[0]["LinkSetDb"]["Link"][94]["Score"]==34950419
        assert record[0]["LinkSetDb"]["Link"][95]["Id"]=="11161560"
        assert record[0]["LinkSetDb"]["Link"][95]["Score"]==34912466
        assert record[0]["LinkSetDb"]["Link"][96]["Id"]=="11967147"
        assert record[0]["LinkSetDb"]["Link"][96]["Score"]==34900540
        assert record[0]["LinkSetDb"]["Link"][97]["Id"]=="14711415"
        assert record[0]["LinkSetDb"]["Link"][97]["Score"]==34883714
        assert record[0]["LinkSetDb"]["Link"][98]["Id"]=="2211824"
        assert record[0]["LinkSetDb"]["Link"][98]["Score"]==34843507
        assert record[0]["LinkSetDb"]["Link"][99]["Id"]=="15737064"
        assert record[0]["LinkSetDb"]["Link"][99]["Score"]==34828187
        assert record[0]["LinkSetDb"]["Link"][100]["Id"]=="7720068"
        assert record[0]["LinkSetDb"]["Link"][100]["Score"]==34811182
        assert record[0]["LinkSetDb"]["Link"][101]["Id"]=="9472001"
        assert record[0]["LinkSetDb"]["Link"][101]["Score"]==34799321
        assert record[0]["LinkSetDb"]["Link"][102]["Id"]=="11792803"
        assert record[0]["LinkSetDb"]["Link"][102]["Score"]==34697393
        assert record[0]["LinkSetDb"]["Link"][103]["Id"]=="11386760"
        assert record[0]["LinkSetDb"]["Link"][103]["Score"]==34684610
        assert record[0]["LinkSetDb"]["Link"][104]["Id"]=="15094189"
        assert record[0]["LinkSetDb"]["Link"][104]["Score"]==34684021
        assert record[0]["LinkSetDb"]["Link"][105]["Id"]=="9763420"
        assert record[0]["LinkSetDb"]["Link"][105]["Score"]==34666950
        assert record[0]["LinkSetDb"]["Link"][106]["Id"]=="10629219"
        assert record[0]["LinkSetDb"]["Link"][106]["Score"]==34422925
        assert record[0]["LinkSetDb"]["Link"][107]["Id"]=="11238410"
        assert record[0]["LinkSetDb"]["Link"][107]["Score"]==34318521
        assert record[0]["LinkSetDb"]["Link"][108]["Id"]=="17199038"
        assert record[0]["LinkSetDb"]["Link"][108]["Score"]==34255594
        assert record[0]["LinkSetDb"]["Link"][109]["Id"]=="12944469"
        assert record[0]["LinkSetDb"]["Link"][109]["Score"]==34249507
        assert record[0]["LinkSetDb"]["Link"][110]["Id"]=="15616192"
        assert record[0]["LinkSetDb"]["Link"][110]["Score"]==34110517
        assert record[0]["LinkSetDb"]["Link"][111]["Id"]=="11146660"
        assert record[0]["LinkSetDb"]["Link"][111]["Score"]==34063257
        assert record[0]["LinkSetDb"]["Link"][112]["Id"]=="11402066"
        assert record[0]["LinkSetDb"]["Link"][112]["Score"]==34012520
        assert record[0]["LinkSetDb"]["Link"][113]["Id"]=="6791901"
        assert record[0]["LinkSetDb"]["Link"][113]["Score"]==33311119
        assert record[0]["LinkSetDb"]["Link"][114]["Id"]=="7172865"
        assert record[0]["LinkSetDb"]["Link"][114]["Score"]==32934223
        assert record[0]["LinkSetDb"]["Link"][115]["Id"]=="8270646"
        assert record[0]["LinkSetDb"]["Link"][115]["Score"]==32898701
        assert record[0]["LinkSetDb"]["Link"][116]["Id"]=="1916263"
        assert record[0]["LinkSetDb"]["Link"][116]["Score"]==32707765
        assert record[0]["LinkSetDb"]["Link"][117]["Id"]=="7588080"
        assert record[0]["LinkSetDb"]["Link"][117]["Score"]==32503526
        assert record[0]["LinkSetDb"]["Link"][118]["Id"]=="7391142"
        assert record[0]["LinkSetDb"]["Link"][118]["Score"]==31632645
        assert record[0]["LinkSetDb"]["Link"][119]["Id"]=="6793236"
        assert record[0]["LinkSetDb"]["Link"][119]["Score"]==31522175
        assert record[0]["LinkSetDb"]["Link"][120]["Id"]=="2512302"
        assert record[0]["LinkSetDb"]["Link"][120]["Score"]==30339372
        assert record[0]["LinkSetDb"]["Link"][121]["Id"]=="7720069"
        assert record[0]["LinkSetDb"]["Link"][121]["Score"]==30024525
        assert record[0]["LinkSetDb"]["Link"][122]["Id"]=="8257792"
        assert record[0]["LinkSetDb"]["Link"][122]["Score"]==29834355
        assert record[0]["LinkSetDb"]["Link"][123]["Id"]=="3417141"
        assert record[0]["LinkSetDb"]["Link"][123]["Score"]==27920818
        assert record[0]["LinkSetDb"]["Link"][124]["Id"]=="3315496"
        assert record[0]["LinkSetDb"]["Link"][124]["Score"]==27422009
        assert record[0]["LinkSetDb"]["Link"][125]["Id"]=="1993311"
        assert record[0]["LinkSetDb"]["Link"][125]["Score"]==26763828
        assert record[0]["LinkSetDb"]["Link"][126]["Id"]=="6185450"
        assert record[0]["LinkSetDb"]["Link"][126]["Score"]==26100420
        assert record[0]["LinkSetDb"]["Link"][127]["Id"]=="1819515"
        assert record[0]["LinkSetDb"]["Link"][127]["Score"]==26036804
        assert record[0]["LinkSetDb"]["Link"][128]["Id"]=="7250964"
        assert record[0]["LinkSetDb"]["Link"][128]["Score"]==25738652
        assert record[0]["LinkSetDb"]["Link"][129]["Id"]=="8489280"
        assert record[0]["LinkSetDb"]["Link"][129]["Score"]==25587858
        assert record[0]["LinkSetDb"]["Link"][130]["Id"]=="7096444"
        assert record[0]["LinkSetDb"]["Link"][130]["Score"]==24642544
        assert record[0]["LinkSetDb"]["Link"][131]["Id"]=="348629"
        assert record[0]["LinkSetDb"]["Link"][131]["Score"]==24432498
        assert record[0]["LinkSetDb"]["Link"][132]["Id"]=="2275018"
        assert record[0]["LinkSetDb"]["Link"][132]["Score"]==23077593
        assert record[0]["LinkSetDb"]["Link"][133]["Id"]=="1747872"
        assert record[0]["LinkSetDb"]["Link"][133]["Score"]==22933494
        assert record[0]["LinkSetDb"]["Link"][134]["Id"]=="3547036"
        assert record[0]["LinkSetDb"]["Link"][134]["Score"]==22925639
        assert record[0]["LinkSetDb"]["Link"][135]["Id"]=="18291669"
        assert record[0]["LinkSetDb"]["Link"][135]["Score"]==22762310
        assert record[0]["LinkSetDb"]["Link"][136]["Id"]=="1576878"
        assert record[0]["LinkSetDb"]["Link"][136]["Score"]==20846041
        assert record[0]["LinkSetDb"]["Link"][137]["Id"]=="6230555"
        assert record[0]["LinkSetDb"]["Link"][137]["Score"]==19354488
        assert record[0]["LinkSetDb"]["Link"][138]["Id"]=="7627547"
        assert record[0]["LinkSetDb"]["Link"][138]["Score"]==18940607
        assert record[0]["LinkSetDb"]["Link"][139]["Id"]=="17678444"
        assert record[0]["LinkSetDb"]["Link"][139]["Score"]==18834135
        assert record[0]["LinkSetDb"]["Link"][140]["Id"]=="3366468"
        assert record[0]["LinkSetDb"]["Link"][140]["Score"]==14831756
        assert record[0]["LinkSetDb"]["Link"][141]["Id"]=="1959920"
        assert record[0]["LinkSetDb"]["Link"][141]["Score"]==14156414
        assert record[0]["LinkSetDb"]["Link"][142]["Id"]=="13242628"
        assert record[0]["LinkSetDb"]["Link"][142]["Score"]==12584732
        assert record[0]["LinkSetDb"]["Link"][143]["Id"]=="17248312"
        assert record[0]["LinkSetDb"]["Link"][143]["Score"]==7610436

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
        assert record[0]["LinkSetDb"]["DbTo"]=="pubmed"
        assert record[0]["LinkSetDb"]["LinkName"]=="pubmed_pubmed"
        assert record[0]["LinkSetDb"]["Link"][0]["Id"]=="11812492"
        assert record[0]["LinkSetDb"]["Link"][0]["Score"]==2147483647
        assert record[0]["LinkSetDb"]["Link"][1]["Id"]=="11774222"
        assert record[0]["LinkSetDb"]["Link"][1]["Score"]==2147483647
        assert record[0]["LinkSetDb"]["Link"][2]["Id"]=="11668631"
        assert record[0]["LinkSetDb"]["Link"][2]["Score"]==86345306
        assert record[0]["LinkSetDb"]["Link"][3]["Id"]=="15111095"
        assert record[0]["LinkSetDb"]["Link"][3]["Score"]==81604359
        assert record[0]["LinkSetDb"]["Link"][4]["Id"]=="10731564"
        assert record[0]["LinkSetDb"]["Link"][4]["Score"]==65665112
        assert record[0]["LinkSetDb"]["Link"][5]["Id"]=="15780005"
        assert record[0]["LinkSetDb"]["Link"][5]["Score"]==62251079
        assert record[0]["LinkSetDb"]["Link"][6]["Id"]=="17885136"
        assert record[0]["LinkSetDb"]["Link"][6]["Score"]==50322134
        assert record[0]["LinkSetDb"]["Link"][7]["Id"]=="17470297"
        assert record[0]["LinkSetDb"]["Link"][7]["Score"]==49148434
        assert record[0]["LinkSetDb"]["Link"][8]["Id"]=="16005284"
        assert record[0]["LinkSetDb"]["Link"][8]["Score"]==49035508
        assert record[0]["LinkSetDb"]["Link"][9]["Id"]=="10856373"
        assert record[0]["LinkSetDb"]["Link"][9]["Score"]==48363137
        assert record[0]["LinkSetDb"]["Link"][10]["Id"]=="15383292"
        assert record[0]["LinkSetDb"]["Link"][10]["Score"]==48347159
        assert record[0]["LinkSetDb"]["Link"][11]["Id"]=="17040125"
        assert record[0]["LinkSetDb"]["Link"][11]["Score"]==48301243
        assert record[0]["LinkSetDb"]["Link"][12]["Id"]=="10770808"
        assert record[0]["LinkSetDb"]["Link"][12]["Score"]==47696325
        assert record[0]["LinkSetDb"]["Link"][13]["Id"]=="11125122"
        assert record[0]["LinkSetDb"]["Link"][13]["Score"]==45889695
        assert record[0]["LinkSetDb"]["Link"][14]["Id"]=="15287587"
        assert record[0]["LinkSetDb"]["Link"][14]["Score"]==45599733
        assert record[0]["LinkSetDb"]["Link"][15]["Id"]=="15839745"
        assert record[0]["LinkSetDb"]["Link"][15]["Score"]==44650620
        assert record[0]["LinkSetDb"]["Link"][16]["Id"]=="10612825"
        assert record[0]["LinkSetDb"]["Link"][16]["Score"]==44445812
        assert record[0]["LinkSetDb"]["Link"][17]["Id"]=="15024419"
        assert record[0]["LinkSetDb"]["Link"][17]["Score"]==44075047
        assert record[0]["LinkSetDb"]["Link"][18]["Id"]=="12743802"
        assert record[0]["LinkSetDb"]["Link"][18]["Score"]==43873158
        assert record[0]["LinkSetDb"]["Link"][19]["Id"]=="15238684"
        assert record[0]["LinkSetDb"]["Link"][19]["Score"]==43856864
        assert record[0]["LinkSetDb"]["Link"][20]["Id"]=="12386340"
        assert record[0]["LinkSetDb"]["Link"][20]["Score"]==43770229
        assert record[0]["LinkSetDb"]["Link"][21]["Id"]=="16269725"
        assert record[0]["LinkSetDb"]["Link"][21]["Score"]==43712594
        assert record[0]["LinkSetDb"]["Link"][22]["Id"]=="10592273"
        assert record[0]["LinkSetDb"]["Link"][22]["Score"]==43640108
        assert record[0]["LinkSetDb"]["Link"][23]["Id"]=="15383308"
        assert record[0]["LinkSetDb"]["Link"][23]["Score"]==42835474
        assert record[0]["LinkSetDb"]["Link"][24]["Id"]=="15676075"
        assert record[0]["LinkSetDb"]["Link"][24]["Score"]==42272663
        assert record[0]["LinkSetDb"]["Link"][25]["Id"]=="11774221"
        assert record[0]["LinkSetDb"]["Link"][25]["Score"]==42058380
        assert record[0]["LinkSetDb"]["Link"][26]["Id"]=="10592272"
        assert record[0]["LinkSetDb"]["Link"][26]["Score"]==41719917
        assert record[0]["LinkSetDb"]["Link"][27]["Id"]=="15997407"
        assert record[0]["LinkSetDb"]["Link"][27]["Score"]==41535461
        assert record[0]["LinkSetDb"]["Link"][28]["Id"]=="15774024"
        assert record[0]["LinkSetDb"]["Link"][28]["Score"]==41351079
        assert record[0]["LinkSetDb"]["Link"][29]["Id"]=="11233160"
        assert record[0]["LinkSetDb"]["Link"][29]["Score"]==41268965
        assert record[0]["LinkSetDb"]["Link"][30]["Id"]=="14702162"
        assert record[0]["LinkSetDb"]["Link"][30]["Score"]==41147661
        assert record[0]["LinkSetDb"]["Link"][31]["Id"]=="16616613"
        assert record[0]["LinkSetDb"]["Link"][31]["Score"]==41073100
        assert record[0]["LinkSetDb"]["Link"][32]["Id"]=="17202370"
        assert record[0]["LinkSetDb"]["Link"][32]["Score"]==40819600
        assert record[0]["LinkSetDb"]["Link"][33]["Id"]=="15478601"
        assert record[0]["LinkSetDb"]["Link"][33]["Score"]==40578911
        assert record[0]["LinkSetDb"]["Link"][34]["Id"]=="15322925"
        assert record[0]["LinkSetDb"]["Link"][34]["Score"]==40548101
        assert record[0]["LinkSetDb"]["Link"][35]["Id"]=="11472559"
        assert record[0]["LinkSetDb"]["Link"][35]["Score"]==40508356
        assert record[0]["LinkSetDb"]["Link"][36]["Id"]=="11925998"
        assert record[0]["LinkSetDb"]["Link"][36]["Score"]==39844751
        assert record[0]["LinkSetDb"]["Link"][37]["Id"]=="12372145"
        assert record[0]["LinkSetDb"]["Link"][37]["Score"]==39809277
        assert record[0]["LinkSetDb"]["Link"][38]["Id"]=="17562224"
        assert record[0]["LinkSetDb"]["Link"][38]["Score"]==38850094
        assert record[0]["LinkSetDb"]["Link"][39]["Id"]=="15037105"
        assert record[0]["LinkSetDb"]["Link"][39]["Score"]==38758229
        assert record[0]["LinkSetDb"]["Link"][40]["Id"]=="14998511"
        assert record[0]["LinkSetDb"]["Link"][40]["Score"]==38608049
        assert record[0]["LinkSetDb"]["Link"][41]["Id"]=="10092480"
        assert record[0]["LinkSetDb"]["Link"][41]["Score"]==38410463
        assert record[0]["LinkSetDb"]["Link"][42]["Id"]=="7729881"
        assert record[0]["LinkSetDb"]["Link"][42]["Score"]==38329800
        assert record[0]["LinkSetDb"]["Link"][43]["Id"]=="12933853"
        assert record[0]["LinkSetDb"]["Link"][43]["Score"]==37881850
        assert record[0]["LinkSetDb"]["Link"][44]["Id"]=="16818783"
        assert record[0]["LinkSetDb"]["Link"][44]["Score"]==37835096
        assert record[0]["LinkSetDb"]["Link"][45]["Id"]=="16406333"
        assert record[0]["LinkSetDb"]["Link"][45]["Score"]==37775136
        assert record[0]["LinkSetDb"]["Link"][46]["Id"]=="11472553"
        assert record[0]["LinkSetDb"]["Link"][46]["Score"]==37750745
        assert record[0]["LinkSetDb"]["Link"][47]["Id"]=="11403387"
        assert record[0]["LinkSetDb"]["Link"][47]["Score"]==37707525
        assert record[0]["LinkSetDb"]["Link"][48]["Id"]=="17306254"
        assert record[0]["LinkSetDb"]["Link"][48]["Score"]==37685833
        assert record[0]["LinkSetDb"]["Link"][49]["Id"]=="11516587"
        assert record[0]["LinkSetDb"]["Link"][49]["Score"]==37620966
        assert record[0]["LinkSetDb"]["Link"][50]["Id"]=="9274032"
        assert record[0]["LinkSetDb"]["Link"][50]["Score"]==37528832
        assert record[0]["LinkSetDb"]["Link"][51]["Id"]=="12856318"
        assert record[0]["LinkSetDb"]["Link"][51]["Score"]==37484650
        assert record[0]["LinkSetDb"]["Link"][52]["Id"]=="14695526"
        assert record[0]["LinkSetDb"]["Link"][52]["Score"]==37429895
        assert record[0]["LinkSetDb"]["Link"][53]["Id"]=="12481045"
        assert record[0]["LinkSetDb"]["Link"][53]["Score"]==37051674
        assert record[0]["LinkSetDb"]["Link"][54]["Id"]=="11752345"
        assert record[0]["LinkSetDb"]["Link"][54]["Score"]==36875760
        assert record[0]["LinkSetDb"]["Link"][55]["Id"]=="12467974"
        assert record[0]["LinkSetDb"]["Link"][55]["Score"]==36787103
        assert record[0]["LinkSetDb"]["Link"][56]["Id"]=="11214099"
        assert record[0]["LinkSetDb"]["Link"][56]["Score"]==36710749
        assert record[0]["LinkSetDb"]["Link"][57]["Id"]=="14638788"
        assert record[0]["LinkSetDb"]["Link"][57]["Score"]==36667774
        assert record[0]["LinkSetDb"]["Link"][58]["Id"]=="16278157"
        assert record[0]["LinkSetDb"]["Link"][58]["Score"]==36598908
        assert record[0]["LinkSetDb"]["Link"][59]["Id"]=="11752242"
        assert record[0]["LinkSetDb"]["Link"][59]["Score"]==36555638
        assert record[0]["LinkSetDb"]["Link"][60]["Id"]=="14681474"
        assert record[0]["LinkSetDb"]["Link"][60]["Score"]==36317853
        assert record[0]["LinkSetDb"]["Link"][61]["Id"]=="15944077"
        assert record[0]["LinkSetDb"]["Link"][61]["Score"]==36264027
        assert record[0]["LinkSetDb"]["Link"][62]["Id"]=="12625936"
        assert record[0]["LinkSetDb"]["Link"][62]["Score"]==36088314
        assert record[0]["LinkSetDb"]["Link"][63]["Id"]=="16672453"
        assert record[0]["LinkSetDb"]["Link"][63]["Score"]==35985060
        assert record[0]["LinkSetDb"]["Link"][64]["Id"]=="14695451"
        assert record[0]["LinkSetDb"]["Link"][64]["Score"]==35971708
        assert record[0]["LinkSetDb"]["Link"][65]["Id"]=="12402526"
        assert record[0]["LinkSetDb"]["Link"][65]["Score"]==35942170
        assert record[0]["LinkSetDb"]["Link"][66]["Id"]=="10592200"
        assert record[0]["LinkSetDb"]["Link"][66]["Score"]==35932875
        assert record[0]["LinkSetDb"]["Link"][67]["Id"]=="17584494"
        assert record[0]["LinkSetDb"]["Link"][67]["Score"]==35869907
        assert record[0]["LinkSetDb"]["Link"][68]["Id"]=="17761848"
        assert record[0]["LinkSetDb"]["Link"][68]["Score"]==35868206
        assert record[0]["LinkSetDb"]["Link"][69]["Id"]=="16697384"
        assert record[0]["LinkSetDb"]["Link"][69]["Score"]==35792791
        assert record[0]["LinkSetDb"]["Link"][70]["Id"]=="8784774"
        assert record[0]["LinkSetDb"]["Link"][70]["Score"]==35787497
        assert record[0]["LinkSetDb"]["Link"][71]["Id"]=="18000556"
        assert record[0]["LinkSetDb"]["Link"][71]["Score"]==35701408
        assert record[0]["LinkSetDb"]["Link"][72]["Id"]=="15828434"
        assert record[0]["LinkSetDb"]["Link"][72]["Score"]==35604052
        assert record[0]["LinkSetDb"]["Link"][73]["Id"]=="10511685"
        assert record[0]["LinkSetDb"]["Link"][73]["Score"]==35598319
        assert record[0]["LinkSetDb"]["Link"][74]["Id"]=="15608284"
        assert record[0]["LinkSetDb"]["Link"][74]["Score"]==35439627
        assert record[0]["LinkSetDb"]["Link"][75]["Id"]=="11125071"
        assert record[0]["LinkSetDb"]["Link"][75]["Score"]==35414962
        assert record[0]["LinkSetDb"]["Link"][76]["Id"]=="11791238"
        assert record[0]["LinkSetDb"]["Link"][76]["Score"]==35411948
        assert record[0]["LinkSetDb"]["Link"][77]["Id"]=="15710433"
        assert record[0]["LinkSetDb"]["Link"][77]["Score"]==35197152
        assert record[0]["LinkSetDb"]["Link"][78]["Id"]=="16164550"
        assert record[0]["LinkSetDb"]["Link"][78]["Score"]==35172458
        assert record[0]["LinkSetDb"]["Link"][79]["Id"]=="17697334"
        assert record[0]["LinkSetDb"]["Link"][79]["Score"]==35121478
        assert record[0]["LinkSetDb"]["Link"][80]["Id"]=="12537121"
        assert record[0]["LinkSetDb"]["Link"][80]["Score"]==35054632
        assert record[0]["LinkSetDb"]["Link"][81]["Id"]=="12860672"
        assert record[0]["LinkSetDb"]["Link"][81]["Score"]==35046651
        assert record[0]["LinkSetDb"]["Link"][82]["Id"]=="15630619"
        assert record[0]["LinkSetDb"]["Link"][82]["Score"]==35034076
        assert record[0]["LinkSetDb"]["Link"][83]["Id"]=="15125639"
        assert record[0]["LinkSetDb"]["Link"][83]["Score"]==35007338
        assert record[0]["LinkSetDb"]["Link"][84]["Id"]=="11443570"
        assert record[0]["LinkSetDb"]["Link"][84]["Score"]==34935553
        assert record[0]["LinkSetDb"]["Link"][85]["Id"]=="12208043"
        assert record[0]["LinkSetDb"]["Link"][85]["Score"]==34923107
        assert record[0]["LinkSetDb"]["Link"][86]["Id"]=="11731507"
        assert record[0]["LinkSetDb"]["Link"][86]["Score"]==34875290
        assert record[0]["LinkSetDb"]["Link"][87]["Id"]=="11988510"
        assert record[0]["LinkSetDb"]["Link"][87]["Score"]==34773036
        assert record[0]["LinkSetDb"]["Link"][88]["Id"]=="11125038"
        assert record[0]["LinkSetDb"]["Link"][88]["Score"]==34754724
        assert record[0]["LinkSetDb"]["Link"][89]["Id"]=="16381944"
        assert record[0]["LinkSetDb"]["Link"][89]["Score"]==34747225
        assert record[0]["LinkSetDb"]["Link"][90]["Id"]=="17135206"
        assert record[0]["LinkSetDb"]["Link"][90]["Score"]==34735015
        assert record[0]["LinkSetDb"]["Link"][91]["Id"]=="17099226"
        assert record[0]["LinkSetDb"]["Link"][91]["Score"]==34698054
        assert record[0]["LinkSetDb"]["Link"][92]["Id"]=="15608233"
        assert record[0]["LinkSetDb"]["Link"][92]["Score"]==34588400
        assert record[0]["LinkSetDb"]["Link"][93]["Id"]=="16672057"
        assert record[0]["LinkSetDb"]["Link"][93]["Score"]==34583177
        assert record[0]["LinkSetDb"]["Link"][94]["Id"]=="15687015"
        assert record[0]["LinkSetDb"]["Link"][94]["Score"]==34357840
        assert record[0]["LinkSetDb"]["Link"][95]["Id"]=="10782070"
        assert record[0]["LinkSetDb"]["Link"][95]["Score"]==34326746
        assert record[0]["LinkSetDb"]["Link"][96]["Id"]=="14970722"
        assert record[0]["LinkSetDb"]["Link"][96]["Score"]==34217911
        assert record[0]["LinkSetDb"]["Link"][97]["Id"]=="18027007"
        assert record[0]["LinkSetDb"]["Link"][97]["Score"]==34185436
        assert record[0]["LinkSetDb"]["Link"][98]["Id"]=="12387845"
        assert record[0]["LinkSetDb"]["Link"][98]["Score"]==34083368
        assert record[0]["LinkSetDb"]["Link"][99]["Id"]=="16237012"
        assert record[0]["LinkSetDb"]["Link"][99]["Score"]==34070163
        assert record[0]["LinkSetDb"]["Link"][100]["Id"]=="16351742"
        assert record[0]["LinkSetDb"]["Link"][100]["Score"]==33775198
        assert record[0]["LinkSetDb"]["Link"][101]["Id"]=="12203989"
        assert record[0]["LinkSetDb"]["Link"][101]["Score"]==33759170
        assert record[0]["LinkSetDb"]["Link"][102]["Id"]=="15474306"
        assert record[0]["LinkSetDb"]["Link"][102]["Score"]==33737675
        assert record[0]["LinkSetDb"]["Link"][103]["Id"]=="15270538"
        assert record[0]["LinkSetDb"]["Link"][103]["Score"]==33697306
        assert record[0]["LinkSetDb"]["Link"][104]["Id"]=="17518759"
        assert record[0]["LinkSetDb"]["Link"][104]["Score"]==33695140
        assert record[0]["LinkSetDb"]["Link"][105]["Id"]=="16085497"
        assert record[0]["LinkSetDb"]["Link"][105]["Score"]==33652537
        assert record[0]["LinkSetDb"]["Link"][106]["Id"]=="16423288"
        assert record[0]["LinkSetDb"]["Link"][106]["Score"]==33564554
        assert record[0]["LinkSetDb"]["Link"][107]["Id"]=="16251775"
        assert record[0]["LinkSetDb"]["Link"][107]["Score"]==33547325
        assert record[0]["LinkSetDb"]["Link"][108]["Id"]=="12632152"
        assert record[0]["LinkSetDb"]["Link"][108]["Score"]==33497998
        assert record[0]["LinkSetDb"]["Link"][109]["Id"]=="11269648"
        assert record[0]["LinkSetDb"]["Link"][109]["Score"]==33493800
        assert record[0]["LinkSetDb"]["Link"][110]["Id"]=="16103603"
        assert record[0]["LinkSetDb"]["Link"][110]["Score"]==33378796
        assert record[0]["LinkSetDb"]["Link"][111]["Id"]=="12816546"
        assert record[0]["LinkSetDb"]["Link"][111]["Score"]==33316167
        assert record[0]["LinkSetDb"]["Link"][112]["Id"]=="10221636"
        assert record[0]["LinkSetDb"]["Link"][112]["Score"]==33310814
        assert record[0]["LinkSetDb"]["Link"][113]["Id"]=="16381973"
        assert record[0]["LinkSetDb"]["Link"][113]["Score"]==33236048
        assert record[0]["LinkSetDb"]["Link"][114]["Id"]=="15977173"
        assert record[0]["LinkSetDb"]["Link"][114]["Score"]==33222497
        assert record[0]["LinkSetDb"]["Link"][115]["Id"]=="16351753"
        assert record[0]["LinkSetDb"]["Link"][115]["Score"]==33205084
        assert record[0]["LinkSetDb"]["Link"][116]["Id"]=="15317790"
        assert record[0]["LinkSetDb"]["Link"][116]["Score"]==33195439
        assert record[0]["LinkSetDb"]["Link"][117]["Id"]=="17135198"
        assert record[0]["LinkSetDb"]["Link"][117]["Score"]==33189951
        assert record[0]["LinkSetDb"]["Link"][118]["Id"]=="12701381"
        assert record[0]["LinkSetDb"]["Link"][118]["Score"]==33172200
        assert record[0]["LinkSetDb"]["Link"][119]["Id"]=="12203988"
        assert record[0]["LinkSetDb"]["Link"][119]["Score"]==33172077
        assert record[0]["LinkSetDb"]["Link"][120]["Id"]=="11456466"
        assert record[0]["LinkSetDb"]["Link"][120]["Score"]==33124900
        assert record[0]["LinkSetDb"]["Link"][121]["Id"]=="16936055"
        assert record[0]["LinkSetDb"]["Link"][121]["Score"]==33081742
        assert record[0]["LinkSetDb"]["Link"][122]["Id"]=="17183477"
        assert record[0]["LinkSetDb"]["Link"][122]["Score"]==33005068
        assert record[0]["LinkSetDb"]["Link"][123]["Id"]=="9455480"
        assert record[0]["LinkSetDb"]["Link"][123]["Score"]==32997067
        assert record[0]["LinkSetDb"]["Link"][124]["Id"]=="12490454"
        assert record[0]["LinkSetDb"]["Link"][124]["Score"]==32995041
        assert record[0]["LinkSetDb"]["Link"][125]["Id"]=="12435493"
        assert record[0]["LinkSetDb"]["Link"][125]["Score"]==32990122
        assert record[0]["LinkSetDb"]["Link"][126]["Id"]=="11038309"
        assert record[0]["LinkSetDb"]["Link"][126]["Score"]==32977663
        assert record[0]["LinkSetDb"]["Link"][127]["Id"]=="10366827"
        assert record[0]["LinkSetDb"]["Link"][127]["Score"]==32903347
        assert record[0]["LinkSetDb"]["Link"][128]["Id"]=="10466136"
        assert record[0]["LinkSetDb"]["Link"][128]["Score"]==32869387
        assert record[0]["LinkSetDb"]["Link"][129]["Id"]=="16381840"
        assert record[0]["LinkSetDb"]["Link"][129]["Score"]==32816923
        assert record[0]["LinkSetDb"]["Link"][130]["Id"]=="11825250"
        assert record[0]["LinkSetDb"]["Link"][130]["Score"]==32776183
        assert record[0]["LinkSetDb"]["Link"][131]["Id"]=="12234534"
        assert record[0]["LinkSetDb"]["Link"][131]["Score"]==32708547
        assert record[0]["LinkSetDb"]["Link"][132]["Id"]=="14624247"
        assert record[0]["LinkSetDb"]["Link"][132]["Score"]==32708542
        assert record[0]["LinkSetDb"]["Link"][133]["Id"]=="12886019"
        assert record[0]["LinkSetDb"]["Link"][133]["Score"]==32653276
        assert record[0]["LinkSetDb"]["Link"][134]["Id"]=="12041732"
        assert record[0]["LinkSetDb"]["Link"][134]["Score"]==32607185
        assert record[0]["LinkSetDb"]["Link"][135]["Id"]=="15336912"
        assert record[0]["LinkSetDb"]["Link"][135]["Score"]==32596453
        assert record[0]["LinkSetDb"]["Link"][136]["Id"]=="12652910"
        assert record[0]["LinkSetDb"]["Link"][136]["Score"]==32567397
        assert record[0]["LinkSetDb"]["Link"][137]["Id"]=="14681353"
        assert record[0]["LinkSetDb"]["Link"][137]["Score"]==32549157
        assert record[0]["LinkSetDb"]["Link"][138]["Id"]=="12586873"
        assert record[0]["LinkSetDb"]["Link"][138]["Score"]==32504063
        assert record[0]["LinkSetDb"]["Link"][139]["Id"]=="11481430"
        assert record[0]["LinkSetDb"]["Link"][139]["Score"]==32462602
        assert record[0]["LinkSetDb"]["Link"][140]["Id"]=="15254259"
        assert record[0]["LinkSetDb"]["Link"][140]["Score"]==32441737
        assert record[0]["LinkSetDb"]["Link"][141]["Id"]=="16873516"
        assert record[0]["LinkSetDb"]["Link"][141]["Score"]==32433603
        assert record[0]["LinkSetDb"]["Link"][142]["Id"]=="17170002"
        assert record[0]["LinkSetDb"]["Link"][142]["Score"]==32425626
        assert record[0]["LinkSetDb"]["Link"][143]["Id"]=="12519941"
        assert record[0]["LinkSetDb"]["Link"][143]["Score"]==32367760
        assert record[0]["LinkSetDb"]["Link"][144]["Id"]=="11197770"
        assert record[0]["LinkSetDb"]["Link"][144]["Score"]==32362623
        assert record[0]["LinkSetDb"]["Link"][145]["Id"]=="11240843"
        assert record[0]["LinkSetDb"]["Link"][145]["Score"]==32347064
        assert record[0]["LinkSetDb"]["Link"][146]["Id"]=="11328780"
        assert record[0]["LinkSetDb"]["Link"][146]["Score"]==32333807
        assert record[0]["LinkSetDb"]["Link"][147]["Id"]=="11875041"
        assert record[0]["LinkSetDb"]["Link"][147]["Score"]==32312036
        assert record[0]["LinkSetDb"]["Link"][148]["Id"]=="11752243"
        assert record[0]["LinkSetDb"]["Link"][148]["Score"]==32268199
        assert record[0]["LinkSetDb"]["Link"][149]["Id"]=="16907992"
        assert record[0]["LinkSetDb"]["Link"][149]["Score"]==32247019
        assert record[0]["LinkSetDb"]["Link"][150]["Id"]=="15046636"
        assert record[0]["LinkSetDb"]["Link"][150]["Score"]==32214942
        assert record[0]["LinkSetDb"]["Link"][151]["Id"]=="10592169"
        assert record[0]["LinkSetDb"]["Link"][151]["Score"]==32137798
        assert record[0]["LinkSetDb"]["Link"][152]["Id"]=="17919582"
        assert record[0]["LinkSetDb"]["Link"][152]["Score"]==32137767
        assert record[0]["LinkSetDb"]["Link"][153]["Id"]=="18025705"
        assert record[0]["LinkSetDb"]["Link"][153]["Score"]==32131322
        assert record[0]["LinkSetDb"]["Link"][154]["Id"]=="11029673"
        assert record[0]["LinkSetDb"]["Link"][154]["Score"]==32126363
        assert record[0]["LinkSetDb"]["Link"][155]["Id"]=="9047337"
        assert record[0]["LinkSetDb"]["Link"][155]["Score"]==32090163
        assert record[0]["LinkSetDb"]["Link"][156]["Id"]=="11080372"
        assert record[0]["LinkSetDb"]["Link"][156]["Score"]==31924475
        assert record[0]["LinkSetDb"]["Link"][157]["Id"]=="18045790"
        assert record[0]["LinkSetDb"]["Link"][157]["Score"]==31834367
        assert record[0]["LinkSetDb"]["Link"][158]["Id"]=="10215019"
        assert record[0]["LinkSetDb"]["Link"][158]["Score"]==31823989
        assert record[0]["LinkSetDb"]["Link"][159]["Id"]=="14706096"
        assert record[0]["LinkSetDb"]["Link"][159]["Score"]==31781977
        assert record[0]["LinkSetDb"]["Link"][160]["Id"]=="17537593"
        assert record[0]["LinkSetDb"]["Link"][160]["Score"]==31771566
        assert record[0]["LinkSetDb"]["Link"][161]["Id"]=="12819149"
        assert record[0]["LinkSetDb"]["Link"][161]["Score"]==31683943
        assert record[0]["LinkSetDb"]["Link"][162]["Id"]=="17880721"
        assert record[0]["LinkSetDb"]["Link"][162]["Score"]==31630816
        assert record[0]["LinkSetDb"]["Link"][163]["Id"]=="14681478"
        assert record[0]["LinkSetDb"]["Link"][163]["Score"]==31620457
        assert record[0]["LinkSetDb"]["Link"][164]["Id"]=="11985867"
        assert record[0]["LinkSetDb"]["Link"][164]["Score"]==31544318
        assert record[0]["LinkSetDb"]["Link"][165]["Id"]=="15608248"
        assert record[0]["LinkSetDb"]["Link"][165]["Score"]==31542256
        assert record[0]["LinkSetDb"]["Link"][166]["Id"]=="17401150"
        assert record[0]["LinkSetDb"]["Link"][166]["Score"]==31497289
        assert record[0]["LinkSetDb"]["Link"][167]["Id"]=="10359795"
        assert record[0]["LinkSetDb"]["Link"][167]["Score"]==31460779
        assert record[0]["LinkSetDb"]["Link"][168]["Id"]=="15608286"
        assert record[0]["LinkSetDb"]["Link"][168]["Score"]==31435112
        assert record[0]["LinkSetDb"]["Link"][169]["Id"]=="15774022"
        assert record[0]["LinkSetDb"]["Link"][169]["Score"]==31425851
        assert record[0]["LinkSetDb"]["Link"][170]["Id"]=="9921679"
        assert record[0]["LinkSetDb"]["Link"][170]["Score"]==31396086
        assert record[0]["LinkSetDb"]["Link"][171]["Id"]=="17038195"
        assert record[0]["LinkSetDb"]["Link"][171]["Score"]==31380822
        assert record[0]["LinkSetDb"]["Link"][172]["Id"]=="15491544"
        assert record[0]["LinkSetDb"]["Link"][172]["Score"]==31294370
        assert record[0]["LinkSetDb"]["Link"][173]["Id"]=="10469257"
        assert record[0]["LinkSetDb"]["Link"][173]["Score"]==31291548
        assert record[0]["LinkSetDb"]["Link"][174]["Id"]=="15487498"
        assert record[0]["LinkSetDb"]["Link"][174]["Score"]==31268351
        assert record[0]["LinkSetDb"]["Link"][175]["Id"]=="15383303"
        assert record[0]["LinkSetDb"]["Link"][175]["Score"]==31264596
        assert record[0]["LinkSetDb"]["Link"][176]["Id"]=="15643605"
        assert record[0]["LinkSetDb"]["Link"][176]["Score"]==31259953
        assert record[0]["LinkSetDb"]["Link"][177]["Id"]=="16418238"
        assert record[0]["LinkSetDb"]["Link"][177]["Score"]==31259003
        assert record[0]["LinkSetDb"]["Link"][178]["Id"]=="15500248"
        assert record[0]["LinkSetDb"]["Link"][178]["Score"]==31252080
        assert record[0]["LinkSetDb"]["Link"][179]["Id"]=="15479945"
        assert record[0]["LinkSetDb"]["Link"][179]["Score"]==31249988
        assert record[0]["LinkSetDb"]["Link"][180]["Id"]=="16962738"
        assert record[0]["LinkSetDb"]["Link"][180]["Score"]==31249405
        assert record[0]["LinkSetDb"]["Link"][181]["Id"]=="15094394"
        assert record[0]["LinkSetDb"]["Link"][181]["Score"]==31200337
        assert record[0]["LinkSetDb"]["Link"][182]["Id"]=="11758285"
        assert record[0]["LinkSetDb"]["Link"][182]["Score"]==31180435
        assert record[0]["LinkSetDb"]["Link"][183]["Id"]=="15723693"
        assert record[0]["LinkSetDb"]["Link"][183]["Score"]==31083464
        assert record[0]["LinkSetDb"]["Link"][184]["Id"]=="16710453"
        assert record[0]["LinkSetDb"]["Link"][184]["Score"]==31083136
        assert record[0]["LinkSetDb"]["Link"][185]["Id"]=="15311460"
        assert record[0]["LinkSetDb"]["Link"][185]["Score"]==31068402
        assert record[0]["LinkSetDb"]["Link"][186]["Id"]=="16549670"
        assert record[0]["LinkSetDb"]["Link"][186]["Score"]==30995148
        assert record[0]["LinkSetDb"]["Link"][187]["Id"]=="18180957"
        assert record[0]["LinkSetDb"]["Link"][187]["Score"]==30973190
        assert record[0]["LinkSetDb"]["Link"][188]["Id"]=="14681351"
        assert record[0]["LinkSetDb"]["Link"][188]["Score"]==30968930
        assert record[0]["LinkSetDb"]["Link"][189]["Id"]=="10902212"
        assert record[0]["LinkSetDb"]["Link"][189]["Score"]==30960861
        assert record[0]["LinkSetDb"]["Link"][190]["Id"]=="15357877"
        assert record[0]["LinkSetDb"]["Link"][190]["Score"]==30947680
        assert record[0]["LinkSetDb"]["Link"][191]["Id"]=="12356773"
        assert record[0]["LinkSetDb"]["Link"][191]["Score"]==30910321
        assert record[0]["LinkSetDb"]["Link"][192]["Id"]=="17537669"
        assert record[0]["LinkSetDb"]["Link"][192]["Score"]==30893205
        assert record[0]["LinkSetDb"]["Link"][193]["Id"]=="16551372"
        assert record[0]["LinkSetDb"]["Link"][193]["Score"]==30889080
        assert record[0]["LinkSetDb"]["Link"][194]["Id"]=="15231810"
        assert record[0]["LinkSetDb"]["Link"][194]["Score"]==30863616
        assert record[0]["LinkSetDb"]["Link"][195]["Id"]=="12819150"
        assert record[0]["LinkSetDb"]["Link"][195]["Score"]==30847027
        assert record[0]["LinkSetDb"]["Link"][196]["Id"]=="15608257"
        assert record[0]["LinkSetDb"]["Link"][196]["Score"]==30840234
        assert record[0]["LinkSetDb"]["Link"][197]["Id"]=="17384426"
        assert record[0]["LinkSetDb"]["Link"][197]["Score"]==30827754
        assert record[0]["LinkSetDb"]["Link"][198]["Id"]=="15811532"
        assert record[0]["LinkSetDb"]["Link"][198]["Score"]==30823185
        assert record[0]["LinkSetDb"]["Link"][199]["Id"]=="10612821"
        assert record[0]["LinkSetDb"]["Link"][199]["Score"]==30822187
        assert record[0]["LinkSetDb"]["Link"][200]["Id"]=="17062145"
        assert record[0]["LinkSetDb"]["Link"][200]["Score"]==30813605
        assert record[0]["LinkSetDb"]["Link"][201]["Id"]=="11355885"
        assert record[0]["LinkSetDb"]["Link"][201]["Score"]==30810648
        assert record[0]["LinkSetDb"]["Link"][202]["Id"]=="15746365"
        assert record[0]["LinkSetDb"]["Link"][202]["Score"]==30784209
        assert record[0]["LinkSetDb"]["Link"][203]["Id"]=="16282300"
        assert record[0]["LinkSetDb"]["Link"][203]["Score"]==30782807
        assert record[0]["LinkSetDb"]["Link"][204]["Id"]=="15546336"
        assert record[0]["LinkSetDb"]["Link"][204]["Score"]==30773578
        assert record[0]["LinkSetDb"]["Link"][205]["Id"]=="11741630"
        assert record[0]["LinkSetDb"]["Link"][205]["Score"]==30764995
        assert record[0]["LinkSetDb"]["Link"][206]["Id"]=="15980532"
        assert record[0]["LinkSetDb"]["Link"][206]["Score"]==30735790
        assert record[0]["LinkSetDb"]["Link"][207]["Id"]=="12519977"
        assert record[0]["LinkSetDb"]["Link"][207]["Score"]==30707395
        assert record[0]["LinkSetDb"]["Link"][208]["Id"]=="12436197"
        assert record[0]["LinkSetDb"]["Link"][208]["Score"]==30705501
        assert record[0]["LinkSetDb"]["Link"][209]["Id"]=="11125059"
        assert record[0]["LinkSetDb"]["Link"][209]["Score"]==30614888
        assert record[0]["LinkSetDb"]["Link"][210]["Id"]=="11163442"
        assert record[0]["LinkSetDb"]["Link"][210]["Score"]==30550965
        assert record[0]["LinkSetDb"]["Link"][211]["Id"]=="12519964"
        assert record[0]["LinkSetDb"]["Link"][211]["Score"]==30518025
        assert record[0]["LinkSetDb"]["Link"][212]["Id"]=="12083398"
        assert record[0]["LinkSetDb"]["Link"][212]["Score"]==30466595
        assert record[0]["LinkSetDb"]["Link"][213]["Id"]=="11908756"
        assert record[0]["LinkSetDb"]["Link"][213]["Score"]==30462080
        assert record[0]["LinkSetDb"]["Link"][214]["Id"]=="15608226"
        assert record[0]["LinkSetDb"]["Link"][214]["Score"]==30335152
        assert record[0]["LinkSetDb"]["Link"][215]["Id"]=="16845091"
        assert record[0]["LinkSetDb"]["Link"][215]["Score"]==30277120
        assert record[0]["LinkSetDb"]["Link"][216]["Id"]=="17338820"
        assert record[0]["LinkSetDb"]["Link"][216]["Score"]==30208452
        assert record[0]["LinkSetDb"]["Link"][217]["Id"]=="10407783"
        assert record[0]["LinkSetDb"]["Link"][217]["Score"]==30171504
        assert record[0]["LinkSetDb"]["Link"][218]["Id"]=="17130148"
        assert record[0]["LinkSetDb"]["Link"][218]["Score"]==30160136
        assert record[0]["LinkSetDb"]["Link"][219]["Id"]=="14681471"
        assert record[0]["LinkSetDb"]["Link"][219]["Score"]==30155757
        assert record[0]["LinkSetDb"]["Link"][220]["Id"]=="17445272"
        assert record[0]["LinkSetDb"]["Link"][220]["Score"]==30015229
        assert record[0]["LinkSetDb"]["Link"][221]["Id"]=="11279516"
        assert record[0]["LinkSetDb"]["Link"][221]["Score"]==29947199
        assert record[0]["LinkSetDb"]["Link"][222]["Id"]=="17221864"
        assert record[0]["LinkSetDb"]["Link"][222]["Score"]==29893674
        assert record[0]["LinkSetDb"]["Link"][223]["Id"]=="15827081"
        assert record[0]["LinkSetDb"]["Link"][223]["Score"]==29891924
        assert record[0]["LinkSetDb"]["Link"][224]["Id"]=="11222582"
        assert record[0]["LinkSetDb"]["Link"][224]["Score"]==29878915
        assert record[0]["LinkSetDb"]["Link"][225]["Id"]=="11384164"
        assert record[0]["LinkSetDb"]["Link"][225]["Score"]==29871698
        assert record[0]["LinkSetDb"]["Link"][226]["Id"]=="17877839"
        assert record[0]["LinkSetDb"]["Link"][226]["Score"]==29843765
        assert record[0]["LinkSetDb"]["Link"][227]["Id"]=="17151077"
        assert record[0]["LinkSetDb"]["Link"][227]["Score"]==29841695
        assert record[0]["LinkSetDb"]["Link"][228]["Id"]=="16381974"
        assert record[0]["LinkSetDb"]["Link"][228]["Score"]==29740312
        assert record[0]["LinkSetDb"]["Link"][229]["Id"]=="10592263"
        assert record[0]["LinkSetDb"]["Link"][229]["Score"]==29633946
        assert record[0]["LinkSetDb"]["Link"][230]["Id"]=="15608212"
        assert record[0]["LinkSetDb"]["Link"][230]["Score"]==29621479
        assert record[0]["LinkSetDb"]["Link"][231]["Id"]=="9847217"
        assert record[0]["LinkSetDb"]["Link"][231]["Score"]==29618439
        assert record[0]["LinkSetDb"]["Link"][232]["Id"]=="17142236"
        assert record[0]["LinkSetDb"]["Link"][232]["Score"]==29577611
        assert record[0]["LinkSetDb"]["Link"][233]["Id"]=="17059604"
        assert record[0]["LinkSetDb"]["Link"][233]["Score"]==29569767
        assert record[0]["LinkSetDb"]["Link"][234]["Id"]=="16845079"
        assert record[0]["LinkSetDb"]["Link"][234]["Score"]==29506663
        assert record[0]["LinkSetDb"]["Link"][235]["Id"]=="14727153"
        assert record[0]["LinkSetDb"]["Link"][235]["Score"]==29368276
        assert record[0]["LinkSetDb"]["Link"][236]["Id"]=="18045498"
        assert record[0]["LinkSetDb"]["Link"][236]["Score"]==29364312
        assert record[0]["LinkSetDb"]["Link"][237]["Id"]=="17185755"
        assert record[0]["LinkSetDb"]["Link"][237]["Score"]==29331905
        assert record[0]["LinkSetDb"]["Link"][238]["Id"]=="18025704"
        assert record[0]["LinkSetDb"]["Link"][238]["Score"]==29323161
        assert record[0]["LinkSetDb"]["Link"][239]["Id"]=="15215374"
        assert record[0]["LinkSetDb"]["Link"][239]["Score"]==29306559
        assert record[0]["LinkSetDb"]["Link"][240]["Id"]=="17135185"
        assert record[0]["LinkSetDb"]["Link"][240]["Score"]==29236297
        assert record[0]["LinkSetDb"]["Link"][241]["Id"]=="10466135"
        assert record[0]["LinkSetDb"]["Link"][241]["Score"]==29231855
        assert record[0]["LinkSetDb"]["Link"][242]["Id"]=="17148475"
        assert record[0]["LinkSetDb"]["Link"][242]["Score"]==29229044
        assert record[0]["LinkSetDb"]["Link"][243]["Id"]=="15657101"
        assert record[0]["LinkSetDb"]["Link"][243]["Score"]==29209567
        assert record[0]["LinkSetDb"]["Link"][244]["Id"]=="14681490"
        assert record[0]["LinkSetDb"]["Link"][244]["Score"]==29189708
        assert record[0]["LinkSetDb"]["Link"][245]["Id"]=="15714328"
        assert record[0]["LinkSetDb"]["Link"][245]["Score"]==29183488
        assert record[0]["LinkSetDb"]["Link"][246]["Id"]=="14960477"
        assert record[0]["LinkSetDb"]["Link"][246]["Score"]==29040531
        assert record[0]["LinkSetDb"]["Link"][247]["Id"]=="11015564"
        assert record[0]["LinkSetDb"]["Link"][247]["Score"]==29011368
        assert record[0]["LinkSetDb"]["Link"][248]["Id"]=="18064491"
        assert record[0]["LinkSetDb"]["Link"][248]["Score"]==28956740
        assert record[0]["LinkSetDb"]["Link"][249]["Id"]=="12734009"
        assert record[0]["LinkSetDb"]["Link"][249]["Score"]==28950064
        assert record[0]["LinkSetDb"]["Link"][250]["Id"]=="17094804"
        assert record[0]["LinkSetDb"]["Link"][250]["Score"]==28906953
        assert record[0]["LinkSetDb"]["Link"][251]["Id"]=="17908294"
        assert record[0]["LinkSetDb"]["Link"][251]["Score"]==28897717
        assert record[0]["LinkSetDb"]["Link"][252]["Id"]=="16176584"
        assert record[0]["LinkSetDb"]["Link"][252]["Score"]==28874470
        assert record[0]["LinkSetDb"]["Link"][253]["Id"]=="14715089"
        assert record[0]["LinkSetDb"]["Link"][253]["Score"]==28763886
        assert record[0]["LinkSetDb"]["Link"][254]["Id"]=="14681408"
        assert record[0]["LinkSetDb"]["Link"][254]["Score"]==28697827
        assert record[0]["LinkSetDb"]["Link"][255]["Id"]=="14594716"
        assert record[0]["LinkSetDb"]["Link"][255]["Score"]==28686075
        assert record[0]["LinkSetDb"]["Link"][256]["Id"]=="16528802"
        assert record[0]["LinkSetDb"]["Link"][256]["Score"]==28644452
        assert record[0]["LinkSetDb"]["Link"][257]["Id"]=="16010002"
        assert record[0]["LinkSetDb"]["Link"][257]["Score"]==28637570
        assert record[0]["LinkSetDb"]["Link"][258]["Id"]=="17430565"
        assert record[0]["LinkSetDb"]["Link"][258]["Score"]==28635513
        assert record[0]["LinkSetDb"]["Link"][259]["Id"]=="16452787"
        assert record[0]["LinkSetDb"]["Link"][259]["Score"]==28631832
        assert record[0]["LinkSetDb"]["Link"][260]["Id"]=="11197127"
        assert record[0]["LinkSetDb"]["Link"][260]["Score"]==28619225
        assert record[0]["LinkSetDb"]["Link"][261]["Id"]=="8682188"
        assert record[0]["LinkSetDb"]["Link"][261]["Score"]==28592521
        assert record[0]["LinkSetDb"]["Link"][262]["Id"]=="12519940"
        assert record[0]["LinkSetDb"]["Link"][262]["Score"]==28573991
        assert record[0]["LinkSetDb"]["Link"][263]["Id"]=="17121775"
        assert record[0]["LinkSetDb"]["Link"][263]["Score"]==28448726
        assert record[0]["LinkSetDb"]["Link"][264]["Id"]=="16371163"
        assert record[0]["LinkSetDb"]["Link"][264]["Score"]==28373394
        assert record[0]["LinkSetDb"]["Link"][265]["Id"]=="15300845"
        assert record[0]["LinkSetDb"]["Link"][265]["Score"]==28338477
        assert record[0]["LinkSetDb"]["Link"][266]["Id"]=="15248903"
        assert record[0]["LinkSetDb"]["Link"][266]["Score"]==28323328
        assert record[0]["LinkSetDb"]["Link"][267]["Id"]=="11319266"
        assert record[0]["LinkSetDb"]["Link"][267]["Score"]==28293166
        assert record[0]["LinkSetDb"]["Link"][268]["Id"]=="16336665"
        assert record[0]["LinkSetDb"]["Link"][268]["Score"]==28231249
        assert record[0]["LinkSetDb"]["Link"][269]["Id"]=="14681350"
        assert record[0]["LinkSetDb"]["Link"][269]["Score"]==28227327
        assert record[0]["LinkSetDb"]["Link"][270]["Id"]=="16216831"
        assert record[0]["LinkSetDb"]["Link"][270]["Score"]==28224610
        assert record[0]["LinkSetDb"]["Link"][271]["Id"]=="15494741"
        assert record[0]["LinkSetDb"]["Link"][271]["Score"]==28190925
        assert record[0]["LinkSetDb"]["Link"][272]["Id"]=="17088289"
        assert record[0]["LinkSetDb"]["Link"][272]["Score"]==28168901
        assert record[0]["LinkSetDb"]["Link"][273]["Id"]=="17099235"
        assert record[0]["LinkSetDb"]["Link"][273]["Score"]==28159766
        assert record[0]["LinkSetDb"]["Link"][274]["Id"]=="15215372"
        assert record[0]["LinkSetDb"]["Link"][274]["Score"]==28129693
        assert record[0]["LinkSetDb"]["Link"][275]["Id"]=="9169870"
        assert record[0]["LinkSetDb"]["Link"][275]["Score"]==28117392
        assert record[0]["LinkSetDb"]["Link"][276]["Id"]=="10077537"
        assert record[0]["LinkSetDb"]["Link"][276]["Score"]==27911205
        assert record[0]["LinkSetDb"]["Link"][277]["Id"]=="18172929"
        assert record[0]["LinkSetDb"]["Link"][277]["Score"]==27885172
        assert record[0]["LinkSetDb"]["Link"][278]["Id"]=="9571806"
        assert record[0]["LinkSetDb"]["Link"][278]["Score"]==27841468
        assert record[0]["LinkSetDb"]["Link"][279]["Id"]=="11752280"
        assert record[0]["LinkSetDb"]["Link"][279]["Score"]==27795833
        assert record[0]["LinkSetDb"]["Link"][280]["Id"]=="11414208"
        assert record[0]["LinkSetDb"]["Link"][280]["Score"]==27725996
        assert record[0]["LinkSetDb"]["Link"][281]["Id"]=="9298642"
        assert record[0]["LinkSetDb"]["Link"][281]["Score"]==27716027
        assert record[0]["LinkSetDb"]["Link"][282]["Id"]=="18073380"
        assert record[0]["LinkSetDb"]["Link"][282]["Score"]==27437383
        assert record[0]["LinkSetDb"]["Link"][283]["Id"]=="14527308"
        assert record[0]["LinkSetDb"]["Link"][283]["Score"]==27332641
        assert record[0]["LinkSetDb"]["Link"][284]["Id"]=="9847220"
        assert record[0]["LinkSetDb"]["Link"][284]["Score"]==27083894
        assert record[0]["LinkSetDb"]["Link"][285]["Id"]=="10413661"
        assert record[0]["LinkSetDb"]["Link"][285]["Score"]==27073030
        assert record[0]["LinkSetDb"]["Link"][286]["Id"]=="10407677"
        assert record[0]["LinkSetDb"]["Link"][286]["Score"]==26907635
        assert record[0]["LinkSetDb"]["Link"][287]["Id"]=="11244060"
        assert record[0]["LinkSetDb"]["Link"][287]["Score"]==26897688
        assert record[0]["LinkSetDb"]["Link"][288]["Id"]=="10227170"
        assert record[0]["LinkSetDb"]["Link"][288]["Score"]==26766431
        assert record[0]["LinkSetDb"]["Link"][289]["Id"]=="8719164"
        assert record[0]["LinkSetDb"]["Link"][289]["Score"]==26515360
        assert record[0]["LinkSetDb"]["Link"][290]["Id"]=="18359019"
        assert record[0]["LinkSetDb"]["Link"][290]["Score"]==26225983
        assert record[0]["LinkSetDb"]["Link"][291]["Id"]=="10511680"
        assert record[0]["LinkSetDb"]["Link"][291]["Score"]==26031196
        assert record[0]["LinkSetDb"]["Link"][292]["Id"]=="9884329"
        assert record[0]["LinkSetDb"]["Link"][292]["Score"]==25992564
        assert record[0]["LinkSetDb"]["Link"][293]["Id"]=="17827295"
        assert record[0]["LinkSetDb"]["Link"][293]["Score"]==25989152
        assert record[0]["LinkSetDb"]["Link"][294]["Id"]=="10899154"
        assert record[0]["LinkSetDb"]["Link"][294]["Score"]==25843128
        assert record[0]["LinkSetDb"]["Link"][295]["Id"]=="11668619"
        assert record[0]["LinkSetDb"]["Link"][295]["Score"]==25822950
        assert record[0]["LinkSetDb"]["Link"][296]["Id"]=="18386064"
        assert record[0]["LinkSetDb"]["Link"][296]["Score"]==25702942
        assert record[0]["LinkSetDb"]["Link"][297]["Id"]=="11092731"
        assert record[0]["LinkSetDb"]["Link"][297]["Score"]==25618899
        assert record[0]["LinkSetDb"]["Link"][298]["Id"]=="9520376"
        assert record[0]["LinkSetDb"]["Link"][298]["Score"]==25549761
        assert record[0]["LinkSetDb"]["Link"][299]["Id"]=="11756688"
        assert record[0]["LinkSetDb"]["Link"][299]["Score"]==25440634
        assert record[0]["LinkSetDb"]["Link"][300]["Id"]=="10737802"
        assert record[0]["LinkSetDb"]["Link"][300]["Score"]==25362744
        assert record[0]["LinkSetDb"]["Link"][301]["Id"]=="9879937"
        assert record[0]["LinkSetDb"]["Link"][301]["Score"]==25277089
        assert record[0]["LinkSetDb"]["Link"][302]["Id"]=="17822801"
        assert record[0]["LinkSetDb"]["Link"][302]["Score"]==25252984
        assert record[0]["LinkSetDb"]["Link"][303]["Id"]=="10965872"
        assert record[0]["LinkSetDb"]["Link"][303]["Score"]==25208185
        assert record[0]["LinkSetDb"]["Link"][304]["Id"]=="10511682"
        assert record[0]["LinkSetDb"]["Link"][304]["Score"]==25183443
        assert record[0]["LinkSetDb"]["Link"][305]["Id"]=="10851186"
        assert record[0]["LinkSetDb"]["Link"][305]["Score"]==25092764
        assert record[0]["LinkSetDb"]["Link"][306]["Id"]=="9775388"
        assert record[0]["LinkSetDb"]["Link"][306]["Score"]==25026910
        assert record[0]["LinkSetDb"]["Link"][307]["Id"]=="10810023"
        assert record[0]["LinkSetDb"]["Link"][307]["Score"]==24904718
        assert record[0]["LinkSetDb"]["Link"][308]["Id"]=="18032438"
        assert record[0]["LinkSetDb"]["Link"][308]["Score"]==24509777
        assert record[0]["LinkSetDb"]["Link"][309]["Id"]=="18377816"
        assert record[0]["LinkSetDb"]["Link"][309]["Score"]==24373788
        assert record[0]["LinkSetDb"]["Link"][310]["Id"]=="11774190"
        assert record[0]["LinkSetDb"]["Link"][310]["Score"]==24185658
        assert record[0]["LinkSetDb"]["Link"][311]["Id"]=="10484179"
        assert record[0]["LinkSetDb"]["Link"][311]["Score"]==24122767
        assert record[0]["LinkSetDb"]["Link"][312]["Id"]=="9625791"
        assert record[0]["LinkSetDb"]["Link"][312]["Score"]==24049917
        assert record[0]["LinkSetDb"]["Link"][313]["Id"]=="11446511"
        assert record[0]["LinkSetDb"]["Link"][313]["Score"]==24048253
        assert record[0]["LinkSetDb"]["Link"][314]["Id"]=="10066467"
        assert record[0]["LinkSetDb"]["Link"][314]["Score"]==23968405
        assert record[0]["LinkSetDb"]["Link"][315]["Id"]=="11783003"
        assert record[0]["LinkSetDb"]["Link"][315]["Score"]==23393870
        assert record[0]["LinkSetDb"]["Link"][316]["Id"]=="10611059"
        assert record[0]["LinkSetDb"]["Link"][316]["Score"]==23255298
        assert record[0]["LinkSetDb"]["Link"][317]["Id"]=="10587943"
        assert record[0]["LinkSetDb"]["Link"][317]["Score"]==23014503
        assert record[0]["LinkSetDb"]["Link"][318]["Id"]=="10612820"
        assert record[0]["LinkSetDb"]["Link"][318]["Score"]==22990878
        assert record[0]["LinkSetDb"]["Link"][319]["Id"]=="9685316"
        assert record[0]["LinkSetDb"]["Link"][319]["Score"]==22771348
        assert record[0]["LinkSetDb"]["Link"][320]["Id"]=="11125121"
        assert record[0]["LinkSetDb"]["Link"][320]["Score"]==22732820
        assert record[0]["LinkSetDb"]["Link"][321]["Id"]=="10075567"
        assert record[0]["LinkSetDb"]["Link"][321]["Score"]==22670427
        assert record[0]["LinkSetDb"]["Link"][322]["Id"]=="11084929"
        assert record[0]["LinkSetDb"]["Link"][322]["Score"]==22397665
        assert record[0]["LinkSetDb"]["Link"][323]["Id"]=="11357826"
        assert record[0]["LinkSetDb"]["Link"][323]["Score"]==22362882
        assert record[0]["LinkSetDb"]["Link"][324]["Id"]=="17983575"
        assert record[0]["LinkSetDb"]["Link"][324]["Score"]==22305320
        assert record[0]["LinkSetDb"]["Link"][325]["Id"]=="11038308"
        assert record[0]["LinkSetDb"]["Link"][325]["Score"]==22115670
        assert record[0]["LinkSetDb"]["Link"][326]["Id"]=="18257289"
        assert record[0]["LinkSetDb"]["Link"][326]["Score"]==22053176
        assert record[0]["LinkSetDb"]["Link"][327]["Id"]=="10419978"
        assert record[0]["LinkSetDb"]["Link"][327]["Score"]==22016184
        assert record[0]["LinkSetDb"]["Link"][328]["Id"]=="9421619"
        assert record[0]["LinkSetDb"]["Link"][328]["Score"]==21957407
        assert record[0]["LinkSetDb"]["Link"][329]["Id"]=="10592198"
        assert record[0]["LinkSetDb"]["Link"][329]["Score"]==21803908
        assert record[0]["LinkSetDb"]["Link"][330]["Id"]=="11483982"
        assert record[0]["LinkSetDb"]["Link"][330]["Score"]==20783817
        assert record[0]["LinkSetDb"]["Link"][331]["Id"]=="11329386"
        assert record[0]["LinkSetDb"]["Link"][331]["Score"]==20223493
        assert record[0]["LinkSetDb"]["Link"][332]["Id"]=="10587942"
        assert record[0]["LinkSetDb"]["Link"][332]["Score"]==20208799
        assert record[0]["LinkSetDb"]["Link"][333]["Id"]=="10810024"
        assert record[0]["LinkSetDb"]["Link"][333]["Score"]==19989188
        assert record[0]["LinkSetDb"]["Link"][334]["Id"]=="11480780"
        assert record[0]["LinkSetDb"]["Link"][334]["Score"]==19974101
        assert record[0]["LinkSetDb"]["Link"][335]["Id"]=="11802378"
        assert record[0]["LinkSetDb"]["Link"][335]["Score"]==19738532
        assert record[0]["LinkSetDb"]["Link"][336]["Id"]=="10610803"
        assert record[0]["LinkSetDb"]["Link"][336]["Score"]==19359100
        assert record[0]["LinkSetDb"]["Link"][337]["Id"]=="10407668"
        assert record[0]["LinkSetDb"]["Link"][337]["Score"]==19070525
        assert record[0]["LinkSetDb"]["Link"][338]["Id"]=="18287701"
        assert record[0]["LinkSetDb"]["Link"][338]["Score"]==19065945
        assert record[0]["LinkSetDb"]["Link"][339]["Id"]=="10963611"
        assert record[0]["LinkSetDb"]["Link"][339]["Score"]==18962273
        assert record[0]["LinkSetDb"]["Link"][340]["Id"]=="10447503"
        assert record[0]["LinkSetDb"]["Link"][340]["Score"]==17406980
        assert record[0]["LinkSetDb"]["Link"][341]["Id"]=="9830540"
        assert record[0]["LinkSetDb"]["Link"][341]["Score"]==17143709
        assert record[0]["LinkSetDb"]["Link"][342]["Id"]=="11462837"
        assert record[0]["LinkSetDb"]["Link"][342]["Score"]==16819799
        assert record[0]["LinkSetDb"]["Link"][343]["Id"]=="10637631"
        assert record[0]["LinkSetDb"]["Link"][343]["Score"]==16390796
        assert record[0]["LinkSetDb"]["Link"][344]["Id"]=="11387032"
        assert record[0]["LinkSetDb"]["Link"][344]["Score"]==15698695
        assert record[0]["LinkSetDb"]["Link"][345]["Id"]=="18365535"
        assert record[0]["LinkSetDb"]["Link"][345]["Score"]==15494816
        assert record[0]["LinkSetDb"]["Link"][346]["Id"]=="15181901"
        assert record[0]["LinkSetDb"]["Link"][346]["Score"]==14385628

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
        assert record[0]["LinkSetDb"]["DbTo"]=="pubmed"
        assert record[0]["LinkSetDb"]["LinkName"]=="pubmed_pubmed"
        assert record[0]["LinkSetDb"]["Link"][0]["Id"]=="12242737"
        assert record[0]["LinkSetDb"]["Link"][0]["Score"]==2147483647
        assert record[0]["LinkSetDb"]["Link"][1]["Id"]=="11218011"
        assert record[0]["LinkSetDb"]["Link"][1]["Score"]==50825961
        assert record[0]["LinkSetDb"]["Link"][2]["Id"]=="11329656"
        assert record[0]["LinkSetDb"]["Link"][2]["Score"]==49822043
        assert record[0]["LinkSetDb"]["Link"][3]["Id"]=="9757294"
        assert record[0]["LinkSetDb"]["Link"][3]["Score"]==42645380
        assert record[0]["LinkSetDb"]["Link"][4]["Id"]=="9456947"
        assert record[0]["LinkSetDb"]["Link"][4]["Score"]==39871666
        assert record[0]["LinkSetDb"]["Link"][5]["Id"]=="17193860"
        assert record[0]["LinkSetDb"]["Link"][5]["Score"]==39717388
        assert record[0]["LinkSetDb"]["Link"][6]["Id"]=="11274884"
        assert record[0]["LinkSetDb"]["Link"][6]["Score"]==39233276
        assert record[0]["LinkSetDb"]["Link"][7]["Id"]=="12878072"
        assert record[0]["LinkSetDb"]["Link"][7]["Score"]==37748327
        assert record[0]["LinkSetDb"]["Link"][8]["Id"]=="11125632"
        assert record[0]["LinkSetDb"]["Link"][8]["Score"]==36227857
        assert record[0]["LinkSetDb"]["Link"][9]["Id"]=="12822521"
        assert record[0]["LinkSetDb"]["Link"][9]["Score"]==36170366
        assert record[0]["LinkSetDb"]["Link"][10]["Id"]=="16999328"
        assert record[0]["LinkSetDb"]["Link"][10]["Score"]==36107139
        assert record[0]["LinkSetDb"]["Link"][11]["Id"]=="17875142"
        assert record[0]["LinkSetDb"]["Link"][11]["Score"]==35736802
        assert record[0]["LinkSetDb"]["Link"][12]["Id"]=="9510579"
        assert record[0]["LinkSetDb"]["Link"][12]["Score"]==35206779
        assert record[0]["LinkSetDb"]["Link"][13]["Id"]=="17354190"
        assert record[0]["LinkSetDb"]["Link"][13]["Score"]==34792954
        assert record[0]["LinkSetDb"]["Link"][14]["Id"]=="11702119"
        assert record[0]["LinkSetDb"]["Link"][14]["Score"]==34618984
        assert record[0]["LinkSetDb"]["Link"][15]["Id"]=="10024396"
        assert record[0]["LinkSetDb"]["Link"][15]["Score"]==33877753
        assert record[0]["LinkSetDb"]["Link"][16]["Id"]=="14650118"
        assert record[0]["LinkSetDb"]["Link"][16]["Score"]==33746160
        assert record[0]["LinkSetDb"]["Link"][17]["Id"]=="17243036"
        assert record[0]["LinkSetDb"]["Link"][17]["Score"]==33198930
        assert record[0]["LinkSetDb"]["Link"][18]["Id"]=="16580806"
        assert record[0]["LinkSetDb"]["Link"][18]["Score"]==33117197
        assert record[0]["LinkSetDb"]["Link"][19]["Id"]=="15278705"
        assert record[0]["LinkSetDb"]["Link"][19]["Score"]==33002826
        assert record[0]["LinkSetDb"]["Link"][20]["Id"]=="15236131"
        assert record[0]["LinkSetDb"]["Link"][20]["Score"]==32808406
        assert record[0]["LinkSetDb"]["Link"][21]["Id"]=="11368937"
        assert record[0]["LinkSetDb"]["Link"][21]["Score"]==32277701
        assert record[0]["LinkSetDb"]["Link"][22]["Id"]=="10688065"
        assert record[0]["LinkSetDb"]["Link"][22]["Score"]==32052850
        assert record[0]["LinkSetDb"]["Link"][23]["Id"]=="15635471"
        assert record[0]["LinkSetDb"]["Link"][23]["Score"]==31938251
        assert record[0]["LinkSetDb"]["Link"][24]["Id"]=="16357381"
        assert record[0]["LinkSetDb"]["Link"][24]["Score"]==31780147
        assert record[0]["LinkSetDb"]["Link"][25]["Id"]=="8153333"
        assert record[0]["LinkSetDb"]["Link"][25]["Score"]==31542202
        assert record[0]["LinkSetDb"]["Link"][26]["Id"]=="16284132"
        assert record[0]["LinkSetDb"]["Link"][26]["Score"]==31290577
        assert record[0]["LinkSetDb"]["Link"][27]["Id"]=="11329162"
        assert record[0]["LinkSetDb"]["Link"][27]["Score"]==31163088
        assert record[0]["LinkSetDb"]["Link"][28]["Id"]=="11973040"
        assert record[0]["LinkSetDb"]["Link"][28]["Score"]==31156707
        assert record[0]["LinkSetDb"]["Link"][29]["Id"]=="15143223"
        assert record[0]["LinkSetDb"]["Link"][29]["Score"]==31025329
        assert record[0]["LinkSetDb"]["Link"][30]["Id"]=="17040637"
        assert record[0]["LinkSetDb"]["Link"][30]["Score"]==30990506
        assert record[0]["LinkSetDb"]["Link"][31]["Id"]=="11016058"
        assert record[0]["LinkSetDb"]["Link"][31]["Score"]==30966482
        assert record[0]["LinkSetDb"]["Link"][32]["Id"]=="9317094"
        assert record[0]["LinkSetDb"]["Link"][32]["Score"]==30935529
        assert record[0]["LinkSetDb"]["Link"][33]["Id"]=="16133609"
        assert record[0]["LinkSetDb"]["Link"][33]["Score"]==30580027
        assert record[0]["LinkSetDb"]["Link"][34]["Id"]=="17325998"
        assert record[0]["LinkSetDb"]["Link"][34]["Score"]==30130533
        assert record[0]["LinkSetDb"]["Link"][35]["Id"]=="15505294"
        assert record[0]["LinkSetDb"]["Link"][35]["Score"]==29430378
        assert record[0]["LinkSetDb"]["Link"][36]["Id"]=="17268692"
        assert record[0]["LinkSetDb"]["Link"][36]["Score"]==29166153
        assert record[0]["LinkSetDb"]["Link"][37]["Id"]=="11329655"
        assert record[0]["LinkSetDb"]["Link"][37]["Score"]==29112282
        assert record[0]["LinkSetDb"]["Link"][38]["Id"]=="11775722"
        assert record[0]["LinkSetDb"]["Link"][38]["Score"]==28940754
        assert record[0]["LinkSetDb"]["Link"][39]["Id"]=="11907356"
        assert record[0]["LinkSetDb"]["Link"][39]["Score"]==28860163
        assert record[0]["LinkSetDb"]["Link"][40]["Id"]=="10222515"
        assert record[0]["LinkSetDb"]["Link"][40]["Score"]==28807143
        assert record[0]["LinkSetDb"]["Link"][41]["Id"]=="17174054"
        assert record[0]["LinkSetDb"]["Link"][41]["Score"]==28790302
        assert record[0]["LinkSetDb"]["Link"][42]["Id"]=="9314960"
        assert record[0]["LinkSetDb"]["Link"][42]["Score"]==28750160
        assert record[0]["LinkSetDb"]["Link"][43]["Id"]=="14661661"
        assert record[0]["LinkSetDb"]["Link"][43]["Score"]==28361423
        assert record[0]["LinkSetDb"]["Link"][44]["Id"]=="17879696"
        assert record[0]["LinkSetDb"]["Link"][44]["Score"]==28120568
        assert record[0]["LinkSetDb"]["Link"][45]["Id"]=="4818442"
        assert record[0]["LinkSetDb"]["Link"][45]["Score"]==28058957
        assert record[0]["LinkSetDb"]["Link"][46]["Id"]=="15141648"
        assert record[0]["LinkSetDb"]["Link"][46]["Score"]==28011681
        assert record[0]["LinkSetDb"]["Link"][47]["Id"]=="8855688"
        assert record[0]["LinkSetDb"]["Link"][47]["Score"]==27711822
        assert record[0]["LinkSetDb"]["Link"][48]["Id"]=="17875143"
        assert record[0]["LinkSetDb"]["Link"][48]["Score"]==27711025
        assert record[0]["LinkSetDb"]["Link"][49]["Id"]=="1481295"
        assert record[0]["LinkSetDb"]["Link"][49]["Score"]==27707751
        assert record[0]["LinkSetDb"]["Link"][50]["Id"]=="8599783"
        assert record[0]["LinkSetDb"]["Link"][50]["Score"]==27683273
        assert record[0]["LinkSetDb"]["Link"][51]["Id"]=="10499696"
        assert record[0]["LinkSetDb"]["Link"][51]["Score"]==27623848
        assert record[0]["LinkSetDb"]["Link"][52]["Id"]=="12733684"
        assert record[0]["LinkSetDb"]["Link"][52]["Score"]==27527242
        assert record[0]["LinkSetDb"]["Link"][53]["Id"]=="18021675"
        assert record[0]["LinkSetDb"]["Link"][53]["Score"]==27495074
        assert record[0]["LinkSetDb"]["Link"][54]["Id"]=="12226761"
        assert record[0]["LinkSetDb"]["Link"][54]["Score"]==27366064
        assert record[0]["LinkSetDb"]["Link"][55]["Id"]=="4808999"
        assert record[0]["LinkSetDb"]["Link"][55]["Score"]==27304472
        assert record[0]["LinkSetDb"]["Link"][56]["Id"]=="16988291"
        assert record[0]["LinkSetDb"]["Link"][56]["Score"]==27295295
        assert record[0]["LinkSetDb"]["Link"][57]["Id"]=="10575758"
        assert record[0]["LinkSetDb"]["Link"][57]["Score"]==27243181
        assert record[0]["LinkSetDb"]["Link"][58]["Id"]=="8903064"
        assert record[0]["LinkSetDb"]["Link"][58]["Score"]==27206664
        assert record[0]["LinkSetDb"]["Link"][59]["Id"]=="10811354"
        assert record[0]["LinkSetDb"]["Link"][59]["Score"]==27088219
        assert record[0]["LinkSetDb"]["Link"][60]["Id"]=="16096604"
        assert record[0]["LinkSetDb"]["Link"][60]["Score"]==26862979
        assert record[0]["LinkSetDb"]["Link"][61]["Id"]=="15788584"
        assert record[0]["LinkSetDb"]["Link"][61]["Score"]==26759584
        assert record[0]["LinkSetDb"]["Link"][62]["Id"]=="17376366"
        assert record[0]["LinkSetDb"]["Link"][62]["Score"]==26743241
        assert record[0]["LinkSetDb"]["Link"][63]["Id"]=="16566645"
        assert record[0]["LinkSetDb"]["Link"][63]["Score"]==26725076
        assert record[0]["LinkSetDb"]["Link"][64]["Id"]=="17259035"
        assert record[0]["LinkSetDb"]["Link"][64]["Score"]==26595433
        assert record[0]["LinkSetDb"]["Link"][65]["Id"]=="9314959"
        assert record[0]["LinkSetDb"]["Link"][65]["Score"]==26445900
        assert record[0]["LinkSetDb"]["Link"][66]["Id"]=="11895298"
        assert record[0]["LinkSetDb"]["Link"][66]["Score"]==26256774
        assert record[0]["LinkSetDb"]["Link"][67]["Id"]=="11740602"
        assert record[0]["LinkSetDb"]["Link"][67]["Score"]==26158189
        assert record[0]["LinkSetDb"]["Link"][68]["Id"]=="15022983"
        assert record[0]["LinkSetDb"]["Link"][68]["Score"]==25889186
        assert record[0]["LinkSetDb"]["Link"][69]["Id"]=="15300544"
        assert record[0]["LinkSetDb"]["Link"][69]["Score"]==25837458
        assert record[0]["LinkSetDb"]["Link"][70]["Id"]=="12719915"
        assert record[0]["LinkSetDb"]["Link"][70]["Score"]==25831232
        assert record[0]["LinkSetDb"]["Link"][71]["Id"]=="14661306"
        assert record[0]["LinkSetDb"]["Link"][71]["Score"]==25788023
        assert record[0]["LinkSetDb"]["Link"][72]["Id"]=="16362812"
        assert record[0]["LinkSetDb"]["Link"][72]["Score"]==25565076
        assert record[0]["LinkSetDb"]["Link"][73]["Id"]=="17320773"
        assert record[0]["LinkSetDb"]["Link"][73]["Score"]==25504305
        assert record[0]["LinkSetDb"]["Link"][74]["Id"]=="11762248"
        assert record[0]["LinkSetDb"]["Link"][74]["Score"]==25504002
        assert record[0]["LinkSetDb"]["Link"][75]["Id"]=="10665303"
        assert record[0]["LinkSetDb"]["Link"][75]["Score"]==25384388
        assert record[0]["LinkSetDb"]["Link"][76]["Id"]=="17453494"
        assert record[0]["LinkSetDb"]["Link"][76]["Score"]==25226372
        assert record[0]["LinkSetDb"]["Link"][77]["Id"]=="9575723"
        assert record[0]["LinkSetDb"]["Link"][77]["Score"]==25174136
        assert record[0]["LinkSetDb"]["Link"][78]["Id"]=="12744498"
        assert record[0]["LinkSetDb"]["Link"][78]["Score"]==24971179
        assert record[0]["LinkSetDb"]["Link"][79]["Id"]=="12352163"
        assert record[0]["LinkSetDb"]["Link"][79]["Score"]==24915990
        assert record[0]["LinkSetDb"]["Link"][80]["Id"]=="8290724"
        assert record[0]["LinkSetDb"]["Link"][80]["Score"]==24909462
        assert record[0]["LinkSetDb"]["Link"][81]["Id"]=="11973504"
        assert record[0]["LinkSetDb"]["Link"][81]["Score"]==24878058
        assert record[0]["LinkSetDb"]["Link"][82]["Id"]=="14661668"
        assert record[0]["LinkSetDb"]["Link"][82]["Score"]==24779779
        assert record[0]["LinkSetDb"]["Link"][83]["Id"]=="16552382"
        assert record[0]["LinkSetDb"]["Link"][83]["Score"]==24760919
        assert record[0]["LinkSetDb"]["Link"][84]["Id"]=="17709829"
        assert record[0]["LinkSetDb"]["Link"][84]["Score"]==24743292
        assert record[0]["LinkSetDb"]["Link"][85]["Id"]=="14528718"
        assert record[0]["LinkSetDb"]["Link"][85]["Score"]==24686212
        assert record[0]["LinkSetDb"]["Link"][86]["Id"]=="15008163"
        assert record[0]["LinkSetDb"]["Link"][86]["Score"]==24612994
        assert record[0]["LinkSetDb"]["Link"][87]["Id"]=="10051883"
        assert record[0]["LinkSetDb"]["Link"][87]["Score"]==24492331
        assert record[0]["LinkSetDb"]["Link"][88]["Id"]=="11027076"
        assert record[0]["LinkSetDb"]["Link"][88]["Score"]==24410525
        assert record[0]["LinkSetDb"]["Link"][89]["Id"]=="17543650"
        assert record[0]["LinkSetDb"]["Link"][89]["Score"]==24371825
        assert record[0]["LinkSetDb"]["Link"][90]["Id"]=="17658095"
        assert record[0]["LinkSetDb"]["Link"][90]["Score"]==24331965
        assert record[0]["LinkSetDb"]["Link"][91]["Id"]=="9193407"
        assert record[0]["LinkSetDb"]["Link"][91]["Score"]==24240252
        assert record[0]["LinkSetDb"]["Link"][92]["Id"]=="10578418"
        assert record[0]["LinkSetDb"]["Link"][92]["Score"]==24091226
        assert record[0]["LinkSetDb"]["Link"][93]["Id"]=="12592155"
        assert record[0]["LinkSetDb"]["Link"][93]["Score"]==24001341
        assert record[0]["LinkSetDb"]["Link"][94]["Id"]=="17157468"
        assert record[0]["LinkSetDb"]["Link"][94]["Score"]==23984321
        assert record[0]["LinkSetDb"]["Link"][95]["Id"]=="15094630"
        assert record[0]["LinkSetDb"]["Link"][95]["Score"]==23912874
        assert record[0]["LinkSetDb"]["Link"][96]["Id"]=="8794574"
        assert record[0]["LinkSetDb"]["Link"][96]["Score"]==23900764
        assert record[0]["LinkSetDb"]["Link"][97]["Id"]=="9125660"
        assert record[0]["LinkSetDb"]["Link"][97]["Score"]==23884352
        assert record[0]["LinkSetDb"]["Link"][98]["Id"]=="8819381"
        assert record[0]["LinkSetDb"]["Link"][98]["Score"]==23839719
        assert record[0]["LinkSetDb"]["Link"][99]["Id"]=="14661666"
        assert record[0]["LinkSetDb"]["Link"][99]["Score"]==23748510
        assert record[0]["LinkSetDb"]["Link"][100]["Id"]=="9658901"
        assert record[0]["LinkSetDb"]["Link"][100]["Score"]==23667126
        assert record[0]["LinkSetDb"]["Link"][101]["Id"]=="12744499"
        assert record[0]["LinkSetDb"]["Link"][101]["Score"]==23647189
        assert record[0]["LinkSetDb"]["Link"][102]["Id"]=="12164574"
        assert record[0]["LinkSetDb"]["Link"][102]["Score"]==23623853
        assert record[0]["LinkSetDb"]["Link"][103]["Id"]=="15136027"
        assert record[0]["LinkSetDb"]["Link"][103]["Score"]==23572558
        assert record[0]["LinkSetDb"]["Link"][104]["Id"]=="14872380"
        assert record[0]["LinkSetDb"]["Link"][104]["Score"]==23460906
        assert record[0]["LinkSetDb"]["Link"][105]["Id"]=="3905087"
        assert record[0]["LinkSetDb"]["Link"][105]["Score"]==23305022
        assert record[0]["LinkSetDb"]["Link"][106]["Id"]=="15642291"
        assert record[0]["LinkSetDb"]["Link"][106]["Score"]==23234831
        assert record[0]["LinkSetDb"]["Link"][107]["Id"]=="16928974"
        assert record[0]["LinkSetDb"]["Link"][107]["Score"]==23223298
        assert record[0]["LinkSetDb"]["Link"][108]["Id"]=="6072516"
        assert record[0]["LinkSetDb"]["Link"][108]["Score"]==23042548
        assert record[0]["LinkSetDb"]["Link"][109]["Id"]=="12949462"
        assert record[0]["LinkSetDb"]["Link"][109]["Score"]==23001441
        assert record[0]["LinkSetDb"]["Link"][110]["Id"]=="10761553"
        assert record[0]["LinkSetDb"]["Link"][110]["Score"]==22995991
        assert record[0]["LinkSetDb"]["Link"][111]["Id"]=="14661663"
        assert record[0]["LinkSetDb"]["Link"][111]["Score"]==22986720
        assert record[0]["LinkSetDb"]["Link"][112]["Id"]=="16338316"
        assert record[0]["LinkSetDb"]["Link"][112]["Score"]==22933288
        assert record[0]["LinkSetDb"]["Link"][113]["Id"]=="17464254"
        assert record[0]["LinkSetDb"]["Link"][113]["Score"]==22912253
        assert record[0]["LinkSetDb"]["Link"][114]["Id"]=="15529836"
        assert record[0]["LinkSetDb"]["Link"][114]["Score"]==22892154
        assert record[0]["LinkSetDb"]["Link"][115]["Id"]=="12361530"
        assert record[0]["LinkSetDb"]["Link"][115]["Score"]==22871698
        assert record[0]["LinkSetDb"]["Link"][116]["Id"]=="12876813"
        assert record[0]["LinkSetDb"]["Link"][116]["Score"]==22822515
        assert record[0]["LinkSetDb"]["Link"][117]["Id"]=="10749221"
        assert record[0]["LinkSetDb"]["Link"][117]["Score"]==22794373
        assert record[0]["LinkSetDb"]["Link"][118]["Id"]=="6482054"
        assert record[0]["LinkSetDb"]["Link"][118]["Score"]==22791927
        assert record[0]["LinkSetDb"]["Link"][119]["Id"]=="9016217"
        assert record[0]["LinkSetDb"]["Link"][119]["Score"]==22738432
        assert record[0]["LinkSetDb"]["Link"][120]["Id"]=="14702442"
        assert record[0]["LinkSetDb"]["Link"][120]["Score"]==22722123
        assert record[0]["LinkSetDb"]["Link"][121]["Id"]=="15279747"
        assert record[0]["LinkSetDb"]["Link"][121]["Score"]==22698787
        assert record[0]["LinkSetDb"]["Link"][122]["Id"]=="7892443"
        assert record[0]["LinkSetDb"]["Link"][122]["Score"]==22642038
        assert record[0]["LinkSetDb"]["Link"][123]["Id"]=="616459"
        assert record[0]["LinkSetDb"]["Link"][123]["Score"]==22591277
        assert record[0]["LinkSetDb"]["Link"][124]["Id"]=="8886718"
        assert record[0]["LinkSetDb"]["Link"][124]["Score"]==22542938
        assert record[0]["LinkSetDb"]["Link"][125]["Id"]=="17245521"
        assert record[0]["LinkSetDb"]["Link"][125]["Score"]==22538649
        assert record[0]["LinkSetDb"]["Link"][126]["Id"]=="1535863"
        assert record[0]["LinkSetDb"]["Link"][126]["Score"]==22468774
        assert record[0]["LinkSetDb"]["Link"][127]["Id"]=="15537403"
        assert record[0]["LinkSetDb"]["Link"][127]["Score"]==22458002
        assert record[0]["LinkSetDb"]["Link"][128]["Id"]=="16040910"
        assert record[0]["LinkSetDb"]["Link"][128]["Score"]==22452119
        assert record[0]["LinkSetDb"]["Link"][129]["Id"]=="16929028"
        assert record[0]["LinkSetDb"]["Link"][129]["Score"]==22433988
        assert record[0]["LinkSetDb"]["Link"][130]["Id"]=="16697589"
        assert record[0]["LinkSetDb"]["Link"][130]["Score"]==22366606
        assert record[0]["LinkSetDb"]["Link"][131]["Id"]=="531835"
        assert record[0]["LinkSetDb"]["Link"][131]["Score"]==22366454
        assert record[0]["LinkSetDb"]["Link"][132]["Id"]=="2308313"
        assert record[0]["LinkSetDb"]["Link"][132]["Score"]==22330898
        assert record[0]["LinkSetDb"]["Link"][133]["Id"]=="12522920"
        assert record[0]["LinkSetDb"]["Link"][133]["Score"]==22178764
        assert record[0]["LinkSetDb"]["Link"][134]["Id"]=="10222521"
        assert record[0]["LinkSetDb"]["Link"][134]["Score"]==22135023
        assert record[0]["LinkSetDb"]["Link"][135]["Id"]=="10499697"
        assert record[0]["LinkSetDb"]["Link"][135]["Score"]==22130302
        assert record[0]["LinkSetDb"]["Link"][136]["Id"]=="8903058"
        assert record[0]["LinkSetDb"]["Link"][136]["Score"]==22113132
        assert record[0]["LinkSetDb"]["Link"][137]["Id"]=="17441569"
        assert record[0]["LinkSetDb"]["Link"][137]["Score"]==22085858
        assert record[0]["LinkSetDb"]["Link"][138]["Id"]=="15284932"
        assert record[0]["LinkSetDb"]["Link"][138]["Score"]==22075791
        assert record[0]["LinkSetDb"]["Link"][139]["Id"]=="15466771"
        assert record[0]["LinkSetDb"]["Link"][139]["Score"]==22075418
        assert record[0]["LinkSetDb"]["Link"][140]["Id"]=="17145267"
        assert record[0]["LinkSetDb"]["Link"][140]["Score"]==22033864
        assert record[0]["LinkSetDb"]["Link"][141]["Id"]=="11329662"
        assert record[0]["LinkSetDb"]["Link"][141]["Score"]==22012948
        assert record[0]["LinkSetDb"]["Link"][142]["Id"]=="10222514"
        assert record[0]["LinkSetDb"]["Link"][142]["Score"]==22009777
        assert record[0]["LinkSetDb"]["Link"][143]["Id"]=="17383530"
        assert record[0]["LinkSetDb"]["Link"][143]["Score"]==22003600
        assert record[0]["LinkSetDb"]["Link"][144]["Id"]=="12455800"
        assert record[0]["LinkSetDb"]["Link"][144]["Score"]==21992674
        assert record[0]["LinkSetDb"]["Link"][145]["Id"]=="15845051"
        assert record[0]["LinkSetDb"]["Link"][145]["Score"]==21946257
        assert record[0]["LinkSetDb"]["Link"][146]["Id"]=="11443295"
        assert record[0]["LinkSetDb"]["Link"][146]["Score"]==21908841
        assert record[0]["LinkSetDb"]["Link"][147]["Id"]=="15162233"
        assert record[0]["LinkSetDb"]["Link"][147]["Score"]==21903624
        assert record[0]["LinkSetDb"]["Link"][148]["Id"]=="16133610"
        assert record[0]["LinkSetDb"]["Link"][148]["Score"]==21872203
        assert record[0]["LinkSetDb"]["Link"][149]["Id"]=="12845461"
        assert record[0]["LinkSetDb"]["Link"][149]["Score"]==21864314
        assert record[0]["LinkSetDb"]["Link"][150]["Id"]=="16947073"
        assert record[0]["LinkSetDb"]["Link"][150]["Score"]==21832153
        assert record[0]["LinkSetDb"]["Link"][151]["Id"]=="7415301"
        assert record[0]["LinkSetDb"]["Link"][151]["Score"]==21822396
        assert record[0]["LinkSetDb"]["Link"][152]["Id"]=="16416239"
        assert record[0]["LinkSetDb"]["Link"][152]["Score"]==21820165
        assert record[0]["LinkSetDb"]["Link"][153]["Id"]=="4848922"
        assert record[0]["LinkSetDb"]["Link"][153]["Score"]==21786194
        assert record[0]["LinkSetDb"]["Link"][154]["Id"]=="12720164"
        assert record[0]["LinkSetDb"]["Link"][154]["Score"]==21785319
        assert record[0]["LinkSetDb"]["Link"][155]["Id"]=="17093987"
        assert record[0]["LinkSetDb"]["Link"][155]["Score"]==21750370
        assert record[0]["LinkSetDb"]["Link"][156]["Id"]=="16769006"
        assert record[0]["LinkSetDb"]["Link"][156]["Score"]==21735873
        assert record[0]["LinkSetDb"]["Link"][157]["Id"]=="17954835"
        assert record[0]["LinkSetDb"]["Link"][157]["Score"]==21733933
        assert record[0]["LinkSetDb"]["Link"][158]["Id"]=="15236134"
        assert record[0]["LinkSetDb"]["Link"][158]["Score"]==21640099
        assert record[0]["LinkSetDb"]["Link"][159]["Id"]=="12524603"
        assert record[0]["LinkSetDb"]["Link"][159]["Score"]==21636724
        assert record[0]["LinkSetDb"]["Link"][160]["Id"]=="16749985"
        assert record[0]["LinkSetDb"]["Link"][160]["Score"]==21628926
        assert record[0]["LinkSetDb"]["Link"][161]["Id"]=="3213296"
        assert record[0]["LinkSetDb"]["Link"][161]["Score"]==21490232
        assert record[0]["LinkSetDb"]["Link"][162]["Id"]=="11409026"
        assert record[0]["LinkSetDb"]["Link"][162]["Score"]==21061296
        assert record[0]["LinkSetDb"]["Link"][163]["Id"]=="9725288"
        assert record[0]["LinkSetDb"]["Link"][163]["Score"]==21053585
        assert record[0]["LinkSetDb"]["Link"][164]["Id"]=="6217136"
        assert record[0]["LinkSetDb"]["Link"][164]["Score"]==21042914
        assert record[0]["LinkSetDb"]["Link"][165]["Id"]=="663071"
        assert record[0]["LinkSetDb"]["Link"][165]["Score"]==20926141
        assert record[0]["LinkSetDb"]["Link"][166]["Id"]=="10341802"
        assert record[0]["LinkSetDb"]["Link"][166]["Score"]==20797282
        assert record[0]["LinkSetDb"]["Link"][167]["Id"]=="6473764"
        assert record[0]["LinkSetDb"]["Link"][167]["Score"]==20757680
        assert record[0]["LinkSetDb"]["Link"][168]["Id"]=="2584497"
        assert record[0]["LinkSetDb"]["Link"][168]["Score"]==20521350
        assert record[0]["LinkSetDb"]["Link"][169]["Id"]=="8338105"
        assert record[0]["LinkSetDb"]["Link"][169]["Score"]==20501334
        assert record[0]["LinkSetDb"]["Link"][170]["Id"]=="18053822"
        assert record[0]["LinkSetDb"]["Link"][170]["Score"]==20275078
        assert record[0]["LinkSetDb"]["Link"][171]["Id"]=="4058411"
        assert record[0]["LinkSetDb"]["Link"][171]["Score"]==20161667
        assert record[0]["LinkSetDb"]["Link"][172]["Id"]=="11669077"
        assert record[0]["LinkSetDb"]["Link"][172]["Score"]==19993282
        assert record[0]["LinkSetDb"]["Link"][173]["Id"]=="11781922"
        assert record[0]["LinkSetDb"]["Link"][173]["Score"]==19969425
        assert record[0]["LinkSetDb"]["Link"][174]["Id"]=="9793138"
        assert record[0]["LinkSetDb"]["Link"][174]["Score"]==19952972
        assert record[0]["LinkSetDb"]["Link"][175]["Id"]=="9391495"
        assert record[0]["LinkSetDb"]["Link"][175]["Score"]==19815538
        assert record[0]["LinkSetDb"]["Link"][176]["Id"]=="10803203"
        assert record[0]["LinkSetDb"]["Link"][176]["Score"]==19495693
        assert record[0]["LinkSetDb"]["Link"][177]["Id"]=="7326186"
        assert record[0]["LinkSetDb"]["Link"][177]["Score"]==19273989
        assert record[0]["LinkSetDb"]["Link"][178]["Id"]=="11868066"
        assert record[0]["LinkSetDb"]["Link"][178]["Score"]==19220137
        assert record[0]["LinkSetDb"]["Link"][179]["Id"]=="10904988"
        assert record[0]["LinkSetDb"]["Link"][179]["Score"]==19203510
        assert record[0]["LinkSetDb"]["Link"][180]["Id"]=="3288780"
        assert record[0]["LinkSetDb"]["Link"][180]["Score"]==18958114
        assert record[0]["LinkSetDb"]["Link"][181]["Id"]=="2047316"
        assert record[0]["LinkSetDb"]["Link"][181]["Score"]==18907473
        assert record[0]["LinkSetDb"]["Link"][182]["Id"]=="12237004"
        assert record[0]["LinkSetDb"]["Link"][182]["Score"]==18751474
        assert record[0]["LinkSetDb"]["Link"][183]["Id"]=="5627987"
        assert record[0]["LinkSetDb"]["Link"][183]["Score"]==18741903
        assert record[0]["LinkSetDb"]["Link"][184]["Id"]=="9269670"
        assert record[0]["LinkSetDb"]["Link"][184]["Score"]==18666426
        assert record[0]["LinkSetDb"]["Link"][185]["Id"]=="8903059"
        assert record[0]["LinkSetDb"]["Link"][185]["Score"]==18653874
        assert record[0]["LinkSetDb"]["Link"][186]["Id"]=="5594242"
        assert record[0]["LinkSetDb"]["Link"][186]["Score"]==18548780
        assert record[0]["LinkSetDb"]["Link"][187]["Id"]=="7068417"
        assert record[0]["LinkSetDb"]["Link"][187]["Score"]==18390022
        assert record[0]["LinkSetDb"]["Link"][188]["Id"]=="7330196"
        assert record[0]["LinkSetDb"]["Link"][188]["Score"]==18371587
        assert record[0]["LinkSetDb"]["Link"][189]["Id"]=="7408592"
        assert record[0]["LinkSetDb"]["Link"][189]["Score"]==18275541
        assert record[0]["LinkSetDb"]["Link"][190]["Id"]=="8835983"
        assert record[0]["LinkSetDb"]["Link"][190]["Score"]==18176923
        assert record[0]["LinkSetDb"]["Link"][191]["Id"]=="6940010"
        assert record[0]["LinkSetDb"]["Link"][191]["Score"]==18011066
        assert record[0]["LinkSetDb"]["Link"][192]["Id"]=="10499712"
        assert record[0]["LinkSetDb"]["Link"][192]["Score"]==17943586
        assert record[0]["LinkSetDb"]["Link"][193]["Id"]=="4539876"
        assert record[0]["LinkSetDb"]["Link"][193]["Score"]==17915154
        assert record[0]["LinkSetDb"]["Link"][194]["Id"]=="1943587"
        assert record[0]["LinkSetDb"]["Link"][194]["Score"]==17752606
        assert record[0]["LinkSetDb"]["Link"][195]["Id"]=="9847909"
        assert record[0]["LinkSetDb"]["Link"][195]["Score"]==17568386
        assert record[0]["LinkSetDb"]["Link"][196]["Id"]=="11578071"
        assert record[0]["LinkSetDb"]["Link"][196]["Score"]==17561413
        assert record[0]["LinkSetDb"]["Link"][197]["Id"]=="11789473"
        assert record[0]["LinkSetDb"]["Link"][197]["Score"]==17435433
        assert record[0]["LinkSetDb"]["Link"][198]["Id"]=="9885599"
        assert record[0]["LinkSetDb"]["Link"][198]["Score"]==17383598
        assert record[0]["LinkSetDb"]["Link"][199]["Id"]=="7423836"
        assert record[0]["LinkSetDb"]["Link"][199]["Score"]==17196872
        assert record[0]["LinkSetDb"]["Link"][200]["Id"]=="10688063"
        assert record[0]["LinkSetDb"]["Link"][200]["Score"]==16453112
        assert record[0]["LinkSetDb"]["Link"][201]["Id"]=="11695100"
        assert record[0]["LinkSetDb"]["Link"][201]["Score"]==16352760
        assert record[0]["LinkSetDb"]["Link"][202]["Id"]=="11329658"
        assert record[0]["LinkSetDb"]["Link"][202]["Score"]==16089885
        assert record[0]["LinkSetDb"]["Link"][203]["Id"]=="11939665"
        assert record[0]["LinkSetDb"]["Link"][203]["Score"]==15947974
        assert record[0]["LinkSetDb"]["Link"][204]["Id"]=="5512349"
        assert record[0]["LinkSetDb"]["Link"][204]["Score"]==15647685
        assert record[0]["LinkSetDb"]["Link"][205]["Id"]=="2222794"
        assert record[0]["LinkSetDb"]["Link"][205]["Score"]==14981157
        assert record[0]["LinkSetDb"]["Link"][206]["Id"]=="5998281"
        assert record[0]["LinkSetDb"]["Link"][206]["Score"]==14226588
        assert record[0]["LinkSetDb"]["Link"][207]["Id"]=="10475937"
        assert record[0]["LinkSetDb"]["Link"][207]["Score"]==13934390
        assert record[0]["LinkSetDb"]["Link"][208]["Id"]=="5046513"
        assert record[0]["LinkSetDb"]["Link"][208]["Score"]==12769605
        assert record[0]["LinkSetDb"]["Link"][209]["Id"]=="1539132"
        assert record[0]["LinkSetDb"]["Link"][209]["Score"]==12395064
        assert record[0]["LinkSetDb"]["Link"][210]["Id"]=="4414214"
        assert record[0]["LinkSetDb"]["Link"][210]["Score"]==10113539

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
        assert record[0]["IdCheckList"][0]["Id"]=="12169658"
        assert record[0]["IdCheckList"][0]["LinkInfo"][0]["DbTo"]=="books"
        assert record[0]["IdCheckList"][0]["LinkInfo"][0]["LinkName"]=="pubmed_books_refs"
        assert record[0]["IdCheckList"][0]["LinkInfo"][0]["MenuTag"]=="Cited in Books"
        assert record[0]["IdCheckList"][0]["LinkInfo"][0]["HtmlTag"]=="Cited in Books"
        assert record[0]["IdCheckList"][0]["LinkInfo"][0]["Priority"]==185
        assert record[0]["IdCheckList"][0]["LinkInfo"][1]["DbTo"]=="gene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][1]["LinkName"]=="pubmed_gene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][1]["MenuTag"]=="Gene Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][1]["HtmlTag"]=="Gene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][1]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][2]["DbTo"]=="geo"
        assert record[0]["IdCheckList"][0]["LinkInfo"][2]["LinkName"]=="pubmed_geo"
        assert record[0]["IdCheckList"][0]["LinkInfo"][2]["MenuTag"]=="GEO Profile Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][2]["HtmlTag"]=="GEO Profiles"
        assert record[0]["IdCheckList"][0]["LinkInfo"][2]["Priority"]==170
        assert record[0]["IdCheckList"][0]["LinkInfo"][3]["DbTo"]=="homologene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][3]["LinkName"]=="pubmed_homologene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][3]["MenuTag"]=="HomoloGene Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][3]["HtmlTag"]=="HomoloGene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][3]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][4]["DbTo"]=="nuccore"
        assert record[0]["IdCheckList"][0]["LinkInfo"][4]["LinkName"]=="pubmed_nuccore"
        assert record[0]["IdCheckList"][0]["LinkInfo"][4]["MenuTag"]=="CoreNucleotide Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][4]["HtmlTag"]=="CoreNucleotide"
        assert record[0]["IdCheckList"][0]["LinkInfo"][4]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][5]["DbTo"]=="nuccore"
        assert record[0]["IdCheckList"][0]["LinkInfo"][5]["LinkName"]=="pubmed_nuccore_refseq"
        assert record[0]["IdCheckList"][0]["LinkInfo"][5]["MenuTag"]=="CoreNucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][5]["HtmlTag"]=="CoreNucleotide (RefSeq)"
        assert record[0]["IdCheckList"][0]["LinkInfo"][5]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][6]["DbTo"]=="nucleotide"
        assert record[0]["IdCheckList"][0]["LinkInfo"][6]["LinkName"]=="pubmed_nucleotide"
        assert record[0]["IdCheckList"][0]["LinkInfo"][6]["MenuTag"]=="Nucleotide Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][6]["HtmlTag"]=="Nucleotide"
        assert record[0]["IdCheckList"][0]["LinkInfo"][6]["Priority"]==135
        assert record[0]["IdCheckList"][0]["LinkInfo"][7]["DbTo"]=="nucleotide"
        assert record[0]["IdCheckList"][0]["LinkInfo"][7]["LinkName"]=="pubmed_nucleotide_refseq"
        assert record[0]["IdCheckList"][0]["LinkInfo"][7]["MenuTag"]=="Nucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][7]["HtmlTag"]=="Nucleotide (RefSeq)"
        assert record[0]["IdCheckList"][0]["LinkInfo"][7]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][8]["DbTo"]=="pcsubstance"
        assert record[0]["IdCheckList"][0]["LinkInfo"][8]["LinkName"]=="pubmed_pcsubstance_mesh"
        assert record[0]["IdCheckList"][0]["LinkInfo"][8]["MenuTag"]=="Substance (MeSH Keyword)"
        assert record[0]["IdCheckList"][0]["LinkInfo"][8]["HtmlTag"]=="Substance (MeSH Keyword)"
        assert record[0]["IdCheckList"][0]["LinkInfo"][8]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][9]["DbTo"]=="pmc"
        assert record[0]["IdCheckList"][0]["LinkInfo"][9]["LinkName"]=="pubmed_pmc_refs"
        assert record[0]["IdCheckList"][0]["LinkInfo"][9]["MenuTag"]=="Cited in PMC"
        assert record[0]["IdCheckList"][0]["LinkInfo"][9]["HtmlTag"]=="Cited in PMC"
        assert record[0]["IdCheckList"][0]["LinkInfo"][9]["Url"]=="http://www.pubmedcentral.gov/tocrender.fcgi?action=cited&tool=pubmed&pubmedid=<@UID@>"
        assert record[0]["IdCheckList"][0]["LinkInfo"][9]["Priority"]==180
        assert record[0]["IdCheckList"][0]["LinkInfo"][10]["DbTo"]=="protein"
        assert record[0]["IdCheckList"][0]["LinkInfo"][10]["LinkName"]=="pubmed_protein"
        assert record[0]["IdCheckList"][0]["LinkInfo"][10]["MenuTag"]=="Protein Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][10]["HtmlTag"]=="Protein"
        assert record[0]["IdCheckList"][0]["LinkInfo"][10]["Priority"]==140
        assert record[0]["IdCheckList"][0]["LinkInfo"][11]["DbTo"]=="protein"
        assert record[0]["IdCheckList"][0]["LinkInfo"][11]["LinkName"]=="pubmed_protein_refseq"
        assert record[0]["IdCheckList"][0]["LinkInfo"][11]["MenuTag"]=="Protein (RefSeq) Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][11]["HtmlTag"]=="Protein (RefSeq)"
        assert record[0]["IdCheckList"][0]["LinkInfo"][11]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][12]["DbTo"]=="pubmed"
        assert record[0]["IdCheckList"][0]["LinkInfo"][12]["LinkName"]=="pubmed_pubmed"
        assert record[0]["IdCheckList"][0]["LinkInfo"][12]["MenuTag"]=="Related Articles"
        assert record[0]["IdCheckList"][0]["LinkInfo"][12]["HtmlTag"]=="Related Articles"
        assert record[0]["IdCheckList"][0]["LinkInfo"][12]["Priority"]==1
        assert record[0]["IdCheckList"][0]["LinkInfo"][13]["DbTo"]=="taxonomy"
        assert record[0]["IdCheckList"][0]["LinkInfo"][13]["LinkName"]=="pubmed_taxonomy_entrez"
        assert record[0]["IdCheckList"][0]["LinkInfo"][13]["MenuTag"]=="Taxonomy via GenBank"
        assert record[0]["IdCheckList"][0]["LinkInfo"][13]["HtmlTag"]=="Taxonomy via GenBank"
        assert record[0]["IdCheckList"][0]["LinkInfo"][13]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][14]["DbTo"]=="unigene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][14]["LinkName"]=="pubmed_unigene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][14]["MenuTag"]=="UniGene Links"
        assert record[0]["IdCheckList"][0]["LinkInfo"][14]["HtmlTag"]=="UniGene"
        assert record[0]["IdCheckList"][0]["LinkInfo"][14]["Priority"]==128
        assert record[0]["IdCheckList"][0]["LinkInfo"][15]["DbTo"]=="LinkOut"
        assert record[0]["IdCheckList"][0]["LinkInfo"][15]["LinkName"]=="ExternalLink"
        assert record[0]["IdCheckList"][0]["LinkInfo"][15]["MenuTag"]=="LinkOut"
        assert record[0]["IdCheckList"][0]["LinkInfo"][15]["HtmlTag"]=="LinkOut"
        assert record[0]["IdCheckList"][0]["LinkInfo"][15]["Priority"]==255
        assert record[0]["IdCheckList"][1]["Id"]=="11748140"
        assert record[0]["IdCheckList"][1]["LinkInfo"][0]["DbTo"]=="books"
        assert record[0]["IdCheckList"][1]["LinkInfo"][0]["LinkName"]=="pubmed_books_refs"
        assert record[0]["IdCheckList"][1]["LinkInfo"][0]["MenuTag"]=="Cited in Books"
        assert record[0]["IdCheckList"][1]["LinkInfo"][0]["HtmlTag"]=="Cited in Books"
        assert record[0]["IdCheckList"][1]["LinkInfo"][0]["Priority"]==185
        assert record[0]["IdCheckList"][1]["LinkInfo"][1]["DbTo"]=="gene"
        assert record[0]["IdCheckList"][1]["LinkInfo"][1]["LinkName"]=="pubmed_gene"
        assert record[0]["IdCheckList"][1]["LinkInfo"][1]["MenuTag"]=="Gene Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][1]["HtmlTag"]=="Gene"
        assert record[0]["IdCheckList"][1]["LinkInfo"][1]["Priority"]==128
        assert record[0]["IdCheckList"][1]["LinkInfo"][2]["DbTo"]=="geo"
        assert record[0]["IdCheckList"][1]["LinkInfo"][2]["LinkName"]=="pubmed_geo"
        assert record[0]["IdCheckList"][1]["LinkInfo"][2]["MenuTag"]=="GEO Profile Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][2]["HtmlTag"]=="GEO Profiles"
        assert record[0]["IdCheckList"][1]["LinkInfo"][2]["Priority"]==170
        assert record[0]["IdCheckList"][1]["LinkInfo"][3]["DbTo"]=="nuccore"
        assert record[0]["IdCheckList"][1]["LinkInfo"][3]["LinkName"]=="pubmed_nuccore"
        assert record[0]["IdCheckList"][1]["LinkInfo"][3]["MenuTag"]=="CoreNucleotide Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][3]["HtmlTag"]=="CoreNucleotide"
        assert record[0]["IdCheckList"][1]["LinkInfo"][3]["Priority"]==128
        assert record[0]["IdCheckList"][1]["LinkInfo"][4]["DbTo"]=="nuccore"
        assert record[0]["IdCheckList"][1]["LinkInfo"][4]["LinkName"]=="pubmed_nuccore_refseq"
        assert record[0]["IdCheckList"][1]["LinkInfo"][4]["MenuTag"]=="CoreNucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][4]["HtmlTag"]=="CoreNucleotide (RefSeq)"
        assert record[0]["IdCheckList"][1]["LinkInfo"][4]["Priority"]==128
        assert record[0]["IdCheckList"][1]["LinkInfo"][5]["DbTo"]=="nucleotide"
        assert record[0]["IdCheckList"][1]["LinkInfo"][5]["LinkName"]=="pubmed_nucleotide"
        assert record[0]["IdCheckList"][1]["LinkInfo"][5]["MenuTag"]=="Nucleotide Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][5]["HtmlTag"]=="Nucleotide"
        assert record[0]["IdCheckList"][1]["LinkInfo"][5]["Priority"]==135
        assert record[0]["IdCheckList"][1]["LinkInfo"][6]["DbTo"]=="nucleotide"
        assert record[0]["IdCheckList"][1]["LinkInfo"][6]["LinkName"]=="pubmed_nucleotide_refseq"
        assert record[0]["IdCheckList"][1]["LinkInfo"][6]["MenuTag"]=="Nucleotide (RefSeq) Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][6]["HtmlTag"]=="Nucleotide (RefSeq)"
        assert record[0]["IdCheckList"][1]["LinkInfo"][6]["Priority"]==128
        assert record[0]["IdCheckList"][1]["LinkInfo"][7]["DbTo"]=="pmc"
        assert record[0]["IdCheckList"][1]["LinkInfo"][7]["LinkName"]=="pubmed_pmc_refs"
        assert record[0]["IdCheckList"][1]["LinkInfo"][7]["MenuTag"]=="Cited in PMC"
        assert record[0]["IdCheckList"][1]["LinkInfo"][7]["HtmlTag"]=="Cited in PMC"
        assert record[0]["IdCheckList"][1]["LinkInfo"][7]["Url"]=="http://www.pubmedcentral.gov/tocrender.fcgi?action=cited&tool=pubmed&pubmedid=<@UID@>"
        assert record[0]["IdCheckList"][1]["LinkInfo"][7]["Priority"]==180
        assert record[0]["IdCheckList"][1]["LinkInfo"][8]["DbTo"]=="protein"
        assert record[0]["IdCheckList"][1]["LinkInfo"][8]["LinkName"]=="pubmed_protein"
        assert record[0]["IdCheckList"][1]["LinkInfo"][8]["MenuTag"]=="Protein Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][8]["HtmlTag"]=="Protein"
        assert record[0]["IdCheckList"][1]["LinkInfo"][8]["Priority"]==140
        assert record[0]["IdCheckList"][1]["LinkInfo"][9]["DbTo"]=="protein"
        assert record[0]["IdCheckList"][1]["LinkInfo"][9]["LinkName"]=="pubmed_protein_refseq"
        assert record[0]["IdCheckList"][1]["LinkInfo"][9]["MenuTag"]=="Protein (RefSeq) Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][9]["HtmlTag"]=="Protein (RefSeq)"
        assert record[0]["IdCheckList"][1]["LinkInfo"][9]["Priority"]==128
        assert record[0]["IdCheckList"][1]["LinkInfo"][10]["DbTo"]=="pubmed"
        assert record[0]["IdCheckList"][1]["LinkInfo"][10]["LinkName"]=="pubmed_pubmed"
        assert record[0]["IdCheckList"][1]["LinkInfo"][10]["MenuTag"]=="Related Articles"
        assert record[0]["IdCheckList"][1]["LinkInfo"][10]["HtmlTag"]=="Related Articles"
        assert record[0]["IdCheckList"][1]["LinkInfo"][10]["Priority"]==1
        assert record[0]["IdCheckList"][1]["LinkInfo"][11]["DbTo"]=="taxonomy"
        assert record[0]["IdCheckList"][1]["LinkInfo"][11]["LinkName"]=="pubmed_taxonomy_entrez"
        assert record[0]["IdCheckList"][1]["LinkInfo"][11]["MenuTag"]=="Taxonomy via GenBank"
        assert record[0]["IdCheckList"][1]["LinkInfo"][11]["HtmlTag"]=="Taxonomy via GenBank"
        assert record[0]["IdCheckList"][1]["LinkInfo"][11]["Priority"]==128
        assert record[0]["IdCheckList"][1]["LinkInfo"][12]["DbTo"]=="unigene"
        assert record[0]["IdCheckList"][1]["LinkInfo"][12]["LinkName"]=="pubmed_unigene"
        assert record[0]["IdCheckList"][1]["LinkInfo"][12]["MenuTag"]=="UniGene Links"
        assert record[0]["IdCheckList"][1]["LinkInfo"][12]["HtmlTag"]=="UniGene"
        assert record[0]["IdCheckList"][1]["LinkInfo"][12]["Priority"]==128
        assert record[0]["IdCheckList"][1]["LinkInfo"][13]["DbTo"]=="LinkOut"
        assert record[0]["IdCheckList"][1]["LinkInfo"][13]["LinkName"]=="ExternalLink"
        assert record[0]["IdCheckList"][1]["LinkInfo"][13]["MenuTag"]=="LinkOut"
        assert record[0]["IdCheckList"][1]["LinkInfo"][13]["HtmlTag"]=="LinkOut"
        assert record[0]["IdCheckList"][1]["LinkInfo"][13]["Priority"]==255


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
        assert len(record[0])==2
        assert record[0]["DbFrom"]=="pubmed"
        assert len(record[0]["IdCheckList"][0])
        assert record[0]["IdCheckList"][0]==["12068369", None, True]

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


if __name__ == '__main__':
    sys.exit(run_tests(sys.argv))
