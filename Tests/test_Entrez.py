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
    tests = [EInfoTest, ESearchTest, EPostTest, ESummaryTest]
    
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



if __name__ == '__main__':
    sys.exit(run_tests(sys.argv))
