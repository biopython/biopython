"""nlmmedline_xml_format.py

A Martel format to parse the NLM's XML format for Medline.

http://www.nlm.nih.gov/databases/dtd/nlmmedline_011101.dtd
http://www.nlm.nih.gov/databases/dtd/nlmmedlinecitation_011101.dtd
http://www.nlm.nih.gov/databases/dtd/nlmcommon_011101.dtd

Formats:
citation_format    Format for one MedlineCitation.
format             Format for a whole file.

"""
import sys

from Martel import *
from Martel import RecordReader

self = sys.modules[__name__]

def _start_elem(element, *attrs):
    if attrs:
        attr_groups = []
        for attr in attrs:
            group = Str(attr) + Str("=") + \
                    Str('"') + Group(attr, Re(r'[^<&"]+')) + Str('"')
            attr_groups.append(group)
        start = Str("<") + Str(element) + \
                Rep(Str(" ") + Alt(*attr_groups)) + \
                Str(">")
    else:
        start = Str("<%s>" % element)
    return start

def _end_elem(element):
    return Str("</%s>" % element)

def simple_elem(element, *attrs):
    """simple_elem(element, *attrs)
    
    Create a Martel Expression in this module's namespace that will
    recognize an XML element in the form of:
    <element>data</element>

    The whole element must be on a single line.  The Expression will
    be created in the module's namespace with the same name as the
    element.

    """
    start, end = _start_elem(element, *attrs), _end_elem(element)
    
    group_name = element
    group_expression = Re(r"[^<]+")
    expr = start + \
           Group(group_name, group_expression) + \
           end + \
           AnyEol()
    setattr(self, element, expr)


# Group expressions.  A group consists of the start and end elements
# with an expression in-between.  The Expression for the group will be
# called "NAME".
def group_elem(element, expr, *attrs):
    start_name, end_name = "%s_start" % element, "%s_end" % element
    start_expr = getattr(self, start_name, None)
    if start_expr is None:
        start_expr = _start_elem(element, *attrs) + AnyEol()
        setattr(self, start_name, start_expr)
    end_expr = getattr(self, end_name, None)
    if end_expr is None:
        end_expr = _end_elem(element) + AnyEol()
        setattr(self, end_name, end_expr)

    group_expr = start_expr + expr + end_expr
    group_expr = Group(element, group_expr)
    setattr(self, element, group_expr)


######################################################################
# Implement Martel expressions that recognize:                       #
# http://www.nlm.nih.gov/databases/dtd/nlmcommon_011101.dtd          #
######################################################################

########################################
# Personal and Author names
elements = [
    "FirstName", "ForeName", "MiddleName", "LastName",
    "Initials", "Suffix",
    "CollectiveName"
    ]
[simple_elem(e) for e in elements]
personal_name = LastName + \
                Opt(Alt(ForeName, FirstName + Opt(MiddleName))) + \
                Opt(Initials) + \
                Opt(Suffix)
author_name = Alt(personal_name, CollectiveName)


########################################
# Dates
elements = [
    "Year", "Month", "Day",
    "Season", "MedlineDate",
    "Hour", "Minute", "Second"
    ]
[simple_elem(e) for e in elements]
normal_date = Year + Month + Day + \
              Opt(Hour + Opt(Minute + Opt(Second)))
pub_date = Alt((Year + Opt(Alt((Month + Opt(Day)), Season))), MedlineDate)


simple_elem("CopyrightInformation")
simple_elem("AbstractText")
group_elem("Abstract", AbstractText + Opt(CopyrightInformation))

########################################
# NCBIArticle

simple_elem("NlmUniqueID")
simple_elem("PMID")
simple_elem("SubHeading", "MajorTopicYN")
simple_elem("QualifierName", "MajorTopicYN")
simple_elem("Descriptor", "MajorTopicYN")
simple_elem("DescriptorName", "MajorTopicYN")
group_elem("MeshHeading",
           Alt(DescriptorName, Descriptor) + \
           Alt(Rep(QualifierName), Rep(SubHeading)))
group_elem("MeshHeadingList", Rep1(MeshHeading))
simple_elem("MedlinePgn")
simple_elem("EndPage")
simple_elem("StartPage")
group_elem("Pagination",
           Alt(StartPage + Opt(EndPage) + Opt(MedlinePgn), MedlinePgn))

simple_elem("Affiliation")
group_elem("Author", author_name + Opt(Affiliation))
group_elem("AuthorList", Rep1(Author), "CompleteYN")
simple_elem("Language")
simple_elem("PublicationType")
group_elem("PublicationTypeList", Rep1(PublicationType))
simple_elem("Title")      # These were moved up, so that the definitions
simple_elem("Volume")     # will be before Book.
simple_elem("VernacularTitle")
simple_elem("CollectionTitle")
simple_elem("ArticleTitle")
simple_elem("Publisher")
group_elem("PubDate", pub_date)
group_elem("Book", PubDate + Publisher + Title +
           Opt(AuthorList) + Opt(CollectionTitle) + Opt(Volume))
simple_elem("Country")
simple_elem("MedlineTA")
simple_elem("MedlineCode")
group_elem("MedlineJournalInfo",
           Opt(Country) + MedlineTA + Opt(MedlineCode) + Opt(NlmUniqueID))
simple_elem("DateOfElectronicPublication")
simple_elem("ISOAbbreviation")
simple_elem("Coden")
simple_elem("Issue")
group_elem("JournalIssue", Opt(Volume) + Opt(Issue) + PubDate)
simple_elem("ISSN")
group_elem("Journal",
           Opt(ISSN) + \
           JournalIssue + \
           Opt(Coden) + \
           Opt(Title) + \
           Opt(ISOAbbreviation)
           )

simple_elem("GrantID")
simple_elem("Acronym")
simple_elem("Agency")
group_elem("Grant", Opt(GrantID) + Opt(Acronym) + Opt(Agency))
group_elem("GrantList", Rep1(Grant), "CompleteYN")
simple_elem("AccessionNumber")
group_elem("AccessionNumberList", Rep1(AccessionNumber))
simple_elem("DataBankName")
group_elem("DataBank", DataBankName + Opt(AccessionNumberList))
group_elem("DataBankList", Rep1(DataBank), "CompleteYN")

group_elem("Article",
           Alt(Journal, Book) + \
           ArticleTitle + \
           Pagination + \
           Opt(Abstract) + \
           Opt(Affiliation) + \
           Opt(AuthorList) + \
           Rep1(Language) + \
           Opt(DataBankList) + \
           Opt(GrantList) + \
           PublicationTypeList + \
           Opt(VernacularTitle) + \
           Opt(DateOfElectronicPublication)
           )
group_elem("NCBIArticle", PMID + Article + Opt(MedlineJournalInfo))






######################################################################
# Implement Martel expressions that recognize:                       #
# http://www.nlm.nih.gov/databases/dtd/nlmmedlinecitation_011101.dtd #
######################################################################


simple_elem("MedlineID")

simple_elem("Note")
simple_elem("RefSource")
Ref_template = RefSource + Opt(Alt(PMID, MedlineID)) + Opt(Note)


########################################
# MedlineCitation

group_elem("OriginalReportIn", Ref_template)
group_elem("SummaryForPatientsIn", Ref_template)
group_elem("CommentOn", Ref_template)
group_elem("CommentIn", Ref_template)
group_elem("ErratumIn", Ref_template)
group_elem("RepublishedFrom", Ref_template)
group_elem("RepublishedIn", Ref_template)
group_elem("RetractionOf", Ref_template)
group_elem("RetractionIn", Ref_template)
group_elem("UpdateIn", Ref_template)
group_elem("UpdateOf", Ref_template)
group_elem("CommentsCorrections",
           Rep(CommentOn) + Rep(CommentIn) + \
           Rep(ErratumIn) + \
           Rep(RepublishedFrom) + Rep(RepublishedIn) + \
           Rep(RetractionOf) + Rep(RetractionIn) + \
           Rep(UpdateIn) + Rep(UpdateOf) + \
           Rep(SummaryForPatientsIn) + Rep(OriginalReportIn)
           )
simple_elem("NumberOfReferences")
group_elem("PersonalNameSubject", personal_name)
group_elem("PersonalNameSubjectList", Rep1(PersonalNameSubject))
simple_elem("GeneSymbol")
group_elem("GeneSymbolList", Rep1(GeneSymbol))
simple_elem("NameOfSubstance")
simple_elem("CASRegistryNumber")
simple_elem("RegistryNumber")
group_elem("Chemical", Alt(CASRegistryNumber, RegistryNumber) + \
           NameOfSubstance)
group_elem("ChemicalList", Rep1(Chemical))
simple_elem("CitationSubset")
simple_elem("GeneralNote", "Owner")
group_elem("Investigator", personal_name + Opt(Affiliation))
group_elem("InvestigatorList", Rep1(Investigator))
simple_elem("OtherID", "Source")
simple_elem("SpaceFlightMission")
simple_elem("Keyword", "MajorTopicYN")
group_elem("KeywordList", Rep1(Keyword), "Owner")
group_elem("OtherAbstract",
           AbstractText + Opt(CopyrightInformation),
           "Type")
group_elem("DateRevised", normal_date)
group_elem("DateCompleted", normal_date)
group_elem("DateCreated", normal_date)
group_elem("MedlineCitation",
           Opt(MedlineID) + \
           Opt(PMID) + \
           DateCreated + \
           Opt(DateCompleted) + \
           Opt(DateRevised) + \
           Article + \
           MedlineJournalInfo + \
           Opt(ChemicalList) + \
           Rep(CitationSubset) + \
           Opt(CommentsCorrections) + \
           Opt(GeneSymbolList) + \
           Opt(MeshHeadingList) + \
           Opt(NumberOfReferences) + \
           Opt(PersonalNameSubjectList) + \
           Rep(OtherID) + \
           Rep(OtherAbstract) + \
           Rep(KeywordList) + \
           Rep(SpaceFlightMission) + \
           Opt(InvestigatorList) + \
           Rep(GeneralNote),
           "Owner", "Status"
           )





           
######################################################################
# Implement Martel expressions that recognize:                       #
# http://www.nlm.nih.gov/databases/dtd/nlmmedline_011101.dtd         #
######################################################################



# The DeleteCitation tags start with spaces, so I have to make a
# special case for it.
space = Any(" \t")
DeleteCitation_start = Rep(space) + Str("<DeleteCitation>") + AnyEol()
DeleteCitation_end = Rep(space) + Str("</DeleteCitation>") + AnyEol()

# The file doesn't always end in a newline, so make MedlineCitationSet
# end in an optional Eol.
MedlineCitationSet_end = Str("</MedlineCitationSet>") + Opt(AnyEol())


group_elem("DeleteCitation", Alt(Rep1(MedlineID), Rep1(PMID)))
group_elem("MedlineCitationSet", Rep(MedlineCitation) + Opt(DeleteCitation))





######################################################################
# Other stuff                                                        #
#                                                                    #
######################################################################


# Should match the proper dtd in here...
DOCTYPE = Str("<!DOCTYPE") + Re(r"[^>]+") + Str(">") + AnyEol()

citation_format = MedlineCitation

# I'm going to use a RecordReader so that I can parse one record at a
# time, instead of sucking the whole XML file into memory.  Each
# citation is going to be a record.  Thus, the header is everything
# before the first citation and the footer is everything after the
# last citation.
header_format = Group("header", DOCTYPE + MedlineCitationSet_start)
footer_format = Opt(DeleteCitation) + MedlineCitationSet_end
format = HeaderFooter(
    None, {},
    # Unfortunately, RecordReader.Until doesn't work because some
    # MedlineCitations have attributes are in the form
    # <MedlineCitation Owner="NLM">.  "<MedlineCitation" by itself
    # won't differentiate between the beginning of a
    # MedlineCitationSet or the beginning of a MedlineCitation.  Thus,
    # I'm just going to read the first 4 lines and hope that's the
    # whole header.
    #header_format, RecordReader.Until, ("<MedlineCitation>",),
    header_format, RecordReader.CountLines, (4,),
    citation_format, RecordReader.EndsWith, ("</MedlineCitation>",),
    footer_format, RecordReader.Everything, (),
    )
