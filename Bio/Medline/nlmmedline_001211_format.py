"""nlmmedline_xml_format.py

A Martel format to parse the NLM's XML format for Medline.

Formats:
citation_format
format

"""
# http://www.nlm.nih.gov/databases/dtd/nlmmedline_001211.dtd
# http://www.nlm.nih.gov/databases/dtd/nlmmedlinecitation_001211.dtd
# http://www.nlm.nih.gov/databases/dtd/nlmcommon_001211.dtd


import string
import sys

from Martel import *
from Martel import RecordReader

self = sys.modules[__name__]

xml_name = Re(r"[\w_:]+")


# Automatically create Expressions that recognize lines that consist
# of either a start or end element.  The Expression that recognizes a
# line with a start element will be named "NAME_start" and respective
# end one will be named "NAME_end".

# e.g.
# <Article>\n     recognized by Article_start
# </Article>\n    recognized by Article_end
separate_start_and_end_elements = [
    ("Abstract", []),
    ("AccessionNumberList", []),
    ("Article", []),
    ("Author", []),
    ("AuthorList", ["CompleteYN"]),
    ("Chemical", []),
    ("ChemicalList", []),
    ("CommentIn", []),
    ("CommentOn", []),
    ("CommentsCorrections", []),
    ("DataBank", []),
    ("DataBankList", ["CompleteYN"]),
    ("DateCompleted", []),
    ("DateCreated", []),
    ("DateRevised", []),
    ("ErratumIn", []),
    ("GeneSymbolList", []),
    ("Grant", []),
    ("GrantList", ["CompleteYN"]),
    ("Journal", []),
    ("JournalIssue", []),
    ("MedlineCitation", []),
    ("MedlineCitationSet", []),
    ("MedlineJournalInfo", []),
    ("MeshHeading", []),
    ("MeshHeadingList", []),
    ("Pagination", []),
    ("PersonalNameSubject", []),
    ("PersonalNameSubjectList", []),
    ("PubDate", []),
    ("PublicationTypeList", []),
    ("RepublishedFrom", []),
    ("RepublishedIn", []),
    ("RetractionIn", []),
    ("RetractionOf", []),
    ("UpdateIn", []),
    ("UpdateOf", []),
    ]

for element, attrs in separate_start_and_end_elements:
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
    start_elem = start + AnyEol()
    setattr(self, "%s_start" % element, start_elem)
    
    end_elem = Str("</%s>" % element) + AnyEol()
    setattr(self, "%s_end" % element, end_elem)

# The file doesn't always end in a newline, so make MedlineCitationSet
# end in an optional Eol.
MedlineCitationSet_end = Str("</MedlineCitationSet>") + Opt(AnyEol())


# Automatically create Expression that recognize lines that consist of
# both the start and end element.  The Expression will be named
# "NAME".
# e.g.
# <PMID>12345</PMID>\n   recognized by PMID
# The data will be in a Group with the same name as the element.

# element name, possible attributes
joined_start_and_end_elements = [
    ("AbstractText", []),
    ("AccessionNumber", []),
    ("Acronym", []),
    ("Affiliation", []),
    ("Agency", []),
    ("ArticleTitle", []),
    ("CASRegistryNumber", []),
    ("CitationSubset", []),
    ("CollectiveName", []),
    ("CopyrightInformation", []),
    ("Country", []),
    ("DataBankName", []),
    ("Day", []),
    ("Descriptor", ["MajorTopicYN"]),
    ("EndPage", []),
    ("FirstName", []),
    ("GeneSymbol", []),
    ("GrantID", []),
    ("ISSN", []),
    ("Initials", []),
    ("Issue", []),
    ("Language", []),
    ("LastName", []),
    ("MedlineCode", []),
    ("MedlineDate", []),
    ("MedlineID", []),
    ("MedlinePgn", []),
    ("MedlineTA", []),
    ("MiddleName", []),
    ("Month", []),
    ("NameOfSubstance", []),
    ("NlmUniqueID", []),
    ("Note", []),
    ("NumberOfReferences", []),
    ("PMID", []),
    ("PublicationType", []),
    ("RefSource", []),
    ("Season", []),
    ("StartPage", []),
    ("SubHeading", ["MajorTopicYN"]),
    ("Suffix", []),
    ("VernacularTitle", []),
    ("Volume", []),
    ("Year", []),
    ]
for element, attrs in joined_start_and_end_elements:
    group_name = element
    group_expression = Re(r"[^<]+")

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
    end = Str("</%s>" % element)
    
    expr = start + \
           Group(group_name, group_expression) + \
           end + \
           AnyEol()
    setattr(self, element, expr)


# Create groups of expressions.  A group consists of the start and end
# elements with an expression in-between.  The Expression for the
# group will be called "NAME_group".

def make_group(element, expr):
    group_expr = getattr(self, "%s_start" % element) + \
                 expr + \
                 getattr(self, "%s_end" % element)
    group_expr = Group(element, group_expr)
    setattr(self, "%s_group" % element, group_expr)

pub_date = Alt((Year + Opt(Alt((Month + Opt(Day)), Season))), MedlineDate)
make_group("DateCreated", pub_date)
make_group("DateCompleted", pub_date)
make_group("DateRevised", pub_date)
make_group("PubDate", pub_date)
make_group("JournalIssue", (Opt(Volume) + \
                            Opt(Issue) + \
                            PubDate_group
                            ))
make_group("Pagination",
           Alt((StartPage + Opt(EndPage) + Opt(MedlinePgn)), MedlinePgn))
make_group("Abstract", AbstractText + Opt(CopyrightInformation))
personal_name = LastName + \
                Opt(FirstName + Opt(MiddleName)) + \
                Opt(Initials) + \
                Opt(Suffix)
author_name = Alt(personal_name, CollectiveName)
make_group("Author", author_name)
make_group("AuthorList", Rep1(Author_group))
make_group("PublicationTypeList", Rep1(PublicationType))
make_group("Journal", (Opt(ISSN) + \
                       JournalIssue_group
                       # Coden?
                       # Title?
                       # ISOAbbreviation?
                       ))
make_group("AccessionNumberList", Rep1(AccessionNumber))
make_group("DataBank", DataBankName + Opt(AccessionNumberList_group))
make_group("DataBankList", Rep1(DataBank_group))
make_group("Grant", (GrantID + \
                     Opt(Acronym) + \
                     Opt(Agency)
                     ))
make_group("GrantList", Rep1(Grant_group))
make_group("Article", (Journal_group + \
                       ArticleTitle + \
                       Pagination_group + \
                       Opt(Abstract_group) + \
                       Opt(Affiliation) + \
                       Opt(AuthorList_group) + \
                       Rep1(Language) + \
                       Opt(DataBankList_group) + \
                       Opt(GrantList_group) + \
                       PublicationTypeList_group + \
                       Opt(VernacularTitle)
                       # DateOfElectronicPublication?
                       ))

make_group("MedlineJournalInfo", (Country + \
                                  MedlineTA + \
                                  MedlineCode + \
                                  Opt(NlmUniqueID)
                                  ))
make_group("Chemical", (CASRegistryNumber + NameOfSubstance))
make_group("ChemicalList", Rep1(Chemical_group))

Ref_template = RefSource + \
               Opt(MedlineID) + \
               Opt(Note)
make_group("CommentOn", Ref_template)
make_group("CommentIn", Ref_template)
make_group("ErratumIn", Ref_template)
make_group("RepublishedFrom", Ref_template)
make_group("RepublishedIn", Ref_template)
make_group("RetractionOf", Ref_template)
make_group("RetractionIn", Ref_template)
make_group("UpdateIn", Ref_template)
make_group("UpdateOf", Ref_template)
make_group("CommentsCorrections", (Rep(CommentOn_group) + \
                                   Rep(CommentIn_group) + \
                                   Rep(ErratumIn_group) + \
                                   Rep(RepublishedFrom_group) + \
                                   Rep(RepublishedIn_group) + \
                                   Rep(RetractionOf_group) + \
                                   Rep(RetractionIn_group) + \
                                   Rep(UpdateIn_group) + \
                                   Rep(UpdateOf_group)
                                   ))
make_group("MeshHeading", Descriptor + Rep(SubHeading))
make_group("MeshHeadingList", Rep1(MeshHeading_group))
make_group("PersonalNameSubject", personal_name)
make_group("PersonalNameSubjectList", Rep1(PersonalNameSubject_group))
make_group("GeneSymbolList", Rep1(GeneSymbol))
make_group("MedlineCitation", (MedlineID + \
                               Opt(PMID) + \
                               DateCreated_group + \
                               Opt(DateCompleted_group) + \
                               Opt(DateRevised_group) + \
                               Article_group + \
                               MedlineJournalInfo_group + \
                               # AdditionalInformation?
                               Opt(ChemicalList_group) + \
                               Rep(CitationSubset) + \
                               Opt(CommentsCorrections_group) + \
                               Opt(GeneSymbolList_group) + \
                               Opt(MeshHeadingList_group) + \
                               Opt(NumberOfReferences) + \
                               Opt(PersonalNameSubjectList_group)
                               ))

# <!ELEMENT AdditionalInformation (OtherAbstract*,
#                                  Keyword*, 
#                                  ProcurementSource*,
#                                  SponsoringAgency*,
#                                  SpaceFlightMission*)>
# <!ELEMENT OtherAbstract (Abstract, AbstractAuthor)>



    
DOCTYPE = Str("<!DOCTYPE") + Re(r"[^>]+") + Str(">") + AnyEol()

citation_format = MedlineCitation_group

header_format = Group("header", DOCTYPE + \
                MedlineCitationSet_start)
footer_format = MedlineCitationSet_end

format = HeaderFooter(
    None,
    header_format, RecordReader.Until, ("<MedlineCitation>",),
    citation_format, RecordReader.EndsWith, ("</MedlineCitation>",),
    footer_format, RecordReader.Everything, (),
    )
