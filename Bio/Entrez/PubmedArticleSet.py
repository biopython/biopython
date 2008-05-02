# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch
# from the PubMed database, as specified by NCBI's DTD file
# pubmed_080101.dtd (2007-11-30 16:19:51)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.


"""
This code is used to parse XML results returned by Entrez's PubMed
Article database.

The code is not meant to be used by itself, but by calling
Bio.Entrez.read() instead.  When parsing a PubMed Article XML file
in this way, the record returned is a list of dictionaries.

For example,

from Bio import Entrez
handle = Entrez.efetch(db="pubmed", id="17238260", retmode="XML")
articles = Entrez.read(handle)
assert len(articles)==1
"""


class DescriptorName(str):
    """Corresponds to one <DescriptorName> block.
       The attribute "MajorTopic" is stored as .attributes["MajorTopic"]
    """
    tag = "DescriptorName"
    def __init__(self, s):
        str.__init__(self, s)
        self.attributes = {}


class QualifierName(str):
    """Corresponds to one <QualifierName> block.
       The attribute "MajorTopic" is stored as .attributes["MajorTopic"]
    """
    tag = "QualifierName"
    def __init__(self, s):
        str.__init__(self, s)
        self.attributes = {}


def startElement(self, name, attrs):
    assert self.element
    if self.element==["PubmedArticleSet"]:
        self.record = []
    elif self.element==["PubmedArticleSet", "PubmedArticle"]:
        self.record.append({"Journal":None,
                            "Volume":None,
                            "Issue":None,
                            "PublicationStatus":None,
                            })
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation"]:
        self.record[-1]["MedlineCitation"] = {}
        if "Owner" in attrs:
            self.record[-1]["MedlineCitation"]["Owner"] = str(attrs["Owner"])
        if "Status" in attrs:
            self.record[-1]["MedlineCitation"]["Status"] = str(attrs["Status"])
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "DateCreated"]:
        self.record[-1]["MedlineCitation"]["DateCreated"] = {}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "DateCompleted"]:
        self.record[-1]["MedlineCitation"]["DateCompleted"] = {}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "DateRevised"]:
        self.record[-1]["MedlineCitation"]["DateRevised"] = {}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article"]:
        self.record[-1]["MedlineCitation"]["Article"] = {}
        if "PubModel" in attrs:
            self.record[-1]["MedlineCitation"]["Article"]["PubModel"] = str(attrs["PubModel"])
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "Journal"]:
        self.record[-1]["MedlineCitation"]["Article"]["Journal"] = {}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "Journal","ISSN"]:
        IssnType = str(attrs["IssnType"])
        self.record[-1]["MedlineCitation"]["Article"]["Journal"]["ISSN"] = [IssnType, None]
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "Journal","JournalIssue"]:
        self.record[-1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"] = {"CitedMedium": str(attrs["CitedMedium"])}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "Journal","JournalIssue", "PubDate"]:
        self.record[-1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"] = {}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "Pagination"] :
        self.record[-1]["MedlineCitation"]["Article"]["Pagination"] = {}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "AuthorList"]:
        CompleteYN = str(attrs["CompleteYN"])
        self.record[-1]["MedlineCitation"]["Article"]["AuthorList"] = [CompleteYN]
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "AuthorList", "Author"]:
        ValidYN = str(attrs["ValidYN"])
        self.record[-1]["MedlineCitation"]["Article"]["AuthorList"].append({"ValidYN": ValidYN})
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "Language"]:
        self.record[-1]["MedlineCitation"]["Article"]["Language"] = []
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "PublicationTypeList"]:
        self.record[-1]["MedlineCitation"]["Article"]["PublicationTypeList"] = []
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "MeshHeadingList"]:
        self.record[-1]["MedlineCitation"]["MeshHeadingList"] = []
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "MeshHeadingList", "MeshHeading"]:
        self.record[-1]["MedlineCitation"]["MeshHeadingList"].append([])
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "MeshHeadingList", "MeshHeading", "DescriptorName"]:
        self._attributes = attrs
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "MeshHeadingList", "MeshHeading", "QualifierName"]:
        self._attributes = attrs
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "OtherID"]:
        if not "OtherID" in self.record[-1]["MedlineCitation"]:
            self.record[-1]["MedlineCitation"]["OtherID"] = []
        value = [str(attrs["Source"]), None]
        self.record[-1]["MedlineCitation"]["OtherID"].append(value)
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "KeywordList"]:
        value = str(attrs["Owner"])
        self.record[-1]["MedlineCitation"]["KeywordList"] = [value]
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "KeywordList", "Keyword"]:
        value = [str(attrs["MajorTopicYN"]), None]
        self.record[-1]["MedlineCitation"]["KeywordList"].append(value)
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "GeneralNote"]:
        if not "GeneralNote" in self.record[-1]["MedlineCitation"]:
            self.record[-1]["MedlineCitation"]["GeneralNote"] = []
        value = [str(attrs["Owner"]), None]
        self.record[-1]["MedlineCitation"]["GeneralNote"].append(value)
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "ChemicalList"]:
        self.record[-1]["MedlineCitation"]["ChemicalList"] = []
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "ChemicalList", "Chemical"]:
        self.record[-1]["MedlineCitation"]["ChemicalList"].append({})
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "MedlineJournalInfo"]:
       self.record[-1]["MedlineCitation"]["MedlineJournalInfo"] = {}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "PubmedData"] :
       self.record[-1]["PubmedData"] = {}
    elif self.element==["PubmedArticleSet", "PubmedArticle", "PubmedData", "ArticleIdList"]:
       self.record[-1]["PubmedData"]["ArticleIdList"] = []
    elif self.element==["PubmedArticleSet", "PubmedArticle", "PubmedData", "ArticleIdList", "ArticleId"]:
        value = [str(attrs["IdType"]), None]
        self.record[-1]["PubmedData"]["ArticleIdList"].append(value)
    elif self.element==["PubmedArticleSet", "PubmedArticle", "PubmedData", "History"]:
        self.record[-1]["PubmedData"]["History"] = []
    elif self.element==["PubmedArticleSet", "PubmedArticle", "PubmedData", "History", "PubMedPubDate"]:
        key = "PubStatus"
        value = str(attrs[key])
        d = {key: value}
        self.record[-1]["PubmedData"]["History"].append(d)
    elif self.element==["PubmedArticleSet", "PubmedArticle", "MedlineCitation", "Article", "Abstract"]:
        self.record[-1]["MedlineCitation"]["Article"]["Abstract"] = {}

def endElement(self, name):
    if name=="Error":
        error = self.content
        raise RuntimeError(error)
    elif self.element[:3]==["PubmedArticleSet", "PubmedArticle", "MedlineCitation"] :
        if self.element[3:]==["PMID"]:
            self.record[-1]["MedlineCitation"]["PMID"] = self.content
        elif self.element[3:]==["CitationSubset"]:
            self.record[-1]["MedlineCitation"]["CitationSubset"] = self.content
        elif self.element[3:5]==["ChemicalList","Chemical"] :
            if len(self.element)==6 :
                self.record[-1]["MedlineCitation"]["ChemicalList"][-1][str(self.element[-1])] = self.content
        elif len(self.element)==5  and str(self.element[3]) in ["DateCreated", "DateCompleted", "DateRevised"] :
            self.record[-1]["MedlineCitation"][str(self.element[3])][str(self.element[4])] = self.content
        elif self.element[3:]==["MedlineJournalInfo", "Country"] :
            self.record[-1]["MedlineCitation"]["MedlineJournalInfo"]["Country"] = self.content
        elif self.element[3:]==["MedlineJournalInfo", "MedlineTA"] :
            self.record[-1]["MedlineCitation"]["MedlineJournalInfo"]["MedlineTA"] = self.content
        elif self.element[3:]==["MedlineJournalInfo", "NlmUniqueID"] :
            self.record[-1]["MedlineCitation"]["MedlineJournalInfo"]["NlmUniqueID"] = self.content
        elif self.element[3:]==["MeshHeadingList", "MeshHeading", "DescriptorName"]:
            d = DescriptorName(self.content)
            d.attributes["MajorTopicYN"] = str(self._attributes["MajorTopicYN"])
            del self._attributes
            self.record[-1]["MedlineCitation"]["MeshHeadingList"][-1].append(d)
        elif self.element[3:]==["MeshHeadingList", "MeshHeading", "QualifierName"]:
            q = QualifierName(self.content)
            q.attributes["MajorTopicYN"] = str(self._attributes["MajorTopicYN"])
            del self._attributes
            self.record[-1]["MedlineCitation"]["MeshHeadingList"][-1].append(q)
        elif self.element[3:]==["NumberOfReferences"]:
            self.record[-1]["MedlineCitation"]["NumberOfReferences"] = int(self.content)
        elif self.element[3:]==["OtherID"]:
            self.record[-1]["MedlineCitation"]["OtherID"][-1][1] = self.content
        elif self.element[3:]==["GeneralNote"]:
            self.record[-1]["MedlineCitation"]["GeneralNote"][-1][1] = self.content
        elif self.element[3:]==["KeywordList", "Keyword"]:
            self.record[-1]["MedlineCitation"]["KeywordList"][-1][1] = self.content
        elif self.element[3:4]==["Article"] :
            if self.element[4:]==["PublicationTypeList", "PublicationType"]:
                self.record[-1]["MedlineCitation"]["Article"]["PublicationTypeList"].append(self.content)
            elif len(self.element)==7 and self.element[4:6]==["AuthorList","Author"] :
                self.record[-1]["MedlineCitation"]["Article"]["AuthorList"][-1][str(name)] = self.content
            elif self.element[4:]==["Journal","ISSN"]:
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["ISSN"][-1] = self.content
            elif self.element[4:]==["Journal","JournalIssue","Volume"] :
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Volume"] = self.content
            elif self.element[4:]==["Journal","JournalIssue","Issue"] :
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["Issue"] = self.content
            elif self.element[4:]==["Journal","JournalIssue","PubDate","Year"]:
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"] = self.content
            elif self.element[4:]==["Journal","JournalIssue","PubDate","Season"]:
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Season"] = self.content
            elif self.element[4:]==["Journal","JournalIssue","PubDate","Month"]:
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Month"] = self.content
            elif self.element[4:]==["Journal","JournalIssue","PubDate","Day"]:
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Day"] = self.content
            elif self.element[4:]==["Journal","Title"] :
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["Title"] = self.content
            elif self.element[4:]==["Journal","ISOAbbreviation"] :
                self.record[-1]["MedlineCitation"]["Article"]["Journal"]["ISOAbbreviation"] = self.content
            elif self.element[4:]==["Abstract","AbstractText"] :
                self.record[-1]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"] = self.content
            elif self.element[4:]==["Abstract","CopyrightInformation"] :
                self.record[-1]["MedlineCitation"]["Article"]["Abstract"]["CopyrightInformation"] = self.content
            elif self.element[4:]==["Pagination", "MedlinePgn"] :
                self.record[-1]["MedlineCitation"]["Article"]["Pagination"]["MedlinePgn"] = self.content
            elif len(self.element) == 5 and self.element[-1]=="Language":
                self.record[-1]["MedlineCitation"]["Article"]["Language"].append(self.content)
            elif self.element[4:]==["ArticleTitle"]:
                self.record[-1]["MedlineCitation"]["Article"]["ArticleTitle"] = self.content
            elif len(self.element) == 5 and self.element[-1]=="Affiliation":
                self.record[-1]["MedlineCitation"]["Article"]["Affiliation"] = self.content
    elif self.element[:3]==["PubmedArticleSet", "PubmedArticle", "PubmedData"] :
        if self.element[3:]==["PublicationStatus"]:
            self.record[-1]["PubmedData"]["PublicationStatus"] = self.content
        elif self.element[3:]==["ArticleIdList", "ArticleId"]:
            self.record[-1]["PubmedData"]["ArticleIdList"][-1][1] = self.content
        elif self.element[3:]==["History", "PubMedPubDate", "Year"]:
            self.record[-1]["PubmedData"]["History"][-1]["Year"] = self.content
        elif self.element[3:]==["History", "PubMedPubDate", "Month"]:
            self.record[-1]["PubmedData"]["History"][-1]["Month"] = self.content
        elif self.element[3:]==["History", "PubMedPubDate", "Day"]:
            self.record[-1]["PubmedData"]["History"][-1]["Day"] = self.content
        elif self.element[3:]==["History", "PubMedPubDate", "Hour"]:
            self.record[-1]["PubmedData"]["History"][-1]["Hour"] = self.content
        elif self.element[3:]==["History", "PubMedPubDate", "Minute"]:
            self.record[-1]["PubmedData"]["History"][-1]["Minute"] = self.content
    self.content = ""

if __name__ == "__main__" :
    print "Quick example/test"
    from Bio import Entrez
    handle = Entrez.efetch(db="pubmed", id=["14630660","17238260"], retmode="XML")
    articles = Entrez.read(handle)
    assert len(articles)==2
    for article in articles :
        print
        for key,value in article.iteritems() :
            if isinstance(value,list) :
                print "%s - %s" % (key,value[0])
                for entry in value[1:] :
                    print " "*len(key) + " - %s" % entry
            else :
                print "%s = %s" % (key,value)
