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
    if name=="PubmedArticleSet":
        object = []
        self.record = object
        self.path = []
    elif name in ("DateCreated",
                  "DateCompleted",
                  "DateRevised",
                  "Journal",
                  "PubDate",
                  "Pagination",
                  "MedlineJournalInfo",
                  "PubmedData" ,
                  "Abstract"):
        object = {}
        self.path[-1][name] = object
    elif name in ("Language",
                  "PublicationTypeList",
                  "MeshHeadingList",
                  "ArticleIdList",
                  "ChemicalList",
                  "History"):
        object = []
        self.path[-1][name] = object
    elif name=="Chemical":
        object = {}
        self.path[-1].append(object)
    elif name=="MeshHeading":
        object = []
        self.path[-1].append(object)
    elif name in ("MedlineCitation",
                  "Article",
                  "JournalIssue",
                  "Author",
                  "PubMedPubDate"):
        object = {}
        keys = attrs.keys()
        for key in keys:
            object[key] = str(attrs[key])
        if type(self.path[-1])==list:
            self.path[-1].append(object)
        elif type(self.path[-1])==dict:
            self.path[-1][name] = object
    elif name=="PubmedArticle":
        object = {"Journal":None,
                  "Volume":None,
                  "Issue":None,
                  "PublicationStatus":None,
                 }
        self.path[-1].append(object)
    elif name=="ISSN":
        IssnType = str(attrs["IssnType"])
        object = [IssnType, None]
        self.path[-1][name] = object
    elif name=="AuthorList":
        CompleteYN = str(attrs["CompleteYN"])
        object = [CompleteYN]
        self.path[-1][name] = object
    elif name=="Keyword":
        object = [str(attrs["MajorTopicYN"]), None]
        self.path[-1].append(object)
    elif name in ("DescriptorName", "QualifierName"):
        self._attributes = attrs
        object = ""
    elif name=="OtherID":
        if not "OtherID" in self.path[-1]:
            self.path[-1]["OtherID"] = []
        object = [str(attrs["Source"]), None]
        self.path[-1][name].append(object)
    elif name=="KeywordList":
        object = [str(attrs["Owner"])]
        self.path[-1][name] = object
    elif name=="GeneralNote":
        if not "GeneralNote" in self.path[-1]:
            self.path[-1]["GeneralNote"] = []
        object = [str(attrs["Owner"]), None]
        self.path[-1][name].append(object)
    elif name=="ArticleId":
        object = [str(attrs["IdType"]), None]
        self.path[-1].append(object)
    else:
        object = ""
    self.path.append(object)

def endElement(self, name):
    object = self.path.pop()
    if name=="Error":
        error = self.content
        raise RuntimeError(error)
    elif name in ("ArticleId",
                  "GeneralNote",
                  "Keyword",
                  "OtherID"):
        object[1] = self.content
    elif name=="ISSN":
        self.path[-1][name][-1] = self.content
    elif name=="PublicationType":
        self.path[-1].append(self.content)
    elif name=="Language":
        self.path[-1][name].append(self.content)
    elif name in ("AbstractText",
                  "Affiliation",
                  "ArticleTitle",
                  "CitationSubset",
                  "CollectiveName",
                  "CopyrightInformation",
                  "Country",
                  "DatesAssociatedWithName",
                  "Day",
                  "FirstName",
                  "ForeName",
                  "Hour",
                  "Initials",
                  "ISOAbbreviation",
                  "Issue",
                  "LastName",
                  "MedlinePgn",
                  "MedlineTA",
                  "MiddleName",
                  "Minute",
                  "Month",
                  "NameOfSubstance",
                  "NameQualifier",
                  "NlmUniqueID",
                  "OtherInformation",
                  "PMID",
                  "PublicationStatus",
                  "RegistryNumber",
                  "Season",
                  "Second",
                  "Suffix",
                  "Title",
                  "TitleAssociatedWithName",
                  "Volume",
                  "Year"):
        self.path[-1][name] = self.content
    elif name=="NumberOfReferences":
        self.path[-1][name] = int(self.content)
    elif name=="DescriptorName":
        d = DescriptorName(self.content)
        d.attributes["MajorTopicYN"] = str(self._attributes["MajorTopicYN"])
        del self._attributes
        self.path[-1].append(d)
    elif name=="QualifierName":
        q = QualifierName(self.content)
        q.attributes["MajorTopicYN"] = str(self._attributes["MajorTopicYN"])
        del self._attributes
        self.path[-1].append(q)
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
