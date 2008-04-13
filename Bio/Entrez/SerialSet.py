# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch
# from the "journals" database. No DTD file seems to be associated with
# these XML files.
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["SerialSet"]:
        self.record = []
    elif self.element==["SerialSet", "Serial"]:
        d = {}
        if "DataCreationMethod" in attrs:
            d["DataCreationMethod"] = str(attrs["DataCreationMethod"])
        self.record.append(d)
    elif self.element==["SerialSet", "Serial", "PublicationInfo"]:
        self.record[-1]["PublicationInfo"] = {}
    elif self.element==["SerialSet", "Serial", "ISSN"]:
        self.record[-1]["ISSN"] = [str(attrs["IssnType"]), None]
    elif self.element==["SerialSet", "Serial", "IndexingHistoryList"]:
        self.record[-1]["IndexingHistoryList"] = []
    elif self.element==["SerialSet", "Serial", "IndexingHistoryList", "IndexingHistory"]:
        d = {}
        for key in attrs.keys():
            d[str(key)] = str(attrs.getValue(key))
        self.record[-1]["IndexingHistoryList"].append(d)
    elif self.element==["SerialSet", "Serial", "IndexingHistoryList", "IndexingHistory", "DateOfAction"]:
        self.record[-1]["IndexingHistoryList"][-1]["DateOfAction"] = {}
    elif self.element==["SerialSet", "Serial", "CurrentlyIndexedForSubset"]:
        d = {}
        for key in attrs.keys():
            d[str(key)] = str(attrs.getValue(key))
        self.record[-1]["CurrentlyIndexedForSubset"] = d
    elif self.element==["SerialSet", "Serial", "BroadJournalHeadingList"]:
        self.record[-1]["BroadJournalHeadingList"] = []
    elif self.element==["SerialSet", "Serial", "CrossReferenceList"]:
        self.record[-1]["CrossReferenceList"] = []
    elif self.element==["SerialSet", "Serial", "CrossReferenceList", "CrossReference"]:
        self.record[-1]["CrossReferenceList"].append([str(attrs["XrType"]), None])
    elif self.element==["SerialSet", "Serial", "IlsCreatedTimestamp"]:
        self.record[-1]["IlsCreatedTimestamp"] = {}
    elif self.element==["SerialSet", "Serial", "IlsUpdatedTimestamp"]:
        self.record[-1]["IlsUpdatedTimestamp"] = {}

def endElement(self, name):
    if self.element==["SerialSet"]:
        self.record = []
    elif self.element==["SerialSet", "Serial", "NlmUniqueID"]:
        self.record[-1]["NlmUniqueID"] = self.content
    elif self.element==["SerialSet", "Serial", "Title"]:
        self.record[-1]["Title"] = self.content
    elif self.element==["SerialSet", "Serial", "MedlineTA"]:
        self.record[-1]["MedlineTA"] = self.content
    elif self.element==["SerialSet", "Serial", "PublicationInfo", "Country"]:
        self.record[-1]["PublicationInfo"]["Country"] = self.content
    elif self.element==["SerialSet", "Serial", "PublicationInfo", "Place"]:
        self.record[-1]["PublicationInfo"]["Place"] = self.content
    elif self.element==["SerialSet", "Serial", "PublicationInfo", "Publisher"]:
        self.record[-1]["PublicationInfo"]["Publisher"] = self.content
    elif self.element==["SerialSet", "Serial", "PublicationInfo", "PublicationFirstYear"]:
        self.record[-1]["PublicationInfo"]["PublicationFirstYear"] = self.content
    elif self.element==["SerialSet", "Serial", "PublicationInfo", "PublicationEndYear"]:
        self.record[-1]["PublicationInfo"]["PublicationEndYear"] = self.content
    elif self.element==["SerialSet", "Serial", "PublicationInfo", "Frequency"]:
        self.record[-1]["PublicationInfo"]["Frequency"] = self.content
    elif self.element==["SerialSet", "Serial", "ISSN"]:
        self.record[-1]["ISSN"][1] = self.content
    elif self.element==["SerialSet", "Serial", "ISOAbbreviation"]:
        self.record[-1]["ISOAbbreviation"] = self.content
    elif self.element==["SerialSet", "Serial", "Language"]:
        if not "Language" in self.record[-1]:
            self.record[-1]["Language"] = []
        self.record[-1]["Language"].append(self.content)
    elif self.element==["SerialSet", "Serial", "ContinuationNotes"]:
        self.record[-1]["ContinuationNotes"] = self.content
    elif self.element==["SerialSet", "Serial", "AcidFreeYN"]:
        self.record[-1]["AcidFreeYN"] = self.content
    elif self.element==["SerialSet", "Serial", "Coden"]:
        self.record[-1]["Coden"] = self.content
    elif self.element==["SerialSet", "Serial", "MinorTitleChangeYN"]:
        self.record[-1]["MinorTitleChangeYN"] = self.content
    elif self.element==["SerialSet", "Serial", "IndexingHistoryList", "IndexingHistory", "DateOfAction", "Year"]:
        self.record[-1]["IndexingHistoryList"][-1]["DateOfAction"]["Year"] = self.content
    elif self.element==["SerialSet", "Serial", "IndexingHistoryList", "IndexingHistory", "DateOfAction", "Month"]:
        self.record[-1]["IndexingHistoryList"][-1]["DateOfAction"]["Month"] = self.content
    elif self.element==["SerialSet", "Serial", "IndexingHistoryList", "IndexingHistory", "DateOfAction", "Day"]:
        self.record[-1]["IndexingHistoryList"][-1]["DateOfAction"]["Day"] = self.content
    elif self.element==["SerialSet", "Serial", "IndexingHistoryList", "IndexingHistory", "Coverage"]:
        self.record[-1]["IndexingHistoryList"][-1]["Coverage"] = self.content
    elif self.element==["SerialSet", "Serial", "CurrentlyIndexedYN"]:
        self.record[-1]["CurrentlyIndexedYN"] = self.content
    elif self.element==["SerialSet", "Serial", "IndexOnlineYN"]:
        self.record[-1]["IndexOnlineYN"] = self.content
    elif self.element==["SerialSet", "Serial", "IndexingSubset"]:
        self.record[-1]["IndexingSubset"] = self.content
    elif self.element==["SerialSet", "Serial", "BroadJournalHeadingList", "BroadJournalHeading"]:
        self.record[-1]["BroadJournalHeadingList"].append(self.content)
    elif self.element==["SerialSet", "Serial", "CrossReferenceList", "CrossReference", "XrTitle"]:
        self.record[-1]["CrossReferenceList"][-1][1] = self.content
    elif self.element==["SerialSet", "Serial", "SortSerialName"]:
        self.record[-1]["SortSerialName"] = self.content
    elif self.element==["SerialSet", "Serial", "IlsCreatedTimestamp", "Year"]:
        self.record[-1]["IlsCreatedTimestamp"]["Year"] = self.content
    elif self.element==["SerialSet", "Serial", "IlsCreatedTimestamp", "Month"]:
        self.record[-1]["IlsCreatedTimestamp"]["Month"] = self.content
    elif self.element==["SerialSet", "Serial", "IlsCreatedTimestamp", "Day"]:
        self.record[-1]["IlsCreatedTimestamp"]["Day"] = self.content
    elif self.element==["SerialSet", "Serial", "IlsUpdatedTimestamp", "Year"]:
        self.record[-1]["IlsUpdatedTimestamp"]["Year"] = self.content
    elif self.element==["SerialSet", "Serial", "IlsUpdatedTimestamp", "Month"]:
        self.record[-1]["IlsUpdatedTimestamp"]["Month"] = self.content
    elif self.element==["SerialSet", "Serial", "IlsUpdatedTimestamp", "Day"]:
        self.record[-1]["IlsUpdatedTimestamp"]["Day"] = self.content
