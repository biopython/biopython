# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch
# from the "journals" database. No DTD file seems to be associated with
# these XML files.
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if name=="SerialSet":
        object = []
        self.path = []
        self.record = object
    else:
        if name in ("IndexingHistoryList",
                    "BroadJournalHeadingList",
                    "CrossReferenceList"):
            object = []
        elif name in ("PublicationInfo",
                      "DateOfAction",
                      "IlsCreatedTimestamp",
                      "IlsUpdatedTimestamp"):
            object = {}
        elif name in ("Serial",
                      "IndexingHistory",
                      "CurrentlyIndexedForSubset"):
            object = {}
            keys = attrs.keys()
            for key in keys:
                object[str(key)] = str(attrs[key])
        elif name=="CrossReference":
            object = [str(attrs["XrType"]), None]
        elif name=="ISSN":
            object = [str(attrs["IssnType"]), None]
        else:
            object = ""
        if object!="":
            current = self.path[-1]
            if type(current)==list:
                current.append(object)
            elif type(current)==dict:
                current[name] = object
    self.path.append(object)

def endElement(self, name):
    self.path.pop()
    if name=="Language":
        if not "Language" in self.path[-1]:
            self.path[-1]["Language"] = []
        self.path[-1]["Language"].append(self.content)
    elif name=="BroadJournalHeading":
        self.path[-1].append(self.content)
    elif name=="XrTitle":
        self.path[-1][1] = self.content
    elif name=="ISSN":
        self.path[-1][name][1] = self.content
    elif name in ("SortSerialName",
                  "Year",
                  "Month",
                  "Day",
                  "IndexingSubset",
                  "IndexOnlineYN",
                  "CurrentlyIndexedYN",
                  "Coverage",
                  "MinorTitleChangeYN",
                  "Coden",
                  "AcidFreeYN",
                  "ContinuationNotes",
                  "NlmUniqueID",
                  "Title",
                  "MedlineTA",
                  "Country",
                  "Place",
                  "Publisher",
                  "PublicationFirstYear",
                  "PublicationEndYear",
                  "Frequency",
                  "ISOAbbreviation"):
        self.path[-1][name] = self.content
