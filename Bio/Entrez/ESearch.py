# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eSearch.
# as specified by NCBI's DTD file eSearch_020511.dtd (2006-06-28 17:35:21)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.



def startElement(self, name, attrs):
    if self.element==["eSearchResult"]:
        self.record = {}
    elif self.element==["eSearchResult", "ErrorList"]:
        self.record["ErrorList"] = {"PhraseNotFound": [],
                                    "FieldNotFound": []}
    elif self.element==["eSearchResult", "WarningList"]:
        self.record["WarningList"] = {"PhraseIgnored": [],
                                      "QuotedPhraseNotFound": [],
                                      "OutputMessage": []}
    elif self.element==["eSearchResult", "IdList"]:
        self.record["IdList"] = []
    elif self.element==["eSearchResult", "TranslationSet"]:
        self.record["TranslationSet"] = []
    elif self.element==["eSearchResult", "TranslationSet", "Translation"]:
        translation = {}
        self.record["TranslationSet"].append(translation)
    elif self.element==["eSearchResult", "TranslationStack"]:
        self.record["TranslationStack"] = []
    elif self.element==["eSearchResult", "TranslationStack", "TermSet"]:
        termset = {}
        self.record["TranslationStack"].append(termset)

def endElement(self, name):
    if self.element==["eSearchResult", "Count"]:
        self.record["Count"] = int(self.content)
    elif self.element==["eSearchResult", "RetMax"]:
        self.record["RetMax"] = int(self.content)
    elif self.element==["eSearchResult", "RetStart"]:
        self.record["RetStart"] = int(self.content)
    elif self.element==["eSearchResult", "QueryKey"]:
        self.record["QueryKey"] = self.content
    elif self.element==["eSearchResult", "QueryTranslation"]:
        self.record["QueryTranslation"] = self.content
    elif self.element==["eSearchResult", "WebEnv"]:
        self.record["WebEnv"] = self.content
    elif self.element==["eSearchResult", "IdList", "Id"]:
        self.record["IdList"].append(self.content)
    elif self.element==["eSearchResult", "TranslationSet", "Translation", "From"]:
        translation = self.record["TranslationSet"][-1]
        translation["From"] = self.content
    elif self.element==["eSearchResult", "TranslationSet", "Translation", "To"]:
        translation = self.record["TranslationSet"][-1]
        translation["To"] = self.content
    elif self.element==["eSearchResult", "TranslationStack", "OP"]:
        self.record["TranslationStack"].append(self.content)
    elif self.element==["eSearchResult", "TranslationStack", "TermSet", "Term"]:
        termset = self.record["TranslationStack"][-1]
        termset["Term"] = self.content
    elif self.element==["eSearchResult", "TranslationStack", "TermSet", "Field"]:
        termset = self.record["TranslationStack"][-1]
        termset["Field"] = self.content
    elif self.element==["eSearchResult", "TranslationStack", "TermSet", "Count"]:
        termset = self.record["TranslationStack"][-1]
        termset["Count"] = int(self.content)
    elif self.element==["eSearchResult", "TranslationStack", "TermSet", "Explode"]:
        termset = self.record["TranslationStack"][-1]
        if self.content=='Y': termset["Explode"] = True
        elif self.content=='N': termset["Explode"] = False
    elif self.element==["eSearchResult", "ErrorList", "PhraseNotFound"]:
        self.record["ErrorList"]["PhraseNotFound"].append(self.content)
    elif self.element==["eSearchResult", "ErrorList", "FieldNotFound"]:
        self.record["ErrorList"]["FieldNotFound"].append(self.content)
    elif self.element==["eSearchResult", "WarningList", "PhraseIgnored"]:
        self.record["WarningList"]["PhraseIgnored"].append(self.content)
    elif self.element==["eSearchResult", "WarningList", "OutputMessage"]:
        self.record["WarningList"]["OutputMessage"].append(self.content)
    elif self.element==["eSearchResult", "WarningList", "QuotedPhraseNotFound"]:
        self.record["WarningList"]["QuotedPhraseNotFound"].append(self.content)
    elif self.element==["eSearchResult", "ERROR"]:
        # Not sure when this occurs. Are we supposed to raise an Exception?
        self.record["ERROR"] = self.content
