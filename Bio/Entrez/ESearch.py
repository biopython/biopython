# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eSearch.
# as specified by NCBI's DTD file eSearch_020511.dtd (2006-06-28 17:35:21)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.



def startElement(self, name, attrs):
    if name=="eSearchResult":
        object = {}
        self.record = object
        self.path = []
    else:
        if name in ("ErrorList", "WarningList", "Translation", "TermSet"):
            object = {}
        elif name in ("PhraseNotFound",
                      "FieldNotFound",
                      "PhraseIgnored",
                      "QuotedPhraseNotFound",
                      "OutputMessage",
                      "IdList",
                      "TranslationSet",
                      "TranslationStack"):
            object = []
        else:
            object = ""
        if object=="":
            pass
        else:
            previous = self.path[-1]
            if type(previous)==list:
                previous.append(object)
            elif type(previous)==dict:
                previous[name] = object
    self.path.append(object)

def endElement(self, name):
    self.path = self.path[:-1]
    if name=="ERROR":
        error = self.content
        raise RuntimeError(error)
    if name in ("PhraseNotFound",
                "FieldNotFound",
                "PhraseIgnored",
                "OutputMessage",
                "QuotedPhraseNotFound"):
        self.path[-1][name].append(self.content)
    else:
        previous = self.path[-1]
        if name in ("Count", "RetMax", "RetStart"):
            value = int(self.content)
        elif name in ("QueryKey",
                      "QueryTranslation",
                      "WebEnv",
                      "From",
                      "To",
                      "Term",
                      "Field",
                      "Id",
                      "OP"):
            value = self.content
        elif name=="Explode":
            if self.content=='Y':
                value = True
            elif self.content=='N':
                value = False
        else:
            return
        if type(previous)==dict:
            self.path[-1][name] = value
        elif type(previous)==list:
            self.path[-1].append(value)
