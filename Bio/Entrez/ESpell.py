# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's ESpell,
# as specified by NCBI's DTD file eSpell.dtd.
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["eSpellResult"]:
        self.record = {}
    elif self.element==["eSpellResult", "SpelledQuery"]:
        self.record["SpelledQuery"] = {"Original": [], "Replaced": []}

def endElement(self, name):
    if self.element==["eSpellResult", "Database"]:
        self.record["Database"] = self.content
    elif self.element==["eSpellResult", "Query"]:
        self.record["Query"] = self.content
    elif self.element==["eSpellResult", "CorrectedQuery"]:
        self.record["CorrectedQuery"] = self.content
    elif self.element==["eSpellResult", "SpelledQuery", "Original"]:
        self.record["SpelledQuery"]["Original"].append(self.content)
    elif self.element==["eSpellResult", "SpelledQuery", "Replaced"]:
        self.record["SpelledQuery"]["Replaced"].append(self.content)
    elif self.element==["eSpellResult", "ERROR"]:
        # Not sure when this occurs. Are we supposed to raise an Exception?
        self.record["ERROR"] = self.content
