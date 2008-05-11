# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's ESpell,
# as specified by NCBI's DTD file eSpell.dtd.
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["eSpellResult"]:
        object = {}
        self.path = []
        self.record = object
    elif name=="SpelledQuery":
        object = {"Original": [], "Replaced": []}
        self.path[-1][name] = object
    else:
        object = ""
    self.path.append(object)

def endElement(self, name):
    self.path.pop()
    if name=="ERROR":
        error = self.content
        if error:
            raise RuntimeError(error)
    if name in ("Database", "Query", "CorrectedQuery", "ERROR"):
        self.path[-1][name] = self.content
    elif name in ("Original", "Replaced"):
        self.path[-1][name].append(self.content)
