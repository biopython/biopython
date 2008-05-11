# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's EGQuery,
# as specified by NCBI's DTD file egquery.dtd (2004-05-03 16:19:48)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if name=="Result":
        object = {}
        self.path = []
        self.record = object
    elif name=="eGQueryResult":
        object = []
        self.path[-1][name] = object
    elif name=="ResultItem":
        object = {}
        self.path[-1].append(object)
    else:
        object = ""
    self.path.append(object)

def endElement(self, name):
    self.path.pop()
    if name in ("Term", "DbName", "MenuName", "Status"):
        value = self.content
    elif name=="Count":
        value = int(self.content)
    else:
        return
    self.path[-1][name] = value
