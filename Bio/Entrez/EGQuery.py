# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's EGQuery,
# as specified by NCBI's DTD file egquery.dtd (2004-05-03 16:19:48)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["Result"]:
        self.record = {}
    elif self.element==["Result", "eGQueryResult"]:
        self.record["eGQueryResult"] = []
    elif self.element==["Result", "eGQueryResult", "ResultItem"]:
        self.record["eGQueryResult"].append({})

def endElement(self, name):
    if self.element==["Result", "Term"]:
        self.record["Term"] = self.content
    elif self.element==["Result", "eGQueryResult", "ResultItem", "DbName"]:
        self.record["eGQueryResult"][-1]["DbName"] = self.content
    elif self.element==["Result", "eGQueryResult", "ResultItem", "MenuName"]:
        self.record["eGQueryResult"][-1]["MenuName"] = self.content
    elif self.element==["Result", "eGQueryResult", "ResultItem", "Count"]:
        self.record["eGQueryResult"][-1]["Count"] = int(self.content)
    elif self.element==["Result", "eGQueryResult", "ResultItem", "Status"]:
        self.record["eGQueryResult"][-1]["Status"] = self.content
