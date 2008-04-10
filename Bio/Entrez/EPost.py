# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's ePost,
# as specified by NCBI's DTD file ePost_020511.dtd (2004-09-28 18:47:37)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.



def startElement(self, name, attrs):
    if self.element==["ePostResult"]:
        self.record = {}
    elif self.element==["ePostResult", "InvalidIdList"]:
        self.record["InvalidIdList"] = []

def endElement(self, name):
    if self.element==["ePostResult", "InvalidIdList", "Id"]:
        self.record["InvalidIdList"].append(self.content)
    elif self.element==["ePostResult", "QueryEnv"]:
        self.record["QueryEnv"] = self.content
    elif self.element==["ePostResult", "WebEnv"]:
        self.record["WebEnv"] = self.content
    elif self.element==["ePostResult", "Error"]:
        self.record["Error"] = self.content
