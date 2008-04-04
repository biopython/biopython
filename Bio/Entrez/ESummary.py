# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eSummary.
# as specified by NCBI's DTD file eSummary_041029.dtd (2004-10-29 15:52:04)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["eSummaryResult"]:
        self.record = {}
    elif self.element[:3]==["eSummaryResult", "DocSum", "Item"]:
        itemname = str(attrs["Name"]) # convert from Unicode
        if attrs["Type"]=="List":
            item = [itemname, []]
        elif attrs["Type"]=="String":
            item = [itemname, ""]
        if not "Item" in self.record:
            self.record["Item"] = [item]
        else:
            n = self.element.count("Item")
            current = self.record["Item"]
            while n > 1:
                current = current[-1][1]
                n-=1
            current.append(item)

def endElement(self, name):
    if self.element==["eSummaryResult", "DocSum", "Id"]:
        self.record["Id"] = self.content
    elif self.element[:3]==["eSummaryResult", "DocSum", "Item"]:
        n = self.element.count("Item")
        current = self.record["Item"]
        while n > 1:
            current = current[-1][1]
            n-=1
        if type(current[-1][1])==str:
            current[-1][1] = self.content
