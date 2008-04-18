# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eSummary.
# as specified by NCBI's DTD file eSummary_041029.dtd (2004-10-29 15:52:04)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["eSummaryResult"]:
        self.record = []
    elif self.element==["eSummaryResult", "DocSum"]:
        self.record.append({})
    elif self.element[:3]==["eSummaryResult", "DocSum", "Item"]:
        docsum = self.record[-1]
        itemname = str(attrs["Name"]) # convert from Unicode
        itemtype = str(attrs["Type"]) # convert from Unicode
        item = [itemname, itemtype, None]
        if itemtype in ["List", "Structure"]: item[2] = []
        if not "Item" in docsum:
            docsum["Item"] = [item]
        else:
            n = self.element.count("Item")
            current = docsum["Item"]
            while n > 1:
                current = current[-1][2]
                n-=1
            current.append(item)

def endElement(self, name):
    if self.element==["eSummaryResult", "ERROR"]:
        raise ValueError(self.content)
    else:
        docsum = self.record[-1]
        if self.element==["eSummaryResult", "DocSum", "Id"]:
            docsum["Id"] = self.content
        elif self.element[:3]==["eSummaryResult", "DocSum", "Item"]:
            n = self.element.count("Item")
            current = docsum["Item"]
            while n > 1:
                current = current[-1][2]
                n-=1
            # Integer|Date|String|Structure|List|Flags|Qualifier|Enumerator|Unknown
            if current[-1][1]=="Integer":
                current[-1][2] = int(self.content)
            elif not current[-1][1] in ("List", "Structure"):
                current[-1][2] = self.content
