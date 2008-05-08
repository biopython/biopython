# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eSummary.
# as specified by NCBI's DTD file eSummary_041029.dtd (2004-10-29 15:52:04)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if name=="eSummaryResult":
        object = []
        self.record = object
        self.path = []
    elif name=="DocSum":
        object = {}
        self.path[-1].append(object)
    elif name=="Item":
        previous = self.path[-1]
        itemname = str(attrs["Name"]) # convert from Unicode
        itemtype = str(attrs["Type"]) # convert from Unicode
        if itemtype=="Structure" or itemname in ("ArticleIds", "History"):
            object = {}
        elif itemtype=="List":
            object = []
        else:
            object = ""
            self.itemname = itemname
            self.itemtype = itemtype
            if itemname in ("pubmed", "medline"):
                previous[itemname] = []
        if object!="":
            if type(previous)==dict:
                previous[itemname] = object
            elif type(previous)==list:
                previous.append(object)
    else:
        object = ""
    self.path.append(object)

def endElement(self, name):
    current = self.path.pop()
    if name=="ERROR":
        raise RuntimeError(self.content)
    # Integer|Date|String|Structure|List|Flags|Qualifier|Enumerator|Unknown
    if current!="": return
    previous = self.path[-1]
    if name=="Id":
        previous[name] = self.content
    elif name=="Item":
        if type(previous)==list:
            previous.append(self.content)
        elif type(previous)==dict:
            itemname = self.itemname
            itemtype = self.itemtype
            value = self.content
            if itemtype=="Integer": value = int(value)
            if itemname in ("pubmed", "medline"):
                previous[itemname].append(value)
            else:
                previous[itemname] = value
