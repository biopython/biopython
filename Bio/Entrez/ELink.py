# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eLink,
# as specified by NCBI's DTD file eLink_020511.dtd (2005-02-18 17:13:40)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

class AttributedString(str): pass

def startElement(self, name, attrs):
    if self.element==["eLinkResult"]:
        object = []
        self.path = []
        self.record = object
    elif name=="LinkSet":
        object = {}
        self.path[-1].append(object)
    elif name in ("IdList", "IdUrlList", "IdCheckList"):
        object = []
        self.path[-1][name] = object
    elif name=="LinkSetDb":
        object = {"Link": []}
        self.path[-1][name] = object
    elif name in ("Link", "LinkInfo"):
        object = {}
        self.path[-1][name].append(object)
    elif name=="IdUrlSet":
        object = {"ObjUrl": []}
        self.path[-1].append(object)
    elif name=="ObjUrl":
        object = {"SubjectType": [], "Attribute": []}
        self.path[-1][name].append(object)
    elif name=="Provider":
        object = {}
        self.path[-1][name] = object
    elif name=="Id":
        object = ""
        attributes = {}
        if "HasLinkOut" in attrs:
            if attrs["HasLinkOut"]=='Y': attributes["HasLinkOut"] = True
            elif attrs["HasLinkOut"]=='N': attributes["HasLinkOut"] = False
        if "HasNeighbor" in attrs:
            if attrs["HasNeighbor"]=='Y': attributes["HasNeighbor"] = True
            elif attrs["HasNeighbor"]=='N': attributes["HasNeighbor"] = False
        self.attributes = attributes
    elif name=="IdLinkSet":
        # This is not in the DTD, but this is what I found when using the
        # following query:
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12169658,11748140",
        #                      cmd="acheck")
        object = {"LinkInfo": []}
        self.path[-1].append(object)
    else:
        object = ""
    self.path.append(object)

def endElement(self, name):
    self.path.pop()
    if name=="Error":
        error = self.content
        raise RuntimeError(error)
    if name=="Id":
        attributes = self.attributes
        value = AttributedString(self.content)
        value.attributes = attributes
        if type(self.path[-1])==list:
            self.path[-1].append(value)
        elif type(self.path[-1])==dict:
            self.path[-1][name] = value
    elif name in ("DbFrom",
                  "DbTo",
                  "Name",
                  "NameAbbr",
                  "IconUrl",
                  "MenuTag",
                  "HtmlTag",
                  "Url",
                  "LinkName"):
        self.path[-1][name] = self.content
    elif name in ("Score", "Priority"):
        self.path[-1][name] = int(self.content)
    elif name in ("SubjectType", "Attribute"):
        self.path[-1][name].append(self.content)
