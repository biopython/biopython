# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eInfo.
# as specified by NCBI's DTD file eInfo_020511.dtd (2006-12-04 21:51:33)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if name=="eInfoResult":
        object = {}
        self.record = object
        self.path = []
    else:
        previous = self.path[-1]
        if name in ("DbList", "FieldList", "LinkList"):
            object = []
        elif name in ("DbInfo", "Field", "Link"):
            object = {}
        else:
            object = ""
        if object=="":
            pass
        elif type(previous)==dict:
            previous[name] = object
        elif type(previous)==list:
            previous.append(object)
    self.path.append(object)

def endElement(self, name):
    self.path = self.path[:-1]
    if name=="ERROR":
        error = self.content
        raise RuntimeError(error)
    if name in ("DbName",
                "Name",
                "Menu",
                "DbTo",
                "LastUpdate",
                "Name",
                "FullName",
                "MenuName",
                "Description"):
        value = self.content
    elif name in ("TermCount", "Count"):
        value = int(self.content)
    elif name in ("IsDate",
                  "IsNumerical",
                  "SingleToken",
                  "Hierarchy",
                  "IsHidden"):
        if self.content=='Y':
            value = True
        elif self.content=='N':
            value = False
    else:
        return
    previous = self.path[-1]
    if type(previous)==list:
        previous.append(value)
    elif type(previous)==dict:
        previous[name] = value
