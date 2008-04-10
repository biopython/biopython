# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eLink,
# as specified by NCBI's DTD file eLink_020511.dtd (2005-02-18 17:13:40)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.



def startElement(self, name, attrs):
    if self.element==["eLinkResult"]:
        self.record = []
    elif self.element==["eLinkResult", "LinkSet"]:
        self.record.append({})
    elif self.element==["eLinkResult", "LinkSet", "IdList"]:
        self.record[-1]["IdList"] = []
    elif self.element==["eLinkResult", "LinkSet", "LinkSetDb"]:
        self.record[-1]["LinkSetDb"] = {"Link": []}
    if self.element==["eLinkResult", "LinkSet", "LinkSetDb", "Link"]:
        self.record[-1]["LinkSetDb"]["Link"].append({})
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList"]:
        self.record[-1]["IdUrlList"] = []
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet"]:
        self.record[-1]["IdUrlList"].append({})
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"] = {"SubjectType": [], "Attribute": []}
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["Provider"] = {}
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList"]:
        self.record[-1]["IdCheckList"] = []

def endElement(self, name):
    if name=="Error":
        error = self.content
        raise RuntimeError(error)
    if self.element==["eLinkResult", "LinkSet", "DbFrom"]:
        self.record[-1]["DbFrom"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdList", "Id"]:
        self.record[-1]["IdList"].append(self.content)
    elif self.element==["eLinkResult", "LinkSet", "LinkSetDb", "DbTo"]:
        self.record[-1]["LinkSetDb"]["DbTo"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "LinkSetDb", "LinkName"]:
        self.record[-1]["LinkSetDb"]["LinkName"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "LinkSetDb", "Link", "Id"]:
        self.record[-1]["LinkSetDb"]["Link"][-1]["Id"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "LinkSetDb", "Link", "Score"]:
        self.record[-1]["LinkSetDb"]["Link"][-1]["Score"] = int(self.content)
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "Id"]:
        self.record[-1]["IdUrlList"][-1]["Id"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Url"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["Url"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "IconUrl"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["IconUrl"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "SubjectType"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["SubjectType"].append(self.content)
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Attribute"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["Attribute"].append(self.content)
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "Name"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["Provider"]["Name"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "NameAbbr"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["Provider"]["NameAbbr"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "Id"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["Provider"]["Id"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "Url"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["Provider"]["Url"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "IconUrl"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"]["Provider"]["IconUrl"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "Id"]:
        self.record[-1]["IdCheckList"].append(self.content)
