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
        self.record[-1]["IdUrlList"].append({"ObjUrl": []})
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl"]:
        d = {"SubjectType": [], "Attribute": []}
        self.record[-1]["IdUrlList"][-1]["ObjUrl"].append(d)
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["Provider"] = {}
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList"]:
        self.record[-1]["IdCheckList"] = []
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "Id"]:
        row = [None, None, None]
        if "HasLinkOut" in attrs:
            if attrs["HasLinkOut"]=='Y': row[1] = True
            elif attrs["HasLinkOut"]=='N': row[1] = False
        if "HasNeighbor" in attrs:
            if attrs["HasNeighbor"]=='Y': row[2] = True
            elif attrs["HasNeighbor"]=='N': row[2] = False
        self.record[-1]["IdCheckList"].append(row)
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet"]:
        # This is not in the DTD, but this is what I found when using the
        # following query:
        # >>> Bio.Entrez.elink(dbfrom="pubmed", id="12169658,11748140",
        #                      cmd="acheck")
        d = {"LinkInfo": []}
        self.record[-1]["IdCheckList"].append(d)
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet", "LinkInfo"]:
        self.record[-1]["IdCheckList"][-1]["LinkInfo"].append({})

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
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["Url"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "IconUrl"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["IconUrl"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "LinkName"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["LinkName"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "SubjectType"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["SubjectType"].append(self.content)
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Attribute"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["Attribute"].append(self.content)
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "Name"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["Provider"]["Name"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "NameAbbr"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["Provider"]["NameAbbr"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "Id"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["Provider"]["Id"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "Url"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["Provider"]["Url"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdUrlList", "IdUrlSet", "ObjUrl", "Provider", "IconUrl"]:
        self.record[-1]["IdUrlList"][-1]["ObjUrl"][-1]["Provider"]["IconUrl"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "Id"]:
        self.record[-1]["IdCheckList"][-1][0] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet", "Id"]:
        self.record[-1]["IdCheckList"][-1]["Id"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet", "LinkInfo", "DbTo"]:
        self.record[-1]["IdCheckList"][-1]["LinkInfo"][-1]["DbTo"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet", "LinkInfo", "LinkName"]:
        self.record[-1]["IdCheckList"][-1]["LinkInfo"][-1]["LinkName"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet", "LinkInfo", "MenuTag"]:
        self.record[-1]["IdCheckList"][-1]["LinkInfo"][-1]["MenuTag"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet", "LinkInfo", "HtmlTag"]:
        self.record[-1]["IdCheckList"][-1]["LinkInfo"][-1]["HtmlTag"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet", "LinkInfo", "Url"]:
        self.record[-1]["IdCheckList"][-1]["LinkInfo"][-1]["Url"] = self.content
    elif self.element==["eLinkResult", "LinkSet", "IdCheckList", "IdLinkSet", "LinkInfo", "Priority"]:
        self.record[-1]["IdCheckList"][-1]["LinkInfo"][-1]["Priority"] = int(self.content)
