# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eInfo.
# as specified by NCBI's DTD file eInfo_020511.dtd (2006-12-04 21:51:33)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["eInfoResult", "DbList"]:
        self.record = []
    elif self.element==["eInfoResult", "DbInfo"]:
        self.record = {}
    elif self.element==["eInfoResult", "DbInfo", "FieldList"]:
        self.record["FieldList"] = []
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field"]:
        field = {}
        self.record["FieldList"].append(field)
    elif self.element==["eInfoResult", "DbInfo", "LinkList"]:
        self.record["LinkList"] = []
    elif self.element==["eInfoResult", "DbInfo", "LinkList", "Link"]:
        link = {}
        self.record["LinkList"].append(link)

def endElement(self, name):
    if self.element==["eInfoResult", "DbList", "DbName"]:
        self.record.append(self.content)
    elif self.element==["eInfoResult", "DbInfo", "DbName"]:
        self.record["DbName"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "MenuName"]:
        self.record["MenuName"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "Description"]:
        self.record["Description"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "Count"]:
        self.record["Count"] = int(self.content)
    elif self.element==["eInfoResult", "DbInfo", "LastUpdate"]:
        self.record["LastUpdate"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "Name"]:
        field = self.record["FieldList"][-1]
        field["Name"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "FullName"]:
        field = self.record["FieldList"][-1]
        field["FullName"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "Description"]:
        field = self.record["FieldList"][-1]
        field["Description"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "TermCount"]:
        field = self.record["FieldList"][-1]
        field["TermCount"] = int(self.content)
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "IsDate"]:
        field = self.record["FieldList"][-1]
        if self.content=='Y': field["IsDate"] = True
        elif self.content=='N': field["IsDate"] = False
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "IsNumerical"]:
        field = self.record["FieldList"][-1]
        if self.content=='Y': field["IsNumerical"] = True
        elif self.content=='N': field["IsNumerical"] = False
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "SingleToken"]:
        field = self.record["FieldList"][-1]
        if self.content=='Y': field["SingleToken"] = True
        elif self.content=='N': field["SingleToken"] = False
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "Hierarchy"]:
        field = self.record["FieldList"][-1]
        if self.content=='Y': field["Hierarchy"] = True
        elif self.content=='N': field["Hierarchy"] = False
    elif self.element==["eInfoResult", "DbInfo", "FieldList", "Field", "IsHidden"]:
        field = self.record["FieldList"][-1]
        if self.content=='Y': field["IsHidden"] = True
        elif self.content=='N': field["IsHidden"] = False
    elif self.element==["eInfoResult", "DbInfo", "LinkList", "Link", "Name"]:
        link = self.record["LinkList"][-1]
        link["Name"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "LinkList", "Link", "Menu"]:
        link = self.record["LinkList"][-1]
        link["Menu"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "LinkList", "Link", "Description"]:
        link = self.record["LinkList"][-1]
        link["Description"] = self.content
    elif self.element==["eInfoResult", "DbInfo", "LinkList", "Link", "DbTo"]:
        link = self.record["LinkList"][-1]
        link["DbTo"] = self.content
