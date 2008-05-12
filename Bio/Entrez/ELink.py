# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eLink,
# as specified by NCBI's DTD file eLink_020511.dtd (2005-02-18 17:13:40)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

class AttributedString(str): pass

def startElement(self, name, attrs):
    if name=="eLinkResult":	# (LinkSet*, ERROR?)
        object = []
        self.path = []
        self.record = object
    elif name=="LinkSet":	# (DbFrom,  ((IdList?, LinkSetDb*) | IdUrlList | IdCheckList | ERROR)
        object = {}
        self.path[-1].append(object)
    elif name in ("IdList",		# (Id*)
                  "IdUrlList",		# (IdUrlSet*,ERROR?)
                  "IdCheckList"):	# (Id*,ERROR?)
        object = []
        self.path[-1][name] = object
    elif name=="LinkSetDb":	# (DbTo, LinkName, (Link*|Info), ERROR?)
        object = {"Link": []}
        self.path[-1][name] = object
    elif name in ("Link",	# (Id, Score?)
                  "LinkInfo"):
        object = {}
        self.path[-1][name].append(object)
    elif name=="IdUrlSet":	# (Id,(ObjUrl+|Info))
        object = {"ObjUrl": []}
        self.path[-1].append(object)
    elif name=="ObjUrl":	# (Url, IconUrl?, LinkName?, SubjectType*, Attribute*, Provider)
        object = {"SubjectType": [], "Attribute": []}
        self.path[-1][name].append(object)
    elif name=="Provider":	# (Name, NameAbbr, Id, Url, IconUrl?)
        object = {}
        self.path[-1][name] = object
    elif name=="Id":	# (#PCDATA)>	<!-- \d+ -->
			# ATTLIST
			# HasLinkOut  (Y|N)	#IMPLIED	
			# HasNeighbor (Y|N)	#IMPLIED
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
    if name=="ERROR":		# (#PCDATA)	<!-- .+ -->
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
    elif name in ("DbFrom",	# (#PCDATA)	<!-- \S+ -->
                  "DbTo",	# (#PCDATA)	<!-- \S+ -->
                  "Name",	# (#PCDATA)	<!-- .+ -->
                  "NameAbbr",	# (#PCDATA)	<!-- \S+ -->
                  "IconUrl",	# (#PCDATA)	<!-- \S+ -->
                  "Info",	# (#PCDATA)>	<!-- .+ -->
                  "Url",	# (#PCDATA)	<!-- \S+ -->
                  "LinkName",	# (#PCDATA)	<!-- \S+ -->
                  "MenuTag",	# Not in the DTD, but found in the XML output
                  "HtmlTag"):	# Not in the DTD, but found in the XML output
        self.path[-1][name] = self.content
    elif name in ("Score",	# (#PCDATA)	<!-- \d+ -->
                  "Priority"):
        self.path[-1][name] = int(self.content)
    elif name in ("SubjectType",	# (#PCDATA)	<!-- .+ -->
                  "Attribute"):		# (#PCDATA)	<!-- .+ -->
        self.path[-1][name].append(self.content)
