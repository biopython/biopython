# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eInfo.
# as specified by NCBI's DTD file eInfo_020511.dtd (2006-12-04 21:51:33)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = "ERROR"		# (#PCDATA)>	<!-- .+ -->

booleans = (
    "IsDate",		# (#PCDATA)	<!-- (Y|N) -->
    "IsNumerical",	# (#PCDATA)	<!-- (Y|N) -->
    "SingleToken",	# (#PCDATA)	<!-- (Y|N) -->
    "Hierarchy",	# (#PCDATA)	<!-- (Y|N) -->
    "IsHidden",		# (#PCDATA)	<!-- (Y|N) -->
)

integers = (
    "TermCount",	# (#PCDATA)	<!-- \d+ -->
    "Count",		# (#PCDATA)	<!-- \d+ -->
)

strings = (
    "DbName",		# (#PCDATA)>    <!-- \S+ -->
    "Name",		# (#PCDATA)>    <!-- .+ -->
    "Menu",		# (#PCDATA)>    <!-- .+ -->
    "DbTo",		# (#PCDATA)>    <!-- \S+ -->
    "LastUpdate",	# (#PCDATA)>    <!-- \d+ -->
    "FullName",		# (#PCDATA)>    <!-- .+ -->
    "MenuName",		# (#PCDATA)>    <!-- .+ -->
    "Description",	# (#PCDATA)>    <!-- .+ -->
)

lists = (
    "DbList",		# (DbName+)
    "FieldList",	# (Field*)
    "LinkList",		# (Link*)
)

dictionaries = (
    "eInfoResult",	# (DbList|DbInfo|ERROR)
    "DbInfo",	# (DbName, MenuName, Description, Count, LastUpdate, FieldList, LinkList?)
    "Field",	# (Name, FullName, Description, TermCount, IsDate, IsNumerical, SingleToken, Hierarchy, IsHidden)
    "Link",	# (Name,Menu,Description,DbTo)
)

structures = {}

items = ()

def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path = self.path[:-1]
