# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eSummary.
# as specified by NCBI's DTD file eSummary_041029.dtd (2004-10-29 15:52:04)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = "ERROR"		# (#PCDATA)	<!-- .+ -->

booleans = (
)

integers = (
)

strings = (
    "Id",		# (#PCDATA)	<!-- \d+ -->
)

lists = (
    "eSummaryResult",	# (DocSum|ERROR)+
)

dictionaries = (
    "DocSum",		# (Id, Item+)
)

structures = {}

items = (
    "Item",	# (#PCDATA|Item)*	<!-- .+ -->
		# ATTLIST Name CDATA #REQUIRED
		#         Type (Integer
		#              |Date
		#              |String
		#              |Structure
		#              |List
		#              |Flags
		#              |Qualifier
		#              |Enumerator
		#              |Unknown) #REQUIRED
)


def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path = self.path[:-1]
