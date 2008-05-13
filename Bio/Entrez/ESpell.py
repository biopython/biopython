# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's ESpell,
# as specified by NCBI's DTD file eSpell.dtd.
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.


error = "ERROR"		# (#PCDATA)	<!-- \d+ -->

booleans = (
)

integers = (
)

strings = (
    "Original",		# (#PCDATA)	<!-- \d+ -->
    "Replaced",		# (#PCDATA)	<!-- \d+ -->
    "Database",		# (#PCDATA)	<!-- \d+ -->
    "Query",		# (#PCDATA)	<!-- \d+ -->
    "CorrectedQuery",	# (#PCDATA)	<!-- \d+ -->
)

lists = (
)

dictionaries = (
    "eSpellResult",	# (Database, Query, CorrectedQuery, SpelledQuery, ERROR)
)

structures = {
    "SpelledQuery": ["Replaced", "Original"],
			# (Replaced|Original)*	<!-- \d+ -->
}

items = ()

def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path = self.path[:-1]
