# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's EGQuery,
# as specified by NCBI's DTD file egquery.dtd (2004-05-03 16:19:48)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = ""

booleans = (
)

integers = (
    "Count",		# (#PCDATA)	<!-- \d+ -->
)

strings = (
    "DbName",		# (#PCDATA)	<!-- .+ -->
    "MenuName",		# (#PCDATA)	<!-- .+ -->
    "Status",		# (#PCDATA)	<!-- .+ -->
    "Term",		# (#PCDATA)	<!-- .+ -->
)

lists = (
    "eGQueryResult",	# (ResultItem+)
)

dictionaries = (
    "Result",		# (Term, eGQueryResult)
    "ResultItem",	# (DbName, MenuName, Count, Status)
)

structures = {}

items = ()
