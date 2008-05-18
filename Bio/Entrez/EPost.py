# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's ePost,
# as specified by NCBI's DTD file ePost_020511.dtd (2004-09-28 18:47:37)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.


error = "ERROR"		# (#PCDATA)>	<!-- .+ -->

booleans = (
)

integers = (
)

strings = (
    "Id",		# (#PCDATA)>	<!-- \d+ -->
    "QueryKey",		# (#PCDATA)>	<!-- \d+ -->
    "WebEnv",		# (#PCDATA)>	<!-- \S+ -->
)

lists = (
    "InvalidIdList",	# (Id+)
)

dictionaries = (
    "ePostResult",	# (InvalidIdList?,(QueryKey,WebEnv)?,ERROR?)
)

structures = {}

items = ()
