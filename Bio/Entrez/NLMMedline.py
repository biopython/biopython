# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results in the format specified by
# NCBI's DTD file nlmmedline_080101.dtd (January 1, 2008).
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.


error = None

booleans = (
)

integers = (
)

strings = (
)

lists = (
    "DeleteCitation",	# (PMID+)
)

dictionaries = (
)

structures = {
    "MedlineCitationSet": ["MedlineCitation"],
			# (MedlineCitation*, DeleteCitation?)
}

items = ()
