# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results in the format specified by
# NCBI's DTD file NCBI_Gene.mod.dtd (04/10/2008 16:04:22).
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = None

booleans = (
)

integers = (
)

strings = (
    "Gene-ref_locus",		# (#PCDATA)
    "Gene-ref_allele",		# (#PCDATA)
    "Gene-ref_desc",		# (#PCDATA)
    "Gene-ref_maploc",		# (#PCDATA)
    "Gene-ref_syn_E",		# (#PCDATA)
    "Gene-ref_locus-tag",	# (#PCDATA)
    "Gene-ref_pseudo",		# EMPTY
				# ATTLIST value ( true | false ) "false"
)

lists = (
    "Gene-ref_db",	# (Dbtag*)
    "Gene-ref_syn",	# (Gene-ref_syn_E*)
)

dictionaries = (
    "Gene-ref",		# (Gene-ref_locus?, 
			#  Gene-ref_allele?, 
			#  Gene-ref_desc?, 
			#  Gene-ref_maploc?, 
			#  Gene-ref_pseudo?, 
			#  Gene-ref_db?, 
			#  Gene-ref_syn?, 
			#  Gene-ref_locus-tag?)
)

structures = {}

items = ()
