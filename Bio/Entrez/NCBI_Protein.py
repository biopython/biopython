# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results in the format specified by
# NCBI's DTD file NCBI_Protein.mod.dtd (04/10/2008 16:04:22).
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = None

booleans = (
)

integers = (
)

strings = (
    "Prot-ref_name_E",		# (#PCDATA)
    "Prot-ref_desc",		# (#PCDATA)
    "Prot-ref_ec_E",		# (#PCDATA)
    "Prot-ref_activity_E",	# (#PCDATA)
    "Prot-ref_processed",	# %ENUM;
				# ATTLIST value (not-set |
				#                preprotein |
				#                mature |
				#                signal-peptide |
				#                transit-peptide) #REQUIRED
)

lists = (
    "Prot-ref_name",		# (Prot-ref_name_E*)
    "Prot-ref_db",		# (Dbtag*)
    "Prot-ref_ec",		# (Prot-ref_ec_E*)
    "Prot-ref_activity",	# (Prot-ref_activity_E*)
)

dictionaries = (
    "Prot-ref",		# ( Prot-ref_name?,
			#   Prot-ref_desc?,
			#   Prot-ref_ec?,
			#   Prot-ref_activity?,
			#   Prot-ref_db?,
			#   Prot-ref_processed?)
)

structures = {}

items = ()
