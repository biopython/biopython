# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results in the format specified by
# NCBI's DTD file NCBI_Organism.mod.dtd (04/10/2008 16:04:22).
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = None

booleans = (
)

integers = (
    "OrgName_gcode",		# (%INTEGER;)>
    "OrgName_mgcode",		# (%INTEGER;)>
    "OrgMod_subtype",		# (%INTEGER;)>
				# ATTLIST value (strain |
				#                substrain |
				#                type |
				#                subtype |
				#                variety |
				#                serotype |
				#                serogroup |
				#                serovar |
				#                cultivar |
				#                pathovar |
				#                chemovar |
				#                biovar |
				#                biotype |
				#                group |
				#                subgroup |
				#                isolate |
				#                common |
				#                acronym |
				#                dosage |
				#                nat-host |
				#                sub-species |
				#                specimen-voucher |
				#                authority |
				#                forma |
				#                forma-specialis |
				#                ecotype |
				#                synonym |
				#                anamorph |
				#                teleomorph |
				#                breed |
				#                gb-acronym |
				#                gb-anamorph |
				#                gb-synonym |
				#                culture-collection |
				#                bio-material |
				#                metagenome-source |
				#                old-lineage |
				#                old-name |
				#                other) #IMPLIED
    "TaxElement_fixed-level",	# (%INTEGER;)
				# ATTLIST value (other |
				#                family |
				#                order |
				#                class) #IMPLIED
)

strings = (
    "TaxElement_level",			# (#PCDATA)
    "BinomialOrgName_genus",		# (#PCDATA)
    "BinomialOrgName_species",		# (#PCDATA)
    "BinomialOrgName_subspecies",	# (#PCDATA)
    "OrgMod_subname",			# (#PCDATA)
    "OrgMod_attrib",			# (#PCDATA)
    "OrgName_div",			# (#PCDATA)
    "OrgName_lineage",			# (#PCDATA)
    "OrgName_attrib",			# (#PCDATA)
    "OrgName_name_virus",		# (#PCDATA)
    "Org-ref_syn_E",			# (#PCDATA)
    "Org-ref_mod_E",			# (#PCDATA)
    "Org-ref_taxname",			# (#PCDATA)
    "Org-ref_common",			# (#PCDATA)
    "TaxElement_name",			# (#PCDATA)
)

lists = (
    "PartialOrgName",		# (TaxElement*)
    "TaxElement",		# (TaxElement_fixed-level, 
				#  TaxElement_level?, 
				#  TaxElement_name)
    "MultiOrgName",		# (OrgName*)
    "OrgName_mod",		# (OrgMod*)
    "Org-ref_mod",		# (Org-ref_mod_E*)
    "Org-ref_db",		# (Dbtag*)
    "Org-ref_syn",		# (Org-ref_syn_E*)
)

dictionaries = (
    "Org-ref",		# (Org-ref_taxname?, 
			#  Org-ref_common?, 
			#  Org-ref_mod?, 
			#  Org-ref_db?, 
			#  Org-ref_syn?, 
			#  Org-ref_orgname?)>
    "Org-ref_orgname",	# (OrgName)>
    "OrgName",		# (OrgName_name?, 
			#  OrgName_attrib?, 
			#  OrgName_mod?, 
			#  OrgName_lineage?, 
			#  OrgName_gcode?, 
			#  OrgName_mgcode?, 
			#  OrgName_div?)
    "OrgName_name",	# (OrgName_name_binomial | 
			#  OrgName_name_virus | 
			#  OrgName_name_hybrid | 
			#  OrgName_name_namedhybrid | 
			#  OrgName_name_partial)
    "OrgName_name_binomial",		# (BinomialOrgName)
    "OrgName_name_hybrid",		# (MultiOrgName)
    "OrgName_name_namedhybrid",		# (BinomialOrgName)
    "OrgName_name_partial",		# (PartialOrgName)
    "OrgMod",		# (OrgMod_subtype, OrgMod_subname, OrgMod_attrib?)
    "BinomialOrgName",	# (BinomialOrgName_genus, 
			#  BinomialOrgName_species?, 
			#  BinomialOrgName_subspecies?)
)

structures = {}

items = ()
