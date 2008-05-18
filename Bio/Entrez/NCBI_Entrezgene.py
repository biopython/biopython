# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results specified by NCBI's DTD file
# NCBI_Entrezgene.mod.dtd (04/10/2008 16:04:22).
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = None

booleans = (
)

integers = (
    "Entrezgene_type",		# (%INTEGER;)
				# ATTLIST value (unknown |
				#                tRNA |
				#                rRNA |
				#                snRNA |
				#                scRNA |
				#                snoRNA |
				#                protein-coding |
				#                pseudo |
				#                transposon |
				#                miscRNA |
				#                other) #IMPLIED
    "Gene-commentary_type",	# (%INTEGER;)
				# ATTLIST value (genomic |
        			#                pre-RNA |
        			#                mRNA |
        			#                rRNA |
        			#                tRNA |
        			#                snRNA |
        			#                scRNA |
        			#                peptide |
        			#                other-genetic |
        			#                genomic-mRNA |
        			#                cRNA |
        			#                mature-peptide |
        			#                pre-protein |
        			#                miscRNA |
        			#                snoRNA |
        			#                property |
        			#                reference |
        			#                generif |
        			#                phenotype |
        			#                complex |
        			#                compound |
        			#                comment |
        			#                other) #IMPLIED
    "Gene-commentary_version",	# (%INTEGER;)
    "Gene-track_geneid",	# (%INTEGER;)
    "Gene-track_status",	# (%INTEGER;)
				# ATTLIST value (live |
				#                secondary |
				#                discontinued |
				#                newentry) #IMPLIED
    "Gene-source_src-int",	# (%INTEGER;)
)

strings = (
    "Entrezgene_summary",		# (#PCDATA)
    "Entrezgene_xtra-index-terms_E",	# (#PCDATA)
    "Gene-commentary_accession",	# (#PCDATA)
    "Gene-commentary_heading",		# (#PCDATA)
    "Gene-commentary_label",		# (#PCDATA)
    "Gene-commentary_text",		# (#PCDATA)
    "Gene-source_extra-terms",		# EMPTY
					# ATTLIST value ( true | false ) "false"
    "Gene-source_gene-display",		# EMPTY
					# ATTLIST value ( true | false ) "false"
    "Gene-source_locus-display",	# EMPTY
					# ATTLIST value ( true | false ) "false"
    "Gene-source_src",			# (#PCDATA)
    "Gene-source_src-str1",		# (#PCDATA)
    "Gene-source_src-str2",		# (#PCDATA)
    "Maps_display-str",			# (#PCDATA)
    "Maps_method_map-type",		# %ENUM;
					# ATTLIST value (cyto |
					#                bp |
					#                cM |
					#                cR |
					#                min) #REQUIRED
    "Maps_method_proxy",		# (#PCDATA)
    "Other-source_anchor",		# (#PCDATA)
    "Other-source_post-text",		# (#PCDATA)
    "Other-source_pre-text",		# (#PCDATA)
    "Other-source_url",			# (#PCDATA)
    "Xtra-Terms_tag",			# (#PCDATA)
    "Xtra-Terms_value",			# (#PCDATA)
)

lists = (
    "Entrezgene_location",	# (Maps*)
    "Entrezgene_locus",		# (Gene-commentary*)
    "Entrezgene_properties",	# (Gene-commentary*)
    "Entrezgene_refgene",	# (Gene-commentary*)
    "Entrezgene_homology",	# (Gene-commentary*)
    "Entrezgene_comments",	# (Gene-commentary*)
    "Entrezgene_unique-keys",	# (Dbtag*)
    "Entrezgene_xtra-index-terms",	# (Entrezgene_xtra-index-terms_E*)
    "Entrezgene_xtra-properties",	# (Xtra-Terms*)>
    "Entrezgene_xtra-iq",		# (Xtra-Terms*)
    "Entrezgene_non-unique-keys",	# (Dbtag*)
    "Entrezgene-Set",		# (Entrezgene*)
    "Gene-commentary_comment",	# (Gene-commentary*)
    "Gene-commentary_genomic-coords",	# (Seq-loc*)
    "Gene-commentary_products",	# (Gene-commentary*)
    "Gene-commentary_properties",	# (Gene-commentary*)
    "Gene-commentary_refs",	# (Pub*)
    "Gene-commentary_seqs",	# (Seq-loc*)
    "Gene-commentary_source",	# (Other-source*)
    "Gene-commentary_xtra-properties",	# (Xtra-Terms*)
    "Gene-track_current-id",	# (Dbtag*)
)

dictionaries = (
    "Entrezgene",	# (Entrezgene_track-info?, 
			#  Entrezgene_type,
			#  Entrezgene_source, 
			#  Entrezgene_gene, 
			#  Entrezgene_prot?, 
			#  Entrezgene_rna?, 
			#  Entrezgene_summary?, 
			#  Entrezgene_location?, 
			#  Entrezgene_gene-source?, 
			#  Entrezgene_locus?, 
			#  Entrezgene_properties?, 
			#  Entrezgene_refgene?, 
			#  Entrezgene_homology?, 
			#  Entrezgene_comments?, 
			#  Entrezgene_unique-keys?, 
			#  Entrezgene_xtra-index-terms?, 
			#  Entrezgene_xtra-properties?, 
			#  Entrezgene_xtra-iq?, 
			#  Entrezgene_non-unique-keys?)
    "Entrezgene_track-info",	# (Gene-track)
    "Entrezgene_source",	# (BioSource)
    "Entrezgene_gene",		# (Gene-ref)
    "Entrezgene_prot",		# (Prot-ref)
    "Entrezgene_rna",		# (RNA-ref)
    "Entrezgene_gene-source",	# (Gene-source)
    "Gene-commentary",		# (Gene-commentary_type, 
				#  Gene-commentary_heading?, 
				#  Gene-commentary_label?, 
				#  Gene-commentary_text?, 
				#  Gene-commentary_accession?, 
				#  Gene-commentary_version?, 
				#  Gene-commentary_xtra-properties?, 
				#  Gene-commentary_refs?, 
				#  Gene-commentary_source?, 
				#  Gene-commentary_genomic-coords?, 
				#  Gene-commentary_seqs?, 
				#  Gene-commentary_products?, 
				#  Gene-commentary_properties?, 
				#  Gene-commentary_comment?, 
				#  Gene-commentary_create-date?, 
				#  Gene-commentary_update-date?)
    "Gene-commentary_create-date",	# (Date)
    "Gene-commentary_update-date",	# (Date)
    "Gene-source",	# (Gene-source_src,
			#  Gene-source_src-int?,
			#  Gene-source_src-str1?,
			#  Gene-source_src-str2?,
			#  Gene-source_gene-display?,
			#  Gene-source_locus-display?,
			#  Gene-source_extra-terms?)
    "Gene-track",	# (Gene-track_geneid, 
			#  Gene-track_status?, 
			#  Gene-track_current-id?, 
			#  Gene-track_create-date, 
			#  Gene-track_update-date, 
			#  Gene-track_discontinue-date?)
    "Gene-track_create-date",		# (Date)
    "Gene-track_update-date",		# (Date)
    "Gene-track_discontinue-date",	# (Date)
    "Maps",		# (Maps_display-str, Maps_method)
    "Maps_method",	# (Maps_method_proxy | Maps_method_map-type)
    "Other-source",	# (Other-source_src?, 
			#  Other-source_pre-text?, 
			#  Other-source_anchor?, 
			#  Other-source_url?, 
			#  Other-source_post-text?)
    "Other-source_src",	# (Dbtag)
    "Xtra-Terms",	# (Xtra-Terms_tag,Xtra-Terms_value)
)

structures = {}

items = ()
