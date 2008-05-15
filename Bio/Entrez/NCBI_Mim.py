# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch from
# the OMIM database, as specified by NCBI's DTD file NCBI_Mim.mod.dtd
# (01/18/2007 23:07:18)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = None

booleans = (
)

integers = (
    "Mim-date_year",			# (%INTEGER;)
    "Mim-date_month",			# (%INTEGER;)
    "Mim-date_day",			# (%INTEGER;)
    "Mim-reference_pubmedUID",		# (%INTEGER;)
    "Mim-reference_citationType",	# (%INTEGER;)
    "Mim-reference_number",		# (%INTEGER;)
    "Mim-cit_year",			# (%INTEGER;)
    "Mim-cit_number",			# (%INTEGER;)
    "Mim-author_index",			# (%INTEGER;)
    "Mim-link_numRelevant",		# (%INTEGER;)
    "Mim-entry_numGeneMaps",		# (%INTEGER;)
    "Mim-entry_mimType",		# (%INTEGER;)
					# ATTLIST value (none |
        				#                star |
        				#                caret |
        				#                pound |
        				#                plus |
        				#                perc) #IMPLIED 
    "Mim-link_num",			# (%INTEGER;)
    "Mim-reference_origNumber",		# (%INTEGER;)
)

strings = (
    "Mim-page_from",		# (#PCDATA)
    "Mim-page_to",		# (#PCDATA)
    "Mim-edit-item_author",	# (#PCDATA)
    "Mim-index-term_terms_E",	# (#PCDATA)
    "Mim-index-term_key",	# (#PCDATA)
    "Mim-reference_miscInfo",	# (#PCDATA)
    "Mim-reference_volume",	# (#PCDATA)
    "Mim-reference_edition",	# (#PCDATA)
    "Mim-reference_journal",	# (#PCDATA)
    "Mim-reference_series",	# (#PCDATA)
    "Mim-reference_publisher",	# (#PCDATA)
    "Mim-reference_place",	# (#PCDATA)
    "Mim-reference_commNote",	# (#PCDATA)
    "Mim-reference_bookTitle",	# (#PCDATA)
    "Mim-reference_primaryAuthor",	# (#PCDATA)
    "Mim-reference_otherAuthors",	# (#PCDATA)
    "Mim-reference_citationTitle",	# (#PCDATA)
    "Mim-cit_others",		# (#PCDATA)
    "Mim-cit_author",		# (#PCDATA)
    "Mim-author_name",		# (#PCDATA)
    "Mim-link_uids",		# (#PCDATA)
    "Mim-allelic-variant_aliases_E",	# (#PCDATA)
    "Mim-allelic-variant_name",		# (#PCDATA)
    "Mim-allelic-variant_number",	# (#PCDATA)
    "Mim-text_label",		# (#PCDATA)
    "Mim-entry_included_E",	# (#PCDATA)
    "Mim-entry_aliases_E",	# (#PCDATA)
    "Mim-entry_synonyms_E",	# (#PCDATA)
    "Mim-entry_locus",		# (#PCDATA)
    "Mim-entry_symbol",		# (#PCDATA)
    "Mim-entry_title",		# (#PCDATA)
    "Mim-entry_mimNumber",	# (#PCDATA)
    "Mim-entry_copyright",	# (#PCDATA)
    "Mim-text_text",		# (#PCDATA)
    "Mim-reference_type",	# %ENUM;
				# ATTLIST value (not-set |
				#                citation |
				#                book |
				#                personal-communication |
				#                book-citation) #REQUIRED
    "Mim-reference_ambiguous",	# EMPTY
				# ATTLIST value ( true | false ) #REQUIRED 
    "Mim-reference_noLink",	# EMPTY
				# ATTLIST value ( true | false ) #REQUIRED
    "Mim-entry_hasSummary",	# EMPTY
				# ATTLIST value ( true | false ) #REQUIRED
    "Mim-entry_hasSynopsis",	# EMPTY
				# ATTLIST value ( true | false ) #REQUIRED
)

lists = (
    "Mim-allelic-variant_aliases",	# (Mim-allelic-variant_aliases_E*)
    "Mim-allelic-variant_description",	# (Mim-text*)
    "Mim-allelic-variant_mutation",	# (Mim-text*)
    "Mim-reference_editors",		# (Mim-author*)
    "Mim-entries",			# (Mim-entry*)
    "Mim-entry_aliases",		# (Mim-entry_aliases_E*)
    "Mim-entry_allelicVariants",	# (Mim-allelic-variant*)
    "Mim-entry_attribution",		# (Mim-edit-item*)
    "Mim-entry_clinicalSynopsis",	# (Mim-index-term*)
    "Mim-entry_editHistory",		# (Mim-edit-item*)
    "Mim-entry_included",		# (Mim-entry_included_E*)
    "Mim-entry_references",		# (Mim-reference*)
    "Mim-entry_seeAlso",		# (Mim-cit*)
    "Mim-entry_summary",		# (Mim-text*)
    "Mim-entry_summaryAttribution",	# (Mim-edit-item*)
    "Mim-entry_summaryEditHistory",	# (Mim-edit-item*)
    "Mim-entry_synonyms",		# (Mim-entry_synonyms_E*)
    "Mim-entry_synopsisAttribution",	# (Mim-edit-item*)
    "Mim-entry_synopsisEditHistory",	# (Mim-edit-item*)
    "Mim-entry_text",			# (Mim-text*)
    "Mim-entry_textfields",		# (Mim-text*)
    "Mim-index-term_terms",		# (Mim-index-term_terms_E*)
    "Mim-reference_authors",		# (Mim-author*)
    "Mim-reference_pages",		# (Mim-page*)
    "Mim-set_mimEntries",		# (Mim-entry*)
)

dictionaries = (
    "Mim-allelic-variant_snpLinks",	# (Mim-link)
    "Mim-entry_genomeLinks",		# (Mim-link)
    "Mim-entry_nucleotideLinks",	# (Mim-link)
    "Mim-entry_medlineLinks",		# (Mim-link)
    "Mim-entry_proteinLinks",		# (Mim-link)
    "Mim-entry_structureLinks",		# (Mim-link)
    "Mim-entry_creationDate",		# (Mim-edit-item)
    "Mim-entry_summaryCreationDate",	# (Mim-edit-item)
    "Mim-entry_synopsisCreationDate",	# (Mim-edit-item)
    "Mim-text_neighbors",		# (Mim-link)
    "Mim-allelic-variant",		#  (Mim-allelic-variant_number, 
					#   Mim-allelic-variant_name, 
					#   Mim-allelic-variant_aliases?, 
					#   Mim-allelic-variant_mutation?, 
					#   Mim-allelic-variant_description?, 
					#   Mim-allelic-variant_snpLinks?)
    "Mim-author",	# (Mim-author_name, Mim-author_index)
    "Mim-cit",		# (Mim-cit_number,
			#  Mim-cit_author,
			#  Mim-cit_others,
			#  Mim-cit_year)
    "Mim-entry",	# (Mim-entry_mimNumber, 
        		#  Mim-entry_mimType, 
        		#  Mim-entry_title, 
        		#  Mim-entry_copyright?, 
        		#  Mim-entry_symbol?, 
        		#  Mim-entry_locus?, 
        		#  Mim-entry_synonyms?, 
        		#  Mim-entry_aliases?, 
        		#  Mim-entry_included?, 
        		#  Mim-entry_seeAlso?, 
        		#  Mim-entry_text?, 
        		#  Mim-entry_textfields?, 
        		#  Mim-entry_hasSummary?, 
        		#  Mim-entry_summary?, 
        		#  Mim-entry_summaryAttribution?, 
        		#  Mim-entry_summaryEditHistory?, 
        		#  Mim-entry_summaryCreationDate?, 
        		#  Mim-entry_allelicVariants?, 
        		#  Mim-entry_hasSynopsis?, 
        		#  Mim-entry_clinicalSynopsis?, 
        		#  Mim-entry_synopsisAttribution?, 
        		#  Mim-entry_synopsisEditHistory?, 
        		#  Mim-entry_synopsisCreationDate?, 
        		#  Mim-entry_editHistory?, 
        		#  Mim-entry_creationDate?, 
        		#  Mim-entry_references?, 
        		#  Mim-entry_attribution?, 
        		#  Mim-entry_numGeneMaps, 
        		#  Mim-entry_medlineLinks?, 
        		#  Mim-entry_proteinLinks?, 
        		#  Mim-entry_nucleotideLinks?, 
        		#  Mim-entry_structureLinks?, 
        		#  Mim-entry_genomeLinks?)
    "Mim-index-term",	# (Mim-index-term_key, Mim-index-term_terms)
    "Mim-link",		# (Mim-link_num, Mim-link_uids, Mim-link_numRelevant?)
    "Mim-reference",	# (Mim-reference_number, 
        		#  Mim-reference_origNumber?, 
        		#  Mim-reference_type?, 
        		#  Mim-reference_authors, 
        		#  Mim-reference_primaryAuthor, 
        		#  Mim-reference_otherAuthors, 
        		#  Mim-reference_citationTitle, 
        		#  Mim-reference_citationType?, 
        		#  Mim-reference_bookTitle?, 
        		#  Mim-reference_editors?, 
        		#  Mim-reference_volume?, 
        		#  Mim-reference_edition?, 
        		#  Mim-reference_journal?, 
        		#  Mim-reference_series?, 
        		#  Mim-reference_publisher?, 
        		#  Mim-reference_place?, 
        		#  Mim-reference_commNote?, 
        		#  Mim-reference_pubDate, 
        		#  Mim-reference_pages?, 
        		#  Mim-reference_miscInfo?, 
        		#  Mim-reference_pubmedUID?, 
        		#  Mim-reference_ambiguous, 
        		#  Mim-reference_noLink?)
    "Mim-date",		#  (Mim-date_year, Mim-date_month?, Mim-date_day?)
    "Mim-edit-item_modDate",	# (Mim-date)
    "Mim-page",		# (Mim-page_from, Mim-page_to?)
    "Mim-reference_pubDate",	# (Mim-date)
    "Mim-edit-item",	# (Mim-edit-item_author, Mim-edit-item_modDate)
    "Mim-text",		# (Mim-text_label, Mim-text_text, Mim-text_neighbors?)
    "Mim-set",		# (Mim-set_releaseDate, Mim-set_mimEntries)
    "Mim-set_releaseDate",	# (Mim-date)
)

structures = {}

items = ()

def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path.pop()
