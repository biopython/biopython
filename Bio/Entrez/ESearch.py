# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eSearch.
# as specified by NCBI's DTD file eSearch_020511.dtd (2006-06-28 17:35:21)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = "ERROR"	# (#PCDATA)	<!-- .+ -->

booleans = (
    "Explode",	# (#PCDATA)	<!-- (Y|N) -->
)

integers = (
    "Count",	# (#PCDATA)	<!-- \d+ -->
    "RetMax",	# (#PCDATA)	<!-- \d+ -->
    "RetStart",	# (#PCDATA)	<!-- \d+ -->
)

strings = (
    "QueryKey",		# (#PCDATA)	<!-- \d+ -->
    "QueryTranslation",	# (#PCDATA)>	<!-- .+ -->
    "WebEnv",		# (#PCDATA)	<!-- \S+ -->
    "From",		# (#PCDATA)	<!-- .+ -->
    "To",		# (#PCDATA)	<!-- .+ -->
    "Term",		# (#PCDATA)	<!-- .+ -->
    "Field",		# (#PCDATA)	<!-- .+ -->
    "Id",		# (#PCDATA)	<!-- \d+ -->
    "OP",		# (#PCDATA)	<!-- (AND|OR|NOT|RANGE|GROUP) -->
    "PhraseNotFound",		# (#PCDATA)	<!-- .+ -->
    "FieldNotFound",		# (#PCDATA)	<!-- .+ -->
    "PhraseIgnored",		# (#PCDATA)	<!-- .+ -->
    "OutputMessage",		# (#PCDATA)	<!-- .+ -->
    "QuotedPhraseNotFound",	# (#PCDATA)	<!-- .+ -->
)

lists = (
    "IdList",			# (Id*)
    "TranslationSet",		# (Translation*)
    "TranslationStack",		# ((TermSet|OP)*)
)

dictionaries = (
    "eSearchResult",	# (((Count,
			#       (RetMax,
			#        RetStart,
			#        QueryKey?,
			#        WebEnv?,
			#        IdList,
			#        TranslationSet,
			#        TranslationStack?,
			#        QueryTranslation
			#       )?
			#   ) | ERROR
			#  ),
			#  ErrorList?,
			#  WarningList?
			# )
    "Translation",	# (From, To)
    "TermSet",		# (Term, Field, Count, Explode)
)

structures = {
    "ErrorList": ("PhraseNotFound", "FieldNotFound"),
	# (PhraseNotFound*,FieldNotFound*)
    "WarningList": ("PhraseIgnored", "QuotedPhraseNotFound", "OutputMessage"),
	# (PhraseIgnored*, QuotedPhraseNotFound*, OutputMessage*)
}

items = ()
