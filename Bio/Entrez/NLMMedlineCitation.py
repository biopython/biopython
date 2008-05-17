# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results in the format specified by
# NCBI's DTD file nlmmedlinecitation_080101.dtd (January 1, 2008).
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.


error = None

booleans = (
)

integers = (
)

strings = (
    "CitationSubset",		# (#PCDATA)
    "GeneSymbol",		# (#PCDATA)
    "NlmDcmsID",		# (#PCDATA)
    "Note",			# (#PCDATA)
    "NumberOfReferences",	# (#PCDATA)
    "RefSource",		# (#PCDATA)
)

lists = (
    "GeneSymbolList",		# (GeneSymbol+)
)

dictionaries = (
    "CommentIn",		# (RefSource, PMID?, Note?)
    "CommentOn",		# (RefSource, PMID?, Note?)
    "OriginalReportIn",		# (RefSource, PMID?, Note?)
    "PartialRetractionIn",	# (RefSource, PMID?, Note?)
    "PartialRetractionOf",	# (RefSource, PMID?, Note?)
    "ErratumFor",		# (RefSource, PMID?, Note?)
    "ErratumIn",		# (RefSource, PMID?, Note?)
    "ReprintIn",		# (RefSource, PMID?, Note?)
    "ReprintOf",		# (RefSource, PMID?, Note?)
    "RepublishedFrom",		# (RefSource, PMID?, Note?)
    "RepublishedIn",		# (RefSource, PMID?, Note?)
    "RetractionIn",		# (RefSource, PMID?, Note?)
    "RetractionOf",		# (RefSource, PMID?, Note?)
    "SummaryForPatientsIn",	# (RefSource, PMID?, Note?)
    "UpdateIn",			# (RefSource, PMID?, Note?)
    "UpdateOf",			# (RefSource, PMID?, Note?)
)

structures = {
    "MedlineCitation": ["CitationSubset", "OtherID", "OtherAbstract",
                        "KeywordList", "SpaceFlightMission", "GeneralNote"],
			# (%NlmDcmsID.Ref;, %PMID.Ref;, %DateCreated.Ref;,
			#  DateCompleted?, DateRevised?, Article,
			#  MedlineJournalInf, ChemicalList?, CitationSubset*,
			#  CommentsCorrections?, GeneSymbolList?,
			#  MeshHeadingList?, NumberOfReferences?,
			#  PersonalNameSubjectList?, OtherID*, OtherAbstract*,
			#  KeywordList*, SpaceFlightMission*,
			#  InvestigatorList?, GeneralNote*)
			# ATTLIST Owner %Owner; "NLM" Status %Status;
    "CommentsCorrections": ["CommentOn", "CommentIn", "ErratumIn", "ErratumFor",
                            "PartialRetractionIn", "PartialRetractionOf",
                            "RepublishedFrom", "RepublishedIn", "RetractionOf",
                            "RetractionIn", "UpdateIn", "UpdateOf", 
                            "SummaryForPatientsIn", "OriginalReportIn",
                            "ReprintOf", "ReprintIn"],
 			# (CommentOn*, CommentIn*, ErratumIn*, ErratumFor*,
			#  PartialRetractionIn*, PartialRetractionOf*,
			#  RepublishedFrom*, RepublishedIn*, RetractionOf*,
			#  RetractionIn*, UpdateIn*, UpdateOf*, 
			#  SummaryForPatientsIn*, OriginalReportIn*,
			#  ReprintOf*, ReprintIn*)
}

items = ()


def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path.pop()
