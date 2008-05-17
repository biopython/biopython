# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results in the format specified by
# NCBI's DTD nlmcommon_080101.dtd (January 1, 2008)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.


error = None

booleans = (
)

integers = (
)

strings = (
    "AbstractText",		# (#PCDATA)
    "AccessionNumber",		# (#PCDATA)
    "Acronym",			# (#PCDATA)
    "Affiliation",		# (#PCDATA)
    "Agency",			# (#PCDATA)
    "ArticleTitle",		# (#PCDATA)
    "BroadJournalHeading",	# (#PCDATA)
    "Coden",			# (#PCDATA)
    "CollectionTitle",		# (#PCDATA)
    "CollectiveName",		# (#PCDATA)
    "CopyrightInformation",	# (#PCDATA)
    "Country",			# (#PCDATA)
    "Coverage",			# (#PCDATA)
    "DataBankName",		# (#PCDATA)
    "DateIssued",		# (#PCDATA)
    "DatesAssociatedWithName",	# (#PCDATA)
    "DatesOfSerialPublication",	# (#PCDATA)
    "Day",			# (#PCDATA)
    "DescriptorName",		# (#PCDATA)
				#  ATTLIST "MajorTopicYN (Y | N) "N"
    "Edition",			# (#PCDATA)>
    "ELocationID",		# (#PCDATA)>
				#  ATTLIST
				#   "EIdType %EIdType;
				#   "ValidYN  (Y | N) "Y"
    "EndPage",			# (#PCDATA)
    "FirstName",		# (#PCDATA)
    "ForeName",			# (#PCDATA)>
    "Frequency",		# (#PCDATA)>
				# !ATTLIST
				# FrequencyType (Current | Former) "Current"
    "GrantID",			# (#PCDATA)>
    "Hour",			# (#PCDATA)>
    "Imprint",			# (#PCDATA)>
    "Initials",			# (#PCDATA)>
    "ISOAbbreviation",		# (#PCDATA)>
    "ISSN",			# (#PCDATA)>
				#  ATTLIST IssnType
				#          (Electronic | Print | Undetermined)
				#           #REQUIRED
    "ISSNLinking",		# (#PCDATA)>
    "Issue",			# (#PCDATA)>
    "Language",			# (#PCDATA)>
    "LastName",			# (#PCDATA)>
    "MedlineDate",		# (#PCDATA)>
    "MedlinePgn",		# (#PCDATA)>
    "MedlineTA",		# (#PCDATA)>
    "MiddleName",		# (#PCDATA)
    "Minute",			# (#PCDATA)
    "Month",			# (#PCDATA)
    "NameQualifier",		# (#PCDATA)
    "NlmUniqueID",		# (#PCDATA)
    "OtherInformation",		# (#PCDATA)
    "Place",			# (#PCDATA)
				#  ATTLIST "ImprintType %ImprintType; "Current"
    "PlaceCode",		# (#PCDATA)>
    "PMID",			# (#PCDATA)>
    "ProjectedPublicationDate",	# (#PCDATA)>
    "PublicationEndYear",	# (#PCDATA)>
    "PublicationFirstYear",	# (#PCDATA)>
    "PublicationType",		# (#PCDATA)>
    "Publisher",		# (#PCDATA)>
    "QualifierName",		# (#PCDATA)>
				#  ATTLIST MajorTopicYN (Y | N) "N"
    "Season",			# (#PCDATA)>
    "Second",			# (#PCDATA)>
    "StartPage",		# (#PCDATA)>
    "Suffix",			# (#PCDATA)>
    "Title",			# (#PCDATA)>
    "TitleAssociatedWithName",	# (#PCDATA)>
    "VernacularTitle",		# (#PCDATA)>
    "Volume",			# (#PCDATA)>
    "Year",			# (#PCDATA)>
)

lists = (
    "AccessionNumberList",	# (AccessionNumber+)
    "AuthorList",		# (Author+) ATTLIST CompleteYN (Y | N) "Y"
    "BroadJournalHeadingList",	# (BroadJournalHeading+)
    "DataBankList",		# (DataBank+) ATTLIST CompleteYN (Y | N) "Y"
    "GrantList",		# (Grant+) ATTLIST CompleteYN (Y | N) "Y"
    "InvestigatorList",		# (Investigator+)
    "MeshHeadingList",		# (MeshHeading+)
    "PublicationTypeList",	# (PublicationType+)
)

dictionaries = (
    "Abstract",		# (%Abstract;)
    "ArticleDate",	# (%normal.date;)
			# ATTLIST "DateType CDATA  #FIXED "Electronic"
    "Author",		# ((%author.name;), Affiliation?,
			#  DatesAssociatedWithName?, NameQualifier?,
			#  OtherInformation?,TitleAssociatedWithName?)
			#   ATTLIST Author ValidYN (Y | N) "Y"
    "Book",		# (%PubDate.Ref;, Publisher, Title, AuthorList?,
			#  CollectionTitle?, Volume?)
    "DataBank",		# (DataBankName, AccessionNumberList?)
    "Grant",		# (%GrantID.Ref;, %Acronym.Ref;, %Agency.Ref;)
    "Investigator",	# (%personal.name;, Affiliation?)
    "Journal",		# (%ISSN.Ref;, JournalIssue, Coden?, Title?,
			#  ISOAbbreviation?)>
    "JournalIssue",	# (Volume?, Issue?, %PubDate.Ref;)
			#  ATTLIST JournalIssue
			#   CitedMedium (Internet | Print) #REQUIRED
    "MedlineJournalInfo",	# (Country?, MedlineTA, %NlmUniqueID.Ref;,
				#  ISSNLinking?)
    "NCBIArticle",	# (PMID, Article, MedlineJournalInfo,InvestigatorList?)
    "Pagination",	# ((StartPage, EndPage?, MedlinePgn?) | MedlinePgn)
    "PubDate",		# (%pub.date;)
)

structures = {
    "Article": ["ELocationID", "ArticleDate", "Language"],
			# ((Journal | Book), %ArticleTitle.Ref;,
			#  ((Pagination, ELocationID*) | ELocationID+),
			#  Abstract?, Affiliation?, AuthorList?, Language+,
			#  DataBankList?, GrantList?, %PublicationType.Ref;,
			#  VernacularTitle?, ArticleDate*)
    "PublicationInfo": ["Imprint", "Place", "Publisher", "DateIssued",
                        "DateofSerialPublication", "Frequency"],
			# (Country?, PlaceCode?, Imprint*, Place*, 
			#  Publisher*, DateIssued*,  ProjectedPublicationDate?,
			#  PublicationFirstYear?, PublicationEndYear?, Edition?,
			#  DatesOfSerialPublication*, Frequency*)
    "MeshHeading": ["QualifierName"],	# (DescriptorName, QualifierName*)
}

items = ()


def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path.pop()
