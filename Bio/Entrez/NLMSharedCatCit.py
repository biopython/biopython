# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results in the format specified by
# NCBI's DTD nlmsharedcatcit_080101.dtd (January 1, 2008)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.


error = None

booleans = (
)

integers = (
)

strings = (
    "GeneralNote",		# (#PCDATA)
    "Keyword",			# (#PCDATA)
				#  ATTLIST MajorTopicYN (Y | N) "N"
    "RegistryNumber",		# (#PCDATA)
    "SpaceFlightMission",	# (#PCDATA)
    "NameOfSubstance",		# (#PCDATA)
    "OtherID",			# (#PCDATA)
				# ATTLIST
				#  Prefix CDATA #IMPLIED
				#  Source %Source;       
)

lists = (
    "PersonalNameSubjectList",	# (PersonalNameSubject+)
    "KeywordList",		# (Keyword+)
				# ATTLIST Owner %Owner; "NLM"
    "ChemicalList",		# (Chemical+)
)

dictionaries = (
    "Chemical",			# (RegistryNumber, NameOfSubstance)
    "DateCompleted",		# (%normal.date;)
    "DateCreated",		# (%normal.date;)
    "DateRevised",		# (%normal.date;)
    "OtherAbstract",		# (%Abstract;)
				#  ATTLIST Type %Type;
    "PersonalNameSubject",	# (%personal.name;, DatesAssociatedWithName?,
				#  NameQualifier?, OtherInformation?,
				#  TitleAssociatedWithName?)
)

structures = {
}

items = ()


def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path.pop()
