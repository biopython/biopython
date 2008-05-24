# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch
# from the "journals" database. No DTD file seems to be associated with
# these XML files.
# The code is not meant to be used by itself, but is used
# from Bio.Entrez.__init__.py.

errors = ()
integers = ()
strings = (
    "SortSerialName",
    "Year",
    "Month",
    "Day",
    "IndexingSubset",
    "IndexOnlineYN",
    "CurrentlyIndexedYN",
    "Coverage",
    "MinorTitleChangeYN",
    "Coden",
    "AcidFreeYN",
    "ContinuationNotes",
    "NlmUniqueID",
    "Title",
    "MedlineTA",
    "Country",
    "Place",
    "Publisher",
    "PublicationFirstYear",
    "PublicationEndYear",
    "Frequency",
    "ISOAbbreviation",
    "ISSN",
    "XrTitle",
    "Language",
    "BroadJournalHeading",
    "CurrentlyIndexedForSubset",
)

lists = (
    "IndexingHistoryList",
    "BroadJournalHeadingList",
    "CrossReferenceList",
    "SerialSet",
)

dictionaries = (
    "CrossReference",
    "PublicationInfo",
    "DateOfAction",
    "IlsCreatedTimestamp",
    "IlsUpdatedTimestamp",
    "IndexingHistory",
)

structures = {
    "Serial": ["Language"],
}

items = ()
