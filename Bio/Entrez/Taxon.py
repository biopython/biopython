# Copyright 2008 by Eric Gibert.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""
This code is used to parse XML results returned by Entrez's Taxonomy
database as specified by NCBI's DTD file taxon.dtd (2005-07-15 21:21:27)

The code is not meant to be used by itself, but by calling
Bio.Entrez.read() instead.  When parsing a Taxonomy XML file
in this way, the record returned is a list of dictionaries.

For example,

from Bio import Entrez
print "Multiple records search..."
search_handle = Entrez.esearch(db="taxonomy",
                               term="orthetrum c*",
                               retmode = "XML")
taxon_list = Entrez.read(search_handle)["IdList"]
print len(taxon_list), "Orthetrum match your search"

handle = Entrez.efetch(db="taxonomy", id=taxon_list, retmode="XML")
orthetrum_list = Entrez.read(handle)
assert len(orthetrum_list) == len(taxon_list)

for orthetrum in orthetrum_list:
    print "%s '%s' has TaxID %s and parent %i" \
        % (orthetrum["Rank"],
           orthetrum["ScientificName"], 
           orthetrum["TaxId"],
           orthetrum["ParentTaxId"])

This should print something like the following on screen:

Multiple records search...
4 Orthetrum match your search
species 'Orthetrum chrysostigma' has TaxID 488050 and parent 126255
species 'Orthetrum chrysis' has TaxID 447868 and parent 126255
species 'Orthetrum coerulescens' has TaxID 333459 and parent 126255
species 'Orthetrum cancellatum' has TaxID 126256 and parent 126255

"""

error = None

booleans = ()

integers = (
    "GCId",		# (#PCDATA)
    "MGCId",		# (#PCDATA)
    "ModId",		# (#PCDATA)
    "ParentTaxId",	# (#PCDATA)
    "RModId",		# (#PCDATA)
    "RTaxId",		# (#PCDATA)
    "TaxId",		# (#PCDATA)
)

strings = (
    "Division",		# (#PCDATA)
    "Rank",		# (#PCDATA)
    "ClassCDE",		# (#PCDATA)
    "DispName",		# (#PCDATA)
    "UniqueName",	# (#PCDATA)
    "GCName",		# (#PCDATA)
    "MGCName",		# (#PCDATA)
    "Lineage",		# (#PCDATA)
    "PropName",		# (#PCDATA)
    "CreateDate",	# (#PCDATA)
    "UpdateDate",	# (#PCDATA)
    "PubDate",		# (#PCDATA)
    "CitId",		# (#PCDATA)
    "CitKey",		# (#PCDATA)
    "CitUrl",		# (#PCDATA)
    "CitText",		# (#PCDATA)
    "CitPubmedId"	# (#PCDATA)
    "CitMedlineId"	# (#PCDATA)
    "ModType",		# (#PCDATA)
    "ModName",		# (#PCDATA)
    "ModGBhidden",	# (#PCDATA)
    "ScientificName",	# (#PCDATA)
    "GenbankCommonName",	# (#PCDATA)
    "GenbankAcronym",		# (#PCDATA)
    "BlastName",	# (#PCDATA)
    "EquivalentName",	# (#PCDATA)
    "Synonym",		# (#PCDATA)
    "Acronym",		# (#PCDATA)
    "Misspelling",	# (#PCDATA)
    "Anamorph",		# (#PCDATA)
    "Includes",		# (#PCDATA)
    "CommonName",	# (#PCDATA)
    "Inpart",		# (#PCDATA)
    "Misnomer",		# (#PCDATA)
    "Teleomorph",	# (#PCDATA)
    "GenbankSynonym",	# (#PCDATA)
    "GenbankAnamorph",	# (#PCDATA)
    "PropValueInt",	# (#PCDATA)
    "PropValueBool",	# (#PCDATA)
    "PropValueString",	# (#PCDATA)
)

lists = (
    "AkaTaxIds",	# ( TaxId* )
    "Citations",	# ( Citation+ )
    "LineageEx",	# ( Taxon* )
    "Modifiers",	# ( Modifier+ )
    "Properties",	# ( Property+ )
    "TaxaSet",		# ( Taxon+ )
)

dictionaries = (
    "Name",		# ( ClassCDE, DispName, UniqueName )
    "GeneticCode",	# ( GCId, GCName )
    "MitoGeneticCode",	# ( MGCId, MGCName )
    "Citation",		# ( CitId,
			#   CitKey,
			#   CitUrl?,
			#   CitText?,
			#   CitPubmedId?,
			#   CitMedlineId?
    "Modifier",		# ( ModId,
			#   ModType,
			#   ModName,
			#   ModGBhidden,
			#   ( RModId | RTaxId )? )
    "Property",		# ( PropName,
			# ( PropValueInt | PropValueBool | PropValueString ) )
    "Taxon",		# ( TaxId,
			#   ScientificName,
			#   OtherNames?, 
			#   ParentTaxId?,
			#   Rank?, 
			#   Division?,
			#   GeneticCode?,
			#   MitoGeneticCode?,
			#   Lineage?,
			#   LineageEx?,
			#   Citations?,
			#   Modifiers?,
			#   Properties?,
			#   CreateDate?,
			#   UpdateDate?,
			#   PubDate?,
			#   AkaTaxIds? )
)

structures = {
    "OtherNames": ["EquivalentName",
                   "Synonym", 
                   "Acronym", 
                   "Misspelling", 
                   "Anamorph", 
                   "Includes", 
                   "CommonName", 
                   "Inpart", 
                   "Misnomer", 
                   "Teleomorph", 
                   "GenbankSynonym", 
                   "GenbankAnamorph"], 
			# ( GenbankCommonName?,
			#   GenbankAcronym?,
			#   BlastName?,
			#   ( EquivalentName |
			#     Synonym        |
			#     Acronym        |
			#     Misspelling    |
			#     Anamorph       |
			#     Includes       |
			#     CommonName     |
			#     Inpart         |
			#     Misnomer       |
			#     Teleomorph     |
			#     GenbankSynonym |
			#     GenbankAnamorph
			#   )*,
			#   Name*
			# )
}

items = ()

def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path = self.path[:-1]
