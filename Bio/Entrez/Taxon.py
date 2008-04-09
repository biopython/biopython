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

def startElement(self, name, attrs):
    if self.element==["TaxaSet"]:
        self.record = []
        
        # Constant for accessing the Triplet
        self.TAXID = 0
        self.SCIENTIFICNAME = 1
        self.RANK = 2

    elif self.element==["TaxaSet", "Taxon"]:
        self._one_taxon_dico = {}

    elif self.element==["TaxaSet", "Taxon", "LineageEx", "Taxon"]:
        self._Triplet = [None, "", ""]
        # will be converted to a tuple (ncbi_taxon_id, Scientific Name, Rank)
                
    elif self.element == ["TaxaSet", "Taxon", "TaxId"]:
        self._one_taxon_dico["TaxId"] = None
        
    elif self.element == ["TaxaSet", "Taxon", "ScientificName"]:
        self._one_taxon_dico["ScientificName"] = ""
 
    elif self.element == ["TaxaSet", "Taxon", "Rank"]:
        self._one_taxon_dico["Rank"] = ""
        
    elif self.element == ["TaxaSet", "Taxon", "OtherNames"]:
        self._one_taxon_dico["OtherNames"] = []
                   
    elif self.element in (["TaxaSet", "Taxon", "Synonym"], \
                          ["TaxaSet", "Taxon", "Division"], \
                          ["TaxaSet", "Taxon", "Lineage"], \
                          ["TaxaSet", "Taxon", "ParentTaxId"], \
                          ["TaxaSet", "Taxon", "CreateDate"], \
                          ["TaxaSet", "Taxon", "UpdateDate"], \
                          ["TaxaSet", "Taxon", "PubDate"]):
        self._one_taxon_dico[self.element[2]] = ""
    
    elif self.element == ["TaxaSet", "Taxon", "GeneticCode"]:
        self._one_taxon_dico["GeneticCode"] = [None, ""] 
        # this will store the [GeneticCode Id, GeneticCode Name]

    elif self.element == ["TaxaSet", "Taxon", "MitoGeneticCode"]:
        self._one_taxon_dico["MitoGeneticCode"] = [None, ""]
        # this will store the [MitoGeneticCode Id, MitoGeneticCode Name]

    # we are now looking at the list of all parents: to be stored in the list LineageEx
    elif self.element == ["TaxaSet", "Taxon", "LineageEx"]:
        self._one_taxon_dico["LineageEx"] = []


def endElement(self, name):
    if self.element==["TaxaSet", "Taxon"]:
        record = (self._one_taxon_dico["TaxId"], \
                  self._one_taxon_dico["ScientificName"], \
                  self._one_taxon_dico["Rank"])
        self._one_taxon_dico["LineageEx"].append( record )
        # add the Source Taxon's TaxId, Sc. Name and Rank at the end of the
        # list so the list is complete
        self.record.append(self._one_taxon_dico)
    
    elif self.element==["TaxaSet", "Taxon", "LineageEx", "Taxon"]:
        self._one_taxon_dico["LineageEx"].append( tuple(self._Triplet) )
        
    elif self.element == ["TaxaSet", "Taxon", "TaxId"]:
        self._one_taxon_dico["TaxId"] = int(self.content)
    elif self.element == ["TaxaSet", "Taxon", "LineageEx", "Taxon", "TaxId"]:
        self._Triplet[self.TAXID] = int(self.content)
            
    elif self.element == ["TaxaSet", "Taxon", "ScientificName"]:
        self._one_taxon_dico["ScientificName"] = self.content
    elif self.element == ["TaxaSet", "Taxon", "LineageEx", "Taxon", "ScientificName"]:
        self._Triplet[self.SCIENTIFICNAME] = self.content

                    
    elif self.element == ["TaxaSet", "Taxon", "Rank"]:
        self._one_taxon_dico["Rank"] = self.content
    elif self.element == ["TaxaSet", "Taxon", "LineageEx", "Taxon", "Rank"]:
            self._Triplet[self.RANK] = self.content

    elif len(self.element) == 4 and self.element[:3] == ["TaxaSet", "Taxon", "OtherNames"]:
        self._one_taxon_dico["OtherNames"].append( (self.element[3], self.content) )        
            
    elif self.element == ["TaxaSet", "Taxon", "ParentTaxId"]:
        self._one_taxon_dico["ParentTaxId"] = int(self.content)
    
    elif self.element == ["TaxaSet", "Taxon", "Division"]:
        self._one_taxon_dico["Division"] = self.content
        
    elif self.element == ["TaxaSet", "Taxon", "GeneticCode", "GCId"]:
        self._one_taxon_dico["GeneticCode"][0] = int(self.content)
    elif self.element == ["TaxaSet", "Taxon", "GeneticCode", "GCName"]:
        self._one_taxon_dico["GeneticCode"][1] = self.content
        
    elif self.element == ["TaxaSet", "Taxon", "MitoGeneticCode", "MGCId"]:
        self._one_taxon_dico["MitoGeneticCode"][0] = int(self.content)
    elif self.element == ["TaxaSet", "Taxon", "MitoGeneticCode", "MGCName"]:
        self._one_taxon_dico["MitoGeneticCode"][1] = self.content
        
    elif self.element == ["TaxaSet", "Taxon", "Lineage"]:
        self._one_taxon_dico["Lineage"] = self.content
        
    elif self.element == ["TaxaSet", "Taxon", "CreateDate"]:
        self._one_taxon_dico["CreateDate"] = self.content
    elif self.element == ["TaxaSet", "Taxon", "UpdateDate"]:
        self._one_taxon_dico["UpdateDate"] = self.content
    elif self.element == ["TaxaSet", "Taxon", "PubDate"]:
        self._one_taxon_dico["PubDate"] = self.content

