# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results returned by Entrez's eFetch from
# the OMIM database, as specified by NCBI's DTD file NCBI_Mim.dtd
# (06/06/2006 23:03:48)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

def startElement(self, name, attrs):
    if self.element==["Mim-entries"]:
        self.record = []
    elif self.element==["Mim-entries", "Mim-entry"]:
        self.record.append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_mimType"]:
        self.record[-1]["Mim-entry_mimType"] = [str(attrs["value"]), None]
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_aliases"]:
        self.record[-1]["Mim-entry_aliases"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synonyms"]:
        self.record[-1]["Mim-entry_synonyms"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_included"]:
        self.record[-1]["Mim-entry_included"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_seeAlso"]:
        self.record[-1]["Mim-entry_seeAlso"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_seeAlso", "Mim-cit"]:
        self.record[-1]["Mim-entry_seeAlso"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_text"]:
        self.record[-1]["Mim-entry_text"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_text", "Mim-text"]:
        self.record[-1]["Mim-entry_text"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields"]:
        self.record[-1]["Mim-entry_textfields"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields", "Mim-text"]:
        self.record[-1]["Mim-entry_textfields"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields", "Mim-text", "Mim-text_neighbors"]:
        self.record[-1]["Mim-entry_textfields"][-1]["Mim-text_neighbors"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields", "Mim-text", "Mim-text_neighbors", "Mim-link"]:
        self.record[-1]["Mim-entry_textfields"][-1]["Mim-text_neighbors"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_hasSummary"]:
        if str(attrs["value"])=="true":
            self.record[-1]["Mim-entry_hasSummary"] = True
        elif str(attrs["value"])=="false":
            self.record[-1]["Mim-entry_hasSummary"] = False
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary"]:
        self.record[-1]["Mim-entry_summary"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary", "Mim-text"]:
        self.record[-1]["Mim-entry_summary"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryAttribution"]:
        self.record[-1]["Mim-entry_summaryAttribution"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryAttribution", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_summaryAttribution"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryAttribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_summaryAttribution"]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryEditHistory"]:
        self.record[-1]["Mim-entry_summaryEditHistory"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryEditHistory", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_summaryEditHistory"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryEditHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_summaryEditHistory"]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryCreationDate"]:
        self.record[-1]["Mim-entry_summaryCreationDate"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryCreationDate", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_summaryCreationDate"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryCreationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_summaryCreationDate"]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary", "Mim-text", "Mim-text_neighbors"]:
        self.record[-1]["Mim-entry_summary"]["Mim-text"]["Mim-text_neighbors"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary", "Mim-text", "Mim-text_neighbors", "Mim-link"]:
        self.record[-1]["Mim-entry_summary"]["Mim-text"]["Mim-text_neighbors"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants"]:
        self.record[-1]["Mim-entry_allelicVariants"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant"]:
        self.record[-1]["Mim-entry_allelicVariants"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant", "Mim-allelic-variant_aliases"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_aliases"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant", "Mim-allelic-variant_snpLinks"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_snpLinks"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant", "Mim-allelic-variant_snpLinks", "Mim-link"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_snpLinks"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_hasSynopsis"]:
        if str(attrs["value"])=="true":
            self.record[-1]["Mim-entry_hasSynopsis"] = True
        elif str(attrs["value"])=="false":
            self.record[-1]["Mim-entry_hasSynopsis"] = False
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis"]:
        self.record[-1]["Mim-entry_synopsis"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis", "Mim-text"]:
        self.record[-1]["Mim-entry_synopsis"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisAttribution"]:
        self.record[-1]["Mim-entry_synopsisAttribution"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisAttribution", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_synopsisAttribution"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisAttribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_synopsisAttribution"]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisEditHistory"]:
        self.record[-1]["Mim-entry_synopsisEditHistory"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisEditHistory", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_synopsisEditHistory"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisEditHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_synopsisEditHistory"]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisCreationDate"]:
        self.record[-1]["Mim-entry_synopsisCreationDate"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisCreationDate", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_synopsisCreationDate"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisCreationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_synopsisCreationDate"]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis", "Mim-text", "Mim-text_neighbors"]:
        self.record[-1]["Mim-entry_synopsis"]["Mim-text"]["Mim-text_neighbors"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis", "Mim-text", "Mim-text_neighbors", "Mim-link"]:
        self.record[-1]["Mim-entry_synopsis"]["Mim-text"]["Mim-text_neighbors"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_clinicalSynopsis"]:
        self.record[-1]["Mim-entry_clinicalSynopsis"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_clinicalSynopsis", "Mim-index-term"]:
        d = {"Mim-index-term_terms": []}
        self.record[-1]["Mim-entry_clinicalSynopsis"].append(d)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_editHistory"]:
        self.record[-1]["Mim-entry_editHistory"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_editHistory", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_editHistory"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_editHistory", "Mim-edit-item", "Mim-edit-item_modDate"]:
        self.record[-1]["Mim-entry_editHistory"][-1]["Mim-edit-item_modDate"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_editHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_editHistory"][-1]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_creationDate"]:
        self.record[-1]["Mim-entry_creationDate"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_creationDate", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_creationDate"]["Mim-edit-item"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_creationDate", "Mim-edit-item", "Mim-edit-item_modDate"]:
        self.record[-1]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_modDate"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_creationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references"]:
        self.record[-1]["Mim-entry_references"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference"]:
        self.record[-1]["Mim-entry_references"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_type"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_type"] = str(attrs["value"])
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_authors"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_authors"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_authors", "Mim-author"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_authors"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pubDate"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pubDate"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pubDate", "Mim-date"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pubDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pages"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pages"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pages", "Mim-page"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pages"]["Mim-page"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_ambiguous"]:
        if str(attrs["value"])=="true":
            self.record[-1]["Mim-entry_references"][-1]["Mim-reference_ambiguous"] = True
        elif str(attrs["value"])=="false":
            self.record[-1]["Mim-entry_references"][-1]["Mim-reference_ambiguous"] = False
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_noLink"]:
        if str(attrs["value"])=="true":
            self.record[-1]["Mim-entry_references"][-1]["Mim-reference_noLink"] = True
        elif str(attrs["value"])=="false":
            self.record[-1]["Mim-entry_references"][-1]["Mim-reference_noLink"] = False
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_attribution"]:
        self.record[-1]["Mim-entry_attribution"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_attribution", "Mim-edit-item"]:
        self.record[-1]["Mim-entry_attribution"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_attribution", "Mim-edit-item", "Mim-edit-item_modDate"]:
        self.record[-1]["Mim-entry_attribution"][-1]["Mim-edit-item_modDate"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_attribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date"]:
        self.record[-1]["Mim-entry_attribution"][-1]["Mim-edit-item_modDate"]["Mim-date"] = {}
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_medlineLinks"]:
        self.record[-1]["Mim-entry_medlineLinks"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_medlineLinks", "Mim-link"]:
        self.record[-1]["Mim-entry_medlineLinks"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_proteinLinks"]:
        self.record[-1]["Mim-entry_proteinLinks"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_proteinLinks", "Mim-link"]:
        self.record[-1]["Mim-entry_proteinLinks"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_nucleotideLinks"]:
        self.record[-1]["Mim-entry_nucleotideLinks"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_nucleotideLinks", "Mim-link"]:
        self.record[-1]["Mim-entry_nucleotideLinks"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_structureLinks"]:
        self.record[-1]["Mim-entry_structureLinks"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_structureLinks", "Mim-link"]:
        self.record[-1]["Mim-entry_structureLinks"].append({})
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_genomeLinks"]:
        self.record[-1]["Mim-entry_genomeLinks"] = []
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_genomeLinks", "Mim-link"]:
        self.record[-1]["Mim-entry_genomeLinks"].append({})

def endElement(self, name):
    if self.element==["Mim-entries", "Mim-entry", "Mim-entry_mimNumber"]:
        self.record[-1]["Mim-entry_mimNumber"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_mimType"]:
        self.record[-1]["Mim-entry_mimType"][1] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_title"]:
        self.record[-1]["Mim-entry_title"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_copyright"]:
        self.record[-1]["Mim-entry_copyright"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_symbol"]:
        self.record[-1]["Mim-entry_symbol"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_locus"]:
        self.record[-1]["Mim-entry_locus"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synonyms", "Mim-entry_synonyms_E"]:
        self.record[-1]["Mim-entry_synonyms"].append(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_aliases", "Mim-entry_aliases_E"]:
        self.record[-1]["Mim-entry_aliases"].append(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_included", "Mim-entry_included_E"]:
        self.record[-1]["Mim-entry_included"].append(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_seeAlso", "Mim-cit", "Mim-cit_number"]:
        self.record[-1]["Mim-entry_seeAlso"]["Mim-cit_number"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_seeAlso", "Mim-cit", "Mim-cit_author"]:
        self.record[-1]["Mim-entry_seeAlso"]["Mim-cit_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_seeAlso", "Mim-cit", "Mim-cit_others"]:
        self.record[-1]["Mim-entry_seeAlso"]["Mim-cit_others"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_seeAlso", "Mim-cit", "Mim-cit_year"]:
        self.record[-1]["Mim-entry_seeAlso"]["Mim-cit_year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_text", "Mim-text", "Mim-text_label"]:
        self.record[-1]["Mim-entry_text"][-1]["Mim-text_label"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_text", "Mim-text", "Mim-text_text"]:
        self.record[-1]["Mim-entry_text"][-1]["Mim-text_text"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields", "Mim-text", "Mim-text_label"]:
        self.record[-1]["Mim-entry_textfields"][-1]["Mim-text_label"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields", "Mim-text", "Mim-text_text"]:
        self.record[-1]["Mim-entry_textfields"][-1]["Mim-text_text"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_textfields"][-1]["Mim-text_neighbors"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_textfields"][-1]["Mim-text_neighbors"][-1]["Mim-link_uids"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_textfields", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_textfields"][-1]["Mim-text_neighbors"][-1]["Mim-link_numRelevant"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary", "Mim-text", "Mim-text_label"]:
        self.record[-1]["Mim-entry_summary"]["Mim-text_label"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary", "Mim-text", "Mim-text_text"]:
        self.record[-1]["Mim-entry_summary"]["Mim-text_text"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_summary"]["Mim-text"]["Mim-text_neighbors"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_summary"]["Mim-text"]["Mim-text_neighbors"][-1]["Mim-link_num"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summary", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_summary"]["Mim-text"]["Mim-text_neighbors"][-1]["Mim-link_numRelevant"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryAttribution", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_summaryAttribution"]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryAttribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_summaryAttribution"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryAttribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_summaryAttribution"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryAttribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_summaryAttribution"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryEditHistory", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_summaryEditHistory"]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryEditHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_summaryEditHistory"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryEditHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_summaryEditHistory"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryEditHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_summaryEditHistory"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryCreationDate", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_summaryCreationDate"]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryCreationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_summaryCreationDate"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryCreationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_summaryCreationDate"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_summaryCreationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_summaryCreationDate"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant_number"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_number"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant_name"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_name"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant_aliases", "Mim-allelic-variant_aliases_E"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_aliases"].append(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant_mutation"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_mutation"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant_description"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_description"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant_snpLinks", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_snpLinks"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant_snpLinks", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_snpLinks"][-1]["Mim-link_uids"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_allelicVariants", "Mim-allelic-variant_snpLinks", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_allelicVariants"][-1]["Mim-allelic-variant_snpLinks"][-1]["Mim-link_numRelevant"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis", "Mim-text", "Mim-text_label"]:
        self.record[-1]["Mim-entry_synopsis"]["Mim-text_label"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis", "Mim-text", "Mim-text_text"]:
        self.record[-1]["Mim-entry_synopsis"]["Mim-text_text"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_synopsis"]["Mim-text"]["Mim-text_neighbors"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_synopsis"]["Mim-text"]["Mim-text_neighbors"][-1]["Mim-link_num"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsis", "Mim-text", "Mim-text_neighbors", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_synopsis"]["Mim-text"]["Mim-text_neighbors"][-1]["Mim-link_numRelevant"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisAttribution", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_synopsisAttribution"]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisAttribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_synopsisAttribution"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisAttribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_synopsisAttribution"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisAttribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_synopsisAttribution"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisEditHistory", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_synopsisEditHistory"]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisEditHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_synopsisEditHistory"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisEditHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_synopsisEditHistory"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisEditHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_synopsisEditHistory"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisCreationDate", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_synopsisCreationDate"]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisCreationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_synopsisCreationDate"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisCreationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_synopsisCreationDate"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_synopsisCreationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_synopsisCreationDate"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date-day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_clinicalSynopsis", "Mim-index-term", "Mim-index-term_key"]:
        self.record[-1]["Mim-entry_clinicalSynopsis"][-1]["Mim-index-term_key"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_clinicalSynopsis", "Mim-index-term", "Mim-index-term_terms", "Mim-index-term_terms_E"]:
        self.record[-1]["Mim-entry_clinicalSynopsis"][-1]["Mim-index-term_terms"].append(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_editHistory", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_editHistory"][-1]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_editHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_editHistory"][-1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_editHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_editHistory"][-1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_editHistory", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_editHistory"][-1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_creationDate", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_creationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_creationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_creationDate", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_creationDate"]["Mim-edit-item"]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_number"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_number"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_origNumber"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_origNumber"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_authors", "Mim-author", "Mim-author_name"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_authors"][-1]["Mim-author_name"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_authors", "Mim-author", "Mim-author_index"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_authors"][-1]["Mim-author_index"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_primaryAuthor"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_primaryAuthor"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_otherAuthors"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_otherAuthors"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_citationTitle"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_citationTitle"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_citationType"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_citationType"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_bookTitle"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_bookTitle"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_editors"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_editors"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_volume"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_volume"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_edition"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_edition"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_journal"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_journal"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_series"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_series"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_publisher"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_publisher"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_place"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_place"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_commNote"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_commNote"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pubDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pubDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pubDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pubDate"]["Mim-date"]["Mim-date_day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pages", "Mim-page", "Mim-page_from"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pages"]["Mim-page"]["Mim-page_from"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pages", "Mim-page", "Mim-page_to"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pages"]["Mim-page"]["Mim-page_to"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_miscInfo"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_miscInfo"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_references", "Mim-reference", "Mim-reference_pubmedUID"]:
        self.record[-1]["Mim-entry_references"][-1]["Mim-reference_pubmedUID"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_attribution", "Mim-edit-item", "Mim-edit-item_author"]:
        self.record[-1]["Mim-entry_attribution"][-1]["Mim-edit-item_author"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_attribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_year"]:
        self.record[-1]["Mim-entry_attribution"][-1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_year"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_attribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_month"]:
        self.record[-1]["Mim-entry_attribution"][-1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_month"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_attribution", "Mim-edit-item", "Mim-edit-item_modDate", "Mim-date", "Mim-date_day"]:
        self.record[-1]["Mim-entry_attribution"][-1]["Mim-edit-item_modDate"]["Mim-date"]["Mim-date_day"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_numGeneMaps"]:
        self.record[-1]["Mim-entry_numGeneMaps"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_medlineLinks", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_medlineLinks"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_medlineLinks", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_medlineLinks"][-1]["Mim-link_uids"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_medlineLinks", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_medlineLinks"][-1]["Mim-link_numRelevant"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_proteinLinks", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_proteinLinks"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_proteinLinks", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_proteinLinks"][-1]["Mim-link_uids"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_proteinLinks", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_proteinLinks"][-1]["Mim-link_numRelevant"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_nucleotideLinks", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_nucleotideLinks"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_nucleotideLinks", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_nucleotideLinks"][-1]["Mim-link_uids"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_nucleotideLinks", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_nucleotideLinks"][-1]["Mim-link_numRelevant"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_structureLinks", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_structureLinks"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_structureLinks", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_structureLinks"][-1]["Mim-link_uids"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_structureLinks", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_structureLinks"][-1]["Mim-link_numRelevant"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_genomeLinks", "Mim-link", "Mim-link_num"]:
        self.record[-1]["Mim-entry_genomeLinks"][-1]["Mim-link_num"] = int(self.content)
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_genomeLinks", "Mim-link", "Mim-link_uids"]:
        self.record[-1]["Mim-entry_genomeLinks"][-1]["Mim-link_uids"] = self.content
    elif self.element==["Mim-entries", "Mim-entry", "Mim-entry_genomeLinks", "Mim-link", "Mim-link_numRelevant"]:
        self.record[-1]["Mim-entry_genomeLinks"][-1]["Mim-link_numRelevant"] = int(self.content)
